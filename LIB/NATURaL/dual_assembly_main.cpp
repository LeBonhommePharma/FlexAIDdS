// dual_assembly_main.cpp — Cotranslational docking CLI entry point.
//
// Usage:
//   dual_assembly --target-pdb <protofibril.pdb>
//                 --sequence   <FASTA-1-letter or @file.fasta>
//                 [--monomer-pdb <monomer.pdb>]
//                 [--checkpoint-interval 10]
//                 [--sim-c-interval 5]
//                 [--monomer-conc-M 1e-6]
//                 [--temperature 310.15]
//                 [--threads 6]
//                 [--no-reciprocal-controls]
//                 [--no-sim-c]
//                 [--output-csv cotranslational_trajectory.csv]
//                 [--nascent-pdb-dir .]
//                 [--synthetic]    # use built-in synthetic GA backend (MVP default)
//                 [--real-ga]      # reserved; exits until FlexAID GA callback lands
//
// The MVP ships a built-in synthetic GA backend (Gaussian pose distributions whose
// width narrows as L_k grows) so the full pipeline can be exercised end-to-end on
// any platform. The `--real-ga` flag intentionally fails closed until the FlexAID
// GA callback wiring lands, so production users cannot mistake synthetic evidence
// for real docking.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include "DualAssemblyRunner.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace {

// ─── argv helpers ────────────────────────────────────────────────────────────
struct Args {
    std::string target_pdb;
    std::string monomer_pdb;
    std::string sequence_or_at_file;
    int    checkpoint_interval = 10;
    int    sim_c_interval      = 5;
    double monomer_conc_M      = 1.0e-6;
    double temperature         = 310.15;
    int    threads             = 6;
    bool   reciprocal_controls = true;
    bool   sim_c_enabled       = true;
    bool   synthetic           = true;   // MVP default
    std::string output_csv     = "cotranslational_trajectory.csv";
    std::string nascent_pdb_dir = ".";
};

std::string load_sequence(const std::string& spec) {
    if (spec.empty() || spec[0] != '@') return spec;
    std::ifstream in(spec.substr(1));
    if (!in) throw std::runtime_error("cannot open FASTA: " + spec.substr(1));
    std::string line, seq;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line)
            if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) seq.push_back(c);
    }
    return seq;
}

void usage() {
    std::cerr <<
        "usage: dual_assembly --target-pdb <pdb> --sequence <FASTA|@file> [opts]\n"
        "see file header for full option list\n";
}

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&](const char* name) -> std::string {
            if (i + 1 >= argc) { usage(); std::exit(2); }
            (void)name;
            return argv[++i];
        };
        if      (k == "--target-pdb")          a.target_pdb = next(k.c_str());
        else if (k == "--monomer-pdb")         a.monomer_pdb = next(k.c_str());
        else if (k == "--sequence")            a.sequence_or_at_file = next(k.c_str());
        else if (k == "--checkpoint-interval") a.checkpoint_interval = std::atoi(next(k.c_str()).c_str());
        else if (k == "--sim-c-interval")      a.sim_c_interval      = std::atoi(next(k.c_str()).c_str());
        else if (k == "--monomer-conc-M")      a.monomer_conc_M      = std::atof(next(k.c_str()).c_str());
        else if (k == "--temperature")         a.temperature         = std::atof(next(k.c_str()).c_str());
        else if (k == "--threads")             a.threads             = std::atoi(next(k.c_str()).c_str());
        else if (k == "--no-reciprocal-controls") a.reciprocal_controls = false;
        else if (k == "--no-sim-c")            a.sim_c_enabled       = false;
        else if (k == "--synthetic")           a.synthetic           = true;
        else if (k == "--real-ga")             a.synthetic           = false;
        else if (k == "--output-csv")          a.output_csv          = next(k.c_str());
        else if (k == "--nascent-pdb-dir")     a.nascent_pdb_dir     = next(k.c_str());
        else if (k == "--help" || k == "-h") { usage(); std::exit(0); }
        else {
            std::cerr << "unknown flag: " << k << "\n";
            usage();
            std::exit(2);
        }
    }
    return a;
}

// ─── Synthetic GA backend ────────────────────────────────────────────────────
// Generates a Gaussian-spread pose ensemble whose width narrows as L_k grows. The
// resulting Shannon entropy crosses through the soft and hard thresholds at
// physiologically plausible chain lengths, exercising every regime of the T/L
// discriminator. This is for end-to-end smoke testing, NOT for production docking.
natural::GAResult synthetic_sim(double temperature_K,
                                 int    n_poses,
                                 double width_A,
                                 double dG_kcal,
                                 unsigned seed)
{
    statmech::StatMechEngine engine(temperature_K);
    std::mt19937 rng(seed);
    std::normal_distribution<double> rmsd_dist(0.0, std::max(1e-3, width_A));
    // We seed the engine with energies drawn from a Gaussian centred at dG_kcal.
    // The Helmholtz F = -kT ln Z then approximates dG_kcal for large N.
    std::normal_distribution<double> e_dist(dG_kcal, 1.0);
    std::vector<double> rmsds;
    rmsds.reserve(n_poses);
    for (int i = 0; i < n_poses; ++i) {
        engine.add_sample(e_dist(rng), 1.0);
        rmsds.push_back(std::abs(rmsd_dist(rng)));
    }
    return {std::move(engine), std::move(rmsds)};
}

natural::SimAFn make_synthetic_sim_a() {
    return [](const std::string&, const std::string&, int L_k, double T_K) {
        // Wider Gaussian when chain is short (still floppy); narrows past L_k ≈ 60.
        const double width = std::max(0.5, 8.0 - 0.07 * L_k);
        const double dG    = -2.0 - 0.04 * L_k;
        return synthetic_sim(T_K, 256, width, dG, /*seed=*/0xA5A5u + L_k);
    };
}

natural::SimBFn make_synthetic_sim_b() {
    return [](const std::string&, const std::string&, int L_k, double T_K) {
        const double width = std::max(0.3, 4.0 - 0.02 * L_k);
        const double dG    = -1.0 - 0.02 * L_k;
        return synthetic_sim(T_K, 128, width, dG, /*seed=*/0xB5B5u + L_k);
    };
}

natural::SimCFn make_synthetic_sim_c() {
    return [](const std::string&, const std::string&, double T_K) {
        // Slightly favourable monomer addition.
        return synthetic_sim(T_K, 128, 1.5, -3.0, /*seed=*/0xC5C5u);
    };
}

// ─── Truncation helper (extended geometry MVP) ───────────────────────────────
// MVP: write a marker PDB with one CA atom per residue along the +X axis at 3.8 Å
// spacing. Real geometry comes later via a Python/ESMFold helper.
natural::TruncateFn make_synthetic_truncate() {
    return [](const std::string& seq, int L_k, const std::string& out_dir) -> std::string {
        const int  L  = std::min(L_k, static_cast<int>(seq.size()));
        fs::create_directories(out_dir);
        const std::string path = (fs::path(out_dir) / ("nascent_L" + std::to_string(L_k) + ".pdb")).string();
        std::ofstream f(path, std::ios::trunc);
        if (!f) throw std::runtime_error("cannot write " + path);
        f << "REMARK 100  Synthetic nascent chain at L_k = " << L_k << "\n";
        for (int i = 0; i < L; ++i) {
            const char aa = seq[i];
            char buf[96];
            std::snprintf(buf, sizeof(buf),
                "ATOM  %5d  CA  %-3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
                i + 1, "ALA", i + 1, 3.8 * i, 0.0, 0.0);
            (void)aa;
            f << buf;
        }
        f << "TER\nEND\n";
        return path;
    };
}

} // namespace

int main(int argc, char** argv) {
    try {
        Args a = parse_args(argc, argv);
        if (!a.synthetic) {
            std::cerr << "[dual_assembly] --real-ga is not implemented yet; refusing to run synthetic backend\n";
            return 2;
        }
        if (a.target_pdb.empty() || a.sequence_or_at_file.empty()) {
            usage();
            return 2;
        }

        natural::DualAssemblyConfig cfg;
        cfg.protofibril_pdb              = a.target_pdb;
        cfg.monomer_pdb                  = a.monomer_pdb;
        cfg.sequence_fasta               = load_sequence(a.sequence_or_at_file);
        cfg.checkpoint_interval          = a.checkpoint_interval;
        cfg.sim_c_interval               = a.sim_c_interval;
        cfg.monomer_conc_M               = a.monomer_conc_M;
        cfg.temperature_K                = a.temperature;
        cfg.n_threads                    = a.threads;
        cfg.include_reciprocal_controls  = a.reciprocal_controls;
        cfg.sim_c_enabled                = a.sim_c_enabled && !a.monomer_pdb.empty();
        cfg.output_csv                   = a.output_csv;
        cfg.nascent_pdb_dir              = a.nascent_pdb_dir;

        natural::SimAFn sim_a = make_synthetic_sim_a();
        natural::SimBFn sim_b = make_synthetic_sim_b();
        natural::SimCFn sim_c = cfg.sim_c_enabled ? make_synthetic_sim_c() : nullptr;
        natural::TruncateFn truncate = make_synthetic_truncate();

        natural::DualAssemblyRunner runner(std::move(cfg),
                                            std::move(sim_a),
                                            std::move(sim_b),
                                            std::move(sim_c),
                                            std::move(truncate));
        auto history = runner.run();
        std::cerr << "[dual_assembly] wrote " << history.size()
                  << " checkpoints to " << a.output_csv << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "dual_assembly: " << e.what() << "\n";
        return 1;
    }
}
