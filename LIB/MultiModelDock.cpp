// MultiModelDock.cpp — CCBM orchestrator implementation
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "MultiModelDock.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <numeric>
#include <stdexcept>

namespace ccbm {

// ─── File type detection ─────────────────────────────────────────────────────

bool MultiModelDock::is_cif_file(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return false;
    std::string ext = path.substr(dot);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return (ext == ".cif" || ext == ".mmcif");
}

// ─── Constructor ─────────────────────────────────────────────────────────────

MultiModelDock::MultiModelDock(const std::string& receptor_path)
    : receptor_path_(receptor_path),
      is_cif_(is_cif_file(receptor_path))
{
}

// ─── Load models ─────────────────────────────────────────────────────────────

int MultiModelDock::load_models(FA_Global* FA, atom** atoms, resid** residue) {
    int n = 0;

    if (is_cif_) {
        n = read_multi_model_cif(FA, atoms, residue, receptor_path_.c_str());
    } else {
        n = read_multi_model_pdb(FA, atoms, residue, receptor_path_.c_str());
    }

    if (n <= 0) {
        // Fall back to single model
        n = 1;
        FA->n_models = 1;
        FA->model_strain.resize(1, 0.0);
    }

    // Build internal ReceptorModel list from FA->model_coords
    models_.resize(FA->n_models);
    for (int i = 0; i < FA->n_models; ++i) {
        models_[i].model_num = i + 1;
        models_[i].sequential_idx = i;
        models_[i].strain_energy = (i < static_cast<int>(FA->model_strain.size()))
                                    ? FA->model_strain[i] : 0.0;
        if (i < static_cast<int>(FA->model_coords.size()))
            models_[i].atom_coords = FA->model_coords[i];
    }

    return FA->n_models;
}

// ─── Set strain energies ─────────────────────────────────────────────────────

void MultiModelDock::set_strain_energies(const std::vector<double>& energies) {
    for (size_t i = 0; i < models_.size() && i < energies.size(); ++i) {
        models_[i].strain_energy = energies[i];
    }
}

void MultiModelDock::set_strain_energies_from_file(const std::string& strain_file) {
    std::ifstream ifs(strain_file);
    if (!ifs.is_open()) {
        fprintf(stderr, "WARNING: Cannot open strain file: %s — using default strains (0.0)\n",
                strain_file.c_str());
        return;
    }

    std::vector<double> energies;
    double val;
    while (ifs >> val) {
        energies.push_back(val);
    }

    set_strain_energies(energies);
}

// ─── Harmonic strain approximation ──────────────────────────────────────────

void MultiModelDock::compute_harmonic_strain(double k) {
    if (models_.empty()) return;

    // Reference: first model
    const auto& ref_coords = models_[0].atom_coords;
    int n_atoms = static_cast<int>(ref_coords.size()) / 3;

    models_[0].strain_energy = 0.0;  // reference has zero strain

    for (size_t m = 1; m < models_.size(); ++m) {
        const auto& coords = models_[m].atom_coords;
        int n = std::min(n_atoms, static_cast<int>(coords.size()) / 3);

        double sum_sq = 0.0;
        for (int a = 0; a < n; ++a) {
            float dx = coords[a*3+0] - ref_coords[a*3+0];
            float dy = coords[a*3+1] - ref_coords[a*3+1];
            float dz = coords[a*3+2] - ref_coords[a*3+2];
            sum_sq += dx*dx + dy*dy + dz*dz;
        }

        // RMSD² = sum_sq / n_atoms
        double rmsd_sq = (n > 0) ? sum_sq / static_cast<double>(n) : 0.0;
        models_[m].strain_energy = 0.5 * k * rmsd_sq;
    }
}

// ─── Configure FA_Global ─────────────────────────────────────────────────────

void MultiModelDock::configure_fa(FA_Global* FA) {
    FA->multi_model = true;
    FA->n_models = static_cast<int>(models_.size());
    FA->model_coords.resize(FA->n_models);
    FA->model_strain.resize(FA->n_models);

    for (int i = 0; i < FA->n_models; ++i) {
        FA->model_coords[i] = models_[i].atom_coords;
        FA->model_strain[i] = models_[i].strain_energy;
    }
}

// ─── Generate report ─────────────────────────────────────────────────────────

CCBMReport MultiModelDock::generate_report(const BindingPopulation& population) const {
    CCBMReport report;
    report.receptor_file = receptor_path_;
    report.n_models = static_cast<int>(models_.size());
    report.n_total_poses = 0;

    // Collect all poses across all binding modes
    std::vector<const Pose*> all_poses;
    const auto& modes = population.get_binding_modes();
    for (const auto& mode : modes) {
        for (int i = 0; i < mode.get_BindingMode_size(); ++i) {
            all_poses.push_back(&mode.get_pose(i));
            report.n_total_poses++;
        }
    }

    // Per-model statistics
    report.per_model.resize(report.n_models);
    for (int m = 0; m < report.n_models; ++m) {
        report.per_model[m].model_index = m;
        report.per_model[m].best_score = 1e30;
        report.per_model[m].mean_score = 0.0;
        report.per_model[m].n_poses = 0;
        report.per_model[m].population = 0.0;
    }

    for (const auto* p : all_poses) {
        int m = p->model_index;
        if (m < 0 || m >= report.n_models) m = 0;
        report.per_model[m].n_poses++;
        double e = p->total_energy();
        if (e < report.per_model[m].best_score)
            report.per_model[m].best_score = e;
    }

    // Compute ensemble thermodynamics from the first (largest) binding mode
    // that has conformer-coupled poses, or aggregate all
    if (!modes.empty()) {
        // Use the top binding mode for entropy decomposition
        const BindingMode& top_mode = modes[0];
        report.entropy_decomp = top_mode.decompose_entropy();

        auto pops = top_mode.conformer_populations();
        for (int m = 0; m < report.n_models && m < static_cast<int>(pops.size()); ++m) {
            report.per_model[m].population = pops[m];
        }

        // Global ensemble from StatMechEngine
        statmech::StatMechEngine global = population.get_global_ensemble();
        auto td = global.compute();
        report.thermo.F = td.free_energy;
        report.thermo.S = td.entropy;
        report.thermo.H = td.mean_energy;
        report.thermo.Cv = td.heat_capacity;
        report.thermo.conformer_populations = pops;
        report.thermo.S_receptor = report.entropy_decomp.S_receptor;
        report.thermo.S_ligand = report.entropy_decomp.S_ligand;
        report.thermo.I_mutual = report.entropy_decomp.I_mutual;
        report.thermo.S_vibrational = report.entropy_decomp.S_vibrational;
    }

    return report;
}

// ─── Print report ────────────────────────────────────────────────────────────

void MultiModelDock::print_report(const CCBMReport& report) {
    printf("\n═══════════════════════════════════════════════════════════════\n");
    printf("  CCBM: Conformer-Coupled Binding Mode Report\n");
    printf("═══════════════════════════════════════════════════════════════\n");
    printf("Receptor:    %s\n", report.receptor_file.c_str());
    printf("Models:      %d\n", report.n_models);
    printf("Total poses: %d\n", report.n_total_poses);
    // Claim firewall + log-parser firewall.
    //
    // These numbers come from a contact-function optimizer ensemble, so they are
    // proxy_only diagnostics; they are never a calibrated kcal/mol thermodynamic
    // quantity (see LIB/statmech.h ScientificProvenance).
    //
    // The label text below is ALSO load-bearing. LIB/DatasetThermoLog.h scrapes
    // this same stdout stream and only gaboom's post-GA ensemble block may reach
    // result.csv. Its predicates are:
    //   free energy : "F-like proxy" AND "F~"   (or "Helmholtz free energy")
    //   mean energy : "Mean CF"                 (or "Mean energy")
    //   entropy     : "S-like diagnostic"       (or "Entropy (conf)")
    // and a matched line is only consumed if it also contains '='.
    // CCBM must stay out of result.csv, so these lines deliberately
    //   (a) omit the "F~" marker, the "Mean CF" stem and the "S-like diagnostic"
    //       spelling, and
    //   (b) use ':' rather than '=' as the separator.
    // Do not "tidy" either property without re-reading DatasetThermoLog.h.
    printf("\n── Global CF-proxy ensemble diagnostics ──\n");
    printf("  claim_validity:      proxy_only\n");
    printf("  F-like proxy:        %10.4f [legacy transform]\n", report.thermo.F);
    printf("  CF-proxy mean <CF>:  %10.4f [CF units]\n", report.thermo.H);
    printf("  S-like proxy:        %10.6f [proxy scale/K]\n", report.thermo.S);
    printf("  C_v-like proxy:      %10.6f [proxy scale/K]\n", report.thermo.Cv);
    printf("\n── Model-scale decomposition (diagnostic) ──\n");
    printf("  S_total:       %10.6f [model scale]\n", report.entropy_decomp.S_total);
    printf("  S_ligand:      %10.6f [model scale]\n", report.entropy_decomp.S_ligand);
    printf("  S_receptor:    %10.6f [model scale]\n", report.entropy_decomp.S_receptor);
    printf("  I(L;R):        %10.6f [model scale]\n", report.entropy_decomp.I_mutual);
    printf("  S_vibrational: %10.6f [model scale]\n", report.entropy_decomp.S_vibrational);
    printf("\n── Per-Model Results ──\n");
    printf("  Model | Poses | Best Score | Population\n");
    printf("  ------|-------|------------|----------\n");
    for (const auto& m : report.per_model) {
        printf("  %5d | %5d | %10.4f | %8.4f\n",
               m.model_index, m.n_poses, m.best_score, m.population);
    }
    printf("═══════════════════════════════════════════════════════════════\n\n");
}

// ─── Write report ────────────────────────────────────────────────────────────

void MultiModelDock::write_report(const CCBMReport& report, const std::string& outpath) {
    FILE* fp = fopen(outpath.c_str(), "w");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot write CCBM report to %s\n", outpath.c_str());
        return;
    }

    fprintf(fp, "# CCBM Conformer-Coupled Binding Mode Report\n");
    fprintf(fp, "receptor: %s\n", report.receptor_file.c_str());
    fprintf(fp, "n_models: %d\n", report.n_models);
    fprintf(fp, "n_total_poses: %d\n", report.n_total_poses);
    // Same firewall as print_report: CF-proxy optimizer samples carry no
    // calibrated energy domain, so the record is stamped proxy_only and the
    // F/H/S/Cv fields below are diagnostics, not kcal/mol thermodynamics.
    fprintf(fp, "\n# Global CF-proxy ensemble diagnostics\n");
    fprintf(fp, "claim_validity: proxy_only\n");
    fprintf(fp, "energy_domain: cf_arbitrary_units\n");
    fprintf(fp, "ensemble_measure: optimizer_samples\n");
    fprintf(fp, "F: %.6f\n", report.thermo.F);
    fprintf(fp, "H: %.6f\n", report.thermo.H);
    fprintf(fp, "S: %.8f\n", report.thermo.S);
    fprintf(fp, "Cv: %.8f\n", report.thermo.Cv);
    fprintf(fp, "\n# Entropy Decomposition\n");
    fprintf(fp, "S_total: %.8f\n", report.entropy_decomp.S_total);
    fprintf(fp, "S_ligand: %.8f\n", report.entropy_decomp.S_ligand);
    fprintf(fp, "S_receptor: %.8f\n", report.entropy_decomp.S_receptor);
    fprintf(fp, "I_mutual: %.8f\n", report.entropy_decomp.I_mutual);
    fprintf(fp, "S_vibrational: %.8f\n", report.entropy_decomp.S_vibrational);
    fprintf(fp, "\n# Per-Model Results (model_index, n_poses, best_score, population)\n");
    for (const auto& m : report.per_model) {
        fprintf(fp, "%d\t%d\t%.6f\t%.6f\n",
                m.model_index, m.n_poses, m.best_score, m.population);
    }

    fclose(fp);
    printf("CCBM report written to %s\n", outpath.c_str());
}

}  // namespace ccbm
