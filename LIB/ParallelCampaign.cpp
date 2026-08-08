// ParallelCampaign.cpp — GPU/SIMD/OpenMP-accelerated virtual screening
//
// C++26 throughout: std::jthread, std::atomic, std::span, std::format (where
// available), std::ranges, constexpr math helpers, structured bindings.
//
// Three-level parallelism:
//   L1  OpenMP task loop over ligands       (#pragma omp parallel for schedule(dynamic))
//   L2  Model-parallel docking per ligand   (nested OpenMP or sequential)
//   L3  GPU batch GA evaluation per dock    (CUDA/ROCm/Metal via HardwareDispatch)
//
// Eigen used for:
//   - Boltzmann weight vectors (Eigen::ArrayXd for vectorized exp/log)
//   - Ensemble consensus log-sum-exp
//   - Reference entropy matrix operations
//
// AVX-512/AVX2 used for:
//   - log-sum-exp in Boltzmann weighting (via HardwareDispatch::log_sum_exp_dispatch)
//   - Distance computations in coordinate builder
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.

#include "ParallelCampaign.h"
#include "LibrarySplitter.h"
#include "ReferenceEntropy.h"
#include "hardware_detect.h"
#include "ProcessLigand/ProcessLigand.h"
#include "ProcessLigand/CoordBuilder.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <ranges>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace campaign {

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

// ─── Eigen-accelerated Boltzmann consensus ───────────────────────────────────

#if 0  // currently unused — reserved for future Eigen-accelerated consensusstatic double boltzmann_consensus_eigen(std::span<const double> dG_values,
                                         double temperature_K) {
    const int N = static_cast<int>(dG_values.size());
    if (N == 0) return 0.0;
    if (N == 1) return dG_values[0];

    const double beta = 1.0 / (reference_entropy::kB_kcal * temperature_K);

    // Eigen vectorized path: exp, log, sum all auto-vectorize to AVX/SSE
    Eigen::Map<const Eigen::ArrayXd> dG(dG_values.data(), N);
    Eigen::ArrayXd exponents = -beta * dG;

    // log-sum-exp with Eigen
    const double max_exp = exponents.maxCoeff();
    const double lse = max_exp + std::log((exponents - max_exp).exp().sum() / N);
    return -lse / beta;
}
#endif  // boltzmann_consensus_eigen

// NON-PHYSICAL SURROGATE — placeholder only.
//
// This is a hardcoded linear function of ligand composition (heavy atoms,
// rotatable bonds, rings) plus a model-index and temperature nudge. It never
// touches the receptor, the pose, or the scoring function, so it carries no
// energy domain at all: it is not a docking score, not a contact-function
// value, and emphatically not a binding free energy in kcal/mol.
//
// It exists so the three-level parallel campaign skeleton can be exercised
// while the GA is still reachable only through top.cpp's globals (see the
// "Level 2" comment in run_campaign). Every consumer of the values it produces
// must present them as `surrogate_*` on a model scale — never as dG, ΔG,
// affinity, Kd, Ki or potency (see LIB/statmech.h ScientificProvenance:
// no calibrated energy domain and no sha256 receipt ⇒ ProxyOnly).
//
// When a real per-model dock is wired in, replace this function AND re-audit
// the labels in run_campaign/write_results_csv/write_top_hits below.
static double surrogate_model_dock_score(const LigandResult& lr,
                                         int model_idx,
                                         double temperature_K) {
    const double base = -0.04 * static_cast<double>(lr.n_atoms)
                      - 0.12 * static_cast<double>(lr.n_rotatable)
                      - 0.03 * static_cast<double>(lr.n_rings);
    const double model_offset = -0.015 * static_cast<double>(model_idx);
    const double thermal = 0.001 * ((temperature_K - 300.0) / 10.0);
    return base + model_offset + thermal;
}

static std::string element_symbol_from_z(int z) {
    switch (z) {
        case 1: return "H";
        case 6: return "C";
        case 7: return "N";
        case 8: return "O";
        case 9: return "F";
        case 15: return "P";
        case 16: return "S";
        case 17: return "Cl";
        case 35: return "Br";
        case 53: return "I";
        default: return "X";
    }
}

// ─── Auto-configure from inputs ──────────────────────────────────────────────

CampaignConfig auto_configure(
    const std::string& receptor_path,
    const std::string& ligand_path,
    const std::string& config_json,
    const std::string& output_prefix,
    bool rigid,
    bool folded)
{
    CampaignConfig cfg;
    cfg.receptor_path = receptor_path;
    cfg.ligand_path   = ligand_path;
    cfg.config_json   = config_json;
    cfg.output_prefix = output_prefix;
    cfg.use_rigid     = rigid;
    cfg.use_folded    = folded;

    // Detect receptor model count
    cfg.receptor_models = library::detect_library_size(receptor_path);
    if (cfg.receptor_models > 1)
        printf("Multi-model receptor: %d conformers detected\n", cfg.receptor_models);

    // Detect ligand library size
    cfg.ligand_count = library::detect_library_size(ligand_path);
    printf("Ligand library: %d compounds detected\n", cfg.ligand_count);

    // Auto-enable reference entropy for multi-model
    cfg.compute_ref_entropy = (cfg.receptor_models > 1);

    // Determine optimal thread counts
    int hw_threads = 1;
#ifdef _OPENMP
    hw_threads = omp_get_max_threads();
#else
    hw_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (hw_threads < 1) hw_threads = 1;
#endif

    // Detect GPU availability
    const auto& hw_info = flexaids::detect_hardware();
    bool has_gpu = hw_info.has_cuda || hw_info.has_metal;

    if (has_gpu) {
        cfg.max_ligand_threads = std::max(1, hw_threads - 2);
        cfg.gpu_per_ligand = false;

        if (hw_info.has_cuda)
            printf("GPU: CUDA (%s)\n", hw_info.cuda_device_name.c_str());
        else if (hw_info.has_metal)
            printf("GPU: Metal (%s)\n", hw_info.metal_gpu_name.c_str());
    } else {
        // CPU-only: all cores for ligand parallelism
        // Each ligand dock uses OpenMP internally for GA evaluation
        // Avoid oversubscription: fewer ligand threads, more per-dock threads
        if (cfg.ligand_count > hw_threads * 2) {
            cfg.max_ligand_threads = std::max(1, hw_threads / 2);
        } else {
            cfg.max_ligand_threads = std::min(cfg.ligand_count, hw_threads);
        }

        if (hw_info.has_avx512)
            printf("CPU: AVX-512 (%d threads)\n", hw_info.openmp_max_threads);
        else if (hw_info.has_avx2)
            printf("CPU: AVX2 (%d threads)\n", hw_info.openmp_max_threads);
    }

    printf("Campaign: %d ligands × %d models = %d docks, %d threads\n",
           cfg.ligand_count, cfg.receptor_models,
           cfg.ligand_count * cfg.receptor_models,
           cfg.max_ligand_threads);

    // Results CSV path
    cfg.results_csv = output_prefix + "_results.csv";

    return cfg;
}

// ─── Run campaign ────────────────────────────────────────────────────────────

CampaignSummary run_campaign(
    const CampaignConfig& config,
    ProgressCallback progress)
{
    auto t_start = Clock::now();

    CampaignSummary summary{};
    summary.total_ligands   = config.ligand_count;
    summary.receptor_models = config.receptor_models;
    summary.total_docks     = config.ligand_count * config.receptor_models;
    summary.cpu_threads_used = config.max_ligand_threads;

    // Split ligand library
    auto lig_lib = library::split_library(config.ligand_path);

    // Split receptor if multi-model
    library::LibraryInfo rec_lib;
    if (config.receptor_models > 1) {
        rec_lib = library::split_library(config.receptor_path);
    } else {
        library::LigandEntry e;
        e.path = config.receptor_path;
        e.name = "receptor";
        e.format = "pdb";
        e.is_temp = false;
        rec_lib.ligands.push_back(e);
        rec_lib.total = 1;
    }

    // Pre-allocate results
    std::vector<LigandResult> results(lig_lib.total);
    std::atomic<int> completed{0};
    std::atomic<int> successful{0};
    std::atomic<int> failed{0};
    std::mutex results_mutex;

    // CSV header (stream as we go)
    std::ofstream csv;
    if (config.stream_results && !config.results_csv.empty()) {
        csv.open(config.results_csv);
        // Score columns are named `surrogate_*`, not dG_*: they are fed by
        // surrogate_model_dock_score(), a composition-only placeholder with no
        // energy domain. Renaming them back to dG_* would re-assert a binding
        // free-energy claim this pipeline cannot support.
        csv << "rank,name,surrogate_consensus,surrogate_best,surrogate_corrected,"
            << "surrogate_mean,surrogate_stddev,"
            << "best_model,n_atoms,n_rotatable,n_rings,mw,time_sec,status\n";
    }

    // ─── Level 1: Parallel over ligands ──────────────────────────────────────
    const int n_lig = lig_lib.total;

#ifdef _OPENMP
    omp_set_nested(1);  // enable nested parallelism for L2
    #pragma omp parallel for \
        num_threads(config.max_ligand_threads) \
        schedule(dynamic, 1) \
        default(none) \
        shared(lig_lib, rec_lib, results, completed, successful, failed, \
               csv, results_mutex, config, progress, n_lig)
#endif
    for (int li = 0; li < n_lig; li++) {
        auto t_lig_start = Clock::now();
        LigandResult& lr = results[li];
        lr.ligand_idx = li;
        lr.name = lig_lib.ligands[li].name;

        // ── Process ligand through ProcessLigand pipeline ────────────────
        bonmol::ProcessOptions pl_opts;
        pl_opts.input  = lig_lib.ligands[li].path;
        if (lig_lib.ligands[li].format == "smiles")
            pl_opts.format = bonmol::InputFormat::SMILES;
        else if (lig_lib.ligands[li].format == "sdf")
            pl_opts.format = bonmol::InputFormat::SDF;
        else
            pl_opts.format = bonmol::InputFormat::MOL2;

        bonmol::ProcessLigand pl;
        auto pl_result = pl.run(pl_opts);

        if (!pl_result.success) {
            lr.success = false;
            lr.error = pl_result.error;
            failed.fetch_add(1, std::memory_order_relaxed);
            completed.fetch_add(1, std::memory_order_relaxed);
            continue;
        }

        // Build 3D if SMILES
        if (lig_lib.ligands[li].format == "smiles") {
            bonmol::CoordBuilderOptions cb_opts;
            if (!bonmol::build_3d_coords(pl_result.mol, cb_opts)) {
                lr.success = false;
                lr.error = "3D coordinate generation failed";
                failed.fetch_add(1, std::memory_order_relaxed);
                completed.fetch_add(1, std::memory_order_relaxed);
                continue;
            }
        }

        lr.n_atoms      = pl_result.num_heavy_atoms;
        lr.n_rotatable  = pl_result.num_rot_bonds;
        lr.n_rings      = pl_result.num_rings;
        lr.molecular_weight = pl_result.molecular_weight;
        lr.pose_xyz.reserve(pl_result.mol.atoms.size());
        lr.pose_atomic_numbers.reserve(pl_result.mol.atoms.size());
        for (int ai = 0; ai < pl_result.mol.num_atoms(); ++ai) {
            lr.pose_xyz.push_back({
                pl_result.mol.coords(0, ai),
                pl_result.mol.coords(1, ai),
                pl_result.mol.coords(2, ai)
            });
            lr.pose_atomic_numbers.push_back(
                static_cast<int>(pl_result.mol.atoms[ai].element));
        }

        // ── Level 2: Dock against each receptor model ────────────────────
        // In a full implementation, each model dock calls the GA engine.
        // For now, we record the infrastructure — the actual GA call requires
        // the full FA/GB/VC setup which is managed by the top.cpp pipeline.
        // This framework enables parallelism when the GA is refactored to
        // accept receptor/ligand as parameters rather than globals.

        const int n_models = rec_lib.total;
        lr.per_model_dG.resize(n_models, 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n_models > 1) num_threads(config.max_model_threads > 0 ? config.max_model_threads : 1)
#endif
        for (int mi = 0; mi < n_models; mi++) {
            lr.per_model_dG[mi] = surrogate_model_dock_score(lr, mi, config.temperature_K);
        }

        // ── Ensemble consensus ───────────────────────────────────────────
        auto consensus = reference_entropy::compute_ensemble_consensus(
            lr.per_model_dG, config.temperature_K);
        lr.dG_consensus = consensus.dG_consensus;
        lr.dG_best      = consensus.dG_best;
        lr.dG_mean      = consensus.dG_mean;
        lr.dG_stddev    = consensus.dG_stddev;
        lr.best_model   = consensus.best_model_idx;

        // ── Reference entropy correction ─────────────────────────────────
        if (config.compute_ref_entropy && n_models > 1) {
            auto ref = reference_entropy::compute_reference_correction(
                lr.per_model_dG, lr.n_rotatable, config.temperature_K);
            lr.ref_entropy_correction = ref.T_dS_total;
            lr.dG_corrected = lr.dG_consensus + ref.T_dS_total;
        } else {
            lr.ref_entropy_correction = 0.0;
            lr.dG_corrected = lr.dG_consensus;
        }

        lr.success = true;
        auto t_lig_end = Clock::now();
        lr.dock_time_sec = std::chrono::duration<double>(t_lig_end - t_lig_start).count();

        successful.fetch_add(1, std::memory_order_relaxed);
        int done = completed.fetch_add(1, std::memory_order_relaxed) + 1;

        // Stream result to CSV
        if (config.stream_results && csv.is_open()) {
            std::lock_guard<std::mutex> lock(results_mutex);
            csv << done << ","
                << lr.name << ","
                << lr.dG_consensus << ","
                << lr.dG_best << ","
                << lr.dG_corrected << ","
                << lr.dG_mean << ","
                << lr.dG_stddev << ","
                << lr.best_model + 1 << ","
                << lr.n_atoms << ","
                << lr.n_rotatable << ","
                << lr.n_rings << ","
                << lr.molecular_weight << ","
                << lr.dock_time_sec << ","
                << (lr.success ? "OK" : lr.error)
                << "\n";
            csv.flush();
        }

        // Progress callback
        if (progress)
            progress(done, n_lig, lr);
    }

    // Finalize
    auto t_end = Clock::now();
    summary.total_time_sec = std::chrono::duration<double>(t_end - t_start).count();
    summary.successful = successful.load();
    summary.failed = failed.load();
    summary.avg_time_per_ligand = summary.total_time_sec / std::max(1, summary.total_ligands);
    summary.throughput_per_hour = 3600.0 / std::max(0.001, summary.avg_time_per_ligand);

    // Sort results by dG_corrected for top hits
    std::vector<LigandResult> sorted_results;
    sorted_results.reserve(results.size());
    for (auto& r : results) {
        if (r.success) sorted_results.push_back(r);
    }
    std::ranges::sort(sorted_results, {},
        [](const LigandResult& r) { return r.dG_corrected; });

    int top_n = std::min(100, static_cast<int>(sorted_results.size()));
    summary.top_hits.assign(sorted_results.begin(), sorted_results.begin() + top_n);

    // Detect GPU backend used
    const auto& hw_final = flexaids::detect_hardware();
    if (hw_final.has_cuda) summary.gpu_backend = "CUDA";
    else if (hw_final.has_metal) summary.gpu_backend = "Metal";
    else if (hw_final.has_avx512) summary.gpu_backend = "AVX-512 (CPU)";
    else if (hw_final.has_avx2) summary.gpu_backend = "AVX2 (CPU)";
    else summary.gpu_backend = "Scalar (CPU)";

    // Print summary
    printf("\n======= Campaign Complete =======\n");
    printf("  Ligands:     %d (%d OK, %d failed)\n",
           summary.total_ligands, summary.successful, summary.failed);
    printf("  Models:      %d\n", summary.receptor_models);
    printf("  Total docks: %d\n", summary.total_docks);
    printf("  Wall time:   %.1f sec (%.1f min)\n",
           summary.total_time_sec, summary.total_time_sec / 60.0);
    printf("  Throughput:  %.0f ligands/hour\n", summary.throughput_per_hour);
    printf("  Backend:     %s, %d CPU threads\n",
           summary.gpu_backend.c_str(), summary.cpu_threads_used);
    if (!summary.top_hits.empty()) {
        // Never "dG = ... kcal/mol": the ranking key is a non-physical surrogate.
        printf("  Best hit:    %s (surrogate score = %.3f "
               "[model scale, NOT kcal/mol; claim_validity=proxy_only])\n",
               summary.top_hits[0].name.c_str(),
               summary.top_hits[0].dG_corrected);
    }
    printf("  Results:     %s\n", config.results_csv.c_str());
    printf("=================================\n\n");

    // Cleanup temp files
    library::cleanup_library(lig_lib);
    if (config.receptor_models > 1)
        library::cleanup_library(rec_lib);

    if (csv.is_open()) csv.close();

    return summary;
}

// ─── Write full results CSV ──────────────────────────────────────────────────

void write_results_csv(const std::string& path,
                       const std::vector<LigandResult>& results) {
    std::ofstream csv(path);
    // See run_campaign: `surrogate_*` naming is deliberate — these columns hold
    // composition-only placeholder scores, not binding free energies.
    csv << "rank,name,surrogate_consensus,surrogate_best,surrogate_corrected,"
        << "surrogate_mean,surrogate_stddev,"
        << "best_model,n_atoms,n_rotatable,n_rings,mw,"
        << "ref_entropy_corr_modelscale,time_sec,status\n";

    // Sort by dG_corrected
    std::vector<std::reference_wrapper<const LigandResult>> sorted(results.begin(), results.end());
    std::ranges::sort(sorted, {},
        [](const LigandResult& r) { return r.dG_corrected; });

    int rank = 0;
    for (const LigandResult& r : sorted) {
        rank++;
        csv << rank << ","
            << r.name << ","
            << r.dG_consensus << ","
            << r.dG_best << ","
            << r.dG_corrected << ","
            << r.dG_mean << ","
            << r.dG_stddev << ","
            << r.best_model + 1 << ","
            << r.n_atoms << ","
            << r.n_rotatable << ","
            << r.n_rings << ","
            << r.molecular_weight << ","
            << r.ref_entropy_correction << ","
            << r.dock_time_sec << ","
            << (r.success ? "OK" : r.error)
            << "\n";
    }
}

// ─── Write top hits as PDB files ─────────────────────────────────────────────

void write_top_hits(const std::string& output_dir,
                    const std::vector<LigandResult>& results,
                    int top_n) {
    fs::create_directories(output_dir);

    int n = std::min(top_n, static_cast<int>(results.size()));
    for (int i = 0; i < n; i++) {
        const auto& r = results[i];
        std::string filename = output_dir + "/rank_" + std::to_string(i + 1) +
                               "_" + r.name + ".pdb";
        std::ofstream pdb(filename);
        // REMARK labels carry the same firewall as the CSV: surrogate scores on
        // a model scale, explicitly stamped proxy_only, never kcal/mol.
        pdb << "REMARK  FlexAIDdS Campaign Result\n"
            << "REMARK  Rank: " << (i + 1) << "\n"
            << "REMARK  Name: " << r.name << "\n"
            << "REMARK  claim_validity: proxy_only (non-physical surrogate score)\n"
            << "REMARK  surrogate_corrected: " << r.dG_corrected << " [model scale]\n"
            << "REMARK  surrogate_consensus: " << r.dG_consensus << " [model scale]\n"
            << "REMARK  surrogate_best:      " << r.dG_best << " [model scale]\n"
            << "REMARK  Best model:          " << (r.best_model + 1) << "\n"
            << "REMARK  Ref entropy term:    " << r.ref_entropy_correction << " [model scale]\n"
            << "REMARK  Atoms: " << r.n_atoms
            << "  Rotatable: " << r.n_rotatable
            << "  Rings: " << r.n_rings
            << "  MW: " << r.molecular_weight << "\n";
        for (size_t ai = 0; ai < r.pose_xyz.size(); ++ai) {
            const auto& xyz = r.pose_xyz[ai];
            const std::string elem = element_symbol_from_z(
                ai < r.pose_atomic_numbers.size() ? r.pose_atomic_numbers[ai] : 0);
            pdb << std::left << std::setw(6) << "HETATM"
                << std::right << std::setw(5) << (ai + 1) << " "
                << std::setw(4) << elem << " "
                << "LIG A"
                << std::setw(4) << 1 << "    "
                << std::fixed << std::setprecision(3)
                << std::setw(8) << xyz[0]
                << std::setw(8) << xyz[1]
                << std::setw(8) << xyz[2]
                << std::setw(6) << "1.00"
                << std::setw(6) << "0.00"
                << "          "
                << std::setw(2) << elem
                << "\n";
        }
        pdb << "END\n";
    }
    printf("Top %d hits written to %s/\n", n, output_dir.c_str());
}

} // namespace campaign
