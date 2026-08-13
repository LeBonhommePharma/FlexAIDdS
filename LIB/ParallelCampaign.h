// ParallelCampaign.h — GPU/SIMD/OpenMP-accelerated virtual screening campaigns
//
// Orchestrates massively parallel docking across:
//   Level 1: Ligand library (N ligands) — OpenMP task parallelism
//   Level 2: Receptor conformers (M models) — OpenMP nested parallelism
//   Level 3: GA population evaluation — CUDA/ROCm/Metal/AVX batch kernels
//
// For a 10,000-ligand × 20-model campaign (200,000 docks):
//   - Level 1: OpenMP distributes ligands across CPU cores
//   - Level 2: Each ligand docks against all models (inner parallel loop)
//   - Level 3: Each GA generation evaluates population on GPU
//   - Ensemble consensus + reference entropy computed per ligand
//   - Results streamed to disk as they complete (no memory blowup)
//
// Hardware utilization:
//   CUDA/ROCm/Metal: GA fitness evaluation (existing cuda_eval/hip_eval/metal_eval)
//   AVX-512/AVX2:    Boltzmann weights, log-sum-exp, distance geometry
//   Eigen:           Coordinate transforms, matrix operations
//   OpenMP:          Ligand-level and model-level parallelism
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "LibrarySplitter.h"
#include "ReferenceEntropy.h"
#include "flexaid.h"

#include <string>
#include <vector>
#include <functional>
#include <atomic>
#include <chrono>
#include <array>

namespace campaign {

// ─── Campaign configuration ─────────────────────────────────────────────────

struct CampaignConfig {
    // Parallelism
    int max_ligand_threads    = 0;   // 0 = auto (OMP_NUM_THREADS or nproc)
    int max_model_threads     = 1;   // models per ligand (usually 1; >1 if many cores)
    bool gpu_per_ligand       = false; // true: each ligand thread gets its own GPU stream
    
    // Receptor
    std::string receptor_path;
    int receptor_models       = 1;   // auto-detected from file
    
    // Ligand library
    std::string ligand_path;
    int ligand_count          = 0;   // auto-detected
    
    // Docking parameters (copied from FA/GB globals)
    double temperature_K      = 300.0;
    bool use_rigid            = false;
    bool use_folded           = false;
    std::string config_json;
    std::string output_prefix = "campaign";
    
    // Reference entropy
    bool compute_ref_entropy  = true;  // auto-enabled for multi-model
    
    // Output
    bool stream_results       = true;  // write results as they complete
    std::string results_csv;           // path to results CSV (auto-generated)

    // P3 grand canonical: default concentration forwarded for TargetServer/GPF when multi-ligand campaign uses grand paths
    double default_conc_M     = 1.0;

    // Optional Stage-1 cube prefilter (default OFF). When true, the ligand
    // library is coarse-screened and only top_n names are docked. Requires a
    // MOL2 target + cleft/site PDB. Does not change scoring of kept ligands.
    bool coarse_prefilter = false;
    int  coarse_prefilter_top_n = 100;
    std::string coarse_target_mol2;
    std::string coarse_cleft_pdb;
};

// ─── Per-ligand result ──────────────────────────────────────────────────────

struct LigandResult {
    std::string name;
    std::string smiles;           // if available
    int    ligand_idx;
    
    // Per-model dG values
    std::vector<double> per_model_dG;
    
    // Ensemble consensus
    double dG_consensus;          // Boltzmann-weighted
    double dG_best;               // best single model
    double dG_mean;
    double dG_stddev;
    int    best_model;
    
    // Reference-corrected
    double dG_corrected;          // dG_consensus + ref entropy correction
    double ref_entropy_correction;
    
    // Diagnostics
    int    n_atoms;
    int    n_rotatable;
    int    n_rings;
    float  molecular_weight;
    double dock_time_sec;
    bool   success;
    std::string error;

    // Optional ligand pose coordinates for top-hit export.
    // These are the coordinates available from preprocessing for now and
    // allow writing coordinate-bearing PDB files.
    std::vector<std::array<float, 3>> pose_xyz;
    std::vector<int> pose_atomic_numbers;
};

// ─── Campaign summary ───────────────────────────────────────────────────────

struct CampaignSummary {
    int total_ligands;
    int successful;
    int failed;
    int receptor_models;
    int total_docks;              // ligands × models
    
    double total_time_sec;
    double avg_time_per_ligand;
    double throughput_per_hour;   // ligands/hour
    
    // Top hits (sorted by dG_corrected)
    std::vector<LigandResult> top_hits;  // top 100 or all if < 100
    
    // Hardware utilization
    std::string gpu_backend;      // "CUDA", "ROCm", "Metal", "CPU"
    int cpu_threads_used;
};

// ─── Progress callback ──────────────────────────────────────────────────────

using ProgressCallback = std::function<void(
    int completed,        // ligands completed so far
    int total,            // total ligands
    const LigandResult& latest  // most recent result
)>;

// ─── Main entry point ───────────────────────────────────────────────────────

/// Run a parallel virtual screening campaign.
/// This is the top-level function that replaces the single-ligand docking
/// path when a library is detected.
CampaignSummary run_campaign(
    const CampaignConfig& config,
    ProgressCallback progress = nullptr);

/// Auto-configure a campaign from command-line arguments.
/// Detects library size, receptor model count, GPU availability,
/// and sets optimal thread counts.
CampaignConfig auto_configure(
    const std::string& receptor_path,
    const std::string& ligand_path,
    const std::string& config_json = "",
    const std::string& output_prefix = "campaign",
    bool rigid = false,
    bool folded = false);

/// Write results CSV with all ligand scores.
void write_results_csv(const std::string& path,
                       const std::vector<LigandResult>& results);

/// Write top-N hits as individual PDB files.
void write_top_hits(const std::string& output_dir,
                    const std::vector<LigandResult>& results,
                    int top_n = 100);

} // namespace campaign
