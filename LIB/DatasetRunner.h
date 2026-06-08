// =============================================================================
// DatasetRunner.h — Benchmark dataset runner for FlexAIDdS
//
// Downloads, prepares, and runs FlexAIDdS against standard docking benchmarks.
// Supports Astex Diverse, Astex Non-Native, HAP2, CASF-2016, PoseBusters,
// DUD-E, BindingDB-ITC, SAMPL6/7 host-guest, PDBbind Refined, and custom sets.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <mutex>
#include <atomic>
#include <vector>

#ifndef _MSC_VER
#include <sys/types.h>   // pid_t
#include <signal.h>      // kill, SIGTERM, SIGKILL
#endif

#include "TargetServer.h"

namespace dataset {

// =============================================================================
// Enums
// =============================================================================

/// Known benchmark dataset identifiers
enum class BenchmarkSet {
    ASTEX_DIVERSE,      // 85 complexes
    ASTEX_NON_NATIVE,   // 65 protein families (table), ~2200 cross-docking pairs
                        // (Verdonk 2008 original: 65 families, 1112 structures)
    HAP2,               // FlexAID JCIM 2015 HAP2 validation set
    CASF_2016,          // 285 complexes (PDBbind core)
    POSEBUSTERS,        // 308 complexes
    DUD_E,              // 102 targets + decoys
    BINDINGDB_ITC,      // ITC-validated subset
    SAMPL6_HG,          // 27 host-guest (OA/TEMOA/CB8)
    SAMPL7_HG,          // ~30 host-guest
    PDBBIND_REFINED,    // 5316 complexes
    CUSTOM_DOI,         // User-provided DOI → parse PDB codes
    CUSTOM_PDB_LIST     // User-provided PDB code list
};

/// Convert BenchmarkSet to string
inline std::string benchmark_set_name(BenchmarkSet s) {
    switch (s) {
        case BenchmarkSet::ASTEX_DIVERSE:    return "Astex Diverse";
        case BenchmarkSet::ASTEX_NON_NATIVE: return "Astex Non-Native";
        case BenchmarkSet::HAP2:             return "HAP2";
        case BenchmarkSet::CASF_2016:        return "CASF-2016";
        case BenchmarkSet::POSEBUSTERS:      return "PoseBusters";
        case BenchmarkSet::DUD_E:            return "DUD-E";
        case BenchmarkSet::BINDINGDB_ITC:    return "BindingDB-ITC";
        case BenchmarkSet::SAMPL6_HG:        return "SAMPL6 Host-Guest";
        case BenchmarkSet::SAMPL7_HG:        return "SAMPL7 Host-Guest";
        case BenchmarkSet::PDBBIND_REFINED:  return "PDBbind Refined";
        case BenchmarkSet::CUSTOM_DOI:       return "Custom DOI";
        case BenchmarkSet::CUSTOM_PDB_LIST:  return "Custom PDB List";
    }
    return "Unknown";
}

/// Parse string to BenchmarkSet
inline std::optional<BenchmarkSet> parse_benchmark_set(const std::string& name) {
    std::string lower = name;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower == "astex" || lower == "astex_diverse")   return BenchmarkSet::ASTEX_DIVERSE;
    if (lower == "astex_nonnative" || lower == "astex_non_native")
        return BenchmarkSet::ASTEX_NON_NATIVE;
    if (lower == "hap2")                                return BenchmarkSet::HAP2;
    if (lower == "casf2016" || lower == "casf_2016")    return BenchmarkSet::CASF_2016;
    if (lower == "posebusters")                         return BenchmarkSet::POSEBUSTERS;
    if (lower == "dude" || lower == "dud_e")            return BenchmarkSet::DUD_E;
    if (lower == "bindingdb_itc" || lower == "bindingdb") return BenchmarkSet::BINDINGDB_ITC;
    if (lower == "sampl6" || lower == "sampl6_hg")      return BenchmarkSet::SAMPL6_HG;
    if (lower == "sampl7" || lower == "sampl7_hg")      return BenchmarkSet::SAMPL7_HG;
    if (lower == "pdbbind" || lower == "pdbbind_refined") return BenchmarkSet::PDBBIND_REFINED;
    return std::nullopt;
}

// =============================================================================
// Data structures
// =============================================================================

/// A single entry in a benchmark dataset
struct DatasetEntry {
    std::string pdb_id;              // PDB code (uppercase)
    std::string receptor_path;       // path to downloaded PDB/CIF
    std::string ligand_path;         // path to extracted ligand SDF
    float experimental_affinity{-1.0f};  // pKd/pKi if available
    float experimental_dH{0.0f};     // ΔH in kcal/mol (ITC)
    float experimental_TdS{0.0f};    // TΔS in kcal/mol (ITC)
    std::string source;              // "Astex Diverse", "CASF-2016", etc.

    bool has_affinity()    const { return experimental_affinity >= 0.0f; }
    bool has_enthalpy()    const { return experimental_dH != 0.0f; }
    bool has_entropy()     const { return experimental_TdS != 0.0f; }
};

/// Result of docking a single entry
struct DockingResult {
    std::string pdb_id;
    float best_score{0.0f};           // FlexAIDdS free energy (kcal/mol)
    float rmsd_to_crystal{999.0f};    // serial-order RMSD to crystal ligand (Å)
    float rmsd_hungarian{999.0f};     // symmetry-corrected (Hungarian) RMSD (Å)
    float predicted_dG{0.0f};         // predicted ΔG (kcal/mol)
    float predicted_dH{0.0f};         // predicted ΔH (kcal/mol)
    float predicted_TdS{0.0f};        // predicted TΔS (kcal/mol)
    float predicted_IEE{0.0f};        // Enthalpy-Entropy Index (Williams 2017) — diagnostic only
    bool  has_IEE{false};             // false when |ΔG| < 1e-6 or not yet computed
    float shannon_entropy{0.0f};      // ensemble Shannon entropy
    int   num_poses{0};               // number of binding modes found
    double wall_time_s{0.0};          // docking wall time
    bool  success{false};             // RMSD < 2.0 Å
    // Clash diagnostics (populated from stdout parsing)
    long  individuals_clashed{0};     // total clashing evaluations
    long  individuals_total{0};       // total evaluations (across all generations)
    float clash_rate{0.0f};           // clashed / total — high (>0.95) = stuck GA
    bool  stuck{false};               // true when clash_rate > 0.95 and F > 0
};

/// Aggregate benchmark report
struct BenchmarkReport {
    std::string dataset_name;
    int total_systems{0};
    int successful{0};               // RMSD < 2.0 Å
    double success_rate{0.0};        // fraction
    double mean_rmsd{0.0};
    double median_rmsd{0.0};
    double pearson_r{0.0};           // predicted vs experimental affinity
    double spearman_rho{0.0};
    double kendall_tau{0.0};
    std::vector<DockingResult> results;

    // ── Cross-ligand competitive binding analysis (per receptor) ─────────
    /// One entry per unique receptor that had ≥2 ligands docked.
    struct CrossLigandResult {
        std::string receptor_id;                       // PDB code of the receptor
        int n_ligands{0};                              // total ligands docked against this receptor
        int n_completed{0};                            // ligands that produced binding modes
        std::vector<target::GrandPartitionFunction::LigandRank> ranked_ligands;  // sorted by ΔG
    };
    std::vector<CrossLigandResult> cross_ligand_results;
};

// =============================================================================
// Lightweight docking config for benchmarks
// =============================================================================

struct DockingConfig {
    // GA parameters — canonical benchmark spec (publication benchmarks):
    //   500 generations × 1000 chromosomes = 510,000 eval_chromosome calls/complex
    //   Matches config_defaults.h num_generations=500, num_chromosomes=1000.
    //   (BENCHMARKING_PLAN.md retired from this slim publication tree; spec lives in
    //    scripts/validate_benchmark_results.py thresholds + paper methods.)
    int    ga_generations{500};
    int    ga_population{1000};
    float  grid_spacing{0.375f};      // Å — 0.5 for coarse pass, 0.375 for full
    float  temperature{300.0f};       // Kelvin
    /// Concurrent FlexAIDdS worker processes (dataset-level parallelism).
    /// Each worker is an independent OS process; they do NOT share OMP threads.
    int    num_threads{1};
    /// OMP threads assigned to each FlexAIDdS subprocess.
    /// 0 = auto: floor(hardware_concurrency / num_threads), minimum 1.
    /// Explicit example: --threads 1 --omp-threads 6  (M3 Pro optimal)
    int    omp_threads_per_worker{0};
    bool   use_gpu{false};
    std::string gpu_backend{"cuda"};  // "cuda" or "metal"
    std::string output_dir{"."};
    std::string clustering_algorithm{"CF"}; // "CF", "FO" (FastOPTICS), or "DP" (DensityPeak)
    /// When true (default), skip targets whose output directory already contains
    /// at least one clustered pose PDB and a non-empty stdout.log.
    bool   skip_completed{true};
    /// Per-job timeout in seconds. 0 = no timeout (block indefinitely).
    /// Default 3600 s (1 h) — covers 8 min/complex with generous headroom.
    int    per_job_timeout_s{3600};
};

// =============================================================================
// SubprocessGuard — RAII process lifecycle manager
// =============================================================================

/// Tracks all spawned child PIDs and kills them on destruction.
/// Creates a dedicated process group so children can be killed as a batch.
/// Thread-safe: all mutations go through an internal mutex.
class SubprocessGuard {
public:
    SubprocessGuard();
    ~SubprocessGuard();

    /// Fork, set process group, exec.  Returns child PID (or -1 on failure).
    /// The child PID is registered for automatic cleanup.
    pid_t fork_exec(const std::string& cmd);

    /// Wait for a specific child with a timeout.
    /// Returns exit code (0-255) on normal exit, -1 on signal/timeout/error.
    /// On timeout, sends SIGTERM then SIGKILL after grace period.
    int wait_with_timeout(pid_t pid, int timeout_s);

    /// Unregister a PID (called after successful wait).
    void forget(pid_t pid);

    /// Kill ALL remaining registered children (SIGTERM, then SIGKILL).
    /// Called automatically by destructor, but safe to call manually.
    void kill_all();

    /// Number of currently tracked (still-running) children.
    size_t active_count() const;

private:
    mutable std::mutex mtx_;
    std::set<pid_t> pids_;
};

// =============================================================================
// Atom structure for ligand extraction
// =============================================================================

struct PDBAtom {
    int    serial{0};
    std::string name;
    std::string altLoc;
    std::string resName;
    std::string chainID;
    int    resSeq{0};
    float  x{0.0f}, y{0.0f}, z{0.0f};
    float  occupancy{1.0f};
    float  tempFactor{0.0f};
    std::string element;
    bool   is_hetatm{false};
};

// =============================================================================
// Statistical helper functions
// =============================================================================

/// Pearson correlation coefficient (computed from scratch)
double compute_pearson_r(const std::vector<double>& x, const std::vector<double>& y);

/// Spearman rank correlation (computed from scratch)
double compute_spearman_rho(const std::vector<double>& x, const std::vector<double>& y);

/// Kendall tau-b rank correlation (computed from scratch)
double compute_kendall_tau(const std::vector<double>& x, const std::vector<double>& y);

/// Compute RMSD between two coordinate sets (3N floats: x1,y1,z1,x2,y2,z2,...)
double compute_rmsd(const std::vector<float>& coords_a,
                    const std::vector<float>& coords_b);

// =============================================================================
// DatasetRunner class
// =============================================================================

class DatasetRunner {
public:
    /// Construct with cache directory (default: ~/.flexaidds/benchmarks/)
    explicit DatasetRunner(const std::string& cache_dir = "");

    /// Download and prepare a standard benchmark dataset.
    /// Returns list of ready-to-dock entries.
    std::vector<DatasetEntry> prepare(BenchmarkSet set);

    /// Download from a DOI: parse paper → extract PDB codes → fetch structures
    std::vector<DatasetEntry> prepare_from_doi(const std::string& doi);

    /// From a plain text file with one PDB code per line
    std::vector<DatasetEntry> prepare_from_pdb_list(const std::string& file_path);

    /// Run FlexAIDdS docking on all entries in the dataset.
    /// Returns per-system results + aggregate statistics.
    BenchmarkReport run(const std::vector<DatasetEntry>& entries,
                        const DockingConfig& config);

    /// Generate publication-ready report (markdown + CSV)
    void write_report(const BenchmarkReport& report,
                      const std::string& output_dir);

    /// Get the cache directory path
    const std::string& cache_dir() const { return cache_dir_; }

    // ── Public utilities for testing ──────────────────────────────────

    /// Download a PDB file from RCSB
    bool download_pdb(const std::string& pdb_id, const std::string& out_path);

    /// Download a CIF file from RCSB
    bool download_cif(const std::string& pdb_id, const std::string& out_path);

    /// Extract the largest non-water/non-ion HETATM ligand from a PDB/mmCIF file
    /// and write it as SDF
    bool extract_ligand(const std::string& structure_path, const std::string& out_sdf);

    /// Write a docking-ready receptor PDB with the cognate ligand removed.
    /// Self-docking requires the native binding site to be empty: the source
    /// PDB still contains the crystal ligand as HETATM, so docking it back
    /// overlaps the embedded copy (r->0 in the r^-12 wall term -> the native
    /// pose is rejected as a clash and the GA is forced into decoy pockets).
    /// Any receptor HETATM within `tol` Å of a ligand SDF atom is dropped.
    /// Returns the cleaned receptor path, or the original on failure.
    std::string write_receptor_without_ligand(const std::string& receptor_path,
                                               const std::string& ligand_sdf,
                                               const std::string& out_receptor,
                                               float tol = 1.3f);

    /// Parse PDB HETATM records into atom structures
    std::vector<PDBAtom> parse_pdb_hetatm(const std::string& pdb_path);

    /// Get the Astex Diverse 85 PDB codes
    static std::vector<std::string> astex_diverse_codes();

    /// Get the CASF-2016 PDB codes (285)
    static std::vector<std::string> casf2016_codes();

    /// Get the DUD-E target list (102)
    static std::vector<std::string> dude_targets();

    /// Get HAP2 target info (59 hardcoded PDB codes)
    static std::vector<std::string> hap2_codes();

private:
    std::string cache_dir_;

    // ── Subprocess lifecycle ─────────────────────────────────────────
    /// RAII guard that tracks all spawned children and kills on destruction.
    /// Shared across all worker threads in run().
    std::unique_ptr<SubprocessGuard> proc_guard_;

    /// Atomic shutdown flag — set by SIGINT/SIGTERM handler or manual stop.
    /// Worker threads check this before launching each new job.
    std::atomic<bool> shutdown_requested_{false};

    // ── TargetServer integration ─────────────────────────────────────
    // One TargetServer per unique receptor. Keyed by receptor_path.
    // Populated in run() before the thread pool starts.
    std::map<std::string, std::unique_ptr<target::TargetServer>> target_servers_;
    std::mutex target_server_mtx_;  // protects target_servers_ map creation

    // ── Dataset-specific fetchers ────────────────────────────────────

    std::vector<DatasetEntry> fetch_astex();
    std::vector<DatasetEntry> fetch_astex_nonnative();
    std::vector<DatasetEntry> fetch_hap2();
    std::vector<DatasetEntry> fetch_casf2016();
    std::vector<DatasetEntry> fetch_posebusters();
    std::vector<DatasetEntry> fetch_bindingdb_itc();
    std::vector<DatasetEntry> fetch_sampl6();
    std::vector<DatasetEntry> fetch_sampl7();
    std::vector<DatasetEntry> fetch_pdbbind_refined();
    std::vector<DatasetEntry> fetch_dud_e();

    // ── PDB structure fetching and preparation ───────────────────────

    /// Generic HTTP download using system curl
    bool http_download(const std::string& url, const std::string& out_path);

    /// Execute a system command and return its exit code
    int exec_cmd(const std::string& cmd);

    /// Execute a system command and capture stdout
    std::string exec_cmd_output(const std::string& cmd);

    /// Execute a docking command with PID tracking and timeout.
    /// Uses SubprocessGuard for process groups and orphan cleanup.
    int exec_dock(const std::string& cmd, int timeout_s);

    /// Ensure a directory exists (create recursively if needed)
    bool ensure_dir(const std::string& path);

    /// DOI parsing: fetch DOI metadata, extract PDB codes from abstract/text
    std::vector<std::string> extract_pdb_codes_from_doi(const std::string& doi);

    /// Expand ~ in paths
    std::string expand_home(const std::string& path);

    /// Download the preferred RCSB coordinate format for a structure.
    /// mmCIF is attempted first because legacy PDB files are lossy or absent
    /// for some benchmark codes; PDB remains only a last-resort fallback.
    bool download_structure(const std::string& pdb_id,
                            const std::string& entry_dir,
                            std::string& out_path);

    /// Prepare a single RCSB entry: download structure + extract ligand
    DatasetEntry prepare_pdb_entry(const std::string& pdb_id,
                                   const std::string& dataset_name,
                                   float affinity = -1.0f,
                                   float dH = 0.0f, float dS = 0.0f);

    /// Common ligand residue exclusion set (water, ions, buffers)
    static const std::set<std::string>& excluded_residues();
};

// =============================================================================
// Astex Non-Native cross-docking mapping
// =============================================================================

/// Astex Non-Native: target PDB → list of alternative conformers for cross-docking
/// Based on Verdonk et al. (2008) J. Chem. Inf. Model.
struct AstexNonNativeTarget {
    std::string target_name;
    std::string native_pdb;
    std::vector<std::string> alternative_pdbs;
};

/// Get the full Astex Non-Native target list (65 targets, 1112 structures)
std::vector<AstexNonNativeTarget> astex_nonnative_targets();

} // namespace dataset
