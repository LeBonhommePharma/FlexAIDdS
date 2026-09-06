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
#include <limits>
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
#include "ProtocolConfig.h"
#include "DatasetRunnerStats.h"  // pure correlation / RMSD helpers (P0 leaf)
#include "DatasetRunnerProvenance.h"  // provenance.json writer (P1 leaf)

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
    std::string rmsd_reference_path; // optional crystal/reference SDF for RMSD/native CF
    std::string binding_site_path;   // oracle binding site PDB (optional; enables LOCCLF mode)
    std::string cleft_sphere_path;   // explicit Get_Cleft sphere PDB for restored multi-cleft runs
    float experimental_affinity{-1.0f};  // pKd/pKi if available
    float experimental_dH{0.0f};     // ΔH in kcal/mol (ITC)
    float experimental_TdS{0.0f};    // TΔS in kcal/mol (ITC)
    std::string source;              // "Astex Diverse", "CASF-2016", etc.

    bool has_affinity()    const { return experimental_affinity >= 0.0f; }
    bool has_enthalpy()    const { return experimental_dH != 0.0f; }
    bool has_entropy()     const { return experimental_TdS != 0.0f; }
    bool has_oracle_site() const { return !binding_site_path.empty(); }
    bool has_cleft_spheres() const { return !cleft_sphere_path.empty(); }
    double conc_M = 1.0;  // P3: ligand concentration for grand canonical
};

/// Result of docking a single entry.
///
/// CSV column names (`best_score`, `predicted_dG`, …) are stable for live
/// campaign compatibility. Semantics (do not conflate):
///   - best_score  ≡ CF/contact-function scoring proxy of the elected pose
///                   (alias concept: cf_score / elected REMARK CF). NOT ΔG.
///   - predicted_dG ≡ ensemble-derived free-energy *estimate* F when the
///                   Post-GA StatMech ledger is available; otherwise may
///                   fall back to CF. NOT experimental binding free energy
///                   unless full Z + vib + solvent/concentration path is
///                   active, validated, and claimed under AGENTS.md rules.
struct DockingResult {
    std::string pdb_id;
    // CF/contact-function scoring proxy of the elected (rank-0) pose.
    // CSV column name kept as `best_score` for campaign compatibility.
    // Prefer language "CF score" / "cf_score" in new docs; see elected_cf.
    float best_score{0.0f};
    float rmsd_to_crystal{-1.0f};     // serial-order RMSD to crystal ligand (Å); -1 = not computed/failed
    float rmsd_hungarian{-1.0f};      // symmetry-corrected (Hungarian) RMSD (Å); -1 = not computed/failed
    // Why rmsd_to_crystal is -1, when it is. "none" iff RMSD >= 0.
    // Values: none | ref_empty | pose_block_empty | count_mismatch |
    //         elem_mismatch | elem_order_mismatch | input_missing.
    // Exists because a bare -1 is ambiguous between a per-pose failure and a
    // wholesale reference-resolution failure (bug 2026-08-22: campaign arms
    // 8/9/10 wrote valid poses with all-RMSD=-1 and 0% success summaries).
    std::string rmsd_fail_reason{"none"};
    // Ensemble free-energy estimate (F = -kT ln Z) when available; else CF
    // fallback. CSV column `predicted_dG` is a historical name — not exp. ΔG.
    float predicted_dG{0.0f};
    // Configurational ΔH ≈ <E> from the same ledger as predicted_dG (kcal/mol proxy units)
    float predicted_dH{0.0f};
    // Configurational TΔS estimate from the ledger (kcal/mol proxy units)
    float predicted_TdS{0.0f};
    float predicted_IEE{0.0f};        // Enthalpy-Entropy Index (Williams 2017) — diagnostic only
    bool  has_IEE{false};             // false when |predicted_dG| < 1e-6 or not yet computed
    float shannon_entropy{0.0f};      // conformational Shannon entropy -Σ p_i ln p_i (nats)
    float search_entropy_proxy{0.0f}; // legacy H_final collapse proxy from GA energy histogram (nats)
    int   num_poses{0};               // number of binding modes found
    // Runtime completion is independent of RMSD/PB success. -1 includes not
    // started or a timeout; the retained child log distinguishes those causes.
    int   docking_exit_code{-1};
    bool  docking_completed{false};
    std::string matrix_md5;          // actual selected matrix input, not an expected default
    double wall_time_s{0.0};          // docking wall time
    // ── Success gates (fixed semantics; never remapped by env) ────────────
    // success_rmsd : ordered direct rmsd_to_crystal <= 2 Å && !seed_echo
    //                (rmsd_hungarian is diagnostic only; never sets success)
    // pb_pass      : PoseBust on elected BindingMode pose (bust_cli preferred;
    //                native_pose_qc / native_pose_qc_fallback when Off/missing).
    //                Never true unless pb_ran (validate_elected_pose contract).
    // success_pb   : success_rmsd && pb_pass
    // claim_ready  : success_pb && pb_backend==bust_cli && tENCoM/Eigen + hashes
    // success      : always == success_rmsd (legacy column)
    bool  success_rmsd{false};
    bool  pb_pass{false};
    bool  success_pb{false};          // == success_rmsd && pb_pass
    bool  claim_ready{false};
    bool  success{false};             // == success_rmsd (stable legacy meaning)
    // Protocol-level native-pose exposure. This remains true even when the
    // elected file is a GA cluster descendant rather than the literal _INI.pdb.
    // A seeded oracle-retention result is useful diagnostically, but cannot be
    // presented as no-seed redocking success.
    bool  native_pose_seeded{false};
    float native_pose_seed_fraction{0.0f};
    bool  protocol_claim_eligible{false};
    // PoseBusters summary (backend-selected: native_pose_qc or bust_cli)
    bool  pb_ran{false};
    int   pb_n_pass{0};
    int   pb_n_fail{0};
    int   pb_n_checks{0};
    std::string pb_failed_keys;
    std::string pb_backend;  // "native_pose_qc" | "bust_cli" | "skipped" | "error"
    // NativePoseQC full-suite report (always filled as a parity diagnostic)
    bool  native_qc_ran{false};
    bool  native_qc_pass{false};
    std::string native_qc_failed_keys;
    float pb_min_lig_prot_dist{std::numeric_limits<float>::quiet_NaN()};
    float pb_volume_overlap{std::numeric_limits<float>::quiet_NaN()};
    // Validator / provenance (must cite same pose hash)
    std::string elected_pose_path;    // absolute or workdir-relative path to elected_pose.pdb
    std::string elected_pose_source;  // original restart path e.g. r1/1G9V_2.pdb
    int         elected_restart{-1};  // -1 = r0 / root
    int         elected_cluster{-1};  // pose index 0–19
    float       elected_cf{std::numeric_limits<float>::quiet_NaN()};
    bool        score_pose_consistent{false}; // emitted coordinates reproduce elected CF
    float       score_pose_delta{std::numeric_limits<float>::quiet_NaN()};
    std::string pose_sha256;          // SHA-256 of elected_pose.pdb
    std::string rmsd_pose_sha256;     // exact elected-pose bytes consumed by RMSD
    std::string posebusters_pose_sha256;  // parent elected PDB consumed by PoseBusters
    std::string posebusters_input_sha256; // derived predicted-ligand SDF consumed by bust
    std::string tencom_status{"not_run"};  // ok | fail | not_run | skipped
    std::string eigen_status{"not_run"};
    std::string tencom_pose_sha256;   // re-hash of the exact pose consumed by tENCoM/Eigen
    int         eigen_n_modes{0};     // positive ligand ANM eigenmodes on elected pose
    float       elected_H_vib{0.0f}; // H(omega) of the exact elected pose (nats)
    std::string eigen_model{"ligand_cartesian_anm"};
    // Clash diagnostics (populated from stdout parsing)
    long  individuals_clashed{0};     // total clashing evaluations
    long  individuals_total{0};       // total evaluations (across all generations)
    float clash_rate{0.0f};           // clashed / total — high (>0.95) = stuck GA
    bool  stuck{false};               // true when clash_rate > 0.95 and F > 0
    // Native-pose CF diagnostic (scored before the GA via FLEXAIDDS_SCORE_NATIVE)
    float cf_native{0.0f};            // CF at crystal pose; 0.0 when not run
    // Conditional scanned-pool ceiling: min ordered direct RMSD over heads/members
    // actually enumerated by the ceiling scan (not guaranteed any-pose / not full
    // emission census unless every head+member file is present). Never mutates
    // elected pose, seed_echo, or pose_source. CSV alias: best_cluster_rmsd.
    float conditional_scanned_pool_ceiling{-1.0f};
    float best_cluster_rmsd{-1.0f};   // == conditional_scanned_pool_ceiling (legacy column)
    int   best_cluster_idx{-1};
    float cf_best_cluster{std::numeric_limits<float>::quiet_NaN()};
    // seed_echo: true when elected path ends in "_INI.pdb". Immutable after set;
    // pool ceiling must never clear it.
    bool  seed_echo{false};
    // Pose provenance: "ini_elitism" | "ga_cluster" | "cf_rank0" | "softbeta" | ""
    std::string pose_source{""};
    // ── Pose-election outcome ────────────────────────────────────────────
    // Which rule actually elected the reported pose, and what the cross-restart
    // consensus vote was.  These were previously emitted only to std::cerr
    // (the "[CONSENSUS]" line), so a successful run -- whose stderr CI discards
    // -- kept the elected RMSD but discarded the rule that produced it.  Two
    // runs with identical configs can elect different poses; without these the
    // artifact cannot say why.  Recorded beside rmsd_to_crystal for that reason.
    //
    // Note the election rule is not independent of FLEXAIDDS_RESTARTS: the
    // consensus veto needs cluster_consensus_k (default 3) distinct restarts to
    // vote, so at restarts < k it can never fire.  Comparing runs across restart
    // counts therefore varies the rule as well as the search budget, which is
    // exactly what these fields make visible.
    // These describe the pose actually REPORTED, so they are set where the
    // election applies its result -- after the v124 guard resolves, not before.
    // Recording them earlier names the incumbent pose the election replaced.
    //
    // "" when no election ran (single candidate / election block skipped).
    // "guard-protected" when the v124 guard VETOED the override and kept the INI
    // seed: naming the gate mode there would credit a rule that did not decide
    // this pose.
    std::string election_mode{""};   // "consensus" | "entropy-midwall" | "entropy-contact"
                                     // | "guard-protected" | ""
    // consensus_count == -1 is AMBIGUOUS ON ITS OWN and must be read together
    // with election_mode:
    //     election_mode ""                -> no election ran (pool < 2 candidates)
    //     election_mode "guard-protected" -> election ran and was VETOED; the
    //                                        elected INI seed is not in the pool,
    //                                        so it has no vote count
    // Reading this column alone cannot distinguish "nothing happened" from
    // "an override was overruled".
    int  consensus_count{-1};        // cross-restart votes for the ELECTED pose
    bool rank0_demoted{false};       // true when the elected pose is not the min-CF one
    // Separate estimands: generator CF top-1 vs entropy/consensus reranked top-1
    std::string cf_top1_pose_path;
    std::string cf_top1_pose_sha256;
    float       cf_top1_score{std::numeric_limits<float>::quiet_NaN()};
    float       cf_top1_rmsd{-1.0f};
    std::string entropy_top1_pose_path;
    std::string entropy_top1_pose_sha256;
    float       entropy_top1_score{std::numeric_limits<float>::quiet_NaN()};
    float       entropy_top1_rmsd{-1.0f};
    // Level-3 H(ω) vibrational-entropy diagnostic. Enabled by default; set
    // FLEXAIDDS_HVIB=0 only for non-claim diagnostic runs. Shannon entropy over ligand ANM eigenvalue spectra of the
    // top-10 emitted cluster reps. See compute_target_hvib() in DatasetRunner.cpp.
    float H_rep_rank0{0.0f};  // H(ω) of rank-0 (best-CF) cluster rep, computed individually
    float H_pop{0.0f};        // pooled population vibrational entropy H(ω)
    float H_rep_mean{0.0f};   // mean per-rep vibrational entropy
    float D_vib{0.0f};        // inter-rep vibrational divergence
    // ThermodynamicEngine decomposition (populated when thermo_engine enabled in config)
    float thermo_G_bind{0.0f};
    float thermo_H_vct{0.0f};        // per-heavy-atom (intensive)
    float thermo_H_vct_raw{0.0f};    // unnormalized ensemble mean CF
    int   thermo_n_heavy{0};         // heavy-atom count used for normalization
    float thermo_TdS_shannon{0.0f};
    float thermo_TdS_vib{0.0f};
    float thermo_D_vib{0.0f};        // H_rep_bound (raw bound-complex tENCoM entropy)
    float thermo_compensation{0.0f};
    bool  has_thermo{false};
    // Reporting-only whiteboard diagnostics (parsed from [THERMO2]); computed
    // at thermo_report_T (default 21.0 = kT_ISMB, ISMB 2017), independent of
    // thermo_T_eff above — never affects thermo_G_bind/CF scoring/GA
    // selection. Whiteboard convention: T=21 defines the LEFT-hand quantity
    // (report_T is that constant, echoed for "(T=21)" labelling downstream).
    float thermo_report_T{21.0f};
    float thermo_I_ES{0.0f};
    float thermo_CF_r2s{0.0f};
    std::string thermo_binding_regime{};
    bool  has_thermo2{false};
// P1: real ensemble log_Z from BindingPopulation.get_log_Z() (preferred over -dG/kT)
    double ensemble_log_Z{0.0};
    // Mid-run H_shannon snapshots at fixed generations (causality test).
    // Populated when FLEXAIDDS_THERMO=1. NaN when that generation was not
    // reached (early exit) or when thermo is disabled.
    float thermo_TdS_shannon_gen500{std::numeric_limits<float>::quiet_NaN()};
    float thermo_TdS_shannon_gen1000{std::numeric_limits<float>::quiet_NaN()};
};

/// Aggregate benchmark report
struct BenchmarkReport {
    std::string dataset_name;
    int total_systems{0};
    int successful{0};               // count of success (== success_rmsd)
    double success_rate{0.0};
    int successful_rmsd{0};
    int successful_pb{0};            // success_rmsd && pb_pass
    int claim_ready_count{0};
    double success_rate_rmsd{0.0};
    double success_rate_pb{0.0};
    double claim_ready_rate{0.0};
    int completed_systems{0};
    int valid_rmsd_count{0};
    double mean_rmsd{std::numeric_limits<double>::quiet_NaN()};
    double median_rmsd{std::numeric_limits<double>::quiet_NaN()};
    // Zero-success plausibility gate (DatasetRunnerStats.h): true when the
    // summary would certify 0% while poses exist and negative RMSDs are
    // dominated by wholesale measurement-side reasons (bug 2026-08-22,
    // arms 8/9/10). Emitted as the trailing `suspect_zero_success` column of
    // <dataset>_summary.csv; never alters any success count.
    bool suspect_zero_success{false};
    int affinity_pairs{0};           // entries with both exp affinity and predicted_dG (F-est / CF fallback)
    double pearson_r{std::numeric_limits<double>::quiet_NaN()};   // predicted_dG-derived pKd vs experimental affinity
    double spearman_rho{std::numeric_limits<double>::quiet_NaN()};
    double kendall_tau{std::numeric_limits<double>::quiet_NaN()};
    std::vector<DockingResult> results;

    // ── Cross-ligand competitive binding analysis (per receptor) ─────────
    /// One entry per unique receptor that had ≥2 ligands docked.
    struct CrossLigandResult {
        std::string receptor_id;                       // PDB code of the receptor
        int n_ligands{0};                              // total ligands docked against this receptor
        int n_completed{0};                            // ligands that produced binding modes
        // Grand-PF ranking by ensemble F estimate (not experimental ΔG_bind)
        std::vector<target::GrandPartitionFunction::LigandRank> ranked_ligands;
    };
    std::vector<CrossLigandResult> cross_ligand_results;
};

// =============================================================================
// Runtime exit status for the CLI; a completed inaccurate pose is not a
// process failure. Preparation/listing intentionally has no docking results.
inline int benchmark_runtime_exit_code(const BenchmarkReport& report,
                                       bool no_docking = false) {
    if (no_docking) return 0;
    if (report.total_systems <= 0 ||
        report.results.size() != static_cast<size_t>(report.total_systems)) return 2;
    for (const auto& result : report.results) {
        if (!result.docking_completed || result.docking_exit_code != 0 ||
            result.num_poses <= 0 || result.stuck) return 2;
    }
    return 0;
}

// Lightweight docking config for benchmarks
// =============================================================================

/// Layer 1: Explicit benchmark protocol mode.
/// Controls both seed_elitism and pose-blinding behavior in DatasetRunner::run().
/// UNSET                 → legacy env-var behavior (FLEXAIDDS_SEED_ELITISM, backward-compat).
/// ORACLE_CEILING        → seed_elitism ON, blinding OFF; ceiling measurement with crystal IC.
/// DEFINED_CLEFT_REDOCK  → seed_elitism OFF, blinding ON; known cleft injected.
/// AUTONOMOUS            → seed_elitism OFF, blinding ON; thesis/publication number.
enum class BenchmarkMode {
    UNSET,            ///< legacy env-var behavior (backward compatible)
    ORACLE_CEILING,   ///< seed_elitism=ON, blinding=OFF (crystal IC anchor)
    DEFINED_CLEFT_REDOCK, ///< known cleft/site, no crystal pose seed
    AUTONOMOUS,       ///< seed_elitism=OFF, blinding=ON (thesis number)
};

struct DockingConfig {
    // GA parameters — canonical benchmark spec (publication benchmarks):
    //   P6: base 2000 generations × 1000 chromosomes, scaled per-target by
    //   ceil(n_genes/4) in DatasetRunner (so a rigid 4-gene ligand gets 2000 gens
    //   and a 12-gene flexible ligand gets 6000) → search budget tracks DoF.
    //   Was 500 (v22) — quadrupled to give the false-minimum-prone targets enough
    //   generations to escape the deepest VCT well via the diversity machinery.
    //   (BENCHMARKING_PLAN.md retired from this slim publication tree; spec lives in
    //    scripts/validate_benchmark_results.py thresholds + paper methods.)
    int    ga_generations{2000};
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
    /// Ablation switch. When true, every generated dock config is written with
    /// flexibility.intramolecular=false, forcing legacy rigid-body docking
    /// (4 genes: translation + 3 rotations, zero torsional DoF). Default false:
    /// ligands dock flexibly (one dihedral gene per perceived rotatable bond).
    bool   force_rigid{false};
    /// Binding-site rotamer pre-relaxation (Option 3 — apo-strain fix).
    /// Greedy Dunbrack rotamer search on pocket sidechains before docking.
    /// Default true (v44+): enabled by default for all benchmark runs.
    /// Set false only for ablation / legacy-baseline comparison.
    bool   receptor_rotamer_prep{true};   // Option 3 apo-strain fix (v44 default)
    /// Layer 1 benchmark protocol selector.
    /// ORACLE_CEILING:       seed_elitism=ON, blinding=OFF (ceiling measurement).
    /// DEFINED_CLEFT_REDOCK: seed_elitism=OFF, blinding=ON, known cleft injected.
    /// AUTONOMOUS:           seed_elitism=OFF, blinding=ON (publication/thesis number).
    /// UNSET (default): preserves legacy env-var-based behavior.
    BenchmarkMode mode{BenchmarkMode::UNSET};
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

    /// Fork, set process group, exec via `/bin/sh -c`.  Returns child PID
    /// (or -1 on failure). The child PID is registered for automatic cleanup.
    /// Prefer fork_exec_argv when no shell metacharacters/redirection are needed.
    pid_t fork_exec(const std::string& cmd);

    /// Fork + execvp(argv) without a shell. argv[0] is the program path;
    /// remaining entries are arguments (no env KEY=VAL prefixes). Use this for
    /// pure binary invocations; use fork_exec(shell string) when the command
    /// must set env vars or redirect stdout/stderr.
    pid_t fork_exec_argv(const std::vector<std::string>& argv);

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
// Restart launch throttling — SCHEDULING ONLY
//
// These helpers decide *when* restart child processes are forked. They never
// touch a child's config, seed, output prefix, or OMP_NUM_THREADS, so they
// cannot move a docked coordinate or a score.
// =============================================================================

/// One restart's exit code, recorded in launch order.
struct RestartResult {
    int ri{0};    ///< restart index (0 = the canonical restart)
    int ret{-1};  ///< exit code; -1 = signal / timeout / fork failure
};

/// Fold per-restart exit codes into a single target-level return code.
///
/// Preserves the historical propagation rule EXACTLY: restart 0's code sets the
/// base, and ANY non-zero code from restarts 1..n-1 overrides it (so the last
/// non-zero failure wins). One failed restart therefore still poisons the whole
/// target, which is what makes `docking_completed` false downstream.
/// `fallback` is returned unchanged when `results` is empty.
int fold_restart_return_codes(const std::vector<RestartResult>& results,
                              int fallback);

/// How many restart children may be alive concurrently.
///
/// `omp_per_worker` already divides the CPU budget by the number of concurrent
/// workers, but the parallel-restart path forks every restart of a target at
/// once and each child inherits `OMP_NUM_THREADS=omp_per_worker`. Real demand
/// was therefore `workers * n_restarts * omp_per_worker` — the restart fan-out
/// was never in the denominator. This bounds it so that
/// `cap * omp_per_worker * num_workers ~= cpu_budget`.
///
/// `env_override` (FLEXAIDDS_MAX_CONCURRENT_RESTARTS):
///   < 0 → auto-derive from the CPU budget
///   = 0 → unlimited (legacy pre-fix fan-out)
///   > 0 → explicit cap
/// Returns 0 for "unlimited"; otherwise >= 1.
int restart_concurrency_cap(int cpu_budget, int omp_per_worker,
                            int num_workers, int env_override);

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

// Statistical helpers (compute_pearson_r / spearman / kendall / compute_rmsd)
// live in DatasetRunnerStats.h — included above for a stable public API.

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

    /// Prepare a single RCSB entry for cognate redocking:
    /// download structure → extract cognate ligand SDF → apo receptor (ligand stripped).
    /// Used by `FlexAIDdS --redock <PDBid>` and by benchmark fetchers.
    /// Cache layout: <cache_dir>/<dataset_name>/<PDBID>/{PDBID.pdb, PDBID_ligand.sdf, PDBID_apo.pdb}
    DatasetEntry prepare_pdb_entry(const std::string& pdb_id,
                                   const std::string& dataset_name,
                                   float affinity = -1.0f,
                                   float dH = 0.0f, float dS = 0.0f);

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

    /// Typed protocol snapshot (seed/restarts/VCT/GA/thermo/data_dir + chunk-2
    /// pose/budget/site knobs). Loaded via ProtocolConfig::from_env() in the
    /// constructor and re-snapshotted at run() entry so long-lived runners pick
    /// up env changes (see docs/implementation/protocol-config.md).
    flexaids::ProtocolConfig protocol_cfg_{};

    // ── Pose selector tuning ─────────────────────────────────────────
    /// CF-window gate (Fix A): when true the freq>1 gate in
    /// select_pose_freq_gated_pooled() also admits singleton clusters whose
    /// CF is within 30 units of the pool minimum. Read once from
    /// FLEXAIDDS_CF_WINDOW_SELECTOR in the constructor (default off).
    bool cf_window_selector_ = false;

    /// Cluster member emission (Fix B): when true, the oracle BCR scan also
    /// folds in cluster *member* poses — not just the cluster representative
    /// (`_N.pdb`) — for "near-miss" clusters whose representative Hungarian
    /// RMSD falls in [2.0, 4.0] Å. A near-native sub-Å member can hide inside
    /// a clash-attractor cluster whose centroid is displaced far from native
    /// (e.g. 1OF6: cluster-0 centroid 14.2 Å, member 2.30 Å); emitting and
    /// scoring members recovers that pose for BCR. Read once from
    /// FLEXAIDDS_CLUSTER_MEMBER_EMIT in the constructor (default off).
    ///
    /// NOTE on coordinate availability: DatasetRunner drives FlexAIDdS as a
    /// subprocess and has no in-process FA atoms[] array. Member Cartesian
    /// coordinates are therefore NOT held in any cluster struct here — the
    /// engine writes only the cluster representative PDB plus a `.mcf` sidecar
    /// of member CF values. Reconstructing member Cartesians from chromosome
    /// IC coordinates is intentionally out of scope (too fragile). This path
    /// consequently recovers members only when member pose PDBs following the
    /// `<head_stem>_member<M>.pdb` convention are present on disk (head_stem is
    /// the CF/DP or FO dual-suffix head without `.pdb`).
    bool cluster_member_emit_ = false;

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
    /// Prefers PDB for fixed-column coordinate parsing; CIF fallback when PDB
    /// is unavailable (large structures / withdrawn formats).
    bool download_structure(const std::string& pdb_id,
                            const std::string& entry_dir,
                            std::string& out_path);

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

// =============================================================================
// Emitting cluster-head enumeration (CF/DP single-suffix + FO dual-suffix)
// =============================================================================

/// One cluster-head PDB under a restart out_prefix.
/// CF/DP: <prefix>_<rank>.pdb              → min_pts = -1
/// FO:    <prefix>_<minPts>_<rank>.pdb     → min_pts >= 0  (BindingMode dual-suffix)
struct EmittedClusterHead {
    std::string path;
    int rank{-1};     ///< trailing rank index (0..N)
    int min_pts{-1};  ///< FO minPts; -1 for single-suffix CF/DP
};

/// Enumerate cluster-head PDBs for one restart prefix.
/// Deduplicated. Single-suffix ranks first, then FO dual-suffix directory scan.
/// Used by election AND oracle BCR so FO packaging does not leave BCR at -1.
std::vector<EmittedClusterHead>
enumerate_emitted_cluster_heads(const std::string& out_prefix);

/// True if any CF/DP or FO dual-suffix cluster head exists for out_prefix.
bool has_emitted_cluster_heads(const std::string& out_prefix);

} // namespace dataset
