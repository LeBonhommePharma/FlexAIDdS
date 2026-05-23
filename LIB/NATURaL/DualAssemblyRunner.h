// DualAssemblyRunner.h — Top-level cotranslational docking driver
//
// Orchestrates Sim A / Sim B / Sim C across NATURaL tracks with the Shannon-driven
// pose-collapse role discriminator. The GA backend is injected as a std::function so
// the runner has no compile-time dependency on the FlexAID GA machinery. Unit tests
// and the current CLI inject synthetic engines; production GA wiring is a separate
// backend integration step.
//
// References: docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md §3, §5.3
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#pragma once

#include "FibrilGrowthOracle.h"
#include "NascentChainScheduler.h"
#include "../statmech.h"

#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace natural {

// ─── Configuration ───────────────────────────────────────────────────────────
struct DualAssemblyConfig {
    std::string protofibril_pdb;            // required
    std::string monomer_pdb;                // required if sim_c_enabled
    std::string sequence_fasta;             // required (1-letter)
    int         transcript_nt              = 0;          // 0 = derived from sequence × 3
    int         checkpoint_interval        = 10;
    int         sim_c_interval             = 5;
    double      monomer_conc_M             = 1.0e-6;     // 1 µM cytosolic free monomer
    double      temperature_K              = 310.15;     // 37 °C
    int         n_threads                  = 6;
    bool        include_reciprocal_controls = true;
    bool        sim_c_enabled               = true;
    std::string output_csv                  = "cotranslational_trajectory.csv";
    std::string nascent_pdb_dir             = ".";
    double      acceptance_threshold        = 0.5;
};

// ─── GA-backend callback signatures ──────────────────────────────────────────
// The runner asks the backend to run one GA and return a populated StatMechEngine
// along with a flat vector of per-pose coordinates used to bin Shannon entropy.
// This entropy is a pose-ensemble collapse metric, not intrinsic conformational
// entropy of either isolated molecule.
struct GAResult {
    statmech::StatMechEngine engine;
    std::vector<double>      pose_rmsds; // any continuous pose coordinate; 1D is enough
};

// Sim A — target = protofibril, ligand = nascent chain at length L_k
using SimAFn = std::function<GAResult(const std::string& target_pdb,
                                       const std::string& ligand_chain_pdb,
                                       int                L_k,
                                       double             temperature_K)>;

// Sim B — target = nascent chain, ligand = monomer
using SimBFn = std::function<GAResult(const std::string& target_chain_pdb,
                                       const std::string& monomer_pdb,
                                       int                L_k,
                                       double             temperature_K)>;

// Sim C — target = protofibril, ligand = free monomer
using SimCFn = std::function<GAResult(const std::string& target_pdb,
                                       const std::string& monomer_pdb,
                                       double             temperature_K)>;

// Nascent-chain truncation helper — given the full sequence and a length L_k, return
// the absolute path of the truncated PDB (the helper writes it to disk under
// cfg.nascent_pdb_dir). MVP: idealised extended geometry.
using TruncateFn = std::function<std::string(const std::string& sequence,
                                              int                L_k,
                                              const std::string& output_dir)>;

// Optional protofibril-state updater. Called only when Sim C accepts elongation.
// Return the PDB path to use as the protofibril target for subsequent checkpoints.
using ProtofibrilAdvanceFn = std::function<std::string(const std::string& current_protofibril_pdb,
                                                       const std::string& monomer_pdb,
                                                       int                next_state_index)>;

// ─── Runner ──────────────────────────────────────────────────────────────────
class DualAssemblyRunner {
public:
    DualAssemblyRunner(DualAssemblyConfig cfg,
                       SimAFn             sim_a,
                       SimBFn             sim_b,
                       SimCFn             sim_c,
                       TruncateFn         truncate,
                       ProtofibrilAdvanceFn advance_protofibril = nullptr);

    // Run the full trajectory. Writes the CSV to cfg.output_csv as it advances and
    // returns the per-checkpoint history.
    std::vector<std::pair<Checkpoint, CheckpointOutcome>> run();

    const NascentChainScheduler& scheduler() const noexcept { return scheduler_; }

    // ── Static pose-entropy role discriminator (exposed for unit tests) ─────
    static void assign_tl(double             H_A_nats,
                          double             H_B_nats,
                          char               prev_tl,
                          CheckpointOutcome& out) noexcept;

private:
    struct TrackState {
        double H_prev_chain_nats = std::numeric_limits<double>::infinity();
        char   tl_prev           = '?';
        int    checkpoints_since_sim_c = 0;
    };

    DualAssemblyConfig                  cfg_;
    NascentChainScheduler               scheduler_;
    FibrilGrowthOracle                  oracle_;
    SimAFn                              sim_a_;
    SimBFn                              sim_b_;
    SimCFn                              sim_c_;
    TruncateFn                          truncate_;
    std::unordered_map<int, TrackState> track_state_;
    ProtofibrilAdvanceFn                advance_protofibril_;
    std::string                         current_protofibril_pdb_;
    int                                 protofibril_state_index_ = 0;

    // Shannon entropy in nats from a population of pose coordinates. Wraps
    // shannon_thermo::compute_shannon_entropy with a default bin count.
    static double shannon_entropy_nats(const std::vector<double>& values);

    // ΔG = compute().free_energy. Returns NaN if engine is empty.
    static double free_energy_kcal(const statmech::StatMechEngine& engine);

    void write_csv_header(const std::string& path) const;
    void write_csv_row(const std::string& path,
                       const Checkpoint& ck,
                       const CheckpointOutcome& out) const;
};

} // namespace natural
