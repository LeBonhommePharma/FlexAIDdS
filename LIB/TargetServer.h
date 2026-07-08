// TargetServer.h — Central owner of all receptor/target-level state
//
// The TargetServer:
//   1. Validates the target structure (steric clashes, completeness, etc.)
//   2. Owns all receptor-level data (FA_Global, atoms, residues, grid)
//   3. Creates independent DockingSession handles for each ligand
//   4. Maintains a grand canonical partition function Ξ across all ligands
//   5. Accumulates cross-ligand knowledge (conformer priors, grid hotspots)
//
// Thread-safety: create_session() and register_result() are thread-safe.
// Multiple ligand sessions can run concurrently.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flexaid.h"
#include "statmech.h"
#include "GrandPartitionFunction.h"
#include "TargetKnowledgeBase.h"
#include "TargetValidation.h"

#include <string>
#include <vector>
#include <memory>
#include <atomic>
#include <optional>

namespace target {

// ── Configuration ──────────────────────────────────────────────────────
struct TargetConfig {
    double temperature_K = 300.0;
    int    n_models = 1;         // set > 1 for multi-model receptor
    double default_conc_M = 1.0; // P3: for grand canonical fugacity z = c/1M
};

// ── Per-ligand docking session ─────────────────────────────────────────
// Owns all ligand-specific state. Thread-safe: each session runs independently.
struct DockingSession {
    int         session_id;
    std::string ligand_name;
    double      log_Z = 0.0;     // set after docking completes
    bool        completed = false;
    int         n_poses = 0;
    double      conc_M = 1.0;  // P3: concentration for fugacity

    // CCBM conformer populations (empty if single model)
    std::vector<double> conformer_populations;

    // Best pose grid coordinates (for knowledge base)
    float best_center[3] = {0, 0, 0};
    double best_energy = 0.0;
};

// ════════════════════════════════════════════════════════════════════════
//  TargetServer
// ════════════════════════════════════════════════════════════════════════
class TargetServer {
public:
    /// Construct a TargetServer with the given configuration.
    explicit TargetServer(const TargetConfig& config = {});
    ~TargetServer() = default;

    // Non-copyable (owns shared state)
    TargetServer(const TargetServer&) = delete;
    TargetServer& operator=(const TargetServer&) = delete;

    // ── Configuration ──────────────────────────────────────────────────
    double temperature() const noexcept { return config_.temperature_K; }
    int    n_models() const noexcept { return config_.n_models; }
    double default_conc_M() const noexcept { return config_.default_conc_M; }  // P3

    // ── Validation ─────────────────────────────────────────────────────
    /// Validate a target structure given its FA_Global and atom arrays.
    /// Does NOT take ownership of any pointers.
    validation::ValidationResult validate(const FA_Global* FA,
                                           const atom* atoms,
                                           const resid* residue,
                                           int num_grd) const;

    // ── Session management (thread-safe) ───────────────────────────────

    /// Create a new docking session for a ligand.
    DockingSession create_session(const std::string& ligand_name);

    /// Register a completed session's results.
    /// Updates grand partition function and knowledge base.
    void register_result(const DockingSession& session);

    // ── Grand Partition Function ───────────────────────────────────────

    const GrandPartitionFunction& grand_partition() const { return grand_xi_; }

    /// p(ligand_i bound) = Z_i / Ξ
    double binding_probability(const std::string& ligand_name) const;

    /// K_a / K_b = Z_a / Z_b
    double selectivity_ratio(const std::string& a, const std::string& b) const;

    /// Rank all ligands by ΔG (ascending).
    std::vector<GrandPartitionFunction::LigandRank> rank_ligands() const;

    /// Number of completed ligand sessions.
    int completed_sessions() const;

    // ── Cross-ligand knowledge ─────────────────────────────────────────

    const TargetKnowledgeBase& knowledge_base() const { return knowledge_; }

    /// Bayesian posterior over receptor conformer populations.
    std::vector<double> conformer_priors() const;

    /// Mean energy per grid point accumulated across all ligands.
    std::vector<float> grid_hotspot_energies() const;

private:
    TargetConfig config_;
    GrandPartitionFunction grand_xi_;
    TargetKnowledgeBase knowledge_;
    mutable std::mutex knowledge_mtx_;    // protects knowledge_ reads/writes
    std::atomic<int> next_session_id_{0};
    std::atomic<int> completed_count_{0};
};

} // namespace target
