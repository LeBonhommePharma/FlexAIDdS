// TargetServer.cpp — Central owner of all receptor/target-level state
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "TargetServer.h"

namespace target {

// ────────────────────────────────────────────────────────────────────────
// Construction
// ────────────────────────────────────────────────────────────────────────

TargetServer::TargetServer(const TargetConfig& config)
    : config_(config)
    , grand_xi_(config.temperature_K)
    , knowledge_(config.n_models, 0)  // grid points unknown until validate()
{
}

// ────────────────────────────────────────────────────────────────────────
// Validation
// ────────────────────────────────────────────────────────────────────────

validation::ValidationResult TargetServer::validate(const FA_Global* FA,
                                                     const atom* atoms,
                                                     const resid* residue,
                                                     int num_grd) const
{
    return validation::run_all_checks(FA, atoms, residue, num_grd);
}

// ────────────────────────────────────────────────────────────────────────
// Session management
// ────────────────────────────────────────────────────────────────────────

DockingSession TargetServer::create_session(const std::string& ligand_name)
{
    DockingSession session;
    session.session_id = next_session_id_.fetch_add(1, std::memory_order_relaxed);
    session.ligand_name = ligand_name;
    session.completed = false;
    session.log_Z = 0.0;
    session.n_poses = 0;
    session.conc_M = config_.default_conc_M;  // P3 from config
    session.best_energy = 0.0;
    session.best_center[0] = session.best_center[1] = session.best_center[2] = 0.0f;
    return session;
}

void TargetServer::register_result(const DockingSession& session)
{
    if (!session.completed) return;

    // Register into grand partition function (atomic insert-or-overwrite)
    // P3: pass conc_M for fugacity
    grand_xi_.add_or_overwrite(session.ligand_name, session.log_Z, session.conc_M);

    // Accumulate conformer knowledge (CCBM) and binding center
    {
        std::lock_guard<std::mutex> lock(knowledge_mtx_);
        if (!session.conformer_populations.empty()) {
            knowledge_.accumulate_conformer_weights(session.conformer_populations);
        }
        knowledge_.accumulate_binding_center(
            session.best_center[0], session.best_center[1], session.best_center[2],
            session.best_energy, session.ligand_name);
    }

    completed_count_.fetch_add(1, std::memory_order_relaxed);
}

// ────────────────────────────────────────────────────────────────────────
// Grand partition function queries (GPF handles its own locking)
// ────────────────────────────────────────────────────────────────────────

double TargetServer::binding_probability(const std::string& ligand_name) const
{
    return grand_xi_.binding_probability(ligand_name);
}

double TargetServer::selectivity_ratio(const std::string& a,
                                        const std::string& b) const
{
    return grand_xi_.selectivity(a, b);
}

std::vector<GrandPartitionFunction::LigandRank> TargetServer::rank_ligands() const
{
    return grand_xi_.rank();
}

int TargetServer::completed_sessions() const
{
    return completed_count_.load(std::memory_order_relaxed);
}

// ────────────────────────────────────────────────────────────────────────
// Knowledge base queries
// ────────────────────────────────────────────────────────────────────────

std::vector<double> TargetServer::conformer_priors() const
{
    std::lock_guard<std::mutex> lock(knowledge_mtx_);
    return knowledge_.conformer_posterior();
}

std::vector<float> TargetServer::grid_hotspot_energies() const
{
    std::lock_guard<std::mutex> lock(knowledge_mtx_);
    return knowledge_.grid_mean_energies();
}

} // namespace target
