// ============================================================================
//  MAKE COMPETITIVE BINDING GREAT AGAIN
// ============================================================================
// GrandCanonicalEngine.h — Grand-canonical (μVT) engine for FlexAIDdS
//
// Thermodynamic justification
// ───────────────────────────
// Two equivalent μVT constructions are supported:
//
// (A) Particle-number channel (classical GC for occupancy N = 0,1,…,N_max):
//
//     Ξ(μ,V,T) = Σ_{N=0}^{N_max} λ^N Z(N,V,T),   λ = exp(βμ)
//
//     ⟨N⟩ = (1/β) ∂lnΞ/∂μ = (1/Ξ) Σ_N N λ^N Z_N
//
//     In log space: ln Ξ = log_sum_exp_N( N·ln λ + ln Z_N )
//     using Shannon-stabilized log-sum-exp (AVX-512 / Metal / scalar).
//
// (B) Multi-species competitive channel (single binding site):
//
//     Ξ = 1 + Σ_i z_i Z_i ,   z_i = c_i / c° = exp(β μ_i)  (μ°=0)
//
//     Implemented by composing target::GrandPartitionFunction — the
//     production single-site competitive engine already in the tree.
//     GrandCanonicalEngine does NOT reimplement that Ξ; it exposes it
//     under the PartitionFunctionBase polymorphic surface and adds the
//     multi-N channel + outer OpenMP summation for large N_max.
//
// Entropy diagnostics (mixing entropy, ligand entropy collapse, isotherms)
// live on GrandPartitionFunction and are forwarded here.
//
// Ranking guardrail: μVT updates concentrations / fugacities only. CF GA
// ranking is unchanged (see VoronoiCFBatch / BindingMode comments).
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "PartitionFunctionBase.h"
#include "GrandPartitionFunction.h"

#include <cstddef>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace flexaids {

/// Grand-canonical engine: multi-N occupancy + multi-species competitive Ξ.
///
/// Inherit PartitionFunctionBase so ledger / audit code can treat NVT adapters
/// and μVT engines uniformly via log_partition() / free_energy().
class GrandCanonicalEngine final : public PartitionFunctionBase {
public:
    explicit GrandCanonicalEngine(double temperature_K = 300.0);

    // ── PartitionFunctionBase ──────────────────────────────────────────
    [[nodiscard]] double temperature() const noexcept override;
    [[nodiscard]] double log_partition() const override;  // ln Ξ (active channel)
    [[nodiscard]] EnsembleKind ensemble() const noexcept override {
        return EnsembleKind::muVT;
    }
    [[nodiscard]] std::string_view engine_name() const noexcept override {
        return "GrandCanonicalEngine";
    }

    // ── (A) Multi-N occupancy channel ──────────────────────────────────
    //
    // Register canonical Z(N,V,T) for occupancy number N.
    // Thermodynamic: Z_N comes from NVT StatMechEngine on the N-particle
    // ensemble (or a synthetic toy value in unit tests).
    void set_canonical_log_Z(int N, double log_Z_N);

    /// Register from a live StatMechEngine (N-particle NVT ensemble).
    void set_canonical_from_engine(int N, const statmech::StatMechEngine& eng);

    /// Ideal-solution chemical potential μ (kcal/mol), μ° = 0 convention.
    /// Fugacity λ = exp(βμ). Invalidates ln Ξ cache.
    void set_chemical_potential(double mu_kcal_mol);

    /// Set fugacity λ = exp(βμ) directly (λ > 0).
    void set_fugacity(double lambda);

    /// λ from concentration: λ = c/c° (ideal, c° = 1 M).
    void set_concentration_M(double concentration_M);

    /// ln Ξ from multi-N sum (OpenMP-friendly outer reduction when N_max large).
    [[nodiscard]] double log_Xi_multiN() const;

    /// ⟨N⟩ = (1/β) ∂lnΞ/∂μ = Σ_N N · p(N).
    [[nodiscard]] double mean_N_multiN() const;

    /// Variance of N: ⟨N²⟩ − ⟨N⟩² (fluctuation formula).
    [[nodiscard]] double var_N_multiN() const;

    /// p(N) = λ^N Z_N / Ξ for registered N.
    [[nodiscard]] double occupancy_probability(int N) const;

    /// Current fugacity λ.
    [[nodiscard]] double fugacity() const;

    /// Current chemical potential μ = kT ln λ (kcal/mol).
    [[nodiscard]] double chemical_potential() const;

    /// Highest registered N (or −1 if none).
    [[nodiscard]] int max_N() const;

    // ── (B) Multi-species competitive channel (wraps GPF) ──────────────
    /// Mutable access to the production single-site competitive engine.
    [[nodiscard]] target::GrandPartitionFunction& competitive() noexcept {
        return competitive_;
    }
    [[nodiscard]] const target::GrandPartitionFunction& competitive() const noexcept {
        return competitive_;
    }

    /// Convenience: register a competitive ligand (forwards to GPF).
    void add_competitive_ligand(const std::string& name, double log_Z,
                                double concentration_M = 1.0);

    /// Forward NRGsuite set_concentration to the competitive channel.
    void set_competitive_concentration(const std::string& name,
                                       double concentration_M);
    void set_competitive_concentrations(const std::vector<std::string>& names,
                                        const std::vector<double>& concentrations_M);

    /// ln Ξ from competitive multi-species channel.
    [[nodiscard]] double log_Xi_competitive() const;

    /// ⟨N⟩ from competitive channel (binary occupancy 0/1).
    [[nodiscard]] double mean_N_competitive() const;

    /// Apparent selectivity A/B (competitive channel).
    [[nodiscard]] double selectivity(const std::string& a,
                                     const std::string& b) const;

    /// Intrinsic log-selectivity ln(Z_A/Z_B).
    [[nodiscard]] double log_intrinsic_selectivity(const std::string& a,
                                                    const std::string& b) const;

    // ── Active channel for PartitionFunctionBase::log_partition ────────
    enum class ActiveChannel {
        MultiN,       ///< Use multi-N occupancy Ξ
        Competitive,  ///< Use multi-species GPF Ξ
    };
    void set_active_channel(ActiveChannel ch) noexcept { active_ = ch; }
    [[nodiscard]] ActiveChannel active_channel() const noexcept { return active_; }

    // ── Entropy / collapse (competitive channel diagnostics) ───────────
    [[nodiscard]] double mixing_entropy() const {
        return competitive_.mixing_entropy();
    }
    [[nodiscard]] double minus_T_mixing_entropy() const {
        return competitive_.minus_T_mixing_entropy();
    }
    [[nodiscard]] double ligand_entropy_collapse() const {
        return competitive_.ligand_entropy_collapse();
    }
    [[nodiscard]] std::vector<target::GrandPartitionFunction::OccupancyPoint>
    occupancy_vs_concentration(const std::string& titrate,
                               const std::vector<double>& concentrations_M) const {
        return competitive_.occupancy_vs_concentration(titrate, concentrations_M);
    }

    // ── Multi-ligand concentration vector (GA-adjacent metadata) ───────
    //
    // Thermodynamic: each GA individual may tag which ligand species it
    // represents and the bath concentration used when post-hoc μVT reweighting
    // is requested. This does NOT change CF scores inside VoronoiCFBatch.
    struct LigandVector {
        std::vector<std::string> names;
        std::vector<double> concentrations_M;  ///< parallel to names
        std::vector<double> log_Z;             ///< optional; empty = unknown
    };

    /// Store a multi-ligand concentration vector for ledger / reweight hooks.
    void set_ligand_vector(LigandVector vec);

    [[nodiscard]] const LigandVector& ligand_vector() const noexcept {
        return ligand_vector_;
    }

    /// Apply ligand_vector_ concentrations into the competitive channel
    /// (requires ligands already registered with matching names).
    void apply_ligand_vector_to_competitive();

private:
    double T_;
    double beta_;
    ActiveChannel active_ = ActiveChannel::Competitive;

    // Multi-N state: sparse map N → log Z_N (index by N, optional empty slots)
    std::vector<std::optional<double>> log_Z_by_N_;
    double log_lambda_ = 0.0;  // ln λ = βμ
    mutable std::optional<double> cached_log_xi_multiN_;
    mutable std::mutex mtx_;

    target::GrandPartitionFunction competitive_;
    LigandVector ligand_vector_;

    double compute_log_Xi_multiN_fresh() const;
};

} // namespace flexaids
