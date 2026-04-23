// GrandPartitionFunction.h — Grand canonical partition function for competitive binding
//
// For a binding site that can be empty or occupied by one of N ligands:
//
//   Ξ = 1 + z_A·Z_A + z_B·Z_B + z_C·Z_C + ...
//
// where Z_i is ligand i's canonical partition function (from StatMechEngine),
// z_i = c_i / c° is the fugacity (activity at concentration c_i, standard
// state c° = 1 M), and the "1" represents the empty (apo) site.
//
// Internally stores ln(z_i · Z_i) = ln(c_i/c°) + ln Z_i for each ligand,
// so all existing log-space arithmetic is preserved.
//
// From Ξ we compute:
//   p(empty)     = 1 / Ξ
//   p(ligand_i)  = z_i·Z_i / Ξ
//   F_bound(i)   = −kT ln Z_i          (concentration-independent)
//   selectivity(A/B) = (z_A·Z_A) / (z_B·Z_B)
//
// All arithmetic uses log-space (stores ln Z_i) for numerical stability
// since partition functions can span hundreds of orders of magnitude.
//
// Thread-safe: all mutations are mutex-protected.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "statmech.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <optional>

namespace target {

class GrandPartitionFunction {
public:
    explicit GrandPartitionFunction(double temperature_K = 300.0);

    // ── Non-copyable, non-movable ────────────────────────────────────
    // The receptor binding state (ligand map + Ξ cache, guarded by mtx_)
    // is mutable. Copying would silently duplicate that state and
    // produce divergent ensemble statistics; moving would invalidate
    // pointers held by any concurrent observer. Construct once per
    // binding site and share by reference/pointer.
    GrandPartitionFunction(const GrandPartitionFunction&) = delete;
    GrandPartitionFunction& operator=(const GrandPartitionFunction&) = delete;
    GrandPartitionFunction(GrandPartitionFunction&&) = delete;
    GrandPartitionFunction& operator=(GrandPartitionFunction&&) = delete;

    // ── Ligand registration ────────────────────────────────────────────

    /// Register a ligand with its partition function and optional concentration.
    /// @param name       Ligand identifier
    /// @param log_Z      ln(Z_i) from StatMechEngine
    /// @param concentration_M  Ligand concentration in MOLAR (M). Standard reference state = 1.0 M.
    ///                         Do NOT pass µM or nM directly — convert to M first.
    ///                         Values <= 0 or > 1000 M are rejected (domain_error).
    void add_ligand(const std::string& name, double log_Z, double concentration_M = 1.0);

    /// Convenience: extract log_Z from a StatMechEngine
    /// @param concentration_M  Ligand concentration in MOLAR (M). Standard reference state = 1.0 M.
    void add_ligand(const std::string& name, const statmech::StatMechEngine& engine,
                    double concentration_M = 1.0);

    /// Atomic insert-or-overwrite: if ligand exists, overwrite; otherwise insert.
    /// Thread-safe — avoids the TOCTOU race between has_ligand() + add/overwrite.
    /// @param concentration_M  Ligand concentration in MOLAR (M). Standard reference state = 1.0 M.
    void add_or_overwrite(const std::string& name, double log_Z, double concentration_M = 1.0);

    /// Overwrite an existing ligand's Z (e.g., after re-docking with a better estimate).
    /// Preserves the existing concentration.
    void overwrite_ligand(const std::string& name, double new_log_Z);

    /// Merge an independent ensemble into an existing ligand's Z (log-sum-exp).
    /// Assumes new_log_Z is at the same concentration as the currently registered entry.
    void merge_ligand(const std::string& name, double new_log_Z);

    /// Remove a ligand from the ensemble.
    void remove_ligand(const std::string& name);

    // ── Thermodynamic queries ──────────────────────────────────────────

    /// ln(Ξ) = ln(1 + Σ_i z_i·Z_i) using log-sum-exp for stability
    [[nodiscard]] double log_Xi() const;

    /// p(ligand_i bound) = z_i·Z_i / Ξ
    [[nodiscard]] double binding_probability(const std::string& name) const;

    /// p(empty) = 1/Ξ = exp(−ln Ξ)
    [[nodiscard]] double empty_probability() const;

    /// Mean site occupancy ⟨n⟩ = 1 − p(empty) = Σ_i p(ligand_i).
    /// For a single binding site n ∈ {0,1}; this equals the probability
    /// that the receptor is bound by ANY ligand.
    [[nodiscard]] double mean_occupancy() const;

    /// Occupancy variance Var(n) = ⟨n⟩·(1 − ⟨n⟩).
    /// Maximised at 0.25 when mean_occupancy = 0.5 (half-saturation).
    /// Useful for detecting near-saturation or near-apo conditions.
    [[nodiscard]] double occupancy_variance() const;

    /// Helmholtz free energy of the bound ensemble: F_bound = −kT · ln Z_i  (kcal/mol).
    /// This is NOT the binding free energy (which requires an unbound reference).
    /// Use delta_G_bind() to get ΔG_bind = F_bound − F_ref.
    [[nodiscard]] double F_bound(const std::string& name) const;

    /// @deprecated Use F_bound() — free_energy() was ambiguous between Helmholtz F
    ///             and binding ΔG_bind, causing silent errors in Kd calculations.
    [[deprecated("Use F_bound() — free_energy() was ambiguous")]]
    double free_energy(const std::string& name) const { return F_bound(name); }

    /// Binding free energy: ΔG_bind = F_bound − F_ref.
    /// @param name    Ligand name
    /// @param F_ref   Reference-state free energy (unbound ligand in solution).
    ///                Typically obtained from a separate StatMechEngine on the
    ///                unbound ensemble. If 0.0, returns F_bound directly.
    [[nodiscard]] double delta_G_bind(const std::string& name, double F_ref = 0.0) const;

    /// Apparent selectivity ratio: (z_A·Z_A) / (z_B·Z_B) — concentration-weighted.
    /// Returns +Inf or 0.0 for extreme ratios (|ΔΔG| > ~700 kT).
    [[nodiscard]] double selectivity(const std::string& a, const std::string& b) const;

    /// ln[(z_A·Z_A) / (z_B·Z_B)] — apparent (concentration-weighted). Overflow-safe.
    [[nodiscard]] double log_selectivity(const std::string& a, const std::string& b) const;

    /// Intrinsic log-selectivity: ln(Z_A / Z_B) = β(ΔG_B − ΔG_A).
    ///
    /// Concentration-INDEPENDENT — cancels the z_i = c_i/c° factors and
    /// depends only on the partition functions, i.e. on the binding
    /// affinities themselves. This is the correct quantity for SAR and
    /// potency-series comparisons where assay concentrations vary.
    ///
    /// Sign convention: positive ⇒ A preferred (Z_A > Z_B); negative ⇒ B preferred.
    /// Compare with log_selectivity(), which is concentration-weighted.
    [[nodiscard]] double log_intrinsic_selectivity(const std::string& a,
                                                    const std::string& b) const;

    // ── Ranking ────────────────────────────────────────────────────────

    struct LigandRank {
        std::string name;
        double log_Z;   // ln Z_i
        double dG;      // −kT ln Z_i (kcal/mol)
        double p_bound; // z_i·Z_i / Ξ
    };

    /// Rank all ligands by ΔG (ascending = most favorable first).
    [[nodiscard]] std::vector<LigandRank> rank() const;

    // ── State queries ──────────────────────────────────────────────────

    [[nodiscard]] int  num_ligands() const;
    [[nodiscard]] bool has_ligand(const std::string& name) const;
    [[nodiscard]] double temperature() const noexcept { return T_; }

    /// Get all ligand names and their intrinsic log(Z) values.
    [[nodiscard]] std::vector<std::pair<std::string, double>> all_log_Z() const;

    /// Get all ligand names and their concentration-weighted log(z·Z) values.
    [[nodiscard]] std::vector<std::pair<std::string, double>> all_log_zZ() const;

private:
    double T_;
    double beta_;    // 1/(kT)

    struct LigandEntry {
        double log_Z;   // intrinsic ln(Z_i)
        double log_c;   // ln(c_i/c°) — stored directly to avoid subtraction drift
        double log_zZ;  // ln(z_i · Z_i) = log_c + log_Z (cached for performance)
    };
    std::unordered_map<std::string, LigandEntry> ligands_;
    mutable std::mutex mtx_;
    mutable std::optional<double> cached_log_xi_;  // nullopt = dirty

    // Compute ln(Ξ) from scratch (caller must hold lock)
    double compute_log_Xi_fresh() const;

    // Return cached ln(Ξ), recomputing if dirty (caller must hold lock)
    double log_Xi_cached() const;
};

} // namespace target
