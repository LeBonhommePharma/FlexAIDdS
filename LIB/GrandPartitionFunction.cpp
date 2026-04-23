// GrandPartitionFunction.cpp — Grand canonical partition function
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "GrandPartitionFunction.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace target {

GrandPartitionFunction::GrandPartitionFunction(double temperature_K)
    : T_(temperature_K)
    , beta_(1.0 / (statmech::kB_kcal * temperature_K))
{
    if (temperature_K <= 0.0)
        throw std::invalid_argument("Temperature must be positive");
}

// ── Ligand registration ────────────────────────────────────────────────

void GrandPartitionFunction::add_ligand(const std::string& name, double log_Z,
                                         double concentration_M)
{
    // concentration_M: molar concentration [M]. Standard reference state = 1.0 M.
    if (concentration_M <= 0.0)
        throw std::domain_error("concentration_M must be > 0 (got "
                                + std::to_string(concentration_M) + " M)");
    if (concentration_M > 1e3)
        throw std::invalid_argument(
            "Concentration > 1000 M — did you pass µM or nM without conversion to M?");
    double log_c = std::log(concentration_M);
    double log_zZ = log_c + log_Z;

    std::scoped_lock lock(mtx_);
    if (ligands_.count(name))
        throw std::invalid_argument("Ligand '" + name + "' already registered");
    ligands_[name] = {log_Z, log_c, log_zZ};
    cached_log_xi_.reset();
}

void GrandPartitionFunction::add_ligand(const std::string& name,
                                         const statmech::StatMechEngine& engine,
                                         double concentration_M)
{
    if (engine.size() == 0)
        throw std::invalid_argument("Cannot add ligand with empty ensemble");
    auto thermo = engine.compute();
    add_ligand(name, thermo.log_Z, concentration_M);
}

void GrandPartitionFunction::add_or_overwrite(const std::string& name, double log_Z,
                                               double concentration_M)
{
    // concentration_M: molar concentration [M]. Standard reference state = 1.0 M.
    if (concentration_M <= 0.0)
        throw std::domain_error("concentration_M must be > 0 (got "
                                + std::to_string(concentration_M) + " M)");
    if (concentration_M > 1e3)
        throw std::invalid_argument(
            "Concentration > 1000 M — did you pass µM or nM without conversion to M?");
    double log_c = std::log(concentration_M);
    double log_zZ = log_c + log_Z;

    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it != ligands_.end()) {
        it->second.log_Z = log_Z;
        it->second.log_c = log_c;
        it->second.log_zZ = log_zZ;
    } else {
        ligands_[name] = {log_Z, log_c, log_zZ};
    }
    cached_log_xi_.reset();
}

void GrandPartitionFunction::overwrite_ligand(const std::string& name, double new_log_Z)
{
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");
    it->second.log_Z = new_log_Z;
    it->second.log_zZ = it->second.log_c + new_log_Z;
    cached_log_xi_.reset();
}

void GrandPartitionFunction::merge_ligand(const std::string& name, double new_log_Z)
{
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");

    double log_c = it->second.log_c;
    double a = it->second.log_zZ;
    double b = log_c + new_log_Z;
    double max_val = std::max(a, b);
    it->second.log_zZ = max_val + std::log(std::exp(a - max_val) + std::exp(b - max_val));
    it->second.log_Z = it->second.log_zZ - log_c;
    cached_log_xi_.reset();
}

void GrandPartitionFunction::remove_ligand(const std::string& name)
{
    std::scoped_lock lock(mtx_);
    if (!ligands_.erase(name))
        throw std::invalid_argument("Ligand '" + name + "' not found");
    cached_log_xi_.reset();
}

// ── Thermodynamic queries ──────────────────────────────────────────────

double GrandPartitionFunction::compute_log_Xi_fresh() const
{
    // Ξ = 1 + Σ_i z_i·Z_i = 1 + Σ_i exp(ln(z_i·Z_i))
    // ln Ξ = log_sum_exp(0, ln(z_1·Z_1), ln(z_2·Z_2), ...)
    //
    // The "0" term is ln(1) for the empty site — it MUST be included
    // in the max-search to anchor log-sum-exp correctly.

    if (ligands_.empty()) return 0.0;

    // Anchor: unoccupied state contributes log(1) = 0 to log-sum-exp.
    // This ensures Ξ ≥ 1 always (receptor can always be empty).
    // See: Hill, T.L. "An Introduction to Statistical Thermodynamics" §15.
    double max_val = 0.0;
    for (const auto& [name, entry] : ligands_)
        max_val = std::max(max_val, entry.log_zZ);

    double sum = std::exp(0.0 - max_val);  // empty site contribution
    for (const auto& [name, entry] : ligands_)
        sum += std::exp(entry.log_zZ - max_val);
    return max_val + std::log(sum);
}

double GrandPartitionFunction::log_Xi_cached() const
{
    if (!cached_log_xi_.has_value())
        cached_log_xi_ = compute_log_Xi_fresh();
    return *cached_log_xi_;
}

double GrandPartitionFunction::log_Xi() const
{
    std::scoped_lock lock(mtx_);
    return log_Xi_cached();
}

double GrandPartitionFunction::binding_probability(const std::string& name) const
{
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");
    double log_xi = log_Xi_cached();
    return std::exp(it->second.log_zZ - log_xi);
}

double GrandPartitionFunction::empty_probability() const
{
    std::scoped_lock lock(mtx_);
    return std::exp(-log_Xi_cached());
}

double GrandPartitionFunction::mean_occupancy() const
{
    // ⟨n⟩ = 1 − p(empty)  (n is binary: 0 = apo, 1 = any ligand bound)
    return 1.0 - empty_probability();
}

double GrandPartitionFunction::occupancy_variance() const
{
    // Var(n) = ⟨n²⟩ − ⟨n⟩²  = ⟨n⟩(1 − ⟨n⟩)  for binary n ∈ {0,1}
    double mu = mean_occupancy();
    return mu * (1.0 - mu);
}

double GrandPartitionFunction::F_bound(const std::string& name) const
{
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");
    return -(1.0 / beta_) * it->second.log_Z;
}

double GrandPartitionFunction::delta_G_bind(const std::string& name, double F_ref) const
{
    return F_bound(name) - F_ref;
}

double GrandPartitionFunction::selectivity(const std::string& a,
                                            const std::string& b) const
{
    double diff = log_selectivity(a, b);
    // Sentinel: DBL_MAX (not infinity — UB under -ffast-math/-ffinite-math-only)
    if (diff > 700.0)  return std::numeric_limits<double>::max();
    if (diff < -700.0) return 0.0;
    return std::exp(diff);
}

double GrandPartitionFunction::log_selectivity(const std::string& a,
                                                const std::string& b) const
{
    std::scoped_lock lock(mtx_);
    auto it_a = ligands_.find(a);
    auto it_b = ligands_.find(b);
    if (it_a == ligands_.end())
        throw std::invalid_argument("Ligand '" + a + "' not found");
    if (it_b == ligands_.end())
        throw std::invalid_argument("Ligand '" + b + "' not found");
    // ln[(z_A·Z_A) / (z_B·Z_B)] — apparent (concentration-weighted)
    return it_a->second.log_zZ - it_b->second.log_zZ;
}

double GrandPartitionFunction::log_intrinsic_selectivity(const std::string& a,
                                                          const std::string& b) const
{
    std::scoped_lock lock(mtx_);
    auto it_a = ligands_.find(a);
    auto it_b = ligands_.find(b);
    if (it_a == ligands_.end())
        throw std::invalid_argument("Ligand '" + a + "' not found");
    if (it_b == ligands_.end())
        throw std::invalid_argument("Ligand '" + b + "' not found");
    // ln(Z_A / Z_B) — intrinsic (concentration-independent)
    return it_a->second.log_Z - it_b->second.log_Z;
}

// ── Ranking ────────────────────────────────────────────────────────────

std::vector<GrandPartitionFunction::LigandRank> GrandPartitionFunction::rank() const
{
    std::scoped_lock lock(mtx_);
    double log_xi = log_Xi_cached();
    double kT = 1.0 / beta_;

    std::vector<LigandRank> ranks;
    ranks.reserve(ligands_.size());
    for (const auto& [name, entry] : ligands_) {
        ranks.push_back({
            name,
            entry.log_Z,
            -kT * entry.log_Z,
            std::exp(entry.log_zZ - log_xi)
        });
    }

    std::sort(ranks.begin(), ranks.end(),
              [](const LigandRank& a, const LigandRank& b) {
                  return a.dG < b.dG;
              });

    return ranks;
}

// ── State queries ──────────────────────────────────────────────────────

int GrandPartitionFunction::num_ligands() const
{
    std::scoped_lock lock(mtx_);
    return static_cast<int>(ligands_.size());
}

bool GrandPartitionFunction::has_ligand(const std::string& name) const
{
    std::scoped_lock lock(mtx_);
    return ligands_.count(name) > 0;
}

std::vector<std::pair<std::string, double>> GrandPartitionFunction::all_log_Z() const
{
    std::scoped_lock lock(mtx_);
    std::vector<std::pair<std::string, double>> result;
    result.reserve(ligands_.size());
    for (const auto& [name, entry] : ligands_)
        result.emplace_back(name, entry.log_Z);
    return result;
}

std::vector<std::pair<std::string, double>> GrandPartitionFunction::all_log_zZ() const
{
    std::scoped_lock lock(mtx_);
    std::vector<std::pair<std::string, double>> result;
    result.reserve(ligands_.size());
    for (const auto& [name, entry] : ligands_)
        result.emplace_back(name, entry.log_zZ);
    return result;
}

} // namespace target
