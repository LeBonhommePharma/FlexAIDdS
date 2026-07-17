// GrandPartitionFunction.cpp — Grand canonical partition function
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "GrandPartitionFunction.h"
#include "UnifiedHardwareDispatch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

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
    double log_c = std::log(concentration_M / c_standard);  // ln(c_i / c°)
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
    double log_c = std::log(concentration_M / c_standard);  // ln(c_i / c°)
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

void GrandPartitionFunction::set_concentration(const std::string& name,
                                               double concentration_M)
{
    // Thermodynamic: only the fugacity z = c/c° changes; Z is fixed from
    // the canonical ensemble. In-place update avoids remove/re-add churn.
    if (concentration_M <= 0.0)
        throw std::domain_error("concentration_M must be > 0 (got "
                                + std::to_string(concentration_M) + " M)");
    if (concentration_M > 1e3)
        throw std::invalid_argument(
            "Concentration > 1000 M — did you pass µM or nM without conversion to M?");

    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");

    it->second.log_c = std::log(concentration_M / c_standard);
    it->second.log_zZ = it->second.log_c + it->second.log_Z;
    cached_log_xi_.reset();
}

void GrandPartitionFunction::set_concentrations(
    const std::vector<std::string>& names,
    const std::vector<double>& concentrations_M)
{
    if (names.size() != concentrations_M.size())
        throw std::invalid_argument(
            "set_concentrations: names and concentrations size mismatch");
    for (std::size_t i = 0; i < names.size(); ++i)
        set_concentration(names[i], concentrations_M[i]);
}

// ── Thermodynamic queries ──────────────────────────────────────────────

double GrandPartitionFunction::compute_log_Xi_fresh() const
{
    // Ξ = 1 + Σ_i z_i·Z_i = 1 + Σ_i exp(ln(z_i·Z_i))
    // ln Ξ = log_sum_exp(0, ln(z_1·Z_1), ln(z_2·Z_2), ...)
    //
    // Thermodynamic: the "0" term is ln(1) for the empty (apo) site — it MUST
    // be included so Ξ ≥ 1 always (Hill, Statistical Thermodynamics §15).
    //
    // Numerics: reuse Shannon-stabilized log-sum-exp (AVX-512 / Metal /
    // OpenMP / scalar) via UnifiedHardwareDispatch — same kernel as
    // StatMechEngine::log_sum_exp for NVT Z.

    if (ligands_.empty()) return 0.0;

    std::vector<double> terms;
    terms.reserve(ligands_.size() + 1);
    terms.push_back(0.0);  // empty site: ln(1) = 0
    for (const auto& [name, entry] : ligands_)
        terms.push_back(entry.log_zZ);

    return flexaids::log_sum_exp_dispatch(terms);
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

double GrandPartitionFunction::chemical_potential(const std::string& name) const
{
    // Ideal solution: μ − μ° = kT ln(c/c°) = kT · log_c
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");
    return (1.0 / beta_) * it->second.log_c;
}

double GrandPartitionFunction::concentration(const std::string& name) const
{
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");
    // c = c° · exp(log_c); c° = 1 M
    return c_standard * std::exp(it->second.log_c);
}

double GrandPartitionFunction::apparent_Ki_M(const std::string& name) const
{
    // Competitive Langmuir (diagnostic): p_i / p_empty = c_i / K_i
    // ⇒ K_i = c_i · p_empty / p_i
    std::scoped_lock lock(mtx_);
    auto it = ligands_.find(name);
    if (it == ligands_.end())
        throw std::invalid_argument("Ligand '" + name + "' not found");

    const double log_xi = log_Xi_cached();
    const double p_i = std::exp(it->second.log_zZ - log_xi);
    const double p_e = std::exp(-log_xi);
    if (p_i < 1e-300)
        return std::numeric_limits<double>::max();
    const double c = c_standard * std::exp(it->second.log_c);
    return c * p_e / p_i;
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

    // ── DO NOT CHANGE TO log_c − log_c ────────────────────────────────
    // This method is *intrinsic* selectivity: ln(Z_A / Z_B) = β(ΔG_B − ΔG_A).
    // It must be independent of the concentrations stored in log_c, otherwise
    // it would just reproduce the concentration ratio and duplicate
    // log_selectivity() at c_A = c_B = 1 M.
    // Concentration-WEIGHTED selectivity lives in log_selectivity(), which
    // returns (log_c_A + log_Z_A) − (log_c_B + log_Z_B).
    // See test ConcentrationInvarianceOfIntrinsicSelectivity.
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

// ── GC entropy diagnostics ─────────────────────────────────────────────

double GrandPartitionFunction::mixing_entropy() const
{
    // S_mix = −k_B Σ_α p_α ln p_α over {empty + all ligands}.
    // Additive μVT extension of the NVT −TΔS ledger; diagnostic only.
    std::scoped_lock lock(mtx_);
    const double log_xi = log_Xi_cached();
    double S_nats = 0.0;

    const double p_empty = std::exp(-log_xi);
    if (p_empty > 1e-300)
        S_nats -= p_empty * std::log(p_empty);

    for (const auto& [name, entry] : ligands_) {
        const double p = std::exp(entry.log_zZ - log_xi);
        if (p > 1e-300)
            S_nats -= p * std::log(p);
    }
    return statmech::kB_kcal * S_nats;
}

double GrandPartitionFunction::minus_T_mixing_entropy() const
{
    return -T_ * mixing_entropy();
}

double GrandPartitionFunction::ligand_entropy_collapse() const
{
    // Collapse on *bound* species only (empty excluded).
    // S_lig / ln(M) → 0 when one ligand dominates → collapse → 1.
    std::scoped_lock lock(mtx_);
    const int M = static_cast<int>(ligands_.size());
    if (M <= 1)
        return 0.0;

    const double log_xi = log_Xi_cached();
    const double p_empty = std::exp(-log_xi);
    const double p_bound = 1.0 - p_empty;
    if (p_bound < 1e-15)
        return 0.0;  // fully apo — no ligand collapse to report

    double S_nats = 0.0;
    for (const auto& [name, entry] : ligands_) {
        const double p_tilde = std::exp(entry.log_zZ - log_xi) / p_bound;
        if (p_tilde > 1e-300)
            S_nats -= p_tilde * std::log(p_tilde);
    }
    const double S_max = std::log(static_cast<double>(M));
    if (S_max < 1e-15)
        return 0.0;
    double collapse = 1.0 - S_nats / S_max;
    if (collapse < 0.0) collapse = 0.0;
    if (collapse > 1.0) collapse = 1.0;
    return collapse;
}

std::vector<GrandPartitionFunction::OccupancyPoint>
GrandPartitionFunction::occupancy_vs_concentration(
    const std::string& titrate,
    const std::vector<double>& concentrations_M) const
{
    // Snapshot state, then sweep titrate concentration on a temporary GPF
    // so this query is const and thread-safe w.r.t. the live ensemble.
    std::vector<std::pair<std::string, LigandEntry>> snapshot;
    {
        std::scoped_lock lock(mtx_);
        if (!ligands_.count(titrate))
            throw std::invalid_argument("Titrate ligand '" + titrate + "' not found");
        snapshot.reserve(ligands_.size());
        for (const auto& [n, e] : ligands_)
            snapshot.emplace_back(n, e);
    }

    std::vector<OccupancyPoint> curve;
    curve.reserve(concentrations_M.size());

    for (double c : concentrations_M) {
        if (c <= 0.0)
            continue;
        GrandPartitionFunction tmp(T_);
        for (const auto& [n, e] : snapshot) {
            const double ci = (n == titrate) ? c : (c_standard * std::exp(e.log_c));
            tmp.add_ligand(n, e.log_Z, ci);
        }
        OccupancyPoint pt;
        pt.concentration_M = c;
        pt.p_bound = tmp.mean_occupancy();
        pt.p_species = tmp.binding_probability(titrate);
        pt.mean_N = tmp.mean_N();
        curve.push_back(pt);
    }
    return curve;
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
