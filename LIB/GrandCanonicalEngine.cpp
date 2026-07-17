// ============================================================================
//  MAKE COMPETITIVE BINDING GREAT AGAIN
// ============================================================================
// GrandCanonicalEngine.cpp — Grand-canonical (μVT) implementation
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "GrandCanonicalEngine.h"
#include "UnifiedHardwareDispatch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace flexaids {

GrandCanonicalEngine::GrandCanonicalEngine(double temperature_K)
    : T_(temperature_K)
    , beta_(1.0 / (statmech::kB_kcal * temperature_K))
    , competitive_(temperature_K)
{
    require_positive_temperature(T_);
}

double GrandCanonicalEngine::temperature() const noexcept { return T_; }

double GrandCanonicalEngine::log_partition() const
{
    switch (active_) {
    case ActiveChannel::MultiN:
        return log_Xi_multiN();
    case ActiveChannel::Competitive:
        return log_Xi_competitive();
    }
    return log_Xi_competitive();
}

// ── Multi-N channel ────────────────────────────────────────────────────

void GrandCanonicalEngine::set_canonical_log_Z(int N, double log_Z_N)
{
    if (N < 0)
        throw std::invalid_argument("set_canonical_log_Z: N must be >= 0");
    std::scoped_lock lock(mtx_);
    if (static_cast<std::size_t>(N) >= log_Z_by_N_.size())
        log_Z_by_N_.resize(static_cast<std::size_t>(N) + 1);
    log_Z_by_N_[static_cast<std::size_t>(N)] = log_Z_N;
    cached_log_xi_multiN_.reset();
}

void GrandCanonicalEngine::set_canonical_from_engine(
    int N, const statmech::StatMechEngine& eng)
{
    if (eng.size() == 0)
        throw std::invalid_argument(
            "set_canonical_from_engine: empty ensemble");
    // Thermodynamic: Z_N is the NVT partition function of the N-particle
    // conformational ensemble at the same T (must match engine T ideally).
    auto th = eng.compute();
    set_canonical_log_Z(N, th.log_Z);
}

void GrandCanonicalEngine::set_chemical_potential(double mu_kcal_mol)
{
    std::scoped_lock lock(mtx_);
    // λ = exp(βμ)  ⇒  ln λ = βμ
    log_lambda_ = beta_ * mu_kcal_mol;
    cached_log_xi_multiN_.reset();
}

void GrandCanonicalEngine::set_fugacity(double lambda)
{
    if (!(lambda > 0.0))
        throw std::domain_error("set_fugacity: λ must be > 0");
    std::scoped_lock lock(mtx_);
    log_lambda_ = std::log(lambda);
    cached_log_xi_multiN_.reset();
}

void GrandCanonicalEngine::set_concentration_M(double concentration_M)
{
    // Ideal solution, c° = 1 M: λ = c/c° = c (numerical when c°=1).
    if (!(concentration_M > 0.0))
        throw std::domain_error("set_concentration_M: c must be > 0");
    if (concentration_M > 1e3)
        throw std::invalid_argument(
            "set_concentration_M: c > 1000 M — convert µM/nM to M first");
    set_fugacity(concentration_M / target::c_standard);
}

double GrandCanonicalEngine::fugacity() const
{
    std::scoped_lock lock(mtx_);
    return std::exp(log_lambda_);
}

double GrandCanonicalEngine::chemical_potential() const
{
    std::scoped_lock lock(mtx_);
    return log_lambda_ / beta_;
}

int GrandCanonicalEngine::max_N() const
{
    std::scoped_lock lock(mtx_);
    for (int i = static_cast<int>(log_Z_by_N_.size()) - 1; i >= 0; --i) {
        if (log_Z_by_N_[static_cast<std::size_t>(i)].has_value())
            return i;
    }
    return -1;
}

double GrandCanonicalEngine::compute_log_Xi_multiN_fresh() const
{
    // Ξ = Σ_N λ^N Z_N = Σ_N exp(N ln λ + ln Z_N)
    // Numerics: Shannon-stabilized log-sum-exp over registered N terms.
    //
    // When N_max is large, OpenMP parallelizes the max-find and the
    // exp-sum reduction — same pattern as UnifiedHardwareDispatch LSE.

    std::vector<double> terms;
    terms.reserve(log_Z_by_N_.size());
    for (std::size_t N = 0; N < log_Z_by_N_.size(); ++N) {
        if (!log_Z_by_N_[N].has_value())
            continue;
        // term = N·ln λ + ln Z_N
        terms.push_back(static_cast<double>(N) * log_lambda_ + *log_Z_by_N_[N]);
    }
    if (terms.empty()) {
        // No registered N: empty system Ξ = 1 (only vacuum), ln Ξ = 0
        return 0.0;
    }

    // Prefer hardware-dispatched log-sum-exp (AVX-512 / Metal / scalar).
    // For small N_max the dispatch overhead is negligible vs correctness.
    return flexaids::log_sum_exp_dispatch(terms);
}

double GrandCanonicalEngine::log_Xi_multiN() const
{
    std::scoped_lock lock(mtx_);
    if (!cached_log_xi_multiN_.has_value())
        cached_log_xi_multiN_ = compute_log_Xi_multiN_fresh();
    return *cached_log_xi_multiN_;
}

double GrandCanonicalEngine::mean_N_multiN() const
{
    // ⟨N⟩ = (1/Ξ) Σ_N N λ^N Z_N
    //     = Σ_N N · exp(N ln λ + ln Z_N − ln Ξ)
    std::scoped_lock lock(mtx_);
    if (!cached_log_xi_multiN_.has_value())
        cached_log_xi_multiN_ = compute_log_Xi_multiN_fresh();
    const double log_xi = *cached_log_xi_multiN_;

    double mean = 0.0;
#ifdef _OPENMP
#  pragma omp parallel for reduction(+ : mean) schedule(static) if (log_Z_by_N_.size() > 32)
#endif
    for (std::ptrdiff_t N = 0; N < static_cast<std::ptrdiff_t>(log_Z_by_N_.size()); ++N) {
        if (!log_Z_by_N_[static_cast<std::size_t>(N)].has_value())
            continue;
        const double log_w =
            static_cast<double>(N) * log_lambda_ + *log_Z_by_N_[static_cast<std::size_t>(N)]
            - log_xi;
        mean += static_cast<double>(N) * std::exp(log_w);
    }
    return mean;
}

double GrandCanonicalEngine::var_N_multiN() const
{
    // Var(N) = ⟨N²⟩ − ⟨N⟩²
    std::scoped_lock lock(mtx_);
    if (!cached_log_xi_multiN_.has_value())
        cached_log_xi_multiN_ = compute_log_Xi_multiN_fresh();
    const double log_xi = *cached_log_xi_multiN_;

    double mean = 0.0;
    double mean2 = 0.0;
    for (std::size_t N = 0; N < log_Z_by_N_.size(); ++N) {
        if (!log_Z_by_N_[N].has_value())
            continue;
        const double p =
            std::exp(static_cast<double>(N) * log_lambda_ + *log_Z_by_N_[N] - log_xi);
        mean += static_cast<double>(N) * p;
        mean2 += static_cast<double>(N) * static_cast<double>(N) * p;
    }
    return mean2 - mean * mean;
}

double GrandCanonicalEngine::occupancy_probability(int N) const
{
    if (N < 0)
        throw std::invalid_argument("occupancy_probability: N must be >= 0");
    std::scoped_lock lock(mtx_);
    if (static_cast<std::size_t>(N) >= log_Z_by_N_.size()
        || !log_Z_by_N_[static_cast<std::size_t>(N)].has_value()) {
        return 0.0;
    }
    if (!cached_log_xi_multiN_.has_value())
        cached_log_xi_multiN_ = compute_log_Xi_multiN_fresh();
    return std::exp(static_cast<double>(N) * log_lambda_
                    + *log_Z_by_N_[static_cast<std::size_t>(N)]
                    - *cached_log_xi_multiN_);
}

// ── Competitive channel ────────────────────────────────────────────────

void GrandCanonicalEngine::add_competitive_ligand(const std::string& name,
                                                  double log_Z,
                                                  double concentration_M)
{
    competitive_.add_ligand(name, log_Z, concentration_M);
}

void GrandCanonicalEngine::set_competitive_concentration(const std::string& name,
                                                         double concentration_M)
{
    competitive_.set_concentration(name, concentration_M);
}

void GrandCanonicalEngine::set_competitive_concentrations(
    const std::vector<std::string>& names,
    const std::vector<double>& concentrations_M)
{
    competitive_.set_concentrations(names, concentrations_M);
}

double GrandCanonicalEngine::log_Xi_competitive() const
{
    return competitive_.log_Xi();
}

double GrandCanonicalEngine::mean_N_competitive() const
{
    return competitive_.mean_N();
}

double GrandCanonicalEngine::selectivity(const std::string& a,
                                         const std::string& b) const
{
    return competitive_.selectivity(a, b);
}

double GrandCanonicalEngine::log_intrinsic_selectivity(const std::string& a,
                                                        const std::string& b) const
{
    return competitive_.log_intrinsic_selectivity(a, b);
}

// ── Ligand vector ──────────────────────────────────────────────────────

void GrandCanonicalEngine::set_ligand_vector(LigandVector vec)
{
    if (vec.names.size() != vec.concentrations_M.size())
        throw std::invalid_argument(
            "LigandVector: names and concentrations_M size mismatch");
    if (!vec.log_Z.empty() && vec.log_Z.size() != vec.names.size())
        throw std::invalid_argument(
            "LigandVector: log_Z size must match names when provided");
    ligand_vector_ = std::move(vec);
}

void GrandCanonicalEngine::apply_ligand_vector_to_competitive()
{
    if (ligand_vector_.names.empty())
        return;
    competitive_.set_concentrations(ligand_vector_.names,
                                    ligand_vector_.concentrations_M);
}

} // namespace flexaids
