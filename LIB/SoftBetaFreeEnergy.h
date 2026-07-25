// SoftBetaFreeEnergy.h — single ranking objective for FlexAIDdS + DatasetRunner
//
// Soft-β free energy on the CF/contact-function scoring proxy (arbitrary units).
// Mathematically identical formulations (local members only):
//
//   ACF  = E_min − T_soft · ln Σ_i exp(−(E_i − E_min)/T_soft)   [legacy name]
//   G̃   = H̃ − T_soft · S̃
//        H̃ = Σ p_i E_i ,  S̃ = −Σ p_i ln p_i
//        p_i = exp(−(E_i − E_min)/T_soft) / Z
//
// Proof: G̃ = E_min − T_soft ln Z = ACF.
//
// **T_soft is dimensionless in CF arbitrary units (CF_AU).** It is NOT physical
// Kelvin and NOT 1/(k_B T). Historical labels "T_K" / "temperature_K" mean the
// soft temperature scale on CF scores. Classic ACF is kept as a legacy alias
// for diagnostics only; production strict re-ranking prefers
// free_energy_strict() (duplicate-invariant).
//
// Used by:
//   - LIB/cluster.cpp          (ACF emission order when TEMPER > 0 — legacy)
//   - LIB/BindingMode.cpp      (mode F_conf; vib may be added on top)
//   - LIB/DatasetRunner.cpp    (S1 election across restarts — **feature-flagged**)
//
// ── Policy (post pilot8 Softβ postmortem) ──────────────────────────────────
// Softβ is an **election / reordering** method over already-sampled modes.
// It is **not** a sampling method and **cannot** create ≤2 Å poses if the
// ensemble has BCR=0 (no near-native among emitted heads). DatasetRunner S1
// Softβ election defaults OFF (`FLEXAIDDS_SOFTBETA_ELECTION=0` /
// `FLEXAIDDS_ELECTION_SHANNON_F=0`). Engine TEMPER>0 ACF emission is a separate
// classic FlexAID product path — do not equate arm-B FO clustering with
// DatasetRunner Softβ S1 rescoring of CF ensembles.
//
// AGENTS.md: CF is a scoring proxy; this is not experimental ΔG unless a full
// validated thermodynamic ledger is active. Prefer "CF soft-β ranking proxy
// G̃" over "true binding free energy ΔG".
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

namespace flexaids {
namespace soft_beta {

struct FreeEnergy {
    double G{0.0};  ///< H̃ − T_soft·S̃  (lower is better) == ACF form
    double H{0.0};  ///< Σ p_i E_i
    double S{0.0};  ///< −Σ p_i ln p_i  (nats)
    double Z{0.0};  ///< partition sum in shifted frame (exp terms only)
    double Emin{0.0};
    int    n{0};
};

/// Soft temperature on CF arbitrary units (CF_AU). Not Kelvin.
using SoftT = double;

/// Soft-β free energy over a list of CF values (cluster members).
/// Emin is taken from **finite** entries only (NaN-first must not poison).
/// Empty / no finite → G = +∞.  Single finite member → G = E, S = 0.
///
/// Parameter name: T_soft (CF_AU). Historical callers may still pass "T_K"
/// variables — the value is always soft temperature, never physical Kelvin.
inline FreeEnergy free_energy(const std::vector<double>& energies,
                              SoftT T_soft) noexcept
{
    FreeEnergy out;
    if (energies.empty()) {
        out.G = std::numeric_limits<double>::infinity();
        return out;
    }
    const double T = (T_soft > 1e-12) ? T_soft : 1e-12;

    // Emin from finite entries only — never seed from energies[0] (may be NaN).
    double Emin = std::numeric_limits<double>::infinity();
    int n_finite = 0;
    for (double e : energies) {
        if (!std::isfinite(e))
            continue;
        ++n_finite;
        if (e < Emin)
            Emin = e;
    }
    if (n_finite == 0 || !std::isfinite(Emin)) {
        out.G = std::numeric_limits<double>::infinity();
        return out;
    }

    double Z = 0.0;
    int n = 0;
    for (double e : energies) {
        if (!std::isfinite(e))
            continue;
        Z += std::exp(-(e - Emin) / T);
        ++n;
    }
    out.n = n;
    out.Emin = Emin;
    out.Z = Z;
    if (!(Z > 0.0) || !std::isfinite(Z) || n == 0) {
        out.G = Emin;
        out.H = Emin;
        out.S = 0.0;
        return out;
    }
    // ACF form (numerically stable) — legacy identity; G̃ = Emin − T ln Z.
    out.G = Emin - T * std::log(Z);
    double H = 0.0;
    double S = 0.0;
    for (double e : energies) {
        if (!std::isfinite(e))
            continue;
        const double p = std::exp(-(e - Emin) / T) / Z;
        if (p <= 0.0)
            continue;
        H += p * e;
        S -= p * std::log(p);
    }
    out.H = H;
    out.S = S;
    return out;
}

/// Convenience: G only (lower better).
inline double free_energy_G(const std::vector<double>& energies,
                            SoftT T_soft) noexcept
{
    return free_energy(energies, T_soft).G;
}

/// Legacy ACF alias (cluster.cpp naming). Diagnostic only — prefer
/// free_energy_strict() for claim re-ranking (duplicate-invariant).
inline double acf(const std::vector<double>& energies, SoftT T_soft) noexcept
{
    return free_energy_G(energies, T_soft);
}

/// Strict re-rank modes for duplicate-invariant Softβ.
enum class StrictRerankMode {
    /// Collapse exact equal CF values then classic free_energy (default claim path).
    /// Near-duplicates (dense basin with slight CF variation) still contribute to Z;
    /// exact clones of the same CF do not inflate Softβ.
    UniqueGeometry,
    /// Log-mean-exp over unique finite energies (multiplicity-agnostic mean exp).
    /// Prefer UniqueGeometry for mode election; LogMeanExp for diagnostics.
    LogMeanExp,
};

/// Duplicate-invariant Softβ free energy.
/// Exact CF duplicates (cloned members / re-emitted heads) must not deepen G̃
/// via multiplicity inflation. Default UniqueGeometry: collapse exact equal CF
/// then classic free_energy. LogMeanExp: G = Emin − T ln(mean exp) over unique.
inline FreeEnergy free_energy_strict(
    const std::vector<double>& energies,
    SoftT T_soft,
    StrictRerankMode mode = StrictRerankMode::UniqueGeometry) noexcept
{
    // Collect unique finite energies (exact equality; CF is discrete proxy).
    std::vector<double> unique;
    unique.reserve(energies.size());
    for (double e : energies) {
        if (!std::isfinite(e))
            continue;
        bool seen = false;
        for (double u : unique) {
            if (u == e) { seen = true; break; }
        }
        if (!seen)
            unique.push_back(e);
    }
    if (unique.empty()) {
        FreeEnergy out;
        out.G = std::numeric_limits<double>::infinity();
        return out;
    }
    if (mode == StrictRerankMode::UniqueGeometry) {
        return free_energy(unique, T_soft);
    }
    // LogMeanExp: multiplicity-free mean of Boltzmann factors.
    FreeEnergy out;
    const double T = (T_soft > 1e-12) ? T_soft : 1e-12;
    double Emin = unique[0];
    for (double e : unique)
        if (e < Emin) Emin = e;
    double sum_exp = 0.0;
    for (double e : unique)
        sum_exp += std::exp(-(e - Emin) / T);
    const double n = static_cast<double>(unique.size());
    const double mean_exp = sum_exp / n;
    out.n = static_cast<int>(unique.size());
    out.Emin = Emin;
    out.Z = mean_exp;  // store mean exp for diagnostics
    out.G = Emin - T * std::log(mean_exp);
    // H, S under uniform unique measure on the same Boltzmann weights (normalized).
    double Z_u = sum_exp;
    double H = 0.0, S = 0.0;
    for (double e : unique) {
        const double p = std::exp(-(e - Emin) / T) / Z_u;
        H += p * e;
        if (p > 0.0) S -= p * std::log(p);
    }
    out.H = H;
    out.S = S;
    return out;
}

// ── Cluster emission basin score (E1b / rank_miss remediation) ──────────────
//
// Between-cluster ranking in cluster.cpp uses this score (lower better).
// free_energy_strict collapses exact-CF duplicates so multiplicity alone
// cannot elect a large wrong basin over a better CF singleton (rank_miss).
//
// Env (product path = free_energy_strict ON; campaign A/B via LEGACY opt-out):
//   unset / FLEXAIDDS_ACF_STRICT=1  → free_energy_strict (default ON, E1b)
//   FLEXAIDDS_ELECT_LEGACY_ACF=1    → legacy multiplicity-inflated acf
//   FLEXAIDDS_ACF_STRICT=0         → legacy acf (explicit opt-out alias)

/// True when cluster emission should use free_energy_strict (default ON).
inline bool cluster_use_free_energy_strict_from_env() noexcept
{
    const char* legacy = std::getenv("FLEXAIDDS_ELECT_LEGACY_ACF");
    if (legacy != nullptr && std::atoi(legacy) != 0)
        return false;
    const char* strict = std::getenv("FLEXAIDDS_ACF_STRICT");
    if (strict != nullptr && std::atoi(strict) == 0)
        return false;
    // unset or ACF_STRICT!=0 → strict (product default after E1b)
    return true;
}

/// Cluster-local basin score for emission order (lower better).
/// Pure function — unit-testable without linking cluster.cpp.
inline double cluster_basin_score(const std::vector<double>& member_energies,
                                  SoftT T_soft,
                                  bool use_strict) noexcept
{
    if (use_strict)
        return free_energy_strict(member_energies, T_soft,
                                  StrictRerankMode::UniqueGeometry)
            .G;
    return acf(member_energies, T_soft);
}

/// Env-gated cluster basin score (production cluster.cpp path).
inline double cluster_basin_score_from_env(
    const std::vector<double>& member_energies,
    SoftT T_soft) noexcept
{
    return cluster_basin_score(member_energies, T_soft,
                               cluster_use_free_energy_strict_from_env());
}

// ── Gated Softβ election (DatasetRunner S1 / offline re-rank) ──────────────
//
// Production path never consults crystal RMSD. Softβ is enabled only by an
// explicit feature flag (ProtocolConfig::election_shannon_free_energy /
// FLEXAIDDS_SOFTBETA_ELECTION / FLEXAIDDS_ELECTION_SHANNON_F). When OFF,
// electors must fall back to CF rank-0 and must not claim Softβ improvement.
//
// Crystal / native-CF oracles are diagnostics only (offline ablation tables).

/// One already-clustered mode / emitted head for Softβ reordering.
struct ModeCandidate {
    double cf{std::numeric_limits<double>::infinity()};  ///< head CF
    std::vector<double> member_cfs;  ///< .mcf / cluster members; empty → {cf}
};

/// Soft-β free energy of one mode (local Z over members only).
/// Default uses free_energy_strict (duplicate-invariant LogMeanExp).
/// Pass use_strict=false only for classic ACF diagnostic comparisons.
inline FreeEnergy mode_free_energy(const ModeCandidate& m,
                                   SoftT T_soft,
                                   bool use_strict = true) noexcept
{
    const std::vector<double>* src = nullptr;
    std::vector<double> singleton;
    if (!m.member_cfs.empty()) {
        src = &m.member_cfs;
    } else if (std::isfinite(m.cf)) {
        singleton = {m.cf};
        src = &singleton;
    } else {
        FreeEnergy bad;
        bad.G = std::numeric_limits<double>::infinity();
        return bad;
    }
    if (use_strict)
        return free_energy_strict(*src, T_soft, StrictRerankMode::UniqueGeometry);
    return free_energy(*src, T_soft);  // legacy ACF diagnostic
}

/// Index of min head CF among finite candidates. Empty → npos.
inline std::size_t elect_cf_rank0(const std::vector<ModeCandidate>& modes) noexcept
{
    std::size_t best = static_cast<std::size_t>(-1);
    double best_cf = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < modes.size(); ++i) {
        if (!std::isfinite(modes[i].cf))
            continue;
        if (modes[i].cf < best_cf) {
            best_cf = modes[i].cf;
            best = i;
        }
    }
    return best;
}

/// Index of min Softβ G̃ among modes at soft temperature T_soft (CF_AU).
/// Uses duplicate-invariant strict reranker. Empty → npos.
inline std::size_t elect_softbeta(const std::vector<ModeCandidate>& modes,
                                  SoftT T_soft) noexcept
{
    std::size_t best = static_cast<std::size_t>(-1);
    double best_G = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < modes.size(); ++i) {
        const double G = mode_free_energy(modes[i], T_soft, /*strict=*/true).G;
        if (!std::isfinite(G))
            continue;
        if (G < best_G) {
            best_G = G;
            best = i;
        }
    }
    return best;
}

/// Gated election: Softβ only when flag is true; else CF rank-0.
/// Returns (index, used_softbeta). Index is npos if no finite candidate.
inline std::pair<std::size_t, bool>
elect_gated(const std::vector<ModeCandidate>& modes,
            SoftT T_soft,
            bool softbeta_election_enabled) noexcept
{
    if (softbeta_election_enabled) {
        const std::size_t i = elect_softbeta(modes, T_soft);
        return {i, i != static_cast<std::size_t>(-1)};
    }
    return {elect_cf_rank0(modes), false};
}

/// Offline / diagnostic only: true if any mode RMSD ≤ threshold (Å).
/// **Never** wire this into production election — crystal-gated ranking is
/// scientifically invalid for blind claims. Use for ablation tables that ask
/// "would Softβ have helped **if** a near-native existed?"
inline bool diagnostic_near_native_present(const std::vector<double>& mode_rmsds,
                                           double threshold_A = 2.0) noexcept
{
    for (double r : mode_rmsds) {
        if (std::isfinite(r) && r <= threshold_A)
            return true;
    }
    return false;
}

/// Offline guidance: Softβ can only reorder; if no near-native is present in
/// the ensemble (BCR=0), Softβ cannot produce ≤2 Å S1 success. Production
/// still uses the feature flag alone (never crystal).
inline bool diagnostic_softbeta_can_help_s1(
    const std::vector<double>& mode_rmsds,
    double threshold_A = 2.0) noexcept
{
    return diagnostic_near_native_present(mode_rmsds, threshold_A);
}

/// Resolve soft-β T_soft (CF_AU): env override > dock TEMPER > fallback.
/// TEMPER / historical "temperature_K" are soft temperature on CF a.u., **not**
/// physical Kelvin and **not** k_B·T in kcal.
inline SoftT resolve_soft_T(double election_soft_T_env,
                            double dock_soft_T,
                            SoftT fallback = 298.0) noexcept
{
    if (election_soft_T_env > 0.0)
        return election_soft_T_env;
    if (dock_soft_T > 0.0)
        return dock_soft_T;
    return (fallback > 1e-12) ? fallback : 298.0;
}

}  // namespace soft_beta
}  // namespace flexaids
