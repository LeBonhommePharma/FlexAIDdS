// SoftBetaFreeEnergy.h — single ranking objective for FlexAIDdS + DatasetRunner
//
// Soft-β free energy on the CF/contact-function scoring proxy (arbitrary units).
// Mathematically identical formulations (local members only):
//
//   ACF  = E_min − T · ln Σ_i exp(−(E_i − E_min)/T)     [cluster.cpp]
//   G̃   = H̃ − T · S̃
//        H̃ = Σ p_i E_i ,  S̃ = −Σ p_i ln p_i
//        p_i = exp(−(E_i − E_min)/T) / Z
//
// Proof: G̃ = E_min − T ln Z = ACF.  β = 1/T (kelvin), NOT 1/(k_B T).
//
// Used by:
//   - LIB/cluster.cpp          (ACF emission order when TEMPER > 0)
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
#include <limits>
#include <utility>
#include <vector>

namespace flexaids {
namespace soft_beta {

struct FreeEnergy {
    double G{0.0};  ///< H̃ − T·S̃  (lower is better) == ACF
    double H{0.0};  ///< Σ p_i E_i
    double S{0.0};  ///< −Σ p_i ln p_i  (nats)
    double Z{0.0};  ///< partition sum in shifted frame (exp terms only)
    double Emin{0.0};
    int    n{0};
};

/// Soft-β free energy over a list of CF values (cluster members).
/// Empty → G = +∞.  Single member → G = E, S = 0.
inline FreeEnergy free_energy(const std::vector<double>& energies, double T_K) noexcept
{
    FreeEnergy out;
    if (energies.empty()) {
        out.G = std::numeric_limits<double>::infinity();
        return out;
    }
    const double T = (T_K > 1e-12) ? T_K : 1e-12;
    double Emin = energies[0];
    for (double e : energies) {
        if (std::isfinite(e) && e < Emin)
            Emin = e;
    }
    if (!std::isfinite(Emin)) {
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
    // ACF form (numerically stable)
    out.G = Emin - T * std::log(Z);
    // Explicit H̃, S̃ for logging / 3Dsig identity checks
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
    // Prefer ACF form for G (identical to H−TS analytically; more stable)
    return out;
}

/// Convenience: G only (lower better).
inline double free_energy_G(const std::vector<double>& energies, double T_K) noexcept
{
    return free_energy(energies, T_K).G;
}

/// ACF alias (cluster.cpp naming).
inline double acf(const std::vector<double>& energies, double T_K) noexcept
{
    return free_energy_G(energies, T_K);
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
inline FreeEnergy mode_free_energy(const ModeCandidate& m, double T_K) noexcept
{
    if (!m.member_cfs.empty())
        return free_energy(m.member_cfs, T_K);
    if (std::isfinite(m.cf))
        return free_energy(std::vector<double>{m.cf}, T_K);
    FreeEnergy bad;
    bad.G = std::numeric_limits<double>::infinity();
    return bad;
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

/// Index of min Softβ G̃ among modes at soft temperature T_K. Empty → npos.
inline std::size_t elect_softbeta(const std::vector<ModeCandidate>& modes,
                                  double T_K) noexcept
{
    std::size_t best = static_cast<std::size_t>(-1);
    double best_G = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < modes.size(); ++i) {
        const double G = mode_free_energy(modes[i], T_K).G;
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
            double T_K,
            bool softbeta_election_enabled) noexcept
{
    if (softbeta_election_enabled) {
        const std::size_t i = elect_softbeta(modes, T_K);
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

/// Resolve soft-β T: env override > dock TEMPER > fallback (default 298).
/// TEMPER is FlexAID soft temperature (β=1/T on CF a.u.), **not** k_B·T in kcal.
inline double resolve_soft_T(double election_soft_T_env,
                             double dock_temperature_K,
                             double fallback_K = 298.0) noexcept
{
    if (election_soft_T_env > 0.0)
        return election_soft_T_env;
    if (dock_temperature_K > 0.0)
        return dock_temperature_K;
    return (fallback_K > 1e-12) ? fallback_K : 298.0;
}

}  // namespace soft_beta
}  // namespace flexaids
