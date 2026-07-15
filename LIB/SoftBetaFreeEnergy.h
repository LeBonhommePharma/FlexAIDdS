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
//   - LIB/cluster.cpp          (ACF emission order)
//   - LIB/BindingMode.cpp      (mode F_conf; vib may be added on top)
//   - LIB/DatasetRunner.cpp    (S1 election across restarts)
//
// AGENTS.md: CF is a scoring proxy; this is not experimental ΔG unless a full
// validated thermodynamic ledger is active.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <limits>
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

}  // namespace soft_beta
}  // namespace flexaids
