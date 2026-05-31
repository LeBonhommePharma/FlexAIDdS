// ThermalExtrapolation.h — Kirchhoff / Robertson-Murphy ΔG(T) extrapolation
//
// Evaluates:
//   ΔG(T) = ΔHm · (1 − T/Tm) − ΔCp · [(Tm − T) + T · ln(T/Tm)]
//
// References:
//   Robertson & Murphy, Chem. Rev. 1997, 97, 1251–1268.  (eq. 7)
//   Murphy & Freire, Adv. Protein Chem. 1992, 43, 313–361.
//
// Inputs (all in kcal/mol and K):
//   Tm    — melting / reference temperature where ΔG(Tm) = 0 by definition
//             (from DSF, TSA, or thermal denaturation ITC)
//   ΔHm   — binding enthalpy at Tm (kcal/mol; from ITC at Tm or extrapolated)
//   ΔCp   — heat capacity change of binding (kcal/mol·K; from compute_delta_Cp())
//             Typically negative for hydrophobic burial; positive for polar burial.
//
// Unit contract: all energies kcal/mol, temperatures K.
//
// This file is header-only. No CMake changes required.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma / NRGlab, Université de Montréal.
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace thermal_extrap {

// ─── Input / output structs ───────────────────────────────────────────────────

struct KirchhoffInput {
    double Tm_K      = 0.0;   // melting / reference temperature (K) — ΔG(Tm) ≡ 0
    double delta_Hm  = 0.0;   // ΔH at Tm (kcal/mol)
    double delta_Cp  = 0.0;   // ΔCp of binding (kcal/mol·K)
};

struct KirchhoffResult {
    double T_K      = 0.0;   // target temperature (K)
    double delta_G  = 0.0;   // ΔG(T) (kcal/mol)
    double delta_H  = 0.0;   // ΔH(T) = ΔHm + ΔCp·(T − Tm)  (kcal/mol)
    double T_delta_S = 0.0;  // T·ΔS(T) = ΔH(T) − ΔG(T)     (kcal/mol)
    double delta_S  = 0.0;   // ΔS(T) = T_delta_S / T         (kcal/mol·K)
};

// ─── Core evaluator ──────────────────────────────────────────────────────────
//
// Derivation from first principles:
//
//   ΔH(T) = ΔHm + ΔCp·(T − Tm)                    [Kirchhoff heat content]
//
//   ΔS(T) = ΔSm + ΔCp·ln(T/Tm)
//         = (ΔHm/Tm) + ΔCp·ln(T/Tm)               [since ΔG(Tm)=0 → ΔSm = ΔHm/Tm]
//
//   ΔG(T) = ΔH(T) − T·ΔS(T)
//          = [ΔHm + ΔCp(T−Tm)] − T·[(ΔHm/Tm) + ΔCp·ln(T/Tm)]
//          = ΔHm·(1 − T/Tm) − ΔCp·[(Tm − T) + T·ln(T/Tm)]
//
// At T = Tm: ΔG = 0 by construction (verified analytically and by test).
// At ΔCp = 0: reduces to van't Hoff linear extrapolation.

inline KirchhoffResult kirchhoff_deltaG(
    const KirchhoffInput& in,
    double T_K)
{
    if (in.Tm_K <= 0.0)
        throw std::invalid_argument(
            "kirchhoff_deltaG: Tm_K must be > 0 K");
    if (T_K <= 0.0)
        throw std::invalid_argument(
            "kirchhoff_deltaG: T_K must be > 0 K");

    KirchhoffResult r;
    r.T_K = T_K;

    // ΔH(T) = ΔHm + ΔCp·(T − Tm)
    r.delta_H = in.delta_Hm + in.delta_Cp * (T_K - in.Tm_K);

    if (T_K == in.Tm_K) {
        // Exact: ΔG(Tm) = 0 by definition
        r.delta_G   = 0.0;
        r.T_delta_S = r.delta_H;
        r.delta_S   = r.delta_H / T_K;
        return r;
    }

    // ΔG(T) = ΔHm·(1 − T/Tm) − ΔCp·[(Tm − T) + T·ln(T/Tm)]
    r.delta_G = in.delta_Hm * (1.0 - T_K / in.Tm_K)
              - in.delta_Cp * ((in.Tm_K - T_K) + T_K * std::log(T_K / in.Tm_K));

    r.T_delta_S = r.delta_H - r.delta_G;
    r.delta_S   = r.T_delta_S / T_K;
    return r;
}

// ─── Temperature scan ────────────────────────────────────────────────────────
// Returns n_steps evenly-spaced KirchhoffResult values from T_lo_K to T_hi_K.
// Useful for plotting ΔG(T) curves and finding the stability window.

inline std::vector<KirchhoffResult> kirchhoff_scan(
    const KirchhoffInput& in,
    double T_lo_K,
    double T_hi_K,
    int    n_steps = 50)
{
    if (n_steps < 2)
        throw std::invalid_argument("kirchhoff_scan: n_steps must be >= 2");
    if (T_lo_K <= 0.0 || T_hi_K <= T_lo_K)
        throw std::invalid_argument(
            "kirchhoff_scan: need 0 < T_lo_K < T_hi_K");

    std::vector<KirchhoffResult> out;
    out.reserve(n_steps);
    const double step = (T_hi_K - T_lo_K) / static_cast<double>(n_steps - 1);
    for (int i = 0; i < n_steps; ++i)
        out.push_back(kirchhoff_deltaG(in, T_lo_K + i * step));
    return out;
}

// ─── Tm finder ───────────────────────────────────────────────────────────────
// Locates the zero-crossing of ΔG(T) in a pre-computed scan via linear
// interpolation. Returns -1.0 if no crossing is found in the scan range.
// Cross-check: find_Tm_crossing(kirchhoff_scan(in, T_lo, T_hi)) ≈ in.Tm_K.

inline double find_Tm_crossing(const std::vector<KirchhoffResult>& scan)
{
    for (std::size_t i = 1; i < scan.size(); ++i) {
        const double g0 = scan[i-1].delta_G;
        const double g1 = scan[i  ].delta_G;
        if (g0 * g1 <= 0.0) {
            // Linear interpolation of zero crossing: T* = T0 − g0·(T1−T0)/(g1−g0)
            const double T0 = scan[i-1].T_K;
            const double T1 = scan[i  ].T_K;
            return T0 - g0 * (T1 - T0) / (g1 - g0);
        }
    }
    return -1.0; // no crossing found
}

// ─── ΔΔG between two ligands at T ────────────────────────────────────────────
// Selectivity: ΔΔG(T) = ΔG_A(T) − ΔG_B(T)
// Negative means A binds more tightly at T.

inline double kirchhoff_selectivity(
    const KirchhoffInput& ligand_A,
    const KirchhoffInput& ligand_B,
    double T_K)
{
    return kirchhoff_deltaG(ligand_A, T_K).delta_G
         - kirchhoff_deltaG(ligand_B, T_K).delta_G;
}

// ─── Stability window ─────────────────────────────────────────────────────────
// Returns the temperature range [T_lo, T_hi] over which ΔG(T) < threshold
// (default threshold = 0.0: the window of spontaneous binding).
// Returns {-1, -1} if the condition is never met in the scan range.

struct StabilityWindow {
    double T_lo_K = -1.0;
    double T_hi_K = -1.0;
    bool   valid  = false;
};

inline StabilityWindow stability_window(
    const KirchhoffInput& in,
    double T_lo_K,
    double T_hi_K,
    double threshold_kcal_mol = 0.0,
    int    n_steps = 200)
{
    const auto scan = kirchhoff_scan(in, T_lo_K, T_hi_K, n_steps);
    StabilityWindow w;
    for (const auto& r : scan) {
        if (r.delta_G < threshold_kcal_mol) {
            if (!w.valid) { w.T_lo_K = r.T_K; w.valid = true; }
            w.T_hi_K = r.T_K;
        }
    }
    return w;
}

} // namespace thermal_extrap
