#pragma once
// Reporting-only thermodynamic formulas (LP whiteboard derivations, ISMB 2017
// calibration). These are ADDITIVE diagnostics computed alongside the existing
// ThermodynamicEngine output — none of them feed back into CF scoring, G_bind,
// or GA selection. They use their own reporting temperature (default T=21 in
// FlexAID internal CF units, the value calibrated for the ISMB 2017 results),
// which is intentionally decoupled from ThermodynamicEngine's T_eff (0.596,
// used by G_bind/TdS_shannon and left untouched here).
//
// T_ISMB naming: per the whiteboard, T=21 is not merely substituted on the
// RHS of these equations — it's baked into the LEFT-hand quantity itself
// (ΔG₂₁, P_i(T=21), ...). kT_ISMB is that calibration constant; every output
// field/column/log tag downstream is labelled with the "_T21"/"(T=21)"
// suffix to keep that distinction visible even when --temperature overrides
// the reporting T away from 21.
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <array>

namespace thermo_whiteboard {

// ISMB 2017 calibrated reporting temperature (FlexAID internal CF units).
// This is the LEFT-hand-side-defining constant (ΔG₂₁, P_i(T=21), ...), not
// just a RHS substitution — see file header. Override via config/
// --temperature / FLEXAIDDS_REPORT_T; never used by CF scoring or the GA.
constexpr float kT_ISMB = 21.0f;
// Back-compat alias.
constexpr float kDefaultReportT = kT_ISMB;

// ── MC_st0r5.2_6 (ISMB 2017 gold-standard matrix) active atom types ────
// Of the 40 pairwise-contact types in MC_st0r5.2_6.txt (820 upper-triangular
// entries), only these carry non-zero contact energies; the rest (8, 16, 17,
// 20, 21, 27-34, 36-39) are all-zero dummy/unused slots. CF_r2s's LS/T,T3
// contact-sum decomposition (Eq. 3) only needs to classify contacts between
// types in this set — a contact where either atom's type is NOT in
// kActiveTypesMC_st0r5 contributes 0 to CF_i already, so it can be skipped
// when/if the per-contact LS/T,T3 accumulators are added to vcfunction.cpp.
// Type 40 has non-zero contacts with most other active types (universal/
// aromatic-boundary type).
constexpr std::array<int, 21> kActiveTypesMC_st0r5 = {
    1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15,
    18, 19, 22, 23, 24, 25, 26 /* , 28, 31, 35, 40 continue below */
};
constexpr std::array<int, 4> kActiveTypesMC_st0r5_ext = { 28, 31, 35, 40 };

inline bool is_active_contact_type(int t) {
    for (int a : kActiveTypesMC_st0r5) if (a == t) return true;
    for (int a : kActiveTypesMC_st0r5_ext) if (a == t) return true;
    return false;
}

// ── Eq. 7: Boltzmann pose probabilities ─────────────────────────────────
// P_i ∝ exp(-E_i / T). Works for both the complex-state form (E_i = ΔG_i)
// and the solution-state form (E_i = E_Δα,i − E_Δ0,i) — caller supplies the
// appropriate energy vector.
inline std::vector<float> boltzmann_probabilities(const std::vector<float>& energies, float T) {
    std::vector<float> P(energies.size(), 0.0f);
    if (energies.empty() || T == 0.0f) return P;

    float e_min = *std::min_element(energies.begin(), energies.end());
    double denom = 0.0;
    std::vector<double> w(energies.size());
    for (size_t i = 0; i < energies.size(); ++i) {
        w[i] = std::exp(-(static_cast<double>(energies[i]) - e_min) / T);
        denom += w[i];
    }
    if (denom <= 0.0) return P;
    for (size_t i = 0; i < energies.size(); ++i)
        P[i] = static_cast<float>(w[i] / denom);
    return P;
}

// Shannon entropy of a probability distribution: ΔS = −Σ P_i ln P_i (nats)
inline float shannon_entropy_nats(const std::vector<float>& P) {
    double H = 0.0;
    for (float p : P) {
        if (p > 0.0f) H -= static_cast<double>(p) * std::log(static_cast<double>(p));
    }
    return static_cast<float>(H);
}

// Probability-weighted mean: ΔH = Σ P_i · x_i
inline float weighted_mean(const std::vector<float>& P, const std::vector<float>& x) {
    double sum = 0.0;
    size_t n = std::min(P.size(), x.size());
    for (size_t i = 0; i < n; ++i) sum += static_cast<double>(P[i]) * static_cast<double>(x[i]);
    return static_cast<float>(sum);
}

// ── Eq. 1: Enthalpy-Entropy Index (Williams et al. 2017, DDT) ──────────
// I_ES = (ΔH + TΔS) / (ΔH − TΔS), bounded [-1, +1] when ΔH, TΔS finite and
// not both zero. Guarded against ΔH == TΔS (division by zero).
inline float enthalpy_entropy_index(float dH, float TdS_signed /* already T·ΔS */) {
    float denom = dH - TdS_signed;
    if (std::abs(denom) < 1e-6f) return 0.0f;
    return (dH + TdS_signed) / denom;
}

// ── Eq. 2: per-conformation entropy ΔS_j (proto-formulation) ───────────
// ΔS_j = (s*_tot − s*_not)·P_j − (H_holo − S_apo)
// s_tot, s_not: total / "not-bound" reference microstate entropies (caller-
// supplied scalars from the pose's contact decomposition); P_j: this pose's
// Boltzmann weight; H_holo: Shannon entropy of the holo (bound) pose
// ensemble; S_apo: Shannon entropy of the apo (unbound ligand) ensemble —
// requires a separate apo conformer sample (see delta_sj_apo_reference()).
inline float delta_Sj(float s_tot, float s_not, float P_j, float H_holo, float S_apo) {
    return (s_tot - s_not) * P_j - (H_holo - S_apo);
}

// Convenience: Shannon entropy of an apo (unbound ligand) gene/conformer
// population, using the same per-gene histogram method as
// ThermodynamicEngine::shannon_entropy, so ΔS_j's two entropy terms are
// computed identically. n_bins matches ThermodynamicEngine (256).
inline float apo_shannon_entropy(const std::vector<std::vector<float>>& apo_pop) {
    if (apo_pop.empty()) return 0.0f;
    const int n_chrom = static_cast<int>(apo_pop.size());
    const int n_genes  = static_cast<int>(apo_pop[0].size());
    float H_total = 0.0f;
    for (int g = 0; g < n_genes; ++g) {
        int hist[256] = {0};
        for (const auto& chrom : apo_pop) {
            int bin = std::clamp(static_cast<int>(chrom[g] * 255.0f), 0, 255);
            ++hist[bin];
        }
        for (int b = 0; b < 256; ++b) {
            if (hist[b] > 0) {
                float p = static_cast<float>(hist[b]) / n_chrom;
                H_total -= p * std::log2(p);
            }
        }
    }
    return H_total;
}

// ── Eq. 3: CF_r2s binding-specific enthalpy (solvation decomposition) ──
// CF_r2s = CF_i − Σ CF_cont,LS − Σ CF_cont,T,T3
// sum_LS: summed ligand-solution contact contributions; sum_TT3: summed
// target-target contact contributions. NOTE: neither sum is currently
// tracked by vcfunction.cpp's per-contact loop (contacts are only tagged
// intramolecular vs not, not classified as LS / T,T3) — callers without
// that instrumentation should pass 0.0f for both, which reduces this to
// CF_r2s == CF_i (documented no-op) until the contact loop is extended.
inline float cf_r2s(float cf_i, float sum_LS, float sum_TT3) {
    return cf_i - sum_LS - sum_TT3;
}

// ── Eq. 4: ligand solvation thermodynamics ──────────────────────────────
// ΔHs = Σ_i P_i · N_i · r2s_i   (N_i: contact count for pose i, r2s_i: that
// pose's CF_r2s per-contact enthalpy)
inline float delta_Hs(const std::vector<float>& P, const std::vector<float>& N,
                       const std::vector<float>& r2s) {
    double sum = 0.0;
    size_t n = std::min({P.size(), N.size(), r2s.size()});
    for (size_t i = 0; i < n; ++i)
        sum += static_cast<double>(P[i]) * static_cast<double>(N[i]) * static_cast<double>(r2s[i]);
    return static_cast<float>(sum);
}

// ΔSs = −Σ_i P_i,s · ln P_i,s  (entropy over the solution-state pose distribution)
inline float delta_Ss(const std::vector<float>& P_solution) {
    return shannon_entropy_nats(P_solution);
}

// ΔΔHc = ΔHc − ΔHp − ΔHL  (complex minus isolated protein minus isolated ligand)
inline float delta_delta_Hc(float dHc, float dHp, float dHL) {
    return dHc - dHp - dHL;
}

// ── Eq. 5: temperature-dependent ΔG via Gibbs-Helmholtz with heat capacity ──
// ΔG(T) = ΔHm·(1 − T/Tm) − ΔCp·[(Tm − T) + T·ln(T/Tm)]
// Tm, T must be > 0 (absolute/reference temperatures in the same units as
// the calibration ΔHm/ΔCp were derived in). Returns NaN on invalid Tm/T.
inline float delta_G_at_temperature(float T, float Tm, float dHm, float dCp) {
    if (T <= 0.0f || Tm <= 0.0f) return std::numeric_limits<float>::quiet_NaN();
    return dHm * (1.0f - T / Tm)
         - dCp * ((Tm - T) + T * std::log(T / Tm));
}

// ── Eq. 6: binding regime classifier ────────────────────────────────────
inline std::string classify_binding_regime(float dH, float dS, float T) {
    const float TdS = T * dS;
    if (dH > 0.0f && dS < 0.0f) return "no_binding";
    if (dH < 0.0f && dS < 0.0f) return (std::abs(dH) > std::abs(TdS)) ? "enthalpy_driven" : "no_binding";
    if (dH < 0.0f && dS > 0.0f) return "both_favorable";
    if (dH > 0.0f && dS > 0.0f) return "entropy_driven";
    return "borderline"; // dH == 0 or dS == 0
}

} // namespace thermo_whiteboard
