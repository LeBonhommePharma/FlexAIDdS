#pragma once
#include <vector>
#include <array>
#include <string>

// ── Thermodynamic impossibility gate (LP whiteboard IMG_3696) ───────────────
// ΔG = ΔH − TΔS. When ΔH > 0 and ΔS < 0 simultaneously, −TΔS > 0 for every
// T > 0, so ΔG is strictly positive at all temperatures: there is no
// temperature at which such a pose binds spontaneously. Such poses are given a
// large positive sentinel.
//
// NOT WIRED IN: the sentinel currently lands only in the reported
// ThermoResult::dG_eff. No clustering, ranking or selection path reads it, so
// the gate cannot presently demote a pose. See the ΔG_eff block below.
//
// Header-only so unit tests need no extra translation unit.
//
// NOTE on the ΔS source: the population Shannon entropy H = −Σ P_i·ln P_i is
// ≥ 0 by construction (0 < P_i ≤ 1 ⇒ −P_i·ln P_i ≥ 0), so feeding H here makes
// the ΔS < 0 arm unreachable and the gate dead code. The vibrational entropy
// term (TdS_vib) is the only ΔS in this engine that takes negative values —
// measured −1.86 … −1.95 on 1SG0/2GBP/1OF1 — so it is the sign source used.
// sign(T·ΔS) == sign(ΔS) for T > 0, so using the T-scaled term is valid here.
namespace thermo_gate {

constexpr float kImpossibleSentinel = 1000.0f;

/// Strict two-sided physics test: positive enthalpy AND negative entropy.
inline bool is_impossible(float dH, float dS) {
    return dH > 0.0f && dS < 0.0f;
}

/// Returns the penalised ΔG_eff when impossible, otherwise the input unchanged.
/// `flagged` (optional) receives the gate verdict.
inline float apply_gate(float dG_eff, float dH, float dS, bool* flagged = nullptr) {
    const bool bad = is_impossible(dH, dS);
    if (flagged) *flagged = bad;
    return bad ? kImpossibleSentinel : dG_eff;
}

}  // namespace thermo_gate

struct ThermoResult {
    float G_bind;        // ΔG_bind = T_eff*H_vct_raw − TdS_shannon + TdS_vib  (kcal/mol); entropy costs binding (+TdS), vib gain stabilizes (−)
    float H_vct;         // <CF_vct>/n_heavy — intensive enthalpic term (kcal/mol per heavy atom)
    float H_vct_raw;     // <CF_vct> ensemble mean — unnormalized (for diagnostics)
    int   n_heavy_atoms; // heavy-atom count used for normalization
    float TdS_shannon;   // T_eff * H_shannon — configurational entropy
    float TdS_vib;       // tencom_scale * (H_rep_bound - H_rep_ref) — vibrational entropy
    float compensation;  // H_vct / (TdS_shannon + TdS_vib)

    // ── ΔG_eff: Boltzmann pose-population free energy (LP derivation board) ──
    //   P_i     = e^(−CF_i/T) / Z          (Boltzmann weights over the pose population)
    //   <CF>    = Σ P_i·CF_i               (enthalpy proxy)
    //   H       = −Σ P_i·ln P_i            (Shannon entropy of the pose population, nats)
    //   ΔG_eff  = <CF> − T·H
    // A broad shallow population (high H) is penalized relative to a sharp deep
    // minimum (low H). Computed at BOTH calibrations because T sets the P_i
    // distribution itself, not just the entropy prefactor: T_eff (scoring
    // temperature, default 0.596) and report_T (ISMB 2017, default 21.0).
    //
// DIAGNOSTIC ONLY — dG_eff does not affect pose selection, at any flag setting.
// FLEXAIDDS_THERMO_SCORE=1 only enables the impossibility gate that rewrites
// this field; it does not promote dG_eff to a ranking criterion. The sole
// consumer is the [THERMO3] printf in gaboom.cpp, which runs after the
// QuickSort that establishes the ranking. Note also that dG_eff is a single
// ensemble-level scalar over the whole population, not a per-pose quantity, so
// using it as a ranking key would require a per-pose reformulation first —
// not merely re-ordering the existing calls.
    float dG_eff;        // <CF> − T_eff·H      at T_eff
    float mean_CF;       // <CF> = Σ P_i·CF_i   at T_eff
    float H_pose;        // H    = −Σ P_i·ln P_i at T_eff (nats)
    float T_eff_used;    // T_eff actually used (echoed for logging)
    float dG_eff_T21;    // <CF> − report_T·H   at report_T
    float mean_CF_T21;   // <CF> = Σ P_i·CF_i   at report_T
    float H_pose_T21;    // H    = −Σ P_i·ln P_i at report_T (nats)

    // ── Impossibility gate (only populated when FLEXAIDDS_THERMO_SCORE=1) ──
    bool  thermo_impossible;   // aggregate verdict; reported dG_eff forced to +1000 when true (diagnostic only — nothing ranks on it)
    int   n_impossible_poses;  // per-pose violations (ΔH_i = CF_i > 0 with ΔS < 0)
    float gate_dS_used;        // ΔS value the gate tested (TdS_vib)

    // ── Reporting-only whiteboard diagnostics (ThermoWhiteboard.h) ──
    // Computed at report_T (default kT_ISMB=21.0, ISMB 2017 calibration), NOT
    // T_eff; none of these feed back into G_bind/compensation/CF scoring or
    // GA. Per the whiteboard, T=21 is baked into the LEFT-hand quantity
    // itself (ΔG₂₁, P_i(T=21)) — report_T is that defining constant, echoed
    // here so every downstream consumer can label output "(T=21)"/"_T21".
    float report_T;      // reporting temperature actually used (default 21.0 = kT_ISMB)
    float I_ES;           // Enthalpy-Entropy Index at T=report_T: (dH+TdS)/(dH-TdS), bounded [-1,+1]
    float CF_r2s;          // CF_i minus LS/T-T3 contact sums (0/0 passed in => CF_r2s == H_vct_raw, documented no-op)
    std::string binding_regime; // "no_binding" | "enthalpy_driven" | "both_favorable" | "entropy_driven" | "borderline"
};

class ThermodynamicEngine {
public:
    explicit ThermodynamicEngine(float T_eff, float tencom_scale = 1.0f);

    // Call once after receptor load, before GA: set unbound vibrational entropy reference
    void set_unbound_reference(float H_rep_receptor_only, float H_rep_ligand_free);

    // Call after GA convergence with final population
    ThermoResult compute(
        const std::vector<std::vector<float>>& final_gene_pop,  // [n_chrom][n_genes], values in [0,1]
        const std::vector<float>& cf_values,                     // CF per chromosome
        float H_rep_bound_complex,                               // tENCoM of bound pose
        int n_heavy_atoms = 0,                                   // 0 = no normalization
        float report_T = 21.0f,                                  // ΔG₂₁/P_i(T=21) reporting T (ISMB 2017 calibration); does NOT affect G_bind
        float sum_LS = 0.0f,                                     // ligand-solution contact CF sum for CF_r2s (0 = untracked, no-op)
        float sum_TT3 = 0.0f                                     // target-target contact CF sum for CF_r2s (0 = untracked, no-op)
    ) const;

private:
    float T_eff_;
    float tencom_scale_;
    float H_rep_ref_;  // H_rep_receptor + H_rep_ligand (set by set_unbound_reference)

    static float shannon_entropy(const std::vector<std::vector<float>>& pop);
    static float ensemble_mean(const std::vector<float>& cf_values);
};
