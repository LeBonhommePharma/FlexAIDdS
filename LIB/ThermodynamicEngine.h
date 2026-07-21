#pragma once
#include <vector>
#include <array>
#include <string>

struct ThermoResult {
    float G_bind;        // ΔG_bind = T_eff*H_vct_raw − TdS_shannon + TdS_vib  (kcal/mol); entropy costs binding (+TdS), vib gain stabilizes (−)
    float H_vct;         // <CF_vct>/n_heavy — intensive enthalpic term (kcal/mol per heavy atom)
    float H_vct_raw;     // <CF_vct> ensemble mean — unnormalized (for diagnostics)
    int   n_heavy_atoms; // heavy-atom count used for normalization
    float TdS_shannon;   // T_eff * H_shannon — configurational entropy
    float TdS_vib;       // tencom_scale * (H_rep_bound - H_rep_ref) — vibrational entropy
    float compensation;  // H_vct / (TdS_shannon + TdS_vib)

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
