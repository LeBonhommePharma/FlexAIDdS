#pragma once
#include <vector>
#include <array>

struct ThermoResult {
    float G_bind;        // ΔG_bind = H_vct + TdS_shannon − TdS_vib  (kcal/mol); entropy costs binding (+), vib gain reduces G (−)
    float H_vct;         // <CF_vct>/n_heavy — intensive enthalpic term (kcal/mol per heavy atom)
    float H_vct_raw;     // <CF_vct> ensemble mean — unnormalized (for diagnostics)
    int   n_heavy_atoms; // heavy-atom count used for normalization
    float TdS_shannon;   // T_eff * H_shannon — configurational entropy
    float TdS_vib;       // tencom_scale * (H_rep_bound - H_rep_ref) — vibrational entropy
    float compensation;  // H_vct / (TdS_shannon + TdS_vib)
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
        int n_heavy_atoms = 0                                    // 0 = no normalization
    ) const;

private:
    float T_eff_;
    float tencom_scale_;
    float H_rep_ref_;  // H_rep_receptor + H_rep_ligand (set by set_unbound_reference)

    static float shannon_entropy(const std::vector<std::vector<float>>& pop);
    static float ensemble_mean(const std::vector<float>& cf_values);
};
