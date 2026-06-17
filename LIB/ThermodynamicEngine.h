#pragma once
#include <vector>
#include <array>

struct ThermoResult {
    float G_bind;       // total free energy estimate
    float H_vct;        // <CF_vct> ensemble mean — enthalpic term
    float TdS_shannon;  // T_eff * H_shannon — configurational entropy
    float TdS_vib;      // T_eff * dH_rep_tencom — vibrational entropy
    float compensation; // H_vct / (TdS_shannon + TdS_vib) — should ~1.0 at calibration
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
        float H_rep_bound_complex                                // tENCoM of bound pose
    ) const;

private:
    float T_eff_;
    float tencom_scale_;
    float H_rep_ref_;  // H_rep_receptor + H_rep_ligand (set by set_unbound_reference)

    static float shannon_entropy(const std::vector<std::vector<float>>& pop);
    static float ensemble_mean(const std::vector<float>& cf_values);
};
