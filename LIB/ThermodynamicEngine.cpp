#include "ThermodynamicEngine.h"
#include <cmath>
#include <numeric>
#include <algorithm>

ThermodynamicEngine::ThermodynamicEngine(float T_eff, float tencom_scale)
    : T_eff_(T_eff), tencom_scale_(tencom_scale), H_rep_ref_(0.0f) {}

void ThermodynamicEngine::set_unbound_reference(float H_rep_receptor, float H_rep_ligand) {
    H_rep_ref_ = H_rep_receptor + H_rep_ligand;
}

float ThermodynamicEngine::shannon_entropy(const std::vector<std::vector<float>>& pop) {
    if (pop.empty()) return 0.0f;
    const int n_chrom = static_cast<int>(pop.size());
    const int n_genes = static_cast<int>(pop[0].size());

    float H_total = 0.0f;
    for (int g = 0; g < n_genes; ++g) {
        std::array<int, 256> hist{};
        for (const auto& chrom : pop) {
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

float ThermodynamicEngine::ensemble_mean(const std::vector<float>& cf) {
    if (cf.empty()) return 0.0f;
    return std::accumulate(cf.begin(), cf.end(), 0.0f) / static_cast<float>(cf.size());
}

ThermoResult ThermodynamicEngine::compute(
        const std::vector<std::vector<float>>& final_pop,
        const std::vector<float>& cf_values,
        float H_rep_bound,
        int n_heavy_atoms) const {

    ThermoResult r{};
    float n_heavy    = (n_heavy_atoms > 0) ? static_cast<float>(n_heavy_atoms) : 1.0f;
    r.H_vct_raw      = ensemble_mean(cf_values);
    r.H_vct          = r.H_vct_raw / n_heavy;   // intensive (ITC-comparable, kcal/mol per heavy atom)
    r.n_heavy_atoms  = n_heavy_atoms;
    // per-gene mean Shannon entropy × T_eff: dividing by n_genes prevents
    // large ligands (many torsion genes) from inflating TdS_shannon linearly.
    {
        const int n_genes = final_pop.empty() ? 1
                          : static_cast<int>(final_pop[0].size());
        r.TdS_shannon = T_eff_ * shannon_entropy(final_pop)
                      / static_cast<float>(n_genes > 0 ? n_genes : 1);
    }
    r.TdS_vib        = tencom_scale_ * (H_rep_bound - H_rep_ref_);
    r.G_bind         = r.H_vct + r.TdS_shannon - r.TdS_vib;  // ΔG = ΔH + TΔS_conf − TΔS_vib: entropy costs (+), vib gain reduces G (−)
    float denom      = r.TdS_shannon + r.TdS_vib;
    r.compensation   = (std::abs(denom) > 1e-6f) ? r.H_vct / denom : 0.0f;
    return r;
}
