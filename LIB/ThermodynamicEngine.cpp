#include "ThermodynamicEngine.h"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

namespace {
constexpr float kThermoVctClashGuard = 1.0e5f;
}

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
    if (cf.empty()) return std::numeric_limits<float>::quiet_NaN();

    double valid_sum  = 0.0;   // non-clash finite values
    double finite_sum = 0.0;   // all finite values
    int finite_count  = 0;
    int valid_count   = 0;
    for (float x : cf) {
        if (!std::isfinite(x)) continue;
        ++finite_count;
        finite_sum += static_cast<double>(x);
        if (std::abs(x) >= kThermoVctClashGuard) continue;
        valid_sum += static_cast<double>(x);
        ++valid_count;
    }

    if (finite_count == 0)
        return std::numeric_limits<float>::quiet_NaN();

    // Raw clash penalties are not binding enthalpy. Prefer the mean of the
    // non-clash subset. If the converged population is dominated by — or made
    // entirely of — catastrophic clash scores, fall back to the overall finite
    // mean rather than returning NaN, which would poison every downstream ΔG/ΔH
    // and silently drop the target from selection.
    if (valid_count > 0)
        return static_cast<float>(valid_sum / static_cast<double>(valid_count));
    return static_cast<float>(finite_sum / static_cast<double>(finite_count));
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
    // Guard: when H_rep_bound ≈ 0, the bound-state vibrational entropy is
    // undefined — the ligand ENM failed to converge (rigid/tiny ligand) or the
    // ligand is displaced outside the pocket (PoseX). In both cases the receptor
    // flexibility change upon binding is physically undefined; TdS_vib = 0 is
    // the only mechanistically correct choice.
    // Without this guard: TdS_vib = scale*(0 − H_rep_ref) → large negative.
    {
        float raw_vib = (std::abs(H_rep_bound) < 1e-6f)
                        ? 0.0f
                        : tencom_scale_ * (H_rep_bound - H_rep_ref_);
        // Secondary safety clamp: |TdS_vib| > 5 nats is unphysical for any
        // real docking scenario. Catches residual explosion when a displaced
        // ligand ENM does build but H_rep_bound << H_rep_ref_.
        r.TdS_vib = std::max(-5.0f, std::min(5.0f, raw_vib));
    }
    // ── G_bind: v100-regression fix (restores v88 effective signal weighting) ──
    //   v88 (working, 91.7% BCD):   G = T_eff * H_vct_raw  − TdS_shannon  + TdS_vib
    //   v100 (broken,  9.4% BCD):   G = H_vct_raw/n_heavy  + TdS_shannon  − TdS_vib
    // Three v100 changes annihilated the native-vs-decoy gap: (1) /n_heavy
    // attenuated the VCT signal ~12×, (2) the TdS_shannon sign flipped (− → +),
    // (3) the TdS_vib sign flipped (+ → −). We restore T_eff weighting on the
    // RAW (extensive) VCT enthalpy and both correct signs, while KEEPING the two
    // v100 structural improvements that are not regressions: per-gene Shannon
    // normalization (baked into r.TdS_shannon) and the ±5-nat TdS_vib clamp.
    // r.H_vct (intensive, per-heavy-atom) is retained as an ITC-comparable
    // diagnostic only; it no longer enters G_bind.
    r.G_bind         = T_eff_ * r.H_vct_raw - r.TdS_shannon + r.TdS_vib;  // ΔG = T·ΔH_vct − TΔS_conf + TΔS_vib: entropy costs (+TdS), vib gain stabilizes (−)
    float denom      = r.TdS_shannon + r.TdS_vib;
    r.compensation   = (std::abs(denom) > 1e-6f) ? r.H_vct / denom : 0.0f;
    return r;
}
