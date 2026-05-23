// disco_natural.hpp — DISCO multimodal diffusion guidance for NATURaL co-translational assembly
//
// Integrates the DISCO (DIffusion for Sequence-structure CO-design) framework's
// inference-time scaling mechanics into NATURaL's co-translational simulation:
//
//   1. Entropy-Adaptive Temperature Scaling
//      Maps DISCO's per-step sampling temperature modulation onto the co-translational
//      genetic algorithm.  High local Shannon entropy (flat / allosteric landscape) →
//      raise effective temperature to prevent kinetic trapping.
//      Reference: DISCO "designable" inference preset (github.com/DISCO-design/DISCO).
//
//   2. Noisy Guidance Penalty
//      Continuous corrective force that steers the search away from non-physical
//      chain orientations by penalising rotations where |det(R) − 1| > ε.
//      Analogous to DISCO's noisy guidance coefficient α.
//
//   3. FunnelScorer
//      Evaluates each co-translational GrowthStep as either a deep enthalpic funnel
//      (orthosteric, easily detected by AI co-folding) or a neutral-frustration
//      landscape (allosteric, invisible to pattern-recognition architectures).
//      Score combines ΔH (CF score), −T_eff·ΔS (Shannon entropy correction), and
//      the noisy-guidance penalty.
//
// Chain orientation is encoded as R_flat — a flattened, row-major 3×3 rotation
// matrix (9 doubles, compatible with Eigen::Map<Eigen::Matrix3d>).
//
// Integration point in DualAssemblyEngine:
//   DiscoGuidance guide = DiscoGuidance::make_designable(cfg.temperature_K);
//   for (auto& step : steps)
//       step.cf_score += guide.score_step(step.shannon_entropy,
//                                         step.cf_score,
//                                         compute_chain_R_flat(step));
//
// References:
//   DISCO: github.com/DISCO-design/DISCO (General Multimodal Protein Design)
//   Ferreiro 2014 PNAS — local frustration analysis, neutral-frustration landscapes
//   Zhao 2011 J. Phys. Chem. B — master equation for co-translational elongation
#pragma once

#include <array>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace natural {

// ─── R_flat ───────────────────────────────────────────────────────────────────
// Flattened 3×3 rotation matrix (row-major, 9 doubles).
// Encodes the spatial orientation of the nascent chain segment or TM helix at
// each co-translational GrowthStep.
//
// Compatible with Eigen via:
//   Eigen::Map<const Eigen::Matrix<double,3,3,Eigen::RowMajor>> R(rf.data());
using R_flat = std::array<double, 9>;

// ─── FunnelParams ─────────────────────────────────────────────────────────────
struct FunnelParams {
    // Base sampling temperature (K).
    double temperature_K = 310.0;

    // Entropy-adaptive weight w: T_eff = T_base × (1 + w × H).
    // Mirrors DISCO's per-step temperature modulation in the "designable" preset.
    // Set to 0.0 to disable adaptive scaling (equivalent to DISCO "diverse" preset).
    double entropy_weight = 1.0;

    // Noisy-guidance coefficient α: scales the rotation-validity penalty.
    // Penalises orientations with |det(R) − 1| > ε (non-physical rotations).
    // Set to 0.0 to disable (DISCO "diverse" preset: guidance off).
    double noisy_guidance_alpha = 0.1;

    // Neutral-frustration fraction threshold for allosteric classification.
    // Allosteric sites: ~68–71 % neutral frustration in apo and holo states.
    // Reference: Ferreiro 2014 PNAS 111, 12186.
    double neutral_frustration_threshold = 0.68;

    // Minimum funnel depth (kcal/mol) to classify a step as deep-funnel
    // (orthosteric).  Steps with score > funnel_depth_min are allosteric.
    double funnel_depth_min = -2.0;
};

// ─── FunnelScorer ─────────────────────────────────────────────────────────────
// Evaluates a co-translational growth step with DISCO-informed thermodynamics.
//
// The composite score is:
//   score = cf_score  −  T_eff · kB · entropy  +  guidance_penalty(R_flat_)
//
// where T_eff = adaptive_temperature(entropy) and kB = 0.001987 kcal/(mol·K).
class FunnelScorer {
public:
    explicit FunnelScorer(FunnelParams params = {}) noexcept
        : params_(params) {}

    // ── Primary entry point ──────────────────────────────────────────────────
    // \param R_flat_   Flattened 3×3 row-major rotation matrix for chain
    //                  orientation at this step.  Passed FIRST so callers
    //                  cannot accidentally omit the orientation argument.
    // \param entropy   Shannon entropy of the current growth ensemble (nats).
    // \param cf_score  Contact Function score (kcal/mol).
    // \returns         Composite funnel score (kcal/mol).
    //                  More negative = deeper funnel (orthosteric character).
    [[nodiscard]] double Score(const R_flat& R_flat_,
                                double        entropy,
                                double        cf_score) const noexcept;

    // ── Entropy-adaptive temperature ─────────────────────────────────────────
    // T_eff = T_base × (1 + entropy_weight × H)
    [[nodiscard]] double adaptive_temperature(double entropy) const noexcept;

    // ── Noisy guidance penalty ────────────────────────────────────────────────
    // penalty = alpha × |det(R) − 1| × |cf_score|
    // Penalises orientations that deviate from SO(3) (det ≠ 1).
    [[nodiscard]] double guidance_penalty(const R_flat& R_flat_,
                                           double        cf_score) const noexcept;

    // ── Frustration classification ────────────────────────────────────────────
    // Returns true when the step exhibits neutral-frustration entropy consistent
    // with an allosteric (non-deep-funnel) landscape.
    [[nodiscard]] bool is_neutral_frustration(double entropy) const noexcept;

    [[nodiscard]] const FunnelParams& params() const noexcept { return params_; }

private:
    FunnelParams params_;
};

// ─── DiscoGuidance ────────────────────────────────────────────────────────────
// High-level wrapper that holds a FunnelScorer (member: funnel_) and exposes
// score_step() for DualAssemblyEngine integration.
//
// Two factory presets mirror DISCO's inference configurations:
//   make_designable() — entropy-adaptive scaling + noisy guidance active
//   make_diverse()    — both mechanisms disabled; free landscape sampling
class DiscoGuidance {
public:
    explicit DiscoGuidance(FunnelParams params = {}) noexcept
        : funnel_(params) {}

    // ── Per-step scoring ─────────────────────────────────────────────────────
    // \param entropy       Shannon entropy at this growth step (nats).
    // \param cf_score      Contact Function score (kcal/mol).
    // \param chain_R_flat  Orientation of the emerging chain segment.
    // \returns             Funnel score (kcal/mol).
    [[nodiscard]] double score_step(double        entropy,
                                     double        cf_score,
                                     const R_flat& chain_R_flat) const noexcept {
        // funnel_.Score: R_flat is FIRST argument (chain orientation).
        return funnel_.Score(chain_R_flat, entropy, cf_score);
    }

    // ── DISCO "designable" preset ─────────────────────────────────────────────
    // Entropy-adaptive temperature scaling + noisy guidance both active.
    // Prioritises structural designability over conformational diversity.
    [[nodiscard]] static DiscoGuidance make_designable(
            double temperature_K = 310.0) noexcept {
        return DiscoGuidance(FunnelParams{
            .temperature_K                = temperature_K,
            .entropy_weight               = 1.0,
            .noisy_guidance_alpha         = 0.1,
            .neutral_frustration_threshold = 0.68,
            .funnel_depth_min             = -2.0,
        });
    }

    // ── DISCO "diverse" preset ────────────────────────────────────────────────
    // Both mechanisms disabled; model samples freely from the landscape.
    // Greater structural variety at the cost of lower average designability.
    [[nodiscard]] static DiscoGuidance make_diverse(
            double temperature_K = 310.0) noexcept {
        return DiscoGuidance(FunnelParams{
            .temperature_K                = temperature_K,
            .entropy_weight               = 0.0,
            .noisy_guidance_alpha         = 0.0,
            .neutral_frustration_threshold = 0.68,
            .funnel_depth_min             = -2.0,
        });
    }

    [[nodiscard]] const FunnelScorer& funnel() const noexcept { return funnel_; }

private:
    FunnelScorer funnel_;
};

// ─── R_flat helpers ───────────────────────────────────────────────────────────

// Returns the 3×3 identity rotation (no orientation bias).
[[nodiscard]] inline R_flat identity_R_flat() noexcept {
    return {1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0};
}

// Convert an Eigen::Matrix3d to a flat row-major array.
[[nodiscard]] inline R_flat R_flat_from_eigen(const Eigen::Matrix3d& R) noexcept {
    R_flat out{};
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(out.data()) = R;
    return out;
}

// ─── FunnelScorer inline implementations ─────────────────────────────────────

inline double FunnelScorer::Score(const R_flat& R_flat_,
                                   double        entropy,
                                   double        cf_score) const noexcept {
    static constexpr double kB = 0.001987; // kcal / (mol · K)
    const double T_eff   = adaptive_temperature(entropy);
    const double penalty = guidance_penalty(R_flat_, cf_score);
    // ΔG_eff = ΔH  −  T_eff · ΔS · kB  +  noisy-guidance penalty
    return cf_score - T_eff * kB * entropy + penalty;
}

inline double FunnelScorer::adaptive_temperature(double entropy) const noexcept {
    // DISCO entropy-adaptive: raise T in high-entropy (flat / allosteric) regions
    // so the GA avoids becoming trapped in neutral-frustration wells.
    return params_.temperature_K * (1.0 + params_.entropy_weight * entropy);
}

inline double FunnelScorer::guidance_penalty(const R_flat& R_flat_,
                                              double        cf_score) const noexcept {
    if (params_.noisy_guidance_alpha == 0.0) return 0.0;
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> R(R_flat_.data());
    const double det_dev = std::abs(R.determinant() - 1.0);
    return params_.noisy_guidance_alpha * det_dev * std::abs(cf_score);
}

inline bool FunnelScorer::is_neutral_frustration(double entropy) const noexcept {
    // Allosteric sites retain persistently high entropy in both apo and holo states;
    // orthosteric sites resolve to low entropy upon ligand binding (deep funnel).
    static constexpr double kB = 0.001987; // kcal / (mol · K)
    const double threshold =
        params_.neutral_frustration_threshold * params_.temperature_K * kB;
    return entropy > threshold;
}

} // namespace natural
