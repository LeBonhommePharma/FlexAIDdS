// tencom_diff.h — Eigenvalue/eigenvector differential engine for tENCoM
//
// Computes vibrational entropy differentials between a reference and target
// TorsionalENM, including per-mode eigenvalue differences, mode overlaps,
// and aggregate ΔS_vib / ΔF_vib.
#pragma once

#include "tencm.h"
#include "../encom.h"

#include <string>
#include <vector>
#include <cmath>
#include <limits>

namespace tencom_diff {

// Per-mode comparison between reference and target
struct ModeComparison {
    int    mode_index;        // 0-based among kept modes
    double eigenvalue_ref;
    double eigenvalue_tgt;
    double delta_eigenvalue;  // target - reference
    double overlap;           // |v_ref · v_tgt| (NaN if dimension mismatch)
};

// Full differential result between reference and one target
struct DifferentialResult {
    std::string ref_name;
    std::string tgt_name;

    // Per-mode differentials
    std::vector<ModeComparison> mode_comparisons;

    // Aggregate vibrational entropy
    encom::VibrationalEntropy svib_ref;
    encom::VibrationalEntropy svib_tgt;
    double delta_S_vib;   // S_vib(target) - S_vib(ref), kcal/mol/K
    double delta_F_vib;   // -T × delta_S_vib, kcal/mol

    // Per-residue B-factor differentials
    std::vector<float> bfactors_ref;
    std::vector<float> bfactors_tgt;
    std::vector<float> delta_bfactors;

    // Per-residue vibrational entropy decomposition
    // S_vib_i = Σ_k (v_ki² × S_mode_k) where v_ki is residue i's component
    std::vector<double> per_residue_svib_ref;
    std::vector<double> per_residue_svib_tgt;
    std::vector<double> per_residue_delta_svib;

    double temperature;
    encom::FrequencyCalibration frequency_calibration;
};

// Convert tencm::NormalMode vector to encom::NormalMode vector
// for use with ENCoMEngine::compute_vibrational_entropy()
std::vector<encom::NormalMode> to_encom_modes(
    const std::vector<tencm::NormalMode>& tencm_modes);

// Compute full differential between reference and target ENMs.
// Both must be built (is_built() == true).
DifferentialResult compute_differential(
    const tencm::TorsionalENM& ref_enm,
    const tencm::TorsionalENM& tgt_enm,
    const std::string& ref_name = "reference",
    const std::string& tgt_name = "target",
    double temperature_K = 300.0,
    double eigenvalue_cutoff = 1e-6);

DifferentialResult compute_differential(
    const tencm::TorsionalENM& ref_enm,
    const tencm::TorsionalENM& tgt_enm,
    const std::string& ref_name,
    const std::string& tgt_name,
    double temperature_K,
    const encom::FrequencyCalibration& calibration,
    double eigenvalue_cutoff = 1e-6);

}  // namespace tencom_diff
