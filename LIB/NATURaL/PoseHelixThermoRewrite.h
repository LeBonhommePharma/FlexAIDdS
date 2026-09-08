// PoseHelixThermoRewrite.h — experimental 2014 NATURAL pose→helix ΔH/ΔS glue
//
// 2014 seminar NATURAL = "Native Assembly of Transcriptionally-Unified RNA And
// Ligand" (LP Morency, master’s seminar R2, 2014-07-25): Perl implementation of
// Zhao, Zhang, Chen, J. Chem. Phys. 135, 245101 (2011) doi:10.1063/1.3671644
// (RNA cotranscriptional folding kinetics) with (1) clustering, (2) population
// inheritance across elongation, (3) slower transcription at nucleotide repeats,
// and (4) docking poses that rewrite helix ΔH/ΔS.
//
// This module is piece (4) only. It is NOT Zhao et al. J. Phys. Chem. B 115,
// 3987 (2011) (ribosome master equation / protein cotranslation — see
// RibosomeElongation.h).
//
// Physics (seminar):
//   • Decision-helix patches only (configurable HelixSegment list).
//   • Ligand H-bond / base-pair-like contacts to loop nts → ΔS < 0
//     (loop entropy reduction).
//   • Aromatic stacking on the helix stem → more negative ΔH (stabilizing).
//   • Mixture over top-k poses with Boltzmann weights from docking scores.
//   • ΔG = ΔH − T·ΔS with ΔH in kcal, ΔS in cal/K (ΔS/1000 before T·ΔS).
//
// Experimental: never elects poses, never mutates CF, never writes
// FA->natural_deltaG. DualAssemblyEngine::run() does not apply these patches.
// See docs/POSE_HELIX_THERMO_REWRITE.md.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#pragma once

#include <string>
#include <vector>

namespace natural {

// Inclusive 0-based nucleotide ranges for one decision helix (hairpin).
// stem_i / stem_j are the two sides of the stem; loop is the intervening loop.
struct HelixSegment {
    int stem_i_begin = 0;
    int stem_i_end   = -1;
    int stem_j_begin = 0;
    int stem_j_end   = -1;
    int loop_begin   = 0;
    int loop_end     = -1;
    std::string label;
};

// One RNA nucleotide used only for contact geometry (no FlexAID atom_struct).
struct ReceptorNtCoord {
    int    nt_index = -1;
    double x = 0.0, y = 0.0, z = 0.0;
    bool   is_hbond_partner = false; // N/O donor/acceptor
    bool   is_aromatic      = false; // nucleobase ring
    // Ring normal. A zero vector skips the stacking-angle test (distance only).
    double nx = 0.0, ny = 0.0, nz = 0.0;
};

struct LigandAtomCoord {
    double x = 0.0, y = 0.0, z = 0.0;
    bool   is_hbond_partner = false;
    bool   is_aromatic      = false;
    double nx = 0.0, ny = 0.0, nz = 0.0;
};

// One docking pose. `score_kcal` is the CF/contact-function scoring proxy used
// only to form Boltzmann mixture weights — it is not a free-energy claim.
struct LigandPose {
    int    rank = 0;           // 0 = best
    double score_kcal = 0.0;  // CF proxy (more negative → higher weight)
    std::vector<LigandAtomCoord> atoms;
};

struct HelixThermoPatch {
    std::string helix_label;
    double delta_H_kcal = 0.0;
    double delta_S_cal_per_K = 0.0;
    double delta_G_kcal = 0.0;
    std::vector<int> touched_nt;
    int    pose_rank = -1;
    double pose_boltzmann_w = 0.0;
};

struct PoseHelixRewriteConfig {
    double T_K = 298.15;
    double hbond_max_A = 3.5;
    double stack_max_A = 4.5;
    double stack_angle_max_deg = 40.0;
    double dH_per_stack_kcal = -1.2;
    double dS_per_hbond_cal_K = -3.0;
    int    top_k_poses = 5;
    bool   require_decision_helix_only = true;
};

struct PoseHelixRewriteResult {
    std::vector<HelixThermoPatch> patches;   // per pose × helix
    std::vector<HelixThermoPatch> ensemble; // Boltzmann mix, one per helix
    bool applied = true; // false when DualAssembly flag is off
};

// kcal mol⁻¹ K⁻¹, same truncated value as LIB/statmech.h kB_kcal.
inline constexpr double kB_pose_helix_kcal = 0.001987206;

// ΔG [kcal] = ΔH[kcal] − T[K] · (ΔS[cal/K] / 1000).
double helix_thermo_delta_g_kcal(double delta_H_kcal,
                                 double delta_S_cal_per_K,
                                 double T_K);

// Log-sum-exp Boltzmann weights from CF scores. Empty input → empty vector.
std::vector<double> boltzmann_weights_from_scores(const std::vector<double>& scores_kcal,
                                                   double T_K);

// Pure function: docking poses → helix ΔH/ΔS patches. No GA, no CF election.
PoseHelixRewriteResult rewrite_helices_from_poses(
    const std::vector<HelixSegment>& helices,
    const std::vector<ReceptorNtCoord>& rna_nts,
    const std::vector<LigandPose>& poses,
    const PoseHelixRewriteConfig& cfg = {});

} // namespace natural
