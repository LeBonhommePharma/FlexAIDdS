// =============================================================================
// DatasetRunnerStats.h — Pure correlation / RMSD helpers for DatasetRunner
//
// Leaf extraction (P0 of the DatasetRunner split plan): no I/O, no process
// control, no DatasetRunner state. Safe to unit-test in isolation.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace dataset {

/// Pearson correlation coefficient (computed from scratch).
/// Returns 0.0 if sizes differ, n < 2, or either series is constant.
double compute_pearson_r(const std::vector<double>& x, const std::vector<double>& y);

/// Spearman rank correlation (Pearson r of average ranks; ties supported).
double compute_spearman_rho(const std::vector<double>& x, const std::vector<double>& y);

/// Kendall tau-b rank correlation (handles ties).
double compute_kendall_tau(const std::vector<double>& x, const std::vector<double>& y);

/// RMSD between two coordinate sets (3N floats: x1,y1,z1,x2,y2,z2,...).
/// Returns -1.0 on size mismatch, empty input, or non-multiple-of-3 length.
double compute_rmsd(const std::vector<float>& coords_a,
                    const std::vector<float>& coords_b);

/// Symmetry-corrected (Hungarian) heavy-atom RMSD between two element-labelled
/// coordinate sets. Atoms are matched by an optimal assignment computed
/// SEPARATELY WITHIN EACH ELEMENT TYPE, so an equivalent relabelling of
/// symmetric atoms is not penalised while a cross-element pairing (e.g. Cl<->C)
/// is impossible. Mirrors FlexAID's calc_Hungarian_RMSD. Returns -1.0f on empty
/// input. This is the single canonical implementation used by the DatasetRunner
/// production RMSD paths (compute_pose_ligand_rmsd / pose_pose_rmsd); it is
/// defined in DatasetRunner.cpp and declared here so it can be unit-tested
/// directly against the real metric.
float hungarian_rmsd(
    const std::vector<std::pair<std::string, std::array<float, 3>>>& crystal,
    const std::vector<std::pair<std::string, std::array<float, 3>>>& docked);

/// Outcome of compute_pose_ligand_rmsd: both RMSD estimands plus a
/// machine-readable failure reason. fail_reason == "none" iff serial >= 0.
/// Values (bug 2026-08-22, arms 8/9/10 wrote 0% summaries over valid poses):
///   none               — serial RMSD computed
///   ref_empty          — crystal reference missing/unparseable at runtime
///                        (the previously SILENT wholesale-failure path)
///   pose_block_empty   — CONECT fingerprint selected no heavy atoms
///   count_mismatch     — pose/crystal heavy-atom counts differ (fail-closed)
///   elem_mismatch      — element-vector length mismatch (fail-closed)
///   elem_order_mismatch— ordered RMSD undefined (element sequences differ);
///                        hungarian may still be valid
struct PoseRmsdOutcome {
    float serial{-1.0f};
    float hungarian{-1.0f};
    std::string fail_reason{"none"};
};

/// Production rank-0 / ceiling RMSD path (see DatasetRunner.cpp). Declared
/// here so unit tests exercise the exact production reason-code contract,
/// mirroring hungarian_rmsd above.
PoseRmsdOutcome compute_pose_ligand_rmsd(
    const std::string& pose_pdb,
    const std::vector<std::array<float, 3>>& crystal_xyz,
    const std::vector<std::string>& crystal_elem,
    const std::string& pdb_id,
    bool warn);

/// ── Zero-success plausibility gate (bug 2026-08-22) ─────────────────────
/// A 0% RMSD summary over a run that produced poses is more plausibly a
/// crystal-reference/measurement failure than a true all-fail: genuine
/// docking misses yield VALID RMSD values above threshold, not -1 sentinels.
/// Campaign arms 8/9/10 certified 0/85 this way while their offline pooled
/// ceilings were 37-42/85.
struct ZeroSuccessGateInput {
    int total_systems{0};
    int successful_rmsd{0};          ///< rows with success_rmsd == true
    int rows_with_any_poses{0};      ///< rows with num_poses > 0 (GA produced modes)
    int rmsd_negative_rows{0};       ///< rows with rmsd_to_crystal < 0
    /// Negative rows whose reason is wholesale/measurement-side:
    /// ref_empty | input_missing | pose_block_empty.
    int rmsd_negative_wholesale{0};
};

/// True when a 0%-success summary must be flagged suspect_zero_success:
/// enough targets for the rate to mean something, the GA produced modes on at
/// least half the rows, and negative RMSDs are dominated by wholesale reasons.
inline bool zero_success_is_suspect(const ZeroSuccessGateInput& in) {
    if (in.total_systems < 8) return false;
    if (in.successful_rmsd != 0) return false;
    if (in.rows_with_any_poses * 2 < in.total_systems) return false;
    if (in.rmsd_negative_rows <= 0) return false;
    return in.rmsd_negative_wholesale * 2 >= in.rmsd_negative_rows;
}

} // namespace dataset
