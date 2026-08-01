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

} // namespace dataset
