// =============================================================================
// DatasetRunnerStats.h — Pure correlation / RMSD helpers for DatasetRunner
//
// Leaf extraction (P0 of the DatasetRunner split plan): no I/O, no process
// control, no DatasetRunner state. Safe to unit-test in isolation.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#pragma once

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

} // namespace dataset
