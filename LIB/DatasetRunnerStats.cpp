// =============================================================================
// DatasetRunnerStats.cpp — Pure correlation / RMSD helpers for DatasetRunner
//
// Extracted from DatasetRunner.cpp (P0 leaf of the split plan). Behaviour is
// intentionally identical to the pre-split free functions.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include "DatasetRunnerStats.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace dataset {

double compute_pearson_r(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    const size_t n = x.size();

    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / static_cast<double>(n);
    double mean_y = sum_y / static_cast<double>(n);

    double cov_xy = 0.0, var_x = 0.0, var_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        cov_xy += dx * dy;
        var_x  += dx * dx;
        var_y  += dy * dy;
    }

    double denom = std::sqrt(var_x * var_y);
    if (denom < 1e-15) return 0.0;
    return cov_xy / denom;
}

/// Helper: compute ranks for a vector (average rank for ties)
static std::vector<double> compute_ranks(const std::vector<double>& vals) {
    const size_t n = vals.size();
    std::vector<size_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&vals](size_t a, size_t b) { return vals[a] < vals[b]; });

    std::vector<double> ranks(n);
    size_t i = 0;
    while (i < n) {
        size_t j = i;
        // Find all tied elements
        while (j < n && vals[indices[j]] == vals[indices[i]]) ++j;
        // Average rank for ties (1-based)
        double avg_rank = 0.5 * (static_cast<double>(i + 1) + static_cast<double>(j));
        for (size_t k = i; k < j; ++k) {
            ranks[indices[k]] = avg_rank;
        }
        i = j;
    }
    return ranks;
}

double compute_spearman_rho(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    std::vector<double> rx = compute_ranks(x);
    std::vector<double> ry = compute_ranks(y);
    return compute_pearson_r(rx, ry);
}

double compute_kendall_tau(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 2) return 0.0;
    const size_t n = x.size();

    int64_t concordant = 0, discordant = 0;
    int64_t ties_x = 0, ties_y = 0, ties_xy = 0;

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double sign_product = dx * dy;

            bool tx = (std::abs(dx) < 1e-12);
            bool ty = (std::abs(dy) < 1e-12);

            if (tx && ty) {
                ties_xy++;
            } else if (tx) {
                ties_x++;
            } else if (ty) {
                ties_y++;
            } else if (sign_product > 0) {
                concordant++;
            } else {
                discordant++;
            }
        }
    }

    int64_t n_pairs = static_cast<int64_t>(n) * (static_cast<int64_t>(n) - 1) / 2;
    double n0 = static_cast<double>(n_pairs);
    double n1 = static_cast<double>(ties_x + ties_xy);
    double n2 = static_cast<double>(ties_y + ties_xy);

    double denom = std::sqrt((n0 - n1) * (n0 - n2));
    if (denom < 1e-15) return 0.0;

    return static_cast<double>(concordant - discordant) / denom;
}

double compute_rmsd(const std::vector<float>& coords_a,
                    const std::vector<float>& coords_b) {
    if (coords_a.size() != coords_b.size() || coords_a.empty()) return -1.0;
    if (coords_a.size() % 3 != 0) return -1.0;

    const size_t n_atoms = coords_a.size() / 3;
    double sum_sq = 0.0;
    for (size_t i = 0; i < coords_a.size(); ++i) {
        double d = static_cast<double>(coords_a[i]) - static_cast<double>(coords_b[i]);
        sum_sq += d * d;
    }
    return std::sqrt(sum_sq / static_cast<double>(n_atoms));
}

} // namespace dataset
