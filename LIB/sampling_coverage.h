// sampling_coverage.h — flag-gated search/coarse-init coverage helpers
//
// Pure helpers for raising near-native head coverage (BCR class failures).
// Wired into LIB/coarse_init.cpp. Unit tests call these entry points directly.
//
// Env (default OFF = baseline bit-behavior):
//   FLEXAIDDS_SAMPLE_COVERAGE_BOOST=1
//     Multiplies coarse-init n_orient and n_seeds (default ×2 / ×2, capped).
//   FLEXAIDDS_COARSE_INIT_FORCE_RANKED=1
//     When every placement is clash-sentinel (CF ≥ CLASH_THRESHOLD), still inject
//     the best CF-ranked seeds so GA gen-0 is not empty (1M2Z-class pathology).
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace flexaids {
namespace sampling {

/// Default clash sentinel (matches flexaid.h CLASH_THRESHOLD).
constexpr double kDefaultClashThreshold = 1e4;

/// True when env var is set to a non-zero integer.
inline bool env_flag_enabled(const char* name) noexcept {
    const char* e = std::getenv(name);
    if (!e || !*e) return false;
    return std::atoi(e) != 0;
}

/// Positive integer env override; returns fallback when unset/invalid.
inline int env_int_or(const char* name, int fallback) noexcept {
    const char* e = std::getenv(name);
    if (!e || !*e) return fallback;
    const int v = std::atoi(e);
    return (v > 0) ? v : fallback;
}

struct CoverageParams {
    int n_orient{64};
    int n_seeds{25};
    bool boost_applied{false};
};

/// Apply FLEXAIDDS_SAMPLE_COVERAGE_BOOST (default OFF).
/// When enabled: n_orient *= orient_mult (default 2, max 256),
///               n_seeds  *= seeds_mult  (default 2, max 100).
inline CoverageParams apply_coverage_boost(int n_orient, int n_seeds) noexcept {
    CoverageParams out;
    out.n_orient = std::max(1, n_orient);
    out.n_seeds  = std::max(1, n_seeds);
    if (!env_flag_enabled("FLEXAIDDS_SAMPLE_COVERAGE_BOOST")) {
        return out;
    }
    const int om = env_int_or("FLEXAIDDS_SAMPLE_ORIENT_MULT", 2);
    const int sm = env_int_or("FLEXAIDDS_SAMPLE_SEEDS_MULT", 2);
    out.n_orient = std::min(256, out.n_orient * std::max(1, om));
    out.n_seeds  = std::min(100, out.n_seeds  * std::max(1, sm));
    out.boost_applied = true;
    return out;
}

/// Indices of ascending CF order (stable for ties via original index).
inline std::vector<std::size_t>
rank_indices_by_cf_asc(const std::vector<double>& cf_vals) {
    std::vector<std::size_t> idx(cf_vals.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        const double ca = cf_vals[a];
        const double cb = cf_vals[b];
        const bool fa = std::isfinite(ca);
        const bool fb = std::isfinite(cb);
        if (fa != fb) return fa;  // finite before non-finite
        if (!fa) return a < b;
        if (ca != cb) return ca < cb;
        return a < b;
    });
    return idx;
}

/// Select up to n_keep ranked seed indices, dropping clash-sentinel CF
/// (cf >= clash_threshold) unless force_ranked is true.
/// Returns indices into cf_vals in ascending CF order.
inline std::vector<std::size_t>
select_ranked_seed_indices(const std::vector<double>& cf_vals,
                           int n_keep,
                           double clash_threshold = kDefaultClashThreshold,
                           bool force_ranked = false) {
    std::vector<std::size_t> out;
    if (cf_vals.empty() || n_keep <= 0) return out;
    const auto ranked = rank_indices_by_cf_asc(cf_vals);
    out.reserve(static_cast<std::size_t>(n_keep));
    for (std::size_t i : ranked) {
        if (static_cast<int>(out.size()) >= n_keep) break;
        const double c = cf_vals[i];
        if (!std::isfinite(c)) continue;
        if (!force_ranked && c >= clash_threshold) {
            // Sorted ascending: remaining are clash-scale or worse.
            break;
        }
        out.push_back(i);
    }
    // 1M2Z-class: every placement is clash-sentinel. Optionally keep the best
    // ranked finite CFs so gen-0 is not empty (search still has a gradient).
    if (out.empty() && force_ranked) {
        for (std::size_t i : ranked) {
            if (static_cast<int>(out.size()) >= n_keep) break;
            if (std::isfinite(cf_vals[i])) out.push_back(i);
        }
    }
    return out;
}

/// Whether force-ranked clash fallback is enabled (env default OFF).
inline bool force_ranked_seeds_enabled() noexcept {
    return env_flag_enabled("FLEXAIDDS_COARSE_INIT_FORCE_RANKED");
}

/// Spatial diversity filter: from ranked indices, greedily keep placements
/// whose 3D positions are at least min_dist_A apart (Å). Empty coords → no-op.
/// coords[i] = {x,y,z} for cf_vals[i]; missing/short coords skip diversity.
inline std::vector<std::size_t>
diversify_by_min_distance(const std::vector<std::size_t>& ranked,
                          const std::vector<std::array<double, 3>>* coords,
                          int n_keep,
                          double min_dist_A) {
    if (!coords || coords->empty() || n_keep <= 0 || min_dist_A <= 0.0) {
        std::vector<std::size_t> out;
        for (std::size_t i : ranked) {
            if (static_cast<int>(out.size()) >= n_keep) break;
            out.push_back(i);
        }
        return out;
    }
    const double min2 = min_dist_A * min_dist_A;
    std::vector<std::size_t> out;
    out.reserve(static_cast<std::size_t>(n_keep));
    for (std::size_t i : ranked) {
        if (static_cast<int>(out.size()) >= n_keep) break;
        if (i >= coords->size()) {
            out.push_back(i);
            continue;
        }
        const auto& p = (*coords)[i];
        bool ok = true;
        for (std::size_t j : out) {
            if (j >= coords->size()) continue;
            const auto& q = (*coords)[j];
            const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
            if (dx * dx + dy * dy + dz * dz < min2) {
                ok = false;
                break;
            }
        }
        if (ok) out.push_back(i);
    }
    // If diversity was too strict, fill remaining from ranked order.
    if (static_cast<int>(out.size()) < n_keep) {
        for (std::size_t i : ranked) {
            if (static_cast<int>(out.size()) >= n_keep) break;
            if (std::find(out.begin(), out.end(), i) == out.end())
                out.push_back(i);
        }
    }
    return out;
}

}  // namespace sampling
}  // namespace flexaids
