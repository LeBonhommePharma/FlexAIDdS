// emission_guards.h — fail-closed emission / RMSD reporting guards
//
// Pure helpers for pathological campaign modes:
//   (a) empty election / sentinel RMSD (1KZK-class)
//   (b) clash-sentinel CF≈10000 coarse-init (1M2Z-class) — see sampling_coverage.h
//   (c) absurd elected RMSD ≫ pocket scale (1J3J/1IGJ-class)
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
#pragma once

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace flexaids {
namespace emission {

/// Sentinel RMSD used when pose/RMSD was not computed (DatasetRunner convention).
constexpr float kSentinelRmsd = -1.0f;

/// Default absurd-RMSD threshold (Å). Pocket-scale poses are ≪ this.
constexpr float kDefaultAbsurdRmsdA = 25.0f;

/// True if RMSD is the not-computed sentinel (< 0).
inline bool is_sentinel_rmsd(float rmsd) noexcept {
    return !std::isfinite(rmsd) || rmsd < 0.0f;
}

/// True if RMSD is finite but absurdly large vs pocket scale.
inline bool is_absurd_rmsd(float rmsd,
                           float threshold_A = kDefaultAbsurdRmsdA) noexcept {
    return std::isfinite(rmsd) && rmsd > threshold_A;
}

/// True if CF is the hard-clash floor (≈10000).
inline bool is_clash_sentinel_cf(double cf,
                                 double clash_threshold = 1e4) noexcept {
    return std::isfinite(cf) && cf >= clash_threshold;
}

/// Fail-closed elected RMSD for reporting: keep finite non-absurd values;
/// map non-finite to sentinel; flag absurd without inventing a success.
struct RmsdReport {
    float rmsd{kSentinelRmsd};
    bool  is_sentinel{true};
    bool  is_absurd{false};
    bool  usable_for_success{false};  ///< true only if finite and in [0, 2] potential band setup
};

inline RmsdReport report_elected_rmsd(float rmsd,
                                      float absurd_threshold_A = kDefaultAbsurdRmsdA) noexcept {
    RmsdReport out;
    if (is_sentinel_rmsd(rmsd)) {
        out.rmsd = kSentinelRmsd;
        out.is_sentinel = true;
        out.is_absurd = false;
        out.usable_for_success = false;
        return out;
    }
    out.rmsd = rmsd;
    out.is_sentinel = false;
    out.is_absurd = is_absurd_rmsd(rmsd, absurd_threshold_A);
    // Success gate still uses caller threshold (2 Å); absurd is never success.
    out.usable_for_success = !out.is_absurd && rmsd >= 0.0f;
    return out;
}

/// Election completeness: docking is claim-ready for pose election only when
/// we have at least one real pose path (or counted head) and not stuck.
inline bool election_inputs_ok(int n_poses,
                               bool stuck,
                               bool has_pose_path) noexcept {
    return n_poses > 0 && !stuck && has_pose_path;
}

/// Recover pose count when root dir lacks _mode_/_cluster_ names but restart
/// prefixes already enumerated emitted cluster heads (1KZK-class).
inline int recover_pose_count(int counted_n_poses, int enumerated_heads) noexcept {
    if (enumerated_heads < 0) enumerated_heads = 0;
    if (counted_n_poses < 0) counted_n_poses = 0;
    return counted_n_poses > enumerated_heads ? counted_n_poses : enumerated_heads;
}

/// Choose elected path: prefer primary election; else first non-empty fallback.
inline std::string choose_elected_path(const std::string& primary,
                                       const std::vector<std::string>& fallbacks) {
    if (!primary.empty()) return primary;
    for (const auto& p : fallbacks) {
        if (!p.empty()) return p;
    }
    return {};
}

/// Human-readable pathology tag for logs/CSV diagnostics.
inline const char* pathology_tag(const RmsdReport& r, double best_cf) noexcept {
    if (r.is_sentinel) return "empty_or_sentinel_rmsd";
    if (r.is_absurd) return "absurd_rmsd";
    if (is_clash_sentinel_cf(best_cf)) return "clash_sentinel_cf";
    return "ok";
}

}  // namespace emission
}  // namespace flexaids
