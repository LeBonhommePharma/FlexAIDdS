// SearchCoverage.h — flag-gated GA search-coverage knobs (S1 remediation)
//
// FLEXAIDDS_SEARCH_COVERAGE=1 (default OFF) densifies BOOM injection and
// widens niche radius so near-native basins are less likely to starve under
// over-burial / false-minima dominance. Targets bulk BCR>2 failures
// (sampling_fail + sampling_near). Does not change product defaults.
//
// Critical: campaign dock_configs often set boom_inject_fraction=0 (BOOM
// fully disabled). Halving boom_inject_interval alone is a no-op when
// fraction==0 (gaboom requires interval>0 AND fraction>0). When coverage is
// enabled, fraction is restored to a positive default (1.0) if currently ≤0.
//
// Optional absolute override: FLEXAIDDS_BOOM_INTERVAL=<N> (gens, >0).
//
// Pure helpers are unit-testable without linking the GA.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace flexaids {
namespace search_coverage {

/// Default boom_inject_fraction restored when coverage enables a dead BOOM path.
constexpr double kDefaultBoomFraction = 1.0;

/// True when FLEXAIDDS_SEARCH_COVERAGE is a non-zero integer.
inline bool enabled_from_env() noexcept
{
	const char* e = std::getenv("FLEXAIDDS_SEARCH_COVERAGE");
	return e != nullptr && std::atoi(e) != 0;
}

/// Absolute boom inject interval override when FLEXAIDDS_BOOM_INTERVAL > 0.
/// Returns 0 when unset / invalid (no override).
inline int boom_interval_override_from_env() noexcept
{
	const char* e = std::getenv("FLEXAIDDS_BOOM_INTERVAL");
	if (e == nullptr)
		return 0;
	const int v = std::atoi(e);
	return v > 0 ? v : 0;
}

/// Chromosomes replaced per BOOM event (mirrors gaboom: fraction × num_chrom/2).
/// Returns 0 when BOOM is disabled (interval≤0 or fraction≤0).
inline int inject_count_per_event(int boom_interval,
                                  double boom_fraction,
                                  int num_chrom) noexcept
{
	if (boom_interval <= 0 || !(boom_fraction > 0.0) || !std::isfinite(boom_fraction))
		return 0;
	if (num_chrom <= 0)
		return 0;
	const int half = num_chrom / 2;
	const int n = static_cast<int>(boom_fraction * static_cast<double>(half));
	return n > 0 ? n : 0;
}

/// Approximate total injects over a run (events × per-event count).
/// Events fire at gens interval, 2*interval, ... while < max_generations.
inline long long total_injects_estimate(int boom_interval,
                                        double boom_fraction,
                                        int num_chrom,
                                        int max_generations) noexcept
{
	const int per = inject_count_per_event(boom_interval, boom_fraction, num_chrom);
	if (per <= 0 || max_generations <= 0 || boom_interval <= 0)
		return 0;
	// gens i+1 where (i+1)%interval==0 and (i+1)<max_generations
	long long events = 0;
	for (int g = boom_interval; g < max_generations; g += boom_interval)
		++events;
	return events * static_cast<long long>(per);
}

/// Apply S1 coverage knobs in-place.
/// - boom_interval: halved, floor 25 (when original > 0 and enabled)
/// - sharing_scale: ×0.5 when enabled (sig_share ∝ 1/scale → wider niches)
/// - boom_fraction: if ≤0 or non-finite, set to kDefaultBoomFraction so BOOM
///   actually fires (campaign configs often leave fraction=0)
/// No-op when enabled is false.
inline void apply(int& boom_interval,
                  double& sharing_scale,
                  double& boom_fraction,
                  bool enabled) noexcept
{
	if (!enabled)
		return;
	if (boom_interval > 0) {
		const int halved = boom_interval / 2;
		boom_interval = std::max(25, halved);
	}
	if (std::isfinite(sharing_scale) && sharing_scale > 0.0) {
		sharing_scale *= 0.5;
	}
	// Restore dead BOOM path: fraction==0 is the campaign inventory failure mode.
	if (!std::isfinite(boom_fraction) || !(boom_fraction > 0.0)) {
		boom_fraction = kDefaultBoomFraction;
	}
}

/// Backward-compatible overload (fraction not adjusted — prefer 3-arg form).
inline void apply(int& boom_interval, double& sharing_scale, bool enabled) noexcept
{
	double frac = kDefaultBoomFraction;  // local dummy; does not restore caller fraction
	apply(boom_interval, sharing_scale, frac, enabled);
	(void)frac;
}

/// Compose env-driven coverage + optional absolute boom override.
/// Returns true if any knob changed vs the inputs.
inline bool apply_from_env(int& boom_interval,
                           double& sharing_scale,
                           double& boom_fraction) noexcept
{
	const int boom0 = boom_interval;
	const double scale0 = sharing_scale;
	const double frac0 = boom_fraction;
	apply(boom_interval, sharing_scale, boom_fraction, enabled_from_env());
	const int ov = boom_interval_override_from_env();
	if (ov > 0)
		boom_interval = ov;
	return boom_interval != boom0 || sharing_scale != scale0 ||
	       boom_fraction != frac0;
}

/// Legacy 2-arg env apply (does not restore boom_fraction — prefer 3-arg).
inline bool apply_from_env(int& boom_interval, double& sharing_scale) noexcept
{
	double frac = 0.0;
	const bool en = enabled_from_env();
	// When coverage on, restore local frac so interval logic is exercised;
	// caller fraction is unchanged (use 3-arg for real BOOM restore).
	const int boom0 = boom_interval;
	const double scale0 = sharing_scale;
	apply(boom_interval, sharing_scale, frac, en);
	const int ov = boom_interval_override_from_env();
	if (ov > 0)
		boom_interval = ov;
	return boom_interval != boom0 || sharing_scale != scale0;
}

}  // namespace search_coverage
}  // namespace flexaids
