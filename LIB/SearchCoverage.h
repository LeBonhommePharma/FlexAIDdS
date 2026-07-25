// SearchCoverage.h — flag-gated GA search-coverage knobs (S1 remediation)
//
// FLEXAIDDS_SEARCH_COVERAGE=1 (default OFF) densifies BOOM injection and
// widens niche radius so near-native basins are less likely to starve under
// over-burial / false-minima dominance. Targets bulk BCR>2 failures
// (sampling_fail + sampling_near). Does not change product defaults.
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

/// Apply S1 coverage knobs in-place.
/// - boom_interval: halved, floor 25 (when original > 0 and enabled)
/// - sharing_scale: ×0.5 when enabled (sig_share ∝ 1/scale → wider niches)
/// No-op when enabled is false.
inline void apply(int& boom_interval, double& sharing_scale, bool enabled) noexcept
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
}

/// Compose env-driven coverage + optional absolute boom override.
/// Returns true if any knob changed vs the inputs.
inline bool apply_from_env(int& boom_interval, double& sharing_scale) noexcept
{
	const int boom0 = boom_interval;
	const double scale0 = sharing_scale;
	apply(boom_interval, sharing_scale, enabled_from_env());
	const int ov = boom_interval_override_from_env();
	if (ov > 0)
		boom_interval = ov;
	return boom_interval != boom0 || sharing_scale != scale0;
}

}  // namespace search_coverage
}  // namespace flexaids
