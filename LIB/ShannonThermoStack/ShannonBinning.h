// ShannonBinning.h — the single definition of the Shannon histogram support
//
// Two independent estimators for the same quantity live in this tree:
//   shannon_thermo::compute_shannon_entropy      (ShannonThermoStack.cpp)
//   UnifiedHardwareDispatch::compute_shannon_entropy (per-backend kernels)
//
// They each derived their own support and their own bin index, so once one
// gained an outlier-robust support the two disagreed by ~2 nats on the same
// input — enough to put them on opposite sides of the GA collapse gate. The
// dispatch shells stay separate (the per-backend kernels exist so
// benchmark_dispatch can time them individually), but the SUPPORT and the BIN
// INDEX are defined once, here, and included by both.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "../EnvFlags.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace shannon_thermo {

// ─── robust support ──────────────────────────────────────────────────────────
//
// The histogram support is normally the sample's own [min, max]. That lets a
// single extreme value rescale every bin: one clash/wall pose (evalue ~1e4) in
// an otherwise diverse population pushes every real sample into bin 0, so H
// reads ~0.08 bits and a caller's collapse gate fires on a population that has
// not collapsed at all.
//
// When the raw range is much wider than the bulk, the support is instead taken
// from a Tukey far-out fence around the interquartile range. Values outside the
// fence are not discarded: bin_index() clamps, so they land in the edge bins and
// keep their probability mass.
//
// Engagement is range/IQR based and the raw range is not itself robust, so which
// samples trigger the fence is tail-dependent rather than merely "pathological".
// Measured firing rates: 0% on Gaussian (n = 1e3 and 1e5), uniform and
// exponential — those stay bit-identical to the pre-fence estimator — but ~100%
// on heavy-tailed samples and on any population carrying clash poses, where H
// rises by 1.5-2.3 nats. That is the intended repair, but it is a real shift in
// the GA collapse gate's operating point, so it is exposed as an A/B arm:
// FLEXAIDDS_SHANNON_ROBUST=0 restores the previous raw min/max support.

inline constexpr std::size_t kRobustMinSamples = 8;    // below this, quartiles are meaningless
inline constexpr double      kRobustTriggerIQR = 20.0; // engage only well outside the bulk
inline constexpr double      kRobustFenceIQR   = 3.0;  // Tukey "far out" fence

inline bool robust_support_enabled() noexcept
{
    static const bool enabled = flexaids::env_bool("FLEXAIDDS_SHANNON_ROBUST", true);
    return enabled;
}

/// Returns true (and fills lo/hi) when a robust support should replace [min,max].
inline bool robust_support(const std::vector<double>& values,
                           double raw_min, double raw_max,
                           double& lo, double& hi)
{
    const std::size_t n = values.size();
    if (n < kRobustMinSamples) return false;
    if (!robust_support_enabled()) return false;

    std::vector<double> scratch(values);
    auto quantile = [&scratch, n](double f) {
        const std::size_t idx =
            static_cast<std::size_t>(f * static_cast<double>(n - 1));
        // nth_element is valid on any permutation, so successive calls on the
        // same (partially reordered) buffer still return the correct element.
        std::nth_element(scratch.begin(),
                         scratch.begin() + static_cast<std::ptrdiff_t>(idx),
                         scratch.end());
        return scratch[idx];
    };

    const double q1  = quantile(0.25);
    const double q3  = quantile(0.75);
    const double iqr = q3 - q1;
    if (!(iqr > 0.0) || !std::isfinite(iqr)) return false;

    if ((raw_max - raw_min) <= kRobustTriggerIQR * iqr) return false;

    lo = q1 - kRobustFenceIQR * iqr;
    hi = q3 + kRobustFenceIQR * iqr;
    // Must clear the +1e-10 bin-width floor applied by the callers, not merely
    // be nonzero: a fence narrower than num_bins*1e-10 is swallowed by that
    // epsilon, so every bulk sample collapses into bin 0 and the fence silently
    // does nothing while still paying for the copy.
    return (hi - lo) > 1e-9;
}

// ─── bin index ───────────────────────────────────────────────────────────────
//
// The clamp is applied in DOUBLE, before the narrowing conversion.
//
// With a raw min/max support every value satisfied 0 <= (v-min)*inv_bw <=
// num_bins, so casting first and clamping the int was safe. Once the support can
// be a fence, out-of-fence values are arbitrarily far outside it: a tight bulk
// plus one clash pose at 1e4 produces a raw index above 5e10, and converting
// that to int is undefined behaviour. In practice x86-64 cvttsd2si yields
// INT_MIN, which then clamps to bin 0 — placing the clash pose ON TOP OF the
// bulk, collapsing the histogram to a single occupied bin and firing the very
// gate the fence exists to protect — while ARM saturates to INT_MAX and lands it
// in the top bin. Clamping in double removes both the undefined behaviour and
// the platform divergence.
inline int bin_index(double v, double min_v, double inv_bw, int num_bins) noexcept
{
    const double t = (v - min_v) * inv_bw;
    if (!(t > 0.0)) return 0;                      // also catches NaN
    const double top = static_cast<double>(num_bins - 1);
    return (t >= top) ? num_bins - 1 : static_cast<int>(t);
}

/// Resolve the histogram support for `values`, applying the robust fence when it
/// engages. Returns false when the sample is degenerate (zero width).
inline bool histogram_support(const std::vector<double>& values,
                              double& min_v, double& max_v)
{
    if (values.empty()) return false;
    const auto [it_min, it_max] = std::minmax_element(values.begin(), values.end());
    min_v = *it_min;
    max_v = *it_max;
    if (max_v - min_v < 1e-12) return false;

    double lo = 0.0, hi = 0.0;
    if (robust_support(values, min_v, max_v, lo, hi)) {
        min_v = lo;
        max_v = hi;
    }
    return true;
}

}  // namespace shannon_thermo
