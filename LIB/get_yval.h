// get_yval.h — piecewise-linear energy-matrix lookup + optional LUT
//
// Default path: get_yval_scan() — binary search + lerp on the flat knot table
// (same knots/slopes as the historical linear scan; LUT still OFF).
// FLEXAIDDS_GET_YVAL_LUT default OFF (METHODOLOGY.md §1): get_yval() is that
// scan. ON: 256-bin lerp of the same scan, resolution knob
// FLEXAIDDS_GET_YVAL_LUT_BINS (default 256, clamp 16..1024). LUT stays opt-in
// until a §1 parity run shows ranking-preserving / bit-identical CF vs LUT-OFF.
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flexaid.h"
#include "EnvFlags.h"

#include <cmath>
#include <cstdlib>

namespace flexaids {

/// Live getenv. Hot-loop code must use get_yval_lut_enabled_cached() instead.
inline bool get_yval_lut_enabled() noexcept
{
    return env_bool("FLEXAIDDS_GET_YVAL_LUT", false);
}

/// Live getenv. Hot-loop code must use get_yval_lut_bins_cached() instead.
inline int get_yval_lut_bins() noexcept
{
    const char* s = std::getenv("FLEXAIDDS_GET_YVAL_LUT_BINS");
    if (!s || !*s) return 256;
    int n = std::atoi(s);
    if (n < 16) n = 16;
    if (n > 1024) n = 1024;
    return n;
}

}  // namespace flexaids

/// Piecewise-linear interpolation on em->flat_* (binary search on knots).
double get_yval_scan(struct energy_matrix* em, double relative_area);

/// Dispatch: scan unless the cached FLEXAIDDS_GET_YVAL_LUT snapshot is on.
double get_yval(struct energy_matrix* em, double relative_area);

/// Same dispatch with an explicit LUT flag (vcfunction hoists the snapshot;
/// tests can exercise LUT ON without mutating the process-wide cache).
double get_yval(struct energy_matrix* em, double relative_area, bool use_lut);
