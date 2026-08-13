// get_yval.h — piecewise-linear energy-matrix lookup + optional LUT
//
// FLEXAIDDS_GET_YVAL_LUT default OFF: get_yval() is the linear scan (bit-identical
// to the pre-LUT vcfunction.cpp implementation). ON: 256-bin lerp of the same
// scan, resolution knob FLEXAIDDS_GET_YVAL_LUT_BINS (default 256, clamp 16..1024).
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "flexaid.h"
#include "EnvFlags.h"

#include <cmath>
#include <cstdlib>

namespace flexaids {

inline bool get_yval_lut_enabled() noexcept
{
    return env_bool("FLEXAIDDS_GET_YVAL_LUT", false);
}

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

/// Linear-scan interpolation (the historical get_yval body).
double get_yval_scan(struct energy_matrix* em, double relative_area);

/// Dispatch: scan unless FLEXAIDDS_GET_YVAL_LUT is on.
double get_yval(struct energy_matrix* em, double relative_area);
