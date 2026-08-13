// pprop.h — NRGRank v2 rank-fraction (pProp) / ΔpProp (default unused)
//
// pProp = rank / N (1-based rank). ΔpProp = pProp_this − pProp_ref.
// FLEXAIDDS_PPROP_MAX unset → no filter. When set to (0,1], drop ligands
// with pProp above that cap (selectivity / vPAINS-style prune).
//
// Method: DesCôteaux, Mailhot, Najmanovich bioRxiv 2025.02.17.638675v2 §2.5.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "EnvFlags.h"

#include <cstdlib>

namespace flexaids {

inline double pprop(int rank_1based, int n)
{
    if (n <= 0) return 1.0;
    if (rank_1based < 1) rank_1based = 1;
    return static_cast<double>(rank_1based) / static_cast<double>(n);
}

inline double delta_pprop(double p_this, double p_ref)
{
    return p_this - p_ref;
}

/// NaN if unset / unparseable (no filter).
inline double pprop_max_cap() noexcept
{
    const char* s = std::getenv("FLEXAIDDS_PPROP_MAX");
    if (!s || !*s) return -1.0;
    char* end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s || v <= 0.0 || v > 1.0) return -1.0;
    return v;
}

inline bool pprop_keep(int rank_1based, int n) noexcept
{
    const double cap = pprop_max_cap();
    if (cap < 0.0) return true;
    return pprop(rank_1based, n) <= cap;
}

}  // namespace flexaids
