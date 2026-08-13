// log_sum_exp.h — sequential (thread-count-independent) log-sum-exp
//
// FLEXAIDDS_FIXED_ORDER_LSE default OFF: callers keep the hardware-dispatched
// reduction. ON: this sequential max+sum, independent of OpenMP/AVX reduction
// order.
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "EnvFlags.h"

#include <cmath>
#include <limits>
#include <span>

namespace flexaids {

inline bool fixed_order_lse_enabled() noexcept
{
    return env_bool("FLEXAIDDS_FIXED_ORDER_LSE", false);
}

inline double log_sum_exp_fixed_order(std::span<const double> values)
{
    if (values.empty())
        return -std::numeric_limits<double>::infinity();
    double xmax = values[0];
    for (std::size_t i = 1; i < values.size(); ++i)
        if (values[i] > xmax) xmax = values[i];
    if (!std::isfinite(xmax))
        return xmax;
    double sum = 0.0;
    for (double v : values)
        sum += std::exp(v - xmax);
    return xmax + std::log(sum);
}

}  // namespace flexaids
