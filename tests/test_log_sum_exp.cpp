// tests/test_log_sum_exp.cpp — sequential LSE used when FIXED_ORDER_LSE is on
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "log_sum_exp.h"

#include <cmath>
#include <cstdlib>
#include <vector>

TEST(FixedOrderLse, DefaultOff) {
#ifdef _WIN32
    _putenv_s("FLEXAIDDS_FIXED_ORDER_LSE", "");
#else
    unsetenv("FLEXAIDDS_FIXED_ORDER_LSE");
#endif
    EXPECT_FALSE(flexaids::fixed_order_lse_enabled());
}

TEST(FixedOrderLse, MatchesClosedFormTwoEqual) {
    const double xs[] = {1.0, 1.0};
    const double got = flexaids::log_sum_exp_fixed_order(xs);
    EXPECT_NEAR(got, 1.0 + std::log(2.0), 1e-12);
}

TEST(FixedOrderLse, OrderIndependent) {
    const double a[] = {3.0, 1.0, -4.0, 2.5};
    const double b[] = {-4.0, 2.5, 3.0, 1.0};
    EXPECT_NEAR(flexaids::log_sum_exp_fixed_order(a),
                flexaids::log_sum_exp_fixed_order(b), 1e-12);
}

TEST(FixedOrderLse, EmptyIsNegInf) {
    std::span<const double> empty{};
    EXPECT_TRUE(std::isinf(flexaids::log_sum_exp_fixed_order(empty)));
}
