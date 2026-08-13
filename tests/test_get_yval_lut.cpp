// tests/test_get_yval_lut.cpp — scan vs LUT dispatch for the shipped get_yval
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "get_yval.h"

#include <cmath>
#include <cstdlib>
#include <vector>

namespace {

void set_env(const char* k, const char* v) {
#ifdef _WIN32
    _putenv_s(k, v);
#else
    setenv(k, v, 1);
#endif
}
void unset_env(const char* k) {
#ifdef _WIN32
    _putenv_s(k, "");
#else
    unsetenv(k);
#endif
}

struct Piecewise {
    energy_values ev{};
    energy_matrix em{};
    std::vector<float> xs, ys, sl;
    Piecewise() {
        ev.x = 0; ev.y = 0; ev.next_value = nullptr;
        xs = {0.f, 0.5f, 1.f};
        ys = {0.f, 10.f, 4.f};
        sl = {(10.f - 0.f) / 0.5f, (4.f - 10.f) / 0.5f};
        em.weight = 0;
        em.energy_values = &ev;
        em.flat_n = 3;
        em.flat_x = xs.data();
        em.flat_y = ys.data();
        em.flat_slope = sl.data();
    }
};

}  // namespace

TEST(GetYval, DefaultOffUsesScan) {
    unset_env("FLEXAIDDS_GET_YVAL_LUT");
    EXPECT_FALSE(flexaids::get_yval_lut_enabled());
    Piecewise p;
    EXPECT_DOUBLE_EQ(get_yval(&p.em, 0.25), get_yval_scan(&p.em, 0.25));
    EXPECT_NEAR(get_yval(&p.em, 0.25), 5.0, 1e-5);
}

TEST(GetYval, LutOnInterpolatesScanSamples) {
    set_env("FLEXAIDDS_GET_YVAL_LUT", "1");
    EXPECT_TRUE(flexaids::get_yval_lut_enabled());
    Piecewise p;
    const double scan = get_yval_scan(&p.em, 0.25);
    const double lut = get_yval(&p.em, 0.25);
    EXPECT_NEAR(lut, scan, 0.1);  // 256-bin lerp of the same scan
    EXPECT_DOUBLE_EQ(get_yval(&p.em, -0.5), get_yval_scan(&p.em, -0.5));
    unset_env("FLEXAIDDS_GET_YVAL_LUT");
}

TEST(GetYval, OutOfRangeFallsBackToScan) {
    set_env("FLEXAIDDS_GET_YVAL_LUT", "1");
    Piecewise p;
    EXPECT_DOUBLE_EQ(get_yval(&p.em, 1.5), get_yval_scan(&p.em, 1.5));
    unset_env("FLEXAIDDS_GET_YVAL_LUT");
}

TEST(GetYval, TwoMatricesLutCopiesUnderLock) {
    set_env("FLEXAIDDS_GET_YVAL_LUT", "1");
    Piecewise a, b;
    b.ys = {0.f, 20.f, 8.f};
    b.sl = {(20.f - 0.f) / 0.5f, (8.f - 20.f) / 0.5f};
    // Inserting a second key must not invalidate interpolation of the first.
    const double a0 = get_yval(&a.em, 0.25);
    const double b0 = get_yval(&b.em, 0.25);
    const double a1 = get_yval(&a.em, 0.25);
    EXPECT_DOUBLE_EQ(a0, a1);
    EXPECT_NEAR(a0, get_yval_scan(&a.em, 0.25), 0.1);
    EXPECT_NEAR(b0, get_yval_scan(&b.em, 0.25), 0.1);
    unset_env("FLEXAIDDS_GET_YVAL_LUT");
}
