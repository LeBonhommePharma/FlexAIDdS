// tests/test_get_yval_lut.cpp — scan vs LUT dispatch for the shipped get_yval
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>
#include "get_yval.h"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <thread>
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

// Historical linear scan — the binary-search path must match this exactly.
double scan_linear(energy_matrix* em, double relative_area)
{
    if (!em || !em->energy_values) return 0.0;
    if (em->weight) return (double)em->energy_values->y;
    const int n = em->flat_n;
    if (n == 0) return 0.0;
    const float ra = (float)relative_area;
    const float* fx = em->flat_x;
    const float* fy = em->flat_y;
    if (!fx || !fy) return 0.0;
    if (ra < fx[0]) return 0.0;
    if (ra >= fx[n - 1]) return (double)fy[n - 1];
    int i = 0;
    while (i < n - 2 && ra >= fx[i + 1]) ++i;
    return (double)(fy[i] + em->flat_slope[i] * (ra - fx[i]));
}

}  // namespace

TEST(GetYval, DefaultOffUsesScan) {
    unset_env("FLEXAIDDS_GET_YVAL_LUT");
    EXPECT_FALSE(flexaids::get_yval_lut_enabled());
    Piecewise p;
    EXPECT_DOUBLE_EQ(get_yval(&p.em, 0.25, false), get_yval_scan(&p.em, 0.25));
    EXPECT_NEAR(get_yval(&p.em, 0.25, false), 5.0, 1e-5);
}

TEST(GetYval, CachedLutFlagIgnoresLaterGetenv) {
    const bool cached = flexaids::get_yval_lut_enabled_cached();
    Piecewise p;
    const double via_cache = get_yval(&p.em, 0.25);
    const double via_flag = get_yval(&p.em, 0.25, cached);
    EXPECT_DOUBLE_EQ(via_cache, via_flag);

    if (cached) {
        unset_env("FLEXAIDDS_GET_YVAL_LUT");
    } else {
        set_env("FLEXAIDDS_GET_YVAL_LUT", "1");
    }
    EXPECT_EQ(flexaids::get_yval_lut_enabled_cached(), cached);
    EXPECT_EQ(flexaids::get_yval_lut_enabled(), !cached)
        << "live getenv must be independent of the hot-loop snapshot";
    EXPECT_DOUBLE_EQ(get_yval(&p.em, 0.25), via_flag);

    if (cached) {
        set_env("FLEXAIDDS_GET_YVAL_LUT", "1");
    } else {
        unset_env("FLEXAIDDS_GET_YVAL_LUT");
    }
}

TEST(GetYval, LutOnInterpolatesScanSamples) {
    set_env("FLEXAIDDS_GET_YVAL_LUT", "1");
    EXPECT_TRUE(flexaids::get_yval_lut_enabled());
    Piecewise p;
    const double scan = get_yval_scan(&p.em, 0.25);
    const double lut = get_yval(&p.em, 0.25, true);
    EXPECT_NEAR(lut, scan, 0.1);  // 256-bin lerp of the same scan
    EXPECT_DOUBLE_EQ(get_yval(&p.em, -0.5, true), get_yval_scan(&p.em, -0.5));
    unset_env("FLEXAIDDS_GET_YVAL_LUT");
}

TEST(GetYval, OutOfRangeFallsBackToScan) {
    Piecewise p;
    EXPECT_DOUBLE_EQ(get_yval(&p.em, 1.5, true), get_yval_scan(&p.em, 1.5));
}

TEST(GetYval, TwoMatricesLutIndependent) {
    Piecewise a, b;
    b.ys = {0.f, 20.f, 8.f};
    b.sl = {(20.f - 0.f) / 0.5f, (8.f - 20.f) / 0.5f};
    const double a0 = get_yval(&a.em, 0.25, true);
    const double b0 = get_yval(&b.em, 0.25, true);
    const double a1 = get_yval(&a.em, 0.25, true);
    EXPECT_DOUBLE_EQ(a0, a1);
    EXPECT_NEAR(a0, get_yval_scan(&a.em, 0.25), 0.1);
    EXPECT_NEAR(b0, get_yval_scan(&b.em, 0.25), 0.1);
}

TEST(GetYval, BinarySearchMatchesLinearScan) {
    energy_values ev{};
    ev.next_value = nullptr;
    std::vector<float> xs{0.f, 0.2f, 0.5f, 0.75f, 1.f};
    std::vector<float> ys{1.f, 3.f, -1.f, 4.f, 2.f};
    std::vector<float> sl(xs.size() - 1);
    for (size_t i = 0; i < sl.size(); ++i)
        sl[i] = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]);
    energy_matrix em{};
    em.weight = 0;
    em.energy_values = &ev;
    em.flat_n = static_cast<int>(xs.size());
    em.flat_x = xs.data();
    em.flat_y = ys.data();
    em.flat_slope = sl.data();

    const double samples[] = {
        -0.5, 0.0, 0.1, 0.2, 0.35, 0.5, 0.625, 0.75, 0.9, 1.0, 1.5,
    };
    for (double ra : samples) {
        EXPECT_DOUBLE_EQ(get_yval_scan(&em, ra), scan_linear(&em, ra)) << "ra=" << ra;
    }
    for (int k = 0; k <= 40; ++k) {
        const double ra = -0.05 + 0.03 * static_cast<double>(k);
        EXPECT_DOUBLE_EQ(get_yval_scan(&em, ra), scan_linear(&em, ra)) << "ra=" << ra;
    }

    Piecewise p;
    EXPECT_DOUBLE_EQ(get_yval_scan(&p.em, 0.0), scan_linear(&p.em, 0.0));
    EXPECT_DOUBLE_EQ(get_yval_scan(&p.em, 0.25), scan_linear(&p.em, 0.25));
    EXPECT_DOUBLE_EQ(get_yval_scan(&p.em, 0.5), scan_linear(&p.em, 0.5));
    EXPECT_DOUBLE_EQ(get_yval_scan(&p.em, 0.75), scan_linear(&p.em, 0.75));
    EXPECT_DOUBLE_EQ(get_yval_scan(&p.em, 1.0), scan_linear(&p.em, 1.0));
}

TEST(GetYval, LutLookupsConcurrentSafe) {
    Piecewise a, b;
    b.ys = {0.f, 20.f, 8.f};
    b.sl = {(20.f - 0.f) / 0.5f, (8.f - 20.f) / 0.5f};
    std::atomic<int> mismatches{0};
    auto worker = [&](int seed) {
        energy_matrix* ems[2] = {&a.em, &b.em};
        for (int k = 0; k < 2000; ++k) {
            energy_matrix* em = ems[(seed + k) & 1];
            const double ra = 0.01 * static_cast<double>((k * 17 + seed) % 100);
            const double lut = get_yval(em, ra, true);
            const double scan = get_yval_scan(em, ra);
            if (std::abs(lut - scan) > 0.15)
                mismatches.fetch_add(1, std::memory_order_relaxed);
        }
    };
    std::vector<std::thread> ts;
    ts.reserve(8);
    for (int t = 0; t < 8; ++t)
        ts.emplace_back(worker, t);
    for (auto& th : ts)
        th.join();
    EXPECT_EQ(mismatches.load(), 0);
}
