// tests/test_search_coverage.cpp — S1 FLEXAIDDS_SEARCH_COVERAGE pure knobs
//
// Drives LIB/SearchCoverage.h (shipped path). Default OFF must leave
// boom/scale/fraction unchanged; ON halves boom interval, halves sharing_scale,
// and restores boom_inject_fraction when campaign configs left it at 0.
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include "SearchCoverage.h"

#include <gtest/gtest.h>

#include <cstdlib>
#include <optional>
#include <string>

namespace {

class ScopedEnv {
public:
    ScopedEnv(const char* key, const char* val) : key_(key)
    {
        const char* prev = std::getenv(key);
        if (prev)
            prev_ = prev;
        if (val)
            ::setenv(key, val, 1);
        else
            ::unsetenv(key);
    }
    ~ScopedEnv()
    {
        if (prev_)
            ::setenv(key_, prev_->c_str(), 1);
        else
            ::unsetenv(key_);
    }

private:
    const char* key_;
    std::optional<std::string> prev_;
};

}  // namespace

TEST(SearchCoverage, DefaultOffIsNoOp)
{
    int boom = 100;
    double scale = 10.0;
    double frac = 0.0;
    flexaids::search_coverage::apply(boom, scale, frac, /*enabled=*/false);
    EXPECT_EQ(boom, 100);
    EXPECT_DOUBLE_EQ(scale, 10.0);
    EXPECT_DOUBLE_EQ(frac, 0.0);
}

TEST(SearchCoverage, EnabledHalvesBoomAndScale)
{
    int boom = 100;
    double scale = 10.0;
    double frac = 1.0;
    flexaids::search_coverage::apply(boom, scale, frac, /*enabled=*/true);
    EXPECT_EQ(boom, 50);           // 100/2
    EXPECT_DOUBLE_EQ(scale, 5.0);  // ×0.5 → wider sig_share
    EXPECT_DOUBLE_EQ(frac, 1.0);   // already positive — unchanged
}

TEST(SearchCoverage, BoomFloorIs25)
{
    int boom = 40;
    double scale = 8.0;
    double frac = 1.0;
    flexaids::search_coverage::apply(boom, scale, frac, true);
    EXPECT_EQ(boom, 25);  // max(25, 40/2=20)
    EXPECT_DOUBLE_EQ(scale, 4.0);
}

TEST(SearchCoverage, ZeroBoomStaysOffCadence)
{
    int boom = 0;
    double scale = 10.0;
    double frac = 1.0;
    flexaids::search_coverage::apply(boom, scale, frac, true);
    EXPECT_EQ(boom, 0);  // 0 means BOOM injection disabled
    EXPECT_DOUBLE_EQ(scale, 5.0);
}

TEST(SearchCoverage, RestoresZeroBoomFraction)
{
    // Campaign inventory mode: interval>0 but fraction=0 → BOOM is a no-op.
    int boom = 100;
    double scale = 10.0;
    double frac = 0.0;
    const int inject_before =
        flexaids::search_coverage::inject_count_per_event(boom, frac, 1000);
    EXPECT_EQ(inject_before, 0);

    flexaids::search_coverage::apply(boom, scale, frac, /*enabled=*/true);
    EXPECT_EQ(boom, 50);
    EXPECT_DOUBLE_EQ(frac, flexaids::search_coverage::kDefaultBoomFraction);
    const int inject_after =
        flexaids::search_coverage::inject_count_per_event(boom, frac, 1000);
    // half of 1000 = 500; fraction 1.0 → 500 injects per event (>0).
    EXPECT_EQ(inject_after, 500);
    EXPECT_GT(inject_after, inject_before);
}

TEST(SearchCoverage, TotalInjectsIncreaseWhenCoverageOn)
{
    // Representative campaign: interval=100, fraction=0, pop=1000, gen=2000.
    const int num_chrom = 1000;
    const int max_gen = 2000;
    int boom = 100;
    double scale = 10.0;
    double frac = 0.0;
    const long long total_before =
        flexaids::search_coverage::total_injects_estimate(
            boom, frac, num_chrom, max_gen);
    EXPECT_EQ(total_before, 0);

    flexaids::search_coverage::apply(boom, scale, frac, true);
    const long long total_after =
        flexaids::search_coverage::total_injects_estimate(
            boom, frac, num_chrom, max_gen);
    // interval 50 → events at 50,100,...,1950 = 39 events × 500 = 19500
    EXPECT_GT(total_after, total_before);
    EXPECT_EQ(total_after, 39LL * 500LL);
}

TEST(SearchCoverage, EnvCoverageChangesKnobsAndFraction)
{
    ScopedEnv cov("FLEXAIDDS_SEARCH_COVERAGE", "1");
    ScopedEnv boom_ov("FLEXAIDDS_BOOM_INTERVAL", nullptr);
    int boom = 100;
    double scale = 10.0;
    double frac = 0.0;
    EXPECT_TRUE(flexaids::search_coverage::apply_from_env(boom, scale, frac));
    EXPECT_EQ(boom, 50);
    EXPECT_DOUBLE_EQ(scale, 5.0);
    EXPECT_DOUBLE_EQ(frac, 1.0);
}

TEST(SearchCoverage, BoomIntervalOverrideWins)
{
    ScopedEnv cov("FLEXAIDDS_SEARCH_COVERAGE", "1");
    ScopedEnv boom_ov("FLEXAIDDS_BOOM_INTERVAL", "33");
    int boom = 100;
    double scale = 10.0;
    double frac = 0.0;
    EXPECT_TRUE(flexaids::search_coverage::apply_from_env(boom, scale, frac));
    EXPECT_EQ(boom, 33);
    EXPECT_DOUBLE_EQ(scale, 5.0);
    EXPECT_DOUBLE_EQ(frac, 1.0);
}

TEST(SearchCoverage, EnvDefaultDoesNotChange)
{
    ScopedEnv cov("FLEXAIDDS_SEARCH_COVERAGE", nullptr);
    ScopedEnv boom_ov("FLEXAIDDS_BOOM_INTERVAL", nullptr);
    int boom = 100;
    double scale = 10.0;
    double frac = 0.0;
    EXPECT_FALSE(flexaids::search_coverage::enabled_from_env());
    EXPECT_FALSE(flexaids::search_coverage::apply_from_env(boom, scale, frac));
    EXPECT_EQ(boom, 100);
    EXPECT_DOUBLE_EQ(scale, 10.0);
    EXPECT_DOUBLE_EQ(frac, 0.0);
}
