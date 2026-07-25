// tests/test_search_coverage.cpp — S1 FLEXAIDDS_SEARCH_COVERAGE pure knobs
//
// Drives LIB/SearchCoverage.h (shipped path). Default OFF must leave
// boom/scale unchanged; ON halves boom (floor 25) and halves sharing_scale.
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include "SearchCoverage.h"

#include <gtest/gtest.h>

#include <cstdlib>
#include <optional>
#include <string>

namespace {

// Set/clear env for one test (restore on destruction).
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
    flexaids::search_coverage::apply(boom, scale, /*enabled=*/false);
    EXPECT_EQ(boom, 100);
    EXPECT_DOUBLE_EQ(scale, 10.0);
}

TEST(SearchCoverage, EnabledHalvesBoomAndScale)
{
    int boom = 100;
    double scale = 10.0;
    flexaids::search_coverage::apply(boom, scale, /*enabled=*/true);
    EXPECT_EQ(boom, 50);           // 100/2
    EXPECT_DOUBLE_EQ(scale, 5.0);  // ×0.5 → wider sig_share
}

TEST(SearchCoverage, BoomFloorIs25)
{
    int boom = 40;
    double scale = 8.0;
    flexaids::search_coverage::apply(boom, scale, true);
    EXPECT_EQ(boom, 25);  // max(25, 40/2=20)
    EXPECT_DOUBLE_EQ(scale, 4.0);
}

TEST(SearchCoverage, ZeroBoomStaysOffCadence)
{
    int boom = 0;
    double scale = 10.0;
    flexaids::search_coverage::apply(boom, scale, true);
    EXPECT_EQ(boom, 0);  // 0 means BOOM injection disabled
    EXPECT_DOUBLE_EQ(scale, 5.0);
}

TEST(SearchCoverage, EnvCoverageChangesKnobs)
{
    ScopedEnv cov("FLEXAIDDS_SEARCH_COVERAGE", "1");
    ScopedEnv boom_ov("FLEXAIDDS_BOOM_INTERVAL", nullptr);
    int boom = 100;
    double scale = 10.0;
    EXPECT_TRUE(flexaids::search_coverage::apply_from_env(boom, scale));
    EXPECT_EQ(boom, 50);
    EXPECT_DOUBLE_EQ(scale, 5.0);
}

TEST(SearchCoverage, BoomIntervalOverrideWins)
{
    ScopedEnv cov("FLEXAIDDS_SEARCH_COVERAGE", "1");
    ScopedEnv boom_ov("FLEXAIDDS_BOOM_INTERVAL", "33");
    int boom = 100;
    double scale = 10.0;
    EXPECT_TRUE(flexaids::search_coverage::apply_from_env(boom, scale));
    EXPECT_EQ(boom, 33);
    EXPECT_DOUBLE_EQ(scale, 5.0);
}

TEST(SearchCoverage, EnvDefaultDoesNotChange)
{
    ScopedEnv cov("FLEXAIDDS_SEARCH_COVERAGE", nullptr);
    ScopedEnv boom_ov("FLEXAIDDS_BOOM_INTERVAL", nullptr);
    int boom = 100;
    double scale = 10.0;
    EXPECT_FALSE(flexaids::search_coverage::enabled_from_env());
    EXPECT_FALSE(flexaids::search_coverage::apply_from_env(boom, scale));
    EXPECT_EQ(boom, 100);
    EXPECT_DOUBLE_EQ(scale, 10.0);
}
