// tests/test_sampling_coverage.cpp — coarse-init / seed-selection pure helpers
//
// Drives LIB/sampling_coverage.h (shipped path used by coarse_init.cpp).
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include "sampling_coverage.h"

#include <gtest/gtest.h>

#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

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

TEST(SamplingCoverage, DefaultBoostOffLeavesParams)
{
    ScopedEnv boost("FLEXAIDDS_SAMPLE_COVERAGE_BOOST", nullptr);
    auto p = flexaids::sampling::apply_coverage_boost(64, 25);
    EXPECT_EQ(p.n_orient, 64);
    EXPECT_EQ(p.n_seeds, 25);
    EXPECT_FALSE(p.boost_applied);
}

TEST(SamplingCoverage, BoostDoublesOrientAndSeeds)
{
    ScopedEnv boost("FLEXAIDDS_SAMPLE_COVERAGE_BOOST", "1");
    ScopedEnv om("FLEXAIDDS_SAMPLE_ORIENT_MULT", nullptr);
    ScopedEnv sm("FLEXAIDDS_SAMPLE_SEEDS_MULT", nullptr);
    auto p = flexaids::sampling::apply_coverage_boost(64, 25);
    EXPECT_TRUE(p.boost_applied);
    EXPECT_EQ(p.n_orient, 128);
    EXPECT_EQ(p.n_seeds, 50);
}

TEST(SamplingCoverage, SelectDropsClashSentinel)
{
    // CF ascending: good, ok, clash
    const std::vector<double> cf = {-40.0, -10.0, 10000.0, 10000.0};
    auto idx = flexaids::sampling::select_ranked_seed_indices(
        cf, /*n_keep=*/4, /*clash=*/1e4, /*force_ranked=*/false);
    ASSERT_EQ(idx.size(), 2u);
    EXPECT_EQ(idx[0], 0u);
    EXPECT_EQ(idx[1], 1u);
}

TEST(SamplingCoverage, ForceRankedInjectsClashScale)
{
    const std::vector<double> cf = {10000.0, 10000.0, 10000.0};
    auto empty = flexaids::sampling::select_ranked_seed_indices(
        cf, 2, 1e4, /*force_ranked=*/false);
    EXPECT_TRUE(empty.empty());

    auto forced = flexaids::sampling::select_ranked_seed_indices(
        cf, 2, 1e4, /*force_ranked=*/true);
    ASSERT_EQ(forced.size(), 2u);
}

TEST(SamplingCoverage, RankOrderStableByCf)
{
    const std::vector<double> cf = {-1.0, -50.0, -20.0};
    auto ranked = flexaids::sampling::rank_indices_by_cf_asc(cf);
    ASSERT_EQ(ranked.size(), 3u);
    EXPECT_EQ(ranked[0], 1u);  // -50
    EXPECT_EQ(ranked[1], 2u);  // -20
    EXPECT_EQ(ranked[2], 0u);  // -1
}

TEST(SamplingCoverage, DiversifyKeepsDistantSeeds)
{
    std::vector<std::array<double, 3>> coords = {
        {0.0, 0.0, 0.0},
        {0.1, 0.0, 0.0},   // near first
        {10.0, 0.0, 0.0},  // far
    };
    // Prefer index order 0,1,2 (as if CF-ranked that way)
    const std::vector<std::size_t> ranked = {0, 1, 2};
    auto kept = flexaids::sampling::diversify_by_min_distance(
        ranked, &coords, /*n_keep=*/2, /*min_dist_A=*/1.0);
    ASSERT_EQ(kept.size(), 2u);
    EXPECT_EQ(kept[0], 0u);
    EXPECT_EQ(kept[1], 2u);  // skip near-duplicate 1
}
