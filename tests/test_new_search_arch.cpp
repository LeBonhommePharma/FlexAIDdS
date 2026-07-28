// Unit tests for S4 new_search_arch A+B (drives shipped LIB/new_search_arch.h).
#include "new_search_arch.h"

#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

namespace {

void clear_env(const char* key) {
#if defined(_WIN32)
    _putenv_s(key, "");
#else
    unsetenv(key);
#endif
}

void set_env(const char* key, const char* val) {
#if defined(_WIN32)
    _putenv_s(key, val);
#else
    setenv(key, val, 1);
#endif
}

class NewSearchArchEnv : public ::testing::Test {
protected:
    void SetUp() override {
        clear_env("FLEXAIDDS_NEW_SEARCH_ARCH");
        clear_env("FLEXAIDDS_PHENOTYPE_UNIQUE");
        clear_env("FLEXAIDDS_BASIN_REINJECT");
        clear_env("FLEXAIDDS_BASIN_SIGMA_ANG");
    }
    void TearDown() override {
        clear_env("FLEXAIDDS_NEW_SEARCH_ARCH");
        clear_env("FLEXAIDDS_PHENOTYPE_UNIQUE");
        clear_env("FLEXAIDDS_BASIN_REINJECT");
        clear_env("FLEXAIDDS_BASIN_SIGMA_ANG");
    }
};

struct FakeGene {
    double to_ic = 0.0;
    int32_t to_int32 = 0;
};
struct FakeLim {
    double min = 0.0;
    double max = 10.0;
    double nbin = 11.0;
};

}  // namespace

TEST_F(NewSearchArchEnv, DefaultOff) {
    EXPECT_FALSE(flexaids::new_search::phenotype_unique_enabled());
    EXPECT_FALSE(flexaids::new_search::basin_reinject_enabled());
}

TEST_F(NewSearchArchEnv, SeparateEnvGates) {
    set_env("FLEXAIDDS_PHENOTYPE_UNIQUE", "1");
    EXPECT_TRUE(flexaids::new_search::phenotype_unique_enabled());
    EXPECT_FALSE(flexaids::new_search::basin_reinject_enabled());
    clear_env("FLEXAIDDS_PHENOTYPE_UNIQUE");
    set_env("FLEXAIDDS_BASIN_REINJECT", "1");
    EXPECT_FALSE(flexaids::new_search::phenotype_unique_enabled());
    EXPECT_TRUE(flexaids::new_search::basin_reinject_enabled());
}

TEST_F(NewSearchArchEnv, NewSearchArchOneEnablesBoth) {
    set_env("FLEXAIDDS_NEW_SEARCH_ARCH", "1");
    EXPECT_TRUE(flexaids::new_search::phenotype_unique_enabled());
    EXPECT_TRUE(flexaids::new_search::basin_reinject_enabled());
}

TEST_F(NewSearchArchEnv, NewSearchArchTokenList) {
    set_env("FLEXAIDDS_NEW_SEARCH_ARCH", "phenotype_unique,basin_reinject");
    EXPECT_TRUE(flexaids::new_search::phenotype_unique_enabled());
    EXPECT_TRUE(flexaids::new_search::basin_reinject_enabled());
    set_env("FLEXAIDDS_NEW_SEARCH_ARCH", "phenotype_unique");
    EXPECT_TRUE(flexaids::new_search::phenotype_unique_enabled());
    EXPECT_FALSE(flexaids::new_search::basin_reinject_enabled());
}

TEST_F(NewSearchArchEnv, PhenotypeBinAndHashStable) {
    FakeGene g[2];
    FakeLim lim[2];
    lim[0] = {0.0, 10.0, 11.0};
    lim[1] = {0.0, 10.0, 11.0};
    g[0].to_ic = 5.0;
    g[1].to_ic = 2.0;
    const auto h1 = flexaids::new_search::hash_phenotype_bins(g, 2, lim);
    g[0].to_ic = 5.01;  // same bin for nbin=11 over [0,10]
    const auto h2 = flexaids::new_search::hash_phenotype_bins(g, 2, lim);
    EXPECT_EQ(h1, h2);
    g[0].to_ic = 9.5;  // different bin
    const auto h3 = flexaids::new_search::hash_phenotype_bins(g, 2, lim);
    EXPECT_NE(h1, h3);
}

TEST_F(NewSearchArchEnv, ApplyPhenotypeBinStepChangesInt) {
    FakeGene g;
    g.to_int32 = 1000000000;
    const bool changed =
        flexaids::new_search::apply_phenotype_bin_step(&g, /*nbin=*/11.0, /*sign=*/1,
                                                       /*k_bins=*/1);
    EXPECT_TRUE(changed);
    EXPECT_GT(g.to_int32, 1000000000);
}

TEST_F(NewSearchArchEnv, OutsideBasinAndSigma) {
    EXPECT_TRUE(flexaids::new_search::outside_basin(2.5, 2.0));
    EXPECT_FALSE(flexaids::new_search::outside_basin(1.5, 2.0));
    EXPECT_FALSE(flexaids::new_search::outside_basin(2.0, 2.0));
    EXPECT_DOUBLE_EQ(flexaids::new_search::basin_sigma_ang(2.0), 2.0);
    set_env("FLEXAIDDS_BASIN_SIGMA_ANG", "3.5");
    EXPECT_DOUBLE_EQ(flexaids::new_search::basin_sigma_ang(2.0), 3.5);
}
