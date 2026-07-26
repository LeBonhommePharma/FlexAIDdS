// Unit tests for G4.2 niche pair distance (drives shipped LIB/niche_distance.h).
#include "niche_distance.h"

#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <vector>

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

class NicheDistanceEnv : public ::testing::Test {
protected:
    void SetUp() override {
        clear_env("FLEXAIDDS_NICHE_CARTESIAN");
        clear_env("FLEXAIDDS_NICHE_SIGMA_ANG");
    }
    void TearDown() override {
        clear_env("FLEXAIDDS_NICHE_CARTESIAN");
        clear_env("FLEXAIDDS_NICHE_SIGMA_ANG");
    }
};

}  // namespace

TEST_F(NicheDistanceEnv, CartesianEnvDefaultOff) {
    EXPECT_FALSE(flexaids::niche_cartesian_env_enabled());
}

TEST_F(NicheDistanceEnv, CartesianEnvOn) {
    set_env("FLEXAIDDS_NICHE_CARTESIAN", "1");
    EXPECT_TRUE(flexaids::niche_cartesian_env_enabled());
}

TEST_F(NicheDistanceEnv, SigmaAngDefaultAndOverride) {
    EXPECT_DOUBLE_EQ(flexaids::niche_cartesian_sigma_ang(2.0), 2.0);
    set_env("FLEXAIDDS_NICHE_SIGMA_ANG", "3.5");
    EXPECT_DOUBLE_EQ(flexaids::niche_cartesian_sigma_ang(2.0), 3.5);
}

// Gene-space RMSP: pure to_ic RMS — gene0 ordinal mixed with angles (defect).
TEST(NicheGeneRmsp, IdenticalZero) {
    const double a[3] = {0.0, 10.0, -5.0};
    const double b[3] = {0.0, 10.0, -5.0};
    EXPECT_NEAR(flexaids::niche_gene_rmsp(a, b, 3), 0.0, 1e-12);
}

TEST(NicheGeneRmsp, AllAngleFlipMatchesPhase4Scale) {
    // PHASE4: 9 angles × 180° + gene0 fixed → rmsp ≈ 170.8 < 204 sigma
    std::vector<double> base(10, 0.0);
    base[0] = 1000.0;  // ordinal
    std::vector<double> flipped = base;
    for (int i = 1; i < 10; ++i) {
        flipped[static_cast<size_t>(i)] = 180.0;
    }
    const double d = flexaids::niche_gene_rmsp(base.data(), flipped.data(), 10);
    EXPECT_NEAR(d, std::sqrt(9.0 * 180.0 * 180.0 / 10.0), 1e-9);
    EXPECT_LT(d, 204.19);
}

TEST(NicheGeneRmsp, LargeOrdinalJumpExceedsAngleFlip) {
    std::vector<double> base(10, 0.0);
    base[0] = 1000.0;
    std::vector<double> ord = base;
    ord[0] = 1000.0 + 5000.0;
    std::vector<double> ang = base;
    for (int i = 1; i < 10; ++i) {
        ang[static_cast<size_t>(i)] = 180.0;
    }
    const double d_ord = flexaids::niche_gene_rmsp(base.data(), ord.data(), 10);
    const double d_ang = flexaids::niche_gene_rmsp(base.data(), ang.data(), 10);
    EXPECT_GT(d_ord, d_ang);
}

// Cartesian path: pure XYZ RMSD — NO gene/ordinal input.
TEST(NicheCartesianRmsd, IdenticalZero) {
    const float a[6] = {0.f, 0.f, 0.f, 1.f, 0.f, 0.f};
    const float b[6] = {0.f, 0.f, 0.f, 1.f, 0.f, 0.f};
    EXPECT_NEAR(flexaids::niche_cartesian_rmsd(a, b, 2), 0.0, 1e-12);
}

TEST(NicheCartesianRmsd, KnownShift) {
    // One atom shifted by 3 Å in x → RMSD = 3
    const float a[3] = {0.f, 0.f, 0.f};
    const float b[3] = {3.f, 0.f, 0.f};
    EXPECT_NEAR(flexaids::niche_cartesian_rmsd(a, b, 1), 3.0, 1e-6);
}

TEST(NicheCartesianRmsd, TwoAtoms) {
    // atom0: (0,0,0)->(0,0,0); atom1: (0,0,0)->(0,4,0) → RMSD = sqrt(16/2)=sqrt(8)
    const float a[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
    const float b[6] = {0.f, 0.f, 0.f, 0.f, 4.f, 0.f};
    EXPECT_NEAR(flexaids::niche_cartesian_rmsd(a, b, 2), std::sqrt(8.0), 1e-6);
}

// Dispatch: use_cartesian chooses path; cart ignores gene ordinal.
TEST(NichePairDispatch, CartesianIgnoresGeneOrdinal) {
    // Gene vectors differ only in ordinal gene0 (huge)
    const double ic_a[3] = {0.0, 0.0, 0.0};
    const double ic_b[3] = {10000.0, 0.0, 0.0};
    // XYZ identical → cart distance 0 regardless of genes
    const float xyz[3] = {1.f, 2.f, 3.f};
    const double d_gene =
        flexaids::niche_pair_distance(false, ic_a, ic_b, 3, xyz, xyz, 1);
    const double d_cart =
        flexaids::niche_pair_distance(true, ic_a, ic_b, 3, xyz, xyz, 1);
    EXPECT_GT(d_gene, 100.0);  // ordinal dominates gene path
    EXPECT_NEAR(d_cart, 0.0, 1e-12);  // cart path never sees gene0
}

TEST(NichePairDispatch, GenePathMatchesNicheGeneRmsp) {
    const double a[4] = {1.0, 2.0, 3.0, 4.0};
    const double b[4] = {1.0, 2.0, 3.0, 7.0};
    const float dummy[3] = {0.f, 0.f, 0.f};
    const double d = flexaids::niche_pair_distance(false, a, b, 4, dummy, dummy, 1);
    EXPECT_NEAR(d, flexaids::niche_gene_rmsp(a, b, 4), 1e-12);
}
