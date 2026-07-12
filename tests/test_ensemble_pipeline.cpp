// tests/test_ensemble_pipeline.cpp
// Reproducibility contract for the 4-layer ensemble pipeline (pure helpers).
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/ensemble_pipeline.h"

#include <cmath>
#include <vector>

using ensemble::FrameChartStatus;

TEST(FrameChart, ClassifyOkWarnFail) {
    EXPECT_EQ(ensemble::classify_frame_chart_rmsd(0.05, false), FrameChartStatus::Ok);
    EXPECT_EQ(ensemble::classify_frame_chart_rmsd(0.05, true), FrameChartStatus::Ok);
    EXPECT_EQ(ensemble::classify_frame_chart_rmsd(0.5, true), FrameChartStatus::Fail);
    EXPECT_EQ(ensemble::classify_frame_chart_rmsd(1.5, false), FrameChartStatus::Warn);
    EXPECT_EQ(ensemble::classify_frame_chart_rmsd(1.5, true), FrameChartStatus::Fail);
}

TEST(FrameChart, GeneChartResidualWithinBin) {
    const double min_ic = 0.0, max_ic = 10.0, del = 0.1;
    auto ictogene = [&](double ic) {
        if (ic < min_ic) ic = min_ic;
        if (ic > max_ic) ic = max_ic;
        return static_cast<int>(std::llround((ic - min_ic) / del));
    };
    auto genetoic = [&](int g) {
        return min_ic + static_cast<double>(g) * del;
    };
    const double resid = ensemble::gene_chart_max_residual(
        min_ic, max_ic, del, 21, genetoic, ictogene);
    EXPECT_LE(resid, del * 0.5 + 1e-9);
}

TEST(PocketSupport, LigandableScorePrefersCompactVolume) {
    // Larger compact cluster beats sparse thin one.
    const double s_dense = ensemble::ligandable_score(30, 2.0, 8.0);
    const double s_sparse = ensemble::ligandable_score(5, 1.0, 20.0);
    EXPECT_GT(s_dense, s_sparse);
}

TEST(PocketSupport, SelectTopKByScore) {
    std::vector<std::pair<int, double>> c = {{1, 10.0}, {2, 50.0}, {3, 30.0}, {4, 5.0}};
    auto keep = ensemble::select_top_k_clefts(c, 2);
    ASSERT_EQ(keep.size(), 2u);
    EXPECT_EQ(keep[0], 2);
    EXPECT_EQ(keep[1], 3);
}

TEST(PocketSupport, ValidSphereRadius) {
    EXPECT_TRUE(ensemble::valid_sphere_radius(1.5f));
    EXPECT_FALSE(ensemble::valid_sphere_radius(0.0f));
    EXPECT_FALSE(ensemble::valid_sphere_radius(0.1f));
    EXPECT_FALSE(ensemble::valid_sphere_radius(100.0f));
}

TEST(PocketSupport, CleftCentroidExtent) {
    sphere a{}, b{};
    a.center[0] = 0; a.center[1] = 0; a.center[2] = 0; a.radius = 1.0f;
    b.center[0] = 4; b.center[1] = 0; b.center[2] = 0; b.radius = 1.0f;
    a.prev = &b;
    b.prev = nullptr;
    ensemble::CleftCentroid g{};
    ASSERT_TRUE(ensemble::cleft_centroid_extent(&a, &g));
    EXPECT_NEAR(g.cx, 2.0, 1e-9);
    EXPECT_EQ(g.n_spheres, 2);
    EXPECT_NEAR(g.extent_A, 2.0 + 1.0, 1e-6);  // |2| + radius
}

TEST(SmfreeSoftBeta, SelectionBetaIsOneOverT) {
    double beta = 0.0;
    ASSERT_TRUE(ensemble::soft_selection_beta(300.0, &beta));
    EXPECT_NEAR(beta, 1.0 / 300.0, 1e-12);
    EXPECT_FALSE(ensemble::soft_selection_beta(0.0, &beta));
    EXPECT_FALSE(ensemble::smfree_requires_temperature(0.0));
    EXPECT_TRUE(ensemble::smfree_requires_temperature(300.0));
}

// Layer 4: classic entropy election (1HNN-class synthetic, reproducibility exhibit)
TEST(ClassicElection, OneHnnClassFlip) {
    // Pre-fix 1HNN: ACF basin freq 29 vs CF champion
    std::vector<double> acf = {-49.3, -83.4, -48.9, -263.4, -221.2};
    std::vector<double> cf  = {-189.9, -120.0, -100.0, -72.1, -80.0};
    EXPECT_EQ(ensemble::elect_rank0_index(acf, cf, false, 300), 3);
    EXPECT_EQ(ensemble::elect_rank0_index(acf, cf, true, 300), 0);
    EXPECT_EQ(ensemble::elect_rank0_index(acf, cf, false, 0), 0);
}
