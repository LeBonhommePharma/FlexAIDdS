// tests/test_pose_helix_thermo_rewrite.cpp
// Experimental 2014 NATURAL pose→helix ΔH/ΔS rewrite (not claim-ready).
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "NATURaL/PoseHelixThermoRewrite.h"

#include <cmath>
#include <vector>

using namespace natural;

static constexpr double kTol = 1e-9;

static HelixSegment decision_helix() {
    HelixSegment h;
    h.stem_i_begin = 0;
    h.stem_i_end   = 1;
    h.loop_begin   = 2;
    h.loop_end     = 3;
    h.stem_j_begin = 4;
    h.stem_j_end   = 5;
    h.label        = "decision";
    return h;
}

static ReceptorNtCoord nt_at(int idx, double x, double y, double z,
                              bool hbond, bool aromatic) {
    ReceptorNtCoord n;
    n.nt_index = idx;
    n.x = x; n.y = y; n.z = z;
    n.is_hbond_partner = hbond;
    n.is_aromatic = aromatic;
    n.nx = 0.0; n.ny = 0.0; n.nz = 1.0;
    return n;
}

static LigandAtomCoord atom_at(double x, double y, double z,
                               bool hbond, bool aromatic) {
    LigandAtomCoord a;
    a.x = x; a.y = y; a.z = z;
    a.is_hbond_partner = hbond;
    a.is_aromatic = aromatic;
    a.nx = 0.0; a.ny = 0.0; a.nz = 1.0;
    return a;
}

static void expect_finite_result(const PoseHelixRewriteResult& r) {
    EXPECT_TRUE(r.applied);
    for (const auto& p : r.patches) {
        EXPECT_TRUE(std::isfinite(p.delta_H_kcal));
        EXPECT_TRUE(std::isfinite(p.delta_S_cal_per_K));
        EXPECT_TRUE(std::isfinite(p.delta_G_kcal));
        EXPECT_TRUE(std::isfinite(p.pose_boltzmann_w));
        EXPECT_FALSE(std::isnan(p.delta_H_kcal));
        EXPECT_FALSE(std::isnan(p.delta_S_cal_per_K));
        EXPECT_FALSE(std::isnan(p.delta_G_kcal));
    }
    for (const auto& p : r.ensemble) {
        EXPECT_TRUE(std::isfinite(p.delta_H_kcal));
        EXPECT_TRUE(std::isfinite(p.delta_S_cal_per_K));
        EXPECT_TRUE(std::isfinite(p.delta_G_kcal));
        EXPECT_FALSE(std::isnan(p.delta_G_kcal));
    }
}

TEST(PoseHelixThermoRewrite, StackingOnStemGivesNegativeDeltaH) {
    PoseHelixRewriteConfig cfg;
    LigandPose pose;
    pose.rank = 0;
    pose.score_kcal = -1.0;
    pose.atoms.push_back(atom_at(0.0, 0.0, 3.5, /*hbond=*/false, /*aromatic=*/true));

    std::vector<ReceptorNtCoord> rna = {
        nt_at(0, 0.0, 0.0, 0.0, /*hbond=*/false, /*aromatic=*/true),
    };

    auto r = rewrite_helices_from_poses({decision_helix()}, rna, {pose}, cfg);
    ASSERT_FALSE(r.patches.empty());
    EXPECT_LT(r.patches[0].delta_H_kcal, 0.0);
    EXPECT_NEAR(r.patches[0].delta_H_kcal, cfg.dH_per_stack_kcal, kTol);
    EXPECT_NEAR(r.patches[0].delta_S_cal_per_K, 0.0, kTol);
    expect_finite_result(r);
}

TEST(PoseHelixThermoRewrite, HBondIntoLoopGivesNegativeDeltaS) {
    PoseHelixRewriteConfig cfg;
    LigandPose pose;
    pose.rank = 0;
    pose.score_kcal = -1.0;
    pose.atoms.push_back(atom_at(10.0, 0.0, 3.0, /*hbond=*/true, /*aromatic=*/false));

    std::vector<ReceptorNtCoord> rna = {
        nt_at(2, 10.0, 0.0, 0.0, /*hbond=*/true, /*aromatic=*/false),
    };

    auto r = rewrite_helices_from_poses({decision_helix()}, rna, {pose}, cfg);
    ASSERT_FALSE(r.patches.empty());
    EXPECT_LT(r.patches[0].delta_S_cal_per_K, 0.0);
    EXPECT_NEAR(r.patches[0].delta_S_cal_per_K, cfg.dS_per_hbond_cal_K, kTol);
    EXPECT_NEAR(r.patches[0].delta_H_kcal, 0.0, kTol);
    const double expect_dG = helix_thermo_delta_g_kcal(
        0.0, cfg.dS_per_hbond_cal_K, cfg.T_K);
    EXPECT_NEAR(r.patches[0].delta_G_kcal, expect_dG, kTol);
    expect_finite_result(r);
}

TEST(PoseHelixThermoRewrite, ContactOutsideDecisionHelixNoPatchWhenRequired) {
    PoseHelixRewriteConfig cfg;
    cfg.require_decision_helix_only = true;

    LigandPose pose;
    pose.rank = 0;
    pose.score_kcal = -2.0;
    pose.atoms.push_back(atom_at(50.0, 0.0, 3.0, /*hbond=*/true, /*aromatic=*/true));

    std::vector<ReceptorNtCoord> rna = {
        nt_at(20, 50.0, 0.0, 0.0, /*hbond=*/true, /*aromatic=*/true),
    };

    auto r = rewrite_helices_from_poses({decision_helix()}, rna, {pose}, cfg);
    EXPECT_TRUE(r.patches.empty());
    EXPECT_TRUE(r.ensemble.empty());
    expect_finite_result(r);
}

TEST(PoseHelixThermoRewrite, ContactOutsideDecisionHelixPatchesWhenNotRequired) {
    PoseHelixRewriteConfig cfg;
    cfg.require_decision_helix_only = false;

    LigandPose pose;
    pose.rank = 0;
    pose.score_kcal = -2.0;
    pose.atoms.push_back(atom_at(50.0, 0.0, 3.0, /*hbond=*/true, /*aromatic=*/false));

    std::vector<ReceptorNtCoord> rna = {
        nt_at(20, 50.0, 0.0, 0.0, /*hbond=*/true, /*aromatic=*/false),
    };

    auto r = rewrite_helices_from_poses({decision_helix()}, rna, {pose}, cfg);
    ASSERT_EQ(r.patches.size(), 1u);
    EXPECT_EQ(r.patches[0].helix_label, "non_decision");
    EXPECT_LT(r.patches[0].delta_S_cal_per_K, 0.0);
    expect_finite_result(r);
}

TEST(PoseHelixThermoRewrite, TwoPoseMixtureWeights) {
    PoseHelixRewriteConfig cfg;
    cfg.T_K = 298.15;
    cfg.top_k_poses = 5;

    LigandPose p0;
    p0.rank = 0;
    p0.score_kcal = -2.0;
    p0.atoms.push_back(atom_at(0.0, 0.0, 3.5, false, true));

    LigandPose p1;
    p1.rank = 1;
    p1.score_kcal = 0.0;
    p1.atoms.push_back(atom_at(0.0, 0.0, 3.5, false, true));

    std::vector<ReceptorNtCoord> rna = {
        nt_at(0, 0.0, 0.0, 0.0, false, true),
    };

    auto r = rewrite_helices_from_poses({decision_helix()}, rna, {p0, p1}, cfg);
    ASSERT_EQ(r.patches.size(), 2u);

    const auto w = boltzmann_weights_from_scores({-2.0, 0.0}, cfg.T_K);
    ASSERT_EQ(w.size(), 2u);
    EXPECT_NEAR(r.patches[0].pose_boltzmann_w, w[0], kTol);
    EXPECT_NEAR(r.patches[1].pose_boltzmann_w, w[1], kTol);
    EXPECT_GT(w[0], w[1]); // more negative CF score → higher weight
    EXPECT_NEAR(w[0] + w[1], 1.0, kTol);

    ASSERT_EQ(r.ensemble.size(), 1u);
    const double expect_dH = w[0] * cfg.dH_per_stack_kcal + w[1] * cfg.dH_per_stack_kcal;
    EXPECT_NEAR(r.ensemble[0].delta_H_kcal, expect_dH, 1e-12);
    expect_finite_result(r);
}

TEST(PoseHelixThermoRewrite, NoNaNOnEmptyAndMixed) {
    PoseHelixRewriteConfig cfg;
    auto empty = rewrite_helices_from_poses({}, {}, {}, cfg);
    EXPECT_TRUE(empty.patches.empty());
    EXPECT_TRUE(empty.ensemble.empty());
    EXPECT_TRUE(empty.applied);

    LigandPose pose;
    pose.rank = 0;
    pose.score_kcal = -0.5;
    pose.atoms.push_back(atom_at(0.0, 0.0, 3.5, false, true));
    pose.atoms.push_back(atom_at(10.0, 0.0, 3.0, true, false));

    std::vector<ReceptorNtCoord> rna = {
        nt_at(0, 0.0, 0.0, 0.0, false, true),
        nt_at(2, 10.0, 0.0, 0.0, true, false),
    };

    auto r = rewrite_helices_from_poses({decision_helix()}, rna, {pose}, cfg);
    ASSERT_FALSE(r.patches.empty());
    EXPECT_LT(r.patches[0].delta_H_kcal, 0.0);
    EXPECT_LT(r.patches[0].delta_S_cal_per_K, 0.0);
    expect_finite_result(r);

    const double dG = helix_thermo_delta_g_kcal(
        r.patches[0].delta_H_kcal, r.patches[0].delta_S_cal_per_K, cfg.T_K);
    EXPECT_NEAR(r.patches[0].delta_G_kcal, dG, kTol);
}

TEST(PoseHelixThermoRewrite, EqualScoresGiveEqualWeights) {
    const auto w = boltzmann_weights_from_scores({-1.0, -1.0, -1.0}, 298.15);
    ASSERT_EQ(w.size(), 3u);
    EXPECT_NEAR(w[0], 1.0 / 3.0, kTol);
    EXPECT_NEAR(w[1], 1.0 / 3.0, kTol);
    EXPECT_NEAR(w[2], 1.0 / 3.0, kTol);
}

TEST(PoseHelixThermoRewrite, UnitsDeltaGMatchesDeltaHMinusTDeltaS) {
    const double dH = -1.2;
    const double dS = -3.0; // cal/K
    const double T = 298.15;
    const double dG = helix_thermo_delta_g_kcal(dH, dS, T);
    EXPECT_NEAR(dG, dH - T * (dS / 1000.0), 1e-12);
    EXPECT_TRUE(std::isfinite(dG));
}
