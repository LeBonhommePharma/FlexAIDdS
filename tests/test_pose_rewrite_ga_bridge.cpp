// tests/test_pose_rewrite_ga_bridge.cpp
// Opt-in DualAssembly / NATURAL GA-coord plumbing for PoseHelix / PoseLocal.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

#include "NATURaL/PoseRewriteGaBridge.h"

#include <cmath>
#include <vector>

using namespace natural;

namespace {

HelixSegment decision_helix()
{
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

ReceptorNtCoord stem_nt()
{
    ReceptorNtCoord n;
    n.nt_index = 0;
    n.x = 0.0;
    n.y = 0.0;
    n.z = 0.0;
    n.is_aromatic = true;
    n.nz = 1.0;
    return n;
}

LigandPose stacking_pose()
{
    LigandAtomCoord a;
    a.x = 0.0;
    a.y = 0.0;
    a.z = 3.5;
    a.is_aromatic = true;
    a.nz = 1.0;
    LigandPose p;
    p.rank = 0;
    p.score_kcal = -1.0;
    p.atoms.push_back(a);
    return p;
}

DecisionElement rna_hp()
{
    DecisionElement el;
    el.ss = SSClass::RNA_STEM_LOOP;
    el.label = "hp";
    el.na_start = 0;
    el.na_end = 10;
    return el;
}

} // namespace

TEST(PoseRewriteGaBridge, CoordsPresentRequiresAtoms)
{
    EXPECT_FALSE(ga_atom_coords_present({}, {}));
    EXPECT_FALSE(ga_atom_coords_present({stem_nt()}, {}));
    LigandPose empty_pose;
    empty_pose.score_kcal = -1.0;
    EXPECT_FALSE(ga_atom_coords_present({stem_nt()}, {empty_pose}));
    EXPECT_TRUE(ga_atom_coords_present({stem_nt()}, {stacking_pose()}));
}

TEST(PoseRewriteGaBridge, HelixDefaultOffDoesNotRewrite)
{
    const auto r = pose_helix_from_ga_or_closed(
        /*enabled=*/false,
        {decision_helix()},
        {stem_nt()},
        {stacking_pose()});
    EXPECT_FALSE(r.applied);
    EXPECT_TRUE(r.patches.empty());
    EXPECT_TRUE(r.ensemble.empty());
}

TEST(PoseRewriteGaBridge, HelixOptInWithCoordsCallsRewrite)
{
    const auto r = pose_helix_from_ga_or_closed(
        /*enabled=*/true,
        {decision_helix()},
        {stem_nt()},
        {stacking_pose()});
    EXPECT_TRUE(r.applied);
    ASSERT_FALSE(r.ensemble.empty());
    EXPECT_LT(r.ensemble.front().delta_H_kcal, 0.0);
    EXPECT_TRUE(std::isfinite(r.ensemble.front().delta_G_kcal));
}

TEST(PoseRewriteGaBridge, HelixMissingCoordsFailClosed)
{
    const auto r = pose_helix_from_ga_or_closed(
        /*enabled=*/true,
        {decision_helix()},
        {},
        {});
    EXPECT_FALSE(r.applied);
    EXPECT_TRUE(r.patches.empty());
    EXPECT_TRUE(r.ensemble.empty());
}

TEST(PoseRewriteGaBridge, HelixMissingHelicesFailClosed)
{
    const auto r = pose_helix_from_ga_or_closed(
        /*enabled=*/true,
        {},
        {stem_nt()},
        {stacking_pose()});
    EXPECT_FALSE(r.applied);
}

TEST(PoseRewriteGaBridge, LocalDefaultOffDoesNotRewrite)
{
    PoseThermoRewriteConfig cfg = default_pose_thermo_rewrite_config();
    const auto mix = pose_local_from_ga_or_closed(
        /*enabled=*/false,
        {rna_hp()},
        {},
        {stem_nt()},
        {stacking_pose()},
        cfg);
    EXPECT_TRUE(mix.empty());
    EXPECT_EQ(mix.n_poses_used, 0);
}

TEST(PoseRewriteGaBridge, LocalOptInWithCoordsCallsRewrite)
{
    PoseThermoRewriteConfig cfg = default_pose_thermo_rewrite_config();
    const auto mix = pose_local_from_ga_or_closed(
        /*enabled=*/true,
        {rna_hp()},
        {},
        {stem_nt()},
        {stacking_pose()},
        cfg);
    EXPECT_FALSE(mix.empty());
    EXPECT_LT(mix.dH_kcal, 0.0);
    EXPECT_NEAR(mix.dH_kcal, kXia1998RnaWcStackMean_dH_kcal, 1e-12);
    EXPECT_TRUE(std::isfinite(mix.dG_kcal));
}

TEST(PoseRewriteGaBridge, LocalMissingCoordsFailClosed)
{
    PoseThermoRewriteConfig cfg = default_pose_thermo_rewrite_config();
    const auto mix = pose_local_from_ga_or_closed(
        /*enabled=*/true,
        {rna_hp()},
        {},
        {},
        {},
        cfg);
    EXPECT_TRUE(mix.empty());
}

TEST(PoseRewriteGaBridge, LocalLabelledFallbackWithoutCoords)
{
    PoseThermoRewriteConfig cfg = default_pose_thermo_rewrite_config();
    PoseView labelled;
    labelled.score_kcal = -1.0;
    PoseContact c;
    c.kind = ContactKind::RNA_STEM_STACK;
    c.geometry_weight = 1.0;
    c.nt_index = 3;
    labelled.contacts.push_back(c);

    const auto mix = pose_local_from_ga_or_closed(
        /*enabled=*/true,
        {rna_hp()},
        {labelled},
        {},
        {},
        cfg);
    EXPECT_FALSE(mix.empty());
    EXPECT_NEAR(mix.dH_kcal, kXia1998RnaWcStackMean_dH_kcal, 1e-12);
}

TEST(PoseRewriteGaBridge, ConverterEmitsRnaStemStack)
{
    const auto views = pose_views_from_ga_atom_coords(
        {stem_nt()}, {stacking_pose()}, {rna_hp()});
    ASSERT_EQ(views.size(), 1u);
    ASSERT_EQ(views[0].contacts.size(), 1u);
    EXPECT_EQ(views[0].contacts[0].kind, ContactKind::RNA_STEM_STACK);
    EXPECT_EQ(views[0].contacts[0].nt_index, 0);
}
