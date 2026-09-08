// tests/test_pose_local_thermo_rewrite.cpp
// Experimental PoseLocalThermoRewrite: class-specific ΔH/ΔS tables share algebra
// only. RNA Turner-shaped ≠ DNA SantaLucia-shaped ≠ protein helix ≠ sheet.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include <gtest/gtest.h>

#include "NATURaL/PoseLocalThermoRewrite.h"

#include <cmath>
#include <string>
#include <vector>

using natural::ContactKind;
using natural::DecisionElement;
using natural::PoseContact;
using natural::PoseView;
using natural::SSClass;
using natural::delta_G_kcal;
using natural::kCalPerKcal;
using natural::kPoseRewrite_kB_kcal_mol_K;
using natural::rewrite_local_thermo_from_poses;
using natural::ss_class_for_kind;
using natural::ss_class_name;
using natural::table_for_ss;

namespace {

PoseView one_contact(ContactKind kind, double w, int nt, int res, int partner = -1)
{
    PoseView p;
    p.score_kcal = 0.0;
    PoseContact c;
    c.kind = kind;
    c.geometry_weight = w;
    c.nt_index = nt;
    c.res_index = res;
    c.partner_res_index = partner;
    p.contacts.push_back(c);
    return p;
}

DecisionElement rna_loop(const std::string& label, int lo, int hi)
{
    DecisionElement el;
    el.ss = SSClass::RNA_STEM_LOOP;
    el.label = label;
    el.na_start = lo;
    el.na_end = hi;
    return el;
}

DecisionElement dna_loop(const std::string& label, int lo, int hi)
{
    DecisionElement el;
    el.ss = SSClass::DNA_STEM_LOOP;
    el.label = label;
    el.na_start = lo;
    el.na_end = hi;
    return el;
}

DecisionElement helix(const std::string& label, int lo, int hi)
{
    DecisionElement el;
    el.ss = SSClass::PROTEIN_HELIX;
    el.label = label;
    el.res_start = lo;
    el.res_end = hi;
    return el;
}

DecisionElement sheet(const std::string& label, int lo, int hi, int plo, int phi)
{
    DecisionElement el;
    el.ss = SSClass::PROTEIN_SHEET;
    el.label = label;
    el.res_start = lo;
    el.res_end = hi;
    el.sheet_partner_start = plo;
    el.sheet_partner_end = phi;
    return el;
}

} // namespace

TEST(PoseLocalThermoRewrite, TablesAreClassSpecific)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const auto& rna = table_for_ss(cfg, SSClass::RNA_STEM_LOOP);
    const auto& dna = table_for_ss(cfg, SSClass::DNA_STEM_LOOP);
    const auto& hx  = table_for_ss(cfg, SSClass::PROTEIN_HELIX);
    const auto& sh  = table_for_ss(cfg, SSClass::PROTEIN_SHEET);

    const int rna_stack = static_cast<int>(ContactKind::RNA_STEM_STACK);
    const int dna_stack = static_cast<int>(ContactKind::DNA_STEM_STACK);
    EXPECT_NEAR(rna.increments[rna_stack].dH_kcal, -1.2, 1e-12);
    EXPECT_NEAR(dna.increments[dna_stack].dH_kcal, -1.0, 1e-12);
    EXPECT_NE(rna.increments[rna_stack].dH_kcal, dna.increments[dna_stack].dH_kcal);
    EXPECT_NEAR(dna.increments[rna_stack].dH_kcal, 0.0, 1e-12);
    EXPECT_NEAR(rna.increments[dna_stack].dH_kcal, 0.0, 1e-12);

    const int hx_hb = static_cast<int>(ContactKind::PROTEIN_HELIX_BB_HBOND);
    const int sh_hb = static_cast<int>(ContactKind::PROTEIN_SHEET_BRIDGE_HBOND);
    EXPECT_NEAR(hx.increments[hx_hb].dH_kcal, -1.5, 1e-12);
    EXPECT_NEAR(sh.increments[sh_hb].dH_kcal, -1.8, 1e-12);
    EXPECT_NE(hx.increments[hx_hb].dH_kcal, sh.increments[sh_hb].dH_kcal);
    EXPECT_NE(hx.increments[hx_hb].dS_cal_per_mol_K, sh.increments[sh_hb].dS_cal_per_mol_K);
}

TEST(PoseLocalThermoRewrite, RnaStemStackDeltaHNegative)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const auto mix = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::RNA_STEM_STACK, 1.0, /*nt=*/3, /*res=*/-1)},
        {rna_loop("hp", 0, 10)},
        cfg);
    ASSERT_EQ(mix.n_poses_used, 1);
    EXPECT_LT(mix.dH_kcal, 0.0);
    EXPECT_NEAR(mix.dH_kcal, -1.2, 1e-12);
    EXPECT_NEAR(mix.dS_cal_per_mol_K, 0.0, 1e-12);
    EXPECT_EQ(ss_class_name(mix.per_element.at(0).ss), std::string("RNA_STEM_LOOP"));
}

TEST(PoseLocalThermoRewrite, RnaLoopHbondDeltaSNegative)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const auto mix = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::RNA_LOOP_HBOND, 1.0, /*nt=*/5, /*res=*/-1)},
        {rna_loop("hp", 0, 10)},
        cfg);
    ASSERT_EQ(mix.n_poses_used, 1);
    EXPECT_LT(mix.dS_cal_per_mol_K, 0.0);
    EXPECT_NEAR(mix.dS_cal_per_mol_K, -3.0, 1e-12);
    EXPECT_NEAR(mix.dH_kcal, 0.0, 1e-12);
}

TEST(PoseLocalThermoRewrite, DnaStemStackUsesDnaTableNotRna)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    auto rna_pose = one_contact(ContactKind::RNA_STEM_STACK, 1.0, 3, -1);
    auto dna_pose = one_contact(ContactKind::DNA_STEM_STACK, 1.0, 3, -1);
    const auto rna = rewrite_local_thermo_from_poses(
        {rna_pose}, {rna_loop("r", 0, 10)}, cfg);
    const auto dna = rewrite_local_thermo_from_poses(
        {dna_pose}, {dna_loop("d", 0, 10)}, cfg);

    ASSERT_EQ(dna.n_poses_used, 1);
    EXPECT_EQ(ss_class_name(dna.per_element.at(0).ss), std::string("DNA_STEM_LOOP"));
    EXPECT_LT(dna.dH_kcal, 0.0);
    EXPECT_NEAR(dna.dH_kcal, -1.0, 1e-12);
    EXPECT_NEAR(rna.dH_kcal, -1.2, 1e-12);
    EXPECT_NE(std::fabs(dna.dH_kcal), std::fabs(rna.dH_kcal));
}

TEST(PoseLocalThermoRewrite, ProteinHelixUsesHelixTable)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const auto mix = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::PROTEIN_HELIX_BB_HBOND, 1.0, -1, /*res=*/12)},
        {helix("H1", 10, 20)},
        cfg);
    ASSERT_EQ(mix.n_poses_used, 1);
    EXPECT_EQ(mix.per_element.at(0).ss, SSClass::PROTEIN_HELIX);
    EXPECT_NEAR(mix.dH_kcal, -1.5, 1e-12);
    EXPECT_NEAR(mix.dS_cal_per_mol_K, -1.0, 1e-12);
    EXPECT_NE(ss_class_for_kind(ContactKind::PROTEIN_HELIX_BB_HBOND), SSClass::RNA_STEM_LOOP);
    EXPECT_NE(ss_class_for_kind(ContactKind::PROTEIN_HELIX_BB_HBOND), SSClass::DNA_STEM_LOOP);
}

TEST(PoseLocalThermoRewrite, ProteinSheetBridgeDiffersFromHelix)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const auto hx = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::PROTEIN_HELIX_BB_HBOND, 1.0, -1, 12)},
        {helix("H1", 10, 20)},
        cfg);
    const auto sh = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::PROTEIN_SHEET_BRIDGE_HBOND, 1.0, -1, 12, /*partner=*/40)},
        {sheet("E1", 10, 20, 35, 45)},
        cfg);
    ASSERT_EQ(sh.n_poses_used, 1);
    EXPECT_EQ(sh.per_element.at(0).ss, SSClass::PROTEIN_SHEET);
    EXPECT_NEAR(sh.dH_kcal, -1.8, 1e-12);
    EXPECT_NEAR(sh.dS_cal_per_mol_K, -1.5, 1e-12);
    EXPECT_NE(sh.dH_kcal, hx.dH_kcal);
    EXPECT_NE(sh.dS_cal_per_mol_K, hx.dS_cal_per_mol_K);
}

TEST(PoseLocalThermoRewrite, OutsideDecisionElementEmptyWhenRequired)
{
    auto cfg = natural::default_pose_thermo_rewrite_config();
    cfg.require_decision_element_only = true;
    const auto mix = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::RNA_STEM_STACK, 1.0, /*nt=*/99, -1)},
        {rna_loop("hp", 0, 10)},
        cfg);
    EXPECT_TRUE(mix.empty());
    EXPECT_EQ(mix.n_poses_used, 0);
    EXPECT_TRUE(mix.per_element.empty());
    EXPECT_NEAR(mix.dH_kcal, 0.0, 1e-12);
    EXPECT_NEAR(mix.dG_kcal, 0.0, 1e-12);
}

TEST(PoseLocalThermoRewrite, TwoPoseBoltzmannMixture)
{
    auto cfg = natural::default_pose_thermo_rewrite_config();
    cfg.temperature_K = 298.15;
    cfg.top_k = 2;
    PoseView a = one_contact(ContactKind::RNA_STEM_STACK, 1.0, 3, -1);
    PoseView b = one_contact(ContactKind::RNA_STEM_STACK, 2.0, 4, -1);
    a.score_kcal = -1.0;
    b.score_kcal = -0.5;
    const auto mix = rewrite_local_thermo_from_poses(
        {a, b}, {rna_loop("hp", 0, 10)}, cfg);

    ASSERT_EQ(mix.n_poses_used, 2);
    ASSERT_EQ(mix.pose_weights.size(), 2u);

    const double dH0 = -1.2;
    const double dH1 = -2.4;
    const double dG0 = delta_G_kcal(dH0, 0.0, cfg.temperature_K);
    const double dG1 = delta_G_kcal(dH1, 0.0, cfg.temperature_K);
    const double kT = kPoseRewrite_kB_kcal_mol_K * cfg.temperature_K;
    const double lw0 = -dG0 / kT;
    const double lw1 = -dG1 / kT;
    const double m = std::max(lw0, lw1);
    const double Z = std::exp(lw0 - m) + std::exp(lw1 - m);
    const double p0 = std::exp(lw0 - m) / Z;
    const double p1 = std::exp(lw1 - m) / Z;
    EXPECT_NEAR(mix.pose_weights[0], p0, 1e-10);
    EXPECT_NEAR(mix.pose_weights[1], p1, 1e-10);
    EXPECT_NEAR(mix.dH_kcal, p0 * dH0 + p1 * dH1, 1e-10);
    EXPECT_NEAR(mix.pose_weights[0] + mix.pose_weights[1], 1.0, 1e-12);
}

TEST(PoseLocalThermoRewrite, NoNaNAndDocumentedUnits)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const std::vector<PoseView> poses = {
        one_contact(ContactKind::RNA_STEM_STACK, 1.0, 2, -1),
        one_contact(ContactKind::DNA_STEM_STACK, 1.0, 2, -1),
        one_contact(ContactKind::PROTEIN_HELIX_PACKING, 1.0, -1, 11),
        one_contact(ContactKind::PROTEIN_SHEET_SHEAR, 1.0, -1, 12, 40),
    };
    const std::vector<DecisionElement> els = {
        rna_loop("r", 0, 10),
        dna_loop("d", 0, 10),
        helix("h", 0, 20),
        sheet("e", 0, 20, 30, 50),
    };
    const auto mix = rewrite_local_thermo_from_poses(poses, els, cfg);
    EXPECT_TRUE(std::isfinite(mix.dH_kcal));
    EXPECT_TRUE(std::isfinite(mix.dS_cal_per_mol_K));
    EXPECT_TRUE(std::isfinite(mix.dG_kcal));
    for (double w : mix.pose_weights)
        EXPECT_TRUE(std::isfinite(w));
    for (const auto& p : mix.per_element) {
        EXPECT_TRUE(std::isfinite(p.dH_kcal));
        EXPECT_TRUE(std::isfinite(p.dS_cal_per_mol_K));
        EXPECT_TRUE(std::isfinite(p.dG_kcal));
        EXPECT_NEAR(p.dG_kcal,
                    p.dH_kcal - cfg.temperature_K * (p.dS_cal_per_mol_K / kCalPerKcal),
                    1e-12);
    }
    EXPECT_NEAR(mix.dG_kcal,
                mix.dH_kcal - cfg.temperature_K * (mix.dS_cal_per_mol_K / kCalPerKcal),
                1e-12);
}

TEST(PoseLocalThermoRewrite, ProteinOtherDoesNotRewriteByDefault)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    DecisionElement other;
    other.ss = SSClass::PROTEIN_OTHER;
    other.label = "loop";
    other.res_start = 0;
    other.res_end = 50;
    PoseView pose = one_contact(ContactKind::PROTEIN_HELIX_BB_HBOND, 1.0, -1, 10);
    const auto mix = rewrite_local_thermo_from_poses({pose}, {other}, cfg);
    EXPECT_TRUE(mix.empty());
}

TEST(PoseLocalThermoRewrite, RnaKindOnDnaElementIsIgnored)
{
    const auto cfg = natural::default_pose_thermo_rewrite_config();
    const auto mix = rewrite_local_thermo_from_poses(
        {one_contact(ContactKind::RNA_STEM_STACK, 1.0, 3, -1)},
        {dna_loop("d", 0, 10)},
        cfg);
    EXPECT_TRUE(mix.empty());
}
