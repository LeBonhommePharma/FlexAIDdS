// tests/test_emission_guards.cpp — pathological emission/RMSD guards
//
// Drives LIB/emission_guards.h (shipped path used by DatasetRunner.cpp).
// Covers 1KZK (sentinel/empty), 1M2Z (clash CF), 1J3J/1IGJ (absurd RMSD).
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include "emission_guards.h"

#include <gtest/gtest.h>

#include <limits>
#include <string>
#include <vector>

TEST(EmissionGuards, SentinelRmsd)
{
    EXPECT_TRUE(flexaids::emission::is_sentinel_rmsd(-1.0f));
    EXPECT_TRUE(flexaids::emission::is_sentinel_rmsd(
        std::numeric_limits<float>::quiet_NaN()));
    EXPECT_FALSE(flexaids::emission::is_sentinel_rmsd(1.5f));
}

TEST(EmissionGuards, AbsurdRmsd)
{
    EXPECT_TRUE(flexaids::emission::is_absurd_rmsd(62.0f));
    EXPECT_TRUE(flexaids::emission::is_absurd_rmsd(25.01f));
    EXPECT_FALSE(flexaids::emission::is_absurd_rmsd(2.0f));
    EXPECT_FALSE(flexaids::emission::is_absurd_rmsd(24.9f));
}

TEST(EmissionGuards, ClashSentinelCf)
{
    EXPECT_TRUE(flexaids::emission::is_clash_sentinel_cf(10000.0));
    EXPECT_TRUE(flexaids::emission::is_clash_sentinel_cf(1e4));
    EXPECT_FALSE(flexaids::emission::is_clash_sentinel_cf(-40.0));
    EXPECT_FALSE(flexaids::emission::is_clash_sentinel_cf(9999.0));
}

TEST(EmissionGuards, ReportElectedRmsdBands)
{
    auto ok = flexaids::emission::report_elected_rmsd(1.2f);
    EXPECT_FALSE(ok.is_sentinel);
    EXPECT_FALSE(ok.is_absurd);
    EXPECT_TRUE(ok.usable_for_success);
    EXPECT_FLOAT_EQ(ok.rmsd, 1.2f);

    auto sent = flexaids::emission::report_elected_rmsd(-1.0f);
    EXPECT_TRUE(sent.is_sentinel);
    EXPECT_FALSE(sent.usable_for_success);
    EXPECT_FLOAT_EQ(sent.rmsd, flexaids::emission::kSentinelRmsd);

    auto absd = flexaids::emission::report_elected_rmsd(74.0f);
    EXPECT_TRUE(absd.is_absurd);
    EXPECT_FALSE(absd.usable_for_success);
    EXPECT_FLOAT_EQ(absd.rmsd, 74.0f);
}

TEST(EmissionGuards, RecoverPoseCountFromHeads)
{
    // 1KZK-class: root count 0, restart heads exist
    EXPECT_EQ(flexaids::emission::recover_pose_count(0, 50), 50);
    EXPECT_EQ(flexaids::emission::recover_pose_count(10, 5), 10);
    EXPECT_EQ(flexaids::emission::recover_pose_count(-1, 3), 3);
}

TEST(EmissionGuards, ChooseElectedPathFallback)
{
    EXPECT_EQ(flexaids::emission::choose_elected_path("a.pdb", {"b.pdb"}), "a.pdb");
    EXPECT_EQ(flexaids::emission::choose_elected_path("", {"", "b.pdb", "c.pdb"}),
              "b.pdb");
    EXPECT_TRUE(flexaids::emission::choose_elected_path("", {"", ""}).empty());
}

TEST(EmissionGuards, ElectionInputsOk)
{
    EXPECT_TRUE(flexaids::emission::election_inputs_ok(5, false, true));
    EXPECT_FALSE(flexaids::emission::election_inputs_ok(0, false, true));
    EXPECT_FALSE(flexaids::emission::election_inputs_ok(5, true, true));
    EXPECT_FALSE(flexaids::emission::election_inputs_ok(5, false, false));
}

TEST(EmissionGuards, PathologyTag)
{
    auto sent = flexaids::emission::report_elected_rmsd(-1.0f);
    EXPECT_STREQ(flexaids::emission::pathology_tag(sent, -10.0),
                 "empty_or_sentinel_rmsd");
    auto absd = flexaids::emission::report_elected_rmsd(60.0f);
    EXPECT_STREQ(flexaids::emission::pathology_tag(absd, -10.0), "absurd_rmsd");
    auto ok = flexaids::emission::report_elected_rmsd(1.0f);
    EXPECT_STREQ(flexaids::emission::pathology_tag(ok, 10000.0),
                 "clash_sentinel_cf");
    EXPECT_STREQ(flexaids::emission::pathology_tag(ok, -50.0), "ok");
}
