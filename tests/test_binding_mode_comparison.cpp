// tests/test_binding_mode_comparison.cpp — Pose and BindingMode comparison operators
// Tests: PoseClassifier, Pose sorting, Pose construction
// Note: Pose::operator< is inline-defined in BindingMode.cpp and not linkable
//       from external TUs — PoseClassifier tests cover the same comparison logic.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <cstring>

#include "flexaid.h"
#include "gaboom.h"
#include "BindingMode.h"

// ═══════════════════════════════════════════════════════════════════════
// Helpers: create lightweight Pose objects for testing
// ═══════════════════════════════════════════════════════════════════════

namespace {

// Create a minimal chromosome with a given energy value
chromosome make_chrom(double evalue) {
    chromosome c;
    std::memset(&c, 0, sizeof(chromosome));
    c.cf.com = evalue;
    c.cf.wal = 0;
    c.cf.sas = 0;
    c.cf.con = 0;
    c.cf.elec = 0;
    c.app_evalue = evalue;
    c.fitnes = 0.0;
    c.genes = nullptr;
    return c;
}

// Pose constructor: Pose(chrom*, chrom_index, order, reachDist, temperature, vPose)
Pose make_pose(chromosome& c, int index, int order, float dist) {
    return Pose(&c, index, order, dist, 300, std::vector<float>{});
}

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════
// PoseClassifier tests
// ═══════════════════════════════════════════════════════════════════════

class PoseClassifierTest : public ::testing::Test {
protected:
    chromosome c1, c2;
    void SetUp() override {
        c1 = make_chrom(-10.0);
        c2 = make_chrom(-5.0);
    }
};

TEST_F(PoseClassifierTest, LowerOrderComesFirst) {
    Pose p1 = make_pose(c1, 0, 1, 1.0f);
    Pose p2 = make_pose(c2, 0, 2, 1.0f);
    PoseClassifier cmp;
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(PoseClassifierTest, SameOrder_LowerDistBreaks) {
    Pose p1 = make_pose(c1, 0, 1, 0.5f);
    Pose p2 = make_pose(c2, 0, 1, 1.5f);
    PoseClassifier cmp;
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(PoseClassifierTest, SameOrderSameDist_LowerIndexBreaks) {
    Pose p1 = make_pose(c1, 3, 1, 1.0f);
    Pose p2 = make_pose(c2, 7, 1, 1.0f);
    PoseClassifier cmp;
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(PoseClassifierTest, IdenticalPoses_NotLess) {
    Pose p1 = make_pose(c1, 5, 1, 2.0f);
    Pose p2 = make_pose(c2, 5, 1, 2.0f);
    PoseClassifier cmp;
    EXPECT_FALSE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(PoseClassifierTest, StrictWeakOrdering) {
    // Verify transitivity: a < b and b < c implies a < c
    Pose pa = make_pose(c1, 0, 1, 1.0f);
    Pose pb = make_pose(c1, 0, 2, 1.0f);
    Pose pc = make_pose(c1, 0, 3, 1.0f);
    PoseClassifier cmp;
    EXPECT_TRUE(cmp(pa, pb));
    EXPECT_TRUE(cmp(pb, pc));
    EXPECT_TRUE(cmp(pa, pc));
}

TEST_F(PoseClassifierTest, EnergyDoesNotAffectOrder) {
    // Energy values differ but comparison ignores them
    Pose p_low_e = make_pose(c1, 0, 1, 1.0f);  // evalue=-10
    Pose p_high_e = make_pose(c2, 0, 1, 1.0f); // evalue=-5
    PoseClassifier cmp;
    EXPECT_FALSE(cmp(p_low_e, p_high_e));
    EXPECT_FALSE(cmp(p_high_e, p_low_e));
}

// ═══════════════════════════════════════════════════════════════════════
// Sorting with PoseClassifier
// ═══════════════════════════════════════════════════════════════════════

TEST(PoseSortTest, SortByOrderThenDistThenIndex) {
    chromosome c = make_chrom(-5.0);
    std::vector<Pose> poses;
    // Insert in non-sorted order
    poses.push_back(make_pose(c, 5, 2, 1.0f));  // order=2
    poses.push_back(make_pose(c, 1, 1, 3.0f));  // order=1, dist=3
    poses.push_back(make_pose(c, 3, 1, 1.0f));  // order=1, dist=1
    poses.push_back(make_pose(c, 2, 1, 1.0f));  // order=1, dist=1, index=2

    std::sort(poses.begin(), poses.end(), PoseClassifier());

    // Expected: (order=1,dist=1,idx=2), (order=1,dist=1,idx=3), (order=1,dist=3,idx=1), (order=2,...)
    EXPECT_EQ(poses[0].chrom_index, 2);
    EXPECT_EQ(poses[1].chrom_index, 3);
    EXPECT_EQ(poses[2].chrom_index, 1);
    EXPECT_EQ(poses[3].chrom_index, 5);
}

TEST(PoseSortTest, AlreadySorted) {
    chromosome c = make_chrom(-5.0);
    std::vector<Pose> poses;
    poses.push_back(make_pose(c, 0, 1, 1.0f));
    poses.push_back(make_pose(c, 1, 1, 2.0f));
    poses.push_back(make_pose(c, 2, 2, 1.0f));

    std::sort(poses.begin(), poses.end(), PoseClassifier());

    EXPECT_EQ(poses[0].chrom_index, 0);
    EXPECT_EQ(poses[1].chrom_index, 1);
    EXPECT_EQ(poses[2].chrom_index, 2);
}

TEST(PoseSortTest, ReverseSorted) {
    chromosome c = make_chrom(-5.0);
    std::vector<Pose> poses;
    poses.push_back(make_pose(c, 2, 3, 1.0f));
    poses.push_back(make_pose(c, 1, 2, 1.0f));
    poses.push_back(make_pose(c, 0, 1, 1.0f));

    std::sort(poses.begin(), poses.end(), PoseClassifier());

    EXPECT_EQ(poses[0].chrom_index, 0);
    EXPECT_EQ(poses[1].chrom_index, 1);
    EXPECT_EQ(poses[2].chrom_index, 2);
}

TEST(PoseSortTest, SingleElement) {
    chromosome c = make_chrom(-5.0);
    std::vector<Pose> poses;
    poses.push_back(make_pose(c, 42, 7, 3.0f));

    std::sort(poses.begin(), poses.end(), PoseClassifier());

    ASSERT_EQ(poses.size(), 1u);
    EXPECT_EQ(poses[0].chrom_index, 42);
}

TEST(PoseSortTest, EmptyVector) {
    std::vector<Pose> poses;
    std::sort(poses.begin(), poses.end(), PoseClassifier());
    EXPECT_TRUE(poses.empty());
}

// ═══════════════════════════════════════════════════════════════════════
// Pose CF field propagation from chromosome
// ═══════════════════════════════════════════════════════════════════════

TEST(PoseConstructionTest, CFMatchesChromEvalue) {
    chromosome c = make_chrom(-42.5);
    Pose p = make_pose(c, 0, 1, 1.0f);
    EXPECT_DOUBLE_EQ(p.CF, -42.5);
}

TEST(PoseConstructionTest, IndexAndOrderPropagated) {
    chromosome c = make_chrom(0.0);
    Pose p = make_pose(c, 7, 3, 2.5f);
    EXPECT_EQ(p.chrom_index, 7);
    EXPECT_EQ(p.order, 3);
    EXPECT_FLOAT_EQ(p.reachDist, 2.5f);
}

TEST(PoseConstructionTest, ZeroEnergy) {
    chromosome c = make_chrom(0.0);
    Pose p = make_pose(c, 0, 1, 0.0f);
    EXPECT_DOUBLE_EQ(p.CF, 0.0);
}

TEST(PoseConstructionTest, NegativeEnergy) {
    chromosome c = make_chrom(-999.99);
    Pose p = make_pose(c, 0, 1, 1.0f);
    EXPECT_DOUBLE_EQ(p.CF, -999.99);
}
