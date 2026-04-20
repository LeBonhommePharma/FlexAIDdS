// tests/test_binding_mode_comparison.cpp
// Unit tests for BindingMode comparison logic, representative election,
// delta_G_relative_to, PoseClassifier, and population-level operations.
//
// NOTE: Pose::operator< and BindingMode::operator< are declared as `inline`
// in BindingMode.h but defined in BindingMode.cpp — they are only callable
// from within that TU.  We test the equivalent logic via PoseClassifier
// (public functor) and compute_energy() (public method).
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <cstring>
#include <vector>
#include <algorithm>

#include "BindingMode.h"
#include "gaboom.h"

// ═══════════════════════════════════════════════════════════════════════
// Test fixture — minimal mock structures
// ═══════════════════════════════════════════════════════════════════════

class ComparisonTest : public ::testing::Test {
protected:
    FA_Global*    fa  = nullptr;
    GB_Global*    gb  = nullptr;
    VC_Global*    vc  = nullptr;
    chromosome*   chroms = nullptr;
    genlim*       glim   = nullptr;
    atom*         atoms  = nullptr;
    resid*        res    = nullptr;
    gridpoint*    grid   = nullptr;
    BindingPopulation* pop = nullptr;

    static constexpr int NCHROM = 8;
    static constexpr int NGENES = 6;
    static constexpr uint TEST_TEMP = 300;

    void SetUp() override {
        fa = new FA_Global();
        std::memset(static_cast<void*>(fa), 0, sizeof(FA_Global));
        fa->temperature = TEST_TEMP;

        gb = new GB_Global();
        std::memset(gb, 0, sizeof(GB_Global));
        gb->num_genes = NGENES;

        vc = new VC_Global();
        std::memset(vc, 0, sizeof(VC_Global));

        chroms = new chromosome[NCHROM];
        for (int i = 0; i < NCHROM; ++i) {
            chroms[i].genes = new gene[NGENES];
            std::memset(chroms[i].genes, 0, sizeof(gene) * NGENES);
            chroms[i].evalue     = 0.0;
            chroms[i].app_evalue = static_cast<double>(-20 + i * 3);
            chroms[i].fitnes     = 0.0;
            chroms[i].status     = 'n';
        }
        glim  = new genlim[NGENES];
        std::memset(glim, 0, sizeof(genlim) * NGENES);
        atoms = new atom[10];
        std::memset(atoms, 0, sizeof(atom) * 10);
        res   = new resid[2];
        std::memset(res, 0, sizeof(resid) * 2);
        grid  = new gridpoint[100];
        std::memset(grid, 0, sizeof(gridpoint) * 100);

        pop = new BindingPopulation(fa, gb, vc, chroms, glim, atoms, res, grid, NCHROM);
    }

    void TearDown() override {
        delete pop;
        delete[] grid;
        delete[] res;
        delete[] atoms;
        delete[] glim;
        for (int i = 0; i < NCHROM; ++i) delete[] chroms[i].genes;
        delete[] chroms;
        delete vc;
        delete gb;
        delete fa;
    }

    Pose make_pose(int idx, double cf, int order = 0, float reachDist = 0.0f) {
        std::vector<float> empty;
        Pose p(&chroms[idx], idx, order, reachDist, TEST_TEMP, empty);
        p.CF = cf;
        return p;
    }
};

// ═══════════════════════════════════════════════════════════════════════
// Pose ordering logic — tested via PoseClassifier (public functor)
// PoseClassifier implements the same logic as Pose::operator<:
//   order → reachDist → chrom_index
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, PoseClassifierLessByOrder) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0, -10.0, /*order=*/1);
    Pose p2 = make_pose(1, -20.0, /*order=*/2);
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(ComparisonTest, PoseClassifierLessByReachDistOnEqualOrder) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0, -10.0, /*order=*/3, /*reachDist=*/1.0f);
    Pose p2 = make_pose(1, -20.0, /*order=*/3, /*reachDist=*/2.0f);
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(ComparisonTest, PoseClassifierLessByChromIndexOnEqualOrderAndDist) {
    PoseClassifier cmp;
    Pose p1 = make_pose(2, -10.0, /*order=*/3, /*reachDist=*/1.5f);
    Pose p2 = make_pose(5, -20.0, /*order=*/3, /*reachDist=*/1.5f);
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(ComparisonTest, PoseClassifierEqualReturnsFalse) {
    PoseClassifier cmp;
    Pose p1 = make_pose(3, -10.0, /*order=*/2, /*reachDist=*/1.0f);
    Pose p2 = make_pose(3, -20.0, /*order=*/2, /*reachDist=*/1.0f);
    EXPECT_FALSE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(ComparisonTest, PoseClassifierTransitive) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0, -5.0,  /*order=*/1, /*reachDist=*/0.5f);
    Pose p2 = make_pose(1, -10.0, /*order=*/2, /*reachDist=*/0.3f);
    Pose p3 = make_pose(2, -15.0, /*order=*/3, /*reachDist=*/0.1f);
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_TRUE(cmp(p2, p3));
    EXPECT_TRUE(cmp(p1, p3));
}

TEST_F(ComparisonTest, PoseClassifierSortVector) {
    std::vector<Pose> poses;
    poses.push_back(make_pose(0, -5.0,  /*order=*/3, /*reachDist=*/2.0f));
    poses.push_back(make_pose(1, -10.0, /*order=*/1, /*reachDist=*/0.5f));
    poses.push_back(make_pose(2, -15.0, /*order=*/1, /*reachDist=*/1.0f));
    poses.push_back(make_pose(3, -20.0, /*order=*/2, /*reachDist=*/0.0f));

    std::sort(poses.begin(), poses.end(), PoseClassifier());

    EXPECT_EQ(poses[0].order, 1);
    EXPECT_FLOAT_EQ(poses[0].reachDist, 0.5f);
    EXPECT_EQ(poses[1].order, 1);
    EXPECT_FLOAT_EQ(poses[1].reachDist, 1.0f);
    EXPECT_EQ(poses[2].order, 2);
    EXPECT_EQ(poses[3].order, 3);
}

TEST_F(ComparisonTest, PoseClassifierSortDescendingInput) {
    std::vector<Pose> poses;
    for (int i = 7; i >= 0; --i) {
        poses.push_back(make_pose(i, -(double)i, i, (float)(7 - i)));
    }
    std::sort(poses.begin(), poses.end(), PoseClassifier());

    for (size_t i = 1; i < poses.size(); ++i) {
        EXPECT_LE(poses[i - 1].order, poses[i].order);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// BindingMode energy comparison — via compute_energy()
// BindingMode::operator< compares compute_energy(), same logic as EnergyComparator
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, ModeEnergyLowerCFIsLess) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -20.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -5.0);
    m2.add_Pose(p2);

    EXPECT_LT(m1.compute_energy(), m2.compute_energy());
}

TEST_F(ComparisonTest, ModeEnergyEqualCFIsEqual) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -10.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -10.0);
    m2.add_Pose(p2);

    EXPECT_NEAR(m1.compute_energy(), m2.compute_energy(), 1e-10);
}

TEST_F(ComparisonTest, ModeEnergyTransitive) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -30.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -20.0);
    m2.add_Pose(p2);

    BindingMode m3(pop);
    Pose p3 = make_pose(2, -10.0);
    m3.add_Pose(p3);

    EXPECT_LT(m1.compute_energy(), m2.compute_energy());
    EXPECT_LT(m2.compute_energy(), m3.compute_energy());
    EXPECT_LT(m1.compute_energy(), m3.compute_energy());
}

// ═══════════════════════════════════════════════════════════════════════
// elect_Representative — CF mode (useOPTICSordering=false)
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, ElectRepresentativeByLowestCF) {
    BindingMode mode(pop);
    Pose p1 = make_pose(0, -5.0);
    Pose p2 = make_pose(1, -20.0);
    Pose p3 = make_pose(2, -10.0);
    mode.add_Pose(p1);
    mode.add_Pose(p2);
    mode.add_Pose(p3);

    auto rep = mode.elect_Representative(false);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_DOUBLE_EQ(rep->CF, -20.0);
    EXPECT_EQ(rep->chrom_index, 1);
}

TEST_F(ComparisonTest, ElectRepresentativeSinglePose) {
    BindingMode mode(pop);
    Pose p = make_pose(0, -15.0);
    mode.add_Pose(p);

    auto rep = mode.elect_Representative(false);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_DOUBLE_EQ(rep->CF, -15.0);
}

TEST_F(ComparisonTest, ElectRepresentativeAllEqualCF) {
    BindingMode mode(pop);
    for (int i = 0; i < 4; ++i) {
        Pose p = make_pose(i, -10.0);
        mode.add_Pose(p);
    }

    auto rep = mode.elect_Representative(false);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_EQ(rep->chrom_index, 0);
}

// ═══════════════════════════════════════════════════════════════════════
// elect_Representative — OPTICS mode (useOPTICSordering=true)
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, ElectRepresentativeByLowestReachDist) {
    BindingMode mode(pop);
    Pose p1 = make_pose(0, -5.0,  1, 3.0f);
    Pose p2 = make_pose(1, -20.0, 1, 0.5f);
    Pose p3 = make_pose(2, -10.0, 1, 1.5f);
    mode.add_Pose(p1);
    mode.add_Pose(p2);
    mode.add_Pose(p3);

    auto rep = mode.elect_Representative(true);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_FLOAT_EQ(rep->reachDist, 0.5f);
}

TEST_F(ComparisonTest, ElectRepresentativeOPTICSOrder) {
    // rep starts as first pose. Any later pose with a defined, lower reachDist wins.
    BindingMode mode(pop);
    Pose p1 = make_pose(0, -5.0,  1, 3.0f);
    Pose p2 = make_pose(1, -20.0, 1, 0.5f);
    Pose p3 = make_pose(2, -10.0, 1, 1.5f);
    mode.add_Pose(p1);
    mode.add_Pose(p2);
    mode.add_Pose(p3);

    auto rep = mode.elect_Representative(true);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_FLOAT_EQ(rep->reachDist, 0.5f);
}

TEST_F(ComparisonTest, ElectRepresentativeOPTICSAllUndefined) {
    BindingMode mode(pop);
    for (int i = 0; i < 3; ++i) {
        Pose p = make_pose(i, -10.0 - i, 1, -0.1f);
        mode.add_Pose(p);
    }

    auto rep = mode.elect_Representative(true);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_EQ(rep->chrom_index, 0);
}

// ═══════════════════════════════════════════════════════════════════════
// delta_G_relative_to — ΔG between binding modes
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, DeltaGRelativeLowerCFIsNegative) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -20.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -5.0);
    m2.add_Pose(p2);

    double dG = m1.delta_G_relative_to(m2);
    EXPECT_LT(dG, 0.0);
}

TEST_F(ComparisonTest, DeltaGRelativeHigherCFIsPositive) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -5.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -20.0);
    m2.add_Pose(p2);

    double dG = m1.delta_G_relative_to(m2);
    EXPECT_GT(dG, 0.0);
}

TEST_F(ComparisonTest, DeltaGRelativeEqualIsZero) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -10.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -10.0);
    m2.add_Pose(p2);

    double dG = m1.delta_G_relative_to(m2);
    EXPECT_NEAR(dG, 0.0, 1e-10);
}

TEST_F(ComparisonTest, DeltaGRelativeIsAntisymmetric) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -15.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -8.0);
    m2.add_Pose(p2);

    double dG12 = m1.delta_G_relative_to(m2);
    double dG21 = m2.delta_G_relative_to(m1);
    EXPECT_NEAR(dG12, -dG21, 1e-10);
}

// ═══════════════════════════════════════════════════════════════════════
// EnergyComparator — tested indirectly via compute_energy() comparison
// The EnergyComparator struct inside BindingPopulation uses compute_energy()
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, EnergyComparatorMatchesComputeEnergy) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -30.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -10.0);
    m2.add_Pose(p2);

    BindingMode m3(pop);
    Pose p3 = make_pose(2, -20.0);
    m3.add_Pose(p3);

    // Verify compute_energy ordering matches expected CF ordering
    EXPECT_LT(m1.compute_energy(), m3.compute_energy());
    EXPECT_LT(m3.compute_energy(), m2.compute_energy());
    EXPECT_LT(m1.compute_energy(), m2.compute_energy());
}

// ═══════════════════════════════════════════════════════════════════════
// Population-level operations with multiple modes
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, PopulationDeltaGMatrix) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -20.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -10.0);
    m2.add_Pose(p2);

    BindingMode m3(pop);
    Pose p3 = make_pose(2, -5.0);
    m3.add_Pose(p3);

    pop->add_BindingMode(m1);
    pop->add_BindingMode(m2);
    pop->add_BindingMode(m3);

    auto dG_matrix = pop->get_deltaG_matrix();
    ASSERT_EQ(dG_matrix.size(), 3u);

    // Diagonal should be ~0
    EXPECT_NEAR(dG_matrix[0][0], 0.0, 1e-10);
    EXPECT_NEAR(dG_matrix[1][1], 0.0, 1e-10);
    EXPECT_NEAR(dG_matrix[2][2], 0.0, 1e-10);

    // Antisymmetric off-diagonal
    EXPECT_NEAR(dG_matrix[0][1], -dG_matrix[1][0], 1e-10);
    EXPECT_LT(dG_matrix[0][1], 0.0);  // mode0 (-20) - mode1 (-10) = -10
}

TEST_F(ComparisonTest, PopulationShannonEntropyPositive) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -10.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -10.01);
    m2.add_Pose(p2);

    pop->add_BindingMode(m1);
    pop->add_BindingMode(m2);

    double S = pop->get_shannon_entropy();
    EXPECT_GT(S, 0.0);
}

TEST_F(ComparisonTest, PopulationComputeDeltaG) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -25.0);
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -15.0);
    m2.add_Pose(p2);

    pop->add_BindingMode(m1);
    pop->add_BindingMode(m2);

    double dG = pop->compute_delta_G(
        pop->get_binding_mode(0),
        pop->get_binding_mode(1)
    );
    EXPECT_LT(dG, 0.0);
}

// ═══════════════════════════════════════════════════════════════════════
// CCBM: total_energy affects comparisons
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, CCBMTotalEnergyWithStrain) {
    Pose p1 = make_pose(0, -15.0);
    p1.receptor_strain = 3.0;
    EXPECT_DOUBLE_EQ(p1.total_energy(), -12.0);

    Pose p2 = make_pose(1, -14.0);
    p2.receptor_strain = 0.0;
    EXPECT_DOUBLE_EQ(p2.total_energy(), -14.0);

    EXPECT_LT(p2.total_energy(), p1.total_energy());
}

TEST_F(ComparisonTest, CCBMStrainChangesModeOrdering) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -20.0);
    p1.receptor_strain = 15.0;  // total = -5.0
    m1.add_Pose(p1);

    BindingMode m2(pop);
    Pose p2 = make_pose(1, -8.0);
    p2.receptor_strain = 0.0;   // total = -8.0
    m2.add_Pose(p2);

    // m2 has lower compute_energy() than m1 despite higher raw CF
    EXPECT_LT(m2.compute_energy(), m1.compute_energy());
}

// ═══════════════════════════════════════════════════════════════════════
// Edge cases
// ═══════════════════════════════════════════════════════════════════════

TEST_F(ComparisonTest, SingleModePopulationDeltaGMatrix) {
    BindingMode m(pop);
    Pose p = make_pose(0, -10.0);
    m.add_Pose(p);
    pop->add_BindingMode(m);

    auto dG_matrix = pop->get_deltaG_matrix();
    ASSERT_EQ(dG_matrix.size(), 1u);
    ASSERT_EQ(dG_matrix[0].size(), 1u);
    EXPECT_NEAR(dG_matrix[0][0], 0.0, 1e-10);
}

TEST_F(ComparisonTest, EmptyModeElectRepresentative) {
    BindingMode mode(pop);
    auto rep = mode.elect_Representative(false);
    EXPECT_EQ(rep, mode.get_poses().end());
}

TEST_F(ComparisonTest, ModeWithManyPosesElectsLowestCF) {
    BindingMode mode(pop);
    double lowest_cf = 0.0;
    int lowest_idx = -1;
    for (int i = 0; i < NCHROM; ++i) {
        double cf = -5.0 * (i + 1);
        if (i == 0 || cf < lowest_cf) {
            lowest_cf = cf;
            lowest_idx = i;
        }
        Pose p = make_pose(i, cf, i, (float)(NCHROM - i));
        mode.add_Pose(p);
    }

    auto rep = mode.elect_Representative(false);
    ASSERT_NE(rep, mode.get_poses().end());
    EXPECT_DOUBLE_EQ(rep->CF, lowest_cf);
    EXPECT_EQ(rep->chrom_index, lowest_idx);
}
