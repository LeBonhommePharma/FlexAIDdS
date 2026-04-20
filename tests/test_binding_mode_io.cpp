// tests/test_binding_mode_io.cpp
// Unit tests for Pose construction, BindingMode pose management,
// and BindingPopulation I/O plumbing.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <cstring>
#include <vector>

#include "BindingMode.h"
#include "gaboom.h"

// ═══════════════════════════════════════════════════════════════════════
// Test fixture — minimal mock structures
// ═══════════════════════════════════════════════════════════════════════

class BindingModeIOTest : public ::testing::Test {
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

    static constexpr int NCHROM = 5;
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
            chroms[i].app_evalue = static_cast<double>(-10 + i * 2);
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

    Pose make_pose(int idx, double cf_override = -999.0) {
        std::vector<float> empty;
        Pose p(&chroms[idx], idx, /*order=*/0, /*dist=*/0.0f, TEST_TEMP, empty);
        if (cf_override != -999.0) p.CF = cf_override;
        return p;
    }
};

// ═══════════════════════════════════════════════════════════════════════
// Pose construction
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeIOTest, PoseStoresChromIndex) {
    Pose p = make_pose(3);
    EXPECT_EQ(p.chrom_index, 3);
}

TEST_F(BindingModeIOTest, PoseStoresChromPointer) {
    Pose p = make_pose(2);
    EXPECT_EQ(p.chrom, &chroms[2]);
}

TEST_F(BindingModeIOTest, PoseCFFromAppEvalue) {
    // Pose constructor sets CF = chrom->app_evalue
    Pose p = make_pose(0);  // chroms[0].app_evalue = -10
    EXPECT_DOUBLE_EQ(p.CF, -10.0);
}

TEST_F(BindingModeIOTest, PoseCFOverride) {
    Pose p = make_pose(1, -42.0);
    EXPECT_DOUBLE_EQ(p.CF, -42.0);
}

TEST_F(BindingModeIOTest, PoseDefaultCCBMFields) {
    Pose p = make_pose(0);
    EXPECT_EQ(p.model_index, 0);
    EXPECT_EQ(p.model_coords, nullptr);
    EXPECT_DOUBLE_EQ(p.receptor_strain, 0.0);
}

TEST_F(BindingModeIOTest, PoseTotalEnergyNoStrain) {
    Pose p = make_pose(0, -15.0);
    EXPECT_DOUBLE_EQ(p.total_energy(), -15.0);
}

TEST_F(BindingModeIOTest, PoseTotalEnergyWithStrain) {
    Pose p = make_pose(0, -15.0);
    p.receptor_strain = 3.5;
    EXPECT_DOUBLE_EQ(p.total_energy(), -11.5);
}

TEST_F(BindingModeIOTest, PoseBoltzmannWeightPositive) {
    Pose p = make_pose(0);  // CF = -10 at T=300
    EXPECT_GT(p.boltzmann_weight, 0.0);
}

TEST_F(BindingModeIOTest, PoseVPoseEmpty) {
    std::vector<float> empty;
    Pose p(&chroms[0], 0, 0, 0.0f, TEST_TEMP, empty);
    EXPECT_TRUE(p.vPose.empty());
}

TEST_F(BindingModeIOTest, PoseVPosePreserved) {
    std::vector<float> coords = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    Pose p(&chroms[0], 0, 0, 0.0f, TEST_TEMP, coords);
    ASSERT_EQ(p.vPose.size(), 6u);
    EXPECT_FLOAT_EQ(p.vPose[0], 1.0f);
    EXPECT_FLOAT_EQ(p.vPose[5], 6.0f);
}

TEST_F(BindingModeIOTest, PoseOrderAndReachDist) {
    Pose p = make_pose(0);
    p.order    = 5;
    p.reachDist = 1.23f;
    EXPECT_EQ(p.order, 5);
    EXPECT_FLOAT_EQ(p.reachDist, 1.23f);
}

// ═══════════════════════════════════════════════════════════════════════
// BindingMode pose management
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeIOTest, EmptyModeHasSizeZero) {
    BindingMode mode(pop);
    EXPECT_EQ(mode.get_BindingMode_size(), 0);
}

TEST_F(BindingModeIOTest, AddPoseIncrementsSize) {
    BindingMode mode(pop);
    Pose p = make_pose(0, -10.0);
    mode.add_Pose(p);
    EXPECT_EQ(mode.get_BindingMode_size(), 1);

    Pose p2 = make_pose(1, -12.0);
    mode.add_Pose(p2);
    EXPECT_EQ(mode.get_BindingMode_size(), 2);
}

TEST_F(BindingModeIOTest, GetPoseReturnsCorrectCF) {
    BindingMode mode(pop);
    Pose p1 = make_pose(0, -10.0);
    Pose p2 = make_pose(1, -20.0);
    mode.add_Pose(p1);
    mode.add_Pose(p2);

    EXPECT_DOUBLE_EQ(mode.get_pose(0).CF, -10.0);
    EXPECT_DOUBLE_EQ(mode.get_pose(1).CF, -20.0);
}

TEST_F(BindingModeIOTest, GetPosesReturnsAll) {
    BindingMode mode(pop);
    for (int i = 0; i < 3; ++i) {
        Pose p = make_pose(i, -10.0 - i * 5.0);
        mode.add_Pose(p);
    }
    const auto& poses = mode.get_poses();
    ASSERT_EQ(poses.size(), 3u);
    EXPECT_DOUBLE_EQ(poses[0].CF, -10.0);
    EXPECT_DOUBLE_EQ(poses[1].CF, -15.0);
    EXPECT_DOUBLE_EQ(poses[2].CF, -20.0);
}

TEST_F(BindingModeIOTest, ClearPosesResetsSize) {
    BindingMode mode(pop);
    for (int i = 0; i < 3; ++i) {
        Pose p = make_pose(i, -10.0);
        mode.add_Pose(p);
    }
    EXPECT_EQ(mode.get_BindingMode_size(), 3);

    mode.clear_Poses();
    EXPECT_EQ(mode.get_BindingMode_size(), 0);
}

TEST_F(BindingModeIOTest, ClearAndReaddWorks) {
    BindingMode mode(pop);
    Pose p1 = make_pose(0, -10.0);
    mode.add_Pose(p1);
    mode.clear_Poses();

    Pose p2 = make_pose(1, -20.0);
    mode.add_Pose(p2);
    ASSERT_EQ(mode.get_BindingMode_size(), 1);
    EXPECT_DOUBLE_EQ(mode.get_pose(0).CF, -20.0);
}

// ═══════════════════════════════════════════════════════════════════════
// BindingPopulation management
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeIOTest, EmptyPopulationHasSizeZero) {
    EXPECT_EQ(pop->get_Population_size(), 0);
}

TEST_F(BindingModeIOTest, AddBindingModeIncrementsSize) {
    BindingMode m1(pop);
    Pose p1 = make_pose(0, -10.0);
    m1.add_Pose(p1);
    BindingMode m2(pop);
    Pose p2 = make_pose(1, -12.0);
    m2.add_Pose(p2);
    pop->add_BindingMode(m1);
    EXPECT_EQ(pop->get_Population_size(), 1);
    pop->add_BindingMode(m2);
    EXPECT_EQ(pop->get_Population_size(), 2);
}

TEST_F(BindingModeIOTest, GetBindingModeByIndex) {
    // Modes are sorted by energy (Entropize) after add_BindingMode.
    // m1 has CF=-5,-6 → higher energy; m2 has CF=-20,-21,-22 → lower energy
    // After sort: index 0 = m2 (3 poses), index 1 = m1 (2 poses)
    BindingMode m1(pop);
    for (int i = 0; i < 2; ++i) {
        Pose p = make_pose(i, -5.0 - i);
        m1.add_Pose(p);
    }
    BindingMode m2(pop);
    for (int i = 0; i < 3; ++i) {
        Pose p = make_pose(i, -20.0 - i);
        m2.add_Pose(p);
    }

    pop->add_BindingMode(m1);
    pop->add_BindingMode(m2);

    // After Entropize sort: m2 (lower energy) comes first
    EXPECT_EQ(pop->get_binding_mode(0).get_BindingMode_size(), 3);
    EXPECT_EQ(pop->get_binding_mode(1).get_BindingMode_size(), 2);
}

TEST_F(BindingModeIOTest, GetBindingModesReturnsAll) {
    BindingMode m1(pop);
    BindingMode m2(pop);
    Pose p1 = make_pose(0, -10.0);
    Pose p2 = make_pose(1, -12.0);
    m1.add_Pose(p1);
    m2.add_Pose(p2);

    pop->add_BindingMode(m1);
    pop->add_BindingMode(m2);

    const auto& modes = pop->get_binding_modes();
    ASSERT_EQ(modes.size(), 2u);
}

TEST_F(BindingModeIOTest, PopulationTemperatureStored) {
    EXPECT_EQ(pop->Temperature, TEST_TEMP);
}

// ═══════════════════════════════════════════════════════════════════════
// CCBM fields on Pose
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeIOTest, PoseModelIndexMutable) {
    Pose p = make_pose(0);
    p.model_index = 3;
    EXPECT_EQ(p.model_index, 3);
}

TEST_F(BindingModeIOTest, PoseReceptorStrainMutable) {
    Pose p = make_pose(0, -12.0);
    p.receptor_strain = 5.0;
    EXPECT_DOUBLE_EQ(p.total_energy(), -7.0);
}

// ═══════════════════════════════════════════════════════════════════════
// Multiple modes with varying pose counts
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeIOTest, ManyModesInPopulation) {
    for (int m = 0; m < 10; ++m) {
        BindingMode mode(pop);
        for (int i = 0; i < m + 1; ++i) {
            Pose p = make_pose(i % NCHROM, -10.0 - m - i);
            mode.add_Pose(p);
        }
        pop->add_BindingMode(mode);
    }
    EXPECT_EQ(pop->get_Population_size(), 10);

    // Modes are sorted by energy (Entropize). The mode added last (m=9)
    // has the most negative CF values, so it has the lowest energy and
    // appears first. Verify total pose count across all modes is preserved.
    int total_poses = 0;
    for (int m = 0; m < 10; ++m) {
        total_poses += pop->get_binding_mode(m).get_BindingMode_size();
    }
    EXPECT_EQ(total_poses, 10 + 9 + 8 + 7 + 6 + 5 + 4 + 3 + 2 + 1);  // = 55
}

// ═══════════════════════════════════════════════════════════════════════
// PoseClassifier ordering
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeIOTest, PoseClassifierByOrder) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0);
    p1.order = 1;
    Pose p2 = make_pose(1);
    p2.order = 2;
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(BindingModeIOTest, PoseClassifierByReachDistOnEqualOrder) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0);
    p1.order = 1;
    p1.reachDist = 0.5f;
    Pose p2 = make_pose(1);
    p2.order = 1;
    p2.reachDist = 1.0f;
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(BindingModeIOTest, PoseClassifierByChromIndexOnEqualOrderAndDist) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0);
    p1.order = 1;
    p1.reachDist = 1.0f;
    Pose p2 = make_pose(1);
    p2.order = 1;
    p2.reachDist = 1.0f;
    EXPECT_TRUE(cmp(p1, p2));
    EXPECT_FALSE(cmp(p2, p1));
}

TEST_F(BindingModeIOTest, PoseClassifierEqualReturnsFalse) {
    PoseClassifier cmp;
    Pose p1 = make_pose(0);
    p1.order = 1;
    p1.reachDist = 1.0f;
    Pose p2 = make_pose(0);  // same chrom_index
    p2.order = 1;
    p2.reachDist = 1.0f;
    EXPECT_FALSE(cmp(p1, p2));
}
