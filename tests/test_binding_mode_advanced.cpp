// tests/test_binding_mode_advanced.cpp
// Tests for BindingMode advanced methods: delta_G_relative_to, elect_Representative,
// free_energy_profile, compute_vibrational_correction
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/BindingMode.h"
#include "../LIB/statmech.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════
// TEST FIXTURE
// ═══════════════════════════════════════════════════════════════════════

class BindingModeAdvancedTest : public ::testing::Test {
protected:
    FA_Global* fa;
    GB_Global* gb;
    VC_Global* vc;
    chromosome* chroms;
    genlim* gene_lim;
    atom* atoms;
    resid* residue;
    gridpoint* cleftgrid;
    BindingPopulation* population;

    const double TEMP = 300.0;
    const double EPS = 1e-6;

    void SetUp() override {
        fa = new FA_Global();
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
        std::memset(fa, 0, sizeof(FA_Global));
#pragma clang diagnostic pop
        fa->temperature = static_cast<uint>(TEMP);
        // Advanced thermo/profile tests use physical per-mode F.
        fa->force_cf_rank_emission = true;

        gb = new GB_Global();
        std::memset(gb, 0, sizeof(GB_Global));
        gb->num_genes = 6;

        vc = new VC_Global();
        std::memset(vc, 0, sizeof(VC_Global));

        chroms = new chromosome[10];
        for (int i = 0; i < 10; ++i) {
            chroms[i].genes = new gene[gb->num_genes];
            std::memset(chroms[i].genes, 0, sizeof(gene) * gb->num_genes);
            chroms[i].evalue = 0.0;
            chroms[i].app_evalue = 0.0;
            chroms[i].fitnes = 0.0;
            chroms[i].status = 'n';
        }
        gene_lim = new genlim[gb->num_genes];
        atoms = new atom[5];
        std::memset(atoms, 0, sizeof(atom) * 5);
        residue = new resid[1];
        std::memset(residue, 0, sizeof(resid) * 1);
        cleftgrid = new gridpoint[100];
        std::memset(cleftgrid, 0, sizeof(gridpoint) * 100);

        population = new BindingPopulation(
            fa, gb, vc, chroms, gene_lim, atoms, residue, cleftgrid, 10
        );
    }

    void TearDown() override {
        delete population;
        delete[] cleftgrid;
        delete[] residue;
        delete[] atoms;
        delete[] gene_lim;
        for (int i = 0; i < 10; ++i) delete[] chroms[i].genes;
        delete[] chroms;
        delete vc;
        delete gb;
        delete fa;
    }

    // Helper: create pose and add it to a BindingMode (avoids rvalue-to-lvalue-ref issue)
    void add_pose(BindingMode& mode, double cf, int idx) {
        std::vector<float> v;
        Pose p(&chroms[idx], idx, 0, 0.0f, static_cast<uint>(TEMP), v);
        p.CF = cf;
        mode.add_Pose(p);
    }
};

// ═══════════════════════════════════════════════════════════════════════
// delta_G_relative_to Tests
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, DeltaGRelativeToSelfIsZero) {
    BindingMode mode(population);
    for (int i = 0; i < 3; ++i)
        add_pose(mode, -10.0 - i, i);

    double dg = mode.delta_G_relative_to(mode);
    EXPECT_NEAR(dg, 0.0, EPS);
}

TEST_F(BindingModeAdvancedTest, DeltaGRelativeToIsAntisymmetric) {
    BindingMode mode_a(population);
    BindingMode mode_b(population);

    for (int i = 0; i < 3; ++i) {
        add_pose(mode_a, -15.0 - i, i);
        add_pose(mode_b, -10.0 - i, i);
    }

    double dg_ab = mode_a.delta_G_relative_to(mode_b);
    double dg_ba = mode_b.delta_G_relative_to(mode_a);

    EXPECT_NEAR(dg_ab, -dg_ba, EPS)
        << "delta_G should be antisymmetric";
}

TEST_F(BindingModeAdvancedTest, DeltaGMatchesFreeEnergyDifference) {
    BindingMode mode_a(population);
    BindingMode mode_b(population);

    for (int i = 0; i < 4; ++i) {
        add_pose(mode_a, -20.0 + i, i);
        add_pose(mode_b, -12.0 + i, i);
    }

    double fa_val = mode_a.get_free_energy();
    double fb_val = mode_b.get_free_energy();
    double dg = mode_a.delta_G_relative_to(mode_b);

    // ΔG(A→B) = F_A - F_B
    EXPECT_NEAR(dg, fa_val - fb_val, EPS);
}

TEST_F(BindingModeAdvancedTest, DeltaGLowerEnergyIsNegative) {
    BindingMode lower(population);
    BindingMode higher(population);

    for (int i = 0; i < 3; ++i) {
        add_pose(lower, -20.0 - i, i);
        add_pose(higher, -5.0 - i, i);
    }

    // Lower energy mode relative to higher should be negative
    double dg = lower.delta_G_relative_to(higher);
    EXPECT_LT(dg, 0.0);
}

// ═══════════════════════════════════════════════════════════════════════
// elect_Representative Tests
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, ElectRepresentativeFindsLowestCF) {
    BindingMode mode(population);

    // Add poses with different CFs; -25 is the best (lowest)
    add_pose(mode, -10.0, 0);
    add_pose(mode, -25.0, 1);
    add_pose(mode, -15.0, 2);

    auto rep = mode.elect_Representative(false);  // not OPTICS ordering
    EXPECT_NEAR(rep->CF, -25.0, EPS)
        << "Representative should be the pose with lowest CF";
}

TEST_F(BindingModeAdvancedTest, ElectRepresentativeSinglePose) {
    BindingMode mode(population);
    add_pose(mode, -8.0, 0);

    auto rep = mode.elect_Representative(false);
    EXPECT_NEAR(rep->CF, -8.0, EPS);
}

// ═══════════════════════════════════════════════════════════════════════
// free_energy_profile Tests (WHAM-based)
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, FreeEnergyProfileReturnsCorrectBinCount) {
    BindingMode mode(population);
    for (int i = 0; i < 5; ++i)
        add_pose(mode, -10.0 - i, i);

    // Provide 1D coordinates matching number of poses
    std::vector<double> coords = {0.0, 1.0, 2.0, 3.0, 4.0};
    int nbins = 10;

    auto profile = mode.free_energy_profile(coords, nbins);

    // Should return bins (might be less than nbins if some are empty)
    EXPECT_LE(profile.size(), static_cast<size_t>(nbins));
}

TEST_F(BindingModeAdvancedTest, FreeEnergyProfileBinsHaveFiniteValues) {
    BindingMode mode(population);
    for (int i = 0; i < 8; ++i)
        add_pose(mode, -10.0 - i * 0.5, i);

    std::vector<double> coords = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5};
    auto profile = mode.free_energy_profile(coords, 5);

    for (const auto& bin : profile) {
        EXPECT_TRUE(std::isfinite(bin.coord_center))
            << "Bin coordinate should be finite";
        EXPECT_TRUE(std::isfinite(bin.free_energy))
            << "Bin free energy should be finite";
    }
}

// ═══════════════════════════════════════════════════════════════════════
// get_heat_capacity Tests
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, HeatCapacityNonNegative) {
    BindingMode mode(population);
    for (int i = 0; i < 5; ++i)
        add_pose(mode, -10.0 - i * 2.0, i);

    double cv = mode.get_heat_capacity();
    EXPECT_GE(cv, 0.0) << "Heat capacity must be non-negative";
}

TEST_F(BindingModeAdvancedTest, SinglePoseZeroHeatCapacity) {
    BindingMode mode(population);
    add_pose(mode, -12.0, 0);

    double cv = mode.get_heat_capacity();
    EXPECT_NEAR(cv, 0.0, EPS) << "Single pose has zero energy variance -> Cv = 0";
}

// ═══════════════════════════════════════════════════════════════════════
// BindingPopulation::compute_delta_G Tests
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, PopulationDeltaGConsistentWithModes) {
    BindingMode mode_a(population);
    BindingMode mode_b(population);

    for (int i = 0; i < 3; ++i) {
        add_pose(mode_a, -18.0 - i, i);
        add_pose(mode_b, -8.0 - i, i);
    }

    double pop_dg = population->compute_delta_G(mode_a, mode_b);
    double mode_dg = mode_a.delta_G_relative_to(mode_b);

    // Both compute F_A - F_B, so should be equal
    EXPECT_NEAR(pop_dg, mode_dg, EPS);
}

// ═══════════════════════════════════════════════════════════════════════
// Global ensemble Tests
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, GlobalEnsembleAggregatesModes) {
    BindingMode mode_a(population);
    BindingMode mode_b(population);

    for (int i = 0; i < 4; ++i) {
        add_pose(mode_a, -15.0 - i, i);
        add_pose(mode_b, -10.0 - i, i + 4);
    }

    population->add_BindingMode(mode_a);
    population->add_BindingMode(mode_b);

    // Population should contain both binding modes
    EXPECT_EQ(population->get_Population_size(), 2);

    // Global ensemble should have aggregated poses
    auto global = population->get_global_ensemble();
    EXPECT_GT(global.size(), 0u);
}

// ═══════════════════════════════════════════════════════════════════════
// FO .mcf sidecar (member CFs for DatasetRunner Shannon G̃)
// ═══════════════════════════════════════════════════════════════════════

TEST_F(BindingModeAdvancedTest, WriteMcfSidecar_HeadFirstFiniteOnly) {
    // Matches cluster.cpp format: one app_evalue per line; head first.
    // Non-finite values must be skipped so DatasetRunner parse_pose stays clean.
    char tmpl[] = "/tmp/flexaidds_fo_mcf_XXXXXX";
    int fd = mkstemp(tmpl);
    ASSERT_GE(fd, 0);
    close(fd);
    std::string pdb_path = std::string(tmpl) + ".pdb";
    // mkstemp created empty file at tmpl; rename conceptually to .pdb basename
    std::remove(tmpl);
    {
        FILE* pf = std::fopen(pdb_path.c_str(), "w");
        ASSERT_NE(pf, nullptr);
        std::fputs("REMARK dummy FO head\n", pf);
        std::fclose(pf);
    }

    std::vector<double> vals = {
        -12.345678,  // head
        -10.0,
        std::numeric_limits<double>::quiet_NaN(),  // skip
        -8.5,
        std::numeric_limits<double>::infinity(),   // skip
        -7.0
    };
    ASSERT_TRUE(write_mcf_sidecar(pdb_path.c_str(), vals));

    std::string mcf_path = pdb_path.substr(0, pdb_path.size() - 4) + ".mcf";
    FILE* mf = std::fopen(mcf_path.c_str(), "r");
    ASSERT_NE(mf, nullptr);
    std::vector<double> read_vals;
    char line[128];
    while (std::fgets(line, sizeof(line), mf)) {
        read_vals.push_back(std::strtod(line, nullptr));
    }
    std::fclose(mf);

    ASSERT_EQ(read_vals.size(), 4u);
    EXPECT_NEAR(read_vals[0], -12.345678, 1e-6);  // head first, %.6f
    EXPECT_NEAR(read_vals[1], -10.0, 1e-6);
    EXPECT_NEAR(read_vals[2], -8.5, 1e-6);
    EXPECT_NEAR(read_vals[3], -7.0, 1e-6);

    // First line formatting: six decimal places like cluster.cpp
    mf = std::fopen(mcf_path.c_str(), "r");
    ASSERT_NE(mf, nullptr);
    ASSERT_NE(std::fgets(line, sizeof(line), mf), nullptr);
    std::fclose(mf);
    EXPECT_STREQ(line, "-12.345678\n");

    std::remove(pdb_path.c_str());
    std::remove(mcf_path.c_str());
}

// ═══════════════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════════════

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
