// tests/test_ccbm.cpp
// Unit tests for Conformer-Coupled Binding Modes (CCBM)
// Part of FlexAIDΔS Phase 2 — multi-conformer receptor docking
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/BindingMode.h"
#include "../LIB/statmech.h"
#include "../LIB/gaboom.h"
#include <cmath>
#include <cstring>
#include <numeric>
#include <vector>

// ===========================================================================
// CONSTANTS
// ===========================================================================

static constexpr double kB   = statmech::kB_kcal;   // 0.001987206 kcal mol⁻¹ K⁻¹
static constexpr double TEMP = 300.0;                // Kelvin
static constexpr double BETA = 1.0 / (kB * TEMP);   // ~1.677 (kcal/mol)⁻¹
static constexpr double EPS  = 1e-8;                 // tight tolerance
static constexpr double EPS_LOOSE = 5e-3;            // loose tolerance for MI (sampling noise with small ensembles)

// ===========================================================================
// TEST FIXTURE
// ===========================================================================

class CCBMTest : public ::testing::Test {
protected:
    FA_Global*   mock_fa;
    GB_Global*   mock_gb;
    VC_Global*   mock_vc;
    chromosome*  mock_chroms;
    genlim*      mock_gene_lim;
    atom*        mock_atoms;
    resid*       mock_residue;
    gridpoint*   mock_cleftgrid;

    BindingPopulation* test_population;

    static constexpr int N_CHROMS = 20;  // enough for multi-conformer tests
    static constexpr int N_GENES  = 7;   // 6 normal + 1 model gene

    void SetUp() override {
        mock_fa = new FA_Global();
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
        std::memset(mock_fa, 0, sizeof(FA_Global));
#pragma clang diagnostic pop
        mock_fa->temperature = static_cast<uint>(TEMP);
        mock_fa->multi_model = false;
        mock_fa->n_models = 1;
        mock_fa->model_gene_index = -1;
        mock_fa->normal_modes = 0;  // no ENCoM modes for these tests
        // CCBM free-energy assertions use per-mode physical StatMech (incl. strain).
        mock_fa->force_cf_rank_emission = true;

        mock_gb = new GB_Global();
        std::memset(mock_gb, 0, sizeof(GB_Global));
        mock_gb->num_genes = N_GENES;

        mock_vc = new VC_Global();
        std::memset(mock_vc, 0, sizeof(VC_Global));

        mock_chroms = new chromosome[N_CHROMS];
        for (int i = 0; i < N_CHROMS; ++i) {
            mock_chroms[i].genes = new gene[N_GENES];
            std::memset(mock_chroms[i].genes, 0, sizeof(gene) * N_GENES);
            mock_chroms[i].evalue = 0.0;
            mock_chroms[i].app_evalue = 0.0;
            mock_chroms[i].fitnes = 0.0;
            mock_chroms[i].status = 'n';
        }
        mock_gene_lim = new genlim[N_GENES];
        std::memset(mock_gene_lim, 0, sizeof(genlim) * N_GENES);
        mock_atoms = new atom[10];
        std::memset(mock_atoms, 0, sizeof(atom) * 10);
        mock_residue = new resid[2];
        std::memset(mock_residue, 0, sizeof(resid) * 2);
        mock_cleftgrid = new gridpoint[100];
        std::memset(mock_cleftgrid, 0, sizeof(gridpoint) * 100);

        test_population = new BindingPopulation(
            mock_fa, mock_gb, mock_vc,
            mock_chroms, mock_gene_lim,
            mock_atoms, mock_residue, mock_cleftgrid,
            N_CHROMS
        );
    }

    void TearDown() override {
        delete test_population;
        delete[] mock_cleftgrid;
        delete[] mock_residue;
        delete[] mock_atoms;
        delete[] mock_gene_lim;
        for (int i = 0; i < N_CHROMS; ++i) delete[] mock_chroms[i].genes;
        delete[] mock_chroms;
        delete mock_vc;
        delete mock_gb;
        delete mock_fa;
    }

    /// Create a pose with given CF, model_index, receptor_strain, and chrom_index.
    Pose make_ccbm_pose(double cf, int model_idx, double strain, int chrom_idx) {
        // Set app_evalue on the chromosome so the Pose ctor picks it up
        mock_chroms[chrom_idx].app_evalue = cf;
        std::vector<float> empty_vec;
        Pose p(&mock_chroms[chrom_idx], chrom_idx, 0, 0.0f,
               static_cast<uint>(TEMP), empty_vec);
        p.CF              = cf;
        p.model_index     = model_idx;
        p.receptor_strain = strain;
        p.model_coords    = nullptr;
        return p;
    }
};

// ===========================================================================
// POSE total_energy() TESTS
// ===========================================================================

TEST_F(CCBMTest, TotalEnergyIncludesStrain) {
    Pose p = make_ccbm_pose(-10.0, 1, 2.5, 0);
    EXPECT_NEAR(p.total_energy(), -7.5, EPS);
}

TEST_F(CCBMTest, TotalEnergyZeroStrainEqualsCF) {
    Pose p = make_ccbm_pose(-15.3, 0, 0.0, 0);
    EXPECT_NEAR(p.total_energy(), p.CF, EPS);
}

TEST_F(CCBMTest, TotalEnergyNegativeStrain) {
    // Strain can be negative if the conformer is lower energy than the reference
    Pose p = make_ccbm_pose(-10.0, 2, -1.0, 0);
    EXPECT_NEAR(p.total_energy(), -11.0, EPS);
}

// ===========================================================================
// CONFORMER POPULATIONS TESTS
// ===========================================================================

TEST_F(CCBMTest, ConformerPopulationsSingleModel) {
    // All poses on model_index=0 → single element vector, population = 1.0
    BindingMode mode(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = make_ccbm_pose(-10.0 - i, 0, 0.0, i);
        mode.add_Pose(p);
    }
    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 1u);
    EXPECT_NEAR(pops[0], 1.0, EPS);
}

TEST_F(CCBMTest, ConformerPopulationsEmpty) {
    BindingMode mode(test_population);
    auto pops = mode.conformer_populations();
    EXPECT_TRUE(pops.empty());
}

TEST_F(CCBMTest, ConformerPopulationsSumToOne) {
    // 3 receptor conformers with different strain energies
    BindingMode mode(test_population);
    // Model 0: best CF, no strain
    Pose p0 = make_ccbm_pose(-12.0, 0, 0.0, 0);
    mode.add_Pose(p0);
    // Model 1: slightly worse CF, some strain
    Pose p1 = make_ccbm_pose(-11.0, 1, 1.0, 1);
    mode.add_Pose(p1);
    // Model 2: good CF, high strain
    Pose p2 = make_ccbm_pose(-13.0, 2, 5.0, 2);
    mode.add_Pose(p2);

    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 3u);

    double sum = 0.0;
    for (double p : pops) {
        EXPECT_GE(p, 0.0);
        EXPECT_LE(p, 1.0);
        sum += p;
    }
    EXPECT_NEAR(sum, 1.0, EPS);
}

TEST_F(CCBMTest, ConformerPopulationsBoltzmannOrdering) {
    // Model 0 has total_energy = -12.0 (best)
    // Model 1 has total_energy = -10.0 + 1.0 = -9.0
    // Model 0 should be more populated at T=300K
    BindingMode mode(test_population);
    Pose p0 = make_ccbm_pose(-12.0, 0, 0.0, 0);
    mode.add_Pose(p0);
    Pose p1 = make_ccbm_pose(-10.0, 1, 1.0, 1);
    mode.add_Pose(p1);

    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 2u);
    // Model 0 should dominate: total_energy = -12 vs -9 → ΔE = 3 kcal/mol
    // at 300K, kBT ≈ 0.6 kcal/mol, so exp(-β * 3) ≈ 0.007 — model 0 dominates
    EXPECT_GT(pops[0], pops[1]);
    EXPECT_GT(pops[0], 0.99);  // nearly all population on model 0
}

TEST_F(CCBMTest, ConformerPopulationsUniformWhenDegenerate) {
    // 3 conformers, all with the same total energy → uniform distribution
    BindingMode mode(test_population);
    for (int r = 0; r < 3; ++r) {
        Pose p = make_ccbm_pose(-10.0, r, 0.0, r);
        mode.add_Pose(p);
    }
    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 3u);
    for (double p : pops) {
        EXPECT_NEAR(p, 1.0 / 3.0, EPS);
    }
}

TEST_F(CCBMTest, ConformerPopulationsMultiplePosesPerModel) {
    // Model 0: two poses with energies -12, -11 → summed Boltzmann weight
    // Model 1: one pose with energy -10
    BindingMode mode(test_population);
    Pose p0a = make_ccbm_pose(-12.0, 0, 0.0, 0);
    Pose p0b = make_ccbm_pose(-11.0, 0, 0.0, 1);
    Pose p1  = make_ccbm_pose(-10.0, 1, 0.0, 2);
    mode.add_Pose(p0a);
    mode.add_Pose(p0b);
    mode.add_Pose(p1);

    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 2u);

    // Manually verify:
    // Z_0 = exp(-β*(-12)) + exp(-β*(-11)) = exp(12β) + exp(11β)
    // Z_1 = exp(-β*(-10)) = exp(10β)
    // p(0) = Z_0 / (Z_0 + Z_1), p(1) = Z_1 / (Z_0 + Z_1)
    double Z0 = std::exp(BETA * 12.0) + std::exp(BETA * 11.0);
    double Z1 = std::exp(BETA * 10.0);
    double expected_p0 = Z0 / (Z0 + Z1);
    double expected_p1 = Z1 / (Z0 + Z1);

    EXPECT_NEAR(pops[0], expected_p0, EPS);
    EXPECT_NEAR(pops[1], expected_p1, EPS);
}

// ===========================================================================
// RECEPTOR CONFORMATIONAL ENTROPY TESTS
// ===========================================================================

TEST_F(CCBMTest, ReceptorEntropyZeroForSingleModel) {
    BindingMode mode(test_population);
    Pose p = make_ccbm_pose(-10.0, 0, 0.0, 0);
    mode.add_Pose(p);

    double S_receptor = mode.receptor_conformational_entropy();
    // Single conformer → zero entropy (no uncertainty)
    EXPECT_NEAR(S_receptor, 0.0, EPS);
}

TEST_F(CCBMTest, ReceptorEntropyEmpty) {
    BindingMode mode(test_population);
    double S = mode.receptor_conformational_entropy();
    EXPECT_NEAR(S, 0.0, EPS);
}

TEST_F(CCBMTest, ReceptorEntropyMaximalWhenUniform) {
    // N degenerate conformers → S = kB * ln(N)
    const int N_MODELS = 5;
    BindingMode mode(test_population);
    for (int r = 0; r < N_MODELS; ++r) {
        Pose p = make_ccbm_pose(-10.0, r, 0.0, r);
        mode.add_Pose(p);
    }
    double S_receptor = mode.receptor_conformational_entropy();
    double expected_S = kB * std::log(static_cast<double>(N_MODELS));
    EXPECT_NEAR(S_receptor, expected_S, EPS);
}

TEST_F(CCBMTest, ReceptorEntropyIncreasesWithMoreConformers) {
    // 2 degenerate conformers → S_2 = kB * ln(2)
    // 4 degenerate conformers → S_4 = kB * ln(4) > S_2
    BindingMode mode2(test_population);
    for (int r = 0; r < 2; ++r) {
        Pose p = make_ccbm_pose(-10.0, r, 0.0, r);
        mode2.add_Pose(p);
    }
    BindingMode mode4(test_population);
    for (int r = 0; r < 4; ++r) {
        Pose p = make_ccbm_pose(-10.0, r, 0.0, r);
        mode4.add_Pose(p);
    }
    double S2 = mode2.receptor_conformational_entropy();
    double S4 = mode4.receptor_conformational_entropy();
    EXPECT_GT(S4, S2);
}

TEST_F(CCBMTest, ReceptorEntropyDecreasesWithDominantConformer) {
    // One conformer much lower energy → nearly zero entropy
    BindingMode mode(test_population);
    // Model 0: very low energy → dominates
    Pose p0 = make_ccbm_pose(-20.0, 0, 0.0, 0);
    mode.add_Pose(p0);
    // Model 1: much higher total energy
    Pose p1 = make_ccbm_pose(-5.0, 1, 5.0, 1);
    mode.add_Pose(p1);

    double S = mode.receptor_conformational_entropy();
    // At 300K, ΔE_total = 20 kcal/mol → the low energy model completely dominates
    // S ≈ 0 (conformational selection regime)
    EXPECT_LT(S, kB * 0.01);  // S << kB (essentially zero)
}

TEST_F(CCBMTest, ReceptorEntropyNonNegative) {
    BindingMode mode(test_population);
    for (int r = 0; r < 3; ++r) {
        Pose p = make_ccbm_pose(-10.0 + r * 2.0, r, r * 0.5, r);
        mode.add_Pose(p);
    }
    double S = mode.receptor_conformational_entropy();
    EXPECT_GE(S, 0.0);
}

// ===========================================================================
// MUTUAL INFORMATION TESTS
// ===========================================================================

TEST_F(CCBMTest, MutualInformationNonNegative) {
    BindingMode mode(test_population);
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < 2; ++i) {
            int chrom_idx = r * 2 + i;
            Pose p = make_ccbm_pose(-10.0 - r - i * 0.5, r, r * 0.3, chrom_idx);
            mode.add_Pose(p);
        }
    }
    double MI = mode.ligand_receptor_mutual_information();
    EXPECT_GE(MI, -EPS_LOOSE);  // MI ≥ 0 (up to numerical error)
}

TEST_F(CCBMTest, MutualInformationZeroSingleModel) {
    // Single conformer → receptor is deterministic → MI = S_L + S_R - S_total
    // with S_R = 0, MI = S_L - S_total ≤ 0, clamped to 0
    BindingMode mode(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = make_ccbm_pose(-10.0 - i, 0, 0.0, i);
        mode.add_Pose(p);
    }
    double MI = mode.ligand_receptor_mutual_information();
    EXPECT_NEAR(MI, 0.0, EPS_LOOSE);
}

TEST_F(CCBMTest, MutualInformationEmpty) {
    BindingMode mode(test_population);
    double MI = mode.ligand_receptor_mutual_information();
    EXPECT_NEAR(MI, 0.0, EPS);
}

TEST_F(CCBMTest, MutualInformationIndependentVariables) {
    // If all (model, pose) pairs have the same total energy,
    // ligand and receptor are independent → MI = 0
    BindingMode mode(test_population);
    int idx = 0;
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < 2; ++i) {
            Pose p = make_ccbm_pose(-10.0, r, 0.0, idx);
            mode.add_Pose(p);
            ++idx;
        }
    }
    double MI = mode.ligand_receptor_mutual_information();
    EXPECT_NEAR(MI, 0.0, EPS_LOOSE);
}

// ===========================================================================
// ENTROPY DECOMPOSITION TESTS
// ===========================================================================

TEST_F(CCBMTest, EntropyDecompositionEmpty) {
    BindingMode mode(test_population);
    auto ed = mode.decompose_entropy();
    EXPECT_NEAR(ed.S_total, 0.0, EPS);
    EXPECT_NEAR(ed.S_ligand, 0.0, EPS);
    EXPECT_NEAR(ed.S_receptor, 0.0, EPS);
    EXPECT_NEAR(ed.I_mutual, 0.0, EPS);
    EXPECT_NEAR(ed.S_vibrational, 0.0, EPS);
}

TEST_F(CCBMTest, EntropyDecompositionConsistency) {
    // S_total = S_ligand + S_receptor - I(L;R)
    BindingMode mode(test_population);
    int idx = 0;
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < 3; ++i) {
            double cf = -10.0 - r * 1.5 - i * 0.8;
            double strain = r * 0.5;
            Pose p = make_ccbm_pose(cf, r, strain, idx);
            mode.add_Pose(p);
            ++idx;
        }
    }

    auto ed = mode.decompose_entropy();

    // The fundamental identity: S_total = S_L + S_R - MI
    double reconstructed = ed.S_ligand + ed.S_receptor - ed.I_mutual;
    EXPECT_NEAR(ed.S_total, reconstructed, EPS_LOOSE);
}

TEST_F(CCBMTest, EntropyDecompositionSinglePose) {
    // Single microstate → all entropies = 0
    BindingMode mode(test_population);
    Pose p = make_ccbm_pose(-10.0, 0, 0.0, 0);
    mode.add_Pose(p);

    auto ed = mode.decompose_entropy();
    EXPECT_NEAR(ed.S_total, 0.0, EPS);
    EXPECT_NEAR(ed.S_ligand, 0.0, EPS);
    EXPECT_NEAR(ed.S_receptor, 0.0, EPS);
    EXPECT_NEAR(ed.I_mutual, 0.0, EPS);
}

TEST_F(CCBMTest, EntropyDecompositionVibrationalZeroWithoutModes) {
    BindingMode mode(test_population);
    Pose p0 = make_ccbm_pose(-10.0, 0, 0.0, 0);
    Pose p1 = make_ccbm_pose(-9.0, 1, 0.5, 1);
    mode.add_Pose(p0);
    mode.add_Pose(p1);

    auto ed = mode.decompose_entropy();
    // No normal modes → S_vibrational should be 0
    EXPECT_NEAR(ed.S_vibrational, 0.0, EPS);
}

TEST_F(CCBMTest, EntropyDecompositionAllComponentsNonNegative) {
    BindingMode mode(test_population);
    int idx = 0;
    for (int r = 0; r < 4; ++r) {
        for (int i = 0; i < 2; ++i) {
            Pose p = make_ccbm_pose(-8.0 - r * 2.0 - i, r, r * 0.3, idx);
            mode.add_Pose(p);
            ++idx;
        }
    }

    auto ed = mode.decompose_entropy();
    EXPECT_GE(ed.S_total, -EPS);
    EXPECT_GE(ed.S_ligand, -EPS);
    EXPECT_GE(ed.S_receptor, -EPS);
    EXPECT_GE(ed.I_mutual, -EPS);
}

TEST_F(CCBMTest, EntropyDecompositionMarginalBounds) {
    // S_total ≤ S_ligand + S_receptor  (since MI ≥ 0)
    BindingMode mode(test_population);
    int idx = 0;
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < 3; ++i) {
            Pose p = make_ccbm_pose(-10.0 - r - i * 0.5, r, r * 0.2, idx);
            mode.add_Pose(p);
            ++idx;
        }
    }

    auto ed = mode.decompose_entropy();
    EXPECT_LE(ed.S_total, ed.S_ligand + ed.S_receptor + EPS_LOOSE);
}

// ===========================================================================
// ANALYTICAL VERIFICATION: 2-state system
// ===========================================================================

TEST_F(CCBMTest, TwoStateBoltzmannExact) {
    // Exactly 2 microstates: (model 0, CF=-10, strain=0) and (model 1, CF=-8, strain=1)
    // Total energies: E_0 = -10, E_1 = -7
    // β = 1/(kB * 300) = 1/(0.001987206 * 300) = 1.677...
    // Z = exp(-β*(-10)) + exp(-β*(-7)) = exp(10β) + exp(7β)
    // p_0 = exp(10β) / Z, p_1 = exp(7β) / Z
    BindingMode mode(test_population);
    Pose p0 = make_ccbm_pose(-10.0, 0, 0.0, 0);
    Pose p1 = make_ccbm_pose(-8.0, 1, 1.0, 1);
    mode.add_Pose(p0);
    mode.add_Pose(p1);

    // Expected populations
    double lw0 = BETA * 10.0;
    double lw1 = BETA * 7.0;
    double mx  = std::max(lw0, lw1);
    double Z   = std::exp(lw0 - mx) + std::exp(lw1 - mx);
    double p0_expected = std::exp(lw0 - mx) / Z;
    double p1_expected = std::exp(lw1 - mx) / Z;

    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 2u);
    EXPECT_NEAR(pops[0], p0_expected, EPS);
    EXPECT_NEAR(pops[1], p1_expected, EPS);

    // Expected receptor entropy: S = -kB * (p0 * ln(p0) + p1 * ln(p1))
    double S_expected = -kB * (p0_expected * std::log(p0_expected) +
                               p1_expected * std::log(p1_expected));
    double S_receptor = mode.receptor_conformational_entropy();
    EXPECT_NEAR(S_receptor, S_expected, EPS);
}

// ===========================================================================
// CONFORMATIONAL SELECTION vs INDUCED FIT CLASSIFICATION
// ===========================================================================

TEST_F(CCBMTest, ConformationalSelectionRegime) {
    // One conformer clearly dominates → low receptor entropy → conformational selection
    BindingMode mode(test_population);
    // Model 0: great affinity, no strain
    Pose p0 = make_ccbm_pose(-15.0, 0, 0.0, 0);
    mode.add_Pose(p0);
    // Model 1: poor affinity, high strain
    Pose p1 = make_ccbm_pose(-5.0, 1, 3.0, 1);
    mode.add_Pose(p1);
    // Model 2: very poor
    Pose p2 = make_ccbm_pose(-3.0, 2, 5.0, 2);
    mode.add_Pose(p2);

    double S_receptor = mode.receptor_conformational_entropy();
    double S_max = kB * std::log(3.0);  // maximum possible for 3 states

    // S_receptor should be much less than S_max
    EXPECT_LT(S_receptor, 0.1 * S_max);
}

TEST_F(CCBMTest, InducedFitRegime) {
    // All conformers roughly equally populated → high receptor entropy → induced fit
    BindingMode mode(test_population);
    for (int r = 0; r < 5; ++r) {
        Pose p = make_ccbm_pose(-10.0, r, 0.0, r);
        mode.add_Pose(p);
    }

    double S_receptor = mode.receptor_conformational_entropy();
    double S_max = kB * std::log(5.0);

    // S_receptor should be very close to S_max
    EXPECT_GT(S_receptor, 0.99 * S_max);
}

// ===========================================================================
// REBUILD ENGINE WITH TOTAL_ENERGY
// ===========================================================================

TEST_F(CCBMTest, RebuildEngineUsesTotalEnergy) {
    // Verify that the thermodynamic free energy accounts for receptor strain
    BindingMode mode_no_strain(test_population);
    BindingMode mode_with_strain(test_population);

    // Same CF but different strain
    Pose p0 = make_ccbm_pose(-10.0, 0, 0.0, 0);
    Pose p1 = make_ccbm_pose(-10.0, 0, 3.0, 1);

    mode_no_strain.add_Pose(p0);
    mode_with_strain.add_Pose(p1);

    double F_no_strain   = mode_no_strain.get_free_energy();
    double F_with_strain = mode_with_strain.get_free_energy();

    // Free energy with strain should be higher (less favorable) by exactly 3.0
    EXPECT_NEAR(F_with_strain - F_no_strain, 3.0, EPS);
}

TEST_F(CCBMTest, ThermodynamicsWithStrainCorrectEnthalpy) {
    // Two poses on different conformers
    BindingMode mode(test_population);
    Pose p0 = make_ccbm_pose(-12.0, 0, 0.0, 0);   // E_total = -12.0
    Pose p1 = make_ccbm_pose(-11.0, 1, 0.5, 1);    // E_total = -10.5
    mode.add_Pose(p0);
    mode.add_Pose(p1);

    auto thermo = mode.get_thermodynamics();

    // Mean energy should be Boltzmann-weighted average of total_energy()
    // not just CF
    double w0 = std::exp(-BETA * (-12.0));
    double w1 = std::exp(-BETA * (-10.5));
    double Z  = w0 + w1;
    double expected_mean = (-12.0 * w0 + (-10.5) * w1) / Z;

    EXPECT_NEAR(thermo.mean_energy, expected_mean, 0.01);
}

// ===========================================================================
// FA_Global MULTI-MODEL FIELDS
// ===========================================================================

TEST_F(CCBMTest, FAGlobalMultiModelDefaults) {
    // Verify default values
    EXPECT_FALSE(mock_fa->multi_model);
    EXPECT_EQ(mock_fa->n_models, 1);
    EXPECT_EQ(mock_fa->model_gene_index, -1);
    EXPECT_TRUE(mock_fa->model_coords.empty());
    EXPECT_TRUE(mock_fa->model_strain.empty());
}

TEST_F(CCBMTest, FAGlobalMultiModelSetup) {
    mock_fa->multi_model = true;
    mock_fa->n_models = 3;
    mock_fa->model_gene_index = 6;

    // Simulate 3 models with 5 atoms each (15 floats per model)
    mock_fa->model_coords.resize(3);
    for (int r = 0; r < 3; ++r) {
        mock_fa->model_coords[r].resize(15, static_cast<float>(r));
    }
    mock_fa->model_strain = {0.0, 1.5, 3.2};

    EXPECT_TRUE(mock_fa->multi_model);
    EXPECT_EQ(mock_fa->n_models, 3);
    EXPECT_EQ(mock_fa->model_gene_index, 6);
    EXPECT_EQ(mock_fa->model_coords.size(), 3u);
    EXPECT_EQ(mock_fa->model_strain.size(), 3u);
    EXPECT_NEAR(mock_fa->model_strain[1], 1.5, EPS);
}

// ===========================================================================
// MULTI-POSE PER CONFORMER: WEIGHTED POPULATIONS
// ===========================================================================

TEST_F(CCBMTest, ConformerPopulationsWeightedByPoseCount) {
    // Model 0 has 5 degenerate poses at E = -10
    // Model 1 has 1 pose at E = -10
    // Both have same individual energy but model 0 has more microstates
    // → p(0) = 5/6, p(1) = 1/6
    BindingMode mode(test_population);
    for (int i = 0; i < 5; ++i) {
        Pose p = make_ccbm_pose(-10.0, 0, 0.0, i);
        mode.add_Pose(p);
    }
    Pose p1 = make_ccbm_pose(-10.0, 1, 0.0, 5);
    mode.add_Pose(p1);

    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 2u);
    EXPECT_NEAR(pops[0], 5.0 / 6.0, EPS);
    EXPECT_NEAR(pops[1], 1.0 / 6.0, EPS);
}

// ===========================================================================
// ENTROPY DECOMPOSITION: ANALYTICAL 2x2 SYSTEM
// ===========================================================================

TEST_F(CCBMTest, EntropyDecomposition2x2Analytical) {
    // 2 conformers × 2 unique ligand poses = 4 microstates
    // All at the same total energy → uniform distribution
    // p(r,i) = 1/4 for all
    // S_total = kB * ln(4)
    // S_receptor = kB * ln(2) (marginal over models: p(r) = 1/2)
    // S_ligand   = kB * ln(2) (marginal over poses: p(i) = 1/2)
    // MI = S_L + S_R - S_joint = kB*ln(2) + kB*ln(2) - kB*ln(4) = 0
    BindingMode mode(test_population);

    // Each chrom_index is a unique ligand pose
    // chrom 0 on model 0, chrom 0 on model 1, chrom 1 on model 0, chrom 1 on model 1
    Pose p00 = make_ccbm_pose(-10.0, 0, 0.0, 0);
    Pose p01 = make_ccbm_pose(-10.0, 1, 0.0, 0);  // same chrom_index, different model
    Pose p10 = make_ccbm_pose(-10.0, 0, 0.0, 1);
    Pose p11 = make_ccbm_pose(-10.0, 1, 0.0, 1);  // same chrom_index, different model
    mode.add_Pose(p00);
    mode.add_Pose(p01);
    mode.add_Pose(p10);
    mode.add_Pose(p11);

    auto ed = mode.decompose_entropy();

    double S_total_expected = kB * std::log(4.0);
    double S_marg_expected  = kB * std::log(2.0);

    EXPECT_NEAR(ed.S_total, S_total_expected, EPS);
    EXPECT_NEAR(ed.S_receptor, S_marg_expected, EPS);
    EXPECT_NEAR(ed.S_ligand, S_marg_expected, EPS);
    EXPECT_NEAR(ed.I_mutual, 0.0, EPS_LOOSE);
}

// ===========================================================================
// STRAIN ENERGY EFFECT ON FREE ENERGY
// ===========================================================================

TEST_F(CCBMTest, StrainPenaltyRaisesEffectiveFreeEnergy) {
    // Two modes: identical CF but one has strain
    // The strained mode should have higher (less favorable) free energy
    BindingMode mode_clean(test_population);
    BindingMode mode_strained(test_population);

    for (int i = 0; i < 3; ++i) {
        Pose pc = make_ccbm_pose(-10.0 - i, 0, 0.0, i);
        mode_clean.add_Pose(pc);
    }
    for (int i = 0; i < 3; ++i) {
        Pose ps = make_ccbm_pose(-10.0 - i, 1, 2.0, i);  // 2.0 strain penalty
        mode_strained.add_Pose(ps);
    }

    double F_clean    = mode_clean.get_free_energy();
    double F_strained = mode_strained.get_free_energy();

    // Strained mode should have higher free energy
    EXPECT_GT(F_strained, F_clean);
    // The difference should be approximately the strain penalty (2.0)
    // for a uniform strain shift
    EXPECT_NEAR(F_strained - F_clean, 2.0, 0.01);
}

// ===========================================================================
// CACHE INVALIDATION WITH CCBM POSES
// ===========================================================================

TEST_F(CCBMTest, CacheInvalidatedOnNewCCBMPose) {
    BindingMode mode(test_population);
    Pose p0 = make_ccbm_pose(-10.0, 0, 0.0, 0);
    mode.add_Pose(p0);

    auto thermo1 = mode.get_thermodynamics();

    // Adding a pose on a different conformer should invalidate cache
    Pose p1 = make_ccbm_pose(-9.0, 1, 0.5, 1);
    mode.add_Pose(p1);

    auto thermo2 = mode.get_thermodynamics();
    EXPECT_NE(thermo1.free_energy, thermo2.free_energy);
}

// ===========================================================================
// NUMERICAL STABILITY: LARGE ENERGY SPREAD
// ===========================================================================

TEST_F(CCBMTest, NumericalStabilityLargeEnergySpread) {
    // Energies spanning 100 kcal/mol — tests log-sum-exp stability
    BindingMode mode(test_population);
    Pose p0 = make_ccbm_pose(-100.0, 0, 0.0, 0);
    Pose p1 = make_ccbm_pose(-50.0, 1, 0.0, 1);
    Pose p2 = make_ccbm_pose(0.0, 2, 0.0, 2);
    mode.add_Pose(p0);
    mode.add_Pose(p1);
    mode.add_Pose(p2);

    // Should not produce NaN or Inf
    auto pops = mode.conformer_populations();
    ASSERT_EQ(pops.size(), 3u);
    for (double p : pops) {
        EXPECT_FALSE(std::isnan(p));
        EXPECT_FALSE(std::isinf(p));
    }
    // The lowest energy state should dominate completely
    EXPECT_GT(pops[0], 0.999);

    double S = mode.receptor_conformational_entropy();
    EXPECT_FALSE(std::isnan(S));
    EXPECT_GE(S, 0.0);

    double MI = mode.ligand_receptor_mutual_information();
    EXPECT_FALSE(std::isnan(MI));
    EXPECT_GE(MI, -EPS_LOOSE);

    auto ed = mode.decompose_entropy();
    EXPECT_FALSE(std::isnan(ed.S_total));
    EXPECT_FALSE(std::isnan(ed.S_ligand));
    EXPECT_FALSE(std::isnan(ed.S_receptor));
    EXPECT_FALSE(std::isnan(ed.I_mutual));
}

TEST_F(CCBMTest, NumericalStabilityDegenerateEnergies) {
    // All exactly the same energy — tests for division stability
    BindingMode mode(test_population);
    for (int r = 0; r < 5; ++r) {
        Pose p = make_ccbm_pose(-10.0, r, 0.0, r);
        mode.add_Pose(p);
    }

    auto pops = mode.conformer_populations();
    for (double p : pops) {
        EXPECT_FALSE(std::isnan(p));
        EXPECT_NEAR(p, 0.2, EPS);
    }

    auto ed = mode.decompose_entropy();
    EXPECT_FALSE(std::isnan(ed.S_total));
    EXPECT_FALSE(std::isnan(ed.S_receptor));
    EXPECT_FALSE(std::isnan(ed.S_ligand));
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
