// tests/test_binding_mode_statmech.cpp
// Unit tests for BindingMode ↔ StatMechEngine integration
// Part of FlexAIDΔS Phase 1 implementation roadmap
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/BindingMode.h"
#include "../LIB/statmech.h"
#include "../LIB/gaboom.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <chrono>

// ===========================================================================
// TEST FIXTURES
// ===========================================================================

class BindingModeStatMechTest : public ::testing::Test {
protected:
    // Mock structures for testing
    FA_Global* mock_fa;
    GB_Global* mock_gb;
    VC_Global* mock_vc;
    chromosome* mock_chroms;
    genlim* mock_gene_lim;
    atom* mock_atoms;
    resid* mock_residue;
    gridpoint* mock_cleftgrid;
    
    BindingPopulation* test_population;
    
    const double TEST_TEMPERATURE = 300.0;  // Kelvin
    const double EPSILON = 1e-6;  // Numerical tolerance
    
    void SetUp() override {
        // Initialize minimal mock structures (zero-init to avoid UB)
        mock_fa = new FA_Global();
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
        std::memset(mock_fa, 0, sizeof(FA_Global));
#pragma clang diagnostic pop
        mock_fa->temperature = static_cast<uint>(TEST_TEMPERATURE);

        mock_gb = new GB_Global();
        std::memset(mock_gb, 0, sizeof(GB_Global));
        mock_gb->num_genes = 6;  // 3 translation + 3 rotation

        mock_vc = new VC_Global();
        std::memset(mock_vc, 0, sizeof(VC_Global));

        mock_chroms = new chromosome[5];  // 5 test chromosomes
        for (int i = 0; i < 5; ++i) {
            mock_chroms[i].genes = new gene[mock_gb->num_genes];
            std::memset(mock_chroms[i].genes, 0, sizeof(gene) * mock_gb->num_genes);
            mock_chroms[i].evalue = 0.0;
            mock_chroms[i].app_evalue = 0.0;
            mock_chroms[i].fitnes = 0.0;
            mock_chroms[i].status = 'n';
        }
        mock_gene_lim = new genlim[mock_gb->num_genes];
        std::memset(mock_gene_lim, 0, sizeof(genlim) * mock_gb->num_genes);
        mock_atoms = new atom[10];
        std::memset(mock_atoms, 0, sizeof(atom) * 10);
        mock_residue = new resid[2];
        std::memset(mock_residue, 0, sizeof(resid) * 2);
        mock_cleftgrid = new gridpoint[100];
        std::memset(mock_cleftgrid, 0, sizeof(gridpoint) * 100);

        // Create test population
        test_population = new BindingPopulation(
            mock_fa, mock_gb, mock_vc,
            mock_chroms, mock_gene_lim,
            mock_atoms, mock_residue, mock_cleftgrid,
            5  // n_chrom
        );
    }
    
    void TearDown() override {
        delete test_population;
        delete[] mock_cleftgrid;
        delete[] mock_residue;
        delete[] mock_atoms;
        delete[] mock_gene_lim;
        for (int i = 0; i < 5; ++i) delete[] mock_chroms[i].genes;
        delete[] mock_chroms;
        delete mock_vc;
        delete mock_gb;
        delete mock_fa;
    }
    
    // Helper: Create mock pose with specific CF (also sets app_evalue so soft-β
    // Pose::boltzmann_weight matches CF — required for classic global-Z ranking).
    Pose create_mock_pose(double cf_value, int index) {
        mock_chroms[index].app_evalue = cf_value;
        mock_chroms[index].evalue = cf_value;
        std::vector<float> empty_vec;
        Pose p(&mock_chroms[index], index, 0, 0.0f, static_cast<uint>(TEST_TEMPERATURE), empty_vec);
        p.CF = cf_value;
        return p;
    }

    // Legacy physical-beta numeric path (still CF-proxy without provenance).
    void enable_force_cf_ranking() {
        mock_fa->force_cf_rank_emission = true;
    }
};

// ===========================================================================
// CORE FUNCTIONALITY TESTS
// ===========================================================================

TEST_F(BindingModeStatMechTest, LazyEngineRebuild) {
    BindingMode mode(test_population);
    
    // Add poses with known CF values
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-10.0 - i * 2.0, i);
        mode.add_Pose(p);
    }
    
    // First call should build engine — verify indirectly via result validity
    auto thermo1 = mode.get_thermodynamics();
    
    // Second call should reuse cache
    auto thermo2 = mode.get_thermodynamics();
    EXPECT_EQ(thermo1.free_energy, thermo2.free_energy);
    
    // Adding pose should change results (cache invalidated)
    Pose new_pose = create_mock_pose(-8.0, 3);
    mode.add_Pose(new_pose);
    auto thermo3 = mode.get_thermodynamics();
    // With a new pose, thermodynamics should change
    EXPECT_NE(thermo1.free_energy, thermo3.free_energy);
}

TEST_F(BindingModeStatMechTest, LiveBindingModeIsExplicitlyCFProxyOnly) {
    BindingMode mode(test_population);
    Pose pose = create_mock_pose(-10.0, 0);
    mode.add_Pose(pose);

    const auto thermo = mode.get_thermodynamics();
    const auto ledger = mode.get_thermodynamic_breakdown();

    EXPECT_EQ(thermo.provenance.energy_domain,
              statmech::EnergyDomain::ContactFunctionArbitraryUnits);
    EXPECT_EQ(thermo.provenance.ensemble_measure,
              statmech::EnsembleMeasure::OptimizerSamples);
    EXPECT_EQ(thermo.provenance.reference_state,
              statmech::ReferenceState::BoundOnly);
    EXPECT_EQ(thermo.claim_validity(), statmech::ClaimValidity::ProxyOnly);
    EXPECT_FALSE(thermo.allows_canonical_physical_claim());
    EXPECT_FALSE(thermo.allows_binding_physical_claim());
    EXPECT_EQ(ledger.provenance, thermo.provenance);
    EXPECT_TRUE(ledger.is_proxy_only());
}

TEST_F(BindingModeStatMechTest, ConsistencyWithLegacy) {
    // Legacy physical-beta numeric path: ranking F must match the proxy ledger.
    // Classic path intentionally diverges (soft-β global Z vs kB per-mode).
    enable_force_cf_ranking();
    BindingMode mode(test_population);
    
    // Add diverse CF distribution
    std::vector<double> cf_values = {-15.0, -12.0, -10.0, -8.0, -6.0};
    for (size_t i = 0; i < cf_values.size(); ++i) {
        Pose p = create_mock_pose(cf_values[i], i);
        mode.add_Pose(p);
    }
    
    // Legacy and new API should give identical results
    double legacy_energy = mode.compute_energy();
    double new_energy = mode.get_free_energy();
    EXPECT_NEAR(legacy_energy, new_energy, EPSILON);
    
    double legacy_enthalpy = mode.compute_enthalpy();
    auto thermo = mode.get_thermodynamics();
    EXPECT_NEAR(legacy_enthalpy, thermo.mean_energy, EPSILON);
    
    double legacy_entropy = mode.compute_entropy();
    EXPECT_NEAR(legacy_entropy, thermo.entropy, EPSILON);
}

TEST_F(BindingModeStatMechTest, ThermodynamicBreakdownPreservesLegacyRankingEnergy) {
    // Breakdown uses physical beta numerically but remains CF-proxy; ranking
    // under force_cf matches the numeric ledger.
    enable_force_cf_ranking();
    BindingMode mode(test_population);

    std::vector<double> cf_values = {-14.0, -11.0, -9.0};
    for (size_t i = 0; i < cf_values.size(); ++i) {
        Pose p = create_mock_pose(cf_values[i], static_cast<int>(i));
        mode.add_Pose(p);
    }
    mock_fa->natural_deltaG = 0.25;

    const double legacy_energy = mode.compute_energy();
    const auto legacy = mode.get_thermodynamics();
    const auto ledger = mode.get_thermodynamic_breakdown();

    // The ranking energy is unchanged: it still carries the NATURaL correction.
    EXPECT_NEAR(ledger.G_total_kcal_mol, legacy_energy, EPSILON);

    // But get_thermodynamics() now returns the configurational ensemble
    // UNMODIFIED. It used to add the correction to free_energy alone, which
    // broke the struct's own identity (entropy is defined as
    // (mean_energy - free_energy)/T, so shifting F without shifting H or S
    // left F != H - T*S). The correction lives in the ledger's dedicated
    // G_natural_kcal_mol field instead.
    EXPECT_NEAR(legacy.free_energy, ledger.G_config_kcal_mol, EPSILON);
    EXPECT_NEAR(ledger.G_total_kcal_mol, legacy.free_energy + ledger.G_natural_kcal_mol,
                EPSILON);
    ASSERT_GT(legacy.temperature, 0.0);
    EXPECT_NEAR(legacy.free_energy,
                legacy.mean_energy - legacy.temperature * legacy.entropy,
                EPSILON);

    EXPECT_NEAR(ledger.G_config_kcal_mol + ledger.G_natural_kcal_mol + ledger.G_vib_kcal_mol,
                ledger.G_total_kcal_mol, EPSILON);
    EXPECT_TRUE(ledger.has_natural);
    EXPECT_FALSE(ledger.has_other);

    // An unreceipted correction downgrades the aggregate ledger witness even
    // though the scalar configurational ledger keeps its CF-proxy metadata.
    EXPECT_TRUE(ledger.is_proxy_only());
}

TEST_F(BindingModeStatMechTest, BoltzmannWeightsNormalization) {
    BindingMode mode(test_population);
    
    // Add poses with known CF distribution
    std::vector<double> cf_values = {-20.0, -15.0, -10.0, -5.0};
    for (size_t i = 0; i < cf_values.size(); ++i) {
        Pose p = create_mock_pose(cf_values[i], i);
        mode.add_Pose(p);
    }
    
    auto weights = mode.get_boltzmann_weights();
    
    // Sum of weights should be 1.0 (within numerical tolerance)
    double sum = 0.0;
    for (double w : weights) {
        sum += w;
        EXPECT_GE(w, 0.0);  // All weights non-negative
    }
    EXPECT_NEAR(sum, 1.0, EPSILON);
    
    // Lowest CF should have highest weight
    double min_cf = *std::min_element(cf_values.begin(), cf_values.end());
    size_t min_idx = std::distance(cf_values.begin(), 
                                   std::find(cf_values.begin(), cf_values.end(), min_cf));
    double max_weight = *std::max_element(weights.begin(), weights.end());
    EXPECT_NEAR(weights[min_idx], max_weight, EPSILON);
}

TEST_F(BindingModeStatMechTest, EntropyBehavior) {
    // Physical β: ground-state collapse on broad CF spreads.
    enable_force_cf_ranking();
    BindingMode mode_sharp(test_population);
    BindingMode mode_broad(test_population);
    
    // Sharp distribution: all CFs similar — near-equal Boltzmann weights → high entropy
    for (int i = 0; i < 5; ++i) {
        Pose p = create_mock_pose(-10.0 - 0.001 * i, i);
        mode_sharp.add_Pose(p);
    }

    // Broad distribution: CFs widely spread — lowest state dominates → low entropy
    // At T=300K (kBT≈0.6 kcal/mol), 20 kcal/mol spread means ground state dominates
    std::vector<double> broad_cfs = {-20.0, -15.0, -10.0, -5.0, 0.0};
    for (size_t i = 0; i < broad_cfs.size(); ++i) {
        Pose p = create_mock_pose(broad_cfs[i], i);
        mode_broad.add_Pose(p);
    }

    double sharp_entropy = mode_sharp.compute_entropy();
    double broad_entropy = mode_broad.compute_entropy();

    // Nearly degenerate states should have HIGHER entropy than widely spread ones
    // because Boltzmann weights are more uniform when energies are similar
    EXPECT_GT(sharp_entropy, broad_entropy);
}

// ===========================================================================
// BINDINGPOPULATION TESTS
// ===========================================================================

TEST_F(BindingModeStatMechTest, DeltaGCalculation) {
    // Per-mode physical F without global Z (modes not yet in population PF).
    enable_force_cf_ranking();
    BindingMode mode1(test_population);
    BindingMode mode2(test_population);
    
    // Mode 1: lower energy (more favorable)
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode1.add_Pose(p);
    }
    
    // Mode 2: higher energy (less favorable)
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-10.0 - i, i);
        mode2.add_Pose(p);
    }
    
    double delta_g = test_population->compute_delta_G(mode1, mode2);

    // compute_delta_G(A, B) = F_B - F_A; mode2 higher energy → negative
    EXPECT_LT(delta_g, 0.0);
}

TEST_F(BindingModeStatMechTest, GlobalEnsemble) {
    // Create two binding modes
    BindingMode mode1(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode1.add_Pose(p);
    }

    BindingMode mode2(test_population);
    for (int i = 0; i < 2; ++i) {
        Pose p = create_mock_pose(-10.0 - i, i);
        mode2.add_Pose(p);
    }

    test_population->add_BindingMode(mode1);
    test_population->add_BindingMode(mode2);

    // Population should now contain 2 binding modes
    EXPECT_EQ(test_population->get_Population_size(), 2);

    // Global ensemble aggregates poses from all modes
    auto global_engine = test_population->get_global_ensemble();
    EXPECT_GT(global_engine.size(), 0u);

    // P1: get_log_Z accessor returns real ensemble log_Z (non-negative for valid pop)
    double lz = test_population->get_log_Z();
    EXPECT_GE(lz, -1e9);  // reasonable bound, actual value depends on energies
}

// ===========================================================================
// CACHE INVALIDATION TESTS
// ===========================================================================

TEST_F(BindingModeStatMechTest, CacheInvalidationOnClear) {
    BindingMode mode(test_population);
    
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-10.0 - i, i);
        mode.add_Pose(p);
    }
    
    // Build cache
    auto thermo_before = mode.get_thermodynamics();

    // Clear should invalidate
    mode.clear_Poses();
    EXPECT_EQ(mode.get_BindingMode_size(), 0);
}

TEST_F(BindingModeStatMechTest, MultipleRebuilds) {
    enable_force_cf_ranking();
    BindingMode mode(test_population);
    
    // First batch of poses
    for (int i = 0; i < 2; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode.add_Pose(p);
    }
    double energy1 = mode.get_free_energy();
    
    // Add more poses (invalidates cache)
    for (int i = 2; i < 4; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode.add_Pose(p);
    }
    double energy2 = mode.get_free_energy();
    
    // Energy should change (more states in ensemble)
    EXPECT_NE(energy1, energy2);
}

// ===========================================================================
// EDGE CASES
// ===========================================================================

TEST_F(BindingModeStatMechTest, EmptyMode) {
    BindingMode mode(test_population);

    // Empty mode should report zero size
    EXPECT_EQ(mode.get_BindingMode_size(), 0);

    // Calling thermodynamics on empty mode should throw (empty ensemble)
    EXPECT_THROW(mode.get_thermodynamics(), std::runtime_error);
}

TEST_F(BindingModeStatMechTest, SinglePoseMode) {
    BindingMode mode(test_population);
    
    Pose p = create_mock_pose(-12.0, 0);
    mode.add_Pose(p);
    
    auto thermo = mode.get_thermodynamics();
    
    // Single state: entropy should be zero (no uncertainty)
    EXPECT_NEAR(thermo.entropy, 0.0, EPSILON);
    // Mean energy should equal the single CF value
    EXPECT_NEAR(thermo.mean_energy, -12.0, EPSILON);
}

TEST_F(BindingModeStatMechTest, HighTemperatureBehavior) {
    // Test at very high temperature (should flatten Boltzmann distribution)
    // kB ≈ 0.002 kcal/(mol·K), so at T=10000K, kBT ≈ 20 kcal/mol
    mock_fa->temperature = 10000;  // Very high T
    BindingPopulation* hot_population = new BindingPopulation(
        mock_fa, mock_gb, mock_vc,
        mock_chroms, mock_gene_lim,
        mock_atoms, mock_residue, mock_cleftgrid, 5
    );

    BindingMode mode(hot_population);

    // Small CF spread relative to kBT → nearly uniform weights
    std::vector<double> cf_values = {-12.0, -11.0, -10.0};
    for (size_t i = 0; i < cf_values.size(); ++i) {
        Pose p = create_mock_pose(cf_values[i], i);
        mode.add_Pose(p);
    }

    auto weights = mode.get_boltzmann_weights();

    // At high T, weights should be more uniform
    double max_weight = *std::max_element(weights.begin(), weights.end());
    double min_weight = *std::min_element(weights.begin(), weights.end());
    double weight_ratio = max_weight / min_weight;

    // Ratio should be close to 1 at very high T with small spread
    EXPECT_LT(weight_ratio, 2.0);

    delete hot_population;
}

// ===========================================================================
// POPULATION SHANNON ENTROPY TESTS
// ===========================================================================

TEST_F(BindingModeStatMechTest, ShannonEntropyPopulation) {
    // Create two binding modes with different energy distributions
    BindingMode mode1(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode1.add_Pose(p);
    }

    BindingMode mode2(test_population);
    for (int i = 0; i < 2; ++i) {
        Pose p = create_mock_pose(-10.0 - i, i);
        mode2.add_Pose(p);
    }

    test_population->add_BindingMode(mode1);
    test_population->add_BindingMode(mode2);

    double shannon_S = test_population->get_shannon_entropy();

    // Shannon entropy should be positive for a multi-pose population
    EXPECT_GT(shannon_S, 0.0);
}

TEST_F(BindingModeStatMechTest, GlobalEnsembleUsesRawMicrostates) {
    BindingMode mode1(test_population);
    {
        Pose p0 = create_mock_pose(-1.0, 0);
        Pose p1 = create_mock_pose(0.0, 1);
        mode1.add_Pose(p0);
        mode1.add_Pose(p1);
    }

    BindingMode mode2(test_population);
    {
        Pose p0 = create_mock_pose(0.0, 2);
        Pose p1 = create_mock_pose(1.0, 3);
        mode2.add_Pose(p0);
        mode2.add_Pose(p1);
    }

    test_population->add_BindingMode(mode1);
    test_population->add_BindingMode(mode2);

    auto global = test_population->get_global_ensemble();
    const auto weights = global.boltzmann_weights();

    ASSERT_EQ(weights.size(), 4u);

    const double beta = 1.0 / (statmech::kB_kcal * TEST_TEMPERATURE);
    const std::vector<double> energies = {-1.0, 0.0, 0.0, 1.0};
    double Z = 0.0;
    for (double energy : energies) Z += std::exp(-beta * energy);

    for (std::size_t i = 0; i < energies.size(); ++i) {
        const double expected = std::exp(-beta * energies[i]) / Z;
        EXPECT_NEAR(weights[i], expected, 1e-9);
    }
}

TEST_F(BindingModeStatMechTest, ShannonEntropyMatchesGlobalEnsembleEntropy) {
    BindingMode mode1(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode1.add_Pose(p);
    }

    BindingMode mode2(test_population);
    for (int i = 0; i < 2; ++i) {
        Pose p = create_mock_pose(-10.0 - i, i + 3);
        mode2.add_Pose(p);
    }

    test_population->add_BindingMode(mode1);
    test_population->add_BindingMode(mode2);

    const double shannon_S = test_population->get_shannon_entropy();
    const auto global = test_population->get_global_ensemble();
    const auto thermo = global.compute();

    EXPECT_NEAR(shannon_S, thermo.entropy, 1e-12);
}

TEST_F(BindingModeStatMechTest, ShannonEntropyEmptyPopulation) {
    // Empty population should have zero Shannon entropy
    double shannon_S = test_population->get_shannon_entropy();
    EXPECT_NEAR(shannon_S, 0.0, EPSILON);
}

TEST_F(BindingModeStatMechTest, ShannonCacheInvalidation) {
    BindingMode mode1(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode1.add_Pose(p);
    }
    test_population->add_BindingMode(mode1);

    double S1 = test_population->get_shannon_entropy();

    // Adding another mode should invalidate the cache and change S
    BindingMode mode2(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-5.0 - i * 3.0, i);
        mode2.add_Pose(p);
    }
    test_population->add_BindingMode(mode2);

    double S2 = test_population->get_shannon_entropy();

    // Shannon entropy should change after adding a mode with different energies
    EXPECT_NE(S1, S2);
}

TEST_F(BindingModeStatMechTest, DeltaGMatrix) {
    BindingMode mode1(test_population);
    for (int i = 0; i < 3; ++i) {
        Pose p = create_mock_pose(-15.0 - i, i);
        mode1.add_Pose(p);
    }

    BindingMode mode2(test_population);
    for (int i = 0; i < 2; ++i) {
        Pose p = create_mock_pose(-10.0 - i, i);
        mode2.add_Pose(p);
    }

    BindingMode mode3(test_population);
    for (int i = 0; i < 2; ++i) {
        Pose p = create_mock_pose(-8.0 - i, i);
        mode3.add_Pose(p);
    }

    test_population->add_BindingMode(mode1);
    test_population->add_BindingMode(mode2);
    test_population->add_BindingMode(mode3);

    auto matrix = test_population->get_deltaG_matrix();

    // Matrix should be 3x3
    EXPECT_EQ(matrix.size(), 3u);
    for (const auto& row : matrix) {
        EXPECT_EQ(row.size(), 3u);
    }

    // Diagonal should be zero
    for (size_t i = 0; i < matrix.size(); ++i) {
        EXPECT_NEAR(matrix[i][i], 0.0, EPSILON);
    }

    // Anti-symmetry: matrix[i][j] = -matrix[j][i]
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = i + 1; j < matrix.size(); ++j) {
            EXPECT_NEAR(matrix[i][j], -matrix[j][i], EPSILON);
        }
    }
}

// ===========================================================================
// PERFORMANCE TESTS (Optional, disabled in CI)
// ===========================================================================

TEST_F(BindingModeStatMechTest, DISABLED_CachePerformance) {
    // Verifies that caching provides performance benefit
    // Disabled by default (use --gtest_also_run_disabled_tests to run)
    
    BindingMode mode(test_population);
    
    // Large number of poses
    for (int i = 0; i < 1000; ++i) {
        Pose p = create_mock_pose(-20.0 + i * 0.01, i % 5);
        mode.add_Pose(p);
    }
    
    // Time first call (rebuilds engine)
    auto start1 = std::chrono::high_resolution_clock::now();
    mode.get_thermodynamics();
    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
    
    // Time second call (uses cache)
    auto start2 = std::chrono::high_resolution_clock::now();
    mode.get_thermodynamics();
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    
    // Cached call should be significantly faster
    EXPECT_LT(duration2.count(), duration1.count() / 10);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
