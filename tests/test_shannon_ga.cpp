/**
 * @file test_shannon_ga.cpp
 * @brief Unit tests for ShannonThermodynamicGA
 */
#include <gtest/gtest.h>
#include "../LIB/shannon_ga.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <memory>

using namespace FlexAID::StatMech;

// ─── helper: allocate a minimal chromosome population ────────────────────────
struct TestPopulation {
    std::vector<chromosome> chroms;
    std::vector<std::vector<gene>> gene_storage;
    int num_genes;

    TestPopulation(int pop_size, int n_genes, const std::vector<double>& energies)
        : num_genes(n_genes)
    {
        chroms.resize(pop_size);
        gene_storage.resize(pop_size, std::vector<gene>(n_genes));

        for (int i = 0; i < pop_size; ++i) {
            chroms[i].genes = gene_storage[i].data();
            chroms[i].evalue = (i < static_cast<int>(energies.size()))
                               ? energies[i] : 0.0;
            chroms[i].app_evalue = chroms[i].evalue;
            chroms[i].fitnes = 0.0;
            chroms[i].boltzmann_weight = 0.0;
            chroms[i].free_energy = 0.0;
            chroms[i].status = 'n';
            std::memset(&chroms[i].cf, 0, sizeof(cfstr));

            // Fill genes with distinct values
            for (int g = 0; g < n_genes; ++g) {
                gene_storage[i][g].to_int32 = static_cast<int32_t>(i * 1000 + g);
                gene_storage[i][g].to_ic = static_cast<double>(i * 10 + g);
            }
        }
    }

    chromosome* data() { return chroms.data(); }
    int size() const { return static_cast<int>(chroms.size()); }
};

// ─── tests ───────────────────────────────────────────────────────────────────

TEST(ShannonThermodynamicGA, ConstructorDefaults) {
    ShannonThermodynamicGA ga(100, 0.01);
    EXPECT_DOUBLE_EQ(ga.temperature(), 298.15);
    EXPECT_DOUBLE_EQ(ga.current_entropy(), 0.0);
    EXPECT_DOUBLE_EQ(ga.current_partition_Z(), 0.0);
}

TEST(ShannonThermodynamicGA, ConstructorCustomTemperature) {
    ShannonThermodynamicGA ga(50, 0.05, 310.0);
    EXPECT_DOUBLE_EQ(ga.temperature(), 310.0);
}

TEST(ShannonThermodynamicGA, EvaluateEnthalpyBatch) {
    std::vector<double> energies = {-10.0, -8.0, -6.0, -4.0, -2.0};
    TestPopulation pop(5, 2, energies);

    ShannonThermodynamicGA ga(5, 0.01);
    ga.evaluate_enthalpy_batch(pop.data(), pop.size());

    // After evaluate, the engine should have 5 samples
    // We verify indirectly via collapse_thermodynamics
    auto snap = ga.collapse_thermodynamics(pop.data(), pop.size(), 1);
    EXPECT_EQ(snap.generation, 1);
    // Mean energy should be between -10 and -2 (Boltzmann-weighted toward -10)
    EXPECT_LT(snap.mean_energy, -2.0);
    EXPECT_GT(snap.mean_energy, -10.0);
}

TEST(ShannonThermodynamicGA, CollapseThermodynamics_BoltzmannWeights) {
    // Diverse energies
    std::vector<double> energies = {-20.0, -15.0, -10.0, -5.0, 0.0};
    TestPopulation pop(5, 2, energies);

    ShannonThermodynamicGA ga(5, 0.01);
    ga.evaluate_enthalpy_batch(pop.data(), pop.size());
    auto snap = ga.collapse_thermodynamics(pop.data(), pop.size(), 1);

    // Best energy (-20) should have highest Boltzmann weight
    double max_bw = -1.0;
    int max_idx = -1;
    for (int i = 0; i < pop.size(); ++i) {
        if (pop.chroms[i].boltzmann_weight > max_bw) {
            max_bw = pop.chroms[i].boltzmann_weight;
            max_idx = i;
        }
    }
    // The chromosome with energy -20 should dominate
    EXPECT_EQ(pop.chroms[max_idx].evalue, -20.0);

    // Free energy should be set on all chromosomes
    for (int i = 0; i < pop.size(); ++i) {
        EXPECT_NE(pop.chroms[i].free_energy, 0.0);
    }
}

TEST(ShannonThermodynamicGA, CollapseThermodynamics_ShannonEntropy) {
    // Uniform-ish energies → higher Shannon entropy
    std::vector<double> uniform = {-5.0, -5.1, -4.9, -5.2, -4.8,
                                    -5.05, -4.95, -5.15, -4.85, -5.0};
    TestPopulation pop_uniform(10, 2, uniform);

    ShannonThermodynamicGA ga1(10, 0.01);
    ga1.evaluate_enthalpy_batch(pop_uniform.data(), pop_uniform.size());
    auto snap1 = ga1.collapse_thermodynamics(pop_uniform.data(),
                                              pop_uniform.size(), 1);

    // Collapsed energies → lower Shannon entropy
    std::vector<double> collapsed = {-10.0, -10.0, -10.0, -10.0, -10.0,
                                      -10.0, -10.0, -10.0, -10.0, -10.0};
    TestPopulation pop_collapsed(10, 2, collapsed);

    ShannonThermodynamicGA ga2(10, 0.01);
    ga2.evaluate_enthalpy_batch(pop_collapsed.data(), pop_collapsed.size());
    auto snap2 = ga2.collapse_thermodynamics(pop_collapsed.data(),
                                              pop_collapsed.size(), 1);

    // Uniform population should have higher Shannon entropy than collapsed
    EXPECT_GT(snap1.shannon_H, snap2.shannon_H);
}

TEST(ShannonThermodynamicGA, ShannonEnergyCollapseThresholdsAreNats) {
    EXPECT_NEAR(shannon_thermo::kHSC_soft_nats,
                2.0 * std::log(2.0),
                1e-12);
    EXPECT_NEAR(shannon_thermo::kHSC_hard_nats,
                std::log(2.0),
                1e-12);

    std::vector<double> collapsed = {-10.0, -10.0, -10.0, -10.0, -10.0};
    TestPopulation pop(5, 2, collapsed);

    ShannonThermodynamicGA ga(5, 0.01);
    ga.evaluate_enthalpy_batch(pop.data(), pop.size());
    auto snap = ga.collapse_thermodynamics(pop.data(), pop.size(), 1);

    EXPECT_LT(snap.shannon_H, shannon_thermo::kHSC_soft_nats);
}

TEST(ShannonThermodynamicGA, SelectParent_PrefersBetterEnergies) {
    std::vector<double> energies = {-20.0, -1.0, -1.0, -1.0, -1.0};
    TestPopulation pop(5, 2, energies);

    ShannonThermodynamicGA ga(5, 0.01);
    ga.evaluate_enthalpy_batch(pop.data(), pop.size());
    ga.collapse_thermodynamics(pop.data(), pop.size(), 1);

    // Select many parents — the one with E=-20 should be selected most often
    int counts[5] = {};
    for (int trial = 0; trial < 1000; ++trial) {
        int idx = ga.select_parent(pop.data(), pop.size());
        ASSERT_GE(idx, 0);
        ASSERT_LT(idx, 5);
        counts[idx]++;
    }
    // Index 0 (E=-20) should dominate at 298 K with 19 kcal/mol gap
    EXPECT_GT(counts[0], 500);
}

TEST(ShannonThermodynamicGA, FreeEnergy_IsNegative) {
    std::vector<double> energies = {-10.0, -8.5, -7.0, -6.0, -5.0};
    TestPopulation pop(5, 2, energies);

    ShannonThermodynamicGA ga(5, 0.01);
    ga.evaluate_enthalpy_batch(pop.data(), pop.size());
    auto snap = ga.collapse_thermodynamics(pop.data(), pop.size(), 1);

    // Free energy should be negative (favorable binding)
    EXPECT_LT(snap.free_energy, 0.0);
    // Free energy <= mean energy (since -TS <= 0 at finite T with S >= 0)
    EXPECT_LE(snap.free_energy, snap.mean_energy + 1e-10);
}

TEST(ShannonThermodynamicGA, GenerationThermoFields) {
    std::vector<double> energies = {-12.0, -10.0, -8.0, -6.0};
    TestPopulation pop(4, 2, energies);

    ShannonThermodynamicGA ga(4, 0.01);
    ga.evaluate_enthalpy_batch(pop.data(), pop.size());
    auto snap = ga.collapse_thermodynamics(pop.data(), pop.size(), 42);

    EXPECT_EQ(snap.generation, 42);
    EXPECT_GT(snap.heat_capacity, 0.0);     // C_v > 0 for non-degenerate ensemble
    EXPECT_GE(snap.entropy, 0.0);           // S >= 0
    EXPECT_GE(snap.shannon_H, 0.0);         // H >= 0
}
