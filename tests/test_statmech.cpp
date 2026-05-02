// tests/test_statmech.cpp
// Unit tests for StatMechEngine (partition function, thermodynamics, WHAM, TI)
// Part of FlexAIDΔS Phase 1 implementation roadmap
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/statmech.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <random>

using namespace statmech;

// ===========================================================================
// CONSTANTS
// ===========================================================================

static constexpr double EPSILON = 1e-6;
static constexpr double TEMPERATURE = 300.0;  // Kelvin

// ===========================================================================
// TEST FIXTURE
// ===========================================================================

class StatMechEngineTest : public ::testing::Test {
protected:
    StatMechEngine engine{TEMPERATURE};
};

// ===========================================================================
// CONSTRUCTION & BASIC STATE
// ===========================================================================

TEST_F(StatMechEngineTest, ConstructorSetsTemperature) {
    EXPECT_DOUBLE_EQ(engine.temperature(), TEMPERATURE);
    EXPECT_NEAR(engine.beta(), 1.0 / (kB_kcal * TEMPERATURE), EPSILON);
}

TEST_F(StatMechEngineTest, DefaultEngineIsEmpty) {
    EXPECT_EQ(engine.size(), 0u);
}

TEST_F(StatMechEngineTest, InvalidTemperatureThrows) {
    EXPECT_THROW(StatMechEngine(0.0), std::invalid_argument);
    EXPECT_THROW(StatMechEngine(-100.0), std::invalid_argument);
}

TEST_F(StatMechEngineTest, ComputeOnEmptyThrows) {
    EXPECT_THROW(engine.compute(), std::runtime_error);
}

TEST_F(StatMechEngineTest, AddSampleIncreasesSize) {
    engine.add_sample(-10.0);
    EXPECT_EQ(engine.size(), 1u);
    engine.add_sample(-8.0);
    EXPECT_EQ(engine.size(), 2u);
}

TEST_F(StatMechEngineTest, ClearResetsSize) {
    engine.add_sample(-10.0);
    engine.add_sample(-8.0);
    engine.clear();
    EXPECT_EQ(engine.size(), 0u);
}

// ===========================================================================
// SINGLE STATE THERMODYNAMICS
// ===========================================================================

TEST_F(StatMechEngineTest, SingleStateFreeEnergy) {
    // For a single state with energy E and multiplicity 1:
    //   Z = exp(-βE), ln Z = -βE
    //   F = -kT ln Z = E
    double E = -12.0;
    engine.add_sample(E);
    auto th = engine.compute();

    EXPECT_NEAR(th.free_energy, E, EPSILON);
    EXPECT_NEAR(th.mean_energy, E, EPSILON);
    EXPECT_NEAR(th.entropy, 0.0, EPSILON);
    EXPECT_NEAR(th.heat_capacity, 0.0, EPSILON);
    EXPECT_NEAR(th.std_energy, 0.0, EPSILON);
}

TEST_F(StatMechEngineTest, SingleStateWithMultiplicity) {
    // Single energy level with degeneracy g:
    //   Z = g * exp(-βE), ln Z = ln(g) - βE
    //   F = -kT(ln g - βE) = E - kT ln(g)
    //   ⟨E⟩ = E, S = k ln(g)
    double E = -10.0;
    int g = 5;
    engine.add_sample(E, g);
    auto th = engine.compute();

    double kT = kB_kcal * TEMPERATURE;
    double expected_F = E - kT * std::log(static_cast<double>(g));

    EXPECT_NEAR(th.free_energy, expected_F, EPSILON);
    EXPECT_NEAR(th.mean_energy, E, EPSILON);
    EXPECT_NEAR(th.entropy, kB_kcal * std::log(static_cast<double>(g)), EPSILON);
}

// ===========================================================================
// TWO-STATE SYSTEM (ANALYTICAL VERIFICATION)
// ===========================================================================

TEST_F(StatMechEngineTest, TwoStatePartitionFunction) {
    // Two states: E1 = -10, E2 = -8 (kcal/mol)
    // Z = exp(-β E1) + exp(-β E2)
    double E1 = -10.0, E2 = -8.0;
    double beta = 1.0 / (kB_kcal * TEMPERATURE);

    engine.add_sample(E1);
    engine.add_sample(E2);
    auto th = engine.compute();

    double Z = std::exp(-beta * E1) + std::exp(-beta * E2);
    double expected_F = -(kB_kcal * TEMPERATURE) * std::log(Z);
    double p1 = std::exp(-beta * E1) / Z;
    double p2 = std::exp(-beta * E2) / Z;
    double expected_E = p1 * E1 + p2 * E2;
    double expected_E2 = p1 * E1 * E1 + p2 * E2 * E2;
    double expected_var = expected_E2 - expected_E * expected_E;
    // Correct formula: C_V = Var(E) / (k_B · T²)
    // The previous expected_Cv used (k_B·T)² in the denominator which matches
    // the wrong implementation and masked the bug.
    double expected_Cv = expected_var / (kB_kcal * TEMPERATURE * TEMPERATURE);

    EXPECT_NEAR(th.free_energy, expected_F, EPSILON);
    EXPECT_NEAR(th.mean_energy, expected_E, EPSILON);
    EXPECT_NEAR(th.heat_capacity, expected_Cv, EPSILON);
    EXPECT_NEAR(th.log_Z, std::log(Z), EPSILON);
}

TEST_F(StatMechEngineTest, TwoStateBoltzmannWeights) {
    double E1 = -10.0, E2 = -8.0;
    double beta = 1.0 / (kB_kcal * TEMPERATURE);

    engine.add_sample(E1);
    engine.add_sample(E2);
    auto weights = engine.boltzmann_weights();

    ASSERT_EQ(weights.size(), 2u);

    double Z = std::exp(-beta * E1) + std::exp(-beta * E2);
    double expected_w1 = std::exp(-beta * E1) / Z;
    double expected_w2 = std::exp(-beta * E2) / Z;

    EXPECT_NEAR(weights[0], expected_w1, EPSILON);
    EXPECT_NEAR(weights[1], expected_w2, EPSILON);

    // Lower energy state should have higher weight
    EXPECT_GT(weights[0], weights[1]);

    // Weights must sum to 1
    EXPECT_NEAR(weights[0] + weights[1], 1.0, EPSILON);
}

// ===========================================================================
// BOLTZMANN WEIGHT PROPERTIES
// ===========================================================================

TEST_F(StatMechEngineTest, BoltzmannWeightsNormalization) {
    std::vector<double> energies = {-20.0, -15.0, -10.0, -5.0, 0.0, 5.0};
    for (double e : energies)
        engine.add_sample(e);

    auto weights = engine.boltzmann_weights();
    ASSERT_EQ(weights.size(), energies.size());

    double sum = 0.0;
    for (double w : weights) {
        EXPECT_GE(w, 0.0);
        sum += w;
    }
    EXPECT_NEAR(sum, 1.0, EPSILON);
}

TEST_F(StatMechEngineTest, BoltzmannWeightsOrderedByEnergy) {
    // Lower energy = higher Boltzmann weight
    std::vector<double> energies = {-20.0, -15.0, -10.0, -5.0};
    for (double e : energies)
        engine.add_sample(e);

    auto weights = engine.boltzmann_weights();
    for (size_t i = 1; i < weights.size(); ++i) {
        EXPECT_GE(weights[i - 1], weights[i])
            << "Weight at index " << i - 1 << " should be >= weight at index " << i;
    }
}

TEST_F(StatMechEngineTest, EmptyBoltzmannWeights) {
    auto weights = engine.boltzmann_weights();
    EXPECT_TRUE(weights.empty());
}

// ===========================================================================
// ENTROPY PROPERTIES
// ===========================================================================

TEST_F(StatMechEngineTest, EntropyNonNegative) {
    std::vector<double> energies = {-15.0, -12.0, -10.0, -8.0, -6.0};
    for (double e : energies)
        engine.add_sample(e);

    auto th = engine.compute();
    EXPECT_GE(th.entropy, 0.0);
}

TEST_F(StatMechEngineTest, EntropyUpperBound) {
    // S <= k_B * ln(N) for N equal-energy states
    int N = 10;
    for (int i = 0; i < N; ++i)
        engine.add_sample(-10.0);  // all same energy

    auto th = engine.compute();
    double max_entropy = kB_kcal * std::log(static_cast<double>(N));
    EXPECT_LE(th.entropy, max_entropy + EPSILON);
}

TEST_F(StatMechEngineTest, EqualEnergyStatesMaxEntropy) {
    // N states at same energy → S = k_B ln(N)
    int N = 8;
    for (int i = 0; i < N; ++i)
        engine.add_sample(-10.0);

    auto th = engine.compute();
    double expected_S = kB_kcal * std::log(static_cast<double>(N));
    EXPECT_NEAR(th.entropy, expected_S, EPSILON);
}

TEST_F(StatMechEngineTest, EntropyIncreasesWithSpread) {
    // At very high T, a broader energy spread yields entropy close to
    // the narrow case (both approach S_max = kB ln N).  At moderate T
    // the narrow (nearly degenerate) set actually has *higher* Boltzmann
    // entropy because all states remain equally accessible.
    // Test: at a high enough temperature the broad set still reaches
    // near-maximum entropy comparable to the narrow set.
    StatMechEngine narrow(100000.0);   // very high T so both are near-uniform
    StatMechEngine broad(100000.0);

    for (int i = 0; i < 5; ++i) {
        narrow.add_sample(-10.0 - 0.001 * i);  // nearly degenerate
        broad.add_sample(-10.0 - 5.0 * i);     // wide spread
    }

    auto th_narrow = narrow.compute();
    auto th_broad  = broad.compute();

    double S_max = kB_kcal * std::log(5.0);
    // Both should be close to S_max at this temperature
    EXPECT_NEAR(th_narrow.entropy, S_max, 1e-6);
    EXPECT_NEAR(th_broad.entropy,  S_max, 1e-4);
}

// ===========================================================================
// TEMPERATURE DEPENDENCE
// ===========================================================================

TEST_F(StatMechEngineTest, HighTemperatureFlattensWeights) {
    // At T → ∞, all Boltzmann weights become equal.
    // Need T high enough so β·ΔE ≪ 1.  With ΔE=30 kcal/mol and
    // kB=0.001987 kcal/(mol·K), T=1e7 gives β·ΔE ≈ 1.5e-3.
    StatMechEngine hot(10000000.0);
    std::vector<double> energies = {-20.0, -10.0, 0.0, 10.0};
    for (double e : energies)
        hot.add_sample(e);

    auto weights = hot.boltzmann_weights();
    double mean_w = 1.0 / static_cast<double>(energies.size());
    for (double w : weights)
        EXPECT_NEAR(w, mean_w, 0.03);
}

TEST_F(StatMechEngineTest, LowTemperatureConcentratesWeight) {
    // At low T, weight concentrates on lowest energy
    StatMechEngine cold(10.0);  // very low T
    cold.add_sample(-20.0);
    cold.add_sample(-10.0);
    cold.add_sample(0.0);

    auto weights = cold.boltzmann_weights();
    EXPECT_GT(weights[0], 0.99);  // nearly all weight on lowest energy
}

TEST_F(StatMechEngineTest, FreeEnergyDecreasesWithTemperature) {
    // F = E - TS, so F decreases as T increases (for S > 0)
    std::vector<double> energies = {-15.0, -10.0, -5.0};

    StatMechEngine low_T(200.0);
    StatMechEngine high_T(500.0);
    for (double e : energies) {
        low_T.add_sample(e);
        high_T.add_sample(e);
    }

    auto th_low = low_T.compute();
    auto th_high = high_T.compute();

    EXPECT_LT(th_high.free_energy, th_low.free_energy);
}

// ===========================================================================
// DELTA_G (RELATIVE FREE ENERGY)
// ===========================================================================

TEST_F(StatMechEngineTest, DeltaGSymmetry) {
    // ΔG(A→B) = -ΔG(B→A)
    StatMechEngine engine_a(TEMPERATURE);
    StatMechEngine engine_b(TEMPERATURE);

    engine_a.add_sample(-15.0);
    engine_a.add_sample(-12.0);
    engine_b.add_sample(-10.0);
    engine_b.add_sample(-8.0);

    double dG_ab = engine_a.delta_G(engine_b);
    double dG_ba = engine_b.delta_G(engine_a);

    EXPECT_NEAR(dG_ab, -dG_ba, EPSILON);
}

TEST_F(StatMechEngineTest, DeltaGSelfIsZero) {
    engine.add_sample(-10.0);
    engine.add_sample(-8.0);

    double dG = engine.delta_G(engine);
    EXPECT_NEAR(dG, 0.0, EPSILON);
}

TEST_F(StatMechEngineTest, DeltaGConsistentWithFreeEnergies) {
    StatMechEngine engine_a(TEMPERATURE);
    StatMechEngine engine_b(TEMPERATURE);

    engine_a.add_sample(-15.0);
    engine_a.add_sample(-12.0);
    engine_b.add_sample(-10.0);
    engine_b.add_sample(-8.0);

    double dG = engine_a.delta_G(engine_b);
    double F_a = engine_a.compute().free_energy;
    double F_b = engine_b.compute().free_energy;

    EXPECT_NEAR(dG, F_a - F_b, EPSILON);
}

// ===========================================================================
// HELMHOLTZ CONVENIENCE FUNCTION
// ===========================================================================

TEST_F(StatMechEngineTest, HelmholtzAgreesWithCompute) {
    std::vector<double> energies = {-15.0, -12.0, -10.0, -8.0};
    for (double e : energies)
        engine.add_sample(e);

    double F_compute = engine.compute().free_energy;
    double F_helmholtz = StatMechEngine::helmholtz(energies, TEMPERATURE);

    EXPECT_NEAR(F_compute, F_helmholtz, EPSILON);
}

TEST_F(StatMechEngineTest, HelmholtzEmptyThrows) {
    std::vector<double> empty;
    EXPECT_THROW(StatMechEngine::helmholtz(empty, TEMPERATURE), std::invalid_argument);
}

TEST_F(StatMechEngineTest, HelmholtzSingleEnergy) {
    std::vector<double> energies = {-10.0};
    double F = StatMechEngine::helmholtz(energies, TEMPERATURE);
    EXPECT_NEAR(F, -10.0, EPSILON);
}

// ===========================================================================
// NUMERICAL STABILITY
// ===========================================================================

TEST_F(StatMechEngineTest, LargeEnergyDifference) {
    // Energy difference >> kT should not cause overflow/NaN
    engine.add_sample(-500.0);
    engine.add_sample(0.0);

    auto th = engine.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.mean_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    EXPECT_TRUE(std::isfinite(th.heat_capacity));

    auto weights = engine.boltzmann_weights();
    for (double w : weights)
        EXPECT_TRUE(std::isfinite(w));
}

TEST_F(StatMechEngineTest, VerySmallEnergyDifferences) {
    // Nearly degenerate states
    for (int i = 0; i < 100; ++i)
        engine.add_sample(-10.0 + i * 1e-10);

    auto th = engine.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    // Nearly degenerate → entropy ≈ k_B ln(100)
    double expected_S = kB_kcal * std::log(100.0);
    EXPECT_NEAR(th.entropy, expected_S, 0.01);
}

// ===========================================================================
// REPLICA EXCHANGE
// ===========================================================================

TEST_F(StatMechEngineTest, InitReplicasCorrectCount) {
    std::vector<double> temps = {200.0, 250.0, 300.0, 350.0, 400.0};
    auto replicas = StatMechEngine::init_replicas(temps);

    ASSERT_EQ(replicas.size(), temps.size());
    for (size_t i = 0; i < temps.size(); ++i) {
        EXPECT_EQ(replicas[i].id, static_cast<int>(i));
        EXPECT_DOUBLE_EQ(replicas[i].temperature, temps[i]);
        EXPECT_NEAR(replicas[i].beta, 1.0 / (kB_kcal * temps[i]), EPSILON);
    }
}

TEST_F(StatMechEngineTest, SwapAcceptedWhenFavorable) {
    // Swap is always accepted when Δ = (β_a - β_b)(E_a - E_b) >= 0
    // β_a > β_b (T_a < T_b) and E_a < E_b → Δ > 0 → swap after: cold gets high E
    // Actually: swap when cold replica has lower energy than hot = favorable
    std::vector<double> temps = {200.0, 400.0};
    auto replicas = StatMechEngine::init_replicas(temps);
    replicas[0].current_energy = -20.0;  // cold replica, low energy
    replicas[1].current_energy = -5.0;   // hot replica, high energy

    // Δ = (β_cold - β_hot)(E_cold - E_hot) = (positive)(negative) = negative
    // This means swap is NOT always accepted. Let's flip to make Δ > 0:
    replicas[0].current_energy = -5.0;   // cold replica, high energy
    replicas[1].current_energy = -20.0;  // hot replica, low energy
    // Δ = (β_cold - β_hot)(E_cold - E_hot) = (positive)(positive) = positive → always accept

    std::mt19937 rng(42);
    bool accepted = StatMechEngine::attempt_swap(replicas[0], replicas[1], rng);
    EXPECT_TRUE(accepted);

    // After swap, energies should be exchanged
    EXPECT_DOUBLE_EQ(replicas[0].current_energy, -20.0);
    EXPECT_DOUBLE_EQ(replicas[1].current_energy, -5.0);
}

TEST_F(StatMechEngineTest, SwapStatisticsPhysical) {
    // Over many trials, acceptance rate should be between 0 and 1
    std::vector<double> temps = {300.0, 350.0};
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> edist(-20.0, 0.0);

    int accepted = 0;
    int trials = 10000;
    for (int i = 0; i < trials; ++i) {
        auto replicas = StatMechEngine::init_replicas(temps);
        replicas[0].current_energy = edist(rng);
        replicas[1].current_energy = edist(rng);
        if (StatMechEngine::attempt_swap(replicas[0], replicas[1], rng))
            accepted++;
    }

    double rate = static_cast<double>(accepted) / trials;
    EXPECT_GT(rate, 0.1);   // not all rejected
    EXPECT_LT(rate, 0.95);  // not all accepted
}

// ===========================================================================
// WHAM (Weighted Histogram Analysis Method)
// ===========================================================================

TEST_F(StatMechEngineTest, BoltzmannPMFBasicOutput) {
    // Simple test: uniform energies, linearly spaced coordinates
    std::vector<double> energies(100);
    std::vector<double> coords(100);
    for (int i = 0; i < 100; ++i) {
        energies[i] = -10.0 + 0.1 * i;
        coords[i] = static_cast<double>(i);
    }

    auto bins = StatMechEngine::boltzmann_pmf(energies, coords, TEMPERATURE, 10);

    ASSERT_EQ(bins.size(), 10u);
    for (const auto& bin : bins) {
        EXPECT_TRUE(std::isfinite(bin.free_energy));
        EXPECT_TRUE(std::isfinite(bin.coord_center));
        EXPECT_GE(bin.count, 0.0);
    }
}

TEST_F(StatMechEngineTest, BoltzmannPMFFreeEnergyMinimumShifted) {
    // All bins should have free_energy >= 0 (shifted so minimum = 0)
    std::vector<double> energies = {-15.0, -12.0, -10.0, -8.0, -6.0};
    std::vector<double> coords = {1.0, 2.0, 3.0, 4.0, 5.0};

    auto bins = StatMechEngine::boltzmann_pmf(energies, coords, TEMPERATURE, 5);
    for (const auto& bin : bins)
        EXPECT_GE(bin.free_energy, -EPSILON);
}

TEST_F(StatMechEngineTest, BoltzmannPMFSizeMismatchThrows) {
    std::vector<double> energies = {-10.0, -8.0};
    std::vector<double> coords = {1.0};

    EXPECT_THROW(
        StatMechEngine::boltzmann_pmf(energies, coords, TEMPERATURE, 5),
        std::invalid_argument
    );
}

TEST_F(StatMechEngineTest, BoltzmannPMFEmptyThrows) {
    std::vector<double> empty;
    EXPECT_THROW(
        StatMechEngine::boltzmann_pmf(empty, empty, TEMPERATURE, 5),
        std::invalid_argument
    );
}

TEST_F(StatMechEngineTest, BoltzmannPMFSingleBin) {
    // Edge case: single bin should produce exactly one result
    std::vector<double> energies = {-10.0, -10.0, -10.0};
    std::vector<double> coords = {0.5, 0.5, 0.5};
    auto result = StatMechEngine::boltzmann_pmf(energies, coords, TEMPERATURE, 1);
    EXPECT_EQ(result.size(), 1u);
    EXPECT_TRUE(std::isfinite(result[0].free_energy));
}

TEST_F(StatMechEngineTest, BoltzmannPMFIdenticalCoordinates) {
    // All samples at same coordinate — all land in one bin
    std::vector<double> energies = {-5.0, -10.0, -15.0};
    std::vector<double> coords = {1.0, 1.0, 1.0};
    auto result = StatMechEngine::boltzmann_pmf(energies, coords, TEMPERATURE, 5);
    // Should not crash; at least one bin populated
    int populated = 0;
    for (const auto& bin : result)
        if (bin.count > 0) ++populated;
    EXPECT_GE(populated, 1);
}

// ===========================================================================
// THERMODYNAMIC INTEGRATION
// ===========================================================================

TEST_F(StatMechEngineTest, TIConstantIntegrand) {
    // ∫₀¹ C dλ = C for constant C
    double C = 5.0;
    std::vector<TIPoint> points = {{0.0, C}, {0.5, C}, {1.0, C}};
    double result = StatMechEngine::thermodynamic_integration(points);
    EXPECT_NEAR(result, C, EPSILON);
}

TEST_F(StatMechEngineTest, TILinearIntegrand) {
    // ∫₀¹ 2λ dλ = 1.0 (trapezoidal is exact for linear)
    int N = 11;
    std::vector<TIPoint> points;
    for (int i = 0; i < N; ++i) {
        double lam = static_cast<double>(i) / (N - 1);
        points.push_back({lam, 2.0 * lam});
    }
    double result = StatMechEngine::thermodynamic_integration(points);
    EXPECT_NEAR(result, 1.0, EPSILON);
}

TEST_F(StatMechEngineTest, TIQuadraticIntegrand) {
    // ∫₀¹ 3λ² dλ = 1.0
    // Trapezoidal rule is approximate for quadratic; use many points
    int N = 1001;
    std::vector<TIPoint> points;
    for (int i = 0; i < N; ++i) {
        double lam = static_cast<double>(i) / (N - 1);
        points.push_back({lam, 3.0 * lam * lam});
    }
    double result = StatMechEngine::thermodynamic_integration(points);
    EXPECT_NEAR(result, 1.0, 1e-4);  // trapezoidal error O(h²)
}

TEST_F(StatMechEngineTest, TITooFewPointsThrows) {
    std::vector<TIPoint> single = {{0.0, 1.0}};
    EXPECT_THROW(StatMechEngine::thermodynamic_integration(single), std::invalid_argument);
}

// ===========================================================================
// BOLTZMANN LOOKUP TABLE
// ===========================================================================

TEST_F(StatMechEngineTest, BoltzmannLUTAccuracy) {
    double beta = 1.0 / (kB_kcal * TEMPERATURE);
    BoltzmannLUT lut(beta, -20.0, 0.0, 10000);

    // Check several energy values within range
    for (double e = -19.0; e <= -1.0; e += 1.0) {
        double exact = std::exp(-beta * e);
        double approx = lut(e);
        double rel_err = std::abs(approx - exact) / exact;
        EXPECT_LT(rel_err, 0.01)  // < 1% relative error
            << "LUT error too large at E=" << e;
    }
}

TEST_F(StatMechEngineTest, BoltzmannLUTBoundary) {
    double beta = 1.0 / (kB_kcal * TEMPERATURE);
    BoltzmannLUT lut(beta, -20.0, 0.0, 1000);

    // Out-of-range values should clamp, not crash
    double below = lut(-100.0);
    double above = lut(100.0);
    EXPECT_TRUE(std::isfinite(below));
    EXPECT_TRUE(std::isfinite(above));
    EXPECT_GT(below, 0.0);
    EXPECT_GT(above, 0.0);
}

// ===========================================================================
// HEAT CAPACITY PROPERTIES
// ===========================================================================

TEST_F(StatMechEngineTest, HeatCapacityNonNegative) {
    std::vector<double> energies = {-20.0, -15.0, -10.0, -5.0, 0.0};
    for (double e : energies)
        engine.add_sample(e);

    auto th = engine.compute();
    EXPECT_GE(th.heat_capacity, 0.0);
}

TEST_F(StatMechEngineTest, HeatCapacityZeroForSingleState) {
    engine.add_sample(-10.0);
    auto th = engine.compute();
    EXPECT_NEAR(th.heat_capacity, 0.0, EPSILON);
}

// Regression for C-1: Boltzmann weights (double in 0..1) were silently truncated
// to int=0 when passed as multiplicity, producing log(0)=-inf → NaN everywhere.
TEST_F(StatMechEngineTest, FractionalMultiplicityNoNaN) {
    StatMechEngine eng(300.0);
    eng.add_sample(-10.0, 0.5);
    eng.add_sample(-8.0,  0.3);
    eng.add_sample(-6.0,  0.2);

    auto th = eng.compute();
    EXPECT_FALSE(std::isnan(th.free_energy));
    EXPECT_FALSE(std::isnan(th.entropy));
    EXPECT_FALSE(std::isnan(th.heat_capacity));
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));

    auto weights = eng.boltzmann_weights();
    for (double w : weights) {
        EXPECT_FALSE(std::isnan(w));
        EXPECT_GE(w, 0.0);
    }
}

// ===========================================================================
// NUMERICAL STABILITY — EXTREME TEMPERATURES
// ===========================================================================

TEST_F(StatMechEngineTest, ExtremelyLowTemperatureFinite) {
    // At T → 0, weight collapses to ground state. No NaN/inf should occur.
    StatMechEngine cold(1.0);  // 1 K
    cold.add_sample(-10.0);
    cold.add_sample(-9.0);
    cold.add_sample(-8.0);

    auto th = cold.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    EXPECT_TRUE(std::isfinite(th.heat_capacity));
    // At 1 K, F ≈ ground state energy
    EXPECT_NEAR(th.free_energy, -10.0, 0.1);
}

TEST_F(StatMechEngineTest, VeryHighTemperatureFinite) {
    // At T → ∞, F → mean energy, S → k_B ln(N)
    StatMechEngine hot(1e8);  // 10^8 K
    int N = 5;
    for (int i = 0; i < N; ++i)
        hot.add_sample(-10.0 - i * 5.0);

    auto th = hot.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    EXPECT_TRUE(std::isfinite(th.heat_capacity));
    // Entropy should approach k_B ln(N)
    double S_max = kB_kcal * std::log(static_cast<double>(N));
    EXPECT_NEAR(th.entropy, S_max, S_max * 0.01); // within 1%
}

TEST_F(StatMechEngineTest, LowTempGroundStateWeightDominates) {
    StatMechEngine cold(1.0);
    cold.add_sample(-100.0);
    for (int i = 0; i < 99; ++i) cold.add_sample(0.0);

    auto w = cold.boltzmann_weights();
    EXPECT_GT(w[0], 0.999); // ground state captures essentially all weight
}

TEST_F(StatMechEngineTest, ExtremeEnergySpreadLogsumexpStable) {
    // If naive exponentiation is used, exp(-β×(-500)) would overflow at T=300.
    // log-sum-exp implementation must avoid this.
    StatMechEngine eng(300.0);
    eng.add_sample(-500.0);
    eng.add_sample(-499.0);
    eng.add_sample(-1.0);
    eng.add_sample(500.0);

    auto th = eng.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    EXPECT_TRUE(std::isfinite(th.heat_capacity));

    // With the extreme spread at T=300, most weight is on the lowest energy.
    // E=-500 vs E=-499: gap is only β×1 ≈ 1.68 kT, so -499 gets ~15.7%.
    auto w = eng.boltzmann_weights();
    EXPECT_GT(w[0], 0.80);  // -500 gets ~84%, -499 gets ~16%
    EXPECT_LT(w[2], 1e-50); // -1 and +500 are effectively zero
    EXPECT_LT(w[3], 1e-50);
}

TEST_F(StatMechEngineTest, AllIdenticalEnergiesNoNan) {
    // N states at same energy: F = E - kT ln N, S = k ln N, Cv = 0
    int N = 1000;
    double E = -7.77;
    StatMechEngine eng(300.0);
    for (int i = 0; i < N; ++i) eng.add_sample(E);

    auto th = eng.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    EXPECT_NEAR(th.heat_capacity, 0.0, 1e-9);
    EXPECT_NEAR(th.entropy, kB_kcal * std::log(static_cast<double>(N)), 1e-9);
}

TEST_F(StatMechEngineTest, LargeEnsembleNumericallyStable) {
    // 10,000 samples spanning a wide energy range
    StatMechEngine eng(300.0);
    std::mt19937 rng(999);
    std::normal_distribution<double> dist(-10.0, 3.0);
    for (int i = 0; i < 10000; ++i) eng.add_sample(dist(rng));

    auto th = eng.compute();
    EXPECT_TRUE(std::isfinite(th.free_energy));
    EXPECT_TRUE(std::isfinite(th.entropy));
    EXPECT_GE(th.entropy, 0.0);
    EXPECT_GE(th.heat_capacity, 0.0);
}

TEST_F(StatMechEngineTest, SingleSampleHighMultiplicity) {
    // Multiplicity M: S = k_B ln(M), F = E - kT ln(M)
    int M = 10000;
    double E = -5.0;
    StatMechEngine eng(300.0);
    eng.add_sample(E, M);

    auto th = eng.compute();
    EXPECT_NEAR(th.entropy, kB_kcal * std::log(static_cast<double>(M)), 1e-9);
    double expected_F = E - kB_kcal * 300.0 * std::log(static_cast<double>(M));
    EXPECT_NEAR(th.free_energy, expected_F, 1e-9);
}

// ===========================================================================
// NUMERICAL STABILITY — PARTITION FUNCTION EDGE CASES
// ===========================================================================

TEST_F(StatMechEngineTest, DeltaGWithSingleStateIsAnalytic) {
    // ΔG(A→B) = F_A − F_B; for single states this is just E_A − E_B
    StatMechEngine eng_a(300.0), eng_b(300.0);
    eng_a.add_sample(-12.0);
    eng_b.add_sample(-8.0);

    double dG = eng_a.delta_G(eng_b);
    EXPECT_NEAR(dG, -12.0 - (-8.0), EPSILON);
}

TEST_F(StatMechEngineTest, HeatCapacityPeakNearTransition) {
    // C_v = Var(E) / kT² peaks at the temperature where the two-state
    // system is half-occupied. At T=300 with ΔE=6.0 kcal/mol:
    //   β·ΔE ≈ 10 → cold side dominates → C_v is near-zero at 300 K.
    // Try ΔE=0.6 kcal/mol (β·ΔE ≈ 1): two-state populations are comparable.
    StatMechEngine eng(300.0);
    eng.add_sample(-10.0);
    eng.add_sample(-9.4);  // ΔE = 0.6 kcal/mol

    auto th = eng.compute();
    EXPECT_GT(th.heat_capacity, 0.0);
}

TEST_F(StatMechEngineTest, EntropyZeroForSingleStateMultiplicity1) {
    StatMechEngine eng(300.0);
    eng.add_sample(-10.0, 1);
    auto th = eng.compute();
    EXPECT_NEAR(th.entropy, 0.0, EPSILON);
}

TEST_F(StatMechEngineTest, FreeEnergyAlwaysLEMeanEnergy) {
    // F = <E> - T*S ≤ <E> because S ≥ 0
    std::vector<double> energies = {-20.0, -15.0, -10.0, -5.0, 0.0, 5.0};
    for (double e : energies) engine.add_sample(e);

    auto th = engine.compute();
    EXPECT_LE(th.free_energy, th.mean_energy + EPSILON);
}

TEST_F(StatMechEngineTest, ComputeTwiceReturnsSameResult) {
    engine.add_sample(-10.0);
    engine.add_sample(-8.0);
    engine.add_sample(-6.0);

    auto th1 = engine.compute();
    auto th2 = engine.compute();
    EXPECT_DOUBLE_EQ(th1.free_energy, th2.free_energy);
    EXPECT_DOUBLE_EQ(th1.entropy, th2.entropy);
    EXPECT_DOUBLE_EQ(th1.heat_capacity, th2.heat_capacity);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
