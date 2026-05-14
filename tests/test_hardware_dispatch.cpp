// tests/test_hardware_dispatch.cpp
// Unit tests for the ShannonThermoStack hardware dispatch layer
// Validates Shannon entropy computation, dispatch backend reporting,
// edge cases, and the ShannonEnergyMatrix singleton.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/ShannonThermoStack/ShannonThermoStack.h"
#include "../LIB/statmech.h"
#include "../LIB/tENCoM/tencm.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <random>
#include <string>
#include <array>

using namespace shannon_thermo;

// ===========================================================================
// CONSTANTS
// ===========================================================================

static constexpr double EPSILON = 1e-6;

static tencm::TorsionalENM make_built_tencm_model(int n_residues = 30) {
    std::vector<std::array<float, 3>> ca;
    ca.reserve(n_residues);
    for (int i = 0; i < n_residues; ++i) {
        ca.push_back({
            2.3f * std::cos(static_cast<float>(i) * 1.74532925f),
            2.3f * std::sin(static_cast<float>(i) * 1.74532925f),
            1.5f * static_cast<float>(i)
        });
    }
    tencm::TorsionalENM model;
    model.build_from_ca(ca);
    return model;
}

// ===========================================================================
// SHANNON ENTROPY — BASIC PROPERTIES
// ===========================================================================

TEST(ShannonEntropy, EmptyInputReturnsZero) {
    std::vector<double> empty;
    EXPECT_DOUBLE_EQ(compute_shannon_entropy(empty), 0.0);
}

TEST(ShannonEntropy, SingleValueReturnsZero) {
    // One sample → one bin with 100% → H = 0
    std::vector<double> single = {5.0};
    EXPECT_DOUBLE_EQ(compute_shannon_entropy(single), 0.0);
}

TEST(ShannonEntropy, IdenticalValuesReturnZero) {
    // All values in one bin → H = 0
    std::vector<double> same(100, 3.14);
    EXPECT_DOUBLE_EQ(compute_shannon_entropy(same), 0.0);
}

TEST(ShannonEntropy, NonNegative) {
    // Shannon entropy is always >= 0
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 5.0);
    std::vector<double> values(500);
    for (auto& v : values) v = dist(rng);

    double H = compute_shannon_entropy(values);
    EXPECT_GE(H, 0.0);
}

TEST(ShannonEntropy, UniformDistributionMaxEntropy) {
    // For a uniform distribution across num_bins bins:
    //   H = ln(num_bins)
    int num_bins = 10;
    int per_bin = 100;
    std::vector<double> values;
    values.reserve(num_bins * per_bin);

    // Create values that land evenly in each bin
    for (int b = 0; b < num_bins; ++b) {
        double center = static_cast<double>(b) + 0.5;
        for (int i = 0; i < per_bin; ++i)
            values.push_back(center);
    }

    double H = compute_shannon_entropy(values, num_bins);
    double H_max = std::log(static_cast<double>(num_bins));
    EXPECT_NEAR(H, H_max, 0.1);  // tolerance for binning edge effects
}

TEST(ShannonEntropy, UpperBound) {
    // H <= ln(num_bins) for any distribution
    int num_bins = 20;
    std::mt19937 rng(123);
    std::normal_distribution<double> dist(0.0, 10.0);
    std::vector<double> values(1000);
    for (auto& v : values) v = dist(rng);

    double H = compute_shannon_entropy(values, num_bins);
    double H_max = std::log(static_cast<double>(num_bins));
    EXPECT_LE(H, H_max + EPSILON);
}

TEST(ShannonEntropy, MoreBinsHigherEntropy) {
    // For the same data, more bins generally gives higher entropy
    std::mt19937 rng(99);
    std::uniform_real_distribution<double> dist(0.0, 100.0);
    std::vector<double> values(5000);
    for (auto& v : values) v = dist(rng);

    double H_5  = compute_shannon_entropy(values, 5);
    double H_20 = compute_shannon_entropy(values, 20);

    EXPECT_GT(H_20, H_5);
}

TEST(ShannonEntropy, DefaultBinsIsValid) {
    // Calling without explicit num_bins should work (DEFAULT_HIST_BINS = 20)
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
    double H = compute_shannon_entropy(values);
    EXPECT_GE(H, 0.0);
    EXPECT_TRUE(std::isfinite(H));
}

TEST(ShannonEntropy, NegativeBinsDefaultsGracefully) {
    // num_bins <= 0 should be corrected to DEFAULT_HIST_BINS
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
    double H = compute_shannon_entropy(values, -1);
    EXPECT_GE(H, 0.0);
    EXPECT_TRUE(std::isfinite(H));
}

// ===========================================================================
// SHANNON ENTROPY — DISCRETE VERSION
// ===========================================================================

TEST(ShannonEntropyDiscrete, EmptyInput) {
    std::vector<int> empty;
    EXPECT_DOUBLE_EQ(compute_shannon_entropy_discrete(empty), 0.0);
}

TEST(ShannonEntropyDiscrete, AllInOneBin) {
    std::vector<int> counts = {100, 0, 0, 0};
    EXPECT_DOUBLE_EQ(compute_shannon_entropy_discrete(counts), 0.0);
}

TEST(ShannonEntropyDiscrete, UniformCounts) {
    // 4 bins, equal counts → H = ln(4) ≈ 1.386 nats
    std::vector<int> counts = {100, 100, 100, 100};
    double H = compute_shannon_entropy_discrete(counts);
    EXPECT_NEAR(H, std::log(4.0), 0.01);
}

TEST(ShannonEntropyDiscrete, TwoBinsEqual) {
    // 2 equal bins → H = ln(2) ≈ 0.693 nats
    std::vector<int> counts = {50, 50};
    double H = compute_shannon_entropy_discrete(counts);
    EXPECT_NEAR(H, std::log(2.0), 0.01);
}

TEST(ShannonEntropyDiscrete, NonNegative) {
    std::vector<int> counts = {10, 20, 30, 5, 1};
    EXPECT_GE(compute_shannon_entropy_discrete(counts), 0.0);
}

// ===========================================================================
// SHANNON ENERGY MATRIX — SINGLETON & INITIALIZATION
// ===========================================================================

TEST(ShannonEnergyMatrix, SingletonReturnsSameInstance) {
    auto& m1 = ShannonEnergyMatrix::instance();
    auto& m2 = ShannonEnergyMatrix::instance();
    EXPECT_EQ(&m1, &m2);
}

TEST(ShannonEnergyMatrix, InitialiseIsIdempotent) {
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise();
    EXPECT_TRUE(mat.is_initialised());

    // Second call should be no-op (early return)
    double v_before = mat.lookup(0, 0);
    mat.initialise();
    double v_after = mat.lookup(0, 0);
    EXPECT_DOUBLE_EQ(v_before, v_after);
}

TEST(ShannonEnergyMatrix, LookupValuesAreFinite) {
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise();

    // Spot-check several entries
    for (int i = 0; i < SHANNON_BINS; i += 32) {
        for (int j = 0; j < SHANNON_BINS; j += 32) {
            double v = mat.lookup(i, j);
            EXPECT_TRUE(std::isfinite(v))
                << "Non-finite value at (" << i << "," << j << ")";
        }
    }
}

TEST(ShannonEnergyMatrix, LookupDeterministic) {
    // Matrix is seeded with 42 → same values every time
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise();

    double v1 = mat.lookup(10, 20);
    double v2 = mat.lookup(10, 20);
    EXPECT_DOUBLE_EQ(v1, v2);
}

TEST(ShannonEnergyMatrix, DiagonalEntriesFinite) {
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise();

    for (int i = 0; i < SHANNON_BINS; i += 16) {
        EXPECT_TRUE(std::isfinite(mat.lookup(i, i)));
    }
}

TEST(ShannonEnergyMatrix, InitialiseFromDataUndersized) {
    // Bug fix test: initialise_from_data with fewer elements than 256*256
    // Previously caused OOB read by clamping count UP to expected
    std::vector<float> small(100, 1.0f);
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise_from_data(small.data(), 100);
    EXPECT_TRUE(mat.is_initialised());
    // First 100 values should be 1.0, rest should be zero-filled
    EXPECT_NEAR(mat.data()[0], 1.0, 1e-6);
    EXPECT_NEAR(mat.data()[99], 1.0, 1e-6);
    EXPECT_NEAR(mat.data()[100], 0.0, 1e-6);
    EXPECT_NEAR(mat.data()[SHANNON_BINS * SHANNON_BINS - 1], 0.0, 1e-6);
}

TEST(ShannonEnergyMatrix, InitialiseFromDataExactSize) {
    // Verify full-size buffer works correctly
    std::vector<float> full(SHANNON_BINS * SHANNON_BINS, 0.5f);
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise_from_data(full.data(), SHANNON_BINS * SHANNON_BINS);
    EXPECT_TRUE(mat.is_initialised());
    EXPECT_NEAR(mat.data()[0], 0.5, 1e-6);
    EXPECT_NEAR(mat.data()[SHANNON_BINS * SHANNON_BINS - 1], 0.5, 1e-6);
}

TEST(ShannonEnergyMatrix, InitialiseFromDataOversized) {
    // count > expected should only read expected elements
    std::vector<float> big(SHANNON_BINS * SHANNON_BINS + 1000, 2.0f);
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise_from_data(big.data(), SHANNON_BINS * SHANNON_BINS + 1000);
    EXPECT_TRUE(mat.is_initialised());
    EXPECT_EQ(mat.size(), SHANNON_BINS * SHANNON_BINS);
}

// ===========================================================================
// TORSIONAL VIBRATIONAL ENTROPY
// ===========================================================================

TEST(TorsionalVibEntropy, EmptyModesReturnsZero) {
    std::vector<tencm::NormalMode> empty;
    EXPECT_DOUBLE_EQ(compute_torsional_vibrational_entropy(empty), 0.0);
}

TEST(TorsionalVibEntropy, SkipsFirstSixModes) {
    // First 6 modes (translation+rotation) should be skipped
    // If we provide exactly 6 modes, result is 0
    std::vector<tencm::NormalMode> modes(6);
    for (int i = 0; i < 6; ++i)
        modes[i].eigenvalue = 1.0;  // non-trivial eigenvalues

    EXPECT_DOUBLE_EQ(compute_torsional_vibrational_entropy(modes), 0.0);
}

TEST(TorsionalVibEntropy, SkipsNearZeroEigenvalues) {
    // Modes with eigenvalue < 1e-6 are skipped
    std::vector<tencm::NormalMode> modes(10);
    for (int i = 0; i < 10; ++i)
        modes[i].eigenvalue = 1e-9;  // effectively zero

    EXPECT_DOUBLE_EQ(compute_torsional_vibrational_entropy(modes), 0.0);
}

TEST(TorsionalVibEntropy, ValidModesProducePositiveEntropy) {
    // Modes 6+ with reasonable eigenvalues should give S > 0
    std::vector<tencm::NormalMode> modes(12);
    for (int i = 0; i < 12; ++i)
        modes[i].eigenvalue = 0.5 + 0.1 * i;

    double S = compute_torsional_vibrational_entropy(modes, 300.0);
    EXPECT_GT(S, 0.0);
    EXPECT_TRUE(std::isfinite(S));
}

TEST(TorsionalVibEntropy, HigherTemperatureHigherEntropy) {
    std::vector<tencm::NormalMode> modes(10);
    for (int i = 0; i < 10; ++i)
        modes[i].eigenvalue = 1.0;

    double S_low  = compute_torsional_vibrational_entropy(modes, 200.0);
    double S_high = compute_torsional_vibrational_entropy(modes, 500.0);

    EXPECT_GT(S_high, S_low);
}

TEST(TorsionalVibEntropy, ResultIsFiniteForLargeModes) {
    std::vector<tencm::NormalMode> modes(100);
    for (int i = 0; i < 100; ++i)
        modes[i].eigenvalue = 0.01 * (i + 1);

    double S = compute_torsional_vibrational_entropy(modes, 298.15);
    EXPECT_TRUE(std::isfinite(S));
}

// ===========================================================================
// RUN_SHANNON_THERMO_STACK — FULL PIPELINE
// ===========================================================================

class ShannonThermoStackTest : public ::testing::Test {
protected:
    statmech::StatMechEngine engine{298.15};
    tencm::TorsionalENM tencm_model;  // default-constructed, not built

    void SetUp() override {
        // Populate with a realistic ensemble of energies
        engine.add_sample(-15.0);
        engine.add_sample(-12.0);
        engine.add_sample(-10.0);
        engine.add_sample(-8.0);
        engine.add_sample(-6.0);
    }
};

TEST_F(ShannonThermoStackTest, ProducesFiniteResults) {
    auto result = run_shannon_thermo_stack(engine, tencm_model, -10.0);

    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_TRUE(std::isfinite(result.shannonEntropy));
    EXPECT_TRUE(std::isfinite(result.torsionalVibEntropy));
    EXPECT_TRUE(std::isfinite(result.entropyContribution));
}

TEST_F(ShannonThermoStackTest, ShannonEntropyNonNegative) {
    auto result = run_shannon_thermo_stack(engine, tencm_model, -10.0);
    EXPECT_GE(result.shannonEntropy, 0.0);
}

TEST_F(ShannonThermoStackTest, TorsionalEntropyZeroWhenNotBuilt) {
    // Default-constructed TorsionalENM is not built → S_vib = 0
    EXPECT_FALSE(tencm_model.is_built());
    auto result = run_shannon_thermo_stack(engine, tencm_model, -10.0);
    EXPECT_DOUBLE_EQ(result.torsionalVibEntropy, 0.0);
}

TEST_F(ShannonThermoStackTest, EntropyContributionIsNegative) {
    // -T*S <= 0 (entropy always contributes favorably to free energy)
    auto result = run_shannon_thermo_stack(engine, tencm_model, -10.0);
    EXPECT_LE(result.entropyContribution, 0.0 + EPSILON);
}

TEST_F(ShannonThermoStackTest, DeltaGIncorporatesEntropy) {
    double base_dG = -10.0;
    auto result = run_shannon_thermo_stack(engine, tencm_model, base_dG);

    // deltaG = base_dG + entropy_contribution
    EXPECT_NEAR(result.deltaG, base_dG + result.entropyContribution, EPSILON);
}

TEST_F(ShannonThermoStackTest, HeuristicTorsionalEntropyIsExcludedFromDeltaG) {
    double base_dG = -10.0;
    tencm::TorsionalENM built_model = make_built_tencm_model();
    ASSERT_TRUE(built_model.is_built());

    auto no_tencm = run_shannon_thermo_stack(engine, tencm_model, base_dG);
    auto with_tencm = run_shannon_thermo_stack(engine, built_model, base_dG);

    EXPECT_GT(with_tencm.torsionalVibEntropy, 0.0);
    EXPECT_NEAR(with_tencm.entropyContribution,
                no_tencm.entropyContribution,
                EPSILON);
    EXPECT_NEAR(with_tencm.deltaG, no_tencm.deltaG, EPSILON);
    EXPECT_NE(with_tencm.report.find("excluded from dG"), std::string::npos)
        << with_tencm.report;
}

TEST_F(ShannonThermoStackTest, ReportContainsBackendName) {
    auto result = run_shannon_thermo_stack(engine, tencm_model, -10.0);

    // Report must contain "ShannonThermoStack["
    EXPECT_NE(result.report.find("ShannonThermoStack["), std::string::npos);

    // Must contain one of the known backend names
    bool has_backend =
        result.report.find("CUDA") != std::string::npos ||
        result.report.find("Metal") != std::string::npos ||
        result.report.find("AVX-512") != std::string::npos ||
        result.report.find("OpenMP") != std::string::npos ||
        result.report.find("scalar") != std::string::npos;
    EXPECT_TRUE(has_backend) << "Report missing backend: " << result.report;
}

TEST_F(ShannonThermoStackTest, ReportContainsMetrics) {
    auto result = run_shannon_thermo_stack(engine, tencm_model, -10.0);

    EXPECT_NE(result.report.find("S_conf="), std::string::npos);
    EXPECT_NE(result.report.find("S_vib_heuristic="), std::string::npos);
    EXPECT_NE(result.report.find("kcal/mol"), std::string::npos);
}

// ===========================================================================
// RUN_SHANNON_THERMO_STACK — EDGE CASES
// ===========================================================================

TEST(ShannonThermoStackEdge, SingleSampleEnsemble) {
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-10.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);
    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_TRUE(std::isfinite(result.shannonEntropy));
}

TEST(ShannonThermoStackEdge, LargeEnsemble) {
    statmech::StatMechEngine eng(298.15);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(-15.0, 5.0);

    for (int i = 0; i < 10000; ++i)
        eng.add_sample(dist(rng));

    tencm::TorsionalENM tencm;
    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);

    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_GT(result.shannonEntropy, 0.0);
}

TEST(ShannonThermoStackEdge, ZeroBaseDeltaG) {
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, 0.0);
    // deltaG should equal just the entropy contribution
    EXPECT_NEAR(result.deltaG, result.entropyContribution, EPSILON);
}

TEST(ShannonThermoStackEdge, CustomTemperature) {
    statmech::StatMechEngine eng(310.0);
    eng.add_sample(-12.0);
    eng.add_sample(-8.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0, 310.0);
    EXPECT_TRUE(std::isfinite(result.deltaG));
}

TEST(ShannonThermoStackEdge, DegenerateEnergyEnsemble) {
    // All samples at same energy → Boltzmann weights are equal (w_i = 1/N)
    // → S = -Σ w_i·ln(w_i) = ln(N) = ln(50) ≈ 3.912 nats
    statmech::StatMechEngine eng(298.15);
    for (int i = 0; i < 50; ++i)
        eng.add_sample(-10.0);

    tencm::TorsionalENM tencm;
    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);

    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_NEAR(result.shannonEntropy, std::log(50.0), 0.01);
}

// ===========================================================================
// NUMERICAL STABILITY
// ===========================================================================

TEST(ShannonEntropyStability, VeryLargeValues) {
    std::vector<double> values(100);
    for (int i = 0; i < 100; ++i)
        values[i] = 1e6 + i * 0.01;

    double H = compute_shannon_entropy(values, 10);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
}

TEST(ShannonEntropyStability, VerySmallRange) {
    std::vector<double> values(100);
    for (int i = 0; i < 100; ++i)
        values[i] = 1.0 + i * 1e-15;

    double H = compute_shannon_entropy(values, 10);
    EXPECT_TRUE(std::isfinite(H));
}

TEST(ShannonEntropyStability, NegativeValues) {
    std::vector<double> values = {-100.0, -50.0, -25.0, -10.0, -5.0};
    double H = compute_shannon_entropy(values, 5);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
}

TEST(ShannonEntropyStability, MixedSignValues) {
    std::vector<double> values = {-50.0, -25.0, 0.0, 25.0, 50.0};
    double H = compute_shannon_entropy(values, 5);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
}

TEST(ShannonEntropyStability, TwoValues) {
    // Minimum non-trivial case
    std::vector<double> values = {0.0, 1.0};
    double H = compute_shannon_entropy(values, 2);
    EXPECT_TRUE(std::isfinite(H));
}

// ===========================================================================
// DISPATCH BACKEND CONSISTENCY
// ===========================================================================

TEST(DispatchConsistency, ReproducibleResults) {
    // Same input → same output (deterministic regardless of backend)
    std::mt19937 rng(77);
    std::normal_distribution<double> dist(0.0, 10.0);
    std::vector<double> values(1000);
    for (auto& v : values) v = dist(rng);

    double H1 = compute_shannon_entropy(values, 20);
    double H2 = compute_shannon_entropy(values, 20);
    EXPECT_DOUBLE_EQ(H1, H2);
}

TEST(DispatchConsistency, EntropyMonotonicWithDataSpread) {
    // Wider spread → more bins occupied → higher entropy
    std::vector<double> narrow(500), wide(500);
    std::mt19937 rng(42);

    std::normal_distribution<double> n_dist(0.0, 1.0);
    std::normal_distribution<double> w_dist(0.0, 50.0);
    for (int i = 0; i < 500; ++i) {
        narrow[i] = n_dist(rng);
        wide[i]   = w_dist(rng);
    }

    double H_narrow = compute_shannon_entropy(narrow, 20);
    double H_wide   = compute_shannon_entropy(wide, 20);

    EXPECT_GT(H_wide, H_narrow);
}

// ===========================================================================
// COMPILE-TIME BACKEND DETECTION
// ===========================================================================

TEST(BackendDetection, ActiveBackendReported) {
    // Verify that the runtime report correctly reflects compilation flags
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);

#if defined(FLEXAIDS_USE_CUDA)
    EXPECT_NE(result.report.find("CUDA"), std::string::npos)
        << "Expected CUDA in report: " << result.report;
#elif defined(ENABLE_METAL_CORE)
    EXPECT_NE(result.report.find("Metal"), std::string::npos)
        << "Expected Metal in report: " << result.report;
#elif defined(__AVX512F__)
    EXPECT_NE(result.report.find("AVX-512"), std::string::npos)
        << "Expected AVX-512 in report: " << result.report;
#elif defined(_OPENMP)
    EXPECT_NE(result.report.find("OpenMP"), std::string::npos)
        << "Expected OpenMP in report: " << result.report;
#else
    EXPECT_NE(result.report.find("scalar"), std::string::npos)
        << "Expected scalar in report: " << result.report;
#endif
}

TEST(BackendDetection, EigenTagInReport) {
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-10.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);
    EXPECT_NE(result.report.find("Eigen"), std::string::npos)
        << "Expected +Eigen in report: " << result.report;
}

// ===========================================================================
// DISPATCH EDGE CASES — UNIFIED DISPATCH LAYER COVERAGE
// ===========================================================================

// --- Shannon entropy: single bin always returns zero ---

TEST(ShannonEntropyEdge, SingleBinReturnsZero) {
    // num_bins=1 → all values in one bin → p=1.0 → H = -1*ln(1) = 0
    std::mt19937 rng(55);
    std::normal_distribution<double> dist(0.0, 10.0);
    std::vector<double> values(200);
    for (auto& v : values) v = dist(rng);

    double H = compute_shannon_entropy(values, 1);
    EXPECT_DOUBLE_EQ(H, 0.0);
}

// --- Shannon entropy: more bins than data points ---

TEST(ShannonEntropyEdge, MoreBinsThanSamples) {
    // 5 distinct values spread into 1000 bins → most bins empty, H still valid
    std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
    double H = compute_shannon_entropy(values, 1000);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
    // With 5 values in 1000 bins, at most 5 bins occupied → H <= ln(5)
    EXPECT_LE(H, std::log(5.0) + 0.01);
}

// --- Shannon entropy: very large number of bins ---

TEST(ShannonEntropyEdge, LargeBinCount) {
    // Stress test all dispatch paths with a large bin count
    std::mt19937 rng(88);
    std::uniform_real_distribution<double> dist(0.0, 100.0);
    std::vector<double> values(2000);
    for (auto& v : values) v = dist(rng);

    double H = compute_shannon_entropy(values, 500);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
    EXPECT_LE(H, std::log(500.0) + EPSILON);
}

// --- Shannon entropy discrete: all-zero counts ---

TEST(ShannonEntropyDiscreteEdge, AllZeroCounts) {
    std::vector<int> counts = {0, 0, 0, 0};
    EXPECT_DOUBLE_EQ(compute_shannon_entropy_discrete(counts), 0.0);
}

// --- Shannon entropy discrete: single non-zero bin ---

TEST(ShannonEntropyDiscreteEdge, SingleNonZeroBin) {
    std::vector<int> counts = {0, 0, 500, 0, 0};
    EXPECT_DOUBLE_EQ(compute_shannon_entropy_discrete(counts), 0.0);
}

// --- Shannon entropy discrete: very large counts ---

TEST(ShannonEntropyDiscreteEdge, VeryLargeCounts) {
    // Ensure no integer overflow in total computation
    std::vector<int> counts = {1000000, 1000000, 1000000};
    double H = compute_shannon_entropy_discrete(counts);
    EXPECT_NEAR(H, std::log(3.0), 0.01);
}

// --- ShannonEnergyMatrix: boundary indices ---

TEST(ShannonEnergyMatrixEdge, BoundaryIndices) {
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise();

    // First and last valid indices
    EXPECT_TRUE(std::isfinite(mat.lookup(0, 0)));
    EXPECT_TRUE(std::isfinite(mat.lookup(SHANNON_BINS - 1, SHANNON_BINS - 1)));
    EXPECT_TRUE(std::isfinite(mat.lookup(0, SHANNON_BINS - 1)));
    EXPECT_TRUE(std::isfinite(mat.lookup(SHANNON_BINS - 1, 0)));
}

// --- ShannonEnergyMatrix: symmetry check ---

TEST(ShannonEnergyMatrixEdge, AsymmetricByDesign) {
    // E[i][j] = -kT * p_i * ln(p_j) ≠ E[j][i] in general
    // because p_i and p_j come from independent distributions
    auto& mat = ShannonEnergyMatrix::instance();
    mat.initialise();

    // Just verify both are finite; they need not be equal
    double v_ij = mat.lookup(10, 50);
    double v_ji = mat.lookup(50, 10);
    EXPECT_TRUE(std::isfinite(v_ij));
    EXPECT_TRUE(std::isfinite(v_ji));
}

// --- Full stack: very low temperature ---

TEST(ShannonThermoStackEdge, VeryLowTemperature) {
    // Near 0 K: entropy contribution should be small (−T·S → 0)
    statmech::StatMechEngine eng(1.0);  // 1 K
    eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    eng.add_sample(-8.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0, 1.0);
    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_TRUE(std::isfinite(result.shannonEntropy));
    // At very low T, entropy contribution magnitude should be small
    EXPECT_LT(std::abs(result.entropyContribution), 1.0);
}

// --- Full stack: high temperature ---

TEST(ShannonThermoStackEdge, HighTemperature) {
    statmech::StatMechEngine eng(1000.0);
    eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    eng.add_sample(-8.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0, 1000.0);
    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_TRUE(std::isfinite(result.shannonEntropy));
    // At high T, entropy contribution should be larger in magnitude
    // compared to standard temperature (298.15 K)
    auto result_std = run_shannon_thermo_stack(eng, tencm, -10.0, 298.15);
    EXPECT_LE(result.entropyContribution, result_std.entropyContribution);
}

// --- Full stack: two-sample ensemble ---

TEST(ShannonThermoStackEdge, TwoSampleEnsemble) {
    // Minimal non-trivial ensemble
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -7.5);
    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_NEAR(result.deltaG, -7.5 + result.entropyContribution, EPSILON);
}

// --- Full stack: highly degenerate large ensemble ---

TEST(ShannonThermoStackEdge, HighlyDegenerateEnsemble) {
    // 1000 samples at -10.0, 1 outlier at -5.0
    // Direct Shannon entropy: w_high ≈ 1000/1001, w_low ≈ 1/1001
    // S = -w_high·ln(w_high) - w_low·ln(w_low) ≈ ln(1001) ≈ 6.91 nats
    statmech::StatMechEngine eng(298.15);
    for (int i = 0; i < 1000; ++i)
        eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);
    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_TRUE(std::isfinite(result.shannonEntropy));
    // Correct direct entropy for 1001 samples with nearly uniform weights
    EXPECT_NEAR(result.shannonEntropy, std::log(1001.0), 0.1);
}

// --- Full stack: wide energy spread ---

TEST(ShannonThermoStackEdge, WideEnergySpread) {
    // Energies spanning a huge range: numerical stability test
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-100.0);
    eng.add_sample(-50.0);
    eng.add_sample(0.0);
    eng.add_sample(50.0);
    eng.add_sample(100.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, -10.0);
    EXPECT_TRUE(std::isfinite(result.deltaG));
    EXPECT_TRUE(std::isfinite(result.shannonEntropy));
    EXPECT_GE(result.shannonEntropy, 0.0);
}

// --- Torsional vibrational entropy: single valid mode ---

TEST(TorsionalVibEntropyEdge, SingleValidMode) {
    // 7 modes: first 6 skipped, only mode[6] counted
    std::vector<tencm::NormalMode> modes(7);
    for (int i = 0; i < 7; ++i)
        modes[i].eigenvalue = 1.0;

    double S = compute_torsional_vibrational_entropy(modes, 298.15);
    EXPECT_GT(S, 0.0);
    EXPECT_TRUE(std::isfinite(S));
}

// --- Torsional vibrational entropy: mixed valid/invalid modes ---

TEST(TorsionalVibEntropyEdge, MixedValidInvalidModes) {
    // Some modes above threshold, some below
    std::vector<tencm::NormalMode> modes(12);
    for (int i = 0; i < 12; ++i)
        modes[i].eigenvalue = (i % 2 == 0) ? 1e-9 : 1.0;

    double S = compute_torsional_vibrational_entropy(modes, 298.15);
    EXPECT_TRUE(std::isfinite(S));
    // Only 3 valid modes (indices 7, 9, 11) should contribute
    EXPECT_GT(S, 0.0);
}

// ===========================================================================
// ENTROPY PLATEAU DETECTION (Improvement 1)
// ===========================================================================

TEST(EntropyPlateau, ConstantHistoryDetected) {
    // All entries identical → plateau
    std::vector<double> history = {2.5, 2.5, 2.5, 2.5, 2.5};
    EXPECT_TRUE(detect_entropy_plateau(history, 5, 0.01));
}

TEST(EntropyPlateau, DecreasingHistoryNotDetected) {
    // Monotonically decreasing → no plateau
    std::vector<double> history = {5.0, 4.0, 3.0, 2.0, 1.0};
    EXPECT_FALSE(detect_entropy_plateau(history, 5, 0.01));
}

TEST(EntropyPlateau, PlateauAfterDecrease) {
    // Decrease then stabilise
    std::vector<double> history = {5.0, 3.0, 1.5, 1.5, 1.5, 1.5};
    EXPECT_TRUE(detect_entropy_plateau(history, 3, 0.01));
    EXPECT_FALSE(detect_entropy_plateau(history, 5, 0.01));
}

TEST(EntropyPlateau, WindowLargerThanHistory) {
    std::vector<double> history = {2.5, 2.5};
    EXPECT_FALSE(detect_entropy_plateau(history, 5, 0.01));
}

TEST(EntropyPlateau, EmptyHistory) {
    std::vector<double> history;
    EXPECT_FALSE(detect_entropy_plateau(history, 5, 0.01));
}

TEST(EntropyPlateau, ZeroEntropyPlateau) {
    // All zeros → plateau (degenerate ensemble)
    std::vector<double> history = {0.0, 0.0, 0.0};
    EXPECT_TRUE(detect_entropy_plateau(history, 3, 0.01));
}

TEST(EntropyPlateau, NearThreshold) {
    // Values within 1% of each other
    std::vector<double> history = {1.000, 1.005, 1.009, 1.003};
    EXPECT_TRUE(detect_entropy_plateau(history, 4, 0.01));
    // But not within 0.1%
    EXPECT_FALSE(detect_entropy_plateau(history, 4, 0.001));
}

// ===========================================================================
// CORRECTED SCALING FORMULA (Improvement 2)
// ===========================================================================

TEST(ScalingFormula, AdditiveDecomposition) {
    // Verify total_S = S_conf + S_vib (no quadratic term)
    statmech::StatMechEngine eng(298.15);
    // Add samples with varied energies to get non-trivial Shannon entropy
    for (int i = 0; i < 100; ++i)
        eng.add_sample(-10.0 + 0.1 * i);

    // Build tencm with some valid modes
    tencm::TorsionalENM tencm;
    // No build → S_vib = 0, so entropyContribution = -T * S_conf_phys

    auto result = run_shannon_thermo_stack(eng, tencm, 0.0);
    // S_conf_phys = H_nats * kB (direct, since H is in nats)
    double expected_S_conf = result.shannonEntropy * kB_kcal;
    // entropyContribution = -T * (S_conf + 0) since no tencm modes
    double expected_contrib = -298.15 * expected_S_conf;
    EXPECT_NEAR(result.entropyContribution, expected_contrib, 1e-8);
}

TEST(ScalingFormula, NatsConversion) {
    // Verify S_conf = k_B * H_nats directly (no log2/ln(2) conversion needed)
    statmech::StatMechEngine eng(298.15);
    eng.add_sample(-10.0);
    eng.add_sample(-5.0);
    tencm::TorsionalENM tencm;

    auto result = run_shannon_thermo_stack(eng, tencm, 0.0);
    double H = result.shannonEntropy;
    double expected_S = H * kB_kcal;
    double expected_contrib = -298.15 * expected_S;
    EXPECT_NEAR(result.entropyContribution, expected_contrib, 1e-8);
}

// ===========================================================================
// GPU DISPATCH THRESHOLD (Improvement 3)
// ===========================================================================

TEST(GPUThreshold, SmallInputUsesScalarPath) {
    // With N=100, results should be identical to scalar regardless of
    // compile-time GPU flags (GPU dispatch is gated by threshold)
    std::mt19937 rng(123);
    std::normal_distribution<double> dist(0.0, 5.0);
    std::vector<double> values(100);
    for (auto& v : values) v = dist(rng);

    double H = compute_shannon_entropy(values, 20);
    EXPECT_GE(H, 0.0);
    EXPECT_TRUE(std::isfinite(H));

    // Compute again — should be deterministic
    double H2 = compute_shannon_entropy(values, 20);
    EXPECT_DOUBLE_EQ(H, H2);
}

TEST(GPUThreshold, ResultConsistencyAcrossSizes) {
    // Same distribution at different sizes should give similar entropy
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);

    std::vector<double> small(500), large(5000);
    for (auto& v : small) v = dist(rng);
    rng.seed(42);
    for (auto& v : large) v = dist(rng);

    double H_small = compute_shannon_entropy(small, 20);
    double H_large = compute_shannon_entropy(large, 20);

    // Both should be positive and in similar range for same distribution
    EXPECT_GT(H_small, 0.0);
    EXPECT_GT(H_large, 0.0);
    // Entropy from same distribution should be reasonably close
    EXPECT_NEAR(H_small, H_large, 1.0);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
