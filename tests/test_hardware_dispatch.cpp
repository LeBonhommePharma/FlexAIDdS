// tests/test_hardware_dispatch.cpp
// Unit tests for the ShannonThermoStack hardware dispatch layer
// Validates Shannon entropy computation, dispatch backend reporting,
// edge cases, and the ShannonEnergyMatrix singleton.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/ShannonThermoStack/ShannonThermoStack.h"
#include "../LIB/UnifiedHardwareDispatch.h"
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

TEST(TorsionalVibEntropy, InternalModesAreNotSkipped) {
    // Internal (torsional) coordinates carry no rigid-body translation/rotation,
    // so there is no 6-fold zero-mode manifold to skip. With 6 modes above the
    // eigenvalue threshold, entropy is POSITIVE — previously this incorrectly
    // returned 0 because the first 6 modes (the softest real torsions) were
    // positionally skipped.
    std::vector<tencm::NormalMode> modes(6);
    for (int i = 0; i < 6; ++i)
        modes[i].eigenvalue = 1.0;  // non-trivial eigenvalues

    EXPECT_GT(compute_torsional_vibrational_entropy(modes, 298.15), 0.0);
}

TEST(TorsionalVibEntropy, SkipsNearZeroEigenvalues) {
    // Modes with eigenvalue < 1e-6 are skipped
    std::vector<tencm::NormalMode> modes(10);
    for (int i = 0; i < 10; ++i)
        modes[i].eigenvalue = 1e-9;  // effectively zero

    EXPECT_DOUBLE_EQ(compute_torsional_vibrational_entropy(modes), 0.0);
}

TEST(TorsionalVibEntropy, ValidModesProducePositiveEntropy) {
    // All modes above the eigenvalue threshold contribute (no positional skip)
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
    // Internal-coordinate torsional modes are NOT positionally skipped: every
    // mode above the eigenvalue threshold contributes (no rigid-body zero-mode
    // manifold in dihedral space).
    std::vector<tencm::NormalMode> modes(7);
    for (int i = 0; i < 7; ++i)
        modes[i].eigenvalue = 1.0;

    double S = compute_torsional_vibrational_entropy(modes, 298.15);
    EXPECT_GT(S, 0.0);
    EXPECT_TRUE(std::isfinite(S));

    // Regression pin for the internal-coordinate fix: with exactly 6 equal modes
    // the old "skip first 6" logic left the eigenvalue buffer empty and returned
    // 0. All 6 must now count, and (all equal) S must be exactly 6x a single mode.
    std::vector<tencm::NormalMode> six(6);
    for (auto& m : six) m.eigenvalue = 1.0;
    std::vector<tencm::NormalMode> one(1, tencm::NormalMode{1.0, {}});
    double S6 = compute_torsional_vibrational_entropy(six, 298.15);
    double S1 = compute_torsional_vibrational_entropy(one, 298.15);
    EXPECT_GT(S6, 0.0);
    EXPECT_NEAR(S6, 6.0 * S1, 1e-9);
}

// --- Torsional vibrational entropy: mixed valid/invalid modes ---

TEST(TorsionalVibEntropyEdge, MixedValidInvalidModes) {
    // Some modes above threshold, some below
    std::vector<tencm::NormalMode> modes(12);
    for (int i = 0; i < 12; ++i)
        modes[i].eigenvalue = (i % 2 == 0) ? 1e-9 : 1.0;

    double S = compute_torsional_vibrational_entropy(modes, 298.15);
    EXPECT_TRUE(std::isfinite(S));
    // No positional skip: all 6 odd-index modes (1,3,5,7,9,11) contribute; the 6
    // sub-threshold even-index modes are dropped by the eigenvalue guard. (Under
    // the old skip-first-6 logic only indices 7,9,11 counted.)
    EXPECT_GT(S, 0.0);
    std::vector<tencm::NormalMode> one(1, tencm::NormalMode{1.0, {}});
    double S1 = compute_torsional_vibrational_entropy(one, 298.15);
    EXPECT_NEAR(S, 6.0 * S1, 1e-9);
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
// OUTLIER-ROBUST HISTOGRAM SUPPORT
//
// The histogram support is normally the sample's own [min, max], which lets a
// single extreme value rescale every bin. These tests pin both halves of the
// contract: the robust fence must neutralise a lone pathological sample, and
// it must stay out of the way for well-behaved data.
// ===========================================================================

TEST(RobustBinning, SingleClashPoseDoesNotFakeCollapse) {
    // A diverse population plus one clash/wall pose. Without a robust support
    // every real sample lands in bin 0 and H reads ~0.056 nats (0.08 bits),
    // far below the hard-collapse line — a GA gate reading this would stop a
    // population that has not collapsed at all.
    std::vector<double> energies;
    energies.reserve(100);
    for (int i = 0; i < 99; ++i)
        energies.push_back(-100.0 + 0.5 * static_cast<double>(i));  // spread over [-100,-51]
    energies.push_back(1.0e4);                                      // clash pose

    const double H = compute_shannon_entropy(energies, DEFAULT_HIST_BINS);

    EXPECT_GT(H, kHSC_soft_nats)
        << "One outlier collapsed the histogram support: H=" << H
        << " nats <= soft gate " << kHSC_soft_nats;
    EXPECT_GT(H, kHSC_hard_nats);
}

TEST(RobustBinning, UniformDataIsUnaffectedByTheFence) {
    // Range = 2*IQR for a uniform sample, far below the trigger, so the fence
    // must not engage and H must still saturate at ln(num_bins).
    std::vector<double> values(2000);
    for (std::size_t i = 0; i < values.size(); ++i)
        values[i] = static_cast<double>(i) / static_cast<double>(values.size());

    const int bins = 20;
    const double H = compute_shannon_entropy(values, bins);
    EXPECT_NEAR(H, std::log(static_cast<double>(bins)), 1e-9);
}

TEST(RobustBinning, GaussianDataIsUnaffectedByTheFence) {
    // A Gaussian's range is only a few IQR, so results must be bit-identical
    // to the pre-fence estimator. Computed here against the same data run
    // through an explicit min/max histogram.
    std::mt19937 rng(12345);
    std::normal_distribution<double> dist(0.0, 1.0);
    std::vector<double> values(1000);
    for (auto& v : values) v = dist(rng);

    const int bins = 20;
    const double H = compute_shannon_entropy(values, bins);

    // Reference: plain min/max binning, no fence.
    const auto [it_min, it_max] = std::minmax_element(values.begin(), values.end());
    const double bw = (*it_max - *it_min) / bins + 1e-10;
    std::vector<int> counts(bins, 0);
    for (double v : values) {
        int b = static_cast<int>((v - *it_min) / bw);
        counts[std::min(std::max(b, 0), bins - 1)]++;
    }
    double H_ref = 0.0;
    for (int c : counts) {
        if (c > 0) {
            const double p = static_cast<double>(c) / static_cast<double>(values.size());
            H_ref -= p * std::log(p);
        }
    }
    EXPECT_NEAR(H, H_ref, 1e-12);
}

TEST(RobustBinning, FarOutlierLandsInTheTopBinNotOnTopOfTheBulk) {
    // Regression for a narrowing-conversion defect the fence introduced.
    //
    // With a fenced support, an outlier's raw bin index is unbounded: a tight
    // bulk plus a clash pose at 1e4 produces an index far above INT_MAX, and
    // converting that to int is undefined behaviour. x86-64 cvttsd2si yields
    // INT_MIN, which clamps to bin 0 — putting the clash pose ON TOP OF the
    // bulk, collapsing the histogram to one occupied bin (H = 0) and firing the
    // very gate the fence exists to protect — while ARM saturates to INT_MAX
    // and lands it in the top bin. n >= kRobustMinSamples so the fence engages.
    std::vector<double> energies;
    energies.reserve(100);
    for (int i = 0; i < 99; ++i)
        energies.push_back(-100.0 + 1.0e-6 * static_cast<double>(i));
    energies.push_back(1.0e4);

    const double H = compute_shannon_entropy(energies, DEFAULT_HIST_BINS);

    // The fence resolves the bulk across several bins and isolates the clash
    // pose, so H must sit well clear of the collapse gates and below the ln N
    // ceiling. If the outlier's index were converted while out of int32 range
    // and clamped to bin 0, it would merge into the bulk and H would fall.
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GT(H, kHSC_soft_nats) << "H=" << H;
    EXPECT_LT(H, std::log(static_cast<double>(DEFAULT_HIST_BINS)) + 1e-12);
}

TEST(RobustBinning, DegenerateAndTinySamplesUnchanged) {
    // All-identical input still has exactly zero entropy.
    EXPECT_DOUBLE_EQ(compute_shannon_entropy(std::vector<double>(50, 7.0), 20), 0.0);
    // Below the quartile-sanity floor the fence never engages; a 2-sample
    // input must not throw or produce NaN.
    const double H = compute_shannon_entropy({0.0, 1.0e6}, 20);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
}

TEST(RobustBinning, CollapseThresholdScalesWithSupport) {
    // The absolute kHSC_* constants were derived for 256 bins but are applied
    // to a 20-bin estimator; the helper expresses the same intent at any bin
    // count. ln(256)/4 = 1.386 nats (= 2 bits) recovers the original line.
    EXPECT_NEAR(collapse_threshold_nats(256, kHSC_soft_frac_of_max),
                kHSC_soft_nats, 1e-12);
    EXPECT_NEAR(collapse_threshold_nats(256, kHSC_hard_frac_of_max),
                kHSC_hard_nats, 1e-12);
    // At the bin count actually used, the scale-correct gate is stricter than
    // the shipped absolute constant.
    EXPECT_LT(collapse_threshold_nats(DEFAULT_HIST_BINS, kHSC_soft_frac_of_max),
              kHSC_soft_nats);
    EXPECT_DOUBLE_EQ(collapse_threshold_nats(1), 0.0);
}

// ===========================================================================
// ONE ESTIMATOR, TWO DISPATCH SHELLS
//
// shannon_thermo::compute_shannon_entropy and
// UnifiedHardwareDispatch::compute_shannon_entropy are separate dispatch
// shells over the same quantity. They each used to derive their own support
// and bin index, so once one gained the robust fence the two disagreed by
// ~2 nats on the same input — enough to land on opposite sides of the collapse
// gate. Both now share ShannonBinning.h; these tests pin that they agree.
// ===========================================================================

TEST(EstimatorParity, AgreesWithUnifiedDispatchOnCleanData) {
    std::mt19937 rng(2024);
    std::normal_distribution<double> dist(-50.0, 8.0);
    std::vector<double> values(1000);
    for (auto& v : values) v = dist(rng);

    const double H_stack = compute_shannon_entropy(values, DEFAULT_HIST_BINS);
    const double H_disp  = hw::UnifiedHardwareDispatch::instance()
                               .compute_shannon_entropy(values, DEFAULT_HIST_BINS,
                                                        hw::Backend::SCALAR);
    EXPECT_NEAR(H_stack, H_disp, 1e-12);
}

TEST(EstimatorParity, AgreesWithUnifiedDispatchWhenTheFenceEngages) {
    // The case that used to split them: the fence engages for one estimator
    // only, so they straddled the collapse gate.
    std::vector<double> energies;
    energies.reserve(100);
    for (int i = 0; i < 99; ++i)
        energies.push_back(-100.0 + 0.5 * static_cast<double>(i));
    energies.push_back(1.0e4);

    const double H_stack = compute_shannon_entropy(energies, DEFAULT_HIST_BINS);
    const double H_disp  = hw::UnifiedHardwareDispatch::instance()
                               .compute_shannon_entropy(energies, DEFAULT_HIST_BINS,
                                                        hw::Backend::SCALAR);
    EXPECT_NEAR(H_stack, H_disp, 1e-12);
    EXPECT_GT(H_stack, kHSC_soft_nats);  // and both clear the gate
}

// ===========================================================================
// TORSIONAL ENM HESSIAN — DIRECTIONAL PROJECTION
//
// The elastic-network potential is a sum of springs ALONG each contact
// vector, V = ½ Σ k_ij [û_ij·(δr_i − δr_j)]². Projecting onto û_ij is what
// makes the model rotationally invariant, so motions that do not change any
// contact LENGTH must cost no energy.
// ===========================================================================

TEST(TorsionalHessian, SpectrumIsPhysicallyAdmissible) {
    tencm::TorsionalENM model = make_built_tencm_model(30);
    ASSERT_TRUE(model.is_built());
    ASSERT_FALSE(model.modes().empty());

    for (const auto& m : model.modes()) {
        EXPECT_TRUE(std::isfinite(m.eigenvalue));
        // A sum of k·(p⊗p) outer products is positive semi-definite; small
        // negative values are numerical noise only.
        EXPECT_GT(m.eigenvalue, -1e-6);
    }
}

TEST(TorsionalHessian, RigidFragmentRotationCostsNothing) {
    // Decisive test of the directional projection.
    //
    // Rotating about the FIRST pseudo-bond turns the entire downstream
    // fragment as a rigid body. If the first residue is placed beyond the
    // contact cutoff from every other residue, that rotation changes no
    // contact LENGTH anywhere in the structure, so its true elastic energy is
    // exactly zero and the Hessian must be singular.
    //
    // Algebraically: for a contact whose two atoms are both downstream of bond
    // k, ΔJ_k = e_k × r_ij, and û_ij·(e_k × r_ij) ≡ 0. The projected form
    // scores exactly zero; the unprojected |ΔJ|² form scores |e_k × r_ij|²,
    // which is large and positive — it charges energy for a pure rotation.
    std::vector<std::array<float, 3>> ca;
    // Residue 0 parked far away: no contacts within the 9 Å default cutoff.
    ca.push_back({0.0f, 0.0f, -500.0f});
    // A compact helical chain for the remaining residues.
    for (int i = 0; i < 24; ++i) {
        ca.push_back({
            2.3f * std::cos(static_cast<float>(i) * 1.74532925f),
            2.3f * std::sin(static_cast<float>(i) * 1.74532925f),
            1.5f * static_cast<float>(i)
        });
    }

    tencm::TorsionalENM model;
    model.build_from_ca(ca);
    ASSERT_TRUE(model.is_built());
    ASSERT_GE(model.modes().size(), 2u);

    // Every contact here is between two atoms that are both downstream of bond
    // 0, and residue 0 has no contacts of its own. With the projection, bond
    // 0's entire row of the Hessian is exactly zero, so its diagonal stiffness
    // H_00 — the energy cost of rotating about bond 0 alone — must vanish.
    // Without the projection, bond 0 instead picks up |e_0 × r_ij|² from every
    // one of those contacts.
    //
    // H is private, but the spectral decomposition recovers any element:
    //     H_kk = Σ_m λ_m · v_m[k]²
    // The null space here is multi-dimensional, so this diagonal form is used
    // rather than inspecting a single (basis-arbitrary) null eigenvector.
    ASSERT_FALSE(model.modes().front().eigenvector.empty())
        << "build_from_ca must populate eigenvectors for this check";

    const int M = model.n_bonds();
    std::vector<double> diag(static_cast<std::size_t>(M), 0.0);
    for (const auto& m : model.modes()) {
        if (static_cast<int>(m.eigenvector.size()) < M) continue;
        for (int k = 0; k < M; ++k)
            diag[static_cast<std::size_t>(k)] +=
                m.eigenvalue * m.eigenvector[static_cast<std::size_t>(k)] *
                m.eigenvector[static_cast<std::size_t>(k)];
    }

    double max_diag = 0.0;
    for (double d : diag) max_diag = std::max(max_diag, std::abs(d));
    ASSERT_GT(max_diag, 0.0) << "degenerate Hessian — test structure is wrong";

    EXPECT_LT(std::abs(diag[0]) / max_diag, 1e-9)
        << "H_00=" << diag[0] << " (max diagonal " << max_diag
        << ") — rotating about bond 0 turns an isolated fragment rigidly and"
           " changes no contact length, so it must cost exactly zero energy."
           " A nonzero H_00 is the signature of a Hessian missing the û_ij"
           " projection.";
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
