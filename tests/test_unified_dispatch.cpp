// tests/test_unified_dispatch.cpp
// Unit tests for the HardwareDispatcher runtime dispatch layer.
// Tests backend detection, kernel routing, edge cases, override, and
// numerical consistency across backends.
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#include <gtest/gtest.h>
#include "../LIB/HardwareDispatch.h"

#include <cmath>
#include <numeric>
#include <random>
#include <vector>

static constexpr double EPSILON = 1e-6;

// ===========================================================================
// HARDWARE DETECTION
// ===========================================================================

TEST(HardwareDetection, DetectDoesNotCrash) {
    auto& d = hw::HardwareDispatcher::instance();
    EXPECT_NO_THROW(d.detect());
}

TEST(HardwareDetection, ScalarAlwaysAvailable) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    EXPECT_TRUE(d.is_available(hw::Backend::SCALAR));
}

TEST(HardwareDetection, AutoAlwaysAvailable) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    EXPECT_TRUE(d.is_available(hw::Backend::AUTO));
}

TEST(HardwareDetection, AvailableBackendsNotEmpty) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    auto avail = d.available_backends();
    ASSERT_FALSE(avail.empty());
    // Last element should always be SCALAR
    EXPECT_EQ(avail.back(), hw::Backend::SCALAR);
}

TEST(HardwareDetection, HardwareReportIsNonEmpty) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    auto report = d.hardware_report();
    EXPECT_GT(report.size(), 50u);
    EXPECT_NE(report.find("CPU:"), std::string::npos);
    EXPECT_NE(report.find("Cores:"), std::string::npos);
}

TEST(HardwareDetection, BackendNameIsValid) {
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::SCALAR), "scalar");
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::OPENMP), "OpenMP");
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::AVX2),   "AVX2");
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::AVX512), "AVX-512");
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::METAL),  "Metal");
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::CUDA),   "CUDA");
    EXPECT_STREQ(hw::HardwareDispatcher::backend_name(hw::Backend::AUTO),   "auto");
}

// ===========================================================================
// BACKEND OVERRIDE
// ===========================================================================

TEST(BackendOverride, OverrideAndClear) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();

    EXPECT_EQ(d.current_override(), hw::Backend::AUTO);
    d.set_override(hw::Backend::SCALAR);
    EXPECT_EQ(d.current_override(), hw::Backend::SCALAR);
    // best_backend should return the override
    EXPECT_EQ(d.best_backend(hw::KernelType::SHANNON_ENTROPY), hw::Backend::SCALAR);
    d.clear_override();
    EXPECT_EQ(d.current_override(), hw::Backend::AUTO);
}

// ===========================================================================
// SHANNON ENTROPY — DISPATCH CORRECTNESS
// ===========================================================================

class ShannonDispatchTest : public ::testing::Test {
protected:
    hw::HardwareDispatcher& d = hw::HardwareDispatcher::instance();
    std::vector<double> data;

    void SetUp() override {
        d.detect();
        d.clear_override();
        std::mt19937 rng(42);
        std::normal_distribution<double> dist(-10.0, 5.0);
        data.resize(5000);
        for (auto& v : data) v = dist(rng);
    }
};

TEST_F(ShannonDispatchTest, EmptyReturnsZero) {
    std::vector<double> empty;
    EXPECT_DOUBLE_EQ(d.compute_shannon_entropy(empty), 0.0);
}

TEST_F(ShannonDispatchTest, ScalarProducesNonNegative) {
    double H = d.compute_shannon_entropy(data, 20, hw::Backend::SCALAR);
    EXPECT_GE(H, 0.0);
    EXPECT_TRUE(std::isfinite(H));
}

TEST_F(ShannonDispatchTest, AutoMatchesScalar) {
    double H_scalar = d.compute_shannon_entropy(data, 20, hw::Backend::SCALAR);
    double H_auto   = d.compute_shannon_entropy(data, 20, hw::Backend::AUTO);
    // They may use different code paths, but should agree within tolerance
    EXPECT_NEAR(H_auto, H_scalar, 0.01);
}

TEST_F(ShannonDispatchTest, OpenMPMatchesScalar) {
    if (!d.is_available(hw::Backend::OPENMP)) GTEST_SKIP();
    double H_scalar = d.compute_shannon_entropy(data, 20, hw::Backend::SCALAR);
    double H_omp    = d.compute_shannon_entropy(data, 20, hw::Backend::OPENMP);
    EXPECT_NEAR(H_omp, H_scalar, EPSILON);
}

TEST_F(ShannonDispatchTest, AVX512MatchesScalar) {
    if (!d.is_available(hw::Backend::AVX512)) GTEST_SKIP();
    double H_scalar = d.compute_shannon_entropy(data, 20, hw::Backend::SCALAR);
    double H_avx    = d.compute_shannon_entropy(data, 20, hw::Backend::AVX512);
    EXPECT_NEAR(H_avx, H_scalar, EPSILON);
}

TEST_F(ShannonDispatchTest, Reproducible) {
    double H1 = d.compute_shannon_entropy(data, 20);
    double H2 = d.compute_shannon_entropy(data, 20);
    EXPECT_DOUBLE_EQ(H1, H2);
}

TEST_F(ShannonDispatchTest, NegativeBinsHandled) {
    double H = d.compute_shannon_entropy(data, -5);
    EXPECT_GE(H, 0.0);
    EXPECT_TRUE(std::isfinite(H));
}

TEST_F(ShannonDispatchTest, SingleValueReturnsZero) {
    std::vector<double> single = {3.14};
    EXPECT_DOUBLE_EQ(d.compute_shannon_entropy(single, 10), 0.0);
}

TEST_F(ShannonDispatchTest, IdenticalValuesReturnZero) {
    std::vector<double> same(100, -5.0);
    EXPECT_DOUBLE_EQ(d.compute_shannon_entropy(same, 10), 0.0);
}

TEST_F(ShannonDispatchTest, UpperBound) {
    int bins = 20;
    double H = d.compute_shannon_entropy(data, bins);
    EXPECT_LE(H, std::log(static_cast<double>(bins)) + EPSILON);
}

// ===========================================================================
// LOG-SUM-EXP — DISPATCH CORRECTNESS
// ===========================================================================

TEST(LogSumExpDispatch, EmptyReturnsLargeNegative) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<double> empty;
    double result = d.log_sum_exp(empty);
    EXPECT_LE(result, -1e300);
}

TEST(LogSumExpDispatch, SingleElement) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<double> single = {5.0};
    double result = d.log_sum_exp(single);
    EXPECT_NEAR(result, 5.0, EPSILON);
}

TEST(LogSumExpDispatch, KnownResult) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    // log(exp(0) + exp(0)) = log(2) ≈ 0.6931
    std::vector<double> vals = {0.0, 0.0};
    double result = d.log_sum_exp(vals);
    EXPECT_NEAR(result, std::log(2.0), EPSILON);
}

TEST(LogSumExpDispatch, NumericalStability) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    // Very large values should not overflow
    std::vector<double> large = {1000.0, 1001.0, 1002.0};
    double result = d.log_sum_exp(large);
    EXPECT_TRUE(std::isfinite(result));
    EXPECT_GT(result, 1001.0);
}

TEST(LogSumExpDispatch, ConsistentAcrossBackends) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(-10.0, 5.0);
    std::vector<double> data(1000);
    for (auto& v : data) v = dist(rng);

    double lse_scalar = d.log_sum_exp(data, hw::Backend::SCALAR);
    double lse_auto   = d.log_sum_exp(data, hw::Backend::AUTO);
    EXPECT_NEAR(lse_auto, lse_scalar, 1e-10);
}

// ===========================================================================
// BOLTZMANN WEIGHTS — DISPATCH CORRECTNESS
// ===========================================================================

TEST(BoltzmannDispatch, EmptyReturnsEmpty) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    auto w = d.compute_boltzmann_weights({}, 1.0);
    EXPECT_TRUE(w.empty());
}

TEST(BoltzmannDispatch, SumToOne) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<double> E = {-15.0, -12.0, -10.0, -8.0, -5.0};
    double beta = 1.0 / (0.001987206 * 298.15);
    auto w = d.compute_boltzmann_weights(E, beta);
    double sum = std::accumulate(w.begin(), w.end(), 0.0);
    EXPECT_NEAR(sum, 1.0, 1e-10);
}

TEST(BoltzmannDispatch, LowestEnergyHighestWeight) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<double> E = {-15.0, -12.0, -10.0};
    double beta = 1.0 / (0.001987206 * 298.15);
    auto w = d.compute_boltzmann_weights(E, beta);
    ASSERT_EQ(w.size(), 3u);
    EXPECT_GT(w[0], w[1]);
    EXPECT_GT(w[1], w[2]);
}

TEST(BoltzmannDispatch, AllWeightsNonNegative) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(-15.0, 5.0);
    std::vector<double> E(1000);
    for (auto& e : E) e = dist(rng);
    double beta = 1.0 / (0.001987206 * 298.15);
    auto w = d.compute_boltzmann_weights(E, beta);
    for (double wi : w) {
        EXPECT_GE(wi, 0.0);
        EXPECT_TRUE(std::isfinite(wi));
    }
}

// ===========================================================================
// RMSD — DISPATCH CORRECTNESS
// ===========================================================================

TEST(RMSDDispatch, ZeroAtoms) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    EXPECT_FLOAT_EQ(d.rmsd(nullptr, nullptr, 0), 0.0f);
}

TEST(RMSDDispatch, IdenticalCoordsZeroRMSD) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<float> c = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    float r = d.rmsd(c.data(), c.data(), 2);
    EXPECT_FLOAT_EQ(r, 0.0f);
}

TEST(RMSDDispatch, KnownValue) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    // Two atoms: A at (0,0,0),(0,0,0); B at (1,0,0),(0,1,0)
    // d^2 = 1 + 1 = 2, RMSD = sqrt(2/2) = 1.0
    std::vector<float> a = {0, 0, 0, 0, 0, 0};
    std::vector<float> b = {1, 0, 0, 0, 1, 0};
    float r = d.rmsd(a.data(), b.data(), 2);
    EXPECT_NEAR(r, 1.0f, 1e-5f);
}

TEST(RMSDDispatch, ConsistentAcrossBackends) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-50.0f, 50.0f);
    int n = 500;
    std::vector<float> a(n * 3), b(n * 3);
    for (auto& x : a) x = dist(rng);
    for (auto& x : b) x = dist(rng);

    float r_scalar = d.rmsd(a.data(), b.data(), n, hw::Backend::SCALAR);
    float r_auto   = d.rmsd(a.data(), b.data(), n, hw::Backend::AUTO);
    EXPECT_NEAR(r_auto, r_scalar, 0.01f);
}

// ===========================================================================
// DISTANCE2 BATCH — DISPATCH CORRECTNESS
// ===========================================================================

TEST(Distance2Dispatch, ScalarCorrectness) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<float> ax = {1, 2, 3, 4};
    std::vector<float> ay = {0, 0, 0, 0};
    std::vector<float> az = {0, 0, 0, 0};
    std::vector<float> out(4);

    d.distance2_batch(ax.data(), ay.data(), az.data(),
                       0.0f, 0.0f, 0.0f, out.data(), 4, hw::Backend::SCALAR);
    EXPECT_NEAR(out[0], 1.0f, 1e-5f);
    EXPECT_NEAR(out[1], 4.0f, 1e-5f);
    EXPECT_NEAR(out[2], 9.0f, 1e-5f);
    EXPECT_NEAR(out[3], 16.0f, 1e-5f);
}

TEST(Distance2Dispatch, ConsistentAcrossBackends) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    int n = 1024;
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-50.0f, 50.0f);
    std::vector<float> ax(n), ay(n), az(n);
    for (int i = 0; i < n; ++i) {
        ax[i] = dist(rng);
        ay[i] = dist(rng);
        az[i] = dist(rng);
    }
    std::vector<float> out_scalar(n), out_auto(n);

    d.distance2_batch(ax.data(), ay.data(), az.data(),
                       1.0f, 2.0f, 3.0f, out_scalar.data(), n, hw::Backend::SCALAR);
    d.distance2_batch(ax.data(), ay.data(), az.data(),
                       1.0f, 2.0f, 3.0f, out_auto.data(), n, hw::Backend::AUTO);

    for (int i = 0; i < n; ++i)
        EXPECT_NEAR(out_auto[i], out_scalar[i], 0.01f) << "Mismatch at i=" << i;
}

// ===========================================================================
// EDGE CASES
// ===========================================================================

TEST(DispatchEdgeCases, LargeDataset) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(-10.0, 5.0);
    std::vector<double> big(100000);
    for (auto& v : big) v = dist(rng);

    double H = d.compute_shannon_entropy(big, 50);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GT(H, 0.0);
}

TEST(DispatchEdgeCases, VeryLargeEnergies) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<double> extreme(100);
    for (int i = 0; i < 100; ++i)
        extreme[i] = 1e6 + i * 0.01;

    double H = d.compute_shannon_entropy(extreme, 10);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
}

TEST(DispatchEdgeCases, MixedSignEnergies) {
    auto& d = hw::HardwareDispatcher::instance();
    d.detect();
    std::vector<double> mixed = {-100.0, -50.0, 0.0, 50.0, 100.0};
    double H = d.compute_shannon_entropy(mixed, 5);
    EXPECT_TRUE(std::isfinite(H));
    EXPECT_GE(H, 0.0);
}

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
