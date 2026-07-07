// tests/test_hardware_detect_dispatch.cpp
// Unit tests for hardware_detect.h and hardware_dispatch.h
//
// Validates:
//   1. Runtime hardware detection (CPUID, OpenMP, Eigen presence)
//   2. Backend selection logic (priority ordering)
//   3. Boltzmann batch computation (correctness across all backends)
//   4. Log-sum-exp dispatch (numerical stability)
//   5. Dispatch report generation
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/hardware_detect.h"
#include "../LIB/UnifiedHardwareDispatch.h"
#include "../LIB/simd_distance.h"

#include <cmath>
#include <vector>
#include <numeric>
#include <random>
#include <limits>
#include <cstring>

using namespace flexaids;

static constexpr double EPSILON = 1e-10;

// ===========================================================================
// HARDWARE DETECTION
// ===========================================================================

TEST(HardwareDetect, DetectReturnsValid) {
    const auto& hw = detect_hardware();
    // Basic sanity: OpenMP threads >= 1
    EXPECT_GE(hw.openmp_max_threads, 1);
}

TEST(HardwareDetect, CachedAcrossCalls) {
    const auto& hw1 = detect_hardware();
    const auto& hw2 = detect_hardware();
    EXPECT_EQ(&hw1, &hw2);  // same static instance
}

TEST(HardwareDetect, SummaryIsNonEmpty) {
    const auto& hw = detect_hardware();
    std::string s = hw.summary();
    EXPECT_FALSE(s.empty());
    EXPECT_NE(s.find("[HW]"), std::string::npos);
}

#if defined(__AVX2__)
TEST(HardwareDetect, DetectsAVX2) {
    const auto& hw = detect_hardware();
    EXPECT_TRUE(hw.has_avx2);
}
#endif

#if defined(__AVX512F__)
TEST(HardwareDetect, DetectsAVX512) {
    const auto& hw = detect_hardware();
    EXPECT_TRUE(hw.has_avx512f);
    // If __AVX512F__ is defined at compile time, CPUID should also report it
}
#endif

#ifdef _OPENMP
TEST(HardwareDetect, DetectsOpenMP) {
    const auto& hw = detect_hardware();
    EXPECT_TRUE(hw.has_openmp);
    EXPECT_GE(hw.openmp_max_threads, 1);
}
#endif

TEST(HardwareDetect, DetectsEigen) {
    const auto& hw = detect_hardware();
    EXPECT_TRUE(hw.has_eigen);
}

// ===========================================================================
// BACKEND SELECTION
// ===========================================================================

TEST(BackendSelection, SelectsValidBackend) {
    HardwareBackend b = select_backend();
    EXPECT_GE(static_cast<int>(b), 0);
    EXPECT_LE(static_cast<int>(b), 6);
}

TEST(BackendSelection, DefaultFitnessBackendIsNotGPU) {
    HardwareBackend b = select_backend();
    EXPECT_NE(b, HardwareBackend::CUDA);
    EXPECT_NE(b, HardwareBackend::ROCM);
    EXPECT_NE(b, HardwareBackend::METAL);
}

TEST(BackendSelection, GpuOverrideDoesNotApplyToFitnessEval) {
    auto& d = hw::UnifiedHardwareDispatch::instance();
    d.detect();

    d.set_override(hw::Backend::CUDA);
    EXPECT_EQ(d.best_backend(hw::KernelType::FITNESS_EVAL), d.select_cpu_backend());

    d.set_override(hw::Backend::ROCM);
    EXPECT_EQ(d.best_backend(hw::KernelType::FITNESS_EVAL), d.select_cpu_backend());

    d.set_override(hw::Backend::METAL);
    EXPECT_EQ(d.best_backend(hw::KernelType::FITNESS_EVAL), d.select_cpu_backend());

    d.clear_override();
}

TEST(BackendSelection, CPUBackendIsNotGPU) {
    HardwareBackend b = select_cpu_backend();
    EXPECT_NE(b, HardwareBackend::CUDA);
    EXPECT_NE(b, HardwareBackend::ROCM);
    EXPECT_NE(b, HardwareBackend::METAL);
}

TEST(BackendSelection, BackendNameValid) {
    for (int i = 0; i <= 6; ++i) {
        const char* name = backend_name(static_cast<HardwareBackend>(i));
        EXPECT_NE(name, nullptr);
        EXPECT_GT(std::strlen(name), 0u);
    }
}

TEST(BackendSelection, BackendNameCoversAll) {
    EXPECT_STREQ(backend_name(HardwareBackend::CUDA), "CUDA");
    EXPECT_STREQ(backend_name(HardwareBackend::ROCM), "ROCm");
    EXPECT_STREQ(backend_name(HardwareBackend::METAL), "Metal");
    EXPECT_STREQ(backend_name(HardwareBackend::AVX512), "AVX-512");
    EXPECT_STREQ(backend_name(HardwareBackend::AVX2), "AVX2");
    EXPECT_STREQ(backend_name(HardwareBackend::OPENMP), "OpenMP");
    EXPECT_STREQ(backend_name(HardwareBackend::SCALAR), "scalar");
}

// ===========================================================================
// BOLTZMANN BATCH COMPUTATION
// ===========================================================================

TEST(BoltzmannBatch, EmptyInput) {
    auto result = compute_boltzmann_batch({}, 1.0);
    EXPECT_TRUE(result.weights.empty());
}

TEST(BoltzmannBatch, SingleEnergy) {
    std::vector<double> E = {-10.0};
    double beta = 1.0 / (0.001987206 * 300.0);
    auto result = compute_boltzmann_batch(E, beta);

    ASSERT_EQ(result.weights.size(), 1u);
    // Single energy: weight = exp(0) = 1.0 (E - E_min = 0)
    EXPECT_NEAR(result.weights[0], 1.0, EPSILON);
    EXPECT_DOUBLE_EQ(result.E_min, -10.0);
}

TEST(BoltzmannBatch, TwoStatesRelativeWeights) {
    // Two states: E1 = -10, E2 = -8
    // w1/w2 = exp(-beta * (E1 - E2)) = exp(-beta * (-2)) = exp(2*beta)
    double T = 300.0;
    double beta = 1.0 / (0.001987206 * T);
    std::vector<double> E = {-10.0, -8.0};

    auto result = compute_boltzmann_batch(E, beta);
    ASSERT_EQ(result.weights.size(), 2u);

    double expected_ratio = std::exp(beta * 2.0);  // w1/w2
    double actual_ratio = result.weights[0] / result.weights[1];
    EXPECT_NEAR(actual_ratio, expected_ratio, expected_ratio * 1e-10);
}

TEST(BoltzmannBatch, WeightsArePositive) {
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(-15.0, 5.0);
    std::vector<double> E(1000);
    for (auto& e : E) e = dist(rng);

    double beta = 1.0 / (0.001987206 * 300.0);
    auto result = compute_boltzmann_batch(E, beta);

    for (double w : result.weights) {
        EXPECT_GT(w, 0.0);
        EXPECT_TRUE(std::isfinite(w));
    }
}

TEST(BoltzmannBatch, PartitionFunctionIsFinite) {
    std::mt19937 rng(77);
    std::normal_distribution<double> dist(-12.0, 3.0);
    std::vector<double> E(500);
    for (auto& e : E) e = dist(rng);

    double beta = 1.0 / (0.001987206 * 300.0);
    auto result = compute_boltzmann_batch(E, beta);

    EXPECT_TRUE(std::isfinite(result.log_Z));
}

TEST(BoltzmannBatch, LargeEnergyDifferences) {
    // Test numerical stability with 100 kT energy range
    double T = 300.0;
    double kT = 0.001987206 * T;
    double beta = 1.0 / kT;

    std::vector<double> E;
    for (int i = 0; i < 100; ++i)
        E.push_back(-10.0 - i * kT);  // span 100 kT

    auto result = compute_boltzmann_batch(E, beta);

    // All weights should be finite
    for (double w : result.weights)
        EXPECT_TRUE(std::isfinite(w));
    EXPECT_TRUE(std::isfinite(result.log_Z));
}

TEST(BoltzmannBatch, DegenerateEnergies) {
    // All same energy → all weights equal
    std::vector<double> E(100, -10.0);
    double beta = 1.0 / (0.001987206 * 300.0);

    auto result = compute_boltzmann_batch(E, beta);

    for (double w : result.weights)
        EXPECT_NEAR(w, 1.0, EPSILON);
}

// ===========================================================================
// LOG-SUM-EXP DISPATCH
// ===========================================================================

TEST(LogSumExp, EmptyReturnsNegInf) {
    double r = log_sum_exp_dispatch({});
    EXPECT_EQ(r, -std::numeric_limits<double>::infinity());
}

TEST(LogSumExp, SingleValue) {
    std::vector<double> v = {5.0};
    EXPECT_NEAR(log_sum_exp_dispatch(v), 5.0, EPSILON);
}

TEST(LogSumExp, TwoEqualValues) {
    // log(exp(3) + exp(3)) = 3 + log(2)
    std::vector<double> v = {3.0, 3.0};
    EXPECT_NEAR(log_sum_exp_dispatch(v), 3.0 + std::log(2.0), EPSILON);
}

TEST(LogSumExp, MatchesNaiveForSmallValues) {
    std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0};
    double naive = 0.0;
    for (double x : v) naive += std::exp(x);
    naive = std::log(naive);

    EXPECT_NEAR(log_sum_exp_dispatch(v), naive, 1e-12);
}

TEST(LogSumExp, NumericalStabilityLargeValues) {
    // exp(1000) overflows, but log-sum-exp should handle it
    std::vector<double> v = {1000.0, 999.0, 998.0};
    double result = log_sum_exp_dispatch(v);
    EXPECT_TRUE(std::isfinite(result));

    // Should be close to 1000 + log(1 + exp(-1) + exp(-2))
    double expected = 1000.0 + std::log(1.0 + std::exp(-1.0) + std::exp(-2.0));
    EXPECT_NEAR(result, expected, 1e-10);
}

TEST(LogSumExp, NumericalStabilityNegativeValues) {
    std::vector<double> v = {-1000.0, -1001.0, -1002.0};
    double result = log_sum_exp_dispatch(v);
    EXPECT_TRUE(std::isfinite(result));

    double expected = -1000.0 + std::log(1.0 + std::exp(-1.0) + std::exp(-2.0));
    EXPECT_NEAR(result, expected, 1e-10);
}

TEST(LogSumExp, LargeArray) {
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 10.0);
    std::vector<double> v(10000);
    for (auto& x : v) x = dist(rng);

    double result = log_sum_exp_dispatch(v);
    EXPECT_TRUE(std::isfinite(result));
}

// ===========================================================================
// DISPATCH REPORT
// ===========================================================================

TEST(DispatchReport, HasValidBackend) {
    auto report = get_dispatch_report();
    EXPECT_GE(static_cast<int>(report.selected), 0);
    EXPECT_LE(static_cast<int>(report.selected), 5);
}

TEST(DispatchReport, ReasonIsNonEmpty) {
    auto report = get_dispatch_report();
    EXPECT_FALSE(report.reason.empty());
}

TEST(DispatchReport, HWSummaryContainsMarker) {
    auto report = get_dispatch_report();
    EXPECT_NE(report.hw_summary.find("[HW]"), std::string::npos);
}

// ===========================================================================
// SIMD PRIMITIVES (AVX-512 when available, else AVX2/scalar)
// ===========================================================================

TEST(SimdDistance, RMSDIdenticalArraysIsZero) {
    std::vector<float> a = {1, 2, 3, 4, 5, 6, 7, 8, 9};  // 3 atoms
    float r = simd::rmsd(a.data(), a.data(), 3);
    EXPECT_NEAR(r, 0.0f, 1e-6f);
}

TEST(SimdDistance, RMSDKnownValue) {
    // 2 atoms: a = (0,0,0), (3,0,0); b = (0,0,0), (0,0,0)
    std::vector<float> a = {0, 0, 0, 3, 0, 0};
    std::vector<float> b = {0, 0, 0, 0, 0, 0};
    float r = simd::rmsd(a.data(), b.data(), 2);
    // sqrt((0 + 9) / 2) = sqrt(4.5) ≈ 2.121
    EXPECT_NEAR(r, std::sqrt(4.5f), 1e-5f);
}

TEST(SimdDistance, Distance2Scalar) {
    float a[3] = {1.0f, 2.0f, 3.0f};
    float b[3] = {4.0f, 6.0f, 3.0f};
    float d2 = simd::distance2_scalar(a, b);
    EXPECT_NEAR(d2, 9.0f + 16.0f + 0.0f, 1e-6f);
}

TEST(SimdDistance, LargeArrayRMSD) {
    // 100 atoms: stress test SIMD tail handling
    const int N = 100;
    std::vector<float> a(N * 3), b(N * 3);
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-10.0f, 10.0f);
    for (int i = 0; i < N * 3; ++i) {
        a[i] = dist(rng);
        b[i] = dist(rng);
    }

    float r = simd::rmsd(a.data(), b.data(), N);
    EXPECT_GT(r, 0.0f);
    EXPECT_TRUE(std::isfinite(r));

    // Verify against scalar
    float sum = 0;
    for (int i = 0; i < N; ++i)
        for (int c = 0; c < 3; ++c)
            sum += simd::sq(a[i*3+c] - b[i*3+c]);
    float expected = std::sqrt(sum / N);
    EXPECT_NEAR(r, expected, 1e-3f);
}

#if FLEXAIDS_HAS_AVX512
TEST(SimdAVX512, Distance2_1x16) {
    float ax[16], ay[16], az[16], out[16];
    for (int i = 0; i < 16; ++i) {
        ax[i] = static_cast<float>(i);
        ay[i] = 0.0f;
        az[i] = 0.0f;
    }
    simd::distance2_1x16(ax, ay, az, 0.0f, 0.0f, 0.0f, out);
    for (int i = 0; i < 16; ++i)
        EXPECT_NEAR(out[i], static_cast<float>(i * i), 1e-4f);
}

TEST(SimdAVX512, LJWall16x) {
    float r2[16], Ewall[16];
    for (int i = 0; i < 16; ++i) r2[i] = 1.0f + 0.1f * i;
    simd::lj_wall_16x(r2, 0.5f, 1.0f, Ewall);
    for (int i = 0; i < 16; ++i) {
        EXPECT_TRUE(std::isfinite(Ewall[i]));
    }
}
#endif

// ===========================================================================
// MAIN
// ===========================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
