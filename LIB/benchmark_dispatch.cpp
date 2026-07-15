// benchmark_dispatch.cpp — Benchmark suite for FlexAIDdS hardware dispatch
//
// Measures throughput and speedup for each available backend across
// key kernels: Shannon entropy, log-sum-exp, Boltzmann weights, RMSD.
//
// Build:
//   cmake -DENABLE_DISPATCH_BENCHMARK=ON -DFLEXAIDS_USE_AVX512=ON \
//         -DFLEXAIDS_USE_OPENMP=ON -DFLEXAIDS_USE_EIGEN=ON ..
//   cmake --build . --target benchmark_dispatch
//
// Run:
//   ./benchmark_dispatch [--size N] [--reps R]
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#include "UnifiedHardwareDispatch.h"  // modern replacement (HardwareDispatch.h is deprecated)
#include "simd_distance.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// ─── timing helper ───────────────────────────────────────────────────────────

struct Timer {
    using clock = std::chrono::high_resolution_clock;
    clock::time_point t0;

    void start() { t0 = clock::now(); }
    double elapsed_ms() const {
        return std::chrono::duration<double, std::milli>(clock::now() - t0).count();
    }
    double elapsed_s() const {
        return std::chrono::duration<double>(clock::now() - t0).count();
    }
};

// ─── test data generation ────────────────────────────────────────────────────

static std::vector<double> random_doubles(int n, double mean, double stddev, int seed = 42) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(mean, stddev);
    std::vector<double> v(n);
    for (auto& x : v) x = dist(rng);
    return v;
}

static std::vector<float> random_coords(int n_atoms, int seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-50.0f, 50.0f);
    std::vector<float> v(n_atoms * 3);
    for (auto& x : v) x = dist(rng);
    return v;
}

// ─── benchmark runners ───────────────────────────────────────────────────────

struct BenchResult {
    std::string kernel;
    std::string backend;
    int         data_size;
    double      time_ms;
    double      throughput;  // items/sec
    double      speedup;     // vs scalar
};

static void print_header() {
    std::cout << std::left
              << std::setw(25) << "Kernel"
              << std::setw(14) << "Backend"
              << std::setw(12) << "Size"
              << std::setw(14) << "Time(ms)"
              << std::setw(16) << "Throughput"
              << std::setw(10) << "Speedup"
              << "\n";
    std::cout << std::string(91, '-') << "\n";
}

static void print_row(const BenchResult& r) {
    std::cout << std::left
              << std::setw(25) << r.kernel
              << std::setw(14) << r.backend
              << std::setw(12) << r.data_size
              << std::fixed << std::setprecision(3)
              << std::setw(14) << r.time_ms
              << std::setprecision(0)
              << std::setw(16) << r.throughput
              << std::setprecision(2)
              << std::setw(10) << r.speedup
              << "\n";
}

static void bench_shannon(hw::UnifiedHardwareDispatch& disp, int n, int reps,
                           std::vector<BenchResult>& results) {
    auto data = random_doubles(n, -10.0, 5.0);
    int bins = 20;

    // Scalar baseline
    Timer t;
    t.start();
    volatile double sink = 0;
    for (int r = 0; r < reps; ++r)
        sink = disp.compute_shannon_entropy(data, bins, hw::Backend::SCALAR);
    double scalar_ms = t.elapsed_ms();
    results.push_back({"ShannonEntropy", "scalar", n, scalar_ms / reps,
                        n * reps / (scalar_ms * 1e-3), 1.0});
    (void)sink;

    // OpenMP
    if (disp.is_available(hw::Backend::OPENMP)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.compute_shannon_entropy(data, bins, hw::Backend::OPENMP);
        double ms = t.elapsed_ms();
        results.push_back({"ShannonEntropy", "OpenMP", n, ms / reps,
                            n * reps / (ms * 1e-3), scalar_ms / ms});
    }

    // AVX-512
    if (disp.is_available(hw::Backend::AVX512)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.compute_shannon_entropy(data, bins, hw::Backend::AVX512);
        double ms = t.elapsed_ms();
        results.push_back({"ShannonEntropy", "AVX-512", n, ms / reps,
                            n * reps / (ms * 1e-3), scalar_ms / ms});
    }

    if (disp.is_available(hw::Backend::METAL)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.compute_shannon_entropy(data, bins, hw::Backend::METAL);
        double ms = t.elapsed_ms();
        results.push_back({"ShannonEntropy", "Metal", n, ms / reps,
                            n * reps / (ms * 1e-3), scalar_ms / ms});
    }
}

static void bench_lse(hw::UnifiedHardwareDispatch& disp, int n, int reps,
                       std::vector<BenchResult>& results) {
    auto data = random_doubles(n, -10.0, 5.0);

    Timer t;
    t.start();
    volatile double sink = 0;
    for (int r = 0; r < reps; ++r)
        sink = disp.log_sum_exp(data, hw::Backend::SCALAR);
    double scalar_ms = t.elapsed_ms();
    results.push_back({"LogSumExp", "scalar", n, scalar_ms / reps,
                        n * reps / (scalar_ms * 1e-3), 1.0});
    (void)sink;

    if (disp.is_available(hw::Backend::OPENMP)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.log_sum_exp(data, hw::Backend::OPENMP);
        double ms = t.elapsed_ms();
        results.push_back({"LogSumExp", "OpenMP", n, ms / reps,
                            n * reps / (ms * 1e-3), scalar_ms / ms});
    }

    if (disp.is_available(hw::Backend::METAL)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.log_sum_exp(data, hw::Backend::METAL);
        double ms = t.elapsed_ms();
        results.push_back({"LogSumExp", "Metal", n, ms / reps,
                            n * reps / (ms * 1e-3), scalar_ms / ms});
    }
}

static void bench_boltzmann(hw::UnifiedHardwareDispatch& disp, int n, int reps,
                             std::vector<BenchResult>& results) {
    auto data = random_doubles(n, -15.0, 5.0);
    double beta = 1.0 / (0.001987206 * 298.15);

    Timer t;
    t.start();
    for (int r = 0; r < reps; ++r) {
        auto w = disp.compute_boltzmann_weights(data, beta, hw::Backend::SCALAR);
    }
    double scalar_ms = t.elapsed_ms();
    results.push_back({"BoltzmannWeights", "scalar", n, scalar_ms / reps,
                        n * reps / (scalar_ms * 1e-3), 1.0});

    // AUTO (best available)
    t.start();
    for (int r = 0; r < reps; ++r) {
        auto w = disp.compute_boltzmann_weights(data, beta, hw::Backend::AUTO);
    }
    double ms = t.elapsed_ms();
    auto best = disp.best_backend(hw::KernelType::BOLTZMANN_WEIGHTS);
    results.push_back({"BoltzmannWeights",
                        hw::UnifiedHardwareDispatch::backend_name(best), n, ms / reps,
                        n * reps / (ms * 1e-3), scalar_ms / ms});

    if (disp.is_available(hw::Backend::METAL)) {
        t.start();
        for (int r = 0; r < reps; ++r) {
            auto w = disp.compute_boltzmann_weights(data, beta, hw::Backend::METAL);
        }
        ms = t.elapsed_ms();
        results.push_back({"BoltzmannWeights", "Metal", n, ms / reps,
                            n * reps / (ms * 1e-3), scalar_ms / ms});
    }
}

static void bench_rmsd(hw::UnifiedHardwareDispatch& disp, int n_atoms, int reps,
                        std::vector<BenchResult>& results) {
    auto a = random_coords(n_atoms, 42);
    auto b = random_coords(n_atoms, 77);

    Timer t;
    volatile float sink = 0;

    // Scalar
    t.start();
    for (int r = 0; r < reps; ++r)
        sink = disp.rmsd(a.data(), b.data(), n_atoms, hw::Backend::SCALAR);
    double scalar_ms = t.elapsed_ms();
    results.push_back({"RMSD", "scalar", n_atoms, scalar_ms / reps,
                        n_atoms * reps / (scalar_ms * 1e-3), 1.0});
    (void)sink;

    if (disp.is_available(hw::Backend::AVX2)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.rmsd(a.data(), b.data(), n_atoms, hw::Backend::AVX2);
        double ms = t.elapsed_ms();
        results.push_back({"RMSD", "AVX2", n_atoms, ms / reps,
                            n_atoms * reps / (ms * 1e-3), scalar_ms / ms});
    }

    if (disp.is_available(hw::Backend::AVX512)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.rmsd(a.data(), b.data(), n_atoms, hw::Backend::AVX512);
        double ms = t.elapsed_ms();
        results.push_back({"RMSD", "AVX-512", n_atoms, ms / reps,
                            n_atoms * reps / (ms * 1e-3), scalar_ms / ms});
    }

    if (disp.is_available(hw::Backend::OPENMP)) {
        t.start();
        for (int r = 0; r < reps; ++r)
            sink = disp.rmsd(a.data(), b.data(), n_atoms, hw::Backend::OPENMP);
        double ms = t.elapsed_ms();
        results.push_back({"RMSD", "OpenMP", n_atoms, ms / reps,
                            n_atoms * reps / (ms * 1e-3), scalar_ms / ms});
    }
}

static void bench_distance2_batch(hw::UnifiedHardwareDispatch& disp, int n, int reps,
                                   std::vector<BenchResult>& results) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-50.0f, 50.0f);
    std::vector<float> ax(n), ay(n), az(n), out(n);
    for (int i = 0; i < n; ++i) {
        ax[i] = dist(rng);
        ay[i] = dist(rng);
        az[i] = dist(rng);
    }
    float bx = 0.0f, by = 0.0f, bz = 0.0f;

    Timer t;

    // Scalar
    t.start();
    for (int r = 0; r < reps; ++r)
        disp.distance2_batch(ax.data(), ay.data(), az.data(), bx, by, bz,
                              out.data(), n, hw::Backend::SCALAR);
    double scalar_ms = t.elapsed_ms();
    results.push_back({"Distance2Batch", "scalar", n, scalar_ms / reps,
                        n * reps / (scalar_ms * 1e-3), 1.0});

    // AUTO
    t.start();
    for (int r = 0; r < reps; ++r)
        disp.distance2_batch(ax.data(), ay.data(), az.data(), bx, by, bz,
                              out.data(), n, hw::Backend::AUTO);
    double ms = t.elapsed_ms();
    auto best = disp.best_backend(hw::KernelType::DISTANCE_BATCH);
    results.push_back({"Distance2Batch",
                        hw::UnifiedHardwareDispatch::backend_name(best), n, ms / reps,
                        n * reps / (ms * 1e-3), scalar_ms / ms});
}

// ─── main ────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    int size = 100000;
    int reps = 100;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--size" && i + 1 < argc)
            size = std::atoi(argv[++i]);
        else if (arg == "--reps" && i + 1 < argc)
            reps = std::atoi(argv[++i]);
        else if (arg == "--help") {
            std::cerr << "Usage: benchmark_dispatch [--size N] [--reps R]\n";
            return 0;
        }
    }

    auto& disp = hw::UnifiedHardwareDispatch::instance();
    disp.detect();

    std::cout << disp.hardware_report() << "\n";
    std::cout << "Benchmark: size=" << size << ", reps=" << reps << "\n\n";

    std::vector<BenchResult> results;

    print_header();

    bench_shannon(disp, size, reps, results);
    for (auto& r : results) print_row(r);
    results.clear();

    bench_lse(disp, size, reps, results);
    for (auto& r : results) print_row(r);
    results.clear();

    bench_boltzmann(disp, size, reps, results);
    for (auto& r : results) print_row(r);
    results.clear();

    bench_rmsd(disp, size / 10, reps, results);
    for (auto& r : results) print_row(r);
    results.clear();

    bench_distance2_batch(disp, size, reps, results);
    for (auto& r : results) print_row(r);
    results.clear();

    std::cout << "\nDone.\n";
    return 0;
}
