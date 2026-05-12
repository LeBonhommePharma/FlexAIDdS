// UnifiedHardwareDispatch.h — Single entry point for all hardware dispatch
//
// Merges the former HardwareDispatch.h (hw:: namespace, Meyers singleton) and
// hardware_dispatch.h (flexaids:: namespace, free functions) into a unified
// API with DispatchResult error propagation.
//
// Backward-compatible: old headers forward-include this file with aliases.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab
#pragma once

#include "hardware_detect.h"
#include <string>
#include <vector>
#include <array>
#include <functional>
#include <chrono>
#include <cstdint>
#include <span>
#include <optional>

// Forward declarations for GPU contexts (opaque handles)
#ifdef FLEXAIDS_USE_CUDA
struct CudaEvalCtx;
#endif
#ifdef FLEXAIDS_USE_METAL
struct MetalEvalCtx;
#endif

namespace hw {

// ─── Backend enumeration (unified) ──────────────────────────────────────────

enum class Backend : uint8_t {
    SCALAR  = 0,
    OPENMP  = 1,
    AVX2    = 2,
    AVX512  = 3,
    METAL   = 4,
    CUDA    = 5,
    ROCM    = 6,
    AUTO    = 255
};

// ─── Kernel types ───────────────────────────────────────────────────────────

enum class KernelType : uint8_t {
    SHANNON_ENTROPY   = 0,
    FITNESS_EVAL      = 1,
    DISTANCE_BATCH    = 2,
    BOLTZMANN_WEIGHTS = 3,
    PARTITION_FUNC    = 4,
    CONTACT_DISC      = 5,
    HESSIAN_ASM       = 6,
    KNN_SEARCH        = 7,
    CAVITY_DET        = 8,
    TURBO_QUANT       = 9,
    RMSD              = 10,
};

// ─── Error handling ─────────────────────────────────────────────────────────

enum class DispatchError : uint8_t {
    OK = 0,
    NO_BACKEND,
    ALLOC_FAILED,
    LAUNCH_FAILED,
    SYNC_FAILED,
    INVALID_ARGS,
    ERR_OVERFLOW,
    DEVICE_LOST,
};

struct DispatchResult {
    DispatchError error        = DispatchError::OK;
    Backend       used_backend = Backend::AUTO;
    double        elapsed_ms   = 0.0;
    std::string   detail;

    explicit operator bool() const { return error == DispatchError::OK; }
};

// ─── Hardware capabilities ──────────────────────────────────────────────────

struct HardwareInfo {
    // CPU
    std::string cpu_name;
    int         logical_cores    = 1;
    bool        has_avx2         = false;
    bool        has_fma          = false;
    bool        has_avx512f      = false;
    bool        has_avx512dq     = false;
    bool        has_avx512bw     = false;
    bool        has_openmp       = false;
    int         omp_max_threads  = 1;

    // GPU
    bool        has_cuda         = false;
    std::string cuda_device_name;
    int         cuda_sm_count    = 0;
    std::size_t cuda_global_mem  = 0;

    bool        has_rocm         = false;
    std::string rocm_device_name;

    bool        has_metal        = false;
    std::string metal_device_name;

    // Eigen
    bool        has_eigen        = false;
};

// ─── Per-call telemetry ─────────────────────────────────────────────────────

struct DispatchTelemetry {
    Backend backend;
    double  wall_time_ms;
    int64_t elements;
    double  throughput_meps;

    std::string summary() const;
};

// ─── Boltzmann batch result ─────────────────────────────────────────────────

struct BoltzmannBatchResult {
    std::vector<double> weights;
    double              log_Z;
    double              E_min;
    DispatchTelemetry   telemetry;
};

// ─── Benchmark result ───────────────────────────────────────────────────────

struct BenchmarkResult {
    Backend     backend;
    KernelType  kernel;
    double      time_seconds;
    double      speedup;
    double      throughput;
    std::string backend_name;
};

// ─── Dispatch report ────────────────────────────────────────────────────────

struct DispatchReport {
    Backend     selected;
    std::string reason;
    std::string hw_summary;
};

// ─── UnifiedHardwareDispatch (Meyers singleton) ─────────────────────────────

class UnifiedHardwareDispatch {
public:
    static UnifiedHardwareDispatch& instance();

    // Detection (idempotent, thread-safe)
    void detect();
    const HardwareInfo& info() const noexcept { return info_; }

    // Backend selection
    Backend best_backend(KernelType kernel = KernelType::SHANNON_ENTROPY) const;
    void set_override(Backend b) noexcept { override_ = b; }
    void clear_override() noexcept { override_ = Backend::AUTO; }
    Backend current_override() const noexcept { return override_; }

    bool is_available(Backend b) const noexcept;
    static const char* backend_name(Backend b) noexcept;
    std::vector<Backend> available_backends() const;
    std::string hardware_report() const;

    // ─── Dispatched compute functions (from former HardwareDispatch.h) ───

    // Returns histogram Shannon entropy in nats. Convert to bits at reporting
    // boundaries with H_bits = H_nats / ln(2).
    double compute_shannon_entropy(const std::vector<double>& values,
                                   int num_bins = 20,
                                   Backend backend = Backend::AUTO);

    std::vector<double> compute_boltzmann_weights(
        const std::vector<double>& energies,
        double beta,
        Backend backend = Backend::AUTO);

    double log_sum_exp(const std::vector<double>& values,
                       Backend backend = Backend::AUTO);

    void distance2_batch(const float* ax, const float* ay, const float* az,
                         float bx, float by, float bz,
                         float* out, int n,
                         Backend backend = Backend::AUTO);

    float rmsd(const float* a_xyz, const float* b_xyz, int n_atoms,
               Backend backend = Backend::AUTO);

    // ─── Dispatched compute with DispatchResult (from former hardware_dispatch.h) ─

    BoltzmannBatchResult compute_boltzmann_batch(
        std::span<const double> energies,
        double beta);

    double log_sum_exp_dispatch(std::span<const double> values);

    DispatchReport get_dispatch_report() const;

    Backend select_cpu_backend() const;

private:
    UnifiedHardwareDispatch() = default;

    HardwareInfo info_;
    bool         detected_  = false;
    Backend      override_  = Backend::AUTO;

    // Internal detection helpers
    void detect_cpu();
    void detect_gpu();
    void detect_libraries();

    // Internal Shannon implementations
    double shannon_scalar(const std::vector<double>& values, int num_bins);
    double shannon_openmp(const std::vector<double>& values, int num_bins);
    double shannon_avx512(const std::vector<double>& values, int num_bins);
    double shannon_avx512_omp(const std::vector<double>& values, int num_bins);

    // Internal log-sum-exp implementations
    double lse_scalar(const std::vector<double>& values);
    double lse_openmp(const std::vector<double>& values);
    double lse_avx512(const std::vector<double>& values);
    double lse_eigen(const std::vector<double>& values);

    // Internal RMSD implementations
    float rmsd_scalar(const float* a, const float* b, int n);
    float rmsd_avx2(const float* a, const float* b, int n);
    float rmsd_avx512(const float* a, const float* b, int n);
    float rmsd_openmp(const float* a, const float* b, int n);
};

}  // namespace hw

// ─── flexaids:: compatibility layer ─────────────────────────────────────────
// Preserves the free-function API from the former hardware_dispatch.h so that
// existing callers (gaboom.cpp, statmech.cpp) continue to compile unchanged.

namespace flexaids {

// Map old HardwareBackend enum values to hw::Backend
enum class HardwareBackend : uint8_t {
    CUDA    = static_cast<uint8_t>(hw::Backend::CUDA),
    ROCM    = static_cast<uint8_t>(hw::Backend::ROCM),
    METAL   = static_cast<uint8_t>(hw::Backend::METAL),
    AVX512  = static_cast<uint8_t>(hw::Backend::AVX512),
    AVX2    = static_cast<uint8_t>(hw::Backend::AVX2),
    OPENMP  = static_cast<uint8_t>(hw::Backend::OPENMP),
    SCALAR  = static_cast<uint8_t>(hw::Backend::SCALAR),
};

const char* backend_name(HardwareBackend b) noexcept;
HardwareBackend select_backend();
HardwareBackend select_cpu_backend();

// Re-export types
using DispatchTelemetry = hw::DispatchTelemetry;
using BoltzmannBatchResult = hw::BoltzmannBatchResult;
using DispatchReport = hw::DispatchReport;

BoltzmannBatchResult compute_boltzmann_batch(
    std::span<const double> energies,
    double beta);

double log_sum_exp_dispatch(std::span<const double> values);

DispatchReport get_dispatch_report();

}  // namespace flexaids

// Legacy alias for backward compatibility
using HardwareDispatcher = hw::UnifiedHardwareDispatch;
