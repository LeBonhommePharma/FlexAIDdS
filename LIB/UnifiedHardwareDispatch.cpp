// UnifiedHardwareDispatch.cpp — Merged implementation of all hardware dispatch
//
// Combines the former HardwareDispatch.cpp (hw:: singleton, Shannon/LSE/RMSD)
// and hardware_dispatch.cpp (flexaids:: free functions, Boltzmann batch) into
// a single translation unit.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#include "UnifiedHardwareDispatch.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <sstream>
#include <thread>
#include <limits>

#ifdef __x86_64__
#  include <cpuid.h>
#endif

#ifdef __AVX2__
#  include <immintrin.h>
#endif

#ifdef __AVX512F__
#  include <immintrin.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <Eigen/Dense>

#if defined(__AVX512F__) && defined(__AVX512DQ__)
#  define HAS_AVX512_RT 1
#else
#  define HAS_AVX512_RT 0
#endif

#ifdef FLEXAIDS_HAS_METAL_SHANNON
#  include "ShannonThermoStack/ShannonMetalBridge.h"
#endif

#ifdef FLEXAIDS_USE_METAL
#  include "metal_eval.h"
#endif

namespace hw {

// ═════════════════════════════════════════════════════════════════════════════
// Singleton
// ═════════════════════════════════════════════════════════════════════════════

UnifiedHardwareDispatch& UnifiedHardwareDispatch::instance() {
    static UnifiedHardwareDispatch inst;
    return inst;
}

// ═════════════════════════════════════════════════════════════════════════════
// Detection
// ═════════════════════════════════════════════════════════════════════════════

void UnifiedHardwareDispatch::detect() {
    if (detected_) return;
    detect_cpu();
    detect_gpu();
    detect_libraries();
    detected_ = true;
}

void UnifiedHardwareDispatch::detect_cpu() {
    info_.logical_cores = static_cast<int>(std::thread::hardware_concurrency());
    if (info_.logical_cores < 1) info_.logical_cores = 1;

#ifdef __x86_64__
    {
        unsigned int eax, ebx, ecx, edx;
        char brand[49] = {};

        for (unsigned int leaf = 0x80000002; leaf <= 0x80000004; ++leaf) {
            if (__get_cpuid(leaf, &eax, &ebx, &ecx, &edx)) {
                int offset = static_cast<int>((leaf - 0x80000002) * 16);
                std::memcpy(brand + offset,      &eax, 4);
                std::memcpy(brand + offset + 4,   &ebx, 4);
                std::memcpy(brand + offset + 8,   &ecx, 4);
                std::memcpy(brand + offset + 12,  &edx, 4);
            }
        }
        brand[48] = '\0';
        info_.cpu_name = brand;
        auto pos = info_.cpu_name.find_first_not_of(' ');
        if (pos != std::string::npos) info_.cpu_name = info_.cpu_name.substr(pos);

        if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
            info_.has_avx2     = (ebx & (1u << 5))  != 0;
            info_.has_avx512f  = (ebx & (1u << 16)) != 0;
            info_.has_avx512dq = (ebx & (1u << 17)) != 0;
            info_.has_avx512bw = (ebx & (1u << 30)) != 0;
        }

        if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
            info_.has_fma = (ecx & (1u << 12)) != 0;
        }
    }
#elif defined(__aarch64__)
    info_.cpu_name = "AArch64";
#endif

#ifdef _OPENMP
    info_.has_openmp      = true;
    info_.omp_max_threads = omp_get_max_threads();
#else
    info_.has_openmp      = false;
    info_.omp_max_threads = 1;
#endif
}

void UnifiedHardwareDispatch::detect_gpu() {
#ifdef FLEXAIDS_USE_CUDA
    info_.has_cuda = true;
    info_.cuda_device_name = "CUDA device (compiled-in)";
#else
    info_.has_cuda = false;
#endif

#ifdef FLEXAIDS_USE_METAL
    MetalCapabilities caps{};
    metal_eval_get_capabilities(&caps);

    info_.has_metal = caps.available;
    if (caps.available) {
        info_.metal_device_name = caps.device_name;
        // Optionally store more in the future (unified mem, cores, etc.)
    } else {
        info_.metal_device_name = "";
    }
#else
    info_.has_metal = false;
#endif

    // Multi-cleft support: allow forcing CPU path to avoid Metal kernel capacity limits (e.g. max_pop 1000)
    // when running many independent per-cleft GAs. Non-behavioral for normal use.
    if (const char* dis = std::getenv("FLEXAIDDS_DISABLE_METAL")) {
        if (*dis) {
            info_.has_metal = false;
            info_.metal_device_name = "disabled-by-env";
        }
    }

#ifdef FLEXAIDS_USE_ROCM
    info_.has_rocm = true;
    info_.rocm_device_name = "ROCm device (compiled-in)";
#else
    info_.has_rocm = false;
#endif
}

void UnifiedHardwareDispatch::detect_libraries() {
    info_.has_eigen = true;
}

// ═════════════════════════════════════════════════════════════════════════════
// Backend queries
// ═════════════════════════════════════════════════════════════════════════════

bool UnifiedHardwareDispatch::is_available(Backend b) const noexcept {
    switch (b) {
        case Backend::SCALAR:  return true;
        case Backend::OPENMP:  return info_.has_openmp;
        case Backend::AVX2:
#ifdef __AVX2__
            return info_.has_avx2 && info_.has_fma;
#else
            return false;
#endif
        case Backend::AVX512:
#ifdef __AVX512F__
            return info_.has_avx512f;
#else
            return false;
#endif
        case Backend::METAL:   return info_.has_metal;
        case Backend::CUDA:    return info_.has_cuda;
        case Backend::ROCM:    return info_.has_rocm;
        case Backend::AUTO:    return true;
    }
    return false;
}

Backend UnifiedHardwareDispatch::best_backend(KernelType kernel) const {
    if (override_ != Backend::AUTO) {
        return override_;
    }

    switch (kernel) {
        case KernelType::SHANNON_ENTROPY:
        case KernelType::FITNESS_EVAL:
        case KernelType::CONTACT_DISC:
        case KernelType::HESSIAN_ASM:
        case KernelType::TURBO_QUANT:
            if (is_available(Backend::CUDA))   return Backend::CUDA;
            if (is_available(Backend::ROCM))   return Backend::ROCM;
            if (is_available(Backend::METAL)) {
                // P1.6 loud dispatch log for Metal on M3 (Shannon for configurational entropy part of thermo ledger / best BindingMode)
                fprintf(stderr, "[HW-DISPATCH] Kernel SHANNON_ENTROPY/FITNESS etc. -> using Metal on Apple (M3 Pro) for entropy/thermo (see .metallib + UnifiedHardwareDispatch)\n");
                return Backend::METAL;
            }
            if (is_available(Backend::AVX512)) return Backend::AVX512;
            if (is_available(Backend::OPENMP)) return Backend::OPENMP;
            return Backend::SCALAR;

        case KernelType::DISTANCE_BATCH:
        case KernelType::RMSD:
            if (is_available(Backend::AVX512)) return Backend::AVX512;
            if (is_available(Backend::AVX2))   return Backend::AVX2;
            if (is_available(Backend::OPENMP)) return Backend::OPENMP;
            return Backend::SCALAR;

        case KernelType::BOLTZMANN_WEIGHTS:
        case KernelType::PARTITION_FUNC:
            if (is_available(Backend::AVX512)) return Backend::AVX512;
            if (is_available(Backend::AVX2))   return Backend::AVX2;
            if (is_available(Backend::OPENMP)) return Backend::OPENMP;
            return Backend::SCALAR;

        case KernelType::KNN_SEARCH:
            if (is_available(Backend::CUDA))   return Backend::CUDA;
            if (is_available(Backend::ROCM))   return Backend::ROCM;
            return Backend::SCALAR;

        case KernelType::CAVITY_DET:
            if (is_available(Backend::METAL))  return Backend::METAL;
            return Backend::SCALAR;
    }
    return Backend::SCALAR;
}

Backend UnifiedHardwareDispatch::select_cpu_backend() const {
#if HAS_AVX512_RT
    if (is_available(Backend::AVX512))
        return Backend::AVX512;
#endif
    if (is_available(Backend::AVX2))
        return Backend::AVX2;
    if (is_available(Backend::OPENMP))
        return Backend::OPENMP;
    return Backend::SCALAR;
}

const char* UnifiedHardwareDispatch::backend_name(Backend b) noexcept {
    switch (b) {
        case Backend::SCALAR:  return "scalar";
        case Backend::OPENMP:  return "OpenMP";
        case Backend::AVX2:    return "AVX2";
        case Backend::AVX512:  return "AVX-512";
        case Backend::METAL:   return "Metal";
        case Backend::CUDA:    return "CUDA";
        case Backend::ROCM:    return "ROCm";
        case Backend::AUTO:    return "auto";
    }
    return "unknown";
}

std::vector<Backend> UnifiedHardwareDispatch::available_backends() const {
    std::vector<Backend> result;
    if (is_available(Backend::CUDA))   result.push_back(Backend::CUDA);
    if (is_available(Backend::ROCM))   result.push_back(Backend::ROCM);
    if (is_available(Backend::METAL))  result.push_back(Backend::METAL);
    if (is_available(Backend::AVX512)) result.push_back(Backend::AVX512);
    if (is_available(Backend::AVX2))   result.push_back(Backend::AVX2);
    if (is_available(Backend::OPENMP)) result.push_back(Backend::OPENMP);
    result.push_back(Backend::SCALAR);
    return result;
}

std::string UnifiedHardwareDispatch::hardware_report() const {
    std::ostringstream os;
    os << "=== FlexAIDdS Hardware Report ===\n";
    os << "CPU: " << info_.cpu_name << "\n";
    os << "Cores: " << info_.logical_cores << "\n";
    os << "AVX2:     " << (info_.has_avx2    ? "YES" : "no") << "\n";
    os << "FMA:      " << (info_.has_fma     ? "YES" : "no") << "\n";
    os << "AVX-512F: " << (info_.has_avx512f ? "YES" : "no") << "\n";
    os << "AVX-512DQ:" << (info_.has_avx512dq? "YES" : "no") << "\n";
    os << "AVX-512BW:" << (info_.has_avx512bw? "YES" : "no") << "\n";
    os << "OpenMP:   " << (info_.has_openmp  ? "YES" : "no");
    if (info_.has_openmp) os << " (" << info_.omp_max_threads << " threads)";
    os << "\n";
    os << "Eigen3:   " << (info_.has_eigen   ? "YES" : "no") << "\n";
    os << "CUDA:     " << (info_.has_cuda    ? "YES" : "no");
    if (info_.has_cuda) os << " (" << info_.cuda_device_name << ")";
    os << "\n";
    os << "ROCm:     " << (info_.has_rocm    ? "YES" : "no");
    if (info_.has_rocm) os << " (" << info_.rocm_device_name << ")";
    os << "\n";
    os << "Metal:    " << (info_.has_metal   ? "YES" : "no");
    if (info_.has_metal) os << " (" << info_.metal_device_name << ")";
    os << "\n";
    os << "Available backends: ";
    auto avail = available_backends();
    for (size_t i = 0; i < avail.size(); ++i) {
        if (i > 0) os << ", ";
        os << backend_name(avail[i]);
    }
    os << "\n";
    os << "Best entropy backend:  " << backend_name(best_backend(KernelType::SHANNON_ENTROPY)) << "\n";
    os << "Best distance backend: " << backend_name(best_backend(KernelType::DISTANCE_BATCH)) << "\n";
    return os.str();
}

// ═════════════════════════════════════════════════════════════════════════════
// Telemetry helper
// ═════════════════════════════════════════════════════════════════════════════

std::string DispatchTelemetry::summary() const {
    std::ostringstream ss;
    ss << UnifiedHardwareDispatch::backend_name(backend) << ": "
       << elements << " elements in "
       << wall_time_ms << " ms ("
       << throughput_meps << " M elem/s)";
    return ss.str();
}

static DispatchTelemetry make_telemetry(
    Backend backend,
    std::chrono::steady_clock::time_point start,
    int64_t elements)
{
    auto end = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();
    double meps = ms > 0.0 ? (elements / 1e6) / (ms / 1000.0) : 0.0;
    return { backend, ms, elements, meps };
}

// ═════════════════════════════════════════════════════════════════════════════
// Shannon entropy dispatch
// ═════════════════════════════════════════════════════════════════════════════

static double entropy_from_hist(const int* counts, int num_bins, int total) {
    if (total == 0) return 0.0;

    Eigen::ArrayXd prob(num_bins);
    for (int b = 0; b < num_bins; ++b)
        prob(b) = static_cast<double>(counts[b]);
    prob /= static_cast<double>(total);
    Eigen::ArrayXd safe_p = (prob > 1e-15).select(prob, Eigen::ArrayXd::Constant(num_bins, 1.0));
    Eigen::ArrayXd lp     = (prob > 1e-15).select(safe_p.log(), Eigen::ArrayXd::Zero(num_bins));
    return -(prob * lp).sum();
}

double UnifiedHardwareDispatch::shannon_scalar(const std::vector<double>& values, int num_bins) {
    double min_v = *std::min_element(values.begin(), values.end());
    double max_v = *std::max_element(values.begin(), values.end());
    if (max_v - min_v < 1e-12) return 0.0;
    double bin_width = (max_v - min_v) / num_bins + 1e-10;
    double inv_bw    = 1.0 / bin_width;
    int    n         = static_cast<int>(values.size());

    std::vector<int> bins(num_bins, 0);
    for (int i = 0; i < n; ++i) {
        int b = static_cast<int>((values[i] - min_v) * inv_bw);
        bins[std::min(std::max(b, 0), num_bins - 1)]++;
    }
    return entropy_from_hist(bins.data(), num_bins, n);
}

double UnifiedHardwareDispatch::shannon_openmp(const std::vector<double>& values, int num_bins) {
#ifdef _OPENMP
    double min_v = *std::min_element(values.begin(), values.end());
    double max_v = *std::max_element(values.begin(), values.end());
    if (max_v - min_v < 1e-12) return 0.0;
    double bin_width = (max_v - min_v) / num_bins + 1e-10;
    double inv_bw    = 1.0 / bin_width;
    int    n         = static_cast<int>(values.size());

    int n_threads = omp_get_max_threads();
    std::vector<std::vector<int>> t_bins(n_threads, std::vector<int>(num_bins, 0));
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        int tid = omp_get_thread_num();
        int b   = static_cast<int>((values[i] - min_v) * inv_bw);
        t_bins[tid][std::min(std::max(b, 0), num_bins - 1)]++;
    }
    std::vector<int> bins(num_bins, 0);
    for (auto& tb : t_bins)
        for (int b = 0; b < num_bins; ++b) bins[b] += tb[b];
    return entropy_from_hist(bins.data(), num_bins, n);
#else
    return shannon_scalar(values, num_bins);
#endif
}

double UnifiedHardwareDispatch::shannon_avx512(const std::vector<double>& values, int num_bins) {
#ifdef __AVX512F__
    double min_v = *std::min_element(values.begin(), values.end());
    double max_v = *std::max_element(values.begin(), values.end());
    if (max_v - min_v < 1e-12) return 0.0;
    double bin_width = (max_v - min_v) / num_bins + 1e-10;
    double inv_bw    = 1.0 / bin_width;
    int    n         = static_cast<int>(values.size());

    std::vector<int> bins(num_bins, 0);

    __m512d vmin   = _mm512_set1_pd(min_v);
    __m512d vinvbw = _mm512_set1_pd(inv_bw);

    int i = 0;
    for (; i + 7 < n; i += 8) {
        __m512d ve   = _mm512_loadu_pd(values.data() + i);
        __m512d vrel = _mm512_mul_pd(_mm512_sub_pd(ve, vmin), vinvbw);
        __m256i v32  = _mm512_cvttpd_epi32(vrel);
        alignas(32) int tmp[8];
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(tmp), v32);
        for (int k = 0; k < 8; ++k) {
            int b = std::min(std::max(tmp[k], 0), num_bins - 1);
            bins[b]++;
        }
    }
    for (; i < n; ++i) {
        int b = static_cast<int>((values[i] - min_v) * inv_bw);
        bins[std::min(std::max(b, 0), num_bins - 1)]++;
    }
    return entropy_from_hist(bins.data(), num_bins, n);
#else
    return shannon_scalar(values, num_bins);
#endif
}

double UnifiedHardwareDispatch::shannon_avx512_omp(const std::vector<double>& values, int num_bins) {
#if defined(__AVX512F__) && defined(_OPENMP)
    double min_v = *std::min_element(values.begin(), values.end());
    double max_v = *std::max_element(values.begin(), values.end());
    if (max_v - min_v < 1e-12) return 0.0;
    double bin_width = (max_v - min_v) / num_bins + 1e-10;
    double inv_bw    = 1.0 / bin_width;
    int    n         = static_cast<int>(values.size());

    int n_threads = omp_get_max_threads();
    std::vector<std::vector<int>> t_bins(n_threads, std::vector<int>(num_bins, 0));

    #pragma omp parallel
    {
        int tid   = omp_get_thread_num();
        int chunk = (n + n_threads - 1) / n_threads;
        int start = tid * chunk;
        int end   = std::min(start + chunk, n);

        __m512d vmin   = _mm512_set1_pd(min_v);
        __m512d vinvbw = _mm512_set1_pd(inv_bw);

        int i = start;
        for (; i + 7 < end; i += 8) {
            __m512d ve   = _mm512_loadu_pd(values.data() + i);
            __m512d vrel = _mm512_mul_pd(_mm512_sub_pd(ve, vmin), vinvbw);
            __m256i v32  = _mm512_cvttpd_epi32(vrel);
            alignas(32) int tmp[8];
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(tmp), v32);
            for (int k = 0; k < 8; ++k) {
                int b = std::min(std::max(tmp[k], 0), num_bins - 1);
                t_bins[tid][b]++;
            }
        }
        for (; i < end; ++i) {
            int b = static_cast<int>((values[i] - min_v) * inv_bw);
            t_bins[tid][std::min(std::max(b, 0), num_bins - 1)]++;
        }
    }
    std::vector<int> bins(num_bins, 0);
    for (auto& tb : t_bins)
        for (int b = 0; b < num_bins; ++b) bins[b] += tb[b];
    return entropy_from_hist(bins.data(), num_bins, n);
#else
    return shannon_avx512(values, num_bins);
#endif
}

double UnifiedHardwareDispatch::compute_shannon_entropy(
    const std::vector<double>& values, int num_bins, Backend backend)
{
    if (values.empty()) return 0.0;
    if (num_bins <= 0)  num_bins = 20;
    if (!detected_)     detect();

    Backend b = (backend == Backend::AUTO) ? best_backend(KernelType::SHANNON_ENTROPY) : backend;

    switch (b) {
        case Backend::CUDA:
        case Backend::ROCM:
        case Backend::METAL:
            if (is_available(Backend::AVX512) && is_available(Backend::OPENMP))
                return shannon_avx512_omp(values, num_bins);
            if (is_available(Backend::AVX512))
                return shannon_avx512(values, num_bins);
            if (is_available(Backend::OPENMP))
                return shannon_openmp(values, num_bins);
            return shannon_scalar(values, num_bins);

        case Backend::AVX512:
            if (is_available(Backend::OPENMP))
                return shannon_avx512_omp(values, num_bins);
            return shannon_avx512(values, num_bins);

        case Backend::AVX2:
        case Backend::OPENMP:
            return shannon_openmp(values, num_bins);

        case Backend::SCALAR:
            return shannon_scalar(values, num_bins);

        default:
            return shannon_scalar(values, num_bins);
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Log-sum-exp dispatch
// ═════════════════════════════════════════════════════════════════════════════

double UnifiedHardwareDispatch::lse_scalar(const std::vector<double>& values) {
    double x_max = *std::max_element(values.begin(), values.end());
    if (x_max <= -1e308) return x_max;
    double sum = 0.0;
    for (double v : values)
        sum += std::exp(v - x_max);
    return x_max + std::log(sum);
}

double UnifiedHardwareDispatch::lse_openmp(const std::vector<double>& values) {
#ifdef _OPENMP
    double x_max = *std::max_element(values.begin(), values.end());
    if (x_max <= -1e308) return x_max;
    double sum = 0.0;
    int n = static_cast<int>(values.size());
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < n; ++i)
        sum += std::exp(values[i] - x_max);
    return x_max + std::log(sum);
#else
    return lse_scalar(values);
#endif
}

double UnifiedHardwareDispatch::lse_avx512(const std::vector<double>& values) {
#ifdef __AVX512F__
    int n = static_cast<int>(values.size());
    double x_max = *std::max_element(values.begin(), values.end());
    if (x_max <= -1e308) return x_max;

    double sum = 0.0;
    for (int i = 0; i < n; ++i)
        sum += std::exp(values[i] - x_max);
    return x_max + std::log(sum);
#else
    return lse_scalar(values);
#endif
}

double UnifiedHardwareDispatch::lse_eigen(const std::vector<double>& values) {
    Eigen::Map<const Eigen::ArrayXd> arr(values.data(),
                                          static_cast<Eigen::Index>(values.size()));
    double x_max = arr.maxCoeff();
    if (x_max <= -1e308) return x_max;
    double sum = (arr - x_max).exp().sum();
    return x_max + std::log(sum);
}

double UnifiedHardwareDispatch::log_sum_exp(const std::vector<double>& values, Backend backend) {
    if (values.empty()) return -1e308;
    if (!detected_) detect();

    Backend b = (backend == Backend::AUTO) ? best_backend(KernelType::PARTITION_FUNC) : backend;

    switch (b) {
        case Backend::AVX512: return lse_avx512(values);
        case Backend::OPENMP: return lse_openmp(values);
        case Backend::SCALAR: return lse_scalar(values);
        default:
            // For AUTO and GPU backends, prefer Eigen if available
            if (info_.has_eigen) return lse_eigen(values);
            return lse_scalar(values);
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Boltzmann weights (vector API)
// ═════════════════════════════════════════════════════════════════════════════

std::vector<double> UnifiedHardwareDispatch::compute_boltzmann_weights(
    const std::vector<double>& energies, double beta, Backend backend)
{
    if (energies.empty()) return {};
    if (!detected_) detect();

    const std::size_t N = energies.size();
    std::vector<double> log_w(N);
    for (std::size_t i = 0; i < N; ++i)
        log_w[i] = -beta * energies[i];

    double lnZ = log_sum_exp(log_w, backend);

    std::vector<double> w(N);

    if (N >= 16) {
        Eigen::Map<const Eigen::ArrayXd> lw(log_w.data(), static_cast<Eigen::Index>(N));
        Eigen::Map<Eigen::ArrayXd> out(w.data(), static_cast<Eigen::Index>(N));
        out = (lw - lnZ).exp();
        return w;
    }

#ifdef _OPENMP
    if (is_available(Backend::OPENMP) && N >= 4096) {
        int n = static_cast<int>(N);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i)
            w[i] = std::exp(log_w[i] - lnZ);
        return w;
    }
#endif

    for (std::size_t i = 0; i < N; ++i)
        w[i] = std::exp(log_w[i] - lnZ);
    return w;
}

// ═════════════════════════════════════════════════════════════════════════════
// Boltzmann batch (span API with telemetry — from hardware_dispatch.cpp)
// ═════════════════════════════════════════════════════════════════════════════

#if 0  // currently unused — scalar boltzmann is inlined in compute_boltzmann_batchstatic BoltzmannBatchResult boltzmann_scalar(
    std::span<const double> energies, double beta)
{
    auto start = std::chrono::steady_clock::now();
    const int n = static_cast<int>(energies.size());
    double E_min = *std::min_element(energies.begin(), energies.end());

    std::vector<double> weights(n);
    for (int i = 0; i < n; ++i)
        weights[i] = std::exp(-beta * (energies[i] - E_min));

    double sum_w = std::accumulate(weights.begin(), weights.end(), 0.0);
    double log_Z = std::log(sum_w) - beta * E_min;

    auto telemetry = make_telemetry(Backend::SCALAR, start, n);
    return { std::move(weights), log_Z, E_min, telemetry };
}
#endif  // boltzmann_scalar

#ifdef _OPENMP
static BoltzmannBatchResult boltzmann_openmp(
    std::span<const double> energies, double beta)
{
    auto start = std::chrono::steady_clock::now();
    const int n = static_cast<int>(energies.size());

    double E_min = std::numeric_limits<double>::max();
    #pragma omp parallel for reduction(min:E_min) schedule(static)
    for (int i = 0; i < n; ++i)
        E_min = std::min(E_min, energies[i]);

    std::vector<double> weights(n);
    double sum_w = 0.0;

    #pragma omp parallel for reduction(+:sum_w) schedule(static)
    for (int i = 0; i < n; ++i) {
        weights[i] = std::exp(-beta * (energies[i] - E_min));
        sum_w += weights[i];
    }

    double log_Z = std::log(sum_w) - beta * E_min;
    auto telemetry = make_telemetry(Backend::OPENMP, start, n);
    return { std::move(weights), log_Z, E_min, telemetry };
}
#endif

#if HAS_AVX512_RT
static BoltzmannBatchResult boltzmann_avx512(
    std::span<const double> energies, double beta)
{
    auto start = std::chrono::steady_clock::now();
    const int n = static_cast<int>(energies.size());

    double E_min = std::numeric_limits<double>::max();
    {
        __m512d vmin = _mm512_set1_pd(E_min);
        int i = 0;
        for (; i + 7 < n; i += 8) {
            __m512d ve = _mm512_loadu_pd(energies.data() + i);
            vmin = _mm512_min_pd(vmin, ve);
        }
        E_min = _mm512_reduce_min_pd(vmin);
        for (; i < n; ++i)
            E_min = std::min(E_min, energies[i]);
    }

    std::vector<double> weights(n);
    __m512d vneg_beta = _mm512_set1_pd(-beta);
    __m512d vEmin     = _mm512_set1_pd(E_min);
    double sum_w = 0.0;

    int i = 0;

#ifdef _OPENMP
    #pragma omp parallel reduction(+:sum_w)
    {
        int tid = omp_get_thread_num();
        int nt  = omp_get_num_threads();
        int chunk = (n + nt - 1) / nt;
        int start_idx = tid * chunk;
        int end   = std::min(start_idx + chunk, n);

        int j = start_idx;
        for (; j + 7 < end; j += 8) {
            __m512d ve  = _mm512_loadu_pd(energies.data() + j);
            __m512d arg = _mm512_mul_pd(_mm512_sub_pd(ve, vEmin), vneg_beta);

            alignas(64) double tmp[8];
            _mm512_storeu_pd(tmp, arg);
            for (int k = 0; k < 8; ++k) tmp[k] = std::exp(tmp[k]);
            __m512d vw = _mm512_loadu_pd(tmp);
            _mm512_storeu_pd(weights.data() + j, vw);
            sum_w += _mm512_reduce_add_pd(vw);
        }
        for (; j < end; ++j) {
            weights[j] = std::exp(-beta * (energies[j] - E_min));
            sum_w += weights[j];
        }
    }
#else
    for (; i + 7 < n; i += 8) {
        __m512d ve  = _mm512_loadu_pd(energies.data() + i);
        __m512d arg = _mm512_mul_pd(_mm512_sub_pd(ve, vEmin), vneg_beta);

        alignas(64) double tmp[8];
        _mm512_storeu_pd(tmp, arg);
        for (int k = 0; k < 8; ++k) tmp[k] = std::exp(tmp[k]);
        __m512d vw = _mm512_loadu_pd(tmp);
        _mm512_storeu_pd(weights.data() + i, vw);
        sum_w += _mm512_reduce_add_pd(vw);
    }
    for (; i < n; ++i) {
        weights[i] = std::exp(-beta * (energies[i] - E_min));
        sum_w += weights[i];
    }
#endif

    double log_Z = std::log(sum_w) - beta * E_min;
    auto telemetry = make_telemetry(Backend::AVX512, start, n);
    return { std::move(weights), log_Z, E_min, telemetry };
}
#endif

static BoltzmannBatchResult boltzmann_eigen(
    std::span<const double> energies, double beta, Backend cpu_backend)
{
    auto start = std::chrono::steady_clock::now();
    const int n = static_cast<int>(energies.size());
    Eigen::Map<const Eigen::ArrayXd> E(energies.data(), n);

    double E_min = E.minCoeff();
    Eigen::ArrayXd w = (-beta * (E - E_min)).exp();
    double sum_w = w.sum();
    double log_Z = std::log(sum_w) - beta * E_min;

    std::vector<double> weights(n);
    Eigen::Map<Eigen::ArrayXd>(weights.data(), n) = w;
    auto telemetry = make_telemetry(cpu_backend, start, n);
    return { std::move(weights), log_Z, E_min, telemetry };
}

#ifdef FLEXAIDS_HAS_METAL_SHANNON
static BoltzmannBatchResult boltzmann_metal(
    std::span<const double> energies, double beta)
{
    auto start = std::chrono::steady_clock::now();
    const int n = static_cast<int>(energies.size());

    std::vector<double> energy_vec(energies.begin(), energies.end());
    double sum_w = 0.0;
    double E_min = 0.0;

    auto weights = ShannonMetalBridge::compute_boltzmann_weights_metal(
        energy_vec, beta, sum_w, E_min);

    double log_Z = std::log(sum_w) - beta * E_min;
    auto telemetry = make_telemetry(Backend::METAL, start, n);
    return { std::move(weights), log_Z, E_min, telemetry };
}
#endif

BoltzmannBatchResult UnifiedHardwareDispatch::compute_boltzmann_batch(
    std::span<const double> energies,
    double beta)
{
    if (energies.empty())
        return { {}, 0.0, 0.0, { Backend::SCALAR, 0.0, 0, 0.0 } };

    if (!detected_) detect();

    [[maybe_unused]] Backend best = best_backend(KernelType::BOLTZMANN_WEIGHTS);
#ifdef FLEXAIDS_HAS_METAL_SHANNON
    if (is_available(Backend::METAL) && energies.size() >= 256) {
        if (ShannonMetalBridge::is_metal_available())
            return boltzmann_metal(energies, beta);
    }
#endif

    Backend cpu = select_cpu_backend();

#if HAS_AVX512_RT
    if (cpu == Backend::AVX512)
        return boltzmann_avx512(energies, beta);
#endif

    return boltzmann_eigen(energies, beta, cpu);
}

// ═════════════════════════════════════════════════════════════════════════════
// Log-sum-exp dispatch (span API — from hardware_dispatch.cpp)
// ═════════════════════════════════════════════════════════════════════════════

double UnifiedHardwareDispatch::log_sum_exp_dispatch(std::span<const double> values) {
    if (values.empty())
        return -std::numeric_limits<double>::infinity();

    if (!detected_) detect();

    const int n = static_cast<int>(values.size());

    double x_max = *std::max_element(values.begin(), values.end());
    if (!std::isfinite(x_max))
        return x_max;

#ifdef FLEXAIDS_HAS_METAL_SHANNON
    if (n >= 1024 && ShannonMetalBridge::is_metal_available()) {
        std::vector<double> vals(values.begin(), values.end());
        return ShannonMetalBridge::log_sum_exp_metal(vals);
    }
#endif

#if HAS_AVX512_RT
    {
        double sum = 0.0;
        int i = 0;
        for (; i + 7 < n; i += 8) {
            __m512d vx    = _mm512_loadu_pd(values.data() + i);
            __m512d vdiff = _mm512_sub_pd(vx, _mm512_set1_pd(x_max));

            alignas(64) double tmp[8];
            _mm512_storeu_pd(tmp, vdiff);
            for (int k = 0; k < 8; ++k) tmp[k] = std::exp(tmp[k]);
            sum += _mm512_reduce_add_pd(_mm512_loadu_pd(tmp));
        }
        for (; i < n; ++i)
            sum += std::exp(values[i] - x_max);
        return x_max + std::log(sum);
    }
#endif

    // Eigen path
    {
        Eigen::Map<const Eigen::ArrayXd> x(values.data(), n);
        return x_max + std::log((x - x_max).exp().sum());
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Distance / RMSD dispatch
// ═════════════════════════════════════════════════════════════════════════════

static inline float sq(float x) { return x * x; }

float UnifiedHardwareDispatch::rmsd_scalar(const float* a, const float* b, int n) {
    float sum = 0.0f;
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < 3; ++c)
            sum += sq(a[i * 3 + c] - b[i * 3 + c]);
    return std::sqrt(sum / static_cast<float>(n));
}

float UnifiedHardwareDispatch::rmsd_avx2(const float* a, const float* b, int n) {
#ifdef __AVX2__
    __m256 acc = _mm256_setzero_ps();
    int i = 0;
    for (; i <= n - 8; i += 8) {
        for (int c = 0; c < 3; ++c) {
            float a8[8], b8[8];
            for (int k = 0; k < 8; ++k) {
                a8[k] = a[(i + k) * 3 + c];
                b8[k] = b[(i + k) * 3 + c];
            }
            __m256 da = _mm256_sub_ps(_mm256_loadu_ps(a8), _mm256_loadu_ps(b8));
            acc = _mm256_fmadd_ps(da, da, acc);
        }
    }
    __m128 lo = _mm256_castps256_ps128(acc);
    __m128 hi = _mm256_extractf128_ps(acc, 1);
    lo = _mm_add_ps(lo, hi);
    lo = _mm_add_ps(lo, _mm_movehl_ps(lo, lo));
    lo = _mm_add_ss(lo, _mm_movehdup_ps(lo));
    float sum = _mm_cvtss_f32(lo);
    for (; i < n; ++i)
        for (int c = 0; c < 3; ++c)
            sum += sq(a[i * 3 + c] - b[i * 3 + c]);
    return std::sqrt(sum / static_cast<float>(n));
#else
    return rmsd_scalar(a, b, n);
#endif
}

float UnifiedHardwareDispatch::rmsd_avx512(const float* a, const float* b, int n) {
#ifdef __AVX512F__
    __m512 acc = _mm512_setzero_ps();
    int i = 0;
    for (; i <= n - 16; i += 16) {
        for (int c = 0; c < 3; ++c) {
            float a16[16], b16[16];
            for (int k = 0; k < 16; ++k) {
                a16[k] = a[(i + k) * 3 + c];
                b16[k] = b[(i + k) * 3 + c];
            }
            __m512 da = _mm512_sub_ps(_mm512_loadu_ps(a16), _mm512_loadu_ps(b16));
            acc = _mm512_fmadd_ps(da, da, acc);
        }
    }
    float sum = _mm512_reduce_add_ps(acc);
    for (; i < n; ++i)
        for (int c = 0; c < 3; ++c)
            sum += sq(a[i * 3 + c] - b[i * 3 + c]);
    return std::sqrt(sum / static_cast<float>(n));
#else
    return rmsd_avx2(a, b, n);
#endif
}

float UnifiedHardwareDispatch::rmsd_openmp(const float* a, const float* b, int n) {
#ifdef _OPENMP
    float sum = 0.0f;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < n; ++i)
        for (int c = 0; c < 3; ++c)
            sum += sq(a[i * 3 + c] - b[i * 3 + c]);
    return std::sqrt(sum / static_cast<float>(n));
#else
    return rmsd_scalar(a, b, n);
#endif
}

void UnifiedHardwareDispatch::distance2_batch(
    const float* ax, const float* ay, const float* az,
    float bx, float by, float bz,
    float* out, int n, Backend backend)
{
    if (!detected_) detect();
    [[maybe_unused]] Backend b = (backend == Backend::AUTO) ? best_backend(KernelType::DISTANCE_BATCH) : backend;
#ifdef __AVX512F__
    if (b == Backend::AVX512 && is_available(Backend::AVX512)) {
        __m512 vbx = _mm512_set1_ps(bx);
        __m512 vby = _mm512_set1_ps(by);
        __m512 vbz = _mm512_set1_ps(bz);
        int i = 0;
        for (; i + 15 < n; i += 16) {
            __m512 dx = _mm512_sub_ps(_mm512_loadu_ps(ax + i), vbx);
            __m512 dy = _mm512_sub_ps(_mm512_loadu_ps(ay + i), vby);
            __m512 dz = _mm512_sub_ps(_mm512_loadu_ps(az + i), vbz);
            __m512 r2 = _mm512_fmadd_ps(dz, dz,
                        _mm512_fmadd_ps(dy, dy,
                        _mm512_mul_ps(dx, dx)));
            _mm512_storeu_ps(out + i, r2);
        }
        for (; i < n; ++i)
            out[i] = sq(ax[i] - bx) + sq(ay[i] - by) + sq(az[i] - bz);
        return;
    }
#endif

#ifdef __AVX2__
    if ((b == Backend::AVX2 || b == Backend::AVX512) && is_available(Backend::AVX2)) {
        __m256 vbx = _mm256_set1_ps(bx);
        __m256 vby = _mm256_set1_ps(by);
        __m256 vbz = _mm256_set1_ps(bz);
        int i = 0;
        for (; i + 7 < n; i += 8) {
            __m256 dx = _mm256_sub_ps(_mm256_loadu_ps(ax + i), vbx);
            __m256 dy = _mm256_sub_ps(_mm256_loadu_ps(ay + i), vby);
            __m256 dz = _mm256_sub_ps(_mm256_loadu_ps(az + i), vbz);
            __m256 r2 = _mm256_fmadd_ps(dz, dz,
                        _mm256_fmadd_ps(dy, dy,
                        _mm256_mul_ps(dx, dx)));
            _mm256_storeu_ps(out + i, r2);
        }
        for (; i < n; ++i)
            out[i] = sq(ax[i] - bx) + sq(ay[i] - by) + sq(az[i] - bz);
        return;
    }
#endif

    for (int i = 0; i < n; ++i)
        out[i] = sq(ax[i] - bx) + sq(ay[i] - by) + sq(az[i] - bz);
}

float UnifiedHardwareDispatch::rmsd(const float* a_xyz, const float* b_xyz,
                                     int n_atoms, Backend backend) {
    if (n_atoms <= 0) return 0.0f;
    if (!detected_) detect();

    Backend b = (backend == Backend::AUTO) ? best_backend(KernelType::DISTANCE_BATCH) : backend;

    switch (b) {
        case Backend::AVX512: return rmsd_avx512(a_xyz, b_xyz, n_atoms);
        case Backend::AVX2:   return rmsd_avx2(a_xyz, b_xyz, n_atoms);
        case Backend::OPENMP: return rmsd_openmp(a_xyz, b_xyz, n_atoms);
        default:              return rmsd_scalar(a_xyz, b_xyz, n_atoms);
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// Dispatch report
// ═════════════════════════════════════════════════════════════════════════════

DispatchReport UnifiedHardwareDispatch::get_dispatch_report() const {
    Backend sel = best_backend(KernelType::FITNESS_EVAL);
    const auto& hwd = flexaids::detect_hardware();

    std::string reason;
    switch (sel) {
        case Backend::CUDA:
            reason = "CUDA GPU detected (" + hwd.cuda_device_name + ", " + hwd.cuda_arch + ")";
            break;
        case Backend::METAL:
            reason = "Metal GPU detected (" + hwd.metal_gpu_name + ")";
#ifdef FLEXAIDS_HAS_METAL_SHANNON
            reason += " [Boltzmann + LogSumExp + Histogram kernels cached]";
#endif
            break;
        case Backend::ROCM:
            reason = "ROCm GPU detected (" + hwd.rocm_device_name + ", " + hwd.rocm_arch + ")";
            break;
        case Backend::AVX512:
            reason = "AVX-512 F+DQ+BW detected on CPU";
            break;
        case Backend::AVX2:
            reason = "AVX2+FMA detected on CPU";
            break;
        case Backend::OPENMP:
            reason = "OpenMP with " + std::to_string(hwd.openmp_max_threads) + " threads";
            break;
        case Backend::SCALAR:
            reason = "No acceleration available; using scalar baseline";
            break;
        default:
            reason = "auto";
            break;
    }

    return { sel, reason, hwd.summary() };
}

}  // namespace hw

// ═════════════════════════════════════════════════════════════════════════════
// flexaids:: compatibility free functions
// ═════════════════════════════════════════════════════════════════════════════

namespace flexaids {

const char* backend_name(HardwareBackend b) noexcept {
    return hw::UnifiedHardwareDispatch::backend_name(
        static_cast<hw::Backend>(static_cast<uint8_t>(b)));
}

HardwareBackend select_backend() {
    auto& d = hw::UnifiedHardwareDispatch::instance();
    d.detect();
    hw::Backend b = d.best_backend(hw::KernelType::FITNESS_EVAL);
    return static_cast<HardwareBackend>(static_cast<uint8_t>(b));
}

HardwareBackend select_cpu_backend() {
    auto& d = hw::UnifiedHardwareDispatch::instance();
    d.detect();
    hw::Backend b = d.select_cpu_backend();
    return static_cast<HardwareBackend>(static_cast<uint8_t>(b));
}

BoltzmannBatchResult compute_boltzmann_batch(
    std::span<const double> energies, double beta)
{
    return hw::UnifiedHardwareDispatch::instance().compute_boltzmann_batch(energies, beta);
}

double log_sum_exp_dispatch(std::span<const double> values) {
    return hw::UnifiedHardwareDispatch::instance().log_sum_exp_dispatch(values);
}

DispatchReport get_dispatch_report() {
    auto& d = hw::UnifiedHardwareDispatch::instance();
    d.detect();
    return d.get_dispatch_report();
}

}  // namespace flexaids
