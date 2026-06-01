// hardware_detect.h — Runtime hardware capability detection for FlexAIDdS
//
// Probes the execution environment at runtime to determine available
// acceleration backends (CUDA, Metal, AVX-512, AVX2, OpenMP) and their
// properties.  The result is cached after the first call.
//
// Usage:
//   const auto& hw = flexaids::detect_hardware();
//   if (hw.has_avx512) { /* use 512-bit SIMD path */ }
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#pragma once

#include <string>
#include <cstdint>

namespace flexaids {

// Collected hardware capabilities of the current execution environment.
struct HardwareCapabilities {
    // ── GPU: CUDA ──
    bool        has_cuda          = false;
    int         cuda_device_count = 0;
    std::string cuda_device_name;       // e.g. "NVIDIA GeForce RTX 4090"
    std::string cuda_arch;              // e.g. "sm_89"
    int         cuda_sm_major     = 0;
    int         cuda_sm_minor     = 0;
    std::size_t cuda_global_mem   = 0;  // bytes

    // ── GPU: ROCm/HIP (AMD) ──
    bool        has_rocm          = false;
    int         rocm_device_count = 0;
    std::string rocm_device_name;       // e.g. "AMD Instinct MI300X"
    std::string rocm_arch;              // e.g. "gfx90a"
    int         rocm_compute_units = 0;
    std::size_t rocm_global_mem    = 0; // bytes
    int         rocm_wavefront     = 0; // typically 64

    // ── GPU: Metal (Apple) ──
    bool        has_metal               = false;
    std::string metal_gpu_name;         // e.g. "Apple M3 Max"
    std::size_t metal_unified_memory    = 0;   // bytes (from MetalCapabilities)
    std::size_t metal_max_buffer        = 0;   // bytes
    int         metal_gpu_core_estimate = 0;   // rough core count (M3/M3 Pro/etc.)

    // ── SIMD: x86-64 ──
    bool has_sse42   = false;
    bool has_avx2    = false;
    bool has_fma     = false;
    bool has_avx512f = false;            // foundation
    bool has_avx512dq = false;           // doubleword/quadword
    bool has_avx512bw = false;           // byte/word
    bool has_avx512vnni = false;         // Vector Neural Network Instructions

    // Convenience composite
    bool has_avx512  = false;            // avx512f && avx512dq && avx512bw

    // ── OpenMP ──
    bool has_openmp       = false;
    int  openmp_max_threads = 1;

    // ── Eigen ──
    bool has_eigen = false;

    // Human-readable summary (one line per backend)
    std::string summary() const;
};

// Detect hardware capabilities (cached after first call).
const HardwareCapabilities& detect_hardware();

}  // namespace flexaids
