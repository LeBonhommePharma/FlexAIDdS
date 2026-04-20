// ROCmDispatch.h — ROCm/HIP backend additions to FlexAIDdS HardwareDispatch
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// This header extends the existing hardware_dispatch.h with AMD ROCm/HIP support.
// It adds:
//   - BackendType::ROCM = 6 (fits between CUDA=5 and existing values)
//   - HardwareInfo::has_rocm, rocm_device_name, rocm_compute_units, rocm_global_mem
//   - Updated backend priority: CUDA > ROCm > Metal > AVX-512 > AVX-2 > OpenMP > Scalar
//   - Runtime detection via hipGetDeviceCount / hipGetDeviceProperties
//   - Stub ROCm dispatch class mirroring the CUDA dispatch interface
//
// HIP code is source-compatible with CUDA. Kernels in hip_eval.hip can be
// compiled with either nvcc (as .cu) or hipcc (as .hip) — same source.
//
// Include this header alongside hardware_dispatch.h. It does NOT modify that
// file; instead it provides an extension class ROCmDispatch and a free function
// detect_rocm() that updates an existing HardwareInfo struct.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Conditional HIP runtime header inclusion
// HIP headers are only available when compiled with hipcc or when ROCm is
// installed. Guard with __HIP_PLATFORM_HCC__ or __HIPCC__.
// ---------------------------------------------------------------------------

#if defined(__HIP_PLATFORM_HCC__) || defined(__HIPCC__) || defined(HAVE_HIP)
#  include <hip/hip_runtime.h>
#  define FLEXAIDDS_HAVE_HIP 1
#else
// Forward-declare minimal HIP types for compilation without ROCm installed.
// These allow the header to be included in CPU-only builds.
struct hipDeviceProp_t {
    char   name[256];
    size_t totalGlobalMem;
    int    multiProcessorCount;
    int    major, minor;
    int    warpSize;
    int    maxThreadsPerBlock;
};
using hipError_t = int;
constexpr hipError_t hipSuccess = 0;
#endif

namespace flexaidds {

// ---------------------------------------------------------------------------
// Extended BackendType enum value
// ---------------------------------------------------------------------------
// In the existing hardware_dispatch.h, BackendType is:
//   SCALAR=0, OPENMP=1, AVX2=2, AVX512=3, METAL=4, CUDA=5
// We add ROCM=6. When merging with the existing header, replace:
//   enum class BackendType { SCALAR=0, OPENMP=1, AVX2=2, AVX512=3, METAL=4, CUDA=5 };
// with:
//   enum class BackendType { SCALAR=0, OPENMP=1, AVX2=2, AVX512=3, METAL=4, CUDA=5, ROCM=6 };

constexpr int ROCM_BACKEND_ID = 6; // BackendType::ROCM

// ---------------------------------------------------------------------------
// ROCm device info (added to HardwareInfo)
// ---------------------------------------------------------------------------

struct ROCmDeviceInfo {
    bool        has_rocm         = false;
    int         device_count     = 0;
    int         active_device    = 0;        // index of the selected device
    std::string device_name;                 // e.g. "AMD Radeon RX 7900 XTX"
    int         compute_units    = 0;        // hipDeviceProp_t::multiProcessorCount
    size_t      global_mem_bytes = 0;        // hipDeviceProp_t::totalGlobalMem
    int         gcn_arch_major   = 0;        // hipDeviceProp_t::major
    int         gcn_arch_minor   = 0;        // hipDeviceProp_t::minor
    int         warp_size        = 64;       // AMD wavefront = 64 threads
    int         max_threads_per_block = 1024;
};

// ---------------------------------------------------------------------------
// Runtime detection
// ---------------------------------------------------------------------------

/// Probe the ROCm runtime and fill a ROCmDeviceInfo.
/// Returns false if ROCm is not available or no devices found.
/// Safe to call on any platform (returns false on non-AMD systems).
bool detect_rocm(ROCmDeviceInfo& info);

/// Select the best ROCm device (highest compute units × global memory product).
/// Returns device index. Sets info.active_device.
int select_rocm_device(ROCmDeviceInfo& info);

// ---------------------------------------------------------------------------
// Backend priority helper
// ---------------------------------------------------------------------------

/// Returns the priority rank of each backend (higher = preferred).
/// Updated priority table:
///   CUDA=100, ROCM=90, METAL=80, AVX512=70, AVX2=60, OPENMP=50, SCALAR=0
constexpr int backend_priority(int backend_id) noexcept {
    switch (backend_id) {
        case 5: return 100; // CUDA
        case 6: return 90;  // ROCM (new)
        case 4: return 80;  // METAL
        case 3: return 70;  // AVX512
        case 2: return 60;  // AVX2
        case 1: return 50;  // OPENMP
        case 0: return 0;   // SCALAR
        default: return -1;
    }
}

// ---------------------------------------------------------------------------
// ROCmDispatch — GPU dispatch class for AMD hardware
// Mirrors the interface of the existing CudaDispatch class.
// ---------------------------------------------------------------------------

class ROCmDispatch {
public:
    /// Initialise the ROCm backend on the given device index.
    /// Throws std::runtime_error if ROCm is unavailable or device invalid.
    explicit ROCmDispatch(int device_index = 0);

    ~ROCmDispatch();

    // Non-copyable, movable
    ROCmDispatch(const ROCmDispatch&)            = delete;
    ROCmDispatch& operator=(const ROCmDispatch&) = delete;
    ROCmDispatch(ROCmDispatch&&) noexcept;
    ROCmDispatch& operator=(ROCmDispatch&&) noexcept;

    // -----------------------------------------------------------------------
    // Memory management
    // -----------------------------------------------------------------------

    /// Allocate device memory (calls hipMalloc).
    void* device_alloc(size_t bytes);

    /// Free device memory (calls hipFree).
    void device_free(void* ptr) noexcept;

    /// Copy host → device.
    void copy_to_device(void* dst, const void* src, size_t bytes);

    /// Copy device → host.
    void copy_from_device(void* dst, const void* src, size_t bytes);

    // -----------------------------------------------------------------------
    // Kernel launch interface
    // -----------------------------------------------------------------------

    /// Evaluate energy for a population of poses.
    /// Mirrors cuda_eval.cu launch_eval_kernel().
    ///   d_coords:       float3* device array of atom coordinates (N_atoms × N_poses)
    ///   d_types:        uint8_t* device array of atom types (N_atoms)
    ///   d_energies:     float* output device array (N_poses)
    ///   n_poses:        number of poses to evaluate
    ///   n_atoms:        atoms per pose
    ///   block_size:     threads per block (default = min(warp_size*2, max_threads))
    void launch_eval_kernel(const void* d_coords, const void* d_types,
                             void* d_energies,
                             int n_poses, int n_atoms,
                             int block_size = 128);

    /// Synchronise device (hipDeviceSynchronize).
    void synchronize();

    /// Return device info.
    const ROCmDeviceInfo& device_info() const noexcept { return info_; }

    /// Return true if the backend is operational.
    bool is_available() const noexcept { return available_; }

private:
    ROCmDeviceInfo info_;
    bool           available_ = false;
    [[maybe_unused]] int            device_id_ = 0;

    void check_hip(hipError_t err, const char* msg);
};

// ---------------------------------------------------------------------------
// HardwareDispatch extension helper
// ---------------------------------------------------------------------------
// Call this in hardware_dispatch.cpp's detect_hardware() to probe ROCm after
// CUDA detection. If CUDA is available, ROCm is skipped (GPU already claimed).

struct HardwareDetectContext {
    bool cuda_available   = false;
    bool rocm_available   = false;
    bool metal_available  = false;
    bool avx512_available = false;
    bool avx2_available   = false;
    bool openmp_available = false;
};

/// Returns the best backend_id given availability flags.
/// Implements the full priority table.
int select_best_backend(const HardwareDetectContext& ctx) noexcept;

} // namespace flexaidds
