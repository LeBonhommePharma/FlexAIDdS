// ROCmDispatch.cpp — ROCm/HIP backend implementation
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Compiled normally (without hipcc) on CPU-only machines: all hipXxx calls
// are guarded and will return graceful no-ops when ROCm is not present.
// When compiled with hipcc, FLEXAIDDS_HAVE_HIP is defined and the real HIP
// runtime API is used.

#include "ROCmDispatch.h"

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <string>

// ---------------------------------------------------------------------------
// Platform stub macros for CPU-only builds
// ---------------------------------------------------------------------------

#ifndef FLEXAIDDS_HAVE_HIP
// These stubs make the file compile on any platform without ROCm.
// They are never actually called (all paths return early if !has_rocm).

[[maybe_unused]] static int  hipGetDeviceCount_stub(int* cnt) { *cnt = 0; return 1; }
[[maybe_unused]] static int  hipSetDevice_stub(int)            { return 1; }
[[maybe_unused]] static int  hipGetDeviceProperties_stub(hipDeviceProp_t*, int) { return 1; }[[maybe_unused]] static int  hipMalloc_stub(void**, size_t)    { return 1; }
[[maybe_unused]] static int  hipFree_stub(void*)               { return 1; }
[[maybe_unused]] static int  hipMemcpy_stub(void*, const void*, size_t, int) { return 1; }
[[maybe_unused]] static int  hipDeviceSynchronize_stub()       { return 1; }

#define hipGetDeviceCount    hipGetDeviceCount_stub
#define hipSetDevice         hipSetDevice_stub
#define hipGetDeviceProperties hipGetDeviceProperties_stub
#define hipMalloc            hipMalloc_stub
#define hipFree              hipFree_stub
#define hipMemcpy            hipMemcpy_stub
#define hipDeviceSynchronize hipDeviceSynchronize_stub

// Stub copy-direction enum[[maybe_unused]] constexpr int hipMemcpyHostToDevice   = 1;
[[maybe_unused]] constexpr int hipMemcpyDeviceToHost   = 2;
[[maybe_unused]] constexpr int hipMemcpyDeviceToDevice = 3;

#endif // FLEXAIDDS_HAVE_HIP

namespace flexaidds {

// ---------------------------------------------------------------------------
// Runtime detection
// ---------------------------------------------------------------------------

bool detect_rocm(ROCmDeviceInfo& info) {
    info = ROCmDeviceInfo{};

    int count = 0;
    auto err = hipGetDeviceCount(&count);
    if (err != hipSuccess || count <= 0) return false;

    // Verify /opt/rocm exists on Linux (belt-and-suspenders for Docker envs)
    // Skip this check on non-Linux platforms.
#if defined(__linux__)
    {
        FILE* f = std::fopen("/opt/rocm/include/hip/hip_version.h", "r");
        if (!f) {
            // /opt/rocm not found, but hipGetDeviceCount succeeded —
            // could be CUDA with HIP shims. Accept anyway.
        } else {
            std::fclose(f);
        }
    }
#endif

    info.has_rocm     = true;
    info.device_count = count;

    // Load info for device 0 initially; select_rocm_device() may change this
    hipDeviceProp_t props;
    std::memset(&props, 0, sizeof(props));
    err = hipGetDeviceProperties(&props, 0);
    if (err == hipSuccess) {
        info.device_name      = std::string(props.name);
        info.compute_units    = props.multiProcessorCount;
        info.global_mem_bytes = props.totalGlobalMem;
        info.gcn_arch_major   = props.major;
        info.gcn_arch_minor   = props.minor;
        info.warp_size        = (props.warpSize > 0) ? props.warpSize : 64;
        info.max_threads_per_block = props.maxThreadsPerBlock;
    }

    return true;
}

int select_rocm_device(ROCmDeviceInfo& info) {
    if (!info.has_rocm || info.device_count <= 0) return 0;

    int best_device = 0;
    long long best_score = -1;

    for (int d = 0; d < info.device_count; ++d) {
        hipDeviceProp_t props;
        std::memset(&props, 0, sizeof(props));
        if (hipGetDeviceProperties(&props, d) != hipSuccess) continue;

        long long score = static_cast<long long>(props.multiProcessorCount)
                        * static_cast<long long>(props.totalGlobalMem >> 20); // MB units
        if (score > best_score) {
            best_score = score;
            best_device = d;
        }
    }

    // Load properties for best device
    hipDeviceProp_t props;
    std::memset(&props, 0, sizeof(props));
    if (hipGetDeviceProperties(&props, best_device) == hipSuccess) {
        info.active_device         = best_device;
        info.device_name           = std::string(props.name);
        info.compute_units         = props.multiProcessorCount;
        info.global_mem_bytes      = props.totalGlobalMem;
        info.gcn_arch_major        = props.major;
        info.gcn_arch_minor        = props.minor;
        info.warp_size             = (props.warpSize > 0) ? props.warpSize : 64;
        info.max_threads_per_block = props.maxThreadsPerBlock;
    }

    return best_device;
}

// ---------------------------------------------------------------------------
// Backend priority selector
// ---------------------------------------------------------------------------

int select_best_backend(const HardwareDetectContext& ctx) noexcept {
    if (ctx.cuda_available)   return 5; // CUDA
    if (ctx.rocm_available)   return 6; // ROCM
    if (ctx.metal_available)  return 4; // METAL
    if (ctx.avx512_available) return 3; // AVX512
    if (ctx.avx2_available)   return 2; // AVX2
    if (ctx.openmp_available) return 1; // OPENMP
    return 0;                           // SCALAR
}

// ---------------------------------------------------------------------------
// ROCmDispatch constructor
// ---------------------------------------------------------------------------

ROCmDispatch::ROCmDispatch(int device_index) : device_id_(device_index) {
    available_ = detect_rocm(info_);
    if (!available_) {
        throw std::runtime_error("ROCmDispatch: ROCm runtime not available or no AMD GPU found");
    }
    if (device_index < 0 || device_index >= info_.device_count) {
        throw std::runtime_error("ROCmDispatch: device index " +
                                 std::to_string(device_index) + " out of range (found " +
                                 std::to_string(info_.device_count) + " device(s))");
    }

#ifdef FLEXAIDDS_HAVE_HIP
    hipError_t err = hipSetDevice(device_index);
    if (err != hipSuccess) {
        throw std::runtime_error("ROCmDispatch: hipSetDevice(" +
                                 std::to_string(device_index) + ") failed");
    }
#endif

    // Re-query properties for the selected device
    hipDeviceProp_t props;
    std::memset(&props, 0, sizeof(props));
    hipGetDeviceProperties(&props, device_index);
    info_.active_device         = device_index;
    info_.device_name           = std::string(props.name);
    info_.compute_units         = props.multiProcessorCount;
    info_.global_mem_bytes      = props.totalGlobalMem;
    info_.gcn_arch_major        = props.major;
    info_.gcn_arch_minor        = props.minor;
    info_.warp_size             = (props.warpSize > 0) ? props.warpSize : 64;
    info_.max_threads_per_block = props.maxThreadsPerBlock;
}

ROCmDispatch::~ROCmDispatch() = default;

ROCmDispatch::ROCmDispatch(ROCmDispatch&&) noexcept = default;
ROCmDispatch& ROCmDispatch::operator=(ROCmDispatch&&) noexcept = default;

// ---------------------------------------------------------------------------
// Error checking helper
// ---------------------------------------------------------------------------

void ROCmDispatch::check_hip(hipError_t err, const char* msg) {
    if (err != hipSuccess) {
        throw std::runtime_error(std::string("ROCmDispatch HIP error in ") +
                                 msg + ": code " + std::to_string(err));
    }
}

// ---------------------------------------------------------------------------
// Memory management
// ---------------------------------------------------------------------------

void* ROCmDispatch::device_alloc(size_t bytes) {
    void* ptr = nullptr;
#ifdef FLEXAIDDS_HAVE_HIP
    check_hip(hipMalloc(&ptr, bytes), "device_alloc/hipMalloc");
#else
    throw std::runtime_error("ROCmDispatch::device_alloc: HIP not compiled in");
#endif
    return ptr;
}

void ROCmDispatch::device_free(void* ptr) noexcept {
    if (!ptr) return;
#ifdef FLEXAIDDS_HAVE_HIP
    hipFree(ptr);
#endif
}

void ROCmDispatch::copy_to_device(void* dst, const void* src, size_t bytes) {
#ifdef FLEXAIDDS_HAVE_HIP
    check_hip(hipMemcpy(dst, src, bytes, hipMemcpyHostToDevice),
              "copy_to_device/hipMemcpy");
#else
    throw std::runtime_error("ROCmDispatch::copy_to_device: HIP not compiled in");
#endif
}

void ROCmDispatch::copy_from_device(void* dst, const void* src, size_t bytes) {
#ifdef FLEXAIDDS_HAVE_HIP
    check_hip(hipMemcpy(dst, src, bytes, hipMemcpyDeviceToHost),
              "copy_from_device/hipMemcpy");
#else
    throw std::runtime_error("ROCmDispatch::copy_from_device: HIP not compiled in");
#endif
}

// ---------------------------------------------------------------------------
// Kernel launch
// ---------------------------------------------------------------------------

// Forward-declare the kernel launcher from hip_eval.hip
// (compiled separately with hipcc when HAVE_HIP is defined)
#ifdef FLEXAIDDS_HAVE_HIP
extern "C" void hip_eval_launch(const float*   d_coords,
                                 const uint8_t* d_types,
                                 float*         d_energies,
                                 int            n_poses,
                                 int            n_atoms,
                                 int            block_size);
#endif

void ROCmDispatch::launch_eval_kernel(const void* d_coords, const void* d_types,
                                       void* d_energies,
                                       int n_poses, int n_atoms,
                                       int block_size) {
    if (!available_) {
        throw std::runtime_error("ROCmDispatch::launch_eval_kernel: ROCm not available");
    }

    // Clamp block_size to hardware limits
    if (block_size <= 0 || block_size > info_.max_threads_per_block)
        block_size = std::min(128, info_.max_threads_per_block);

#ifdef FLEXAIDDS_HAVE_HIP
    hip_eval_launch(static_cast<const float*>(d_coords),
                    static_cast<const uint8_t*>(d_types),
                    static_cast<float*>(d_energies),
                    n_poses, n_atoms, block_size);
#else
    throw std::runtime_error("ROCmDispatch::launch_eval_kernel: HIP not compiled in; "
                             "recompile with hipcc and -DHAVE_HIP");
#endif
}

void ROCmDispatch::synchronize() {
#ifdef FLEXAIDDS_HAVE_HIP
    check_hip(hipDeviceSynchronize(), "synchronize");
#endif
}

} // namespace flexaidds
