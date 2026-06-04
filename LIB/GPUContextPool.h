// GPUContextPool.h — Thread-safe singleton pool for GPU evaluation contexts
//
// GPU contexts (CudaEvalCtx, MetalEvalCtx) are expensive to create (~100ms).
// This pool manages a single context per GPU backend, allowing multiple
// concurrent GA instances to share the same device context safely.
//
// Usage:
//   auto& pool = GPUContextPool::instance();
//   auto handle = pool.acquire_cuda(n_atoms, n_types, ...);
//   // ... use handle.ctx ...
//   pool.release_cuda(handle);
//
// For ParallelDock (same system, different regions): all threads share one
// context since n_atoms and n_types are identical.
//
// For ParallelCampaign (different ligands): acquire() waits for exclusive
// access when dimensions change, serializing GPU rebuilds naturally.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <mutex>
#include <condition_variable>
#include <functional>
#include <cstddef>

#ifdef FLEXAIDS_USE_CUDA
#include "cuda_eval.cuh"
#endif

#ifdef FLEXAIDS_USE_METAL
#include "metal_eval.h"
#endif

#ifdef FLEXAIDS_USE_ROCM
#include "hip_eval.h"
#endif

class GPUContextPool {
public:
    static GPUContextPool& instance() {
        static GPUContextPool pool;
        return pool;
    }

#ifdef FLEXAIDS_USE_CUDA
    struct CudaHandle {
        CudaEvalCtx* ctx;
    };

    // Acquire the shared CUDA context, rebuilding if dimensions changed.
    // Blocks if another thread is rebuilding.
    // init_fn is called to create/recreate the context when dimensions change.
    CudaHandle acquire_cuda(
        int n_atoms, int n_types,
        std::function<CudaEvalCtx*(void)> init_fn
    ) {
        std::unique_lock<std::mutex> lock(cuda_mtx_);
        // Wait if someone is rebuilding with different dimensions
        cuda_cv_.wait(lock, [&] {
            return !cuda_rebuilding_ || (cuda_natom_ == n_atoms && cuda_ntype_ == n_types);
        });

        if (!cuda_ctx_ || cuda_natom_ != n_atoms || cuda_ntype_ != n_types) {
            cuda_rebuilding_ = true;
            // Wait for all users to release before rebuilding
            cuda_cv_.wait(lock, [&] { return cuda_ref_count_ == 0; });

            if (cuda_ctx_) {
                cuda_eval_shutdown(cuda_ctx_);
                cuda_ctx_ = nullptr;
            }
            lock.unlock();
            CudaEvalCtx* new_ctx = nullptr;
            try {
                new_ctx = init_fn();
            } catch (...) {
                lock.lock();
                cuda_rebuilding_ = false;
                cuda_cv_.notify_all();
                throw;
            }
            lock.lock();

            cuda_ctx_   = new_ctx;
            cuda_natom_ = n_atoms;
            cuda_ntype_ = n_types;
            cuda_rebuilding_ = false;
            cuda_cv_.notify_all();
        }

        cuda_ref_count_++;
        return CudaHandle{cuda_ctx_};
    }

    void release_cuda(const CudaHandle& /*handle*/) {
        std::lock_guard<std::mutex> lock(cuda_mtx_);
        cuda_ref_count_--;
        if (cuda_ref_count_ == 0) {
            cuda_cv_.notify_all();
        }
    }
#endif

#ifdef FLEXAIDS_USE_METAL
    struct MetalHandle {
        MetalEvalCtx* ctx;
    };

    MetalHandle acquire_metal(
        int n_atoms, int n_types, int pop_size,
        std::function<MetalEvalCtx*(void)> init_fn
    ) {
        std::unique_lock<std::mutex> lock(metal_mtx_);
        metal_cv_.wait(lock, [&] {
            return !metal_rebuilding_ || (metal_natom_ == n_atoms && metal_ntype_ == n_types
                                          && metal_max_pop_ >= pop_size);
        });

        if (!metal_ctx_ || metal_natom_ != n_atoms || metal_ntype_ != n_types
                        || metal_max_pop_ < pop_size) {
            metal_rebuilding_ = true;
            metal_cv_.wait(lock, [&] { return metal_ref_count_ == 0; });

            if (metal_ctx_) {
                metal_eval_shutdown(metal_ctx_);
                metal_ctx_ = nullptr;
            }
            lock.unlock();
            MetalEvalCtx* new_ctx = nullptr;
            try {
                new_ctx = init_fn();
            } catch (...) {
                lock.lock();
                metal_rebuilding_ = false;
                metal_cv_.notify_all();
                throw;
            }
            lock.lock();

            metal_ctx_     = new_ctx;
            metal_natom_   = n_atoms;
            metal_ntype_   = n_types;
            metal_max_pop_ = pop_size;
            metal_rebuilding_ = false;
            metal_cv_.notify_all();
        }

        metal_ref_count_++;
        return MetalHandle{metal_ctx_};
    }

    void release_metal(const MetalHandle& /*handle*/) {
        std::lock_guard<std::mutex> lock(metal_mtx_);
        metal_ref_count_--;
        if (metal_ref_count_ == 0) {
            metal_cv_.notify_all();
        }
    }
#endif

#ifdef FLEXAIDS_USE_ROCM
    // ─── ROCm/HIP context management ────────────────────────────────────────
    //
    // The hip_eval backend uses module-level state (hip_eval::init/shutdown)
    // rather than an opaque context pointer.  The RocmHandle carries the
    // device_id so callers can route work to the correct device.
    //
    // Dimension tracking uses max_pop + num_genes (the two axes that trigger
    // reallocation inside hip_eval).

    struct RocmHandle {
        int device_id;
    };

    // Acquire the shared ROCm context, re-initialising if the device or
    // population dimensions changed.  Thread-safe: blocks if another thread
    // is rebuilding.
    //
    // init_fn is called to (re-)initialise the ROCm backend.  It receives
    // the device_id and max_pop and should call hip_eval::init().
    RocmHandle acquire_rocm(
        int device_id, int max_pop,
        std::function<void(int /*device_id*/, int /*max_pop*/)> init_fn
    ) {
        std::unique_lock<std::mutex> lock(rocm_mtx_);
        rocm_cv_.wait(lock, [&] {
            return !rocm_rebuilding_ ||
                   (rocm_device_id_ == device_id && rocm_max_pop_ >= max_pop);
        });

        if (!rocm_initialised_ ||
            rocm_device_id_ != device_id ||
            rocm_max_pop_ < max_pop)
        {
            rocm_rebuilding_ = true;
            // Wait for all users to release before rebuilding
            rocm_cv_.wait(lock, [&] { return rocm_ref_count_ == 0; });

            if (rocm_initialised_) {
                hip_eval::shutdown();
                rocm_initialised_ = false;
            }
            lock.unlock();
            try {
                init_fn(device_id, max_pop);
            } catch (...) {
                lock.lock();
                rocm_rebuilding_ = false;
                rocm_cv_.notify_all();
                throw;
            }
            lock.lock();

            rocm_device_id_   = device_id;
            rocm_max_pop_     = max_pop;
            rocm_initialised_ = true;
            rocm_rebuilding_  = false;
            rocm_cv_.notify_all();
        }

        rocm_ref_count_++;
        return RocmHandle{rocm_device_id_};
    }

    void release_rocm(const RocmHandle& /*handle*/) {
        std::lock_guard<std::mutex> lock(rocm_mtx_);
        rocm_ref_count_--;
        if (rocm_ref_count_ == 0) {
            rocm_cv_.notify_all();
        }
    }

    // ─── ROCm device enumeration helpers ─────────────────────────────────────

    // Return the number of available HIP devices, or 0 if ROCm is unavailable.
    static int rocm_device_count() {
        int count = 0;
        if (hipGetDeviceCount(&count) != hipSuccess) return 0;
        return count;
    }

    // Query free and total memory on a HIP device.
    // Returns false if the query fails.
    static bool rocm_mem_info(int device_id, std::size_t& free_bytes, std::size_t& total_bytes) {
        if (hipSetDevice(device_id) != hipSuccess) return false;
        if (hipMemGetInfo(&free_bytes, &total_bytes) != hipSuccess) return false;
        return true;
    }

    // Return the device index with the most free memory (best default for
    // compute workloads).  Returns -1 if no device is available.
    static int rocm_best_device() {
        int count = rocm_device_count();
        if (count == 0) return -1;

        int         best_idx  = 0;
        std::size_t best_free = 0;

        for (int i = 0; i < count; ++i) {
            std::size_t free_b = 0, total_b = 0;
            if (rocm_mem_info(i, free_b, total_b) && free_b > best_free) {
                best_free = free_b;
                best_idx  = i;
            }
        }
        return best_idx;
    }
#endif

private:
    GPUContextPool() = default;
    ~GPUContextPool() {
#ifdef FLEXAIDS_USE_CUDA
        if (cuda_ctx_) cuda_eval_shutdown(cuda_ctx_);
#endif
#ifdef FLEXAIDS_USE_METAL
        if (metal_ctx_) metal_eval_shutdown(metal_ctx_);
#endif
#ifdef FLEXAIDS_USE_ROCM
        if (rocm_initialised_) hip_eval::shutdown();
#endif
    }

    GPUContextPool(const GPUContextPool&) = delete;
    GPUContextPool& operator=(const GPUContextPool&) = delete;

#ifdef FLEXAIDS_USE_CUDA
    std::mutex cuda_mtx_;
    std::condition_variable cuda_cv_;
    CudaEvalCtx* cuda_ctx_ = nullptr;
    int cuda_natom_ = 0;
    int cuda_ntype_ = 0;
    int cuda_ref_count_ = 0;
    bool cuda_rebuilding_ = false;
#endif

#ifdef FLEXAIDS_USE_METAL
    std::mutex metal_mtx_;
    std::condition_variable metal_cv_;
    MetalEvalCtx* metal_ctx_ = nullptr;
    int metal_natom_   = 0;
    int metal_ntype_   = 0;
    int metal_max_pop_ = 0;
    int metal_ref_count_ = 0;
    bool metal_rebuilding_ = false;
#endif

#ifdef FLEXAIDS_USE_ROCM
    std::mutex rocm_mtx_;
    std::condition_variable rocm_cv_;
    bool rocm_initialised_ = false;
    int  rocm_device_id_   = 0;
    int  rocm_max_pop_     = 0;
    int  rocm_ref_count_   = 0;
    bool rocm_rebuilding_  = false;
#endif
};
