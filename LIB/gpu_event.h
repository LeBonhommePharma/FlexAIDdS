// gpu_event.h — Lightweight GPU completion event for async dispatch
//
// Wraps cudaEvent_t / hipEvent_t / MTLSharedEvent into a uniform RAII
// token that callers can query, block on, or attach a callback to.
//
// Usage:
//   GPUEvent ev(Backend::CUDA);
//   ev.record(stream);          // stamp the event on a stream
//   // ... do other work ...
//   if (!ev.is_complete()) ev.wait();   // block until GPU is done
//
// For CPU-only work, GPUEvent acts as an immediately-complete token.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <atomic>
#include <functional>
#include <cstdint>

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#include <immintrin.h>
#endif

#ifdef FLEXAIDS_USE_CUDA
#include <cuda_runtime.h>
#endif

#ifdef FLEXAIDS_USE_ROCM
#include <hip/hip_runtime.h>
#endif

#ifdef __APPLE__
#  ifdef __OBJC__
#    import <Metal/Metal.h>
#  else
     // Forward-declare opaque pointer for non-ObjC++ translation units
     using MTLSharedEventPtr = void*;
#  endif
#endif

namespace flexaids {

// Backend tag (matches the plan's Backend enum values)
enum class EventBackend : uint8_t {
    CPU   = 0,
    CUDA  = 1,
    ROCM  = 2,
    METAL = 3
};

// ─── GPUEvent ──────────────────────────────────────────────────────────────

class GPUEvent {
public:
    // Construct an event for the given backend.
    // CPU events are born "complete" (no GPU work to wait on).
    explicit GPUEvent(EventBackend backend = EventBackend::CPU)
        : backend_(backend)
    {
        switch (backend_) {
#ifdef FLEXAIDS_USE_CUDA
        case EventBackend::CUDA: {
            cudaError_t err = cudaEventCreateWithFlags(&cuda_event_, cudaEventDisableTiming);
            if (err != cudaSuccess) {
                backend_ = EventBackend::CPU;
                cpu_complete_.store(true, std::memory_order_release);
            }
            break;
        }
#endif
#ifdef FLEXAIDS_USE_ROCM
        case EventBackend::ROCM: {
            hipError_t err = hipEventCreateWithFlags(&hip_event_, hipEventDisableTiming);
            if (err != hipSuccess) {
                backend_ = EventBackend::CPU;
                cpu_complete_.store(true, std::memory_order_release);
            }
            break;
        }
#endif
        case EventBackend::METAL:
            // Metal events are signalled via command-buffer completion handlers
            // and tracked by the atomic flag below.
            break;
        case EventBackend::CPU:
        default:
            cpu_complete_.store(true, std::memory_order_release);
            break;
        }
    }

    ~GPUEvent() { destroy(); }

    // Move-only
    GPUEvent(GPUEvent&& o) noexcept { swap(o); }
    GPUEvent& operator=(GPUEvent&& o) noexcept {
        if (this != &o) { destroy(); swap(o); }
        return *this;
    }
    GPUEvent(const GPUEvent&) = delete;
    GPUEvent& operator=(const GPUEvent&) = delete;

    // ── Record ──────────────────────────────────────────────────────────

#ifdef FLEXAIDS_USE_CUDA
    // Record the event on a CUDA stream.
    void record(cudaStream_t stream) {
        cudaEventRecord(cuda_event_, stream);
    }
#endif

#ifdef FLEXAIDS_USE_ROCM
    // Record the event on a HIP stream.
    void record(hipStream_t stream) {
        hipEventRecord(hip_event_, stream);
    }
#endif

    // Mark a CPU/Metal event as complete (called from completion handler
    // or from the submitting thread once CPU work finishes).
    void signal() noexcept {
        cpu_complete_.store(true, std::memory_order_release);
    }

    // ── Query ───────────────────────────────────────────────────────────

    // Non-blocking: has the GPU/CPU work finished?
    bool is_complete() const noexcept {
        switch (backend_) {
#ifdef FLEXAIDS_USE_CUDA
        case EventBackend::CUDA:
            return cudaEventQuery(cuda_event_) == cudaSuccess;
#endif
#ifdef FLEXAIDS_USE_ROCM
        case EventBackend::ROCM:
            return hipEventQuery(hip_event_) == hipSuccess;
#endif
        case EventBackend::METAL:
        case EventBackend::CPU:
        default:
            return cpu_complete_.load(std::memory_order_acquire);
        }
    }

    // ── Wait ────────────────────────────────────────────────────────────

    // Blocking wait until the event completes.
    void wait() const {
        switch (backend_) {
#ifdef FLEXAIDS_USE_CUDA
        case EventBackend::CUDA:
            cudaEventSynchronize(cuda_event_);
            return;
#endif
#ifdef FLEXAIDS_USE_ROCM
        case EventBackend::ROCM:
            hipEventSynchronize(hip_event_);
            return;
#endif
        case EventBackend::METAL:
        case EventBackend::CPU:
        default:
            // Spin-wait (Metal/CPU events are typically fast)
            while (!cpu_complete_.load(std::memory_order_acquire)) {
                // Yield to avoid burning a core
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
                _mm_pause();
#elif defined(__x86_64__) || defined(_M_X64)
                __builtin_ia32_pause();
#elif defined(__aarch64__)
                asm volatile("yield");
#endif
            }
            return;
        }
    }

    // NOTE: elapsed_ms() is NOT provided because events are created with
    // cudaEventDisableTiming for lower overhead.  Use wall-clock timing
    // (e.g. std::chrono) for profiling instead.

    EventBackend backend() const noexcept { return backend_; }

private:
    EventBackend backend_ = EventBackend::CPU;
    std::atomic<bool> cpu_complete_{false};

#ifdef FLEXAIDS_USE_CUDA
    cudaEvent_t cuda_event_{};
#endif
#ifdef FLEXAIDS_USE_ROCM
    hipEvent_t hip_event_{};
#endif

    void destroy() noexcept {
        switch (backend_) {
#ifdef FLEXAIDS_USE_CUDA
        case EventBackend::CUDA:
            if (cuda_event_) { cudaEventDestroy(cuda_event_); cuda_event_ = {}; }
            break;
#endif
#ifdef FLEXAIDS_USE_ROCM
        case EventBackend::ROCM:
            if (hip_event_) { hipEventDestroy(hip_event_); hip_event_ = {}; }
            break;
#endif
        default:
            break;
        }
    }

    void swap(GPUEvent& o) noexcept {
        std::swap(backend_, o.backend_);
        bool tmp = cpu_complete_.load();
        cpu_complete_.store(o.cpu_complete_.load());
        o.cpu_complete_.store(tmp);
#ifdef FLEXAIDS_USE_CUDA
        std::swap(cuda_event_, o.cuda_event_);
#endif
#ifdef FLEXAIDS_USE_ROCM
        std::swap(hip_event_, o.hip_event_);
#endif
    }
};

}  // namespace flexaids
