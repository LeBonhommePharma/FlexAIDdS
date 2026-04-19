// gpu_batch_distance.h — GPU batch pairwise distance computation
//
// Provides hardware-dispatched batch distance kernels for FastOPTICS
// and other clustering algorithms.  Supports CUDA, ROCm/HIP, Metal,
// and CPU SIMD (AVX-512/AVX2/scalar) backends.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab
#pragma once

#include <vector>
#include <cstdint>

namespace gpu_batch {

/// Result of a batch pairwise distance computation.
/// distances[i * N + j] = Euclidean distance between points i and j.
struct BatchDistanceResult {
    std::vector<float> distances;   // N x N row-major
    int    N = 0;
    int    D = 0;                   // dimensionality
    double elapsed_ms = 0.0;
    int    backend_used = 0;        // 0=CPU, 1=CUDA, 2=ROCm, 3=Metal
};

/// CUDA batch pairwise distance (N x N matrix).
/// Returns result with backend_used=1 on success, or falls back to CPU.
BatchDistanceResult cuda_batch_distance(const float* points, int N, int D);

/// ROCm/HIP batch pairwise distance (N x N matrix).
/// Returns result with backend_used=2 on success, or falls back to CPU.
BatchDistanceResult rocm_batch_distance(const float* points, int N, int D);

/// Metal batch pairwise distance (N x N matrix, macOS only).
/// Returns result with backend_used=3 on success, or falls back to CPU.
BatchDistanceResult metal_batch_distance(const float* points, int N, int D);

/// CPU SIMD fallback using flexaids::sum_sq_distances_f().
/// Dispatches to AVX-512 / AVX2 / scalar automatically.
BatchDistanceResult cpu_batch_distance(const float* points, int N, int D);

/// Auto-select best available backend and compute.
/// Priority: CUDA > ROCm > Metal > CPU-SIMD.
BatchDistanceResult auto_batch_distance(const float* points, int N, int D);

}  // namespace gpu_batch
