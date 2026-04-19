// gpu_batch_distance.cpp — CPU SIMD fallback for batch pairwise distance
//
// Uses flexaids::sum_sq_distances_f() which auto-dispatches to
// AVX-512 / AVX2 / scalar.  OpenMP parallelises the outer loop.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#include "gpu_batch_distance.h"
#include "simd_distance.h"
#include <cmath>
#include <chrono>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace gpu_batch {

BatchDistanceResult cpu_batch_distance(const float* points, int N, int D) {
    BatchDistanceResult result;
    if (N <= 0 || D <= 0 || !points) return result;

    auto t0 = std::chrono::high_resolution_clock::now();

    result.distances.resize(static_cast<size_t>(N) * N);
    result.N = N;
    result.D = D;

    #pragma omp parallel for schedule(dynamic) if(N > 32)
    for (int i = 0; i < N; ++i) {
        const float* pi = points + i * D;
        // Diagonal
        result.distances[static_cast<size_t>(i) * N + i] = 0.0f;

        // Upper triangle (j > i), mirror to lower
        for (int j = i + 1; j < N; ++j) {
            const float* pj = points + j * D;
            float d2 = flexaids::sum_sq_distances_f(pi, pj, D);
            float dist = std::sqrt(d2);
            result.distances[static_cast<size_t>(i) * N + j] = dist;
            result.distances[static_cast<size_t>(j) * N + i] = dist;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    result.backend_used = 0;

    return result;
}

}  // namespace gpu_batch
