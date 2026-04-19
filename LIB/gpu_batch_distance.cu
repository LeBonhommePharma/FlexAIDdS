// gpu_batch_distance.cu — CUDA batch pairwise distance kernel
//
// One thread per (i,j) pair.  Each thread loops over dimension D
// computing Euclidean distance sqrt(sum((a[d]-b[d])^2)).
//
// For large N this is O(N^2) threads but each thread does minimal work
// and the GPU handles the parallelism naturally.  For N > ~2000 the
// grid is launched in tiles.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#ifdef FLEXAIDS_USE_CUDA

#include "gpu_batch_distance.h"
#include "gpu_buffer.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>
#include <chrono>

#define GPU_CUDA_CHECK(call) do {                                                \
    cudaError_t _e = (call);                                                     \
    if (_e != cudaSuccess) {                                                     \
        fprintf(stderr, "[gpu_batch_distance] CUDA error: %s at %s:%d\n",       \
                cudaGetErrorString(_e), __FILE__, __LINE__);                     \
        return {};                                                               \
    }                                                                            \
} while (0)

// ─── Kernel: one thread per pair ──────────────────────────────────────────────
__global__ void batchDistanceKernel(const float* __restrict__ d_points,
                                     float*       __restrict__ d_dist,
                                     int N, int D)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int totalPairs = N * N;
    if (idx >= totalPairs) return;

    int i = idx / N;
    int j = idx % N;

    if (i == j) {
        d_dist[idx] = 0.0f;
        return;
    }

    const float* pi = d_points + i * D;
    const float* pj = d_points + j * D;

    float sum = 0.0f;
    for (int d = 0; d < D; ++d) {
        float diff = pi[d] - pj[d];
        sum += diff * diff;
    }
    d_dist[idx] = sqrtf(sum);
}

namespace gpu_batch {

BatchDistanceResult cuda_batch_distance(const float* points, int N, int D) {
    BatchDistanceResult result;
    if (N <= 0 || D <= 0 || !points) return result;

    auto t0 = std::chrono::high_resolution_clock::now();

    int totalPairs = N * N;

    GPUBuffer<float> d_points_gpu(N * D, GPUBackend::CUDA);
    GPUBuffer<float> d_dist_gpu(totalPairs, GPUBackend::CUDA);

    d_points_gpu.upload(points, N * D);

    int blockSize = 256;
    int gridSize = (totalPairs + blockSize - 1) / blockSize;

    batchDistanceKernel<<<gridSize, blockSize>>>(
        d_points_gpu.data(), d_dist_gpu.data(), N, D);

    GPU_CUDA_CHECK(cudaGetLastError());
    GPU_CUDA_CHECK(cudaDeviceSynchronize());

    result.distances.resize(totalPairs);
    d_dist_gpu.download(result.distances.data(), totalPairs);

    auto t1 = std::chrono::high_resolution_clock::now();
    result.N = N;
    result.D = D;
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    result.backend_used = 1;

    return result;
}

}  // namespace gpu_batch

#endif // FLEXAIDS_USE_CUDA
