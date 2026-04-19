// gpu_batch_distance_hip.hip — ROCm/HIP batch pairwise distance kernel
//
// Direct port of gpu_batch_distance.cu for AMD ROCm/HIP runtime.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#ifdef FLEXAIDS_USE_ROCM

#include "gpu_batch_distance.h"
#include "gpu_buffer.h"
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cmath>
#include <chrono>

#define GPU_ROCM_CHECK(call) do {                                                \
    hipError_t _e = (call);                                                      \
    if (_e != hipSuccess) {                                                      \
        fprintf(stderr, "[gpu_batch_distance_hip] ROCm error: %s at %s:%d\n",   \
                hipGetErrorString(_e), __FILE__, __LINE__);                      \
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

BatchDistanceResult rocm_batch_distance(const float* points, int N, int D) {
    BatchDistanceResult result;
    if (N <= 0 || D <= 0 || !points) return result;

    auto t0 = std::chrono::high_resolution_clock::now();

    int totalPairs = N * N;

    GPUBuffer<float> d_points_gpu(N * D, GPUBackend::ROCm);
    GPUBuffer<float> d_dist_gpu(totalPairs, GPUBackend::ROCm);

    d_points_gpu.upload(points, N * D);

    int blockSize = 256;
    int gridSize = (totalPairs + blockSize - 1) / blockSize;

    hipLaunchKernelGGL(batchDistanceKernel,
                       dim3(gridSize), dim3(blockSize), 0, 0,
                       d_points_gpu.data(), d_dist_gpu.data(), N, D);

    GPU_ROCM_CHECK(hipGetLastError());
    GPU_ROCM_CHECK(hipDeviceSynchronize());

    result.distances.resize(totalPairs);
    d_dist_gpu.download(result.distances.data(), totalPairs);

    auto t1 = std::chrono::high_resolution_clock::now();
    result.N = N;
    result.D = D;
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    result.backend_used = 2;

    return result;
}

}  // namespace gpu_batch

#endif // FLEXAIDS_USE_ROCM
