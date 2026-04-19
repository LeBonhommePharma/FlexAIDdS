// gpu_fast_optics_hip.hip — ROCm/HIP GPU kNN for FastOPTICS
//
// When FLEXAIDS_USE_ROCM is defined, this TU provides a GPU kernel that
// accelerates the pairwise distance computation and neighbour identification
// used by the FastOPTICS clustering algorithm.
//
// Direct port of gpu_fast_optics.cu for AMD ROCm/HIP runtime.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#ifdef FLEXAIDS_USE_ROCM

#include "FOPTICS.h"
#include "gpu_buffer.h"
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

#include "flexaid_exception.h"
#include <string>

#define GPU_ROCM_CHECK(call) do {                                                \
    hipError_t _e = (call);                                                      \
    if (_e != hipSuccess) {                                                      \
        throw FlexAIDException(std::string("[gpu_fast_optics_hip] ROCm error at ")\
            + __FILE__ + ":" + std::to_string(__LINE__) + " — "                  \
            + hipGetErrorString(_e));                                            \
    }                                                                            \
} while (0)

static constexpr int BLOCK_SIZE = 256;

// ─── kernel: pairwise Euclidean distance + k-nearest neighbour ────────────────
__global__ void gpuFastOPTICSKernel(const float* __restrict__ d_points,
                                     int*   __restrict__ d_knn_idx,
                                     float* __restrict__ d_knn_dist,
                                     int N, int D, int k)
{
    const int qid = blockIdx.x;
    if (qid >= N) return;

    const float* q = d_points + qid * D;

    extern __shared__ float shared_q[];

    for (int i = threadIdx.x; i < D; i += blockDim.x) {
        shared_q[i] = q[i];
    }
    __syncthreads();

    constexpr int LOCAL_K_MAX = 20;
    float local_dist[LOCAL_K_MAX];
    int   local_idx[LOCAL_K_MAX];
    int   local_count = 0;
    int   actual_k = (k < LOCAL_K_MAX) ? k : LOCAL_K_MAX;

    for (int i = 0; i < actual_k; ++i) {
        local_dist[i] = 1e30f;
        local_idx[i]  = -1;
    }

    for (int pid = threadIdx.x; pid < N; pid += blockDim.x) {
        if (pid == qid) continue;

        const float* p = d_points + pid * D;
        float dist2 = 0.0f;
        for (int d = 0; d < D; ++d) {
            float diff = shared_q[d] - p[d];
            dist2 += diff * diff;
        }
        float dist = sqrtf(dist2);

        if (local_count < actual_k) {
            local_dist[local_count] = dist;
            local_idx[local_count]  = pid;
            local_count++;
        } else {
            int worst = 0;
            for (int j = 1; j < actual_k; ++j) {
                if (local_dist[j] > local_dist[worst]) worst = j;
            }
            if (dist < local_dist[worst]) {
                local_dist[worst] = dist;
                local_idx[worst]  = pid;
            }
        }
    }

    __shared__ float  s_dist[BLOCK_SIZE * 20];
    __shared__ int    s_idx[BLOCK_SIZE * 20];
    __shared__ int    s_count[BLOCK_SIZE];

    int base = threadIdx.x * LOCAL_K_MAX;
    s_count[threadIdx.x] = (local_count < actual_k) ? local_count : actual_k;
    for (int i = 0; i < actual_k && i < local_count; ++i) {
        s_dist[base + i] = local_dist[i];
        s_idx[base + i]  = local_idx[i];
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        float  final_dist[LOCAL_K_MAX];
        int    final_idx[LOCAL_K_MAX];
        int    final_count = 0;

        for (int i = 0; i < actual_k; ++i) {
            final_dist[i] = 1e30f;
            final_idx[i]  = -1;
        }

        for (int t = 0; t < blockDim.x; ++t) {
            int tbase = t * LOCAL_K_MAX;
            int tcount = s_count[t];
            for (int i = 0; i < tcount; ++i) {
                float d = s_dist[tbase + i];
                int   idx = s_idx[tbase + i];
                if (final_count < actual_k) {
                    final_dist[final_count] = d;
                    final_idx[final_count]  = idx;
                    final_count++;
                } else {
                    int worst = 0;
                    for (int j = 1; j < actual_k; ++j) {
                        if (final_dist[j] > final_dist[worst]) worst = j;
                    }
                    if (d < final_dist[worst]) {
                        final_dist[worst] = d;
                        final_idx[worst]  = idx;
                    }
                }
            }
        }

        for (int i = 1; i < actual_k; ++i) {
            float key_d = final_dist[i];
            int   key_i = final_idx[i];
            int j = i - 1;
            while (j >= 0 && final_dist[j] > key_d) {
                final_dist[j + 1] = final_dist[j];
                final_idx[j + 1]  = final_idx[j];
                --j;
            }
            final_dist[j + 1] = key_d;
            final_idx[j + 1]  = key_i;
        }

        int out_base = qid * k;
        for (int i = 0; i < k; ++i) {
            if (i < actual_k) {
                d_knn_idx[out_base + i]  = final_idx[i];
                d_knn_dist[out_base + i] = final_dist[i];
            } else {
                d_knn_idx[out_base + i]  = -1;
                d_knn_dist[out_base + i] = 1e30f;
            }
        }
    }
}

// ─── host wrapper ─────────────────────────────────────────────────────────────
void gpu_rocm_foptics_knn(const std::vector<std::pair<chromosome*, std::vector<float>>>& points,
                           int k, int nDim,
                           std::vector<std::vector<int>>& out_neighbors,
                           std::vector<std::vector<float>>& out_distances)
{
    int N = static_cast<int>(points.size());
    if (N == 0 || k <= 0) return;

    std::vector<float> h_points(N * nDim);
    for (int i = 0; i < N; ++i) {
        const auto& coords = points[i].second;
        for (int d = 0; d < nDim; ++d) {
            h_points[i * nDim + d] = (d < static_cast<int>(coords.size())) ? coords[d] : 0.0f;
        }
    }

    GPUBuffer<float> d_points(N * nDim, GPUBackend::ROCm);
    GPUBuffer<int>   d_knn_idx(N * k, GPUBackend::ROCm);
    GPUBuffer<float> d_knn_dist(N * k, GPUBackend::ROCm);

    d_points.upload(h_points.data(), N * nDim);

    int shared_mem = nDim * sizeof(float);
    hipLaunchKernelGGL(gpuFastOPTICSKernel,
                       dim3(N), dim3(BLOCK_SIZE), shared_mem, 0,
                       d_points.data(), d_knn_idx.data(), d_knn_dist.data(),
                       N, nDim, k);
    GPU_ROCM_CHECK(hipGetLastError());
    GPU_ROCM_CHECK(hipDeviceSynchronize());

    std::vector<int>   h_knn_idx(N * k);
    std::vector<float> h_knn_dist(N * k);
    d_knn_idx.download(h_knn_idx.data(), N * k);
    d_knn_dist.download(h_knn_dist.data(), N * k);

    out_neighbors.resize(N);
    out_distances.resize(N);
    for (int i = 0; i < N; ++i) {
        out_neighbors[i].clear();
        out_distances[i].clear();
        int base = i * k;
        for (int j = 0; j < k; ++j) {
            if (h_knn_idx[base + j] >= 0) {
                out_neighbors[i].push_back(h_knn_idx[base + j]);
                out_distances[i].push_back(h_knn_dist[base + j]);
            }
        }
    }
}

#endif // FLEXAIDS_USE_ROCM
