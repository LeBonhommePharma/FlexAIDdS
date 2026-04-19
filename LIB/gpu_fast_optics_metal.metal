// gpu_fast_optics_metal.metal — Metal compute shader for kNN search
//
// Kernel: each threadgroup handles one query point.
// Threads cooperatively compute distances and maintain a shared top-k list.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#include <metal_stdlib>
using namespace metal;

// Kernel: pairwise Euclidean distance + k-nearest neighbour
// Each threadgroup processes one query point (group_id = query index).
// 256 threads cooperatively scan all N points, keeping the k nearest.
kernel void gpuFastOPTICSKernel(
    device const float* d_points   [[buffer(0)]],
    device int*         d_knn_idx  [[buffer(1)]],
    device float*       d_knn_dist [[buffer(2)]],
    constant int&       N          [[buffer(3)]],
    constant int&       D          [[buffer(4)]],
    constant int&       k          [[buffer(5)]],
    uint                gid        [[thread_index_in_threadgroup]],
    uint                group_id   [[threadgroup_position_in_grid]])
{
    // Shared memory for query point (max 256 dimensions)
    threadgroup float shared_q[256];

    const uint qid = group_id;
    if (qid >= (uint)N) return;

    device const float* q = d_points + qid * D;

    // Cooperative load of query point into shared memory
    for (uint i = gid; i < (uint)D; i += 256) {
        shared_q[i] = q[i];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Thread-local top-k storage (LOCAL_K_MAX = 20)
    constexpr int LOCAL_K_MAX = 20;
    const int actual_k = (k < LOCAL_K_MAX) ? k : LOCAL_K_MAX;

    float local_dist[LOCAL_K_MAX];
    int   local_idx[LOCAL_K_MAX];
    int   local_count = 0;

    for (int i = 0; i < actual_k; ++i) {
        local_dist[i] = 1e30f;
        local_idx[i]  = -1;
    }

    // Strided scan over all candidate points
    for (uint pid = gid; pid < (uint)N; pid += 256) {
        if (pid == qid) continue;

        device const float* p = d_points + pid * D;
        float dist2 = 0.0f;
        for (int d = 0; d < D; ++d) {
            float diff = shared_q[d] - p[d];
            dist2 += diff * diff;
        }
        float dist = metal::sqrt(dist2);

        // Insert into local top-k if closer than current worst
        if (local_count < actual_k) {
            local_dist[local_count] = dist;
            local_idx[local_count]  = (int)pid;
            local_count++;
        } else {
            int worst = 0;
            for (int j = 1; j < actual_k; ++j) {
                if (local_dist[j] > local_dist[worst]) worst = j;
            }
            if (dist < local_dist[worst]) {
                local_dist[worst] = dist;
                local_idx[worst]  = (int)pid;
            }
        }
    }

    // Shared memory for cross-thread merge
    threadgroup float s_dist[256 * LOCAL_K_MAX];
    threadgroup int   s_idx[256 * LOCAL_K_MAX];
    threadgroup int   s_count[256];

    int base = (int)gid * LOCAL_K_MAX;
    s_count[gid] = (local_count < actual_k) ? local_count : actual_k;
    for (int i = 0; i < actual_k && i < local_count; ++i) {
        s_dist[base + i] = local_dist[i];
        s_idx[base + i]  = local_idx[i];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Thread 0 merges all thread-local lists into final top-k
    if (gid == 0) {
        float final_dist[LOCAL_K_MAX];
        int   final_idx[LOCAL_K_MAX];
        int   final_count = 0;

        for (int i = 0; i < actual_k; ++i) {
            final_dist[i] = 1e30f;
            final_idx[i]  = -1;
        }

        for (uint t = 0; t < 256; ++t) {
            int tbase = (int)t * LOCAL_K_MAX;
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

        // Insertion sort by distance
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

        // Write output
        int out_base = (int)qid * k;
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
