// TurboQuant.metal — Metal compute kernels for TurboQuant batch quantize/dequantize
// Apple Silicon GPU acceleration for FlexAIDdS vector quantization.
//
// Architecture:
//   - One threadgroup per vector (threads_per_threadgroup = d, up to 256)
//   - Threadgroup memory for codebook boundaries/centroids + scratch
//   - SIMD group (warp-equivalent) reductions for L2 norm
//   - Coalesced reads via threadgroup-stride access patterns
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.

#include <metal_stdlib>
using namespace metal;

// ─── GPU-side parameter structs ─────────────────────────────────────────────

struct QuantizeParams {
    uint N;               // number of vectors
    uint d;               // dimension per vector
    uint bit_width;       // quantization bits
    uint num_boundaries;  // 2^bit_width - 1
    uint num_centroids;   // 2^bit_width
};

// ─── SIMD group reduction ───────────────────────────────────────────────────

inline float simd_sum(float val) {
    // Metal SIMD group operations for warp-level reduction
    return simd_sum(val);
}

// =============================================================================
// QUANTIZE KERNEL
// =============================================================================
// Each threadgroup processes one input vector.
// Thread j handles coordinate j.
//
// Pipeline:
//   1. Compute y_j = dot(Pi[j,:], x) via sequential accumulation
//   2. Compute L2 norm of x via SIMD reduction
//   3. Find nearest centroid via boundary scan
//   4. Write index to output

kernel void turboquant_quantize(
    device const float*  Pi          [[buffer(0)]],  // [d × d] rotation, row-major
    device const float*  input       [[buffer(1)]],  // [N × d] input vectors
    device uint8_t*      indices     [[buffer(2)]],  // [N × d] output indices
    device float*        norms       [[buffer(3)]],  // [N] output L2 norms
    device const float*  boundaries  [[buffer(4)]],  // [num_boundaries]
    constant QuantizeParams& params  [[buffer(5)]],
    uint2 tg_pos     [[threadgroup_position_in_grid]],
    uint  tid        [[thread_index_in_threadgroup]],
    uint  tg_size    [[threads_per_threadgroup]],
    uint  simd_lane  [[thread_index_in_simdgroup]],
    uint  simd_id    [[simdgroup_index_in_threadgroup]])
{
    const uint vec_id = tg_pos.x;
    const uint j      = tid;
    const uint d      = params.d;
    const uint num_bd = params.num_boundaries;

    if (vec_id >= params.N || j >= d) return;

    // ── Threadgroup memory for boundaries + warp scratch ────────────────
    threadgroup float tg_bounds[256];   // max 255 boundaries (b=8)
    threadgroup float tg_warp[8];       // partial sums from SIMD groups

    // Load boundaries collaboratively
    if (j < num_bd)
        tg_bounds[j] = boundaries[j];
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // ── Step 1: Compute y_j = dot(Pi[j,:], x) ──────────────────────────
    device const float* x_vec  = input + vec_id * d;
    device const float* pi_row = Pi + j * d;

    float y_j = 0.0f;
    for (uint c = 0; c < d; ++c)
        y_j += pi_row[c] * x_vec[c];

    // ── Step 2: Compute L2 norm of x ────────────────────────────────────
    float x_j = x_vec[j];
    float x_sq = x_j * x_j;

    // SIMD group partial sum
    float simd_partial = simd_sum(x_sq);

    // Write per-SIMD-group partial to threadgroup memory
    if (simd_lane == 0)
        tg_warp[simd_id] = simd_partial;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // First SIMD group finalises the reduction
    if (simd_id == 0) {
        uint n_simd_groups = (tg_size + 31) / 32;
        float total = (simd_lane < n_simd_groups) ? tg_warp[simd_lane] : 0.0f;
        total = simd_sum(total);
        if (simd_lane == 0)
            norms[vec_id] = sqrt(total);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // ── Step 3: Scalar quantize y_j → nearest centroid index ────────────
    int idx = 0;
    for (uint bi = 0; bi < num_bd; ++bi) {
        if (y_j >= tg_bounds[bi]) idx = bi + 1;
        else break;
    }

    indices[vec_id * d + j] = static_cast<uint8_t>(idx);
}

// =============================================================================
// DEQUANTIZE KERNEL
// =============================================================================
// Each threadgroup processes one quantized vector.
// Thread j handles coordinate j.
//
// Pipeline:
//   1. Look up centroid from codebook for each coordinate
//   2. Compute x̃_j = dot(PiT[j,:], centroids_vec) (inverse rotation)

kernel void turboquant_dequantize(
    device const float*   PiT        [[buffer(0)]],  // [d × d] transpose rotation
    device const uint8_t* indices    [[buffer(1)]],  // [N × d] centroid indices
    device const float*   centroids  [[buffer(2)]],  // [k] codebook centroids
    device float*         output     [[buffer(3)]],  // [N × d] output vectors
    constant QuantizeParams& params  [[buffer(4)]],
    uint2 tg_pos     [[threadgroup_position_in_grid]],
    uint  tid        [[thread_index_in_threadgroup]])
{
    const uint vec_id = tg_pos.x;
    const uint j      = tid;
    const uint d      = params.d;
    const uint k      = params.num_centroids;

    if (vec_id >= params.N || j >= d) return;

    // ── Threadgroup memory for centroids + reconstructed y_hat ──────────
    threadgroup float tg_centroids[16];  // max 16 centroids (b=4)
    threadgroup float tg_y_hat[256];     // max 256-dim

    // Load centroids collaboratively
    if (j < k)
        tg_centroids[j] = centroids[j];
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Look up centroid for this coordinate
    uint8_t q_idx = indices[vec_id * d + j];
    tg_y_hat[j] = tg_centroids[q_idx];
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // ── Inverse rotation: x̃_j = dot(PiT[j,:], y_hat) ──────────────────
    device const float* pit_row = PiT + j * d;

    float x_hat_j = 0.0f;
    for (uint c = 0; c < d; ++c)
        x_hat_j += pit_row[c] * tg_y_hat[c];

    output[vec_id * d + j] = x_hat_j;
}
