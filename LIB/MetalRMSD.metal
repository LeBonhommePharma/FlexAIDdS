// MetalRMSD.metal — GPU pairwise RMSD kernel for FOPTICS clustering
//
// Computes N*(N-1)/2 pairwise squared-RMSD values in a single dispatch.
// Each GPU thread handles one unique (i,j) pair from the upper triangle.
// Output is a flat upper-triangular distance matrix suitable for OPTICS.
//
// Performance design:
//   – One thread per pair → embarrassingly parallel, no atomics needed
//   – Coalesced reads via contiguous flat buffer layout [N × 3M]
//   – Optimal for M3 Pro 18-core GPU: saturates with N > 256
//
// Tiling strategy for N > 4096:
//   When N is large, the full N*(N-1)/2 output may exceed GPU memory.
//   The host bridge tiles the computation by dispatching batches of rows,
//   each producing a partial output.  The tile size is chosen so that
//   (tile_rows * N * sizeof(float)) + (N * 3M * sizeof(float)) fits in
//   the recommended working set (~1 GB for Apple Silicon unified memory).
//   Typical tile: 512 rows × N columns × 4 bytes.
//
// Apache-2.0 (C) 2026 Le Bonhomme Pharma

#include <metal_stdlib>
using namespace metal;

// ===========================================================================
// PAIRWISE SQUARED-RMSD KERNEL
// ===========================================================================
// Computes RMSD^2 * n_atoms for each unique pair (i,j), i < j.
// RMSD^2 = (1/n_atoms) * sum_k sum_dim (coords[i][k*3+dim] - coords[j][k*3+dim])^2
// We store (RMSD^2 * n_atoms) = sum of squared differences, so the host can
// divide by n_atoms to get RMSD^2, or take sqrt(RMSD^2) for actual RMSD.
// This avoids redundant division inside the kernel.
//
// Output layout: flat upper-triangular, row-major.
//   For pair (i,j) with i < j, the linear index is: i*N - i*(i+1)/2 + (j-i-1)
//   Total output size: N*(N-1)/2 floats.

kernel void pairwise_rmsd_squared(
    device const float* coords   [[buffer(0)]],   // flat [N * 3*M]: M atoms × 3 coords per conformation
    device float*       output   [[buffer(1)]],   // flat [N*(N-1)/2]: upper-triangular squared-distance
    constant uint&      n_conf   [[buffer(2)]],   // N: number of conformations
    constant uint&      n_atoms  [[buffer(3)]],   // M: number of atoms (coords per conf = 3*M)
    uint                gid      [[thread_position_in_grid]])
{
    // Total number of unique pairs
    ulong total_pairs = static_cast<ulong>(n_conf) * (n_conf - 1) / 2;

    // Bounds check
    if (static_cast<ulong>(gid) >= total_pairs) return;

    // ── Decode linear pair index → (i, j) with i < j ──────────────────────
    // Inverse of the upper-triangular indexing:
    //   linear = i*N - i*(i+1)/2 + (j - i - 1)
    // We solve for i first using the quadratic formula:
    //   i^2 - (2N-1)*i + 2*linear = 0  →  i = floor((2N-1 - sqrt(D)) / 2)
    // where D = (2N-1)^2 - 8*linear.
    // Then j = linear - i*N + i*(i+1)/2 + i + 1.

    ulong N = static_cast<ulong>(n_conf);
    ulong lin = static_cast<ulong>(gid);

    // Solve for row index i
    // D = (2N - 1)^2 - 8 * lin
    ulong two_N_minus_1 = 2 * N - 1;
    ulong D = two_N_minus_1 * two_N_minus_1 - 8 * lin;

    // fast inverse square root not needed; Metal has precise sqrt for ulong
    ulong sqrt_D = static_cast<ulong>(metal::sqrt(static_cast<float>(D)));
    // Adjust for floating-point rounding
    while ((sqrt_D + 1) * (sqrt_D + 1) <= D) sqrt_D++;
    while (sqrt_D * sqrt_D > D) sqrt_D--;

    ulong i = (two_N_minus_1 - sqrt_D) / 2;
    // Clamp i to valid range
    if (i >= N) i = N - 1;

    // Compute j from i and linear index
    // j = lin - (i * N - i * (i + 1) / 2) + i + 1
    ulong j = lin - (i * N - i * (i + 1) / 2) + i + 1;

    // Safety: ensure i < j < N
    if (j >= N) return;

    // ── Compute sum of squared coordinate differences ──────────────────────
    ulong offset_i = i * 3 * static_cast<ulong>(n_atoms);
    ulong offset_j = j * 3 * static_cast<ulong>(n_atoms);

    float sum_sq = 0.0f;
    ulong ncoords = 3 * static_cast<ulong>(n_atoms);
    for (ulong k = 0; k < ncoords; ++k) {
        float diff = coords[offset_i + k] - coords[offset_j + k];
        sum_sq += diff * diff;
    }

    // Store raw sum of squared differences (host divides by n_atoms for RMSD^2)
    output[lin] = sum_sq;
}

// ===========================================================================
// TILED PAIRWISE SQUARED-RMSD KERNEL
// ===========================================================================
// For large N, computes RMSD for a contiguous range of row indices
// [row_start, row_end).  Each thread handles one (i,j) pair where
// row_start <= i < row_end and i < j < N.
// Output is written starting at out_offset so host can stitch tiles.

constant uint kTileThreadsPerGroup [[function_constant(0)]]; // not used as constant, just doc

kernel void pairwise_rmsd_squared_tiled(
    device const float* coords     [[buffer(0)]],
    device float*       output     [[buffer(1)]],
    constant uint&      n_conf     [[buffer(2)]],
    constant uint&      n_atoms    [[buffer(3)]],
    constant uint&      row_start  [[buffer(4)]],   // first row index in this tile
    constant uint&      row_end    [[buffer(5)]],   // one past last row index
    constant ulong&     out_offset [[buffer(6)]],   // byte offset into output buffer
    uint                gid        [[thread_position_in_grid]])
{
    ulong N = static_cast<ulong>(n_conf);

    // Each thread handles one pair (i,j) where row_start <= i < row_end, i < j < N
    // Total pairs in this tile: sum_{i=row_start}^{row_end-1} (N - 1 - i)
    ulong total_tile_pairs = 0;
    for (ulong r = static_cast<ulong>(row_start); r < static_cast<ulong>(row_end); ++r) {
        total_tile_pairs += (N - 1 - r);
    }

    if (static_cast<ulong>(gid) >= total_tile_pairs) return;

    // Decode gid → (i, j) within tile
    ulong remaining = static_cast<ulong>(gid);
    ulong i = static_cast<ulong>(row_start);
    while (i < static_cast<ulong>(row_end)) {
        ulong pairs_in_row = N - 1 - i;
        if (remaining < pairs_in_row) {
            break;
        }
        remaining -= pairs_in_row;
        i++;
    }
    ulong j = i + 1 + remaining;

    if (i >= static_cast<ulong>(row_end) || j >= N) return;

    // Compute sum of squared coordinate differences
    ulong offset_i = i * 3 * static_cast<ulong>(n_atoms);
    ulong offset_j = j * 3 * static_cast<ulong>(n_atoms);

    float sum_sq = 0.0f;
    ulong ncoords = 3 * static_cast<ulong>(n_atoms);
    for (ulong k = 0; k < ncoords; ++k) {
        float diff = coords[offset_i + k] - coords[offset_j + k];
        sum_sq += diff * diff;
    }

    // Write to output at the tile's linear position
    output[gid] = sum_sq;
}
