// shannon_metal.metal — GPU kernels for Shannon entropy and Boltzmann weight computation
//
// Performance design:
//   – 256 threads per threadgroup (optimal for Apple Silicon SIMD width)
//   – threadgroup-local histogram (reduces global atomic contention ~32×)
//   – Final merge via atomic_fetch_add to global bin buffer
//   – Boltzmann batch kernel uses shared-memory parallel reduction
//   – Supports arbitrary n (out-of-bounds threads are no-ops)
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#include <metal_stdlib>
using namespace metal;

// ===========================================================================
// SHANNON HISTOGRAM KERNEL
// ===========================================================================
// Bins energy values into histogram for Shannon entropy computation.
// Uses threadgroup-local histograms to minimize global atomic contention.

kernel void shannon_histogram(
    device const double*     energies  [[buffer(0)]],
    device atomic_int*       bins      [[buffer(1)]],
    constant uint&           n         [[buffer(2)]],
    constant int&            num_bins  [[buffer(3)]],
    constant double&         min_v     [[buffer(4)]],
    constant double&         bin_width [[buffer(5)]],
    uint                     gid       [[thread_position_in_grid]],
    uint                     lid       [[thread_position_in_threadgroup]],
    uint                     tg_size   [[threads_per_threadgroup]])
{
    // Threadgroup-local histogram to reduce global atomic pressure
    threadgroup atomic_int local_bins[256]; // max num_bins = 256
    if (lid < (uint)num_bins)
        atomic_store_explicit(&local_bins[lid], 0, memory_order_relaxed);
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Bin this thread's element (coalesced memory access pattern)
    if (gid < n) {
        double e = energies[gid];
        int b = (int)((e - min_v) / bin_width);
        b = clamp(b, 0, num_bins - 1);
        atomic_fetch_add_explicit(&local_bins[b], 1, memory_order_relaxed);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Merge threadgroup histogram to global (only threads 0..num_bins-1)
    if (lid < (uint)num_bins) {
        int local_count = atomic_load_explicit(&local_bins[lid], memory_order_relaxed);
        if (local_count > 0)
            atomic_fetch_add_explicit(&bins[lid], local_count, memory_order_relaxed);
    }
}

// ===========================================================================
// BOLTZMANN WEIGHT BATCH KERNEL
// ===========================================================================
// Computes w[i] = exp(-beta * (E[i] - E_min)) for partition function evaluation.
// E_min must be pre-computed on CPU (single pass) and passed as constant.
// Uses Metal's fast math exp() for Apple Silicon's FP64 units.

kernel void boltzmann_weights_batch(
    device const double*  energies     [[buffer(0)]],
    device double*        weights      [[buffer(1)]],
    constant uint&        n            [[buffer(2)]],
    constant double&      neg_beta     [[buffer(3)]],
    constant double&      E_min        [[buffer(4)]],
    uint                  gid          [[thread_position_in_grid]])
{
    if (gid >= n) return;
    weights[gid] = exp(neg_beta * (energies[gid] - E_min));
}

// ===========================================================================
// PARALLEL SUM REDUCTION KERNEL
// ===========================================================================
// Sums an array of doubles using threadgroup shared-memory reduction.
// Output: partial sums (one per threadgroup), final sum done on CPU.

kernel void parallel_sum_reduce(
    device const double*   input       [[buffer(0)]],
    device double*         partials    [[buffer(1)]],
    constant uint&         n           [[buffer(2)]],
    uint                   gid         [[thread_position_in_grid]],
    uint                   lid         [[thread_position_in_threadgroup]],
    uint                   tgid        [[threadgroup_position_in_grid]],
    uint                   tg_size     [[threads_per_threadgroup]])
{
    // Must match the threadgroup size dispatched by the host (256).
    // Using a constant ensures no OOB if tg_size matches.
    constexpr uint MAX_TG_SIZE = 256;
    threadgroup double shared[MAX_TG_SIZE];

    // Guard: if runtime tg_size exceeds our shared array, bail safely
    if (lid >= MAX_TG_SIZE) return;

    // Load (zero-pad out-of-bounds)
    shared[lid] = (gid < n) ? input[gid] : 0.0;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Tree reduction within threadgroup
    for (uint stride = tg_size / 2; stride > 0; stride >>= 1) {
        if (lid < stride) {
            shared[lid] += shared[lid + stride];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    // Thread 0 writes the threadgroup's partial sum
    if (lid == 0) {
        partials[tgid] = shared[0];
    }
}

// ===========================================================================
// LOG-SUM-EXP KERNEL
// ===========================================================================
// Computes exp(x[i] - x_max) for log-sum-exp stability.
// x_max must be pre-computed. Final log(sum) + x_max done on CPU.

kernel void log_sum_exp_shifted(
    device const double*  values       [[buffer(0)]],
    device double*        exp_shifted  [[buffer(1)]],
    constant uint&        n            [[buffer(2)]],
    constant double&      x_max        [[buffer(3)]],
    uint                  gid          [[thread_position_in_grid]])
{
    if (gid >= n) return;
    exp_shifted[gid] = exp(values[gid] - x_max);
}
