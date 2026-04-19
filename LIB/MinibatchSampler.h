// MinibatchSampler.h — Farthest-point sampling pre-filter for clustering
//
// Reduces the chromosome population from ~500K to ~5K diverse representatives
// before FO/DP/OPTICS clustering, yielding a ~100x speedup in the O(N^2)
// distance phase while preserving cluster structure.
//
// Key features:
//   * farthest_point_sample() — O(k*n) deterministic diversity maximisation
//   * quality_weighted_sample() — blends distance with energy (alpha param)
//   * SIMD-accelerated batch distance via simd_distance.h (AVX2/512/SSE4.2)
//   * OpenMP parallel distance-to-nearest-selected accumulation
//   * Thread-safe (all const operations, no mutable state)
//
// Integration point: called from FastOPTICS_cluster() and DensityPeak_cluster()
// before the BindingPopulation/FOPTICS constructors.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#ifndef MINIBATCHSAMPLER_H
#define MINIBATCHSAMPLER_H

#include "simd_distance.h"

#include <vector>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cassert>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace minibatch {

// ─── Configuration ───────────────────────────────────────────────────────────

struct SamplerConfig {
    int   target_k      = 5000;    // number of representatives to select
    float alpha         = 0.7f;    // diversity weight (1.0 = pure FPS, 0.0 = pure energy)
    bool  enable        = true;    // master switch (disabled when n <= target_k)
    bool  verbose       = true;    // print summary to stdout
};

// ─── Result ──────────────────────────────────────────────────────────────────

struct SampleResult {
    std::vector<int>  selected_indices;  // indices into the original chrom array
    std::vector<int>  assignment;        // assignment[i] = index of nearest selected for chrom i
    double            elapsed_ms;        // wall-clock time
    int               n_original;        // input population size
    int               n_selected;        // output representative count
};

// ─── Coordinate cache ────────────────────────────────────────────────────────
//
// Pre-compute Cartesian coordinates for all chromosomes once, then reuse
// for all FPS distance queries.  This mirrors the coord_cache pattern in
// cluster.cpp but is packaged as a self-contained struct.

struct CoordCache {
    std::vector<float> data;     // (n_chrom * stride) contiguous floats
    int                n_chrom;
    int                stride;   // nAtoms * 3

    const float* operator[](int i) const noexcept {
        return data.data() + static_cast<std::ptrdiff_t>(i) * stride;
    }
};

// ─── Core sampler class ──────────────────────────────────────────────────────

class MinibatchSampler {
public:
    // ── Farthest-point sampling ──────────────────────────────────────────────
    //
    // Deterministic O(k*n) algorithm:
    //   1. Seed with the lowest-energy chromosome
    //   2. Maintain min_distance[i] = min RMSD from chrom[i] to any selected
    //   3. Iteratively pick argmax min_distance
    //
    // coord_cache:   pre-computed coordinate cache (from build_coord_cache)
    // energies:      app_evalue for each chromosome (lower = better)
    // n_atoms:       number of ligand atoms (for RMSD normalisation)
    // k:             target number of representatives
    // verbose:       print timing/summary

    static SampleResult farthest_point_sample(
        const CoordCache& coord_cache,
        const double*     energies,
        int               n_atoms,
        int               k,
        bool              verbose = true);

    // ── Quality-weighted sampling ─────────────────────────────────────────────
    //
    // Blends spatial diversity with energy quality:
    //   score[i] = alpha * normalised_min_distance[i]
    //           + (1 - alpha) * normalised_quality[i]
    //
    // where quality[i] = 1 - (E[i] - E_min) / (E_max - E_min)
    // Picks argmax score at each iteration.

    static SampleResult quality_weighted_sample(
        const CoordCache& coord_cache,
        const double*     energies,
        int               n_atoms,
        int               k,
        float             alpha = 0.7f,
        bool              verbose = true);

    // ── Assign remaining chromosomes to nearest selected ──────────────────────
    //
    // After FPS selects k representatives, assign every original chromosome
    // to its nearest selected representative.  Returns assignment vector
    // where assignment[i] = index into selected_indices.
    // Uses OpenMP for parallel distance computation.

    static std::vector<int> assign_to_nearest(
        const CoordCache&       coord_cache,
        const std::vector<int>& selected_indices,
        int                     n_atoms);

private:
    // Internal: batch RMSD between one reference and all others
    static float rmsd_to_ref(
        const CoordCache& cache,
        int               ref_idx,
        int               other_idx,
        int               n_atoms) noexcept {
        const float* a = cache[ref_idx];
        const float* b = cache[other_idx];
        float ssd = flexaids::sum_sq_distances_f(a, b, n_atoms * 3);
        return std::sqrt(ssd / static_cast<float>(n_atoms));
    }
};

} // namespace minibatch

#endif // MINIBATCHSAMPLER_H
