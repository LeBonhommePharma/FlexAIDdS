// MinibatchSampler.cpp — Farthest-point sampling pre-filter for clustering
//
// Implements O(k*n) farthest-point sampling and quality-weighted sampling
// with SIMD-accelerated batch distance computation and OpenMP parallelism.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#include "MinibatchSampler.h"
#include "gaboom.h"        // chromosome

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace minibatch {

// ═══════════════════════════════════════════════════════════════════════════════
//  Farthest-point sampling
// ═══════════════════════════════════════════════════════════════════════════════

SampleResult MinibatchSampler::farthest_point_sample(
    const CoordCache& coord_cache,
    const double*     energies,
    int               n_atoms,
    int               k,
    bool              verbose)
{
    const int n = coord_cache.n_chrom;

    SampleResult result;
    result.n_original  = n;
    result.elapsed_ms  = 0.0;

    // Edge cases
    if (n <= 0 || k <= 0 || energies == nullptr) {
        result.n_selected = 0;
        return result;
    }
    if (k >= n) {
        // Nothing to subsample; return identity mapping
        result.selected_indices.resize(static_cast<std::size_t>(n));
        result.assignment.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            result.selected_indices[i] = i;
            result.assignment[i] = i;
        }
        result.n_selected = n;
        return result;
    }

    auto t0 = std::chrono::steady_clock::now();

    // min_dist[i] = min RMSD from chrom[i] to any already-selected point
    std::vector<float> min_dist(static_cast<std::size_t>(n),
                                 std::numeric_limits<float>::max());

    // selected: indices of selected representatives
    std::vector<int> selected;
    selected.reserve(static_cast<std::size_t>(k));

    // Step 1: seed with lowest-energy chromosome
    int seed_idx = 0;
    double best_energy = energies[0];
    for (int i = 1; i < n; ++i) {
        if (energies[i] < best_energy) {
            best_energy = energies[i];
            seed_idx = i;
        }
    }
    selected.push_back(seed_idx);
    min_dist[static_cast<std::size_t>(seed_idx)] = 0.0f;

    // Update min_dist after seeding
    const int stride = coord_cache.stride;
    const float* seed_coor = coord_cache[seed_idx];
    const float inv_n_atoms = 1.0f / static_cast<float>(n_atoms);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < n; ++i) {
        if (i == seed_idx) continue;
        const float* coor_i = coord_cache[i];
        float ssd = flexaids::sum_sq_distances_f(coor_i, seed_coor, stride);
        min_dist[i] = std::sqrt(ssd * inv_n_atoms);
    }

    // Step 2: iteratively pick the farthest point
    for (int step = 1; step < k; ++step) {
        // Find argmax min_dist among non-selected
        int farthest = -1;
        float max_dist = -1.0f;
        for (int i = 0; i < n; ++i) {
            if (min_dist[i] > max_dist) {
                max_dist = min_dist[i];
                farthest = i;
            }
        }

        if (farthest < 0 || max_dist <= 0.0f) {
            // All remaining points are duplicates or degenerate
            break;
        }

        selected.push_back(farthest);
        min_dist[static_cast<std::size_t>(farthest)] = 0.0f;

        // Update min_dist: parallel distance computation from newly selected
        const float* far_coor = coord_cache[farthest];

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int i = 0; i < n; ++i) {
            if (min_dist[i] <= 0.0f) continue;  // already selected
            const float* coor_i = coord_cache[i];
            float ssd = flexaids::sum_sq_distances_f(coor_i, far_coor, stride);
            float rmsd = std::sqrt(ssd * inv_n_atoms);
            if (rmsd < min_dist[i]) {
                min_dist[i] = rmsd;
            }
        }
    }

    result.selected_indices = std::move(selected);
    result.n_selected = static_cast<int>(result.selected_indices.size());

    // Build assignment: each original chrom -> nearest selected
    result.assignment = assign_to_nearest(coord_cache, result.selected_indices, n_atoms);

    auto t1 = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (verbose) {
        std::printf(
            "╔══ MinibatchSampler (FPS) ═══════════════════════════════════╗\n"
            "║  Input poses     : %8d                                  ║\n"
            "║  Selected (k)    : %8d                                  ║\n"
            "║  Reduction ratio : %8.1f : 1                             ║\n"
            "║  Elapsed         : %8.1f ms                              ║\n"
            "╚═══════════════════════════════════════════════════════════════╝\n",
            result.n_original,
            result.n_selected,
            static_cast<double>(result.n_original) / std::max(1, result.n_selected),
            result.elapsed_ms);
    }

    return result;
}


// ═══════════════════════════════════════════════════════════════════════════════
//  Quality-weighted sampling
// ═══════════════════════════════════════════════════════════════════════════════

SampleResult MinibatchSampler::quality_weighted_sample(
    const CoordCache& coord_cache,
    const double*     energies,
    int               n_atoms,
    int               k,
    float             alpha,
    bool              verbose)
{
    const int n = coord_cache.n_chrom;

    SampleResult result;
    result.n_original  = n;
    result.elapsed_ms  = 0.0;

    // Edge cases
    if (n <= 0 || k <= 0 || energies == nullptr) {
        result.n_selected = 0;
        return result;
    }
    if (k >= n) {
        result.selected_indices.resize(static_cast<std::size_t>(n));
        result.assignment.resize(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            result.selected_indices[i] = i;
            result.assignment[i] = i;
        }
        result.n_selected = n;
        return result;
    }

    auto t0 = std::chrono::steady_clock::now();

    // Pre-compute energy range for quality normalisation
    double e_min = energies[0], e_max = energies[0];
    for (int i = 1; i < n; ++i) {
        if (energies[i] < e_min) e_min = energies[i];
        if (energies[i] > e_max) e_max = energies[i];
    }
    const double e_range = (e_max - e_min) > 1e-12 ? (e_max - e_min) : 1.0;

    // min_dist[i] = min RMSD to any selected
    std::vector<float> min_dist(static_cast<std::size_t>(n),
                                 std::numeric_limits<float>::max());

    // quality[i] = 1 - (E[i] - E_min) / E_range  (higher = better)
    std::vector<float> quality(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        quality[i] = 1.0f - static_cast<float>((energies[i] - e_min) / e_range);
    }

    std::vector<int> selected;
    selected.reserve(static_cast<std::size_t>(k));

    // Seed with lowest-energy chromosome
    int seed_idx = 0;
    double best_energy = energies[0];
    for (int i = 1; i < n; ++i) {
        if (energies[i] < best_energy) {
            best_energy = energies[i];
            seed_idx = i;
        }
    }
    selected.push_back(seed_idx);
    min_dist[static_cast<std::size_t>(seed_idx)] = 0.0f;

    // Track max min_dist for normalisation
    float max_min_dist = 0.0f;
    const int stride = coord_cache.stride;
    const float inv_n_atoms = 1.0f / static_cast<float>(n_atoms);

    // Update min_dist after seeding
    const float* seed_coor = coord_cache[seed_idx];
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(max:max_min_dist)
    #endif
    for (int i = 0; i < n; ++i) {
        if (i == seed_idx) continue;
        const float* coor_i = coord_cache[i];
        float ssd = flexaids::sum_sq_distances_f(coor_i, seed_coor, stride);
        float rmsd = std::sqrt(ssd * inv_n_atoms);
        min_dist[i] = rmsd;
        if (rmsd > max_min_dist) max_min_dist = rmsd;
    }

    // Iteratively pick the best blend of diversity and quality
    for (int step = 1; step < k; ++step) {
        const float norm_dist = max_min_dist > 1e-12f ? max_min_dist : 1.0f;

        int best = -1;
        float best_score = -1.0f;
        for (int i = 0; i < n; ++i) {
            if (min_dist[i] <= 0.0f) continue;  // already selected
            float score = alpha * (min_dist[i] / norm_dist) + (1.0f - alpha) * quality[i];
            if (score > best_score) {
                best_score = score;
                best = i;
            }
        }

        if (best < 0) break;

        selected.push_back(best);
        min_dist[static_cast<std::size_t>(best)] = 0.0f;

        // Update min_dist from newly selected point
        const float* best_coor = coord_cache[best];
        float local_max = 0.0f;

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(max:local_max)
        #endif
        for (int i = 0; i < n; ++i) {
            if (min_dist[i] <= 0.0f) continue;
            const float* coor_i = coord_cache[i];
            float ssd = flexaids::sum_sq_distances_f(coor_i, best_coor, stride);
            float rmsd = std::sqrt(ssd * inv_n_atoms);
            if (rmsd < min_dist[i]) {
                min_dist[i] = rmsd;
            }
            if (min_dist[i] > local_max) local_max = min_dist[i];
        }
        max_min_dist = local_max;
    }

    result.selected_indices = std::move(selected);
    result.n_selected = static_cast<int>(result.selected_indices.size());
    result.assignment = assign_to_nearest(coord_cache, result.selected_indices, n_atoms);

    auto t1 = std::chrono::steady_clock::now();
    result.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (verbose) {
        std::printf(
            "╔══ MinibatchSampler (quality-weighted, alpha=%.2f) ═════════╗\n"
            "║  Input poses     : %8d                                  ║\n"
            "║  Selected (k)    : %8d                                  ║\n"
            "║  Reduction ratio : %8.1f : 1                             ║\n"
            "║  Elapsed         : %8.1f ms                              ║\n"
            "╚═══════════════════════════════════════════════════════════════╝\n",
            static_cast<double>(alpha),
            result.n_original,
            result.n_selected,
            static_cast<double>(result.n_original) / std::max(1, result.n_selected),
            result.elapsed_ms);
    }

    return result;
}


// ═══════════════════════════════════════════════════════════════════════════════
//  Assign remaining chromosomes to nearest selected
// ═══════════════════════════════════════════════════════════════════════════════

std::vector<int> MinibatchSampler::assign_to_nearest(
    const CoordCache&       coord_cache,
    const std::vector<int>& selected_indices,
    int                     n_atoms)
{
    const int n = coord_cache.n_chrom;
    const int k = static_cast<int>(selected_indices.size());
    const int stride = coord_cache.stride;
    const float inv_n_atoms = 1.0f / static_cast<float>(n_atoms);

    std::vector<int> assignment(static_cast<std::size_t>(n), 0);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (int i = 0; i < n; ++i) {
        const float* coor_i = coord_cache[i];
        float best_rmsd = std::numeric_limits<float>::max();
        int   best_idx  = 0;

        for (int s = 0; s < k; ++s) {
            const float* coor_s = coord_cache[selected_indices[s]];
            float ssd = flexaids::sum_sq_distances_f(coor_i, coor_s, stride);
            float rmsd = std::sqrt(ssd * inv_n_atoms);
            if (rmsd < best_rmsd) {
                best_rmsd = rmsd;
                best_idx  = s;
            }
        }
        assignment[i] = best_idx;
    }

    return assignment;
}

} // namespace minibatch
