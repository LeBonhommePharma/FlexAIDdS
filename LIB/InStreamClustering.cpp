// InStreamClustering.cpp — Online medoid-based clustering for FlexAIDdS GA
//
// Implementation of InStreamCluster: bounded medoid maintenance during
// GA evolution.  Replaces O(n log n) post-hoc clustering on ~500K snapshots
// with O(K * M) per merge, where K = elites per batch, M = current medoids.
//
// Apache-2.0 — see LICENSE.

#include "InStreamClustering.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <cassert>

namespace flexaids {

// ─── Construction / Reset ─────────────────────────────────────────────────

InStreamCluster::InStreamCluster(float rmsd_threshold,
                                 int max_representatives,
                                 int num_genes)
    : rmsd_threshold_(rmsd_threshold)
    , max_representatives_(max_representatives > 0 ? max_representatives : 5000)
    , num_genes_(num_genes)
    , total_merged_(0)
{
    if (num_genes_ > 0) {
        medoids_.reserve(static_cast<size_t>(max_representatives_));
    }
}

void InStreamCluster::reset()
{
    medoids_.clear();
    total_merged_ = 0;
}

// ─── Gene-space RMSD ─────────────────────────────────────────────────────

float InStreamCluster::gene_rmsd(const float* a, const float* b, int n)
{
    if (n <= 0) return 0.0f;

    // First 3 genes are typically translation (x, y, z in Angstroms),
    // which dominate the RMSD signal.  The remaining genes are angles
    // (degrees) and indices — we weight translations more heavily.
    //
    // For a fast approximation, we compute a plain Euclidean distance
    // on the translation genes (first 3) and add a small contribution
    // from the angular genes.  This avoids the need for full coordinate
    // reconstruction via buildcc().
    //
    // Translation contribution: sqrt(mean(dx^2 + dy^2 + dz^2))
    // where n_trans = min(3, n).
    const int n_trans = (n >= 3) ? 3 : n;
    float trans_sq = 0.0f;
    for (int i = 0; i < n_trans; ++i) {
        float d = a[i] - b[i];
        trans_sq += d * d;
    }

    // Angular contribution (small weight: degrees ~ 0.01 Angstroms impact)
    // Only include if there are more than 3 genes.
    float angle_sq = 0.0f;
    for (int i = n_trans; i < n; ++i) {
        float d = a[i] - b[i];
        // Wrap angular differences to [-180, 180]
        // (genes may be dihedral angles in degrees)
        if (d > 180.0f)  d -= 360.0f;
        if (d < -180.0f) d += 360.0f;
        angle_sq += d * d;
    }

    // Weighted combination: translation RMSD dominates
    // Angular genes contribute at reduced weight (0.05 per gene)
    float rmsd_sq = trans_sq / static_cast<float>(n_trans);
    if (n > n_trans) {
        float angle_weight = 0.05f;  // each angular degree contributes ~0.05 Ang equiv
        rmsd_sq += angle_weight * angle_sq / static_cast<float>(n - n_trans);
    }

    return std::sqrt(std::max(rmsd_sq, 0.0f));
}

// ─── Core merge algorithm ────────────────────────────────────────────────

void InStreamCluster::merge_elites(const float* genes_ic,
                                   const double* scores,
                                   int n_elites,
                                   int generation,
                                   int num_genes)
{
    // Guard: empty input
    if (n_elites <= 0 || num_genes <= 0) return;
    if (!genes_ic || !scores) return;

    // Update num_genes_ on first call if not set at construction
    if (num_genes_ <= 0) {
        num_genes_ = num_genes;
        medoids_.reserve(static_cast<size_t>(max_representatives_));
    }

    std::lock_guard<std::mutex> lock(mutex_);

    for (int e = 0; e < n_elites; ++e) {
        const float* elite_genes = genes_ic + static_cast<ptrdiff_t>(e) * num_genes;
        double elite_score = scores[e];

        // Skip individuals with invalid scores
        if (!std::isfinite(elite_score)) continue;

        // Find closest existing medoid
        float best_dist = std::numeric_limits<float>::max();
        int   best_idx  = -1;

        const int n_med = static_cast<int>(medoids_.size());
        for (int m = 0; m < n_med; ++m) {
            if (static_cast<int>(medoids_[m].genes_ic.size()) != num_genes) continue;
            float d = gene_rmsd(elite_genes, medoids_[m].genes_ic.data(), num_genes);
            if (d < best_dist) {
                best_dist = d;
                best_idx = m;
            }
        }

        if (best_idx >= 0 && best_dist < rmsd_threshold_) {
            // Merge: absorb into existing cluster
            ClusterMedoid& med = medoids_[best_idx];
            med.member_count++;
            med.last_updated_gen = generation;
            if (elite_score < med.best_score) {
                med.best_score = elite_score;
                // Update medoid center to the better-scoring individual
                // (greedy: track the actual best member)
                std::copy(elite_genes, elite_genes + num_genes,
                          med.genes_ic.begin());
            }
        } else {
            // Create new medoid
            ClusterMedoid new_med;
            new_med.genes_ic.assign(elite_genes, elite_genes + num_genes);
            new_med.best_score = elite_score;
            new_med.first_seen_gen = generation;
            new_med.last_updated_gen = generation;
            new_med.member_count = 1;
            medoids_.push_back(std::move(new_med));
        }

        total_merged_++;

        // If over capacity, merge the closest pair
        while (static_cast<int>(medoids_.size()) > max_representatives_) {
            merge_closest_pair(generation);
        }
    }
}

void InStreamCluster::merge_closest_pair(int generation)
{
    const int n = static_cast<int>(medoids_.size());
    if (n < 2) return;

    // Find the pair of medoids with smallest gene RMSD
    float min_dist = std::numeric_limits<float>::max();
    int   mi = 0, mj = 1;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int ni = static_cast<int>(medoids_[i].genes_ic.size());
            int nj = static_cast<int>(medoids_[j].genes_ic.size());
            int ng = std::min(ni, nj);
            if (ng <= 0) continue;
            float d = gene_rmsd(medoids_[i].genes_ic.data(),
                                medoids_[j].genes_ic.data(), ng);
            if (d < min_dist) {
                min_dist = d;
                mi = i;
                mj = j;
            }
        }
    }

    // Merge mj into mi: keep the one with better (lower) score
    ClusterMedoid& a = medoids_[mi];
    ClusterMedoid& b = medoids_[mj];

    if (b.best_score < a.best_score) {
        // b is better: make b the medoid center, merge a into b
        a.genes_ic = std::move(b.genes_ic);
        a.best_score = b.best_score;
    }
    a.member_count += b.member_count;
    a.last_updated_gen = std::max(a.last_updated_gen, b.last_updated_gen);

    // Remove mj (swap with last, pop back)
    if (mj < n - 1) {
        medoids_[mj] = std::move(medoids_[n - 1]);
    }
    medoids_.pop_back();
}

// ─── Finalize ─────────────────────────────────────────────────────────────

std::vector<ClusterMedoid> InStreamCluster::finalize()
{
    std::lock_guard<std::mutex> lock(mutex_);

    // Sort medoids by best_score ascending (best clusters first)
    std::sort(medoids_.begin(), medoids_.end(),
              [](const ClusterMedoid& a, const ClusterMedoid& b) {
                  return a.best_score < b.best_score;
              });

    return medoids_;  // return a copy
}

// ─── Accessors ────────────────────────────────────────────────────────────

int InStreamCluster::cluster_count() const
{
    std::lock_guard<std::mutex> lock(mutex_);
    return static_cast<int>(medoids_.size());
}

} // namespace flexaids
