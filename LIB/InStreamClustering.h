// InStreamClustering.h — Online medoid-based clustering for FlexAIDdS GA
//
// Maintains a bounded set of cluster medoids during GA evolution.
// Every N generations, elite chromosomes are merged into the medoid set:
//   - If an elite is within rmsd_threshold of an existing medoid, it is
//     absorbed (member_count++, score updated if better).
//   - Otherwise, a new medoid is created.
//   - If the medoid set exceeds max_representatives, the two closest
//     medoids are merged.
//
// Thread-safe for OpenMP parallel GA via mutex.
//
// Apache-2.0 — see LICENSE.

#ifndef INSTREAM_CLUSTERING_H
#define INSTREAM_CLUSTERING_H

#include <vector>
#include <mutex>
#include <cstddef>
#include <cstdint>

namespace flexaids {

/// A single cluster medoid: its gene IC values, best score, and bookkeeping.
struct ClusterMedoid {
    std::vector<float> genes_ic;  ///< Gene IC values (num_genes floats)
    double best_score;            ///< Best app_evalue seen in this cluster
    int    first_seen_gen;        ///< Generation when this medoid was created
    int    last_updated_gen;      ///< Last generation that updated this medoid
    int    member_count;          ///< Number of elites absorbed into this cluster

    ClusterMedoid() : best_score(1e30), first_seen_gen(0),
                      last_updated_gen(0), member_count(0) {}

    explicit ClusterMedoid(int num_genes)
        : genes_ic(static_cast<size_t>(num_genes), 0.0f),
          best_score(1e30), first_seen_gen(0),
          last_updated_gen(0), member_count(0) {}
};

/// Online medoid clustering that runs during GA evolution.
///
/// Usage:
///   InStreamCluster isc(rmsd_threshold, max_medoids, num_genes);
///   // every K generations:
///   isc.merge_elites(genes_ic, scores, generation, K);
///   // after GA finishes:
///   auto medoids = isc.finalize();
class InStreamCluster {
public:
    /// Constructor.
    /// @param rmsd_threshold  RMSD below which two conformers are "same cluster" (Angstrom)
    /// @param max_representatives  Upper bound on medoid count; triggers merge when exceeded
    /// @param num_genes  Number of GA genes per chromosome
    InStreamCluster(float rmsd_threshold = 2.0f,
                    int max_representatives = 5000,
                    int num_genes = 0);

    /// Merge a batch of elite individuals into the medoid set.
    /// @param genes_ic  Flat array: [elite_idx * num_genes + gene_idx] gene IC values
    /// @param scores    Array: app_evalue for each elite (lower = better)
    /// @param n_elites  Number of elites in this batch
    /// @param generation  Current GA generation number
    /// @param num_genes  Number of genes per chromosome (must match constructor or each call)
    void merge_elites(const float* genes_ic,
                      const double* scores,
                      int n_elites,
                      int generation,
                      int num_genes);

    /// Convenience overload: merge elites from chromosome arrays.
    /// Extracts gene IC values from chromosome::genes[i].to_ic.
    /// Defined in the .cpp to avoid including gaboom.h here.
    void merge_elites_from_chrom(const void* chroms, int n_elites,
                                 int generation, int num_genes,
                                 int chrom_stride_bytes);

    /// Finalize: return the medoid set and perform any last merges.
    /// @return Vector of cluster medoids sorted by best_score (ascending = best first).
    std::vector<ClusterMedoid> finalize();

    /// Current number of clusters (medoids).
    int cluster_count() const;

    /// Total number of individuals merged across all calls.
    int64_t total_merged() const { return total_merged_; }

    /// Get the RMSD threshold.
    float rmsd_threshold() const { return rmsd_threshold_; }

    /// Reset to empty state (reuses allocated capacity).
    void reset();

private:
    float rmsd_threshold_;
    int   max_representatives_;
    int   num_genes_;
    int64_t total_merged_;

    std::vector<ClusterMedoid> medoids_;
    mutable std::mutex mutex_;

    /// Compute RMSD between two gene IC vectors.
    /// Returns sqrt(mean(squared_difference)), same formula as coordinate RMSD
    /// but operating in gene-space.  This is an approximation of the Cartesian
    /// RMSD that avoids the expensive buildcc() coordinate reconstruction.
    static float gene_rmsd(const float* a, const float* b, int n);

    /// Merge the two closest medoids (by gene RMSD) into one.
    void merge_closest_pair(int generation);
};

} // namespace flexaids

#endif // INSTREAM_CLUSTERING_H
