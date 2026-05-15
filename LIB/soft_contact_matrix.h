// soft_contact_matrix.h — 256×256 soft contact interaction matrix
//
// Header-only implementation providing:
//   - 256 KB cache-aligned matrix with O(1) lookup
//   - AVX2 8-wide gather for batch scoring
//   - FastOPTICS super-cluster detection on row vectors
//   - Gaussian supercluster bias modulation
//   - 256→40 SYBYL projection via base_to_sybyl_parent()
//   - Binary I/O with SHNN magic header
//
// The matrix stores pairwise soft-contact interaction energies between
// 256 atom types encoded by atom_typing_256.h.
#pragma once

#include "atom_typing_256.h"
#include "nrgrank_matrix.h"
#include <vector>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>

#ifdef __AVX512F__
#include <immintrin.h>
#elif defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSE4_2__)
#include <nmmintrin.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Dense>

namespace scm {

// ─── constants ──────────────────────────────────────────────────────────────
inline constexpr int    MATRIX_DIM   = 256;
inline constexpr int    MATRIX_SIZE  = MATRIX_DIM * MATRIX_DIM;  // 65536
inline constexpr size_t MATRIX_BYTES = MATRIX_SIZE * sizeof(float); // 256 KB
inline constexpr int    LEGACY_SYBYL_DIM = 40;

inline constexpr const char* ATOM256_SCHEMA_ID =
    "flexaidds.atom256.v1.base64.charge2.hbond1";
inline constexpr const char* LEGACY_40_FALLBACK_SCHEMA_ID =
    "flexaidds.nrgrank.v1.sybyl40";

// SHNN binary blob magic
inline constexpr uint32_t SHNN_MAGIC   = 0x4E4E4853; // "SHNN" little-endian
inline constexpr uint32_t SHNN_VERSION = 1;

// ─── aligned storage ────────────────────────────────────────────────────────

struct alignas(64) SoftContactMatrix {
    float data[MATRIX_SIZE];

    // O(1) lookup
    float lookup(int type_i, int type_j) const noexcept {
        return data[type_i * MATRIX_DIM + type_j];
    }

    void set(int type_i, int type_j, float value) noexcept {
        data[type_i * MATRIX_DIM + type_j] = value;
    }

    // Row pointer for batch operations
    const float* row(int type_i) const noexcept {
        return &data[type_i * MATRIX_DIM];
    }

    float* row(int type_i) noexcept {
        return &data[type_i * MATRIX_DIM];
    }

    void zero() noexcept { std::memset(data, 0, MATRIX_BYTES); }

    void symmetrise() noexcept {
        for (int i = 0; i < MATRIX_DIM; ++i)
            for (int j = i + 1; j < MATRIX_DIM; ++j) {
                float avg = (data[i * MATRIX_DIM + j] +
                             data[j * MATRIX_DIM + i]) * 0.5f;
                data[i * MATRIX_DIM + j] = avg;
                data[j * MATRIX_DIM + i] = avg;
            }
    }

    static SoftContactMatrix from_nrgrank40_fallback() noexcept {
        SoftContactMatrix mat;
        mat.zero();
        for (int ci = 0; ci < MATRIX_DIM; ++ci) {
            const int si =
                atom256::base_to_sybyl_parent(atom256::get_base(ci));
            for (int cj = 0; cj < MATRIX_DIM; ++cj) {
                const int sj =
                    atom256::base_to_sybyl_parent(atom256::get_base(cj));
                if (si >= 1 && si <= LEGACY_SYBYL_DIM &&
                    sj >= 1 && sj <= LEGACY_SYBYL_DIM) {
                    mat.set(ci, cj, static_cast<float>(
                        nrgrank::kEnergyMatrix[si][sj]));
                }
            }
        }
        return mat;
    }

    bool load_or_nrgrank40_fallback(const char* path,
                                    bool* used_fallback = nullptr,
                                    bool fallback_on_load_error = false) noexcept {
        if (path != nullptr && path[0] != '\0' && load(path)) {
            if (used_fallback) *used_fallback = false;
            return true;
        }
        if (path != nullptr && path[0] != '\0' && !fallback_on_load_error) {
            if (used_fallback) *used_fallback = false;
            return false;
        }
        *this = from_nrgrank40_fallback();
        if (used_fallback) *used_fallback = true;
        return true;
    }

    // ── AVX2 batch scoring ──────────────────────────────────────────────

#ifdef __AVX2__
    // Score N contacts: sum of matrix[type_a[k]][type_b[k]] * area[k]
    // Uses vgatherdps for 8-wide gather.
    float score_contacts_avx2(const uint8_t* type_a, const uint8_t* type_b,
                               const float* areas, int n) const noexcept {
        __m256 acc = _mm256_setzero_ps();
        int k = 0;
        for (; k + 7 < n; k += 8) {
            // Build indices: type_a[k]*256 + type_b[k]
            alignas(32) int32_t idx[8];
            for (int q = 0; q < 8; ++q)
                idx[q] = static_cast<int32_t>(type_a[k + q]) * MATRIX_DIM +
                         static_cast<int32_t>(type_b[k + q]);
            __m256i vidx = _mm256_load_si256((__m256i*)idx);
            __m256  vals = _mm256_i32gather_ps(data, vidx, sizeof(float));
            __m256  area = _mm256_loadu_ps(areas + k);
            acc = _mm256_fmadd_ps(vals, area, acc);
        }
        // Horizontal sum
        alignas(32) float tmp[8];
        _mm256_store_ps(tmp, acc);
        float sum = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7];
        // Scalar tail
        for (; k < n; ++k)
            sum += data[type_a[k] * MATRIX_DIM + type_b[k]] * areas[k];
        return sum;
    }
#endif

#ifdef __AVX512F__
    // Score N contacts with AVX-512: 16-wide gather
    float score_contacts_avx512(const uint8_t* type_a, const uint8_t* type_b,
                                 const float* areas, int n) const noexcept {
        __m512 acc = _mm512_setzero_ps();
        int k = 0;
        for (; k + 15 < n; k += 16) {
            alignas(64) int32_t idx[16];
            for (int q = 0; q < 16; ++q)
                idx[q] = static_cast<int32_t>(type_a[k + q]) * MATRIX_DIM +
                         static_cast<int32_t>(type_b[k + q]);
            __m512i vidx = _mm512_load_epi32(idx);
            __m512  vals = _mm512_i32gather_ps(vidx, data, sizeof(float));
            __m512  area = _mm512_loadu_ps(areas + k);
            acc = _mm512_fmadd_ps(vals, area, acc);
        }
        float sum = _mm512_reduce_add_ps(acc);
        for (; k < n; ++k)
            sum += data[type_a[k] * MATRIX_DIM + type_b[k]] * areas[k];
        return sum;
    }
#endif

#ifdef __SSE4_2__
    // Score N contacts with SSE4.2: 4-wide manual gather (no vgatherdps).
    // SSE lacks gather, so we load each matrix element individually and
    // pack into a __m128 for the FMA-equivalent multiply-add.
    float score_contacts_sse42(const uint8_t* type_a, const uint8_t* type_b,
                                const float* areas, int n) const noexcept {
        __m128 acc = _mm_setzero_ps();
        int k = 0;
        for (; k + 3 < n; k += 4) {
            // Manual 4-wide gather: load 4 matrix values individually
            alignas(16) float vals[4];
            vals[0] = data[static_cast<int>(type_a[k])     * MATRIX_DIM + type_b[k]];
            vals[1] = data[static_cast<int>(type_a[k + 1]) * MATRIX_DIM + type_b[k + 1]];
            vals[2] = data[static_cast<int>(type_a[k + 2]) * MATRIX_DIM + type_b[k + 2]];
            vals[3] = data[static_cast<int>(type_a[k + 3]) * MATRIX_DIM + type_b[k + 3]];
            __m128 v    = _mm_load_ps(vals);
            __m128 area = _mm_loadu_ps(areas + k);
            // acc += v * area  (no FMA on SSE4.2, use mul + add)
            acc = _mm_add_ps(acc, _mm_mul_ps(v, area));
        }
        // Horizontal sum
        __m128 hi = _mm_movehl_ps(acc, acc);
        __m128 s  = _mm_add_ps(acc, hi);
        __m128 s2 = _mm_movehdup_ps(s);
        float sum = _mm_cvtss_f32(_mm_add_ss(s, s2));
        // Scalar tail
        for (; k < n; ++k)
            sum += data[type_a[k] * MATRIX_DIM + type_b[k]] * areas[k];
        return sum;
    }
#endif

    // Dispatch: AVX-512 → AVX2 → SSE4.2 → scalar
    float score_contacts(const uint8_t* type_a, const uint8_t* type_b,
                          const float* areas, int n) const noexcept {
#ifdef __AVX512F__
        return score_contacts_avx512(type_a, type_b, areas, n);
#elif defined(__AVX2__)
        return score_contacts_avx2(type_a, type_b, areas, n);
#elif defined(__SSE4_2__)
        return score_contacts_sse42(type_a, type_b, areas, n);
#else
        float sum = 0.0f;
        for (int k = 0; k < n; ++k)
            sum += data[type_a[k] * MATRIX_DIM + type_b[k]] * areas[k];
        return sum;
#endif
    }

    // ── binary I/O ──────────────────────────────────────────────────────

    bool save(const char* path) const {
        FILE* fp = std::fopen(path, "wb");
        if (!fp) return false;
        uint32_t magic = SHNN_MAGIC;
        uint32_t ver   = SHNN_VERSION;
        uint32_t dim   = MATRIX_DIM;
        std::fwrite(&magic, 4, 1, fp);
        std::fwrite(&ver,   4, 1, fp);
        std::fwrite(&dim,   4, 1, fp);
        size_t written = std::fwrite(data, sizeof(float), MATRIX_SIZE, fp);
        std::fclose(fp);
        return written == MATRIX_SIZE;
    }

    bool load(const char* path) {
        FILE* fp = std::fopen(path, "rb");
        if (!fp) return false;
        uint32_t magic, ver, dim;
        if (std::fread(&magic, 4, 1, fp) != 1 ||
            std::fread(&ver,   4, 1, fp) != 1 ||
            std::fread(&dim,   4, 1, fp) != 1) {
            std::fclose(fp); return false;
        }
        if (magic != SHNN_MAGIC || dim != MATRIX_DIM) {
            std::fclose(fp); return false;
        }
        size_t nread = std::fread(data, sizeof(float), MATRIX_SIZE, fp);
        std::fclose(fp);
        return nread == MATRIX_SIZE;
    }

    // ── 256→40 projection ───────────────────────────────────────────────

    // Projects this 256×256 matrix to a 40×40 SYBYL matrix by averaging
    // all 256-type codes that map to each SYBYL parent.
    // Uses OpenMP for outer-loop parallelism when available.
    std::array<float, 40 * 40> project_to_40x40() const {
        // Precompute 256→40 mapping to avoid repeated lookups
        std::array<int, MATRIX_DIM> sybyl_map;
        for (int c = 0; c < MATRIX_DIM; ++c)
            sybyl_map[c] = atom256::base_to_sybyl_parent(atom256::get_base(c)) - 1;

        std::array<float, 40 * 40> out{};
        std::array<int, 40 * 40> counts{};
        out.fill(0.0f);
        counts.fill(0);

#ifdef _OPENMP
        // Thread-private accumulators
        int n_threads = omp_get_max_threads();
        std::vector<std::array<float, 40*40>> t_out(n_threads);
        std::vector<std::array<int, 40*40>> t_cnt(n_threads);
        for (auto& a : t_out) a.fill(0.0f);
        for (auto& a : t_cnt) a.fill(0);

        #pragma omp parallel for schedule(static)
        for (int ci = 0; ci < MATRIX_DIM; ++ci) {
            int tid = omp_get_thread_num();
            int si = sybyl_map[ci];
            if (si < 0 || si >= 40) continue;
            for (int cj = 0; cj < MATRIX_DIM; ++cj) {
                int sj = sybyl_map[cj];
                if (sj < 0 || sj >= 40) continue;
                t_out[tid][si * 40 + sj] += data[ci * MATRIX_DIM + cj];
                t_cnt[tid][si * 40 + sj]++;
            }
        }
        for (int t = 0; t < n_threads; ++t)
            for (int k = 0; k < 40 * 40; ++k) {
                out[k] += t_out[t][k];
                counts[k] += t_cnt[t][k];
            }
#else
        for (int ci = 0; ci < MATRIX_DIM; ++ci) {
            int si = sybyl_map[ci];
            if (si < 0 || si >= 40) continue;
            for (int cj = 0; cj < MATRIX_DIM; ++cj) {
                int sj = sybyl_map[cj];
                if (sj < 0 || sj >= 40) continue;
                out[si * 40 + sj] += data[ci * MATRIX_DIM + cj];
                counts[si * 40 + sj]++;
            }
        }
#endif
        for (int k = 0; k < 40 * 40; ++k)
            if (counts[k] > 0) out[k] /= counts[k];
        return out;
    }
};

// ═══════════════════════════════════════════════════════════════════════════════
// FastOPTICS super-cluster detection on matrix row vectors
// ═══════════════════════════════════════════════════════════════════════════════
//
// Lightweight reimplementation of FastOPTICS (Fast Ordering Points To Identify
// the Clustering Structure) specialised for finding "super-clusters" among
// 256 row vectors of the contact matrix.  When many rows have nearly identical
// interaction profiles, they form a super-cluster — indicating that the atom
// types are functionally interchangeable at the scoring level (entropy collapse).
//
// This uses random projections + split-based neighbor estimation, matching the
// O(N log N) expected complexity of the original Kriegel et al. algorithm, but
// without the FlexAID GA coupling.

struct FOPTICSResult {
    std::vector<int>   order;           // OPTICS ordering of row indices
    std::vector<float> reachability;    // reachability distances
    std::vector<int>   cluster_labels; // cluster assignment (-1 = noise)
    int                n_clusters;
};

// A reduced matrix over atom-contact superclusters. Each cluster groups 256-type
// rows with similar interaction profiles; lookup_type() maps a concrete 256
// type pair to its cluster-pair mean. Noise rows are isolated by default so
// projection never creates artificial equivalence classes.
struct ContactSuperclusterMatrix {
    int n_clusters = 0;
    std::vector<int> type_to_cluster;
    std::vector<int> cluster_sizes;
    std::vector<float> data;

    float lookup_cluster(int ci, int cj) const noexcept {
        return data[static_cast<size_t>(ci) * n_clusters + cj];
    }

    float lookup_type(uint8_t type_i, uint8_t type_j) const noexcept {
        return lookup_cluster(type_to_cluster[type_i], type_to_cluster[type_j]);
    }

    double compression_ratio() const noexcept {
        if (n_clusters <= 0) return 0.0;
        return static_cast<double>(MATRIX_DIM) /
               static_cast<double>(n_clusters);
    }
};

inline ContactSuperclusterMatrix project_to_contact_superclusters(
        const SoftContactMatrix& mat,
        const FOPTICSResult& clusters,
        bool isolate_noise = true) {
    ContactSuperclusterMatrix out;
    out.type_to_cluster.assign(MATRIX_DIM, -1);

    int max_label = -1;
    for (int label : clusters.cluster_labels)
        if (label > max_label) max_label = label;
    std::vector<int> label_to_cluster(static_cast<size_t>(max_label + 1), -1);

    auto cluster_for_label = [&](int label) -> int {
        if (label < 0) {
            if (isolate_noise) return out.n_clusters++;
            return 0;
        }
        int& mapped = label_to_cluster[static_cast<size_t>(label)];
        if (mapped < 0) mapped = out.n_clusters++;
        return mapped;
    };

    for (int t = 0; t < MATRIX_DIM; ++t) {
        const int label = (t < static_cast<int>(clusters.cluster_labels.size()))
            ? clusters.cluster_labels[static_cast<size_t>(t)]
            : -1;
        out.type_to_cluster[static_cast<size_t>(t)] = cluster_for_label(label);
    }

    if (out.n_clusters == 0) {
        out.n_clusters = 1;
        std::fill(out.type_to_cluster.begin(), out.type_to_cluster.end(), 0);
    }

    out.cluster_sizes.assign(static_cast<size_t>(out.n_clusters), 0);
    for (int cid : out.type_to_cluster)
        out.cluster_sizes[static_cast<size_t>(cid)]++;

    const size_t reduced_size =
        static_cast<size_t>(out.n_clusters) * out.n_clusters;
    out.data.assign(reduced_size, 0.0f);
    std::vector<int> counts(reduced_size, 0);

    for (int ti = 0; ti < MATRIX_DIM; ++ti) {
        const int ci = out.type_to_cluster[static_cast<size_t>(ti)];
        for (int tj = 0; tj < MATRIX_DIM; ++tj) {
            const int cj = out.type_to_cluster[static_cast<size_t>(tj)];
            const size_t idx = static_cast<size_t>(ci) * out.n_clusters + cj;
            out.data[idx] += mat.lookup(ti, tj);
            counts[idx]++;
        }
    }

    for (size_t idx = 0; idx < reduced_size; ++idx)
        if (counts[idx] > 0)
            out.data[idx] /= static_cast<float>(counts[idx]);

    return out;
}

namespace detail {

// L2 distance between two 256-dimensional row vectors
// Dispatch: Eigen → AVX2 → scalar
inline float row_distance(const SoftContactMatrix& mat, int a, int b) {
    const float* ra = mat.row(a);
    const float* rb = mat.row(b);

    Eigen::Map<const Eigen::Array<float, MATRIX_DIM, 1>> va(ra);
    Eigen::Map<const Eigen::Array<float, MATRIX_DIM, 1>> vb(rb);
    return (va - vb).matrix().norm();
}

// Random projection: project each row onto a random unit vector
inline std::vector<float> random_projection(const SoftContactMatrix& mat,
                                             std::mt19937& rng) {
    std::normal_distribution<float> gauss(0.0f, 1.0f);
    std::vector<float> proj_vec(MATRIX_DIM);
    float norm = 0.0f;
    for (int d = 0; d < MATRIX_DIM; ++d) {
        proj_vec[d] = gauss(rng);
        norm += proj_vec[d] * proj_vec[d];
    }
    norm = 1.0f / std::sqrt(norm + 1e-12f);
    for (int d = 0; d < MATRIX_DIM; ++d) proj_vec[d] *= norm;

    std::vector<float> projected(MATRIX_DIM);
    for (int i = 0; i < MATRIX_DIM; ++i) {
        float dot = 0.0f;
        const float* r = mat.row(i);
        for (int d = 0; d < MATRIX_DIM; ++d) dot += r[d] * proj_vec[d];
        projected[i] = dot;
    }
    return projected;
}

} // namespace detail

// Main FastOPTICS entry point for row-vector clustering
inline FOPTICSResult find_super_clusters(const SoftContactMatrix& mat,
                                          int min_pts = 5,
                                          int n_projections = 20,
                                          uint32_t seed = 42) {
    constexpr int N = MATRIX_DIM;
    constexpr float UNDEFINED = std::numeric_limits<float>::infinity();

    std::mt19937 rng(seed);

    // Phase 1: Random projections to estimate densities and neighbors
    std::vector<float> inverse_density(N, 0.0f);
    std::vector<std::vector<int>> all_neighbors(N);

    for (int proj = 0; proj < n_projections; ++proj) {
        auto projected = detail::random_projection(mat, rng);

        // Sort indices by projected value
        std::vector<int> sorted_idx(N);
        std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
        std::sort(sorted_idx.begin(), sorted_idx.end(),
                  [&](int a, int b) { return projected[a] < projected[b]; });

        // Recursive split to find neighbor sets
        // Split sets: contiguous runs in sorted order of size ~2*min_pts
        int split_size = std::max(2 * min_pts, 4);
        for (int start = 0; start < N; start += split_size / 2) {
            int end = std::min(start + split_size, N);
            int count = end - start;

            // All pairs within this set are candidate neighbors
            for (int a = start; a < end; ++a) {
                int ia = sorted_idx[a];
                for (int b = a + 1; b < end; ++b) {
                    int ib = sorted_idx[b];
                    all_neighbors[ia].push_back(ib);
                    all_neighbors[ib].push_back(ia);
                }
                // Inverse density: distance to min_pts-th neighbor in projection
                if (count > min_pts) {
                    int far_idx = std::min(a + min_pts, end - 1);
                    float dist_1d = std::fabs(projected[sorted_idx[far_idx]] -
                                              projected[ia]);
                    inverse_density[ia] = std::max(inverse_density[ia], dist_1d);
                }
            }
        }
    }

    // Deduplicate neighbor lists
    for (int i = 0; i < N; ++i) {
        auto& nb = all_neighbors[i];
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    }

    // Phase 2: OPTICS ordering with priority queue
    std::vector<float> reach_dist(N, UNDEFINED);
    std::vector<bool>  processed(N, false);
    std::vector<int>   order;
    order.reserve(N);

    // Start from the point with lowest inverse density
    int start_idx = 0;
    float min_id = inverse_density[0];
    for (int i = 1; i < N; ++i) {
        if (inverse_density[i] < min_id) {
            min_id = inverse_density[i];
            start_idx = i;
        }
    }

    // Simple priority queue (N=256, so O(N²) is fine)
    auto process_point = [&](int p) {
        processed[p] = true;
        order.push_back(p);

        // Core distance: distance to min_pts-th nearest neighbor
        std::vector<float> nb_dists;
        nb_dists.reserve(all_neighbors[p].size());
        for (int nb : all_neighbors[p]) {
            if (!processed[nb])
                nb_dists.push_back(detail::row_distance(mat, p, nb));
        }
        std::sort(nb_dists.begin(), nb_dists.end());
        float core_dist = (static_cast<int>(nb_dists.size()) >= min_pts)
                          ? nb_dists[min_pts - 1]
                          : UNDEFINED;

        if (core_dist >= UNDEFINED) return;

        // Update reachability of unprocessed neighbors
        for (int nb : all_neighbors[p]) {
            if (processed[nb]) continue;
            float dist = detail::row_distance(mat, p, nb);
            float new_reach = std::max(core_dist, dist);
            if (new_reach < reach_dist[nb])
                reach_dist[nb] = new_reach;
        }
    };

    process_point(start_idx);

    for (int step = 1; step < N; ++step) {
        // Find unprocessed point with smallest reachability
        int next = -1;
        float best_reach = UNDEFINED;
        for (int i = 0; i < N; ++i) {
            if (!processed[i] && reach_dist[i] < best_reach) {
                best_reach = reach_dist[i];
                next = i;
            }
        }
        if (next < 0) {
            // No reachable points; pick any unprocessed
            for (int i = 0; i < N; ++i) {
                if (!processed[i]) { next = i; break; }
            }
        }
        if (next < 0) break;
        process_point(next);
    }

    // Phase 3: Extract clusters from reachability plot
    // Xi-steep-downward method (simplified): clusters are contiguous regions
    // below a threshold derived from the median reachability
    std::vector<float> valid_reach;
    for (int i : order) {
        if (reach_dist[i] < UNDEFINED)
            valid_reach.push_back(reach_dist[i]);
    }

    float threshold = UNDEFINED;
    if (!valid_reach.empty()) {
        std::sort(valid_reach.begin(), valid_reach.end());
        float median = valid_reach[valid_reach.size() / 2];
        threshold = median * 1.5f;  // 1.5× median as cluster boundary
    }

    std::vector<int> labels(N, -1);
    int current_label = 0;
    for (int idx = 0; idx < static_cast<int>(order.size()); ++idx) {
        int pt = order[idx];
        if (reach_dist[pt] <= threshold) {
            labels[pt] = current_label;
        } else {
            if (idx > 0 && labels[order[idx - 1]] >= 0) {
                current_label++;  // start new cluster after gap
            }
            // This point is noise or starts a new cluster
            // Check if next points are dense
            if (idx + 1 < static_cast<int>(order.size()) &&
                reach_dist[order[idx + 1]] <= threshold) {
                labels[pt] = current_label;
            }
        }
    }

    return { std::move(order), std::move(reach_dist),
             std::move(labels), current_label + 1 };
}

// ─── Gaussian supercluster bias modulation ──────────────────────────────────
// Applies a Gaussian mask that smooths interaction energies within each
// super-cluster: E'[i][j] = E[i][j] * (1 - alpha * G(cluster_i, cluster_j))
// where G is 1 when i,j are in the same cluster, decaying with inter-cluster
// distance.

inline void apply_supercluster_bias(SoftContactMatrix& mat,
                                     const FOPTICSResult& clusters,
                                     float alpha = 0.3f,
                                     float sigma = 2.0f) {
    // Compute cluster centroids (mean row vector per cluster)
    std::vector<std::vector<float>> centroids(clusters.n_clusters,
                                               std::vector<float>(MATRIX_DIM, 0.0f));
    std::vector<int> counts(clusters.n_clusters, 0);

    for (int i = 0; i < MATRIX_DIM; ++i) {
        int c = clusters.cluster_labels[i];
        if (c < 0) continue;
        const float* r = mat.row(i);
        for (int d = 0; d < MATRIX_DIM; ++d)
            centroids[c][d] += r[d];
        counts[c]++;
    }
    for (int c = 0; c < clusters.n_clusters; ++c) {
        if (counts[c] > 0)
            for (int d = 0; d < MATRIX_DIM; ++d)
                centroids[c][d] /= counts[c];
    }

    // Compute inter-cluster distances
    auto cluster_dist = [&](int c1, int c2) -> float {
        if (c1 == c2) return 0.0f;
        if (c1 < 0 || c2 < 0) return 1e6f;
        float sum = 0.0f;
        for (int d = 0; d < MATRIX_DIM; ++d) {
            float diff = centroids[c1][d] - centroids[c2][d];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    };

    // Apply Gaussian modulation
    float inv_2sigma2 = 1.0f / (2.0f * sigma * sigma);
    for (int i = 0; i < MATRIX_DIM; ++i) {
        int ci = clusters.cluster_labels[i];
        for (int j = 0; j < MATRIX_DIM; ++j) {
            int cj = clusters.cluster_labels[j];
            float d = cluster_dist(ci, cj);
            float gauss = std::exp(-d * d * inv_2sigma2);
            mat.data[i * MATRIX_DIM + j] *= (1.0f - alpha * gauss);
        }
    }
}

} // namespace scm
