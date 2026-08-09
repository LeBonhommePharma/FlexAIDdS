#include "CleftDetector.h"
#include "ensemble_pipeline.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * SURFNET / GetCleft gap-sphere algorithm
 * ----------------------------------------
 * For every pair of surface atoms (i, j) within max_pair_dist:
 *   - place a probe sphere centred at the midpoint
 *   - set its radius = half the inter-atom distance
 *   - shrink until no other atom k (k != i, k != j) overlaps
 *   - keep if radius >= probe_radius_min
 *
 * The surviving probes are clustered (single-linkage) and the
 * largest cluster is returned as the binding cleft.
 */

// ── helpers ──────────────────────────────────────────────────────────────

static inline float sq(float x) { return x * x; }

static float sqdist3(const float a[3], const float b[3]) {
    return sq(a[0] - b[0]) + sq(a[1] - b[1]) + sq(a[2] - b[2]);
}

// ── probe generation ────────────────────────────────────────────────────

struct Probe { float center[3]; float radius; };

static std::vector<Probe> generate_probes(
    const atom* atoms, int atm_cnt,
    const CleftDetectorParams& p)
{
    const float max_pair_sq = p.max_pair_dist * p.max_pair_dist;
    std::vector<Probe> probes;

    // Collect protein (non-HET) atom indices that have coordinates.
    // If oracle_radius > 0, pre-filter to atoms within that radius of
    // oracle_center.  This eliminates the O(N^3) blowup for multimeric
    // receptors: 1OF6 (20826 atoms, 8 chains) drops to ~200-400 atoms
    // within 15 A of the cognate site, giving a ~10,000x speedup.
    // The downstream site-confinement step in top.cpp trims the grid to
    // the cognate centroid sphere anyway, so excluding far-away atoms
    // from SURFNET is correct and safe.
    const bool has_spatial_filter = (p.oracle_radius > 0.0f);
    const float filter_r2 = has_spatial_filter
        ? p.oracle_radius * p.oracle_radius : 0.0f;

    std::vector<int> idx;
    idx.reserve(atm_cnt);
    for (int i = 0; i < atm_cnt; ++i) {
        // skip atoms with zero coordinates (uninitialised/padding)
        if (atoms[i].coor[0] == 0.0f &&
            atoms[i].coor[1] == 0.0f &&
            atoms[i].coor[2] == 0.0f &&
            atoms[i].radius  == 0.0f) continue;
        // spatial pre-filter: skip atoms outside oracle sphere
        if (has_spatial_filter) {
            float dx = atoms[i].coor[0] - p.oracle_center[0];
            float dy = atoms[i].coor[1] - p.oracle_center[1];
            float dz = atoms[i].coor[2] - p.oracle_center[2];
            if (dx*dx + dy*dy + dz*dz > filter_r2) continue;
        }
        idx.push_back(i);
    }
    if (has_spatial_filter) {
        printf("CleftDetector: oracle spatial filter (%.1f A) reduced atoms %d -> %d\n",
               p.oracle_radius, atm_cnt, static_cast<int>(idx.size()));
    }

    const int n = static_cast<int>(idx.size());

#ifdef _OPENMP
    // Per-iteration buckets, concatenated in ascending ii order after the
    // parallel region. The merge order must not depend on thread arrival:
    // probe order fixes generate_grid()'s cleftgrid index assignment, and
    // that index is GA gene 0, so a thread-arrival-ordered merge decorrelates
    // the whole search run-to-run and across thread counts (see
    // determinism_cleft/FINDING_cleft_grid_nondeterminism.md). Bucketing by
    // ii reproduces the serial (ascending ii, ascending jj) order bit-exactly
    // at any thread count and any schedule.
    std::vector<std::vector<Probe>> buckets(static_cast<size_t>(n));
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 64) nowait
        for (int ii = 0; ii < n; ++ii) {
            int i = idx[ii];
            for (int jj = ii + 1; jj < n; ++jj) {
                int j = idx[jj];
                float d2 = sqdist3(atoms[i].coor, atoms[j].coor);
                if (d2 > max_pair_sq || d2 < 1.0f) continue;

                float d = std::sqrt(d2);
                Probe pr;
                pr.center[0] = 0.5f * (atoms[i].coor[0] + atoms[j].coor[0]);
                pr.center[1] = 0.5f * (atoms[i].coor[1] + atoms[j].coor[1]);
                pr.center[2] = 0.5f * (atoms[i].coor[2] + atoms[j].coor[2]);
                pr.radius    = 0.5f * d;

                // Clamp initial radius
                if (pr.radius > p.probe_radius_max)
                    pr.radius = p.probe_radius_max;

                // Shrink until no other atom overlaps (atom radius + probe radius)
                bool keep = true;
                while (pr.radius >= p.probe_radius_min) {
                    bool clash = false;
                    for (int kk = 0; kk < n && !clash; ++kk) {
                        int k = idx[kk];
                        if (k == i || k == j) continue;
                        float dk2 = sqdist3(pr.center, atoms[k].coor);
                        float overlap = atoms[k].radius + pr.radius;
                        if (dk2 < overlap * overlap)
                            clash = true;
                    }
                    if (!clash) break;
                    pr.radius -= p.probe_shrink_step;
                }
                if (pr.radius < p.probe_radius_min) keep = false;
                if (keep) buckets[ii].push_back(pr);
            }
        }
    }
    for (int ii = 0; ii < n; ++ii)
        probes.insert(probes.end(), buckets[ii].begin(), buckets[ii].end());
#else
    for (int ii = 0; ii < n; ++ii) {
        int i = idx[ii];
        for (int jj = ii + 1; jj < n; ++jj) {
            int j = idx[jj];
            float d2 = sqdist3(atoms[i].coor, atoms[j].coor);
            if (d2 > max_pair_sq || d2 < 1.0f) continue;

            float d = std::sqrt(d2);
            Probe pr;
            pr.center[0] = 0.5f * (atoms[i].coor[0] + atoms[j].coor[0]);
            pr.center[1] = 0.5f * (atoms[i].coor[1] + atoms[j].coor[1]);
            pr.center[2] = 0.5f * (atoms[i].coor[2] + atoms[j].coor[2]);
            pr.radius    = 0.5f * d;

            if (pr.radius > p.probe_radius_max)
                pr.radius = p.probe_radius_max;

            bool keep = true;
            while (pr.radius >= p.probe_radius_min) {
                bool clash = false;
                for (int kk = 0; kk < n && !clash; ++kk) {
                    int k = idx[kk];
                    if (k == i || k == j) continue;
                    float dk2 = sqdist3(pr.center, atoms[k].coor);
                    float overlap = atoms[k].radius + pr.radius;
                    if (dk2 < overlap * overlap)
                        clash = true;
                }
                if (!clash) break;
                pr.radius -= p.probe_shrink_step;
            }
            if (pr.radius < p.probe_radius_min) keep = false;
            if (keep) probes.push_back(pr);
        }
    }
#endif

    // CANONICAL GEOMETRIC ORDER (opt-in). Since the ascending-ii bucket merge
    // above, probe order is already deterministic and thread-count-invariant
    // (identical to the serial branch), so this sort is NOT needed for
    // reproducibility. It remains as an opt-in alternative canonical order:
    // it reorders probe tie-breaks relative to the serial order, which changes
    // poses (1G9V 5.13 -> 7.54 A single-thread) and would invalidate prior
    // benchmark numbers if flipped on. See CLEFTSORT_AB_VERDICT.md.
    {
        const char* cs = std::getenv("FLEXAIDDS_CLEFT_SORT");
        if (cs && cs[0] == '1') {
            std::sort(probes.begin(), probes.end(),
                      [](const Probe& a, const Probe& b) {
                          if (a.center[0] != b.center[0]) return a.center[0] < b.center[0];
                          if (a.center[1] != b.center[1]) return a.center[1] < b.center[1];
                          if (a.center[2] != b.center[2]) return a.center[2] < b.center[2];
                          return a.radius < b.radius;
                      });
        }
    }

    return probes;
}

// ── single-linkage clustering ───────────────────────────────────────────

// Union-Find
static int uf_find(std::vector<int>& parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}
static void uf_union(std::vector<int>& parent, std::vector<int>& rank, int a, int b) {
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a == b) return;
    if (rank[a] < rank[b]) std::swap(a, b);
    parent[b] = a;
    if (rank[a] == rank[b]) ++rank[a];
}

static std::vector<int> cluster_probes(const std::vector<Probe>& probes, float cutoff) {
    int n = static_cast<int>(probes.size());
    std::vector<int> parent(n), rank(n, 0);
    std::iota(parent.begin(), parent.end(), 0);

    float cutoff_sq = cutoff * cutoff;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (sqdist3(probes[i].center, probes[j].center) < cutoff_sq)
                uf_union(parent, rank, i, j);

    // canonical labels
    for (int i = 0; i < n; ++i) parent[i] = uf_find(parent, i);
    return parent;
}

// ── public API ──────────────────────────────────────────────────────────

sphere* detect_cleft(const atom* atoms, const resid* /*residue*/,
                     int atm_cnt, int /*res_cnt*/,
                     const CleftDetectorParams& params)
{
    printf("CleftDetector: scanning %d atoms for binding cavities ...\n", atm_cnt);

    // 1. generate gap-spheres
    std::vector<Probe> probes = generate_probes(atoms, atm_cnt, params);
    printf("CleftDetector: %d gap-spheres survived shrinking\n",
           static_cast<int>(probes.size()));

    if (probes.empty()) {
        fprintf(stderr, "CleftDetector WARNING: no cavities found — "
                "try increasing max_pair_dist or decreasing probe_radius_min\n");
        return nullptr;
    }

    // 2. cluster
    std::vector<int> labels = cluster_probes(probes, params.cluster_cutoff);

    // tally cluster sizes
    std::map<int, int> freq;
    for (int l : labels) freq[l]++;

    int best_label = -1, best_count = 0;
    for (auto& kv : freq) {
        if (kv.second > best_count) { best_label = kv.first; best_count = kv.second; }
    }

    // ── Layer 2: ligandable top-K (not largest void only) ───────────────────
    // Score each cluster meeting min_cluster_size by volume×enclosure proxy
    // (ensemble::ligandable_score). Keep top_k_clefts; fallback to largest.
    //
    // Cleft-ranking audit sidecar (FLEXAIDDS_CLEFT_DUMP=<path>): write every probe
    // sphere of every scored cluster as label/x/y/z/radius TSV. Summary statistics
    // cannot answer "is the cognate cleft in the candidate set" -- a bounding-box
    // diagonal is not a pocket radius (observed 1.8-92.5 A, some larger than the
    // protein) and a centroid can sit far outside an elongated cleft that does
    // enclose the ligand. The only sound test is min sphere-to-ligand-atom distance,
    // which needs the spheres themselves. Diagnostic only: opt-in via env, writes to
    // a separate file, and feeds nothing back into scoring or selection.
    // NOTE: detect_clefts() runs MORE THAN ONCE per target (per receptor conformer /
    // restart) and cluster labels are not stable across invocations. Truncating here
    // would leave only the last pass, so the sidecar would disagree with the
    // accumulated log (observed: 32 scored vs 32 dumped labels, only 5 shared).
    // Append instead, and stamp each row with an invocation id so an analysis can
    // pair spheres with the matching log block.
    static int cleft_dump_invocation = 0;
    const int this_invocation = cleft_dump_invocation++;
    std::FILE* cleft_dump = nullptr;
    if (const char* dump_path = std::getenv("FLEXAIDDS_CLEFT_DUMP")) {
        if (dump_path[0] != '\0') {
            const bool fresh = (this_invocation == 0);
            cleft_dump = std::fopen(dump_path, fresh ? "w" : "a");
            if (cleft_dump) {
                if (fresh)
                    std::fprintf(cleft_dump,
                                 "invocation\tlabel\tligandable_score\tx\ty\tz\tradius\n");
            } else {
                std::fprintf(stderr,
                             "WARNING: [CLEFT_DUMP] could not open %s for write\n",
                             dump_path);
            }
        }
    }
    std::vector<std::pair<int, double>> scored;
    scored.reserve(freq.size());
    for (auto& kv : freq) {
        if (kv.second < params.min_cluster_size) continue;
        double sum_r = 0.0, minx = 1e30, miny = 1e30, minz = 1e30;
        double maxx = -1e30, maxy = -1e30, maxz = -1e30;
        // Cleft-ranking audit: accumulate the cluster centroid so the log can be
        // matched against the crystal-ligand centroid offline. Diagnostic only --
        // sum_x/y/z feed nothing but the printf below, so scoring, cleft selection
        // and CF remain bit-identical.
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        int n = 0;
        for (int i = 0; i < static_cast<int>(probes.size()); ++i) {
            if (labels[i] != kv.first) continue;
            sum_r += probes[i].radius;
            sum_x += (double)probes[i].center[0];
            sum_y += (double)probes[i].center[1];
            sum_z += (double)probes[i].center[2];
            minx = std::min(minx, (double)probes[i].center[0]);
            miny = std::min(miny, (double)probes[i].center[1]);
            minz = std::min(minz, (double)probes[i].center[2]);
            maxx = std::max(maxx, (double)probes[i].center[0]);
            maxy = std::max(maxy, (double)probes[i].center[1]);
            maxz = std::max(maxz, (double)probes[i].center[2]);
            ++n;
        }
        if (n <= 0) continue;
        const double mean_r = sum_r / static_cast<double>(n);
        const double dx = maxx - minx, dy = maxy - miny, dz = maxz - minz;
        const double bbox_diag = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double s = ensemble::ligandable_score(n, mean_r, bbox_diag);
        scored.emplace_back(kv.first, s);
        printf("CleftDetector: cluster label=%d n=%d ligandable_score=%.3f "
               "centroid=%.3f,%.3f,%.3f bbox_diag=%.3f mean_r=%.3f\n",
               kv.first, n, s,
               sum_x / static_cast<double>(n),
               sum_y / static_cast<double>(n),
               sum_z / static_cast<double>(n),
               bbox_diag, mean_r);
        // Emit this cluster's spheres WITH its ligandable_score, so the sidecar is
        // self-contained. Cluster labels are NOT comparable across engine processes
        // (the runner invokes the engine more than once per target and labels are
        // reassigned), so an analysis must never join the sidecar to stdout.log by
        // label -- doing so silently pairs unrelated clusters. Written here rather
        // than in the accumulation loop above because `s` is only known now.
        if (cleft_dump) {
            for (int i = 0; i < static_cast<int>(probes.size()); ++i) {
                if (labels[i] != kv.first) continue;
                std::fprintf(cleft_dump, "%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
                             this_invocation, kv.first, s,
                             (double)probes[i].center[0],
                             (double)probes[i].center[1],
                             (double)probes[i].center[2],
                             (double)probes[i].radius);
            }
        }
    }
    if (cleft_dump) {
        std::fclose(cleft_dump);
        cleft_dump = nullptr;
    }

    std::set<int> kept_labels;
    int kept_clusters = 0, kept_spheres = 0;
    if (scored.empty()) {
        fprintf(stderr, "CleftDetector WARNING: no cluster reached min_cluster_size "
                "(largest has %d, min is %d) — falling back to largest cluster\n",
                best_count, params.min_cluster_size);
        if (best_label >= 0) {
            kept_labels.insert(best_label);
            kept_clusters = 1;
            kept_spheres  = best_count;
        }
    } else {
        const int k = (params.top_k_clefts > 0)
            ? params.top_k_clefts
            : static_cast<int>(scored.size());
        const auto top = ensemble::select_top_k_clefts(scored, k);
        for (int lab : top) {
            kept_labels.insert(lab);
            kept_spheres += freq[lab];
        }
        kept_clusters = static_cast<int>(kept_labels.size());
    }

    printf("CleftDetector: keeping %d ligandable cluster(s) (top_k=%d, min_size=%d; "
           "%d spheres total; largest cluster %d)\n",
           kept_clusters, params.top_k_clefts, params.min_cluster_size,
           kept_spheres, best_count);

    // 3. build linked list (same format as read_spheres) from all kept clusters
    sphere* head = nullptr;
    for (int i = 0; i < static_cast<int>(probes.size()); ++i) {
        if (kept_labels.find(labels[i]) == kept_labels.end()) continue;
        sphere* s = (sphere*)malloc(sizeof(sphere));
        if (!s) { fprintf(stderr, "CleftDetector: out of memory\n"); break; }
        s->center[0] = probes[i].center[0];
        s->center[1] = probes[i].center[1];
        s->center[2] = probes[i].center[2];
        s->radius    = probes[i].radius;
        s->prev      = head;
        head          = s;
    }

    return head;
}

void write_cleft_spheres(const sphere* spheres, const char* filename) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "CleftDetector: cannot write %s\n", filename);
        return;
    }
    int n = 1;
    for (const sphere* s = spheres; s; s = s->prev, ++n) {
        fprintf(fp,
            "ATOM  %5d  C   SPH Z   1      %8.3f%8.3f%8.3f  1.00%6.2f\n",
            n, s->center[0], s->center[1], s->center[2], s->radius);
    }
    fclose(fp);
    printf("CleftDetector: wrote %d spheres to %s\n", n - 1, filename);
}

void free_sphere_list(sphere* head) {
    while (head) {
        sphere* tmp = head->prev;
        free(head);
        head = tmp;
    }
}

// ─── Task 8: Flexible Residue Selection (preprocessing only) ────────────────

std::vector<int> select_flexible_residues(
    const atom* atoms,
    const resid* residue,
    int atm_cnt,
    int res_cnt,
    const std::vector<int>& cleft_sphere_residues,
    double distance_shell_A,
    const std::vector<int>& active_site_residues,
    const std::vector<int>& user_fixed_residues,
    const std::vector<int>& user_forced_flexible)
{
    std::set<int> flexible;
    std::set<int> fixed(user_fixed_residues.begin(), user_fixed_residues.end());

    // 1. Add forced-flexible residues first (highest priority)
    for (int r : user_forced_flexible) {
        if (r >= 0 && r < res_cnt) {
            // Skip Gly/Ala only if we had backbone support — for now we allow
            // because the caller is responsible for having a backbone module.
            flexible.insert(r);
        }
    }

    // 2. Collect candidate residues near cleft spheres or active site
    auto add_nearby = [&](int res_index) {
        if (res_index < 0 || res_index >= res_cnt) return;
        if (fixed.count(res_index)) return; // respect fixed
        // Skip Gly/Ala only if we have structure data (residue != null)
        if (residue) {
            const char* resname = residue[res_index].name;
            if (std::strcmp(resname, "GLY") == 0 || std::strcmp(resname, "ALA") == 0) {
                return;
            }
        }
        flexible.insert(res_index);
    };

    // 3. Distance-based inclusion around cleft-related residues
    // (simplified: we treat the passed cleft_sphere_residues as seeds)
    for (int seed : cleft_sphere_residues) {
        add_nearby(seed);
    }
    for (int act : active_site_residues) {
        add_nearby(act);
    }

    // 4. Build deterministic sorted output
    std::vector<int> result(flexible.begin(), flexible.end());
    // Already sorted because we used std::set
    return result;
}
