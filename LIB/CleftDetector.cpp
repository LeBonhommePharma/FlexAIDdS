#include "CleftDetector.h"
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

    // Collect protein (non-HET) atom indices that have coordinates
    std::vector<int> idx;
    idx.reserve(atm_cnt);
    for (int i = 0; i < atm_cnt; ++i) {
        // skip atoms with zero coordinates (uninitialised/padding)
        if (atoms[i].coor[0] == 0.0f &&
            atoms[i].coor[1] == 0.0f &&
            atoms[i].coor[2] == 0.0f &&
            atoms[i].radius  == 0.0f) continue;
        idx.push_back(i);
    }

    const int n = static_cast<int>(idx.size());

#ifdef _OPENMP
    // Each thread collects into a local vector, merged later
    #pragma omp parallel
    {
        std::vector<Probe> local;
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
                if (keep) local.push_back(pr);
            }
        }
        #pragma omp critical
        probes.insert(probes.end(), local.begin(), local.end());
    }
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

    // ── Cognate-aware coverage ──────────────────────────────────────────────
    // SURFNET historically returned ONLY the largest gap-sphere cluster as the
    // binding cleft.  For ~40% of Astex Diverse the largest cavity is NOT the
    // cognate ligand site, so those complexes were starved of grid points at the
    // true pocket and scored 0 successes.  Instead, keep EVERY cluster meeting
    // min_cluster_size (all genuine pockets, not just the biggest).  The
    // downstream site-confinement filter in top.cpp restricts the final grid to
    // the cognate centroid sphere, so handing it all pockets simply guarantees
    // the right one is present.  Tiny (< min_cluster_size) singleton/noise
    // clusters are still dropped to keep the grid build bounded.
    std::set<int> kept_labels;
    int kept_clusters = 0, kept_spheres = 0;
    for (auto& kv : freq) {
        if (kv.second >= params.min_cluster_size) {
            kept_labels.insert(kv.first);
            ++kept_clusters;
            kept_spheres += kv.second;
        }
    }

    if (kept_labels.empty()) {
        // No cluster reached min_cluster_size — fall back to the single largest
        // so detection still yields a usable cleft rather than nothing.
        fprintf(stderr, "CleftDetector WARNING: no cluster reached min_cluster_size "
                "(largest has %d, min is %d) — falling back to largest cluster\n",
                best_count, params.min_cluster_size);
        if (best_label >= 0) {
            kept_labels.insert(best_label);
            kept_clusters = 1;
            kept_spheres  = best_count;
        }
    }

    printf("CleftDetector: keeping %d cluster(s) >= %d spheres (%d spheres total; "
           "largest cluster %d)\n",
           kept_clusters, params.min_cluster_size, kept_spheres, best_count);

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
