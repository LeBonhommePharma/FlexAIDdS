// coarse_init.cpp — Coarse pocket-scan seeding for gen-0 (autonomous docking)
//
// Apache-2.0 © 2026 Le Bonhomme Pharma / NRGlab, Université de Montréal

#include "coarse_init.h"
#include "Vcontacts.h"   // ic2cf declaration

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <unordered_set>
#include <vector>

// ── Helpers ──────────────────────────────────────────────────────────────────

static inline double rand_in_range(double lo, double hi,
                                   std::function<int32_t()>& dice) {
    double u = static_cast<double>(dice()) / (2147483647.0 + 1.0); // [0,1)
    return lo + u * (hi - lo);
}

// Euclidean distance squared between two gridpoints.
static inline float dist_sq(const gridpoint& a, const gridpoint& b) {
    float dx = a.coor[0] - b.coor[0];
    float dy = a.coor[1] - b.coor[1];
    float dz = a.coor[2] - b.coor[2];
    return dx*dx + dy*dy + dz*dz;
}

// ── Public API ────────────────────────────────────────────────────────────────

void run_coarse_pocket_scan(
    FA_Global*                 FA,
    VC_Global*                 VC,
    GB_Global*                 GB,
    atom*                      atoms,
    resid*                     residue,
    gridpoint*                 cleftgrid,
    const genlim*              gene_lim,
    std::function<int32_t()>&  dice)
{
    if (!FA || !VC || !GB || !atoms || !residue || !cleftgrid || !gene_lim)
        return;
    if (FA->num_grd <= 1 || GB->num_genes <= 0)
        return;

    const int n_genes   = GB->num_genes;
    const int n_orient  = std::max(1, FA->coarse_init_n_orient);
    const int n_seeds   = std::max(1, FA->coarse_init_n_seeds);
    const float step    = (FA->coarse_init_grid_step > 0.0f)
                          ? FA->coarse_init_grid_step : 3.0f;
    const float thresh2 = (step * 1.5f) * (step * 1.5f);

    // ── 1. Build candidate grid-point index set ───────────────────────────
    // Seed from reflig_nearest_grid, then expand one neighbourhood shell.
    std::unordered_set<int> cand_set;

    if (FA->reflig_nearest_count > 0 && FA->reflig_nearest_grid) {
        for (int k = 0; k < FA->reflig_nearest_count; k++) {
            int idx = FA->reflig_nearest_grid[k];
            if (idx >= 1 && idx < FA->num_grd)
                cand_set.insert(idx);
        }
    } else {
        // Fallback: use the first MIN(50, num_grd-1) grid points
        int fallback_n = std::min(50, FA->num_grd - 1);
        for (int k = 1; k <= fallback_n; k++)
            cand_set.insert(k);
    }

    // Expand: for each seed, add all grid points within thresh2
    std::vector<int> seeds_snapshot(cand_set.begin(), cand_set.end());
    for (int si : seeds_snapshot) {
        for (int gi = 1; gi < FA->num_grd; gi++) {
            if (cand_set.count(gi)) continue;
            if (dist_sq(cleftgrid[si], cleftgrid[gi]) <= thresh2)
                cand_set.insert(gi);
        }
    }

    std::vector<int> candidates(cand_set.begin(), cand_set.end());
    if (candidates.empty()) {
        fprintf(stderr, "[COARSE-INIT] WARN: empty candidate set, skipping scan\n");
        return;
    }

    printf("[COARSE-INIT] Scanning %d candidate grid points × %d orientations\n",
           static_cast<int>(candidates.size()), n_orient);

    // ── 2. Disable native-anchor seeding so coarse seeds can win gen-0 ──────
    // Without this, gaboom's populate_chromosomes forces gene[0]=0 (the blinded
    // anchor) for ~90% of the population via the native_direct_seed path, completely
    // overriding the coarse seeds. ic2cf already saves/restores FA->ori and the
    // atoms it modifies, so we can call it on the real buffers directly.
    FA->reflig_seed_fraction = 0.0f;

    // ── 3. Evaluate each (candidate, orientation) pair ───────────────────────
    // IMPORTANT: ligand atoms in FlexAID live at special high-index slots
    // (map_par[i].atm ≈ 90001+), far outside a natm-sized copy.  Using a
    // shallow working copy caused out-of-bounds reads producing false CF values.
    // ic2cf already saves/restores FA->ori (and the atoms it moves via saved_atoms
    // on early-exit paths); for normal exits atoms are left modified but the
    // next ic2cf call rebuilds from IC values, so using the real buffers is safe.
    struct ScanResult {
        int    grid_idx;
        double cf_val;
        // IC values for genes 1..n_genes-1
        std::vector<double> ics;
    };
    std::vector<ScanResult> results;
    results.reserve(static_cast<std::size_t>(candidates.size()) * n_orient);

    double icv[MAX_NUM_GENES] = {};

    for (int ci : candidates) {
        // Gene 0: grid index (translation)
        const double grid_ic = static_cast<double>(ci);
        if (grid_ic < gene_lim[0].min || grid_ic > gene_lim[0].max)
            continue;
        icv[0] = grid_ic;

        for (int oi = 0; oi < n_orient; oi++) {
            // Genes 1..N-1: random IC in [min, max]
            for (int g = 1; g < n_genes; g++) {
                icv[g] = rand_in_range(gene_lim[g].min, gene_lim[g].max, dice);
            }

            cfstr cf{};
            try {
                cf = ic2cf(FA, VC, atoms, residue, cleftgrid, n_genes, icv);
            } catch (...) {
                // Scoring errors are non-fatal — skip this orientation.
                continue;
            }

            // Capture the IC values (genes 1..N-1) for storage.
            ScanResult r;
            r.grid_idx = ci;
            r.cf_val   = cf.com;
            r.ics.assign(icv + 1, icv + n_genes);
            results.push_back(std::move(r));
        }
    }

    if (results.empty()) {
        fprintf(stderr, "[COARSE-INIT] WARN: all evaluations failed, skipping seed injection\n");
        return;
    }

    // ── 4. Sort by CF ascending (most negative = best contacts first) ─────
    std::partial_sort(results.begin(),
                      results.begin() + std::min(n_seeds, static_cast<int>(results.size())),
                      results.end(),
                      [](const ScanResult& a, const ScanResult& b) {
                          return a.cf_val < b.cf_val;
                      });

    const int actual_n = std::min(n_seeds, static_cast<int>(results.size()));

    // Filter: only keep results with CF < 0 (at least some contacts).
    // CF≥0 means the placement has no binding contacts — not useful as a seed.
    int keep_n = 0;
    for (int k = 0; k < actual_n; k++) {
        if (results[static_cast<std::size_t>(k)].cf_val < 0.0)
            keep_n++;
        else
            break; // sorted, so all subsequent are ≥0
    }

    if (keep_n == 0) {
        printf("[COARSE-INIT] No contact-forming placements found (best CF=%.2f), "
               "skipping seed injection\n",
               results[0].cf_val);
        return;
    }

    // ── 5. Allocate and fill FA coarse-seed arrays ────────────────────────
    free(FA->coarse_seeds_grid);
    free(FA->coarse_seeds_genes);
    FA->coarse_seeds_grid  = nullptr;
    FA->coarse_seeds_genes = nullptr;
    FA->coarse_seeds_count = 0;

    FA->coarse_seeds_grid  = static_cast<int*>  (malloc(keep_n * sizeof(int)));
    FA->coarse_seeds_genes = static_cast<float*>(malloc(keep_n * (n_genes - 1) * sizeof(float)));

    if (!FA->coarse_seeds_grid || !FA->coarse_seeds_genes) {
        free(FA->coarse_seeds_grid);
        free(FA->coarse_seeds_genes);
        FA->coarse_seeds_grid  = nullptr;
        FA->coarse_seeds_genes = nullptr;
        fprintf(stderr, "[COARSE-INIT] ERROR: allocation failed, skipping seed injection\n");
        return;
    }

    for (int k = 0; k < keep_n; k++) {
        const auto& r = results[static_cast<std::size_t>(k)];
        FA->coarse_seeds_grid[k] = r.grid_idx;
        const int ic_stride = k * (n_genes - 1);
        for (int g = 0; g < n_genes - 1; g++)
            FA->coarse_seeds_genes[ic_stride + g] = static_cast<float>(r.ics[static_cast<std::size_t>(g)]);
    }
    FA->coarse_seeds_count = keep_n;

    printf("[COARSE-INIT] %d contact-forming seeds ready "
           "(best CF=%.2f, worst CF=%.2f)\n",
           keep_n,
           results[0].cf_val,
           results[static_cast<std::size_t>(keep_n - 1)].cf_val);
}
