// coarse_init.cpp — Coarse pocket-scan seeding for gen-0 (no-seed docking)
//
// Apache-2.0 © 2026 Le Bonhomme Pharma / NRGlab, Université de Montréal

#include "coarse_init.h"
#include "sampling_coverage.h"
#include "Vcontacts.h"   // ic2cf declaration

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <tuple>
#include <unordered_set>
#include <vector>

// ── Helpers ──────────────────────────────────────────────────────────────────

static inline double unit_fraction(double value) {
    return value - std::floor(value);
}

static int prime_for_dimension(int dimension) {
    int found = 0;
    for (int candidate = 2; ; ++candidate) {
        bool prime = true;
        for (int divisor = 2; divisor * divisor <= candidate; ++divisor) {
            if (candidate % divisor == 0) {
                prime = false;
                break;
            }
        }
        if (prime && found++ == dimension) return candidate;
    }
}

static double radical_inverse(unsigned long long index, int base) {
    double result = 0.0;
    double factor = 1.0 / static_cast<double>(base);
    while (index > 0) {
        result += static_cast<double>(index % static_cast<unsigned long long>(base)) * factor;
        index /= static_cast<unsigned long long>(base);
        factor /= static_cast<double>(base);
    }
    return result;
}

static double low_discrepancy_ic(const FA_Global* FA,
                                 const genlim* gene_lim,
                                 int gene_index,
                                 int prime_base,
                                 unsigned long long sample_index,
                                 double shift) {
    const double u = unit_fraction(
        radical_inverse(sample_index, prime_base) + shift);

    // The first rotational IC is the polar angle of the ligand GPA frame.
    // Uniform theta would over-sample the poles; uniform cos(theta), combined
    // with the two linearly sampled azimuth/roll ICs, gives isotropic rigid-body
    // coverage without using the native orientation.
    if (gene_index == 1 && FA->map_par && FA->map_par[gene_index].typ == 1) {
        constexpr double kRadToDeg = 57.2957795130823208768;
        const double theta = std::acos(std::clamp(1.0 - 2.0 * u, -1.0, 1.0)) * kRadToDeg;
        return std::clamp(theta, gene_lim[gene_index].min, gene_lim[gene_index].max);
    }

    return gene_lim[gene_index].min
           + u * (gene_lim[gene_index].max - gene_lim[gene_index].min);
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
    // Flag-gated coverage boost (FLEXAIDDS_SAMPLE_COVERAGE_BOOST=1, default OFF).
    const auto cov = flexaids::sampling::apply_coverage_boost(
        FA->coarse_init_n_orient, FA->coarse_init_n_seeds);
    const int n_orient  = cov.n_orient;
    const int n_seeds   = cov.n_seeds;
    if (cov.boost_applied) {
        printf("[COARSE-INIT] SAMPLE_COVERAGE_BOOST on → n_orient=%d n_seeds=%d\n",
               n_orient, n_seeds);
    }
    const float step    = (FA->coarse_init_grid_step > 0.0f)
                          ? FA->coarse_init_grid_step : 3.0f;
    const float thresh2 = (step * 1.5f) * (step * 1.5f);

    // ── 1. Build candidate grid-point index set ───────────────────────────
    // Seed from reflig_nearest_grid, then expand one neighbourhood shell.
    std::unordered_set<int> cand_set;

    const bool has_reflig_candidates =
        FA->reflig_nearest_count > 0 && FA->reflig_nearest_grid;
    if (has_reflig_candidates) {
        for (int k = 0; k < FA->reflig_nearest_count; k++) {
            int idx = FA->reflig_nearest_grid[k];
            if (idx >= 1 && idx < FA->num_grd)
                cand_set.insert(idx);
        }
    } else {
        // No reference ligand is allowed in a defined-cleft benchmark. Cover
        // the complete cleft deterministically with one representative per
        // coarse spatial voxel instead of taking the first 50 storage-order
        // points, which can miss most of an elongated cavity.
        const float voxel_step = std::max(step, 0.5f);
        std::map<std::tuple<int,int,int>, int> voxel_representatives;
        for (int gi = 1; gi < FA->num_grd; ++gi) {
            const auto key = std::make_tuple(
                static_cast<int>(std::floor(cleftgrid[gi].coor[0] / voxel_step)),
                static_cast<int>(std::floor(cleftgrid[gi].coor[1] / voxel_step)),
                static_cast<int>(std::floor(cleftgrid[gi].coor[2] / voxel_step)));
            voxel_representatives.emplace(key, gi);
        }
        for (const auto& [key, gi] : voxel_representatives) {
            (void)key;
            cand_set.insert(gi);
        }
    }

    // Expand: for each seed, add all grid points within thresh2
    if (has_reflig_candidates) {
        std::vector<int> seeds_snapshot(cand_set.begin(), cand_set.end());
        std::sort(seeds_snapshot.begin(), seeds_snapshot.end());
        for (int si : seeds_snapshot) {
            for (int gi = 1; gi < FA->num_grd; gi++) {
                if (cand_set.count(gi)) continue;
                if (dist_sq(cleftgrid[si], cleftgrid[gi]) <= thresh2)
                    cand_set.insert(gi);
            }
        }
    }

    std::vector<int> candidates(cand_set.begin(), cand_set.end());
    std::sort(candidates.begin(), candidates.end());
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
    // Score against the live buffers because ic2cf expects the complete FlexAID
    // atom/residue graph. Normal ic2cf exits intentionally leave Cartesian atom
    // state materialised for pose output, so preserve the pre-scan state and
    // restore it before the actual population is evaluated. Without this reset,
    // gen-0 decodes the screened genes against the final trial pose and the
    // accepted score is not reproducible.
    const std::vector<atom> atom_baseline(
        atoms, atoms + static_cast<std::size_t>(FA->atm_cnt) + 1);
    const std::vector<resid> residue_baseline(
        residue, residue + static_cast<std::size_t>(FA->res_cnt) + 1);
    std::vector<OptRes> optres_baseline;
    if (FA->optres && FA->num_optres > 0) {
        optres_baseline.assign(FA->optres, FA->optres + FA->num_optres);
    }
    const int numcarec_baseline = VC->numcarec;
    const auto restore_scan_state = [&]() {
        std::copy(atom_baseline.begin(), atom_baseline.end(), atoms);
        std::copy(residue_baseline.begin(), residue_baseline.end(), residue);
        if (!optres_baseline.empty()) {
            std::copy(optres_baseline.begin(), optres_baseline.end(), FA->optres);
        }
        VC->numcarec = numcarec_baseline;
    };
    struct ScanResult {
        int    grid_idx;
        double cf_val;
        // IC values for genes 1..n_genes-1
        std::vector<double> ics;
    };
    std::vector<ScanResult> results;
    results.reserve(static_cast<std::size_t>(candidates.size()) * n_orient);

    double icv[MAX_NUM_GENES] = {};
    std::vector<double> sequence_shifts(static_cast<std::size_t>(n_genes), 0.0);
    std::vector<int> sequence_bases(static_cast<std::size_t>(n_genes), 2);
    for (int g = 1; g < n_genes; ++g) {
        sequence_bases[static_cast<std::size_t>(g)] = prime_for_dimension(g - 1);
        sequence_shifts[static_cast<std::size_t>(g)] =
            static_cast<double>(dice()) / (2147483647.0 + 1.0);
    }

    for (std::size_t candidate_pos = 0; candidate_pos < candidates.size(); ++candidate_pos) {
        const int ci = candidates[candidate_pos];
        // Gene 0: grid index (translation)
        const double grid_ic = static_cast<double>(ci);
        if (grid_ic < gene_lim[0].min || grid_ic > gene_lim[0].max)
            continue;
        icv[0] = grid_ic;

        for (int oi = 0; oi < n_orient; oi++) {
            restore_scan_state();
            // Genes 1..N-1: randomized low-discrepancy IC coverage. The
            // per-restart shifts come from the seeded RNG, but the sequence
            // itself stratifies every dimension instead of relying on a joint
            // random draw whose coverage collapses as torsions are added.
            const unsigned long long sample_index =
                static_cast<unsigned long long>(candidate_pos)
                    * static_cast<unsigned long long>(n_orient)
                + static_cast<unsigned long long>(oi) + 1ULL;
            for (int g = 1; g < n_genes; g++) {
                icv[g] = low_discrepancy_ic(
                    FA, gene_lim, g,
                    sequence_bases[static_cast<std::size_t>(g)], sample_index,
                    sequence_shifts[static_cast<std::size_t>(g)]);
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
            // Rank with the same apparent CF used by the GA. Sorting on the
            // attractive contact component alone preferentially selected
            // interpenetrating poses with very negative cf.com but a dominant
            // positive steric wall; those seeds became poor gen-0 chromosomes.
            r.cf_val = get_apparent_cf_evalue(&cf);
            if (!std::isfinite(r.cf_val)) continue;
            r.ics.assign(icv + 1, icv + n_genes);
            results.push_back(std::move(r));
        }
    }

    restore_scan_state();

    if (results.empty()) {
        fprintf(stderr, "[COARSE-INIT] WARN: all evaluations failed, skipping seed injection\n");
        return;
    }

    // ── 4. Rank by CF and select seeds (shared pure helper) ───────────────
    // Absolute CF < 0 is NOT a valid "has contacts" test: SAS weight /
    // FLEXAIDDS_POLAR_DESOLV_WEIGHT can shift the zero-point positive while the
    // placement still has real VCT contacts. Hard-clash sentinels
    // (CF ≥ CLASH_THRESHOLD) are dropped unless FLEXAIDDS_COARSE_INIT_FORCE_RANKED=1
    // (1M2Z-class: inject least-bad ranks so gen-0 is not empty).
    std::vector<double> cf_vals;
    cf_vals.reserve(results.size());
    for (const auto& r : results) cf_vals.push_back(r.cf_val);

    const bool force_ranked = flexaids::sampling::force_ranked_seeds_enabled();
    auto keep_idx = flexaids::sampling::select_ranked_seed_indices(
        cf_vals, n_seeds, static_cast<double>(CLASH_THRESHOLD), force_ranked);

    // Optional spatial diversity when coverage boost is on (spread seeds across
    // the cleft rather than stacking near-identical grid points).
    if (cov.boost_applied && !keep_idx.empty()) {
        std::vector<std::array<double, 3>> coords(results.size());
        for (std::size_t i = 0; i < results.size(); ++i) {
            const int gi = results[i].grid_idx;
            if (gi >= 1 && gi < FA->num_grd) {
                coords[i] = {cleftgrid[gi].coor[0],
                             cleftgrid[gi].coor[1],
                             cleftgrid[gi].coor[2]};
            }
        }
        const auto ranked_all = flexaids::sampling::rank_indices_by_cf_asc(cf_vals);
        // Prefer diverse picks among the non-clash (or force-ranked) pool.
        keep_idx = flexaids::sampling::diversify_by_min_distance(
            keep_idx.empty() ? ranked_all : keep_idx,
            &coords, n_seeds, /*min_dist_A=*/step);
        // Re-filter clash unless force_ranked.
        if (!force_ranked) {
            std::vector<std::size_t> filtered;
            for (std::size_t i : keep_idx) {
                if (cf_vals[i] < static_cast<double>(CLASH_THRESHOLD))
                    filtered.push_back(i);
            }
            keep_idx.swap(filtered);
        }
    }

    if (keep_idx.empty()) {
        printf("[COARSE-INIT] No non-clash placements found (best CF=%.2f ≥ CLASH_THRESHOLD), "
               "skipping seed injection%s\n",
               results[flexaids::sampling::rank_indices_by_cf_asc(cf_vals).front()].cf_val,
               force_ranked ? " (FORCE_RANKED on but no finite CF)" : " (set FLEXAIDDS_COARSE_INIT_FORCE_RANKED=1 to inject least-bad ranks)");
        return;
    }

    if (force_ranked && cf_vals[keep_idx.front()] >= static_cast<double>(CLASH_THRESHOLD)) {
        printf("[COARSE-INIT] FORCE_RANKED: injecting %zu clash-scale seeds "
               "(best CF=%.2f) so gen-0 is non-empty\n",
               keep_idx.size(), cf_vals[keep_idx.front()]);
    }

    const int keep_n = static_cast<int>(keep_idx.size());

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
        const auto& r = results[keep_idx[static_cast<std::size_t>(k)]];
        FA->coarse_seeds_grid[k] = r.grid_idx;
        const int ic_stride = k * (n_genes - 1);
        for (int g = 0; g < n_genes - 1; g++)
            FA->coarse_seeds_genes[ic_stride + g] = static_cast<float>(r.ics[static_cast<std::size_t>(g)]);
    }
    FA->coarse_seeds_count = keep_n;

    printf("[COARSE-INIT] %d ranked seeds ready "
           "(best CF=%.2f, worst CF=%.2f; absolute CF sign not required)\n",
           keep_n,
           cf_vals[keep_idx.front()],
           cf_vals[keep_idx.back()]);
}
