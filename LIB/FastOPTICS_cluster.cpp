#include "FOPTICS.h"
#include "fast_optics.hpp"
#include "MinibatchSampler.h"
#include "ga_constants.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <vector>

namespace {

inline int clampi(int x, int lo, int hi) {
    return std::max(lo, std::min(hi, x));
}

// Ligand flexible dihedral count lives on the ligand residue (not FA_Global).
int fo_ligand_fdih(const FA_Global* FA)
{
    if (!FA || !FA->resligand) return 0;
    return (FA->resligand->fdih > 0) ? FA->resligand->fdih : 0;
}

// Effective dimensionality for pose density clustering (Sander et al. 1998).
// Rigid-body placement ≈ 6 continuous DoF; each flexible dihedral adds one.
// Cap at 20 so 2·dim stays inside Ankerst's practical MinPts band for typical
// drug-like ligands. IC gene count is a fallback when fdih unset.
int fo_effective_dim(const FA_Global* FA, const GB_Global* GB)
{
    const int fdih = fo_ligand_fdih(FA);
    const int npar = (FA && FA->npar > 0) ? FA->npar
                     : (GB && GB->num_genes > 0) ? GB->num_genes : 0;
    int dim = 6 + fdih;  // SE(3) + torsions
    if (fdih <= 0 && npar > 0)
        dim = std::max(2, npar);
    return clampi(dim, 2, 20);
}

// Single MinPts for one FastOPTICS pass (production: run the algorithm once).
//
//   1. Sander et al. 1998: MinPts ≈ 2 · dim  (dim > 2)
//   2. Ankerst et al. 1999: MinPts ∈ [10, 20] “always good results”
//   3. Ester et al. 1996:   MinPts = 4 floor for low-d / small-N
//   4. Ensemble feasibility: MinPts < N and MinPts ≤ N/3 so cores can form
//
// CF diversity only softens MinPts toward Ester's floor on near-degenerate
// landscapes — it does not replace Sander/Ankerst. Do not re-run FO at multiple
// MinPts (legacy triple ladder was a testing-only artifact from the old repo).
int fo_choose_minpts(const FA_Global* FA, const GB_Global* GB,
                     int nChrom, double diversity_ratio, int* dim_out)
{
    const int dim = fo_effective_dim(FA, GB);
    if (dim_out) *dim_out = dim;

    const int sander = 2 * dim;
    const int n_cap = std::max(GA_FOPTICS_MIN_POINTS,
                               std::min(GA_FOPTICS_MAX_MINPTS,
                                        std::max(GA_FOPTICS_MIN_POINTS, nChrom / 3)));

    int minPts;
    if (nChrom < 2 * GA_FOPTICS_ANKERST_LO) {
        // Small snapshot: cannot honour Ankerst 10–20; Ester floor when N allows.
        minPts = clampi(sander, GA_FOPTICS_MIN_POINTS, n_cap);
        if (nChrom < 2 * GA_FOPTICS_MIN_POINTS)
            minPts = clampi(2, 2, std::max(2, nChrom - 1));
    } else {
        minPts = clampi(std::max(sander, GA_FOPTICS_ANKERST_LO),
                        GA_FOPTICS_MIN_POINTS, n_cap);
        minPts = clampi(minPts, GA_FOPTICS_MIN_POINTS, GA_FOPTICS_ANKERST_HI);
        // Very large N + high dim: mild climb (Ankerst: larger MinPts reduces single-link).
        if (nChrom >= 500 && dim >= 12)
            minPts = clampi(std::max(minPts, sander), GA_FOPTICS_ANKERST_LO,
                            std::min(n_cap, GA_FOPTICS_MAX_MINPTS));
    }

    if (diversity_ratio > 0.0 && diversity_ratio < 0.05 && minPts > GA_FOPTICS_MIN_POINTS) {
        minPts = std::max(GA_FOPTICS_MIN_POINTS,
                          minPts - std::max(1, (minPts - GA_FOPTICS_MIN_POINTS) / 2));
    }

    minPts = clampi(minPts, 2, std::max(2, nChrom - 1));

    printf("[FO-MINPTS] literature=Ankerst1999[10-20]+Sander1998(2*dim)+Ester1996(floor4) "
           "nChrom=%d dim_eff=%d fdih=%d npar=%d diversity=%.4f minPts=%d "
           "(single FO pass; Ankerst band [%d,%d])\n",
           nChrom, dim, fo_ligand_fdih(FA),
           (FA && FA->npar > 0) ? FA->npar : (GB ? GB->num_genes : 0),
           diversity_ratio, minPts,
           GA_FOPTICS_ANKERST_LO, GA_FOPTICS_ANKERST_HI);

    return minPts;
}

} // namespace

void FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
    double diversity_ratio = 0.0;
    if (nChrom > 0) {
        std::set<double> distinct_cf;
        for (int i = 0; i < nChrom; ++i)
            distinct_cf.insert(chrom[i].evalue);
        diversity_ratio = static_cast<double>(distinct_cf.size()) /
                          static_cast<double>(nChrom);
    }

    // Optional super-cluster pre-filter (energy 1-D; not a second full FO pose cluster).
    if (FA->use_super_cluster && nChrom > 4) {
        std::vector<fast_optics::Point> energy_pts(static_cast<size_t>(nChrom));
        for (int i = 0; i < nChrom; ++i)
            energy_pts[static_cast<size_t>(i)].coords = { chrom[i].evalue };

        fast_optics::FastOPTICS sc_optics(
            energy_pts, std::max(GA_FOPTICS_MIN_POINTS, nChrom / GA_FOPTICS_DIVISOR));
        auto sc_indices = sc_optics.extractSuperCluster(fast_optics::ClusterMode::SUPER_CLUSTER_ONLY);

        if (!sc_indices.empty() && sc_indices.size() < static_cast<size_t>(nChrom)) {
            std::vector<bool> in_sc(static_cast<size_t>(nChrom), false);
            for (size_t idx : sc_indices)
                in_sc[idx] = true;

            int write_pos = 0;
            for (int i = 0; i < nChrom; ++i) {
                if (in_sc[static_cast<size_t>(i)]) {
                    if (i != write_pos)
                        std::swap(chrom[write_pos], chrom[i]);
                    ++write_pos;
                }
            }

            printf("--- SuperCluster pre-filter: %zu / %d poses in dominant basin ---\n",
                   sc_indices.size(), nChrom);
            nChrom = static_cast<int>(sc_indices.size());

            std::set<double> distinct_cf;
            for (int i = 0; i < nChrom; ++i)
                distinct_cf.insert(chrom[i].evalue);
            diversity_ratio = nChrom > 0
                ? static_cast<double>(distinct_cf.size()) / static_cast<double>(nChrom)
                : 0.0;
        }
    }

    // Minibatch pre-filter for very large ensembles (not a second FO clustering).
    {
        const int MINIBATCH_THRESHOLD = 10000;
        const int MINIBATCH_TARGET    = 5000;

        if (nChrom > MINIBATCH_THRESHOLD) {
            const int nAtoms_mb = residue[atoms[FA->map_par[0].atm].ofres].latm[0]
                                - residue[atoms[FA->map_par[0].atm].ofres].fatm[0] + 1;
            const int stride_mb = nAtoms_mb * 3;

            minibatch::CoordCache coord_cache;
            coord_cache.n_chrom = nChrom;
            coord_cache.stride  = stride_mb;
            coord_cache.data.resize(static_cast<std::size_t>(nChrom) * static_cast<std::size_t>(stride_mb));

            for (int c = 0; c < nChrom; ++c) {
                if (c + 1 < nChrom) {
                    calc_rmsd_chrom(FA, GB, chrom, gene_lim, atoms, residue, cleftgrid,
                                    GB->num_genes, c, c + 1,
                                    &coord_cache.data[static_cast<std::size_t>(c) * stride_mb],
                                    &coord_cache.data[static_cast<std::size_t>(c + 1) * stride_mb], false);
                } else {
                    calc_rmsd_chrom(FA, GB, chrom, gene_lim, atoms, residue, cleftgrid,
                                    GB->num_genes, c, c,
                                    &coord_cache.data[static_cast<std::size_t>(c) * stride_mb], NULL, false);
                }
            }

            std::vector<double> energies(static_cast<std::size_t>(nChrom));
            for (int i = 0; i < nChrom; ++i)
                energies[static_cast<size_t>(i)] = chrom[i].app_evalue;

            auto sample = minibatch::MinibatchSampler::farthest_point_sample(
                coord_cache, energies.data(), nAtoms_mb, MINIBATCH_TARGET, /*verbose=*/true);

            if (sample.n_selected > 0 && sample.n_selected < nChrom) {
                const auto& sel = sample.selected_indices;
                std::vector<bool> is_selected(static_cast<std::size_t>(nChrom), false);
                for (int idx : sel)
                    is_selected[static_cast<size_t>(idx)] = true;

                int write_pos = 0;
                for (int i = 0; i < nChrom; ++i) {
                    if (is_selected[static_cast<size_t>(i)]) {
                        if (i != write_pos)
                            std::swap(chrom[write_pos], chrom[i]);
                        ++write_pos;
                    }
                }
                nChrom = sample.n_selected;

                std::set<double> distinct_cf;
                for (int i = 0; i < nChrom; ++i)
                    distinct_cf.insert(chrom[i].evalue);
                diversity_ratio = nChrom > 0
                    ? static_cast<double>(distinct_cf.size()) / static_cast<double>(nChrom)
                    : 0.0;
            }
        }
    }

    if (nChrom < 2) {
        fprintf(stderr, "[FO-MINPTS] WARN: nChrom=%d too small for FO clustering\n", nChrom);
        printf("-- end of FastOPTICS_cluster --\n");
        return;
    }

    int dim = 0;
    const int minPts = fo_choose_minpts(FA, GB, nChrom, diversity_ratio, &dim);
    (void)dim;

    // Single FastOPTICS pass + single BindingPopulation emission.
    // Output: <prefix>_<minPts>_<rank>.pdb (BindingMode.cpp).
    BindingPopulation population(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom);
    FastOPTICS algo(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom,
                    population, minPts);
    algo.Execute_FastOPTICS(end_strfile, tmp_end_strfile);

    std::cout << "Size of Population is " << population.get_Population_size()
              << " Binding Modes (minPts=" << minPts << ")." << std::endl;

    population.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp,
                                 algo.get_minPoints());
    printf("-- end of FastOPTICS_cluster --\n");
}
