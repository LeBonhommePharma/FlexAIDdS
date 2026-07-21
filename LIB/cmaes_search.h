// cmaes_search.h — FlexAIDdS CMA-ES search backend (chunk 1 adapter)
// Apache-2.0 clean-room CMA-ES (Hansen 2006 style). No GPL code.
//
// Engine seam (implemented in gaboom.cpp / ic2cf.cpp; NOT modified here):
//   set_gene_lim, set_bins, eval_chromosome, get_cf_evalue, get_apparent_cf_evalue
//
// Namespace / prefix: flexaids_cmaes (types) + cmaes_* free functions.

#ifndef FLEXAIDS_CMAES_SEARCH_H
#define FLEXAIDS_CMAES_SEARCH_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Engine types (include-guarded). Prefer including only in .cpp when possible;
// dock / snapshot APIs need complete types for cfstr-by-value and chromosome.
#include "gaboom.h"

namespace flexaids_cmaes {

// ── Configuration ────────────────────────────────────────────────────────────
struct CmaesConfig {
    int          population = 1000;     // λ (offspring per generation)
    std::int64_t max_evals  = 2000000;  // hard evaluation budget
    std::uint32_t seed      = 1;
    double       sigma0     = 0.3;      // initial step-size
    std::string  write_trace;           // empty = no file write from run
    bool         enable_entropy_trace = false;  // fill optional_trace when true
    int          archive_size = 32;     // top-K elites retained for snapshot
};

// ── Result ───────────────────────────────────────────────────────────────────
// status: 0 = completed (budget exhausted);
//        -1 = bad args; -2 = internal error (exception caught)
struct CmaesResult {
    double              best_cf     = 0.0;
    double              best_app_cf = 0.0;
    std::vector<double> best_genes;  // continuous IC values (to_ic), length = num_genes
    int                 n_evals     = 0;
    int                 n_gens      = 0;
    int                 status      = 0;

    // Ranked elite archive (best first) for cmaes_fill_chromosomes top-K
    std::vector<std::vector<double>> archive_genes;
    std::vector<double>              archive_cfs;
    std::vector<double>              archive_app_cfs;
};

// ── Entropy / free-energy proxy sample (diagnostic only) ─────────────────────
// H_search : Shannon entropy of normalized rank-μ selection weights (nats)
// H_energy : Shannon entropy of Boltzmann weights over the λ sample (nats)
// F        : free-energy proxy = best_cf - T * H_energy
//            T = kB_kcal * T_K with T_K = FA->temperature if >0 else 300 K
//            (kB_kcal*300 ≈ 0.596 kcal/mol)
struct EntropyTraceSample {
    int    gen      = 0;
    double H_search = 0.0;
    double H_energy = 0.0;
    double F        = 0.0;
    double best_cf  = 0.0;
    int    n_evals  = 0;
};

}  // namespace flexaids_cmaes

// ── Public free functions (prefix cmaes_) ────────────────────────────────────

using flexaids_cmaes::CmaesConfig;
using flexaids_cmaes::CmaesResult;
using flexaids_cmaes::EntropyTraceSample;

// Target function type matching eval_chromosome's scoring callback.
using CmaesTargetFn = cfstr (*)(FA_Global*, VC_Global*, atom*, resid*,
                                gridpoint*, int, double*);

// Run CMA-ES docking search against the live FlexAID engine seam.
// gene_lim must be pre-allocated (GB->num_genes entries). Calls set_gene_lim
// and set_bins to refresh limits/bins from FA before sampling.
// Minimizes get_cf_evalue (CF / contact-function scoring proxy).
// Returns 0 on success, negative on error. Does not throw across the boundary.
int cmaes_run_dock(
    FA_Global* FA,
    GB_Global* GB,
    VC_Global* VC,
    genlim* gene_lim,
    atom* atoms,
    resid* residue,
    gridpoint* cleftgrid,
    CmaesTargetFn target,
    const CmaesConfig& config,
    CmaesResult* result,
    std::vector<EntropyTraceSample>* optional_trace = nullptr);

// Mock multi-dimensional well: quadratic + mild nearest-neighbor coupling.
// Minimum at origin with f(0)=0. Seed 12345 + dim 5–10 + few thousand evals
// should yield best_cf on the order of 1e-8 … 1e-2.
double cmaes_mock_objective(const double* x, int n);

// Pure CMA-ES on the mock objective (no engine). Bounds default to [-5,5]^dim.
int cmaes_run_mock(
    int dim,
    const CmaesConfig& config,
    CmaesResult* result,
    std::vector<EntropyTraceSample>* optional_trace = nullptr);

// Fill top-K gene vectors from result into a chromosome array for clustering.
// gene_storage must hold at least max_chrom * num_genes genes (contiguous).
// Sets status='n', evalue/app_evalue from archive CF. Returns count filled.
int cmaes_fill_chromosomes(
    const CmaesResult& result,
    int num_genes,
    chromosome* chrom_out,
    int max_chrom,
    gene* gene_storage);

// Write entropy CSV: gen,H_search,H_energy,F,best_cf,n_evals
void cmaes_write_trace_csv(
    const std::string& path,
    const std::vector<EntropyTraceSample>& samples);

#endif  // FLEXAIDS_CMAES_SEARCH_H
