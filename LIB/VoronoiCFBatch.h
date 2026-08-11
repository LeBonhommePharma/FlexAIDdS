// VoronoiCFBatch.h — C++20 std::span-based batch CF evaluation interface
//
// Wraps the existing ic2cf / vcfunction / Vcontacts pipeline with zero-copy
// std::span views for chromosome gene arrays and atom buffers.  Designed for
// the GA's inner evaluation loop where we need to score a large population of
// chromosomes without redundant copies.
//
// Key features
//   • std::span<const gene>   — read-only view of one chromosome's genes
//   • std::span<atom>         — per-thread mutable atom view (no extra alloc)
//   • batch_eval()            — parallel evaluation of chromosome population
//   • benchmark()             — wall-clock comparison: serial vs OpenMP
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#pragma once

#include <span>
#include <vector>
#include <functional>
#include <chrono>
#include <cstring>
#include <cassert>
#include <cmath>

#include "flexaid.h"   // cfstr, FA_Global, VC_Global, atom, resid, gene, genlim
#include "gaboom.h"    // chromosome, eval_chromosome, get_apparent_cf_evalue
#include "Vcontacts.h" // atomsas, ca_struct, contactlist, ptindex, vertex, plane, edgevector

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace voronoi_cf {

// ─── Span-based eval_chromosome overload ─────────────────────────────────────
//
// Identical semantics to the C-style eval_chromosome() in gaboom.cpp, but
// accepts std::span for gene and genlim parameters, removing the reliance on
// implicit size from GB->num_genes.
//
// The caller owns all pointed-to data; spans are non-owning views.

inline cfstr eval_span(
    FA_Global*             FA,
    GB_Global*             GB,
    VC_Global*             VC,
    std::span<const genlim> gene_lim,
    std::span<atom>         atoms,
    std::span<resid>        residue,
    gridpoint*             cleftgrid,
    std::span<const gene>   john,
    cfstr (*function)(FA_Global*, VC_Global*, atom*, resid*, gridpoint*, int, double*))
{
    assert(john.size() == gene_lim.size());
    assert(john.size() <= MAX_NUM_GENES);

    double icv[MAX_NUM_GENES] = {};
    for (std::size_t i = 0; i < john.size(); ++i) {
        double val = john[i].to_ic;
        if (val > gene_lim[i].max) val = gene_lim[i].max;
        if (val < gene_lim[i].min) val = gene_lim[i].min;
        icv[i] = val;
    }
    return (*function)(FA, VC,
                       atoms.data(), residue.data(),
                       cleftgrid,
                       static_cast<int>(john.size()),
                       icv);
}

// ─── Per-thread workspace ─────────────────────────────────────────────────────
//
// Holds private copies of every mutable buffer that ic2cf/vcfunction/Vcontacts
// write to, allowing lock-free parallel evaluation.

struct ThreadWorkspace {
    FA_Global         fa;           // shallow copy; pointer fields redirected
    VC_Global         vc;           // shallow copy; pointer fields redirected

    std::vector<atom>        atoms;
    std::vector<resid>       residue;
    std::vector<int>         contacts;
    std::vector<float>       contributions;
    std::vector<OptRes>      optres;
    std::vector<atomsas>     calc;
    std::vector<int>         calclist;
    std::vector<int>         ca_index;
    std::vector<ca_struct>   ca_rec;
    std::vector<int>         seed;
    std::vector<contactlist> contlist;
    std::vector<ptindex>     ptorder;
    std::vector<vertex>      centerpt;
    std::vector<vertex>      poly;
    std::vector<plane>       cont;
    std::vector<edgevector>  vedge;

    // Construct and wire all buffers for the given FA/VC reference state.
    ThreadWorkspace(const FA_Global* FA, const VC_Global* VC,
                    const atom* ref_atoms, const resid* ref_residue)
    {
        const int natm  = FA->atm_cnt;
        const int natmr = FA->atm_cnt_real;
        const int nres  = FA->res_cnt;
        const int nopt  = FA->num_optres;
        const int nctb  = FA->ntypes * FA->ntypes;

        atoms        .assign(ref_atoms,   ref_atoms   + natm);
        residue      .assign(ref_residue, ref_residue + nres);
        // CONTACTS_BUFFER_SIZE: stamps + the trailing epoch slot (flexaid.h).
        // Allocated together with `fa` below so the epoch and the stamps it
        // guards share this workspace's lifetime.
        contacts     .assign(CONTACTS_BUFFER_SIZE, 0);
        contributions.assign(static_cast<std::size_t>(nctb), 0.0f);
        optres       .assign(FA->optres,  FA->optres  + nopt);
        calc         .resize(static_cast<std::size_t>(natmr));
        calclist     .resize(static_cast<std::size_t>(natmr));
        ca_index     .assign(static_cast<std::size_t>(natmr), -1);
        ca_rec       .resize(static_cast<std::size_t>(VC->ca_recsize));
        seed         .resize(static_cast<std::size_t>(3 * natmr));
        contlist     .resize(10000);
        ptorder      .resize(MAX_PT);
        centerpt     .resize(MAX_PT);
        poly         .resize(MAX_POLY);
        cont         .resize(MAX_PT);
        vedge        .resize(MAX_POLY);

        fa = *FA;
        fa.contacts      = contacts     .data();
        fa.contributions = contributions.data();
        fa.optres        = optres       .data();

        vc = *VC;
        vc.Calc      = calc    .data();
        vc.Calclist  = calclist.data();
        vc.ca_index  = ca_index.data();
        vc.ca_rec    = ca_rec  .data();
        vc.seed      = seed    .data();
        vc.contlist  = contlist.data();
        vc.ptorder   = ptorder .data();
        vc.centerpt  = centerpt.data();
        vc.poly      = poly    .data();
        vc.cont      = cont    .data();
        vc.vedge     = vedge   .data();
    }

    // Reset per-chromosome mutable state without re-allocating.
    void reset(const FA_Global* FA, const atom* ref_atoms, const resid* ref_residue)
    {
        const int nopt = FA->num_optres;
        std::copy(ref_atoms,   ref_atoms   + FA->atm_cnt, atoms  .begin());
        std::copy(ref_residue, ref_residue + FA->res_cnt, residue.begin());
        for (int o = 0; o < nopt; ++o) {
            optres[o].cf.com         = 0.0;
            optres[o].cf.wal         = 0.0;
            optres[o].cf.sas         = 0.0;
            optres[o].cf.totsas      = 0.0;
            optres[o].cf.con         = 0.0;
            optres[o].cf.gist        = 0.0;
            optres[o].cf.hbond       = 0.0;
            optres[o].cf.elec        = 0.0;
            optres[o].cf.metal_coord = 0.0;
            optres[o].cf.gist_desolv = 0.0;
            optres[o].cf.rclash      = 0;
        }
        vc.numcarec = 0;
    }
};

// ─── Batch result ─────────────────────────────────────────────────────────────

struct BatchResult {
    std::vector<cfstr>  cf;          // one cfstr per chromosome
    std::vector<double> app_evalue;  // apparent CF evalue for quick ranking
    double              wall_ms = 0.0; // wall-clock evaluation time (ms)
};

// ─── batch_eval ──────────────────────────────────────────────────────────────
//
// Evaluate a population of chromosomes using OpenMP + std::span zero-copy
// views.  Results are written into BatchResult (does NOT modify chrom[]).
//
// chroms:    read-only view of the chromosome array (status != 'n' are skipped)
// gene_lim:  gene boundary array (std::span, bounds-safe)
// function:  the CF scoring function (ic2cf or cffunction)

inline BatchResult batch_eval(
    std::span<const chromosome> chroms,
    FA_Global*                  FA,
    GB_Global*                  GB,
    VC_Global*                  VC,
    std::span<const genlim>     gene_lim,
    const atom*                 ref_atoms,
    const resid*                ref_residue,
    gridpoint*                  cleftgrid,
    cfstr (*function)(FA_Global*, VC_Global*, atom*, resid*, gridpoint*, int, double*))
{
    const int N = static_cast<int>(chroms.size());
    BatchResult result;
    result.cf       .resize(static_cast<std::size_t>(N));
    result.app_evalue.resize(static_cast<std::size_t>(N));

    const int n_thr =
#ifdef _OPENMP
        omp_get_max_threads();
#else
        1;
#endif

    // Build per-thread workspaces (heap-allocated; one per thread).
    std::vector<ThreadWorkspace> ws;
    ws.reserve(static_cast<std::size_t>(n_thr));
    for (int t = 0; t < n_thr; ++t)
        ws.emplace_back(FA, VC, ref_atoms, ref_residue);

    auto t0 = std::chrono::steady_clock::now();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) \
    shared(chroms, result, ws, gene_lim, cleftgrid, function, FA, GB, \
           ref_atoms, ref_residue, N)
#endif
    for (int i = 0; i < N; ++i) {
        if (chroms[i].status == 'n') continue;
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        ws[tid].reset(FA, ref_atoms, ref_residue);

        std::span<const gene>   gene_span(chroms[i].genes,
                                          static_cast<std::size_t>(GB->num_genes));
        std::span<atom>         atom_span(ws[tid].atoms);
        std::span<resid>        res_span (ws[tid].residue);

        result.cf[static_cast<std::size_t>(i)] = eval_span(
            &ws[tid].fa, GB, &ws[tid].vc,
            gene_lim, atom_span, res_span,
            cleftgrid, gene_span, function);

        result.app_evalue[static_cast<std::size_t>(i)] =
            get_apparent_cf_evalue(&result.cf[static_cast<std::size_t>(i)]);
    }

    auto t1 = std::chrono::steady_clock::now();
    result.wall_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    return result;
}

// ─── benchmark ────────────────────────────────────────────────────────────────
//
// Times serial vs OpenMP batch_eval over `n_reps` repetitions and prints a
// summary table to stdout.  Returns the parallel result for correctness checks.

struct BenchmarkReport {
    double serial_ms;
    double parallel_ms;
    double speedup;
    int    n_threads;
    int    n_chroms;
    BatchResult parallel_result;  // kept for correctness validation
};

inline BenchmarkReport benchmark(
    std::span<const chromosome> chroms,
    FA_Global*                  FA,
    GB_Global*                  GB,
    VC_Global*                  VC,
    std::span<const genlim>     gene_lim,
    const atom*                 ref_atoms,
    const resid*                ref_residue,
    gridpoint*                  cleftgrid,
    cfstr (*function)(FA_Global*, VC_Global*, atom*, resid*, gridpoint*, int, double*),
    int n_reps = 3)
{
    BenchmarkReport report{};
    report.n_chroms  = static_cast<int>(chroms.size());

#ifdef _OPENMP
    report.n_threads = omp_get_max_threads();
#else
    report.n_threads = 1;
#endif

    // Serial baseline (force 1 thread)
    double serial_total = 0.0;
    for (int r = 0; r < n_reps; ++r) {
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        auto res = batch_eval(chroms, FA, GB, VC, gene_lim,
                              ref_atoms, ref_residue, cleftgrid, function);
        serial_total += res.wall_ms;
    }
    report.serial_ms = serial_total / n_reps;

    // Parallel (all cores)
    double par_total = 0.0;
    for (int r = 0; r < n_reps; ++r) {
#ifdef _OPENMP
        omp_set_num_threads(report.n_threads);
#endif
        report.parallel_result = batch_eval(chroms, FA, GB, VC, gene_lim,
                                            ref_atoms, ref_residue,
                                            cleftgrid, function);
        par_total += report.parallel_result.wall_ms;
    }
    report.parallel_ms = par_total / n_reps;
    report.speedup     = report.serial_ms / report.parallel_ms;

    // Print benchmark table to stdout
    std::printf(
        "\n╔══ VoronoiCFBatch benchmark ══════════════════════════════════════╗\n"
        "║  Chromosomes : %6d                                           ║\n"
        "║  Threads     : %6d   (OpenMP)                               ║\n"
        "║  Repetitions : %6d                                           ║\n"
        "║──────────────────────────────────────────────────────────────────║\n"
        "║  Serial      : %8.2f ms                                     ║\n"
        "║  Parallel    : %8.2f ms                                     ║\n"
        "║  Speedup     : %8.2f ×                                      ║\n"
        "╚══════════════════════════════════════════════════════════════════╝\n",
        report.n_chroms, report.n_threads, n_reps,
        report.serial_ms, report.parallel_ms, report.speedup);

    return report;
}

} // namespace voronoi_cf
