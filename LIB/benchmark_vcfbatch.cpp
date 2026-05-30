// benchmark_vcfbatch.cpp — VoronoiCFBatch std::span batch-eval benchmark
//
// Exercises the C++20 std::span-based batch evaluation interface from
// VoronoiCFBatch.h.  Builds a synthetic chromosome population (random genes
// within gene_lim bounds), runs serial (1-thread) and parallel (all-core)
// batch_eval(), and prints the speedup table.
//
// Build:  cmake -B build -DENABLE_VCFBATCH_BENCHMARK=ON && cmake --build build
// Run:    ./build/benchmark_vcfbatch [pop_size] [n_genes]
//         (defaults: pop_size=200, n_genes=20)
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include "VoronoiCFBatch.h"
#include "gaboom.h"
#include "flexaid.h"
#include "Vcontacts.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <vector>
#include <span>

// ── Minimal stub CF function (no real Vcontacts geometry; just returns
//    a synthetic cfstr proportional to the sum of gene values).
//    Replace with ic2cf for real benchmarking on an actual receptor.
static cfstr stub_cf(FA_Global* FA, VC_Global* /*VC*/,
                     atom* /*atoms*/, resid* /*residue*/,
                     gridpoint* /*cleftgrid*/, int n_genes, double* icv)
{
    cfstr cf{};
    double sum = 0.0;
    for (int i = 0; i < n_genes; ++i) sum += icv[i];
    cf.com = -std::abs(sum) * 0.01;   // synthetic complementarity
    cf.wal =  0.0;
    cf.sas =  0.0;
    cf.con =  0.0;
    cf.totsas = 1.0;
    cf.rclash = 0;

    // Simulate realistic compute cost (~10 µs per chromosome) - safe bounded access
    volatile double dummy = 0.0;
    for (int k = 0; k < 200; ++k) {  // reduced from 5000 to avoid any edge issues in synthetic run
        int idx = k % n_genes;
        dummy += std::sin(icv[idx] + static_cast<double>(k));
    }
    (void)dummy;
    return cf;
}

double get_apparent_cf_evalue(cfstr* cf)
{
    return cf ? cf->com + cf->wal + cf->sas + cf->elec : 0.0;
}

int main(int argc, char* argv[])
{
    const int pop_size = (argc > 1) ? std::atoi(argv[1]) : 200;
    const int n_genes  = (argc > 2) ? std::atoi(argv[2]) : 20;

    if (n_genes < 1 || n_genes > MAX_NUM_GENES || pop_size < 1) {
        std::fprintf(stderr, "Usage: benchmark_vcfbatch [pop_size 1-2000] [n_genes 1-%d]\n",
                     MAX_NUM_GENES);
        return 1;
    }

    std::printf("VoronoiCFBatch — std::span batch-eval benchmark\n");
    std::printf("  pop_size = %d   n_genes = %d\n\n", pop_size, n_genes);

    // ── Synthetic GA structures ──────────────────────────────────────────────
    // These are minimal stubs; a real benchmark would load an actual receptor.

    FA_Global FA{};
    FA.atm_cnt      = 1;
    FA.atm_cnt_real = 1;
    FA.res_cnt      = 1;
    FA.num_optres   = 0;
    FA.ntypes       = 1;

    GB_Global GB{};
    GB.num_genes  = n_genes;
    GB.num_chrom  = pop_size;

    VC_Global VC{};
    VC.ca_recsize = 1;

    // Stub atom / residue / optres / contacts / contributions
    // Provide at least one element so ThreadWorkspace assign() calls are safe
    // even when counts are minimal (prevents potential bad pointer arithmetic in synthetic benchmark).
    static atom  stub_atoms[4]{};
    static resid stub_residues[4]{};
    static OptRes stub_optres[1]{};

    atom*  stub_atom   = &stub_atoms[0];
    resid* stub_res    = &stub_residues[0];
    gridpoint stub_gp{};

    FA.contacts      = nullptr;
    FA.contributions = nullptr;
    FA.optres        = &stub_optres[0];   // non-null even if num_optres==0
    FA.num_optres    = 0;                 // keep count at 0 for the test

    // Gene limits: each gene in [−180, 180] (dihedral-like)
    std::vector<genlim> gene_lim(static_cast<std::size_t>(n_genes));
    for (int i = 0; i < n_genes; ++i) {
        gene_lim[i].min = -180.0;
        gene_lim[i].max =  180.0;
    }

    // Build random chromosome population
    std::mt19937 rng{42};
    std::uniform_real_distribution<double> dist(-180.0, 180.0);

    std::vector<chromosome> chroms(static_cast<std::size_t>(pop_size));
    for (auto& c : chroms) {
        c.status = 's';  // 's' = needs scoring (not 'n' = already done)
        for (int j = 0; j < n_genes; ++j) {
            c.genes[j].to_ic = dist(rng);
        }
    }

    // ── std::span views ──────────────────────────────────────────────────────
    std::span<const chromosome> chrom_span(chroms);
    std::span<const genlim>     lim_span(gene_lim);

    // ── Run benchmark ────────────────────────────────────────────────────────
    auto report = voronoi_cf::benchmark(
        chrom_span, &FA, &GB, &VC, lim_span,
        stub_atom, stub_res, &stub_gp,
        stub_cf,
        /*n_reps=*/5);

    // Correctness spot-check: all parallel app_evalue should be ≤ 0
    int ok = 0;
    for (double v : report.parallel_result.app_evalue)
        if (v <= 0.0) ++ok;
    std::printf("Correctness check: %d / %d chromosomes have app_evalue ≤ 0\n",
                ok, pop_size);

    return (report.speedup >= 1.0) ? 0 : 1;
}
