// test_cmaes_search.cpp — GoogleTest unit tests for CMA-ES mock path (chunk 6)
// Apache-2.0 © 2026 Le Bonhomme Pharma
//
// Links only:
//   LIB/cmaes_search.cpp  (+ include path to LIB/)
//   tests/cmaes_mock_seams_stub.cpp  (five engine seams; no flexaid_core)
//   GTest
//
// Does NOT call cmaes_run_dock or require a live FA/GB context.

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <vector>

#include "cmaes_search.h"

namespace {

// Shared mock config matching VALIDATION.md P5 / chunk1 README smoke.
// max_evals >= 5000 (chunk6 acceptance); seed 12345; dim 8.
CmaesConfig MakeMockConfig(std::int64_t max_evals = 5000) {
    CmaesConfig cfg;
    cfg.seed = 12345;
    cfg.population = 20;  // small λ for unit-test wall time
    cfg.max_evals = max_evals;
    cfg.sigma0 = 1.0;
    cfg.enable_entropy_trace = false;
    cfg.write_trace.clear();
    cfg.archive_size = 32;
    return cfg;
}

}  // namespace

// ── MockSeed12345Converges ───────────────────────────────────────────────────
// Seed 12345, dim=8, max_evals>=5000 → best_cf < 1e-2 (prefer < 1e-4).
// Sandbox P5 / chunk1 DONE: best_cf typically ~1e-8 … 1e-20 on this objective.
TEST(CmaesSearch, MockSeed12345Converges) {
    constexpr int kDim = 8;
    CmaesConfig cfg = MakeMockConfig(/*max_evals=*/5000);

    CmaesResult res;
    const int rc = cmaes_run_mock(kDim, cfg, &res, /*optional_trace=*/nullptr);

    ASSERT_EQ(rc, 0);
    ASSERT_EQ(res.status, 0);
    EXPECT_GE(res.n_evals, 5000);
    EXPECT_EQ(static_cast<int>(res.best_genes.size()), kDim);
    ASSERT_TRUE(std::isfinite(res.best_cf));
    EXPECT_GE(res.best_cf, 0.0);  // mock well is non-negative, min at origin
    EXPECT_LT(res.best_cf, 1e-2);
    // Prefer tighter tolerance when the budget is healthy (chunk1 smoke ~1e-8+).
    EXPECT_LT(res.best_cf, 1e-4);
}

// ── EntropyTraceNonEmpty ─────────────────────────────────────────────────────
// enable_entropy_trace → samples non-empty; H_energy and F finite.
TEST(CmaesSearch, EntropyTraceNonEmpty) {
    constexpr int kDim = 8;
    CmaesConfig cfg = MakeMockConfig(/*max_evals=*/5000);
    cfg.enable_entropy_trace = true;

    CmaesResult res;
    std::vector<EntropyTraceSample> trace;
    const int rc = cmaes_run_mock(kDim, cfg, &res, &trace);

    ASSERT_EQ(rc, 0);
    ASSERT_FALSE(trace.empty()) << "enable_entropy_trace must produce ≥1 sample";

    for (const auto& s : trace) {
        EXPECT_TRUE(std::isfinite(s.H_energy)) << "gen=" << s.gen;
        EXPECT_TRUE(std::isfinite(s.F)) << "gen=" << s.gen;
        EXPECT_TRUE(std::isfinite(s.H_search)) << "gen=" << s.gen;
        EXPECT_TRUE(std::isfinite(s.best_cf)) << "gen=" << s.gen;
        EXPECT_GT(s.n_evals, 0);
    }

    // Final sample should track the run best CF.
    EXPECT_NEAR(trace.back().best_cf, res.best_cf, 1e-12);
    EXPECT_EQ(trace.back().n_evals, res.n_evals);
}

// ── SnapshotDims ─────────────────────────────────────────────────────────────
// cmaes_fill_chromosomes fills num_genes correctly for each chromosome slot.
TEST(CmaesSearch, SnapshotDims) {
    constexpr int kDim = 8;
    constexpr int kMaxChrom = 16;

    CmaesConfig cfg = MakeMockConfig(/*max_evals=*/5000);
    CmaesResult res;
    ASSERT_EQ(cmaes_run_mock(kDim, cfg, &res, nullptr), 0);
    ASSERT_FALSE(res.best_genes.empty());
    ASSERT_EQ(static_cast<int>(res.best_genes.size()), kDim);

    std::vector<chromosome> chroms(static_cast<std::size_t>(kMaxChrom));
    std::vector<gene> gene_storage(
        static_cast<std::size_t>(kMaxChrom) * static_cast<std::size_t>(kDim));
    // Poison gene storage so we can detect unwritten slots.
    for (auto& g : gene_storage) {
        g.to_ic = 9999.0;
        g.to_int32 = -1;
    }

    const int filled =
        cmaes_fill_chromosomes(res, kDim, chroms.data(), kMaxChrom, gene_storage.data());

    ASSERT_GT(filled, 0);
    EXPECT_LE(filled, kMaxChrom);
    // Archive was populated by the run (default archive_size=32).
    EXPECT_LE(filled, static_cast<int>(res.archive_genes.size()) > 0
                          ? static_cast<int>(res.archive_genes.size())
                          : 1);

    for (int k = 0; k < filled; ++k) {
        EXPECT_EQ(chroms[static_cast<std::size_t>(k)].status, 'n');
        ASSERT_NE(chroms[static_cast<std::size_t>(k)].genes, nullptr);
        // Contiguous layout: gene_storage[k * num_genes + g]
        EXPECT_EQ(chroms[static_cast<std::size_t>(k)].genes,
                  gene_storage.data() +
                      static_cast<std::size_t>(k) * static_cast<std::size_t>(kDim));

        const std::vector<double>* expected_genes = nullptr;
        double expected_cf = 0.0;
        if (!res.archive_genes.empty()) {
            expected_genes = &res.archive_genes[static_cast<std::size_t>(k)];
            expected_cf = res.archive_cfs[static_cast<std::size_t>(k)];
        } else {
            expected_genes = &res.best_genes;
            expected_cf = res.best_cf;
        }
        ASSERT_EQ(static_cast<int>(expected_genes->size()), kDim);

        EXPECT_DOUBLE_EQ(chroms[static_cast<std::size_t>(k)].evalue, expected_cf);

        for (int g = 0; g < kDim; ++g) {
            EXPECT_DOUBLE_EQ(chroms[static_cast<std::size_t>(k)].genes[g].to_ic,
                             (*expected_genes)[static_cast<std::size_t>(g)])
                << "chrom=" << k << " gene=" << g;
            EXPECT_EQ(chroms[static_cast<std::size_t>(k)].genes[g].to_int32, 0);
        }
    }
}

// ── Extra: invalid args on mock API ──────────────────────────────────────────
TEST(CmaesSearch, MockRejectsBadDim) {
    CmaesConfig cfg = MakeMockConfig(100);
    CmaesResult res;
    EXPECT_EQ(cmaes_run_mock(0, cfg, &res, nullptr), -1);
    EXPECT_EQ(res.status, -1);
}

// ── H_search varies as the population collapses on a smooth well ─────────────
// Rank-only log-weights are generation-invariant (bug fixed). Softmax + allele
// histograms must produce a non-constant H_search series on this objective.
TEST(CmaesSearch, HSearchVariesWithCollapse) {
    constexpr int kDim = 8;
    CmaesConfig cfg = MakeMockConfig(/*max_evals=*/8000);
    cfg.population = 30;
    cfg.enable_entropy_trace = true;

    CmaesResult res;
    std::vector<EntropyTraceSample> trace;
    ASSERT_EQ(cmaes_run_mock(kDim, cfg, &res, &trace), 0);
    ASSERT_GE(trace.size(), 5u);

    double h_min = trace.front().H_search;
    double h_max = trace.front().H_search;
    for (const auto& s : trace) {
        ASSERT_TRUE(std::isfinite(s.H_search));
        h_min = std::min(h_min, s.H_search);
        h_max = std::max(h_max, s.H_search);
    }
    // Must not be a flat constant (the pre-fix pathology).
    EXPECT_GT(h_max - h_min, 1e-6) << "H_search range=[" << h_min << "," << h_max << "]";
}
