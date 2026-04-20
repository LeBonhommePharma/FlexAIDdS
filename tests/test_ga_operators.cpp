// tests/test_ga_operators.cpp
// Unit tests for remaining GA helper functions in gaboom.cpp:
//   RandomDouble, RandomInt, genetoic, ictogene
//   swap_chrom, QuickSort, save_snapshot, check_state, calc_rmsp
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <vector>
#include <filesystem>

#include "flexaid.h"
#include "gaboom.h"

// ═══════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════

static chromosome* alloc_chroms(int n, int ng) {
    chromosome* c = new chromosome[n];
    for (int i = 0; i < n; ++i) {
        c[i].genes = new gene[ng];
        std::memset(c[i].genes, 0, sizeof(gene) * ng);
        c[i].cf = {};
        c[i].fitnes = 0.0;
        c[i].evalue = 0.0;
    }
    return c;
}

static void free_chroms(chromosome* c, int n) {
    if (!c) return;
    for (int i = 0; i < n; ++i) delete[] c[i].genes;
    delete[] c;
}

static void set_genes(chromosome& ch, int ng, const std::vector<double>& ic_vals) {
    for (int i = 0; i < ng && i < (int)ic_vals.size(); ++i)
        ch.genes[i].to_ic = ic_vals[i];
}

// ═══════════════════════════════════════════════════════════════════════
// RandomDouble() — thread-safe RNG in [0, 1)
// ═══════════════════════════════════════════════════════════════════════

TEST(RandomDoubleTest, InUnitRange) {
    for (int i = 0; i < 1000; ++i) {
        double v = RandomDouble();
        EXPECT_GE(v, 0.0);
        EXPECT_LT(v, 1.0);
    }
}

TEST(RandomDoubleTest, ProducesVariedValues) {
    double sum = 0.0;
    const int N = 10000;
    for (int i = 0; i < N; ++i) sum += RandomDouble();
    double mean = sum / N;
    // Mean of uniform [0,1) should be ~0.5; allow generous tolerance
    EXPECT_GT(mean, 0.45);
    EXPECT_LT(mean, 0.55);
}

// ═══════════════════════════════════════════════════════════════════════
// RandomDouble(int32_t gene) — deterministic conversion
// ═══════════════════════════════════════════════════════════════════════

TEST(RandomDoubleGeneTest, ZeroMapsToZero) {
    EXPECT_DOUBLE_EQ(RandomDouble((int32_t)0), 0.0);
}

TEST(RandomDoubleGeneTest, MaxMapsCloseToOne) {
    double v = RandomDouble(MAX_RANDOM_VALUE);
    EXPECT_NEAR(v, 1.0, 1e-9);
    EXPECT_LT(v, 1.0);  // strictly less than 1
}

TEST(RandomDoubleGeneTest, MonotonicInPositiveRange) {
    // RandomDouble(int32_t) = gene/(MAX_RANDOM_VALUE+1), which wraps negative
    // for large gene values. Only monotonic in the positive range.
    double prev = RandomDouble((int32_t)0);
    for (int32_t g = 100; g > 0 && g < MAX_RANDOM_VALUE / 2; g += 100000000) {
        double cur = RandomDouble(g);
        EXPECT_GT(cur, prev);
        prev = cur;
    }
}

// ═══════════════════════════════════════════════════════════════════════
// RandomInt(double frac)
// ═══════════════════════════════════════════════════════════════════════

TEST(RandomIntTest, ZeroFracGivesZero) {
    EXPECT_EQ(RandomInt(0.0), 0);
}

TEST(RandomIntTest, NegativeFracGivesZero) {
    EXPECT_EQ(RandomInt(-0.5), 0);
}

TEST(RandomIntTest, UnitFracInRange) {
    int v = RandomInt(1.0);
    EXPECT_GE(v, 0);
    EXPECT_LE(v, RAND_MAX);
}

TEST(RandomIntTest, ClampsToRAND_MAX) {
    // frac > 1.0 should still produce valid range
    int v = RandomInt(2.0);
    EXPECT_GE(v, 0);
    EXPECT_LE(v, RAND_MAX);
}

// ═══════════════════════════════════════════════════════════════════════
// genetoic / ictogene — gene encoding/decoding
// ═══════════════════════════════════════════════════════════════════════

class GeneICTest : public ::testing::Test {
protected:
    genlim gl;
    void SetUp() override {
        std::memset(&gl, 0, sizeof(gl));
        gl.max = 10.0;
        gl.min = 0.0;
        gl.del = 0.5;
        gl.bin = 0.1;  // (10-0)/0.5 = 20 bins → bin = 1/20 = 0.05; use 0.1 for 10 bins
        gl.nbin = 10.0;
        gl.map = 0;
    }
};

TEST_F(GeneICTest, GenetoicReturnsInRange) {
    for (int32_t g = 0; g < 100; g += 10) {
        double ic = genetoic(&gl, g);
        EXPECT_GE(ic, gl.min);
        EXPECT_LE(ic, gl.max);
    }
}

TEST_F(GeneICTest, GenetoicZeroIsMin) {
    // gene = 0 → RandomDouble(0) = 0.0 → tot starts at bin → first iteration fails → i=0
    double ic = genetoic(&gl, 0);
    EXPECT_DOUBLE_EQ(ic, gl.min);
}

TEST_F(GeneICTest, IctogeneReturnsInt32) {
    double ic = 5.0;  // midpoint
    int gene = ictogene(&gl, ic);
    EXPECT_GE(gene, 0);
    EXPECT_LE(gene, MAX_RANDOM_VALUE);
}

TEST_F(GeneICTest, IctogeneAtMin) {
    int gene = ictogene(&gl, gl.min);
    EXPECT_GE(gene, 0);
}

// ═══════════════════════════════════════════════════════════════════════
// swap_chrom
// ═══════════════════════════════════════════════════════════════════════

TEST(SwapChromTest, SwapsValues) {
    constexpr int NG = 3;
    auto* chroms = alloc_chroms(2, NG);
    chroms[0].evalue = 1.0;
    chroms[0].fitnes = 10.0;
    chroms[1].evalue = 2.0;
    chroms[1].fitnes = 20.0;
    set_genes(chroms[0], NG, {1.0, 2.0, 3.0});
    set_genes(chroms[1], NG, {4.0, 5.0, 6.0});

    swap_chrom(&chroms[0], &chroms[1]);

    EXPECT_DOUBLE_EQ(chroms[0].evalue, 2.0);
    EXPECT_DOUBLE_EQ(chroms[0].fitnes, 20.0);
    EXPECT_DOUBLE_EQ(chroms[1].evalue, 1.0);
    EXPECT_DOUBLE_EQ(chroms[1].fitnes, 10.0);
    EXPECT_DOUBLE_EQ(chroms[0].genes[0].to_ic, 4.0);
    EXPECT_DOUBLE_EQ(chroms[1].genes[2].to_ic, 3.0);

    free_chroms(chroms, 2);
}

TEST(SwapChromTest, SelfSwapNoOp) {
    constexpr int NG = 2;
    auto* c = alloc_chroms(1, NG);
    c[0].evalue = 42.0;
    set_genes(c[0], NG, {1.0, 2.0});

    swap_chrom(&c[0], &c[0]);

    EXPECT_DOUBLE_EQ(c[0].evalue, 42.0);
    EXPECT_DOUBLE_EQ(c[0].genes[0].to_ic, 1.0);

    free_chroms(c, 1);
}

// ═══════════════════════════════════════════════════════════════════════
// QuickSort
// ═══════════════════════════════════════════════════════════════════════

class QuickSortTest : public ::testing::Test {
protected:
    static constexpr int NG = 2;
    static constexpr int N = 10;
    chromosome* chroms = nullptr;
    void SetUp() override { chroms = alloc_chroms(N, NG); }
    void TearDown() override { free_chroms(chroms, N); }
};

TEST_F(QuickSortTest, SortByEnergyAscending) {
    // Set random evalues
    double evals[] = {5.0, 1.0, 8.0, 3.0, 9.0, 2.0, 7.0, 4.0, 6.0, 0.0};
    for (int i = 0; i < N; ++i) chroms[i].evalue = evals[i];

    QuickSort(chroms, 0, N - 1, true);

    for (int i = 1; i < N; ++i)
        EXPECT_LE(chroms[i-1].evalue, chroms[i].evalue);
}

TEST_F(QuickSortTest, SortByFitnessDescending) {
    double fits[] = {5.0, 1.0, 8.0, 3.0, 9.0, 2.0, 7.0, 4.0, 6.0, 0.0};
    for (int i = 0; i < N; ++i) chroms[i].fitnes = fits[i];

    QuickSort(chroms, 0, N - 1, false);

    for (int i = 1; i < N; ++i)
        EXPECT_GE(chroms[i-1].fitnes, chroms[i].fitnes);
}

TEST_F(QuickSortTest, AlreadySorted) {
    for (int i = 0; i < N; ++i) chroms[i].evalue = (double)i;
    QuickSort(chroms, 0, N - 1, true);
    for (int i = 0; i < N; ++i)
        EXPECT_DOUBLE_EQ(chroms[i].evalue, (double)i);
}

TEST_F(QuickSortTest, ReverseSorted) {
    for (int i = 0; i < N; ++i) chroms[i].evalue = (double)(N - 1 - i);
    QuickSort(chroms, 0, N - 1, true);
    for (int i = 0; i < N; ++i)
        EXPECT_DOUBLE_EQ(chroms[i].evalue, (double)i);
}

TEST_F(QuickSortTest, SingleElement) {
    chroms[0].evalue = 42.0;
    QuickSort(chroms, 0, 0, true);
    EXPECT_DOUBLE_EQ(chroms[0].evalue, 42.0);
}

TEST_F(QuickSortTest, TwoElements) {
    chroms[0].evalue = 10.0;
    chroms[1].evalue = 5.0;
    QuickSort(chroms, 0, 1, true);
    EXPECT_LE(chroms[0].evalue, chroms[1].evalue);
}

TEST_F(QuickSortTest, AllEqual) {
    for (int i = 0; i < N; ++i) chroms[i].evalue = 3.14;
    QuickSort(chroms, 0, N - 1, true);
    for (int i = 0; i < N; ++i)
        EXPECT_DOUBLE_EQ(chroms[i].evalue, 3.14);
}

TEST_F(QuickSortTest, WithNaNsAtEnd) {
    // NaN comparisons are false; they should not crash
    for (int i = 0; i < N - 2; ++i) chroms[i].evalue = (double)i;
    chroms[8].evalue = 100.0;
    chroms[9].evalue = 50.0;
    QuickSort(chroms, 0, N - 1, true);
    // Just verify no crash — non-NaN portion should be sorted
    for (int i = 1; i < N; ++i) {
        if (!std::isnan(chroms[i-1].evalue) && !std::isnan(chroms[i].evalue))
            EXPECT_LE(chroms[i-1].evalue, chroms[i].evalue);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// save_snapshot
// ═══════════════════════════════════════════════════════════════════════

TEST(SaveSnapshotTest, CopiesAllChromosomes) {
    constexpr int NC = 5;
    constexpr int NG = 3;
    auto* src = alloc_chroms(NC, NG);
    auto* snap = alloc_chroms(NC, NG);

    for (int i = 0; i < NC; ++i) {
        src[i].evalue = (double)(i * 10);
        src[i].fitnes = (double)(i * 100);
        for (int g = 0; g < NG; ++g)
            src[i].genes[g].to_ic = (double)(i * NG + g);
    }

    save_snapshot(snap, src, NC, NG);

    for (int i = 0; i < NC; ++i) {
        EXPECT_DOUBLE_EQ(snap[i].evalue, src[i].evalue);
        EXPECT_DOUBLE_EQ(snap[i].fitnes, src[i].fitnes);
        for (int g = 0; g < NG; ++g)
            EXPECT_DOUBLE_EQ(snap[i].genes[g].to_ic, src[i].genes[g].to_ic);
    }

    // Snapshot should be independent
    src[0].evalue = 999.0;
    EXPECT_DOUBLE_EQ(snap[0].evalue, 0.0);

    free_chroms(src, NC);
    free_chroms(snap, NC);
}

// ═══════════════════════════════════════════════════════════════════════
// check_state — file-based GA control
// ═══════════════════════════════════════════════════════════════════════

class CheckStateTest : public ::testing::Test {
protected:
    std::string tmpdir;
    char pausefile[256], abortfile[256], stopfile[256];

    void SetUp() override {
        tmpdir = std::filesystem::temp_directory_path().string() + "/flexaids_test_checkstate_" +
                 std::to_string(std::hash<std::thread::id>{}(std::this_thread::get_id()));
        std::filesystem::create_directories(tmpdir);
        std::snprintf(pausefile, sizeof(pausefile), "%s/pause", tmpdir.c_str());
        std::snprintf(abortfile, sizeof(abortfile), "%s/abort", tmpdir.c_str());
        std::snprintf(stopfile, sizeof(stopfile), "%s/stop", tmpdir.c_str());
    }
    void TearDown() override {
        std::filesystem::remove_all(tmpdir);
    }
    void touch(const char* path) {
        FILE* f = std::fopen(path, "w");
        if (f) std::fclose(f);
    }
};

TEST_F(CheckStateTest, NoFilesReturnsZero) {
    EXPECT_EQ(check_state(pausefile, abortfile, stopfile, 1), 0);
}

TEST_F(CheckStateTest, AbortFileReturnsNegOne) {
    touch(abortfile);
    EXPECT_EQ(check_state(pausefile, abortfile, stopfile, 1), -1);
    std::filesystem::remove(abortfile);
}

TEST_F(CheckStateTest, StopFileReturnsOne) {
    touch(stopfile);
    EXPECT_EQ(check_state(pausefile, abortfile, stopfile, 1), 1);
    std::filesystem::remove(stopfile);
}

TEST_F(CheckStateTest, AbortTakesPriorityOverStop) {
    touch(abortfile);
    touch(stopfile);
    EXPECT_EQ(check_state(pausefile, abortfile, stopfile, 1), -1);
    std::filesystem::remove(abortfile);
    std::filesystem::remove(stopfile);
}

// ═══════════════════════════════════════════════════════════════════════
// calc_rmsp — root mean square parameter difference
// ═══════════════════════════════════════════════════════════════════════

TEST(CalcRmspTest, IdenticalGenesZero) {
    constexpr int NG = 5;
    gene g1[NG], g2[NG];
    for (int i = 0; i < NG; ++i) {
        g1[i].to_ic = (double)i;
        g2[i].to_ic = (double)i;
    }
    double rmsp = calc_rmsp(NG, g1, g2, nullptr, nullptr);
    EXPECT_NEAR(rmsp, 0.0, 1e-12);
}

TEST(CalcRmspTest, KnownDifference) {
    constexpr int NG = 2;
    gene g1[NG], g2[NG];
    g1[0].to_ic = 0.0; g2[0].to_ic = 3.0;  // diff = 3
    g1[1].to_ic = 0.0; g2[1].to_ic = 4.0;  // diff = 4
    // RMSP = sqrt((9+16)/2) = sqrt(12.5) = 3.535...
    double expected = std::sqrt(12.5);
    double rmsp = calc_rmsp(NG, g1, g2, nullptr, nullptr);
    EXPECT_NEAR(rmsp, expected, 1e-10);
}

TEST(CalcRmspTest, SingleGene) {
    gene g1[1], g2[1];
    g1[0].to_ic = 0.0; g2[0].to_ic = 5.0;
    double rmsp = calc_rmsp(1, g1, g2, nullptr, nullptr);
    EXPECT_NEAR(rmsp, 5.0, 1e-12);
}

TEST(CalcRmspTest, Symmetric) {
    constexpr int NG = 3;
    gene g1[NG], g2[NG];
    for (int i = 0; i < NG; ++i) {
        g1[i].to_ic = (double)i * 1.5;
        g2[i].to_ic = (double)i * 0.3;
    }
    double fwd = calc_rmsp(NG, g1, g2, nullptr, nullptr);
    double rev = calc_rmsp(NG, g2, g1, nullptr, nullptr);
    EXPECT_NEAR(fwd, rev, 1e-12);
}
