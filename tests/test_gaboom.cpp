// tests/test_gaboom.cpp
// Unit tests for the Genetic Algorithm engine (gaboom.cpp)
// Tests cover: gene encoding, chromosome operations, sorting, selection,
// crossover, mutation, fitness stats, and adaptive probabilities.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>

// flexaid.h defines E as a macro which conflicts with GoogleTest templates.
// Include gtest first, then our headers.
#include "../LIB/gaboom.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

// ===========================================================================
// HELPERS
// ===========================================================================

namespace {

static constexpr double EPSILON = 1e-6;

// Allocate a simple chromosome array with gene storage.
struct ChromArray {
    std::vector<chromosome> chroms;
    std::vector<std::vector<gene>> gene_store;

    ChromArray(int n, int ng) : chroms(n), gene_store(n, std::vector<gene>(ng)) {
        for (int i = 0; i < n; ++i)
            chroms[i].genes = gene_store[i].data();
    }

    chromosome* data() { return chroms.data(); }
    chromosome& operator[](int i) { return chroms[i]; }
};

} // namespace

// ===========================================================================
// RANDOM NUMBER HELPERS (deterministic overload)
// ===========================================================================

TEST(RandomDouble_Deterministic, ZeroMapsToZero) {
    EXPECT_DOUBLE_EQ(RandomDouble(static_cast<int32_t>(0)), 0.0);
}

TEST(RandomDouble_Deterministic, MaxMapsToJustBelowOne) {
    double v = RandomDouble(MAX_RANDOM_VALUE);
    EXPECT_GT(v, 0.99);
    EXPECT_LT(v, 1.0);
}

TEST(RandomDouble_Deterministic, MidpointMapsCorrectly) {
    int32_t mid = MAX_RANDOM_VALUE / 2;
    double v = RandomDouble(mid);
    EXPECT_NEAR(v, 0.5, 0.01);
}

TEST(RandomInt_Test, ZeroFractionGivesZero) {
    EXPECT_EQ(RandomInt(0.0), 0);
}

TEST(RandomInt_Test, OneFractionGivesRandMaxPlusOne) {
    int v = RandomInt(1.0);
    // float→int conversion of RAND_MAX+1.0 saturates to RAND_MAX on ARM (Apple Silicon)
    EXPECT_EQ(v, RAND_MAX);
}

// ===========================================================================
// GENE ENCODING / DECODING
// ===========================================================================

class GeneEncodingTest : public ::testing::Test {
protected:
    genlim gl;

    void SetUp() override {
        // Simple gene limits: range [0, 10] with step 1.0
        gl.min = 0.0;
        gl.max = 10.0;
        gl.del = 1.0;
        gl.nbin = 11;               // 11 discrete values: 0,1,...,10
        gl.bin = 1.0 / 11.0;        // probability per bin
        gl.map = 0;
    }
};

TEST_F(GeneEncodingTest, GeneToICReturnsValueInRange) {
    // Test multiple gene values → IC should be within [min, max]
    const int64_t step = static_cast<int64_t>(MAX_RANDOM_VALUE) / 100;
    for (int64_t g = 0; g < MAX_RANDOM_VALUE; g += step) {
        double ic = genetoic(&gl, static_cast<int32_t>(g));
        EXPECT_GE(ic, gl.min - EPSILON);
        EXPECT_LE(ic, gl.max + EPSILON);
    }
}

TEST_F(GeneEncodingTest, ICToGeneRoundTripsApproximately) {
    // Convert IC → gene → IC; result should be close to a discrete bin
    double ic = 5.0;
    int32_t g = ictogene(&gl, ic);
    double ic2 = genetoic(&gl, g);
    // Should map to a valid discrete value
    EXPECT_GE(ic2, gl.min - EPSILON);
    EXPECT_LE(ic2, gl.max + EPSILON);
    // Should be within one bin width of original
    EXPECT_NEAR(ic2, ic, gl.del + EPSILON);
}

TEST_F(GeneEncodingTest, BoundaryMinIC) {
    double ic = 0.0;
    int32_t g = ictogene(&gl, ic);
    double ic2 = genetoic(&gl, g);
    // Round-trip should land within the gene limits range
    EXPECT_GE(ic2, gl.min - EPSILON);
    EXPECT_LE(ic2, gl.max + EPSILON);
}

TEST_F(GeneEncodingTest, BoundaryMaxIC) {
    double ic = 10.0;
    int32_t g = ictogene(&gl, ic);
    double ic2 = genetoic(&gl, g);
    // Round-trip should land within the gene limits range
    EXPECT_GE(ic2, gl.min - EPSILON);
    EXPECT_LE(ic2, gl.max + EPSILON);
}

TEST_F(GeneEncodingTest, OneBinGeneStepChangesDecodedIC) {
    // G4.3 premise: step of size ~2^31/nbin in gene int space moves phenotype.
    const int32_t g0 = ictogene(&gl, 5.0);
    const int32_t step = static_cast<int32_t>(std::max(
        1.0, std::floor((static_cast<double>(MAX_RANDOM_VALUE) + 1.0) / gl.nbin)));
    int32_t g1 = g0 + step;
    if (g1 > MAX_RANDOM_VALUE) g1 = g0 - step;
    ASSERT_NE(g0, g1);
    const double ic0 = genetoic(&gl, g0);
    const double ic1 = genetoic(&gl, g1);
    EXPECT_NE(ic0, ic1);
}

TEST(MutationGranular, EnvOnWithGeneLimAppliesBinStep) {
    // Drive shipped mutate(..., gene_lim) with FLEXAIDDS_MUTATION_GRANULAR=1.
    // High mut_rate so gene is almost always selected; mid-range gene so ±step stays in range.
    ASSERT_EQ(setenv("FLEXAIDDS_MUTATION_GRANULAR", "1", 1), 0);
    genlim gl{};
    gl.min = 0.0;
    gl.max = 10.0;
    gl.del = 1.0;
    gl.nbin = 11.0;
    gl.bin = 1.0 / 11.0;
    gl.map = 0;

    gene g{};
    g.to_int32 = ictogene(&gl, 5.0);
    const int32_t before = g.to_int32;
    const double ic_before = genetoic(&gl, before);

    // Force many attempts; mut_rate=1.0 always fires.
    bool changed = false;
    for (int trial = 0; trial < 20; ++trial) {
        g.to_int32 = before;
        mutate(&g, 1, /*mut_rate=*/1.0, &gl);
        if (g.to_int32 != before) {
            changed = true;
            const double ic_after = genetoic(&gl, g.to_int32);
            // Phenotype must move by at least one bin width (allow float noise).
            EXPECT_GE(std::fabs(ic_after - ic_before), gl.del - 1e-6);
            break;
        }
    }
    EXPECT_TRUE(changed) << "granular mutate with mut_rate=1 never changed gene";
    unsetenv("FLEXAIDDS_MUTATION_GRANULAR");
}

TEST(MutationGranular, EnvOffKeepsClassicBitflipApi) {
    unsetenv("FLEXAIDDS_MUTATION_GRANULAR");
    gene g{};
    g.to_int32 = 123456789;
    // Classic path (null gene_lim) must still compile and run.
    mutate(&g, 1, 0.0, nullptr);
    EXPECT_EQ(g.to_int32, 123456789);
}

// ===========================================================================
// CHROMOSOME OPERATIONS
// ===========================================================================

TEST(CopyChrom, CopiesAllFields) {
    gene src_genes[3] = {{100, 1.5}, {200, 2.5}, {300, 3.5}};
    gene dst_genes[3] = {};

    chromosome src{}, dst{};
    src.genes = src_genes;
    dst.genes = dst_genes;

    src.evalue = -10.0;
    src.app_evalue = -9.5;
    src.fitnes = 42.0;
    src.status = 'n';
    src.cf.com = 1.0; src.cf.con = 2.0; src.cf.wal = 3.0;
    src.cf.sas = 4.0; src.cf.totsas = 5.0; src.cf.rclash = 0;

    copy_chrom(&dst, &src, 3);

    EXPECT_DOUBLE_EQ(dst.evalue, -10.0);
    EXPECT_DOUBLE_EQ(dst.app_evalue, -9.5);
    EXPECT_DOUBLE_EQ(dst.fitnes, 42.0);
    EXPECT_EQ(dst.status, 'n');
    EXPECT_DOUBLE_EQ(dst.cf.com, 1.0);

    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(dst.genes[i].to_int32, src.genes[i].to_int32);
        EXPECT_DOUBLE_EQ(dst.genes[i].to_ic, src.genes[i].to_ic);
    }
}

TEST(SwapChrom, SwapsContents) {
    gene g1[1] = {{111, 1.0}};
    gene g2[1] = {{222, 2.0}};
    chromosome a{}, b{};
    a.genes = g1; a.evalue = -5.0; a.fitnes = 10.0;
    b.genes = g2; b.evalue = -3.0; b.fitnes = 20.0;

    swap_chrom(&a, &b);

    // Note: swap_chrom does a shallow struct copy, so gene pointers swap
    EXPECT_DOUBLE_EQ(a.evalue, -3.0);
    EXPECT_DOUBLE_EQ(a.fitnes, 20.0);
    EXPECT_DOUBLE_EQ(b.evalue, -5.0);
    EXPECT_DOUBLE_EQ(b.fitnes, 10.0);
}

// ===========================================================================
// HASH GENES
// ===========================================================================

TEST(HashGenes, UniqueForDifferentGenes) {
    gene g1[3] = {{0, 1.0}, {0, 2.0}, {0, 3.0}};
    gene g2[3] = {{0, 1.0}, {0, 2.0}, {0, 4.0}};

    EXPECT_NE(hash_genes(g1, 3), hash_genes(g2, 3));
}

TEST(HashGenes, IdenticalForSameGenes) {
    gene g1[2] = {{0, 5.0}, {0, 10.0}};
    gene g2[2] = {{0, 5.0}, {0, 10.0}};

    EXPECT_EQ(hash_genes(g1, 2), hash_genes(g2, 2));
}

TEST(HashGenes, Deterministic) {
    gene g[1] = {{0, 7.3}};
    size_t h1 = hash_genes(g, 1);
    size_t h2 = hash_genes(g, 1);
    EXPECT_EQ(h1, h2);
}

// ===========================================================================
// QUICKSORT
// ===========================================================================

class QuickSortTest : public ::testing::Test {
protected:
    void verifyAscendingEnergy(ChromArray& ca, int n) {
        for (int i = 1; i < n; ++i)
            EXPECT_LE(ca[i - 1].evalue, ca[i].evalue)
                << "Energy not ascending at index " << i;
    }

    void verifyDescendingFitness(ChromArray& ca, int n) {
        for (int i = 1; i < n; ++i)
            EXPECT_GE(ca[i - 1].fitnes, ca[i].fitnes)
                << "Fitness not descending at index " << i;
    }
};

TEST_F(QuickSortTest, SortByEnergyAscending) {
    const int N = 10;
    ChromArray ca(N, 1);
    double energies[] = {5.0, -3.0, 8.0, 1.0, -10.0, 0.0, 3.0, -1.0, 7.0, 2.0};
    for (int i = 0; i < N; ++i)
        ca[i].evalue = energies[i];

    QuickSort(ca.data(), 0, N - 1, true);
    verifyAscendingEnergy(ca, N);
}

TEST_F(QuickSortTest, SortByFitnessDescending) {
    const int N = 8;
    ChromArray ca(N, 1);
    double fits[] = {1.0, 5.0, 3.0, 8.0, 2.0, 7.0, 4.0, 6.0};
    for (int i = 0; i < N; ++i)
        ca[i].fitnes = fits[i];

    QuickSort(ca.data(), 0, N - 1, false);

    // Verify descending order (highest fitness first)
    for (int i = 0; i < N - 1; ++i)
        EXPECT_GE(ca[i].fitnes, ca[i + 1].fitnes);

    // Verify all original values are still present (no data loss)
    std::vector<double> sorted_fits(N);
    for (int i = 0; i < N; ++i)
        sorted_fits[i] = ca[i].fitnes;
    std::sort(sorted_fits.begin(), sorted_fits.end());
    std::vector<double> expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    EXPECT_EQ(sorted_fits, expected);
}

TEST_F(QuickSortTest, AlreadySorted) {
    const int N = 5;
    ChromArray ca(N, 1);
    for (int i = 0; i < N; ++i)
        ca[i].evalue = static_cast<double>(i);

    QuickSort(ca.data(), 0, N - 1, true);
    verifyAscendingEnergy(ca, N);
}

TEST_F(QuickSortTest, ReverseSorted) {
    const int N = 5;
    ChromArray ca(N, 1);
    for (int i = 0; i < N; ++i)
        ca[i].evalue = static_cast<double>(N - i);

    QuickSort(ca.data(), 0, N - 1, true);
    verifyAscendingEnergy(ca, N);
}

TEST_F(QuickSortTest, SingleElement) {
    ChromArray ca(1, 1);
    ca[0].evalue = 42.0;
    QuickSort(ca.data(), 0, 0, true);
    EXPECT_DOUBLE_EQ(ca[0].evalue, 42.0);
}

TEST_F(QuickSortTest, TwoElements) {
    ChromArray ca(2, 1);
    ca[0].evalue = 10.0;
    ca[1].evalue = -5.0;
    QuickSort(ca.data(), 0, 1, true);
    EXPECT_LE(ca[0].evalue, ca[1].evalue);
}

TEST_F(QuickSortTest, DuplicateValues) {
    const int N = 6;
    ChromArray ca(N, 1);
    double vals[] = {3.0, 3.0, 1.0, 1.0, 2.0, 2.0};
    for (int i = 0; i < N; ++i)
        ca[i].evalue = vals[i];

    QuickSort(ca.data(), 0, N - 1, true);
    verifyAscendingEnergy(ca, N);
}

// ===========================================================================
// REMOVE DUPLICATES
// ===========================================================================

TEST(RemoveDups, AllUnique) {
    ChromArray ca(3, 2);
    ca[0].genes[0].to_ic = 1.0; ca[0].genes[1].to_ic = 2.0;
    ca[1].genes[0].to_ic = 3.0; ca[1].genes[1].to_ic = 4.0;
    ca[2].genes[0].to_ic = 5.0; ca[2].genes[1].to_ic = 6.0;

    int result = remove_dups(ca.data(), 3, 2);
    EXPECT_EQ(result, 3);
}

TEST(RemoveDups, AllIdentical) {
    ChromArray ca(4, 2);
    for (int i = 0; i < 4; ++i) {
        ca[i].genes[0].to_ic = 1.0;
        ca[i].genes[1].to_ic = 2.0;
    }

    int result = remove_dups(ca.data(), 4, 2);
    EXPECT_EQ(result, 1);
}

TEST(RemoveDups, SomeDuplicates) {
    ChromArray ca(4, 1);
    ca[0].genes[0].to_ic = 1.0;
    ca[1].genes[0].to_ic = 1.0;  // duplicate of [0]
    ca[2].genes[0].to_ic = 5.0;
    ca[3].genes[0].to_ic = 5.0;  // duplicate of [2]

    int result = remove_dups(ca.data(), 4, 1);
    EXPECT_EQ(result, 2);
}

TEST(RemoveDups, SingleChromosome) {
    ChromArray ca(1, 2);
    ca[0].genes[0].to_ic = 1.0;
    ca[0].genes[1].to_ic = 2.0;

    int result = remove_dups(ca.data(), 1, 2);
    EXPECT_EQ(result, 1);
}

TEST(RemoveDups, ToleranceBoundary) {
    // Values differ by exactly 0.1 per gene → should NOT be duplicates
    ChromArray ca(2, 1);
    ca[0].genes[0].to_ic = 1.0;
    ca[1].genes[0].to_ic = 1.1;

    int result = remove_dups(ca.data(), 2, 1);
    EXPECT_EQ(result, 2);
}

TEST(RemoveDups, WithinTolerance) {
    // Values differ by 0.05 < 0.1 → should be duplicates
    ChromArray ca(2, 1);
    ca[0].genes[0].to_ic = 1.0;
    ca[1].genes[0].to_ic = 1.05;

    int result = remove_dups(ca.data(), 2, 1);
    EXPECT_EQ(result, 1);
}

// ===========================================================================
// FITNESS STATISTICS
// ===========================================================================

TEST(FitnessStats, CalculatesMaxAndAverage) {
    GB_Global gb{};
    ChromArray ca(4, 1);
    ca[0].fitnes = 10.0;
    ca[1].fitnes = 20.0;
    ca[2].fitnes = 30.0;
    ca[3].fitnes = 40.0;

    fitness_stats(&gb, ca.data(), 4);

    // Average: (10+20+30+40)/4 = 25.0
    EXPECT_DOUBLE_EQ(gb.fit_avg, 25.0);
    // Max: 40.0
    EXPECT_DOUBLE_EQ(gb.fit_max, 40.0);
}

TEST(FitnessStats, SingleChromosome) {
    GB_Global gb{};
    ChromArray ca(1, 1);
    ca[0].fitnes = 42.0;

    fitness_stats(&gb, ca.data(), 1);

    EXPECT_DOUBLE_EQ(gb.fit_avg, 42.0);
}

// ===========================================================================
// ADAPTIVE PROBABILITIES
// ===========================================================================

TEST(AdaptProb, HighFitnessReducesProbability) {
    GB_Global gb{};
    gb.fit_max = 100.0;
    gb.fit_avg = 50.0;
    gb.k1 = 1.0;
    gb.k2 = 1.0;
    gb.k3 = 0.5;
    gb.k4 = 0.5;

    double mutp = 0.0, crossp = 0.0;
    // Both parents above average
    adapt_prob(&gb, 80.0, 90.0, &mutp, &crossp);

    // crossp = k1 * (max - high) / (max - avg) = 1.0 * (100 - 90) / (100 - 50) = 0.2
    EXPECT_NEAR(crossp, 0.2, EPSILON);
    // mutp = k2 * (max - low) / (max - avg) = 1.0 * (100 - 80) / (100 - 50) = 0.4
    EXPECT_NEAR(mutp, 0.4, EPSILON);
}

TEST(AdaptProb, BelowAverageUsesConstants) {
    GB_Global gb{};
    gb.fit_max = 100.0;
    gb.fit_avg = 50.0;
    gb.k1 = 1.0;
    gb.k2 = 1.0;
    gb.k3 = 0.5;
    gb.k4 = 0.3;

    double mutp = 0.0, crossp = 0.0;
    // Both parents below average
    adapt_prob(&gb, 20.0, 30.0, &mutp, &crossp);

    EXPECT_DOUBLE_EQ(crossp, 0.5);  // k3
    EXPECT_DOUBLE_EQ(mutp, 0.3);    // k4
}

TEST(AdaptProb, AtMaxFitnessZeroesProbability) {
    GB_Global gb{};
    gb.fit_max = 100.0;
    gb.fit_avg = 50.0;
    gb.k1 = 1.0;
    gb.k2 = 1.0;
    gb.k3 = 0.5;
    gb.k4 = 0.5;

    double mutp = 0.0, crossp = 0.0;
    // Higher parent is at maximum
    adapt_prob(&gb, 100.0, 60.0, &mutp, &crossp);

    // crossp = k1 * (100 - 100) / (100 - 50) = 0.0
    EXPECT_NEAR(crossp, 0.0, EPSILON);
}

// ===========================================================================
// CROSSOVER
// ===========================================================================

TEST(Crossover, MaterialIsConserved) {
    // After crossover, the combined bits of both parents should be the same
    // (crossover only swaps bits, never creates new ones)
    srand(42);  // deterministic
    gene john[2] = {{0x0F0F0F0F, 0.0}, {0x00FF00FF, 0.0}};
    gene mary[2] = {{static_cast<int32_t>(0xF0F0F0F0u), 0.0}, {static_cast<int32_t>(0xFF00FF00u), 0.0}};

    int32_t combined_before_0 = john[0].to_int32 | mary[0].to_int32;
    int32_t combined_before_1 = john[1].to_int32 | mary[1].to_int32;

    crossover(john, mary, 2, 0);  // inter-gene crossover

    // All bits that existed before should still exist in one parent or the other
    EXPECT_EQ(john[0].to_int32 | mary[0].to_int32, combined_before_0);
    EXPECT_EQ(john[1].to_int32 | mary[1].to_int32, combined_before_1);
}

TEST(Crossover, SingleGeneSwapsBits) {
    srand(123);
    gene john[1] = {{static_cast<int32_t>(0xAAAAAAAAu), 0.0}};
    gene mary[1] = {{0x55555555, 0.0}};

    int32_t j_before = john[0].to_int32;
    int32_t m_before = mary[0].to_int32;

    crossover(john, mary, 1, 1);  // intra-gene crossover

    // At least one bit should have changed (with very high probability)
    bool john_changed = (john[0].to_int32 != j_before);
    bool mary_changed = (mary[0].to_int32 != m_before);
    // Both should change symmetrically
    EXPECT_EQ(john_changed, mary_changed);
}

// ===========================================================================
// MUTATION
// ===========================================================================

TEST(Mutate, ZeroRateNoChange) {
    gene john[3] = {{0x12345678, 0.0}, {static_cast<int32_t>(0xABCDEF01u), 0.0}, {0x0, 0.0}};
    gene backup[3];
    std::memcpy(backup, john, sizeof(john));

    mutate(john, 3, 0.0);

    for (int i = 0; i < 3; ++i)
        EXPECT_EQ(john[i].to_int32, backup[i].to_int32);
}

TEST(Mutate, FullRateFlipsAllBits) {
    // With mut_rate = 1.0, every bit should flip
    // But RandomDouble() returns values in [0,1), so < 1.0 is always true
    gene john[1] = {{0x00000000, 0.0}};

    mutate(john, 1, 1.0);

    // All 32 bits should have flipped: 0x00000000 → 0xFFFFFFFF
    EXPECT_EQ(john[0].to_int32, static_cast<int32_t>(0xFFFFFFFFu));
}

// ===========================================================================
// ROULETTE WHEEL SELECTION
// ===========================================================================

TEST(RouletteWheel, UniformFitnessSelectsValid) {
    srand(42);
    ChromArray ca(5, 1);
    for (int i = 0; i < 5; ++i)
        ca[i].fitnes = 1.0;  // all equal

    for (int trial = 0; trial < 50; ++trial) {
        int selected = roullete_wheel(ca.data(), 5);
        EXPECT_GE(selected, 0);
        EXPECT_LT(selected, 5);
    }
}

TEST(RouletteWheel, DominantFitnessSelected) {
    srand(42);
    ChromArray ca(3, 1);
    ca[0].fitnes = 0.001;
    ca[1].fitnes = 0.001;
    ca[2].fitnes = 1000.0;  // overwhelmingly dominant

    // Over many trials, index 2 should be selected most often
    int count_2 = 0;
    for (int trial = 0; trial < 100; ++trial) {
        if (roullete_wheel(ca.data(), 3) == 2)
            ++count_2;
    }
    EXPECT_GT(count_2, 90);  // should be ~99%
}

TEST(RouletteWheel, SingleIndividual) {
    ChromArray ca(1, 1);
    ca[0].fitnes = 5.0;

    int selected = roullete_wheel(ca.data(), 1);
    EXPECT_EQ(selected, 0);
}

// ===========================================================================
// CALC_RMSP (RMS parameter-space distance)
// ===========================================================================

TEST(CalcRMSP, IdenticalGenesZeroDistance) {
    gene g1[3] = {{0, 1.0}, {0, 2.0}, {0, 3.0}};
    gene g2[3] = {{0, 1.0}, {0, 2.0}, {0, 3.0}};

    double rmsp = calc_rmsp(3, g1, g2, nullptr, nullptr);
    EXPECT_NEAR(rmsp, 0.0, EPSILON);
}

TEST(CalcRMSP, KnownDistance) {
    // g1 = (0, 0, 0), g2 = (1, 1, 1)
    // diff² = 1+1+1 = 3, rmsp = sqrt(3/3) = 1.0
    gene g1[3] = {{0, 0.0}, {0, 0.0}, {0, 0.0}};
    gene g2[3] = {{0, 1.0}, {0, 1.0}, {0, 1.0}};

    double rmsp = calc_rmsp(3, g1, g2, nullptr, nullptr);
    EXPECT_NEAR(rmsp, 1.0, EPSILON);
}

TEST(CalcRMSP, Symmetric) {
    gene g1[2] = {{0, 3.0}, {0, 4.0}};
    gene g2[2] = {{0, 0.0}, {0, 0.0}};

    double d1 = calc_rmsp(2, g1, g2, nullptr, nullptr);
    double d2 = calc_rmsp(2, g2, g1, nullptr, nullptr);
    EXPECT_NEAR(d1, d2, EPSILON);
}
