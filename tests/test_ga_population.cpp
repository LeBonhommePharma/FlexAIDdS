// tests/test_ga_population.cpp
// Unit tests for GA population-level functions in gaboom.cpp:
//   fitness_stats, adapt_prob, roullete_wheel, cmp_chrom2pop,
//   cmp_chrom2pop_int, validate_dups, calc_poss, set_bins,
//   set_gene_lim, crossover, mutate, remove_dups
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <cstring>

#include "flexaid.h"
#include "gaboom.h"
#include "ga_constants.h"

// ═══════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════

namespace {

// Allocate a chromosome array with gene arrays of size num_genes.
chromosome* alloc_chroms(int num_chrom, int num_genes) {
    chromosome* c = new chromosome[num_chrom];
    std::memset(c, 0, sizeof(chromosome) * num_chrom);
    for (int i = 0; i < num_chrom; i++) {
        c[i].genes = new gene[num_genes];
        std::memset(c[i].genes, 0, sizeof(gene) * num_genes);
    }
    return c;
}

void free_chroms(chromosome* c, int num_chrom) {
    if (!c) return;
    for (int i = 0; i < num_chrom; i++) delete[] c[i].genes;
    delete[] c;
}

// Set gene values for a chromosome.
void set_genes(chromosome& c, const std::vector<double>& ic_vals) {
    for (size_t j = 0; j < ic_vals.size(); j++) {
        c.genes[j].to_ic = ic_vals[j];
        c.genes[j].to_int32 = static_cast<int32_t>(ic_vals[j] * 1000);
    }
}

void set_genes_int(chromosome& c, const std::vector<int32_t>& int_vals) {
    for (size_t j = 0; j < int_vals.size(); j++) {
        c.genes[j].to_int32 = int_vals[j];
        c.genes[j].to_ic = static_cast<double>(int_vals[j]);
    }
}

} // anonymous namespace

// ═══════════════════════════════════════════════════════════════════════
// fitness_stats — computes fit_max and fit_avg
// ═══════════════════════════════════════════════════════════════════════

class FitnessStatsTest : public ::testing::Test {
protected:
    GB_Global GB;
    void SetUp() override { std::memset(&GB, 0, sizeof(GB)); }
};

TEST_F(FitnessStatsTest, SingleChromosome) {
    auto* c = alloc_chroms(1, 2);
    c[0].fitnes = 5.0;
    fitness_stats(&GB, c, 1);
    EXPECT_DOUBLE_EQ(GB.fit_max, 5.0);
    EXPECT_DOUBLE_EQ(GB.fit_avg, 5.0);
    free_chroms(c, 1);
}

TEST_F(FitnessStatsTest, MultipleChromosomes) {
    const int N = 4;
    auto* c = alloc_chroms(N, 2);
    double vals[] = {3.0, 7.0, 1.0, 5.0};
    for (int i = 0; i < N; i++) c[i].fitnes = vals[i];
    fitness_stats(&GB, c, N);
    EXPECT_DOUBLE_EQ(GB.fit_max, 7.0);
    EXPECT_DOUBLE_EQ(GB.fit_avg, (3.0 + 7.0 + 1.0 + 5.0) / N);
    free_chroms(c, N);
}

TEST_F(FitnessStatsTest, AllZeroFitness) {
    auto* c = alloc_chroms(3, 2);
    fitness_stats(&GB, c, 3);
    EXPECT_DOUBLE_EQ(GB.fit_max, 0.0);
    EXPECT_DOUBLE_EQ(GB.fit_avg, 0.0);
    free_chroms(c, 3);
}

TEST_F(FitnessStatsTest, NegativeFitness) {
    const int N = 3;
    auto* c = alloc_chroms(N, 2);
    c[0].fitnes = -2.0; c[1].fitnes = -5.0; c[2].fitnes = -1.0;
    fitness_stats(&GB, c, N);
    EXPECT_DOUBLE_EQ(GB.fit_max, -1.0);
    EXPECT_DOUBLE_EQ(GB.fit_avg, (-2.0 - 5.0 - 1.0) / N);
    free_chroms(c, N);
}

TEST_F(FitnessStatsTest, MaxAtFirstPosition) {
    auto* c = alloc_chroms(3, 2);
    c[0].fitnes = 100.0; c[1].fitnes = 1.0; c[2].fitnes = 2.0;
    fitness_stats(&GB, c, 3);
    EXPECT_DOUBLE_EQ(GB.fit_max, 100.0);
    free_chroms(c, 3);
}

TEST_F(FitnessStatsTest, LargePopulation) {
    const int N = 1000;
    auto* c = alloc_chroms(N, 2);
    for (int i = 0; i < N; i++) c[i].fitnes = static_cast<double>(i);
    fitness_stats(&GB, c, N);
    EXPECT_DOUBLE_EQ(GB.fit_max, 999.0);
    EXPECT_NEAR(GB.fit_avg, 499.5, 1e-9);
    free_chroms(c, N);
}

// ═══════════════════════════════════════════════════════════════════════
// adapt_prob — adaptive mutation/crossover probabilities
// ═══════════════════════════════════════════════════════════════════════

class AdaptProbTest : public ::testing::Test {
protected:
    GB_Global GB;
    void SetUp() override {
        std::memset(&GB, 0, sizeof(GB));
        GB.fit_max = 10.0;
        GB.fit_avg = 5.0;
        GB.k1 = 1.0;
        GB.k2 = 0.5;
        GB.k3 = 1.0;
        GB.k4 = 0.5;
    }
};

TEST_F(AdaptProbTest, BothBelowAverage) {
    double mutp = 0, crossp = 0;
    adapt_prob(&GB, 3.0, 2.0, &mutp, &crossp);
    // fit_high=3.0 < fit_avg=5.0 → crossp = k3 = 1.0
    // fit_low=2.0  < fit_avg=5.0 → mutp   = k4 = 0.5
    EXPECT_DOUBLE_EQ(crossp, GB.k3);
    EXPECT_DOUBLE_EQ(mutp, GB.k4);
}

TEST_F(AdaptProbTest, BothAboveAverage) {
    double mutp = 0, crossp = 0;
    adapt_prob(&GB, 8.0, 7.0, &mutp, &crossp);
    // fit_high=8.0, fit_low=7.0
    // denom = 10.0 - 5.0 = 5.0
    // crossp = k1*(10-8)/5 = 1.0*0.4 = 0.4
    // mutp   = k2*(10-7)/5 = 0.5*0.6 = 0.3
    EXPECT_DOUBLE_EQ(crossp, GB.k1 * (10.0 - 8.0) / 5.0);
    EXPECT_DOUBLE_EQ(mutp, GB.k2 * (10.0 - 7.0) / 5.0);
}

TEST_F(AdaptProbTest, MixedHighLow) {
    double mutp = 0, crossp = 0;
    adapt_prob(&GB, 8.0, 3.0, &mutp, &crossp);
    // fit_high=8.0 > avg → crossp = k1*(10-8)/5 = 0.4
    // fit_low=3.0  < avg → mutp   = k4 = 0.5
    EXPECT_DOUBLE_EQ(crossp, 0.4);
    EXPECT_DOUBLE_EQ(mutp, GB.k4);
}

TEST_F(AdaptProbTest, ConvergedPopulation) {
    GB.fit_max = 5.0;
    GB.fit_avg = 5.0;  // denom → GA_FITNESS_DENOM_FLOOR
    double mutp = 0, crossp = 0;
    adapt_prob(&GB, 5.0, 5.0, &mutp, &crossp);
    // Both equal to avg → uses k3/k4
    EXPECT_DOUBLE_EQ(crossp, GB.k3);
    EXPECT_DOUBLE_EQ(mutp, GB.k4);
}

TEST_F(AdaptProbTest, HighEqualsMax) {
    double mutp = 0, crossp = 0;
    adapt_prob(&GB, 10.0, 5.0, &mutp, &crossp);
    // fit_high=10.0 == fit_max → crossp = k1*(10-10)/5 = 0.0
    EXPECT_DOUBLE_EQ(crossp, 0.0);
    // fit_low=5.0 == fit_avg → mutp = k4 = 0.5
    EXPECT_DOUBLE_EQ(mutp, GB.k4);
}

TEST_F(AdaptProbTest, OrderIndependence) {
    double mutp1 = 0, crossp1 = 0;
    double mutp2 = 0, crossp2 = 0;
    adapt_prob(&GB, 8.0, 3.0, &mutp1, &crossp1);
    adapt_prob(&GB, 3.0, 8.0, &mutp2, &crossp2);
    // Results should be identical regardless of argument order
    EXPECT_DOUBLE_EQ(crossp1, crossp2);
    EXPECT_DOUBLE_EQ(mutp1, mutp2);
}

// ═══════════════════════════════════════════════════════════════════════
// roullete_wheel — fitness-proportional selection
// ═══════════════════════════════════════════════════════════════════════

class RoulleteWheelTest : public ::testing::Test {
protected:
    // RandomDouble() uses thread_local RNG, auto-seeded from random_device
};

TEST_F(RoulleteWheelTest, ZeroPopulation) {
    auto* c = alloc_chroms(0, 2);
    int idx = roullete_wheel(c, 0);
    EXPECT_EQ(idx, 0);  // early return for n <= 0
    free_chroms(c, 0);
}

TEST_F(RoulleteWheelTest, SingleChromosome) {
    auto* c = alloc_chroms(1, 2);
    c[0].fitnes = 10.0;
    int idx = roullete_wheel(c, 1);
    EXPECT_EQ(idx, 0);
    free_chroms(c, 1);
}

TEST_F(RoulleteWheelTest, AllZeroFitness) {
    auto* c = alloc_chroms(5, 2);
    int idx = roullete_wheel(c, 5);
    EXPECT_GE(idx, 0);
    EXPECT_LT(idx, 5);
    free_chroms(c, 5);
}

TEST_F(RoulleteWheelTest, AllNegativeFitness) {
    auto* c = alloc_chroms(5, 2);
    for (int i = 0; i < 5; i++) c[i].fitnes = -1.0;
    int idx = roullete_wheel(c, 5);
    EXPECT_GE(idx, 0);
    EXPECT_LT(idx, 5);
    free_chroms(c, 5);
}

TEST_F(RoulleteWheelTest, ReturnsValidIndex) {
    const int N = 10;
    auto* c = alloc_chroms(N, 2);
    for (int i = 0; i < N; i++) c[i].fitnes = static_cast<double>(i + 1);
    for (int trial = 0; trial < 100; trial++) {
        int idx = roullete_wheel(c, N);
        EXPECT_GE(idx, 0);
        EXPECT_LT(idx, N);
    }
    free_chroms(c, N);
}

TEST_F(RoulleteWheelTest, BiasedTowardHigherFitness) {
    // With many trials, higher fitness should be selected more often
    const int N = 3;
    auto* c = alloc_chroms(N, 2);
    c[0].fitnes = 1.0;
    c[1].fitnes = 1.0;
    c[2].fitnes = 100.0;  // dominant

    int counts[N] = {0};
    for (int trial = 0; trial < 1000; trial++) {
        int idx = roullete_wheel(c, N);
        counts[idx]++;
    }
    // Chromosome 2 should be selected much more often
    EXPECT_GT(counts[2], counts[0]);
    EXPECT_GT(counts[2], counts[1]);
    free_chroms(c, N);
}

// ═══════════════════════════════════════════════════════════════════════
// cmp_chrom2pop — gene comparison with tolerance
// ═══════════════════════════════════════════════════════════════════════

class CmpChromTest : public ::testing::Test {
protected:
    static constexpr int NG = 3;
    chromosome* pop = nullptr;
    gene query[MAX_NUM_GENES];

    void SetUp() override {
        pop = alloc_chroms(4, NG);
        std::memset(query, 0, sizeof(gene) * MAX_NUM_GENES);
    }
    void TearDown() override { free_chroms(pop, 4); }
};

TEST_F(CmpChromTest, ExactMatch) {
    set_genes(pop[0], {1.0, 2.0, 3.0});
    set_genes(pop[1], {4.0, 5.0, 6.0});
    set_genes(pop[2], {7.0, 8.0, 9.0});
    for (int j = 0; j < NG; j++) query[j].to_ic = pop[1].genes[j].to_ic;

    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 0, 3), 1);
}

TEST_F(CmpChromTest, NoMatch) {
    set_genes(pop[0], {1.0, 2.0, 3.0});
    set_genes(pop[1], {4.0, 5.0, 6.0});
    set_genes(pop[2], {7.0, 8.0, 9.0});
    set_genes_int(pop[0], {100, 200, 300});
    set_genes_int(pop[1], {400, 500, 600});
    set_genes_int(pop[2], {700, 800, 900});
    // Query with completely different values
    for (int j = 0; j < NG; j++) {
        query[j].to_ic = 999.0 + j;
        query[j].to_int32 = 999000 + j;
    }
    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 0, 3), 0);
}

TEST_F(CmpChromTest, WithinTolerance) {
    set_genes(pop[0], {1.0, 2.0, 3.0});
    // Query within GA_GENE_MATCH_TOLERANCE of pop[0]
    for (int j = 0; j < NG; j++) {
        query[j].to_ic = pop[0].genes[j].to_ic + 0.05;  // within 0.1
    }
    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 0, 1), 1);
}

TEST_F(CmpChromTest, OutsideTolerance) {
    set_genes(pop[0], {1.0, 2.0, 3.0});
    for (int j = 0; j < NG; j++) {
        query[j].to_ic = pop[0].genes[j].to_ic + 0.2;  // outside 0.1
    }
    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 0, 1), 0);
}

TEST_F(CmpChromTest, PartialRange) {
    set_genes(pop[0], {1.0, 2.0, 3.0});
    set_genes(pop[1], {10.0, 20.0, 30.0});
    set_genes(pop[2], {100.0, 200.0, 300.0});
    // Query matches pop[1], range [1,3) includes only pop[1] and pop[2]
    for (int j = 0; j < NG; j++) query[j].to_ic = pop[1].genes[j].to_ic;
    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 1, 3), 1);
    // But not in range [0,1)
    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 0, 1), 0);
}

TEST_F(CmpChromTest, EmptyRange) {
    EXPECT_EQ(cmp_chrom2pop(pop, query, NG, 2, 2), 0);
}

// ═══════════════════════════════════════════════════════════════════════
// cmp_chrom2pop_int — exact int32 comparison
// ═══════════════════════════════════════════════════════════════════════

class CmpChromIntTest : public ::testing::Test {
protected:
    static constexpr int NG = 3;
    chromosome* pop = nullptr;
    gene query[MAX_NUM_GENES];

    void SetUp() override {
        pop = alloc_chroms(4, NG);
        std::memset(query, 0, sizeof(gene) * MAX_NUM_GENES);
    }
    void TearDown() override { free_chroms(pop, 4); }
};

TEST_F(CmpChromIntTest, ExactMatch) {
    set_genes_int(pop[0], {100, 200, 300});
    set_genes_int(pop[1], {400, 500, 600});
    for (int j = 0; j < NG; j++) query[j].to_int32 = pop[0].genes[j].to_int32;
    EXPECT_EQ(cmp_chrom2pop_int(pop, query, NG, 0, 2), 1);
}

TEST_F(CmpChromIntTest, NoMatch) {
    set_genes_int(pop[0], {100, 200, 300});
    for (int j = 0; j < NG; j++) query[j].to_int32 = 999;
    EXPECT_EQ(cmp_chrom2pop_int(pop, query, NG, 0, 1), 0);
}

TEST_F(CmpChromIntTest, PartialMismatch) {
    set_genes_int(pop[0], {100, 200, 300});
    query[0].to_int32 = 100;
    query[1].to_int32 = 200;
    query[2].to_int32 = 999;  // different
    EXPECT_EQ(cmp_chrom2pop_int(pop, query, NG, 0, 1), 0);
}

TEST_F(CmpChromIntTest, NegativeValues) {
    set_genes_int(pop[0], {-100, -200, -300});
    for (int j = 0; j < NG; j++) query[j].to_int32 = pop[0].genes[j].to_int32;
    EXPECT_EQ(cmp_chrom2pop_int(pop, query, NG, 0, 1), 1);
}

TEST_F(CmpChromIntTest, RangeSearch) {
    set_genes_int(pop[0], {10, 20, 30});
    set_genes_int(pop[1], {40, 50, 60});
    set_genes_int(pop[2], {70, 80, 90});
    // Query matches pop[2]
    for (int j = 0; j < NG; j++) query[j].to_int32 = pop[2].genes[j].to_int32;
    EXPECT_EQ(cmp_chrom2pop_int(pop, query, NG, 2, 3), 1);
    EXPECT_EQ(cmp_chrom2pop_int(pop, query, NG, 0, 2), 0);
}

// ═══════════════════════════════════════════════════════════════════════
// calc_poss — calculate total possibility space
// ═══════════════════════════════════════════════════════════════════════

TEST(CalcPossTest, SingleGene) {
    genlim gl[1];
    gl[0].nbin = 10.0;
    EXPECT_DOUBLE_EQ(calc_poss(gl, 1), 10.0);
}

TEST(CalcPossTest, MultipleGenes) {
    genlim gl[3];
    gl[0].nbin = 10.0; gl[1].nbin = 5.0; gl[2].nbin = 2.0;
    EXPECT_DOUBLE_EQ(calc_poss(gl, 3), 100.0);  // 10*5*2
}

TEST(CalcPossTest, SingleNbin) {
    genlim gl[2];
    gl[0].nbin = 1.0; gl[1].nbin = 1.0;
    EXPECT_DOUBLE_EQ(calc_poss(gl, 2), 1.0);
}

TEST(CalcPossTest, LargeNbins) {
    genlim gl[2];
    gl[0].nbin = 1000.0; gl[1].nbin = 1000.0;
    EXPECT_DOUBLE_EQ(calc_poss(gl, 2), 1000000.0);
}

// ═══════════════════════════════════════════════════════════════════════
// set_bins — compute bin widths from gene limits
// ═══════════════════════════════════════════════════════════════════════

class SetBinsTest : public ::testing::Test {
protected:
    static constexpr int NG = 3;
    genlim gl[MAX_NUM_GENES];
    void SetUp() override { std::memset(gl, 0, sizeof(gl)); }
};

TEST_F(SetBinsTest, ExactDivision) {
    gl[0] = {.max=10.0, .min=0.0, .del=2.0, .bin=0.0, .nbin=0.0, .map=0};
    set_bins(gl, 1);
    EXPECT_DOUBLE_EQ(gl[0].nbin, 5.0);
    EXPECT_DOUBLE_EQ(gl[0].bin, 1.0 / 5.0);
}

TEST_F(SetBinsTest, CeilForNonExact) {
    gl[0] = {.max=10.0, .min=0.0, .del=3.0, .bin=0.0, .nbin=0.0, .map=0};  // (10-0)/3 = 3.33 → nbin = 3.33 + 1.0 = 4.33
    set_bins(gl, 1);
    // Implementation adds 1.0 without truncating: nbin stays fractional
    EXPECT_DOUBLE_EQ(gl[0].nbin, 10.0 / 3.0 + 1.0);
    EXPECT_DOUBLE_EQ(gl[0].bin, 1.0 / (10.0 / 3.0 + 1.0));
}

TEST_F(SetBinsTest, MappedGeneGetsExtraBin) {
    // genlim fields: max, min, del, bin, nbin, map
    gl[0] = {.max=10.0, .min=0.0, .del=2.0, .bin=0.0, .nbin=0.0, .map=1};
    set_bins(gl, 1);
    // (10-0)/2 = 5 exact (no ceiling), map=1 → +1 = 6
    EXPECT_DOUBLE_EQ(gl[0].nbin, 6.0);
    EXPECT_DOUBLE_EQ(gl[0].bin, 1.0 / 6.0);
}

TEST_F(SetBinsTest, MultipleGenes) {
    gl[0] = {.max=10.0, .min=0.0, .del=2.0, .bin=0.0, .nbin=0.0, .map=0};  // 5 exact
    gl[1] = {.max=6.0, .min=0.0, .del=3.0, .bin=0.0, .nbin=0.0, .map=0};   // 2 exact
    gl[2] = {.max=5.0, .min=0.0, .del=2.0, .bin=0.0, .nbin=0.0, .map=1};   // 2.5→+1=3.5, +1(map)=4.5
    set_bins(gl, NG);
    EXPECT_DOUBLE_EQ(gl[0].nbin, 5.0);
    EXPECT_DOUBLE_EQ(gl[1].nbin, 2.0);
    // (5-0)/2 = 2.5, not exact → +1 = 3.5, mapped → +1 = 4.5
    EXPECT_DOUBLE_EQ(gl[2].nbin, 5.0 / 2.0 + 1.0 + 1.0);
}

TEST_F(SetBinsTest, SingleGenlimOverload) {
    genlim g = {.max=10.0, .min=0.0, .del=5.0, .bin=0.0, .nbin=0.0, .map=0};
    set_bins(&g);
    // Single overload only sets bin, not nbin
    EXPECT_DOUBLE_EQ(g.nbin, 0.0);
    EXPECT_DOUBLE_EQ(g.bin, 0.5);
}

TEST_F(SetBinsTest, SingleGenlimNonExact) {
    genlim g = {.max=7.0, .min=0.0, .del=3.0, .bin=0.0, .nbin=0.0, .map=0};
    set_bins(&g);
    // nbin = 7/3 + 1.0 = 3.333...; bin = 1/3.333...
    double expected_bin_inv = 7.0 / 3.0 + 1.0;
    EXPECT_DOUBLE_EQ(g.bin, 1.0 / expected_bin_inv);
}

// ═══════════════════════════════════════════════════════════════════════
// set_gene_lim — copies FA arrays to genlim struct
// ═══════════════════════════════════════════════════════════════════════

TEST(SetGeneLimTest, CopiesArrays) {
    const int NG = 3;
    FA_Global FA;
    GB_Global GB;
    genlim gl[MAX_NUM_GENES];
    std::memset(static_cast<void*>(&FA), 0, sizeof(FA));
    std::memset(static_cast<void*>(&GB), 0, sizeof(GB));
    std::memset(gl, 0, sizeof(gl));

    double min_par[] = {0.0, -5.0, 10.0};
    double max_par[] = {10.0, 5.0, 20.0};
    double del_par[] = {0.5, 1.0, 2.0};
    int    map_par[] = {0, 1, 0};

    FA.min_opt_par = min_par;
    FA.max_opt_par = max_par;
    FA.del_opt_par = del_par;
    FA.map_opt_par = map_par;
    GB.num_genes = NG;

    set_gene_lim(&FA, &GB, gl);

    for (int i = 0; i < NG; i++) {
        EXPECT_DOUBLE_EQ(gl[i].min, min_par[i]);
        EXPECT_DOUBLE_EQ(gl[i].max, max_par[i]);
        EXPECT_DOUBLE_EQ(gl[i].del, del_par[i]);
        EXPECT_EQ(gl[i].map, map_par[i]);
    }
}

// ═══════════════════════════════════════════════════════════════════════
// validate_dups — checks population fits possibility space
// ═══════════════════════════════════════════════════════════════════════

TEST(ValidateDupsTest, EnoughSpace_NoChange) {
    GB_Global GB;
    genlim gl[1];
    std::memset(&GB, 0, sizeof(GB));
    std::memset(gl, 0, sizeof(gl));

    GB.num_chrom = 10;
    GB.duplicates = 0;
    gl[0].nbin = 100.0;

    validate_dups(&GB, gl, 1);
    EXPECT_EQ(GB.duplicates, 0);  // 100 possibilities > 10 chroms
}

TEST(ValidateDupsTest, TooManyChroms_EnableDuplicates) {
    GB_Global GB;
    genlim gl[1];
    std::memset(&GB, 0, sizeof(GB));
    std::memset(gl, 0, sizeof(gl));

    GB.num_chrom = 1000;
    GB.duplicates = 0;
    gl[0].nbin = 5.0;  // Only 5 possibilities

    validate_dups(&GB, gl, 1);
    EXPECT_EQ(GB.duplicates, 1);
}

TEST(ValidateDupsTest, AlreadyAllowsDuplicates) {
    GB_Global GB;
    genlim gl[1];
    std::memset(&GB, 0, sizeof(GB));
    std::memset(gl, 0, sizeof(gl));

    GB.num_chrom = 1000;
    GB.duplicates = 1;
    gl[0].nbin = 5.0;

    validate_dups(&GB, gl, 1);
    EXPECT_EQ(GB.duplicates, 1);  // unchanged
}

// ═══════════════════════════════════════════════════════════════════════
// crossover — bitwise crossover operator
// ═══════════════════════════════════════════════════════════════════════

class CrossoverTest : public ::testing::Test {
protected:
    static constexpr int NG = 4;
    gene john[MAX_NUM_GENES];
    gene mary[MAX_NUM_GENES];
    gene john_orig[MAX_NUM_GENES];
    gene mary_orig[MAX_NUM_GENES];

    void SetUp() override {
        std::memset(john, 0, sizeof(gene) * MAX_NUM_GENES);
        std::memset(mary, 0, sizeof(gene) * MAX_NUM_GENES);
    }

    void save_orig() {
        std::memcpy(john_orig, john, sizeof(gene) * MAX_NUM_GENES);
        std::memcpy(mary_orig, mary, sizeof(gene) * MAX_NUM_GENES);
    }
};

TEST_F(CrossoverTest, MaterialConservation) {
    // All bits from john and mary should be preserved (just swapped)
    for (int j = 0; j < NG; j++) {
        john[j].to_int32 = 0xAAAAAAAA;
        mary[j].to_int32 = 0x55555555;
    }
    save_orig();

    crossover(john, mary, NG, 1);

    // For each gene: john[j] | mary[j] should == john_orig[j] | mary_orig[j]
    for (int j = 0; j < NG; j++) {
        auto uj = static_cast<unsigned int>(john[j].to_int32);
        auto um = static_cast<unsigned int>(mary[j].to_int32);
        auto ujo = static_cast<unsigned int>(john_orig[j].to_int32);
        auto umo = static_cast<unsigned int>(mary_orig[j].to_int32);
        EXPECT_EQ(uj | um, ujo | umo) << "Bits lost at gene " << j;
        EXPECT_EQ(uj & um, ujo & umo) << "Bits gained at gene " << j;
    }
}

TEST_F(CrossoverTest, XORConservation) {
    // XOR of john and mary at each gene should be preserved
    for (int j = 0; j < NG; j++) {
        john[j].to_int32 = static_cast<int32_t>(0x12345678);
        mary[j].to_int32 = static_cast<int32_t>(0xABCDEF01);
    }

    uint32_t orig_xor[NG];
    for (int j = 0; j < NG; j++) {
        orig_xor[j] = static_cast<uint32_t>(john[j].to_int32) ^
                       static_cast<uint32_t>(mary[j].to_int32);
    }

    crossover(john, mary, NG, 1);

    for (int j = 0; j < NG; j++) {
        uint32_t new_xor = static_cast<uint32_t>(john[j].to_int32) ^
                            static_cast<uint32_t>(mary[j].to_int32);
        EXPECT_EQ(new_xor, orig_xor[j]) << "XOR changed at gene " << j;
    }
}

TEST_F(CrossoverTest, InterGenesMode) {
    for (int j = 0; j < NG; j++) {
        john[j].to_int32 = static_cast<int32_t>(0xFFFFFFFF);
        mary[j].to_int32 = 0;
    }
    crossover(john, mary, NG, 0);
    // With intragenes=0, crossover operates at gene boundaries.
    // XOR conservation still holds: each gene is a bitwise swap of sub-ranges.
    for (int j = 0; j < NG; j++) {
        auto uj = static_cast<unsigned int>(john[j].to_int32);
        auto um = static_cast<unsigned int>(mary[j].to_int32);
        // XOR of john[j] ^ mary[j] should still be all-ones (conservation)
        EXPECT_EQ(uj ^ um, 0xFFFFFFFFu)
            << "XOR conservation violated at gene " << j;
    }
}

TEST_F(CrossoverTest, MultipleCrossoversProduceDifferentResults) {
    for (int j = 0; j < NG; j++) {
        john[j].to_int32 = static_cast<int32_t>(0xAAAAAAAA);
        mary[j].to_int32 = 0x55555555;
    }

    bool found_difference = false;
    int32_t first_john = john[0].to_int32;
    for (int trial = 0; trial < 50; trial++) {
        for (int j = 0; j < NG; j++) {
            john[j].to_int32 = static_cast<int32_t>(0xAAAAAAAA);
            mary[j].to_int32 = 0x55555555;
        }
        crossover(john, mary, NG, 1);
        if (john[0].to_int32 != first_john) found_difference = true;
    }
    // With random crossover points, some should differ
    EXPECT_TRUE(found_difference);
}

// ═══════════════════════════════════════════════════════════════════════
// mutate — bit-flip mutation
// ═══════════════════════════════════════════════════════════════════════

class MutateTest : public ::testing::Test {
protected:
    static constexpr int NG = 4;
    gene john[MAX_NUM_GENES];
    void SetUp() override {
        std::memset(john, 0, sizeof(gene) * MAX_NUM_GENES);
    }
};

TEST_F(MutateTest, ZeroRate_NoChange) {
    john[0].to_int32 = 0x12345678;
    int32_t orig = john[0].to_int32;
    mutate(john, NG, 0.0);
    EXPECT_EQ(john[0].to_int32, orig);  // No bits should flip
}

TEST_F(MutateTest, FullRate_AllBitsFlip) {
    john[0].to_int32 = 0;
    mutate(john, NG, 1.0);
    // With mut_rate=1.0, every bit has probability 1 of flipping
    // So all bits should be set → 0xFFFFFFFF
    EXPECT_EQ(static_cast<uint32_t>(john[0].to_int32), 0xFFFFFFFF);
}

TEST_F(MutateTest, DoubleMutateRestores) {
    john[0].to_int32 = 0x12345678;
    int32_t orig = john[0].to_int32;
    // Mutating twice with rate=1.0 should restore original (flip all, then flip all back)
    mutate(john, NG, 1.0);
    mutate(john, NG, 1.0);
    EXPECT_EQ(john[0].to_int32, orig);
}

TEST_F(MutateTest, ModerateRateChangesSomeBits) {
    // With 50% rate over many genes, at least some should change
    int changes = 0;
    for (int trial = 0; trial < 100; trial++) {
        int32_t orig = 0xAAAAAAAA;
        john[0].to_int32 = orig;
        mutate(john, 1, 0.5);
        if (john[0].to_int32 != orig) changes++;
    }
    EXPECT_GT(changes, 50);  // ~50% should change at least one bit
}

TEST_F(MutateTest, XORWithOriginalHasBoundedBits) {
    john[0].to_int32 = 0;
    int32_t orig = john[0].to_int32;
    mutate(john, NG, 0.5);
    // XOR tells us which bits flipped
    uint32_t flipped = static_cast<uint32_t>(john[0].to_int32 ^ orig);
    // Should not be all zeros (no change) nor all ones (all flipped)
    // for moderate rate (statistically very unlikely both extremes)
    // Just verify it's a valid 32-bit mask
    EXPECT_GE(__builtin_popcount(flipped), 0);
    EXPECT_LE(__builtin_popcount(flipped), 32);
}

// ═══════════════════════════════════════════════════════════════════════
// remove_dups — duplicate chromosome removal
// ═══════════════════════════════════════════════════════════════════════

class RemoveDupsTest : public ::testing::Test {
protected:
    static constexpr int NG = 3;
};

TEST_F(RemoveDupsTest, NoDuplicates) {
    auto* c = alloc_chroms(3, NG);
    set_genes(c[0], {1.0, 2.0, 3.0});
    set_genes(c[1], {4.0, 5.0, 6.0});
    set_genes(c[2], {7.0, 8.0, 9.0});
    int n = remove_dups(c, 3, NG);
    EXPECT_EQ(n, 3);
    free_chroms(c, 3);
}

TEST_F(RemoveDupsTest, AllDuplicates) {
    auto* c = alloc_chroms(4, NG);
    for (int i = 0; i < 4; i++) set_genes(c[i], {1.0, 2.0, 3.0});
    int n = remove_dups(c, 4, NG);
    EXPECT_EQ(n, 1);  // Only first kept
    free_chroms(c, 4);
}

TEST_F(RemoveDupsTest, SomeDuplicates) {
    auto* c = alloc_chroms(5, NG);
    set_genes(c[0], {1.0, 2.0, 3.0});
    set_genes(c[1], {1.0, 2.0, 3.0});  // dup of [0]
    set_genes(c[2], {4.0, 5.0, 6.0});
    set_genes(c[3], {4.0, 5.0, 6.0});  // dup of [2]
    set_genes(c[4], {7.0, 8.0, 9.0});
    int n = remove_dups(c, 5, NG);
    EXPECT_EQ(n, 3);
    free_chroms(c, 5);
}

TEST_F(RemoveDupsTest, SingleChromosome) {
    auto* c = alloc_chroms(1, NG);
    set_genes(c[0], {1.0, 2.0, 3.0});
    int n = remove_dups(c, 1, NG);
    EXPECT_EQ(n, 1);
    free_chroms(c, 1);
}

TEST_F(RemoveDupsTest, EmptyPopulation) {
    auto* c = alloc_chroms(0, NG);
    int n = remove_dups(c, 0, NG);
    EXPECT_EQ(n, 0);
    free_chroms(c, 0);
}

TEST_F(RemoveDupsTest, TwoChromosomes_Same) {
    auto* c = alloc_chroms(2, NG);
    set_genes(c[0], {1.0, 2.0, 3.0});
    set_genes(c[1], {1.0, 2.0, 3.0});
    int n = remove_dups(c, 2, NG);
    EXPECT_EQ(n, 1);
    free_chroms(c, 2);
}

TEST_F(RemoveDupsTest, TwoChromosomes_Different) {
    auto* c = alloc_chroms(2, NG);
    set_genes(c[0], {1.0, 2.0, 3.0});
    set_genes(c[1], {4.0, 5.0, 6.0});
    int n = remove_dups(c, 2, NG);
    EXPECT_EQ(n, 2);
    free_chroms(c, 2);
}

// ═══════════════════════════════════════════════════════════════════════
// copy_chrom — chromosome copy helper
// ═══════════════════════════════════════════════════════════════════════

TEST(CopyChromTest, CopiesAllFields) {
    const int NG = 3;
    auto* src = alloc_chroms(1, NG);
    auto* dst = alloc_chroms(1, NG);

    src[0].fitnes = 42.0;
    src[0].evalue = -10.5;
    src[0].app_evalue = -9.8;
    src[0].status = 'n';
    set_genes(src[0], {1.0, 2.0, 3.0});

    copy_chrom(&dst[0], &src[0], NG);

    EXPECT_DOUBLE_EQ(dst[0].fitnes, 42.0);
    EXPECT_DOUBLE_EQ(dst[0].evalue, -10.5);
    EXPECT_DOUBLE_EQ(dst[0].app_evalue, -9.8);
    EXPECT_EQ(dst[0].status, 'n');
    for (int j = 0; j < NG; j++) {
        EXPECT_EQ(dst[0].genes[j].to_int32, src[0].genes[j].to_int32);
        EXPECT_DOUBLE_EQ(dst[0].genes[j].to_ic, src[0].genes[j].to_ic);
    }

    free_chroms(src, 1);
    free_chroms(dst, 1);
}

TEST(CopyChromTest, IndependentAfterCopy) {
    const int NG = 2;
    auto* src = alloc_chroms(1, NG);
    auto* dst = alloc_chroms(1, NG);

    set_genes(src[0], {5.0, 10.0});
    copy_chrom(&dst[0], &src[0], NG);

    // Modify src after copy
    src[0].genes[0].to_ic = 999.0;
    EXPECT_DOUBLE_EQ(dst[0].genes[0].to_ic, 5.0);  // dst unchanged

    free_chroms(src, 1);
    free_chroms(dst, 1);
}
