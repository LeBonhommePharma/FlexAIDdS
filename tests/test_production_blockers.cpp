// tests/test_production_blockers.cpp
// Focused regression tests for production-blocker fixes and guardrails.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/gaboom.h"
#include "GPUContextPool.h"

#include <cmath>
#include <vector>

namespace {

struct ChromArray {
    std::vector<chromosome> chroms;
    std::vector<std::vector<gene>> gene_store;

    ChromArray(int n, int ng) : chroms(n), gene_store(n, std::vector<gene>(ng)) {
        for (int i = 0; i < n; ++i) chroms[i].genes = gene_store[i].data();
    }

    chromosome* data() { return chroms.data(); }
    chromosome& operator[](int i) { return chroms[i]; }
};

} // namespace

TEST(ProductionBlockers, RouletteWheelAlwaysReturnsInBounds) {
    ChromArray ca(4, 1);
    ca[0].fitnes = 1.0;
    ca[1].fitnes = 2.0;
    ca[2].fitnes = 3.0;
    ca[3].fitnes = 4.0;

    for (int i = 0; i < 1000; ++i) {
        int idx = roullete_wheel(ca.data(), 4);
        EXPECT_GE(idx, 0);
        EXPECT_LT(idx, 4);
    }
}

TEST(ProductionBlockers, RouletteWheelHandlesNonPositiveTotalFitness) {
    ChromArray ca(3, 1);
    ca[0].fitnes = 0.0;
    ca[1].fitnes = 0.0;
    ca[2].fitnes = 0.0;

    for (int i = 0; i < 100; ++i) {
        int idx = roullete_wheel(ca.data(), 3);
        EXPECT_GE(idx, 0);
        EXPECT_LT(idx, 3);
    }
}

TEST(ProductionBlockers, AdaptProbRemainsFiniteWhenPopulationConverged) {
    GB_Global gb{};
    gb.fit_max = 10.0;
    gb.fit_avg = 10.0;
    gb.k1 = 1.0;
    gb.k2 = 1.0;
    gb.k3 = 0.5;
    gb.k4 = 0.25;

    double mutp = 0.0;
    double crossp = 0.0;
    adapt_prob(&gb, 10.0, 10.0, &mutp, &crossp);

    EXPECT_TRUE(std::isfinite(mutp));
    EXPECT_TRUE(std::isfinite(crossp));
    EXPECT_GE(mutp, 0.0);
    EXPECT_GE(crossp, 0.0);
}

TEST(ProductionBlockers, GPUContextPoolSingletonConstructs) {
    auto& pool1 = GPUContextPool::instance();
    auto& pool2 = GPUContextPool::instance();
    EXPECT_EQ(&pool1, &pool2);
}
