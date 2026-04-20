// tests/test_ga_validation.cpp — GA scoring validation tests
// Tests: cfstr scoring functions, batch result structure, scoring invariants
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <cstring>
#include <numeric>
#include <algorithm>

#include "flexaid.h"
#include "gaboom.h"
#include "VoronoiCFBatch.h"

namespace {

// ═══════════════════════════════════════════════════════════════════════
// cfstr + scoring function tests
// ═══════════════════════════════════════════════════════════════════════

class CFScoreTest : public ::testing::Test {
protected:
    cfstr cf;

    void SetUp() override {
        std::memset(&cf, 0, sizeof(cfstr));
    }
};

// get_apparent_cf_evalue returns com + wal + sas + elec
TEST_F(CFScoreTest, ApparentEvalueFormula) {
    cf.com = -10.5;
    cf.wal = 2.0;
    cf.sas = 1.5;
    cf.elec = 0.3;
    cf.con = 99.0;  // should NOT contribute

    double app = get_apparent_cf_evalue(&cf);
    EXPECT_DOUBLE_EQ(app, -10.5 + 2.0 + 1.5 + 0.3);
}

// get_cf_evalue returns com + wal + sas + con + elec
TEST_F(CFScoreTest, FullEvalueFormula) {
    cf.com = -10.5;
    cf.wal = 2.0;
    cf.sas = 1.5;
    cf.con = 3.0;
    cf.elec = 0.3;

    double full = get_cf_evalue(&cf);
    EXPECT_DOUBLE_EQ(full, -10.5 + 2.0 + 1.5 + 3.0 + 0.3);
}

// Constraint difference: get_cf_evalue - get_apparent_cf_evalue == con
TEST_F(CFScoreTest, ConstraintDifference) {
    cf.com = -8.0;
    cf.wal = 1.0;
    cf.sas = 0.5;
    cf.con = 4.0;

    double diff = get_cf_evalue(&cf) - get_apparent_cf_evalue(&cf);
    EXPECT_DOUBLE_EQ(diff, cf.con);
}

// Zero-initialized cfstr gives zero scores
TEST_F(CFScoreTest, ZeroCfstrGivesZero) {
    EXPECT_DOUBLE_EQ(get_apparent_cf_evalue(&cf), 0.0);
    EXPECT_DOUBLE_EQ(get_cf_evalue(&cf), 0.0);
}

// Both functions return finite results for extreme values
TEST_F(CFScoreTest, FiniteForExtremeValues) {
    cf.com = -1e6;
    cf.wal = 1e6;
    cf.sas = 1e3;
    cf.con = -1e3;

    EXPECT_TRUE(std::isfinite(get_apparent_cf_evalue(&cf)));
    EXPECT_TRUE(std::isfinite(get_cf_evalue(&cf)));
}

// No wall/constraint/elec penalty: apparent == com + sas, full == com + sas
TEST_F(CFScoreTest, NoWallPenalty) {
    cf.com = -15.0;
    cf.wal = 0.0;
    cf.sas = 2.5;
    cf.con = 0.0;
    cf.elec = 0.0;

    EXPECT_DOUBLE_EQ(get_apparent_cf_evalue(&cf), cf.com + cf.sas);
    EXPECT_DOUBLE_EQ(get_cf_evalue(&cf), cf.com + cf.sas);
}

// ═══════════════════════════════════════════════════════════════════════
// cfstr determinism: same input → same output
// ═══════════════════════════════════════════════════════════════════════

TEST_F(CFScoreTest, Deterministic) {
    cf.com = -12.3;
    cf.wal = 0.7;
    cf.sas = 1.1;
    cf.con = 0.2;

    std::vector<double> results;
    for (int i = 0; i < 5; ++i)
        results.push_back(get_apparent_cf_evalue(&cf));

    for (size_t i = 1; i < results.size(); ++i)
        EXPECT_EQ(results[0], results[i])
            << "Scoring must be deterministic";
}

// ═══════════════════════════════════════════════════════════════════════
// BatchResult structure tests
// ═══════════════════════════════════════════════════════════════════════

TEST(BatchResultTest, DefaultConstructionValid) {
    voronoi_cf::BatchResult result;
    EXPECT_TRUE(result.cf.empty());
    EXPECT_TRUE(result.app_evalue.empty());
    EXPECT_EQ(result.wall_ms, 0.0);
}

TEST(BatchResultTest, ResizePreservesSize) {
    voronoi_cf::BatchResult result;
    const size_t N = 100;
    result.cf.resize(N);
    result.app_evalue.resize(N);

    EXPECT_EQ(result.cf.size(), N);
    EXPECT_EQ(result.app_evalue.size(), N);
}

TEST(BatchResultTest, CfAndAppEvalueConsistency) {
    // Verify that batch result cfstr → app_evalue is consistent
    voronoi_cf::BatchResult result;
    const int N = 10;
    result.cf.resize(N);
    result.app_evalue.resize(N);

    for (int i = 0; i < N; ++i) {
        result.cf[i].com = -(double)(i + 1);
        result.cf[i].wal = 0.1 * i;
        result.cf[i].sas = 0.05 * i;
        result.cf[i].con = 0.01 * i;
        result.app_evalue[i] = get_apparent_cf_evalue(&result.cf[i]);
    }

    // Verify consistency: apparent = com + wal + sas + elec
    for (int i = 0; i < N; ++i) {
        double expected = result.cf[i].com + result.cf[i].wal + result.cf[i].sas + result.cf[i].elec;
        EXPECT_DOUBLE_EQ(result.app_evalue[i], expected)
            << "Apparent evalue mismatch at index " << i;
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Scoring ordering invariant: more negative com → better score
// ═══════════════════════════════════════════════════════════════════════

TEST(ScoringInvariant, MoreNegativeComIsBetter) {
    cfstr good, bad;
    std::memset(&good, 0, sizeof(cfstr));
    std::memset(&bad, 0, sizeof(cfstr));

    good.com = -20.0;  // strong complementarity
    bad.com  = -5.0;   // weak complementarity

    EXPECT_LT(get_apparent_cf_evalue(&good), get_apparent_cf_evalue(&bad))
        << "More negative com should give lower (better) apparent evalue";
    EXPECT_LT(get_cf_evalue(&good), get_cf_evalue(&bad))
        << "More negative com should give lower (better) full evalue";
}

TEST(ScoringInvariant, WallPenaltyWorsensScore) {
    cfstr no_clash, clash;
    std::memset(&no_clash, 0, sizeof(cfstr));
    std::memset(&clash, 0, sizeof(cfstr));

    no_clash.com = -10.0;
    no_clash.wal = 0.0;

    clash.com = -10.0;
    clash.wal = 5.0;  // steric clash penalty

    EXPECT_LT(get_apparent_cf_evalue(&no_clash), get_apparent_cf_evalue(&clash))
        << "Wall penalty should worsen the score";
}

// ═══════════════════════════════════════════════════════════════════════
// eval_span: gene clamping correctness
// ═══════════════════════════════════════════════════════════════════════

// eval_span clamps gene IC values to [gene_lim.min, gene_lim.max].
// We verify this indirectly by constructing a trivial scoring function
// that returns the sum of IC values in the cf.com field.

static cfstr test_sum_function(
    FA_Global* /*FA*/, VC_Global* /*VC*/,
    atom* /*atoms*/, resid* /*residue*/, gridpoint* /*cleftgrid*/,
    int n_genes, double* icv)
{
    cfstr result;
    std::memset(&result, 0, sizeof(cfstr));
    double sum = 0.0;
    for (int i = 0; i < n_genes; ++i)
        sum += icv[i];
    result.com = sum;
    return result;
}

TEST(EvalSpanTest, ClampsBeyondMax) {
    FA_Global fa;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
    std::memset(&fa, 0, sizeof(FA_Global));
    GB_Global gb;
    std::memset(&gb, 0, sizeof(GB_Global));
    gb.num_genes = 2;
    VC_Global vc;
    std::memset(&vc, 0, sizeof(VC_Global));
#pragma clang diagnostic pop

    genlim gl[2];
    gl[0].min = -10.0; gl[0].max = 10.0; gl[0].del = 1.0; gl[0].nbin = 20; gl[0].bin = 0.05;
    gl[1].min = 0.0;   gl[1].max = 5.0;  gl[1].del = 0.5; gl[1].nbin = 10; gl[1].bin = 0.1;

    // Gene values beyond max
    gene genes[2];
    genes[0].to_ic = 999.0;   // well above max of 10
    genes[1].to_ic = 100.0;   // well above max of 5

    atom a;
    std::memset(&a, 0, sizeof(atom));
    resid r;
    std::memset(&r, 0, sizeof(resid));

    cfstr result = voronoi_cf::eval_span(
        &fa, &gb, &vc,
        std::span<const genlim>(gl, 2),
        std::span<atom>(&a, 1),
        std::span<resid>(&r, 1),
        nullptr,
        std::span<const gene>(genes, 2),
        test_sum_function
    );

    // IC values should be clamped to max: 10.0 + 5.0 = 15.0
    EXPECT_DOUBLE_EQ(result.com, 15.0);
}

TEST(EvalSpanTest, ClampsBelowMin) {
    FA_Global fa;
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnontrivial-memcall"
    std::memset(&fa, 0, sizeof(FA_Global));
    GB_Global gb;
    std::memset(&gb, 0, sizeof(GB_Global));
    gb.num_genes = 2;
    VC_Global vc;
    std::memset(&vc, 0, sizeof(VC_Global));
#pragma clang diagnostic pop

    genlim gl[2];
    gl[0].min = -10.0; gl[0].max = 10.0; gl[0].del = 1.0; gl[0].nbin = 20; gl[0].bin = 0.05;
    gl[1].min = 0.0;   gl[1].max = 5.0;  gl[1].del = 0.5; gl[1].nbin = 10; gl[1].bin = 0.1;

    gene genes[2];
    genes[0].to_ic = -999.0;  // below min of -10
    genes[1].to_ic = -999.0;  // below min of 0

    atom a;
    std::memset(&a, 0, sizeof(atom));
    resid r;
    std::memset(&r, 0, sizeof(resid));

    cfstr result = voronoi_cf::eval_span(
        &fa, &gb, &vc,
        std::span<const genlim>(gl, 2),
        std::span<atom>(&a, 1),
        std::span<resid>(&r, 1),
        nullptr,
        std::span<const gene>(genes, 2),
        test_sum_function
    );

    // IC values should be clamped to min: -10.0 + 0.0 = -10.0
    EXPECT_DOUBLE_EQ(result.com, -10.0);
}

TEST(EvalSpanTest, WithinBoundsPassthrough) {
    FA_Global fa;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wnontrivial-memcall"
    std::memset(&fa, 0, sizeof(FA_Global));
    GB_Global gb;
    std::memset(&gb, 0, sizeof(GB_Global));
    #pragma clang diagnostic pop
    gb.num_genes = 3;
    VC_Global vc;
    std::memset(&vc, 0, sizeof(VC_Global));

    genlim gl[3];
    for (int i = 0; i < 3; ++i) {
        gl[i].min = -180.0; gl[i].max = 180.0;
        gl[i].del = 1.0; gl[i].nbin = 360; gl[i].bin = 1.0/360.0;
    }

    gene genes[3];
    genes[0].to_ic = 45.0;
    genes[1].to_ic = -90.0;
    genes[2].to_ic = 0.0;

    atom a;
    std::memset(&a, 0, sizeof(atom));
    resid r;
    std::memset(&r, 0, sizeof(resid));

    cfstr result = voronoi_cf::eval_span(
        &fa, &gb, &vc,
        std::span<const genlim>(gl, 3),
        std::span<atom>(&a, 1),
        std::span<resid>(&r, 1),
        nullptr,
        std::span<const gene>(genes, 3),
        test_sum_function
    );

    // Values within bounds: 45 + (-90) + 0 = -45.0
    EXPECT_DOUBLE_EQ(result.com, -45.0);
}

}  // namespace anonymous
