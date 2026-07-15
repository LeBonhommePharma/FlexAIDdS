// tests/test_cf_aggregator.cpp — pins the REAL get_cf_evalue() term set.
//
// Background: unit tests normally link tests/stubs.cpp, which reimplements
// get_cf_evalue() as a reduced sum (com+wal+sas+con+elec). That stub silently
// omits hbond, gist_desolv, metal_coord and entropy, so a bug that dropped a
// term from the production aggregator (LIB/ic2cf.cpp) would pass unnoticed.
//
// This target instead links the real LIB/ic2cf.cpp and asserts the exact
// additive term set. Each cfstr field is probed in isolation, and the aggregate
// is checked for linearity, so adding/removing/reweighting a term forces this
// test to be updated deliberately rather than regressing silently.
//
// Term set as of this pin (FA == nullptr, i.e. default gating):
//   included, weight 1.0 : com, con, wal, sas, elec, hbond, gist_desolv,
//                          metal_coord, entropy
//   gated off, weight 0.0: gist, h_rep, totsas, nor
// gist and h_rep are weighted by FA-supplied coefficients that default to zero
// when no FA_Global is provided.

#include <gtest/gtest.h>
#include "flexaid.h"

extern double get_cf_evalue(cfstr* cf, FA_Global* FA);

namespace {

double eval(const cfstr& c) {
    cfstr tmp = c;
    return get_cf_evalue(&tmp, nullptr);
}

// Helper: a zeroed cfstr with exactly one field set to `v`.
template <typename Setter>
double one_term(double v, Setter set) {
    cfstr c{};
    set(c, v);
    return eval(c);
}

} // namespace

TEST(CFAggregator, ZeroIsZero) {
    cfstr c{};
    EXPECT_DOUBLE_EQ(eval(c), 0.0);
}

TEST(CFAggregator, IncludedTermsContributeUnitWeight) {
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.com = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.con = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.wal = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.sas = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.elec = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.hbond = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.gist_desolv = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.metal_coord = v; }), 1.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.entropy = v; }), 1.0);
}

TEST(CFAggregator, GatedTermsAreExcludedWithoutFA) {
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.gist = v; }), 0.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.h_rep = v; }), 0.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.totsas = v; }), 0.0);
    EXPECT_DOUBLE_EQ(one_term(1.0, [](cfstr& c, double v){ c.nor = v; }), 0.0);
}

TEST(CFAggregator, IsLinearSumOfIncludedTerms) {
    // Distinct values per included term so a dropped or double-counted term
    // changes the total. Sum of these nine values is the expected aggregate.
    cfstr c{};
    c.com         = 2.0;
    c.con         = 3.0;
    c.wal         = 5.0;
    c.sas         = 7.0;
    c.elec        = 11.0;
    c.hbond       = 13.0;
    c.gist_desolv = 17.0;
    c.metal_coord = 19.0;
    c.entropy     = 23.0;
    // Gated-off terms set nonzero to confirm they stay excluded.
    c.gist  = 100.0;
    c.h_rep = 200.0;

    const double expected = 2.0 + 3.0 + 5.0 + 7.0 + 11.0 + 13.0 + 17.0 + 19.0 + 23.0;
    EXPECT_DOUBLE_EQ(eval(c), expected);
}
