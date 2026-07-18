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
#include "Vcontacts.h"

#include <vector>

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

TEST(CFAggregator, ElecAndGistDesolvAreNotSilentZeros) {
    // Regression: production ic2cf previously zeroed elec/gist_desolv then
    // failed to sum them from optres. get_cf_evalue must still count them.
    cfstr only_elec{};
    only_elec.elec = 42.0;
    EXPECT_DOUBLE_EQ(eval(only_elec), 42.0);

    cfstr only_gist{};
    only_gist.gist_desolv = -7.5;
    EXPECT_DOUBLE_EQ(eval(only_gist), -7.5);

    cfstr both{};
    both.elec = 10.0;
    both.gist_desolv = 5.0;
    both.com = 1.0;
    EXPECT_DOUBLE_EQ(eval(both), 16.0);
}

// Mirrors LIB/ic2cf.cpp per-residue reduction (sum optres[i].cf.* into cf).
static cfstr reduce_optres_like_ic2cf(const std::vector<cfstr>& per_res) {
    cfstr cf{};  // value-init like production ic2cf
    for (const auto& r : per_res) {
        cf.com += r.com;
        cf.wal += r.wal;
        cf.sas += r.sas;
        cf.con += r.con;
        cf.elec += r.elec;
        cf.gist_desolv += r.gist_desolv;
        cf.metal_coord += r.metal_coord;
        cf.hbond += r.hbond;
        cf.entropy += r.entropy;
    }
    return cf;
}

TEST(CFAggregator, PerResidueReductionMatchesGetCfEvalue) {
    // Exercise the actual ic2cf accumulation pattern, not only get_cf_evalue
    // on a hand-built aggregate. Distinct residues contribute distinct terms.
    std::vector<cfstr> residues(3);
    residues[0].com = 2.0;
    residues[0].elec = 1.5;
    residues[1].wal = 3.0;
    residues[1].hbond = -0.5;
    residues[2].sas = 4.0;
    residues[2].gist_desolv = 0.25;
    residues[2].metal_coord = 0.1;
    residues[2].entropy = 0.05;
    // Gated fields present but must not affect get_cf_evalue without FA.
    residues[0].gist = 99.0;
    residues[1].h_rep = 88.0;

    cfstr reduced = reduce_optres_like_ic2cf(residues);
    EXPECT_DOUBLE_EQ(reduced.com, 2.0);
    EXPECT_DOUBLE_EQ(reduced.wal, 3.0);
    EXPECT_DOUBLE_EQ(reduced.sas, 4.0);
    EXPECT_DOUBLE_EQ(reduced.elec, 1.5);
    EXPECT_DOUBLE_EQ(reduced.hbond, -0.5);
    EXPECT_DOUBLE_EQ(reduced.gist_desolv, 0.25);
    EXPECT_DOUBLE_EQ(reduced.metal_coord, 0.1);
    EXPECT_DOUBLE_EQ(reduced.entropy, 0.05);
    // Value-init: untouched fields stay 0.
    EXPECT_DOUBLE_EQ(reduced.con, 0.0);
    EXPECT_DOUBLE_EQ(reduced.gist, 0.0);

    const double expected =
        2.0 + 3.0 + 4.0 + 1.5 + (-0.5) + 0.25 + 0.1 + 0.05;
    EXPECT_DOUBLE_EQ(eval(reduced), expected);
}

TEST(CFAggregator, ValueInitIsZeroNotGarbage) {
    cfstr c{};
    EXPECT_DOUBLE_EQ(c.com, 0.0);
    EXPECT_DOUBLE_EQ(c.elec, 0.0);
    EXPECT_DOUBLE_EQ(c.gist_desolv, 0.0);
    EXPECT_DOUBLE_EQ(c.metal_coord, 0.0);
    EXPECT_DOUBLE_EQ(c.hbond, 0.0);
    EXPECT_DOUBLE_EQ(c.entropy, 0.0);
    EXPECT_EQ(c.rclash, 0);
}

// Controllable stubs (tests/cf_aggregator_stubs.cpp)
namespace ic2cf_test_hooks {
extern int g_buildcc_fail_next;
extern int g_vcfunction_error_next;
extern int g_buildcc_calls;
extern int g_vcfunction_calls;
extern bool g_buildcc_scribble_coords;
extern float g_scribble_value;
}
extern cfstr ic2cf(FA_Global*, VC_Global*, atom*, resid*, gridpoint*, int, double*);

namespace {

struct MinimalIc2cfFixture {
    FA_Global FA{};
    VC_Global VC{};
    atom atoms[4]{};
    resid residue[2]{};
    OptRes optres[1]{};
    int mov_buf[1]{1};
    double icv[1]{0.0};
    // map_par needs at least npar entries; keep one dummy gene
    // (FA.map_par is pointer — allocate statically).
    static constexpr int kNpar = 0;

    MinimalIc2cfFixture() {
        FA.nors = 1;
        FA.nmov[0] = 1;
        FA.mov[0] = mov_buf;
        FA.num_optres = 1;
        FA.optres = optres;
        optres[0].rnum = 1;
        optres[0].type = 1;  // ligand
        optres[0].cf = cfstr{};
        FA.res_cnt = 1;
        FA.useflexdee = 0;
        FA.ring_flex_active = 0;
        FA.ori[0] = 1.0f;
        FA.ori[1] = 2.0f;
        FA.ori[2] = 3.0f;
        FA.globalmin[0] = FA.globalmin[1] = FA.globalmin[2] = -10.0f;
        FA.globalmax[0] = FA.globalmax[1] = FA.globalmax[2] = 10.0f;
        atoms[1].coor[0] = 0.5f;
        atoms[1].coor[1] = 0.5f;
        atoms[1].coor[2] = 0.5f;
        atoms[1].ofres = 1;
        residue[1].rot = 0;
        FA.map_par = nullptr;
    }
};

}  // namespace

TEST(Ic2cfContamination, BuildccFailureRestoresBaselineThenNextEvalClean) {
    ic2cf_test_hooks::g_buildcc_fail_next = 0;
    ic2cf_test_hooks::g_vcfunction_error_next = 0;
    ic2cf_test_hooks::g_buildcc_scribble_coords = false;
    ic2cf_test_hooks::g_buildcc_calls = 0;
    ic2cf_test_hooks::g_vcfunction_calls = 0;

    MinimalIc2cfFixture fx;
    const float x0 = fx.atoms[1].coor[0];
    const float ori0 = fx.FA.ori[0];
    fx.residue[1].rot = 3;

    // Eval 1: buildcc fails → must restore ori/atoms/rot and return penalty.
    ic2cf_test_hooks::g_buildcc_fail_next = 1;
    cfstr bad = ic2cf(&fx.FA, &fx.VC, fx.atoms, fx.residue, nullptr, 0, fx.icv);
    EXPECT_EQ(bad.rclash, 1);
    EXPECT_DOUBLE_EQ(bad.wal, 1.0e12);
    EXPECT_FLOAT_EQ(fx.atoms[1].coor[0], x0);
    EXPECT_FLOAT_EQ(fx.FA.ori[0], ori0);
    EXPECT_EQ(fx.residue[1].rot, 3);

    // Eval 2: consecutive call must see clean baseline (no contamination).
    ic2cf_test_hooks::g_buildcc_fail_next = 0;
    ic2cf_test_hooks::g_vcfunction_error_next = 0;
    cfstr ok = ic2cf(&fx.FA, &fx.VC, fx.atoms, fx.residue, nullptr, 0, fx.icv);
    EXPECT_EQ(ok.rclash, 0);
    EXPECT_FLOAT_EQ(fx.atoms[1].coor[0], x0);
    EXPECT_FLOAT_EQ(fx.FA.ori[0], ori0);
    EXPECT_GE(ic2cf_test_hooks::g_vcfunction_calls, 1);
}

TEST(Ic2cfContamination, VcfunctionErrorRestoresAfterBuildccScribble) {
    ic2cf_test_hooks::g_buildcc_fail_next = 0;
    ic2cf_test_hooks::g_vcfunction_error_next = 0;
    ic2cf_test_hooks::g_buildcc_calls = 0;
    ic2cf_test_hooks::g_vcfunction_calls = 0;

    MinimalIc2cfFixture fx;
    const float x0 = fx.atoms[1].coor[0];
    const float y0 = fx.atoms[1].coor[1];
    const float z0 = fx.atoms[1].coor[2];
    const float ori0 = fx.FA.ori[0];
    fx.residue[1].rot = 7;

    // buildcc succeeds but scribbles coords (stay inside OOB margin so
    // vcfunction is reached); vcfunction errors → restore baseline.
    ic2cf_test_hooks::g_buildcc_scribble_coords = true;
    ic2cf_test_hooks::g_scribble_value = 1.25f;  // ≠ baseline, in-bounds
    ic2cf_test_hooks::g_vcfunction_error_next = 1;
    cfstr bad = ic2cf(&fx.FA, &fx.VC, fx.atoms, fx.residue, nullptr, 0, fx.icv);
    EXPECT_EQ(bad.rclash, 1);
    EXPECT_FLOAT_EQ(fx.atoms[1].coor[0], x0);
    EXPECT_FLOAT_EQ(fx.atoms[1].coor[1], y0);
    EXPECT_FLOAT_EQ(fx.atoms[1].coor[2], z0);
    EXPECT_FLOAT_EQ(fx.FA.ori[0], ori0);
    EXPECT_EQ(fx.residue[1].rot, 7);

    // Second evaluation: clean success, no leftover scribble.
    ic2cf_test_hooks::g_buildcc_scribble_coords = false;
    ic2cf_test_hooks::g_vcfunction_error_next = 0;
    cfstr ok = ic2cf(&fx.FA, &fx.VC, fx.atoms, fx.residue, nullptr, 0, fx.icv);
    EXPECT_EQ(ok.rclash, 0);
    EXPECT_FLOAT_EQ(fx.atoms[1].coor[0], x0);
}
