// tests/test_soft_wall_c1.cpp -- the C1-matched physical wall (c1_match=true).
//
// WHAT THIS GUARDS, AND WHY EACH ASSERTION EXISTS. Every number below is
// COMPUTED from soft_wall.h rather than hardcoded. An earlier draft of this
// test carried a literal 0.4056 taken from a finite-difference slope evaluated
// 1e-4 past the join; the analytic slope gives 0.40785 at the same point, a
// 0.5% gap of pure provenance. A hardcoded expectation would have gone red on
// first run and invited widening the tolerance until it passed.
//
// GATE NON-VACUITY IS ASSERTED FIRST. With c1_match defaulting to false, a test
// that forgets to enable it exercises the legacy branch, reports green, and
// proves nothing about the form it was written for -- the same shape as a unit
// test that links a stub instead of the code under test, which this project has
// hit twice. So the first assertion is that the two branches DIFFER.
//
// THE cr SWEEP IS THE POINT. A single-cr test at 3.50 passes every assertion
// here while O-O (cr 3.10) carries an ATTRACTIVE WELL, because the admissibility
// condition s <= 3*k_wal*o_soft is violated below cr 3.1827 and the unclamped
// cubic dips to -5.05 at cr 2.80. Parameter-space coverage and code-path
// coverage fail identically and neither shows up as a red test.
// TWO KINDS OF ASSERTION HERE, AND BOTH ARE NECESSARY. DO NOT PRUNE EITHER.
//
//   ORACLE tests (6, 6b, 6c) re-derive s, c2, c3 and the admissibility
//   threshold from KWALL_D with their own pow(), so they can disagree with the
//   implementation. They catch an IMPLEMENTATION error: wrong exponent, wrong
//   factor, wrong sign in wall_r12_overlap_slope.
//
//   PROPERTY tests (2, 4, 5, 7) assert non-negativity, monotonicity, the
//   never-softer-than-legacy_uncapped relation and divergence. These are true
//   of ANY correct wall regardless of whether our formula is the right one, so
//   they are the only assertions that can catch a SPECIFICATION error -- us
//   being wrong about the physics rather than wrong about the code. If the
//   derivation of 12*KWALL_D/d_s^13 were itself wrong, the oracle and the
//   implementation would be wrong together and agree perfectly.
//
// Deleting the property tests because "the oracle already pins the values"
// removes the only coverage of the case where we are the ones who are wrong.

#include <gtest/gtest.h>
#include "soft_wall.h"
#include <cmath>

namespace {
constexpr float  OS_F = 0.40f;
const     double OS   = static_cast<double>(OS_F);   // the FLOAT join, not 0.40
constexpr double K    = WAL_CONTACT_CAP;             // reused as stiffness here

double legacy(double d, double cr) { return wall_energy_fitness_physical(d, cr, OS_F, 0.0, false); }
double c1    (double d, double cr) { return wall_energy_fitness_physical(d, cr, OS_F, 0.0, true ); }
}

// ── 0. NON-VACUITY: am I actually on the new branch? ──────────────────
TEST(SoftWallC1, GateIsNotVacuous)
{
	const double cr = 3.50, o = 0.10, d = cr - o;
	const double a = legacy(d, cr), b = c1(d, cr);
	ASSERT_GT(a, 0.0) << "legacy branch returned no wall at all";
	EXPECT_LT(b, a * 0.9) << "c1_match produced the legacy value -- gate not engaged";
}

// ── 1. THE DEFAULT MUST NOT MOVE ──────────────────────────────────────
TEST(SoftWallC1, DefaultReproducesTheShippedForm)
{
	for (double cr : {2.80, 3.10, 3.50, 4.20}) {
		for (double o = 0.01; o <= 2.0; o += 0.01) {
			const double d = cr - o;
			double expect;
			if (o <= OS) {
				const double t = o / OS;
				expect = K * OS * OS * t * t * (3.0 - 2.0 * t);
			} else {
				const double de = (d < 0.10) ? 0.10 : d;
				expect = K * OS * OS
				       + (wall_energy_raw_r12(de, cr) - wall_energy_raw_r12(cr - OS, cr));
			}
			ASSERT_NEAR(legacy(d, cr), expect, 1e-9 * std::fabs(expect) + 1e-12)
			    << "default branch moved at cr=" << cr << " o=" << o;
		}
	}
}

// ── 2. NON-NEGATIVE AND MONOTONE ACROSS THE WHOLE RADIUS TABLE ────────
TEST(SoftWallC1, NonNegativeAndMonotoneOverCrSweep)
{
	for (double cr = 2.40; cr <= 4.80001; cr += 0.05) {
		double prev = -1.0;
		for (double o = 0.0; o <= 2.5; o += 0.005) {
			const double v = c1(cr - o, cr);
			ASSERT_GE(v, -1e-9) << "ATTRACTIVE WELL at cr=" << cr << " o=" << o
			                    << " -> " << v;
			ASSERT_GE(v, prev - 1e-9) << "non-monotone at cr=" << cr << " o=" << o;
			prev = v;
		}
	}
}

// ── 3. C0 AT THE JOIN, and the join value is V ────────────────────────
TEST(SoftWallC1, ContinuousAtTheJoin)
{
	for (double cr : {2.80, 3.10, 3.1827, 3.50, 4.20}) {
		const double V = K * OS * OS;
		EXPECT_NEAR(c1(cr - (OS - 1e-9), cr), V, 1e-6) << "cr=" << cr;
		EXPECT_NEAR(c1(cr - (OS + 1e-9), cr), V, 1e-6) << "cr=" << cr;
	}
}

// ── 4. THE ALGEBRAIC DEFICIT BOUND: h - smoothstep = -S_t*t^2*(1-t) ───
// Never exceeds the shipped smoothstep in the soft band, and the deficit
// cannot exceed 4V/9 for any pair at any cr (since S_t <= 3V).
TEST(SoftWallC1, SoftBandNeverExceedsShippedFormAndDeficitIsBounded)
{
	const double V = K * OS * OS, bound = 4.0 * V / 9.0;
	for (double cr = 2.40; cr <= 4.80001; cr += 0.05) {
		for (double o = 0.0; o <= OS; o += OS / 400.0) {
			const double a = legacy(cr - o, cr), b = c1(cr - o, cr);
			ASSERT_LE(b, a + 1e-9) << "soft branch EXCEEDS shipped form, cr=" << cr;
			ASSERT_LE(a - b, bound + 1e-9) << "deficit exceeds 4V/9, cr=" << cr;
		}
	}
}

// ── 5. ABOVE THE JOIN, NEVER SOFTER THAN THE UNCAPPED LEGACY QUADRATIC ─
TEST(SoftWallC1, AboveJoinNeverSofterThanLegacyUncapped)
{
	for (double cr = 2.40; cr <= 4.80001; cr += 0.05)
		for (double o = OS + 1e-3; o <= 2.5; o += 0.005)
			ASSERT_GE(c1(cr - o, cr),
			          soft_wall_fitness_energy(cr - o, cr, OS_F, true) - 1e-6)
			    << "softer than legacy_uncapped at cr=" << cr << " o=" << o;
}

// ── 6. THE CLAMP BOUNDARY, AGAINST AN INDEPENDENT ORACLE ───────────
//
// An earlier draft computed the expected slope with wall_r12_overlap_slope --
// the SAME helper the implementation calls -- and then asserted the predicate
// agreed with it. That is h == h in a costume: it passes for any
// implementation, including a wrong one. An expectation is only an oracle if it
// COULD disagree with the code.
//
// So this re-derives the threshold from KWALL_D alone, with its own pow():
//   admissible  <=>  12*KWALL_D/d_s^13 <= k_wal*o_soft
//               <=>  d_s >= (12*KWALL_D/(k_wal*o_soft))^(1/13)
// At k_wal = 50, o_soft = 0.40 that is d_s >= 2.7827, i.e. cr >= 3.1827 A --
// a value derived analytically on 2026-09-12 and confirmed by a 1e-4 sweep at
// 3.1826, independently of this code.
TEST(SoftWallC1, ClampBoundaryMatchesAnIndependentOracle)
{
	constexpr double KWALL_D = 1.0e6;          // re-declared, not imported
	const double d_s_min  = std::pow(12.0 * KWALL_D / (K * OS), 1.0 / 13.0);
	const double cr_thresh = d_s_min + OS;

	// the oracle must agree with the value derived by hand, or one of them is wrong
	EXPECT_NEAR(cr_thresh, 3.1827, 5e-4) << "oracle drifted from the derived threshold";

	for (double cr = 2.40; cr <= 4.80001; cr += 0.01) {
		const bool expect_admissible = (cr >= cr_thresh);
		EXPECT_EQ(wall_c1_slope_is_admissible(cr, OS, K), expect_admissible)
		    << "predicate disagrees with the independent oracle at cr=" << cr;
	}
}

// ── 6b. AT THE CLAMP, THE FORM IS THE PURE CUBIC V*t^3 ────────────
// Because the clamp sets s_eff = 3*k_wal*o_soft exactly, c2 = 3V - 3V = 0 and
// c3 = 3V - 2V = V. Asserted against V*t^3 computed here, not against the code.
TEST(SoftWallC1, ClampedFormIsThePureCubic)
{
	constexpr double KWALL_D = 1.0e6;
	const double cr_thresh = std::pow(12.0 * KWALL_D / (K * OS), 1.0 / 13.0) + OS;
	const double V = K * OS * OS;
	int checked = 0;
	for (double cr = 2.60; cr < cr_thresh - 1e-3; cr += 0.02) {
		for (double t = 0.1; t <= 0.9; t += 0.1)
			ASSERT_NEAR(c1(cr - OS * t, cr), V * t * t * t, 1e-9 * V + 1e-12)
			    << "clamped form is not V*t^3 at cr=" << cr << " t=" << t;
		++checked;
	}
	ASSERT_GT(checked, 10) << "no clamped radii exercised -- the sweep is vacuous";
}

// ── 6c. THE UNCLAMPED FORM AGAINST THE ANALYTIC HERMITE ───────────
// s, c2, c3 re-derived here from KWALL_D with an independent pow(), so this can
// disagree with the implementation.
TEST(SoftWallC1, UnclampedFormMatchesTheAnalyticHermite)
{
	constexpr double KWALL_D = 1.0e6;
	const double V = K * OS * OS;
	int checked = 0;
	for (double cr : {3.30, 3.50, 3.80, 4.20, 4.60}) {
		const double d_s = cr - OS;
		double p = 1.0; for (int i = 0; i < 13; ++i) p *= d_s;
		const double s   = 2.0 * K * OS + 12.0 * KWALL_D / p;
		if (s > 3.0 * K * OS) continue;              // clamped, covered by 6b
		const double S_t = s * OS, c2 = 3.0 * V - S_t, c3 = S_t - 2.0 * V;
		ASSERT_GE(c2, 0.0) << "oracle says inadmissible but the clamp did not fire";
		for (double t = 0.05; t <= 1.0; t += 0.05) {
			ASSERT_NEAR(c1(cr - OS * t, cr), c2 * t * t + c3 * t * t * t,
			            1e-9 * V + 1e-12) << "cr=" << cr << " t=" << t;
			// and the closed-form deficit: h - smoothstep == -S_t*t^2*(1-t)
			ASSERT_NEAR(c1(cr - OS * t, cr) - legacy(cr - OS * t, cr),
			            -S_t * t * t * (1.0 - t), 1e-9 * V + 1e-12)
			    << "closed-form deficit violated at cr=" << cr << " t=" << t;
		}
		++checked;
	}
	ASSERT_GT(checked, 2) << "no unclamped radii exercised -- the sweep is vacuous";
}

// ── 7. THE DIVERGENCE SURVIVES ────────────────────────────────────────
TEST(SoftWallC1, DivergesAtDeepInterpenetration)
{
	const double cr = 3.50;
	EXPECT_GT(c1(cr - 2.0, cr), 1.0e3);
	EXPECT_GT(c1(cr - 3.0, cr), c1(cr - 2.0, cr) * 10.0);
}
