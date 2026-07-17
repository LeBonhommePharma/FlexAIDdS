// tests/test_soft_wall.cpp — Unit tests for overlap-based soft-core clash potential
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/soft_wall.h"
#include "../LIB/flexaid.h"

#include <cmath>
#include <cstdlib>
#include <string>

static constexpr double EPSILON = 1e-6;

// Finite-difference slope dE/do ≈ (E(o+h)-E(o-h))/(2h) with o = cr-d.
static double slope_at_overlap(double o, double cr, float cutoff, double k_wal, double h = 1e-5)
{
	const double d0 = cr - (o - h);
	const double d1 = cr - (o + h);
	const double e0 = soft_wall_fitness_energy(d0, cr, cutoff, k_wal);
	const double e1 = soft_wall_fitness_energy(d1, cr, cutoff, k_wal);
	return (e1 - e0) / (2.0 * h);
}

TEST(SoftWall, ModerateOverlapBelowClashThreshold) {
	// Typical crystal-coordinate micro-overlap: o = 0.30 Å inside cr = 3.0 Å
	const double cr = 3.0;
	const double d  = cr - 0.30;
	const float cutoff = 0.40f;

	const double e_soft = soft_wall_fitness_energy(d, cr, cutoff, K_WAL_STIFF_DEFAULT);
	// Pure quadratic: E = 50 * 0.30² = 4.5
	EXPECT_NEAR(e_soft, K_WAL_STIFF_DEFAULT * 0.30 * 0.30, 1e-9);
	// ~200 clashing contacts must stay under CLASH_THRESHOLD.
	const double pose_clash_tally = 200.0 * e_soft;
	EXPECT_LT(e_soft, CLASH_THRESHOLD);
	EXPECT_LT(pose_clash_tally, CLASH_THRESHOLD);
}

TEST(SoftWall, DeepOverlapIsUncappedQuadratic) {
	// Regression guard: soft-core must grow as k_wal * o^2 with NO ceiling.
	const double cr = 3.0;
	const double d  = 0.5;
	const float cutoff = 0.40f;

	const double e_soft = soft_wall_fitness_energy(d, cr, cutoff, K_WAL_STIFF_DEFAULT);
	const double e_hard = KWALL * (std::pow(d, -12.0) - std::pow(cr, -12.0));
	const double o = cr - d;
	const double e_expected = K_WAL_STIFF_DEFAULT * o * o;

	EXPECT_GT(e_soft, WAL_CONTACT_CAP);       // must exceed the legacy cap
	EXPECT_LT(e_soft, e_hard);                // but stay below raw r^-12
	EXPECT_NEAR(e_soft, e_expected, 1e-6);    // pure quadratic, uncapped
}

TEST(SoftWall, ZeroAtNoOverlap) {
	const double cr = 2.0;
	const double d  = 2.5;   // no overlap: d > cr
	const float cutoff = 0.40f;
	EXPECT_NEAR(soft_wall_fitness_energy(d, cr, cutoff), 0.0, EPSILON);
}

TEST(SoftWall, PureQuadraticIsC1Everywhere) {
	// Soft-core is E = k o² for all o>0 → C∞; slope at any o is 2 k o.
	const double cr = 3.0;
	const float cutoff = 0.40f;
	const double k = K_WAL_STIFF_DEFAULT;
	for (double o : {0.05, 0.20, 0.40, 0.80, 1.50, 2.50}) {
		const double slope = slope_at_overlap(o, cr, cutoff, k);
		EXPECT_NEAR(slope, 2.0 * k * o, 1e-3) << "o=" << o;
	}
	// At o→0+, slope → 0
	const double slope0 = slope_at_overlap(1e-4, cr, cutoff, k, 1e-5);
	EXPECT_NEAR(slope0, 0.0, 0.05);
}

TEST(SoftWall, MonotonicWithOverlapDepth) {
	const double cr = 3.0;
	const float cutoff = 0.40f;
	double prev = -1.0;
	for (double d = 3.0; d >= 0.2; d -= 0.1) {
		const double e = soft_wall_fitness_energy(d, cr, cutoff, K_WAL_STIFF_DEFAULT);
		EXPECT_GE(e, prev);
		prev = e;
	}
}

TEST(SoftWall, LegacyPathMatchesCappedR12) {
	const double cr = 3.0;
	const double d  = 2.5;
	const double raw = wall_energy_raw_r12(d, cr);
	const double capped = (raw > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : raw;
	EXPECT_NEAR(soft_wall_fitness_energy(d, cr, 0.0f), capped, 1e-3);
}

TEST(SoftWall, ContactAtZeroDistanceIsFinite) {
	// d=0 soft-core: o=cr → E = k cr² (finite). Legacy path caps.
	const double cr = 3.0;
	const double e_soft = soft_wall_fitness_energy(0.0, cr, 0.40f, K_WAL_STIFF_DEFAULT);
	EXPECT_TRUE(std::isfinite(e_soft));
	EXPECT_NEAR(e_soft, K_WAL_STIFF_DEFAULT * cr * cr, 1e-6);

	const double e_leg = soft_wall_fitness_energy(0.01, cr, 0.0f);
	EXPECT_TRUE(std::isfinite(e_leg));
	EXPECT_LE(e_leg, WAL_CONTACT_CAP + EPSILON);
}

TEST(SoftWall, HugeOverlapUncappedQuadratic) {
	const double cr = 3.0;
	const double d  = 0.01;  // o ≈ 2.99
	const double o  = cr - d;
	const double e  = soft_wall_fitness_energy(d, cr, 0.40f, K_WAL_STIFF_DEFAULT);
	EXPECT_NEAR(e, K_WAL_STIFF_DEFAULT * o * o, 1e-4);
	EXPECT_GT(e, WAL_CONTACT_CAP * 5.0);
}

TEST(SoftWall, ExplicitKWalScalesEnergy) {
	const double cr = 3.0;
	const double d  = cr - 1.0;
	const double e50  = soft_wall_fitness_energy(d, cr, 0.40f, 50.0);
	const double e100 = soft_wall_fitness_energy(d, cr, 0.40f, 100.0);
	EXPECT_NEAR(e100, 2.0 * e50, 1e-9);
	EXPECT_NEAR(e50, 50.0 * 1.0 * 1.0, 1e-9);
}

TEST(SoftWall, MultiContactDeepClashExceedsLegacyCapSum) {
	// N contacts at o=1.5: old flat cap → N*50; uncap → N*k*o²
	const double cr = 3.0;
	const double d  = cr - 1.5;
	const double e  = soft_wall_fitness_energy(d, cr, 0.40f, 50.0);
	const int N = 20;
	const double sum_uncap = N * e;
	const double sum_old_cap = N * WAL_CONTACT_CAP;
	EXPECT_NEAR(e, 50.0 * 1.5 * 1.5, 1e-9);  // 112.5
	EXPECT_GT(sum_uncap, sum_old_cap);
	EXPECT_LT(sum_uncap, CLASH_THRESHOLD);  // still a finite pre-filter signal
}

TEST(SoftWall, EnvKWalOverrideViaResolve) {
	// resolve_k_wal reads env each call (no process-static cache).
	// Only run if we can setenv safely.
#if defined(_WIN32)
	GTEST_SKIP() << "setenv not used on Windows in this test";
#else
	const char* prev_a = std::getenv("FLEXAIDDS_K_WAL");
	const char* prev_b = std::getenv("FLEXAID_KWAL");
	std::string save_a = prev_a ? prev_a : "";
	std::string save_b = prev_b ? prev_b : "";
	unsetenv("FLEXAIDDS_K_WAL");
	unsetenv("FLEXAID_KWAL");
	EXPECT_NEAR(resolve_k_wal(0.0), K_WAL_STIFF_DEFAULT, 1e-12);
	setenv("FLEXAID_KWAL", "100", 1);
	EXPECT_NEAR(resolve_k_wal(0.0), 100.0, 1e-12);
	// Explicit arg wins over env
	EXPECT_NEAR(resolve_k_wal(25.0), 25.0, 1e-12);
	setenv("FLEXAID_KWAL", "0", 1);  // invalid → default
	EXPECT_NEAR(resolve_k_wal(0.0), K_WAL_STIFF_DEFAULT, 1e-12);
	// restore
	if (save_a.empty()) unsetenv("FLEXAIDDS_K_WAL"); else setenv("FLEXAIDDS_K_WAL", save_a.c_str(), 1);
	if (save_b.empty()) unsetenv("FLEXAID_KWAL"); else setenv("FLEXAID_KWAL", save_b.c_str(), 1);
#endif
}

TEST(SoftWall, FloatKernelParityWithDouble) {
	const float cr = 3.0f;
	const float d  = 1.5f;
	const float e_f = soft_wall_fitness_energy_f(d, cr, 0.40f, 50.0f);
	const double e_d = soft_wall_fitness_energy(d, cr, 0.40f, 50.0);
	EXPECT_NEAR(static_cast<double>(e_f), e_d, 1e-4);
}

TEST(SoftWall, RelativeVdwCutoffMatchesPoseBustersBoundary) {
	const double oxygen_radius = posebusters_vdw_radius("O", 0.0);
	const double radius_sum = 2.0 * oxygen_radius;
	const float cutoff = 0.75f;
	const double boundary = static_cast<double>(cutoff) * radius_sum;

	EXPECT_NEAR(oxygen_radius, 1.55, EPSILON);
	EXPECT_NEAR(posebusters_vdw_radius("Cl", 0.0), 1.80, EPSILON);
	EXPECT_NEAR(posebusters_vdw_radius("unknown", 1.42), 1.42, EPSILON);
	EXPECT_FALSE(violates_relative_vdw_cutoff(boundary - 1.0, radius_sum, 0.0f));
	EXPECT_TRUE(violates_relative_vdw_cutoff(1.904304, radius_sum, cutoff));
	EXPECT_TRUE(violates_relative_vdw_cutoff(boundary - 1e-6, radius_sum, cutoff));
	EXPECT_FALSE(violates_relative_vdw_cutoff(boundary, radius_sum, cutoff));
	EXPECT_FALSE(violates_relative_vdw_cutoff(boundary + 1e-6, radius_sum, cutoff));
}
