// tests/test_soft_wall.cpp — Unit tests for overlap-based soft-core clash potential
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/soft_wall.h"
#include "../LIB/flexaid.h"

#include <cmath>

static constexpr double EPSILON = 1e-6;

TEST(SoftWall, ModerateOverlapBelowClashThreshold) {
	// Typical crystal-coordinate micro-overlap: o = 0.30 Å inside cr = 3.0 Å
	const double cr = 3.0;
	const double d  = cr - 0.30;
	const float cutoff = 0.40f;

	const double e_soft = soft_wall_fitness_energy(d, cr, cutoff);
	// ~200 clashing contacts is realistic for a near-native pose; a shallow
	// ramp-region overlap must still keep the Vcontacts pre-filter sum
	// under CLASH_THRESHOLD.
	const double pose_clash_tally = 200.0 * e_soft;

	EXPECT_LT(e_soft, CLASH_THRESHOLD);
	EXPECT_LT(pose_clash_tally, CLASH_THRESHOLD);
}

TEST(SoftWall, DeepOverlapIsUncappedQuadratic) {
	// Regression guard for the flat-cap bug: min(E, WAL_CONTACT_CAP) on the
	// soft-core branch flattened the wall past ~1 A overlap, zeroing the
	// GA's gradient away from buried poses. The deep region must now grow
	// as k_wal * o^2 with NO ceiling, while still staying below the raw
	// (uncapped) r^-12 potential at the same depth.
	const double cr = 3.0;
	const double d  = 0.5;
	const float cutoff = 0.40f;

	const double e_soft = soft_wall_fitness_energy(d, cr, cutoff);
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

TEST(SoftWall, MonotonicWithOverlapDepth) {
	const double cr = 3.0;
	const float cutoff = 0.40f;
	double prev = -1.0;
	for (double d = 3.0; d >= 0.2; d -= 0.1) {
		const double e = soft_wall_fitness_energy(d, cr, cutoff);
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

TEST(SoftWall, HermiteRampIsContinuousAtCutoff) {
	const double cr = 3.0;
	const float cutoff = 0.40f;
	const double d_at  = cr - static_cast<double>(cutoff);
	const double d_eps = d_at - 1e-6;

	const double e_at  = soft_wall_fitness_energy(d_at,  cr, cutoff);
	const double e_eps = soft_wall_fitness_energy(d_eps, cr, cutoff);

	EXPECT_NEAR(e_at, e_eps, 0.05);
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
