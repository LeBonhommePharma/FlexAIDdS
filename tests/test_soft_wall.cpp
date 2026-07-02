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
	// ~200 clashing contacts is realistic for a near-native pose; soft-core
	// keeps the Vcontacts pre-filter sum under CLASH_THRESHOLD.
	const double pose_clash_tally = 200.0 * e_soft;

	EXPECT_LT(e_soft, CLASH_THRESHOLD);
	EXPECT_LT(pose_clash_tally, CLASH_THRESHOLD);
	EXPECT_LE(e_soft, WAL_CONTACT_CAP);
}

TEST(SoftWall, DeepOverlapAvoidsR12Explosion) {
	const double cr = 3.0;
	const double d  = 0.5;
	const float cutoff = 0.40f;

	const double e_soft = soft_wall_fitness_energy(d, cr, cutoff);
	const double e_hard = KWALL * (std::pow(d, -12.0) - std::pow(cr, -12.0));

	EXPECT_LT(e_soft, e_hard);
	EXPECT_NEAR(e_soft, WAL_CONTACT_CAP, EPSILON);
}

TEST(SoftWall, PerContactCap) {
	const double cr = 2.0;
	const double d  = 0.5;   // deep overlap
	const float cutoff = 0.40f;

	const double e = soft_wall_fitness_energy(d, cr, cutoff);
	EXPECT_NEAR(e, WAL_CONTACT_CAP, EPSILON);
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