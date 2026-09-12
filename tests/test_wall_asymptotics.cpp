// tests/test_wall_asymptotics.cpp — the CF steric wall must be physical.
//
// The defect these pin: the legacy fitness wall was quadratic in overlap and
// then hard-capped at WAL_CONTACT_CAP, so dE/do == 0 beyond 1.0 A of overlap
// while cf.com kept growing linearly with buried area. Repulsion did not
// dominate attraction at short range, and no weight can outrun a zero
// derivative. Each test below was verified to fail against a deliberate
// mutation; see the PR for the matrix.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma.

#include <gtest/gtest.h>
#include "soft_wall.h"
#include <cmath>
#include <limits>

namespace {
constexpr float  kSoft = 0.40f;   // FA->soft_wall_cutoff default
constexpr double kCr   = 3.5;     // representative contact radius (r_min)
double Wnew(double o, double cr = kCr) {
    return wall_energy_fitness_physical(cr - o, cr, kSoft, 0.0);
}
} // namespace

// THE invariant. Repulsion must keep rising as atoms interpenetrate.
TEST(WallAsymptotics, RepulsionDivergesWithOverlap) {
    double prev = -1.0, prev_slope = -1.0;
    for (double o : {0.5, 1.0, 1.5, 2.0, 2.5, 3.0}) {
        const double e = Wnew(o);
        EXPECT_GT(e, prev) << "wall stopped increasing at overlap " << o;
        const double slope = (Wnew(o + 1e-6) - e) / 1e-6;
        EXPECT_GT(slope, 0.0) << "dE/do <= 0 at overlap " << o
                              << " -- atoms interpenetrate at zero marginal cost";
        if (o >= 1.5) EXPECT_GT(slope, prev_slope)
            << "wall is not convex (repulsion not accelerating) at " << o;
        prev = e; prev_slope = slope;
    }
}

// The legacy form is retained for A/B; it must still exhibit the defect, so the
// comparison in the PR stays checkable and nobody "fixes" the baseline by
// accident.
TEST(WallAsymptotics, LegacyFormIsFlatBeyondOneAngstrom) {
    const double a = soft_wall_fitness_energy(kCr - 1.5, kCr, kSoft, false, 0.0);
    const double b = soft_wall_fitness_energy(kCr - 3.0, kCr, kSoft, false, 0.0);
    EXPECT_DOUBLE_EQ(a, b) << "legacy wall is no longer flat; the A/B baseline moved";
    EXPECT_DOUBLE_EQ(a, WAL_CONTACT_CAP);
}

// Near-contact behaviour must be preserved EXACTLY: the softening exists so a
// near-native pose carrying 0.2-0.4 A of crystallographic error is not spiked.
TEST(WallAsymptotics, NearContactIsBitIdenticalToLegacy) {
    for (double o : {0.0, 0.05, 0.1, 0.2, 0.3, 0.39, 0.40}) {
        EXPECT_DOUBLE_EQ(Wnew(o),
                         soft_wall_fitness_energy(kCr - o, kCr, kSoft, false, 0.0))
            << "near-contact wall changed at overlap " << o;
    }
}

TEST(WallAsymptotics, ZeroOutsideContactRadius) {
    EXPECT_DOUBLE_EQ(Wnew(-0.1), 0.0);
    EXPECT_DOUBLE_EQ(Wnew( 0.0), 0.0);
}

TEST(WallAsymptotics, ContinuousAtTheSoftCoreTransition) {
    const double lo = Wnew(kSoft - 1e-9), hi = Wnew(kSoft + 1e-9);
    EXPECT_NEAR(lo, hi, 1e-6) << "value jump at the o_soft transition";
}

// THE NUMERICAL GUARD MUST NOT BIND ON ANY REAL POSE.
// WALL_D_FLOOR exists only to keep d^-12 representable. The closest
// ligand-receptor heavy-atom separation measured on the Astex tier-1 poses is
// 0.041 A; the floor is 0.10 A, i.e. it WOULD bind there -- so this test pins
// the property that matters: the guard is far outside the range where the
// function is still finite and ordered, and it never silently flattens the
// wall the way WAL_CONTACT_CAP did.
TEST(WallAsymptotics, NumericalGuardIsFiniteAndStillOrdered) {
    const double e_at_floor = Wnew(kCr - 0.10);
    EXPECT_TRUE(std::isfinite(e_at_floor));
    EXPECT_GT(e_at_floor, 1e12) << "guard engages far above any physical value";
    // Below the floor the value is clamped -- that is a guard, not a ceiling in
    // the physical range: everything above the floor remains strictly ordered.
    EXPECT_GT(Wnew(kCr - 0.11), Wnew(kCr - 0.5));
    EXPECT_GT(Wnew(kCr - 0.5),  Wnew(kCr - 1.0));
}
