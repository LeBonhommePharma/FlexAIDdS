// soft_wall.h — Shared overlap-based soft-core clash potential (v43)
//
// Used by Vcontacts pre-filter (get_contlist4 clash_value) and vcfunction
// fitness WAL accumulation so both paths agree on clash tallies.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <string_view>

constexpr double WAL_CONTACT_CAP = 50.0;

inline bool violates_relative_vdw_cutoff(double distance,
                                         double radius_sum,
                                         float cutoff_ratio)
{
	return cutoff_ratio > 0.0f && radius_sum > 0.0 &&
	       distance < static_cast<double>(cutoff_ratio) * radius_sum;
}

// RDKit periodic-table vdW radii used by PoseBusters' intermolecular-distance
// check. Keep this separate from FlexAID's NRG contact radii: the two serve
// different purposes and are not numerically interchangeable.
inline double posebusters_vdw_radius(std::string_view element, double fallback)
{
	if (element == "H")  return 1.20;
	if (element == "C")  return 1.70;
	if (element == "N")  return 1.60;
	if (element == "O")  return 1.55;
	if (element == "F")  return 1.50;
	if (element == "P")  return 1.95;
	if (element == "S")  return 1.80;
	if (element == "Cl") return 1.80;
	if (element == "Br") return 1.90;
	if (element == "I")  return 2.10;
	if (element == "Se") return 1.90;
	if (element == "Mg") return 2.20;
	if (element == "Sr") return 2.55;
	if (element == "Cu") return 2.00;
	if (element == "Mn") return 2.05;
	if (element == "Hg") return 2.05;
	if (element == "Cd") return 2.20;
	if (element == "Ni") return 2.00;
	if (element == "Zn") return 2.10;
	if (element == "Ca") return 2.40;
	if (element == "Fe") return 2.05;
	if (element == "Co") return 2.00;
	if (element == "Na") return 2.40;
	if (element == "K")  return 2.80;
	return fallback;
}

// Coordinating-metal carve-out set for the pb_clash intermolecular scan.
//
// A genuine ligand->metal coordination bond sits at ~1.9-2.3 A, which is far
// below pb_clash_ratio*(vdw_lig + vdw_metal) (e.g. 0.75*(1.55+2.10) = 2.74 A for
// an O donor to Zn). Without a carve-out the pb_clash penalty scores every real
// coordination bond as a hard clash and fights the cf.metal_coord Morse reward
// that is enabled in the generated benchmark configs.
//
// The set is the catalytic/structural transition metals and alkaline earths that
// actually appear as coordinating centres in the benchmark receptors. Alkali
// ions (Na, K) are deliberately NOT included: they are almost always spectator
// ions at crystallographic distances that a clash term should still see.
inline bool posebusters_is_coordinating_metal(std::string_view element)
{
	return element == "Zn" || element == "Fe" || element == "Mg" ||
	       element == "Ca" || element == "Mn" || element == "Co" ||
	       element == "Ni" || element == "Cu";
}

// Raw r^-12 wall energy: KWALL * (d^-12 - cr^-12).  Unbounded as d -> 0.
inline double wall_energy_raw_r12(double d, double cr)
{
	constexpr double KWALL_D = 1.0e6;
	const double d2  = d  * d;
	const double d4  = d2 * d2;
	const double d6  = d4 * d2;
	const double inv_d12  = 1.0 / (d6 * d6);
	const double cr2 = cr * cr;
	const double cr4 = cr2 * cr2;
	const double cr6 = cr4 * cr2;
	const double inv_cr12 = 1.0 / (cr6 * cr6);
	return KWALL_D * (inv_d12 - inv_cr12);
}

// Analytic slope of wall_energy_raw_r12 with respect to OVERLAP o = cr - d:
//   d/do [ KWALL_D * ((cr-o)^-12 - cr^-12) ] = 12 * KWALL_D / (cr-o)^13
// Closed form deliberately, NOT a finite difference. The continuation's second
// derivative at the transition is ~120 per A^2, so sampling the slope 1e-3 A
// past the join reports it 0.27% high with no warning -- measured 2026-09-12,
// where two probes disagreed by exactly C_o * offset to six figures.
inline double wall_r12_overlap_slope(double d)
{
	constexpr double KWALL_D = 1.0e6;
	const double d2 = d * d, d4 = d2 * d2, d8 = d4 * d4;
	return 12.0 * KWALL_D / (d8 * d4 * d);      // 12*KWALL_D / d^13
}

// ADMISSIBILITY OF A C1-MATCHED SOFT BRANCH.
//
// A cubic Hermite on [0, o_soft] with h(0)=h'(0)=0, h(o_soft)=V, h'(o_soft)=s
// is, in t = o/o_soft space,
//   h(t) = c2*t^2 + c3*t^3,  c2 = 3V - S_t,  c3 = S_t - 2V,  S_t = s*o_soft
// and h'(t) = t*(2*c2 + 3*c3*t). With c2 >= 0 the bracket decreases in t to a
// minimum 2*c2 + 3*c3 = S_t >= 0, so monotonicity binds on c2 ALONE:
//
//   ADMISSIBLE  <=>  c2 >= 0  <=>  S_t <= 3V  <=>  s <= 3*k_wal*o_soft
//
// c3 < 0 is harmless: the shipped smoothstep runs c3 = -2V and is monotone. A
// LOWER bound on S_t is NOT required, and asserting one is over-tight.
//
// WHY IT BINDS. Under a quadratic-continued upper branch the matched slope is
// s = 2*k_wal*o_soft + 12*KWALL_D/d_s^13, whose second term grows as d_s^-13.
// Measured at k_wal=50, o_soft=0.40: admissible only for cr >= 3.1827 A (closed
// form d_s^13 >= 12*KWALL_D/(k_wal*o_soft) = 6.0e5; swept value 3.1826). Below
// that c2 goes NEGATIVE and the cubic dips below zero -- an ATTRACTIVE well
// inside the repulsion term, min -5.05 at cr=2.80 and -83.30 at cr=2.40. With
// vdW-sum radii that band contains O-O (3.10) and N-O (3.15), i.e. the
// hydrogen-bonding pairs. A single-cr test at 3.50 sees none of it.
//
// RESCUE OPTIONS, ALL MEASURED 2026-09-12:
//   o_soft  DEAD for O-O at ANY width. The condition is I <= k_wal*o_soft with
//           I = 12*KWALL_D/d_s^13. Shrinking o_soft sends the budget to zero
//           while I stays bounded below by 12*KWALL_D/cr^13; growing it blows I
//           up as d_s^-13. So there is an interior optimum, and for O-O it peaks
//           at -1.68 -- k_wal*o_soft = I has NO ROOT for that pair. (N-O does
//           have one, at o_soft 0.152.)
//   k_wal   Works arithmetically: the exact threshold is
//           12*KWALL_D/(o_soft*d_s^13) = 74.027 for O-O, so k_wal >= 75 clears
//           it and 74 does NOT (margin -0.0109). Available today via
//           FLEXAIDDS_WAL_STIFF: k_wal_override and the WAL_CONTACT_CAP ceiling
//           are already separate roles (:139/:172 stiffness, :184/:189 ceiling),
//           so no constant split is needed. NOT a substitute for the clamp --
//           I scales as d_s^-13, so a 0.03 A radius-convention change (O 1.52
//           vs 1.55) swings I by 34% against a 1.3% margin at k_wal=75, and a
//           test built on the same radius table cannot detect the regression.
//   clamp   s_eff = min(s, 3*k_wal*o_soft). THE BOUNDARY OF THE ADMISSIBLE SET,
//           not a fallback: the steepest arrival achievable without an
//           attractive well. At the bound c2 == 0 exactly and h(t) = V*t^3, a
//           pure cubic -- monotone, non-negative, exact at the endpoint.
//           Admissible for every cr, including radii not yet in the table, by
//           construction rather than by arithmetic. THIS IS WHAT SHIPS.
inline bool wall_c1_slope_is_admissible(double cr, double o_soft, double k_wal)
{
	if (!(o_soft > 0.0) || !(cr > o_soft)) return true;
	return wall_r12_overlap_slope(cr - o_soft) <= k_wal * o_soft;
}

// ── PHYSICAL FITNESS WALL (CF scoring path) ─────────────────────────────────
//
// THE DEFECT THIS REPLACES. soft_wall_fitness_energy() below is quadratic in
// overlap and then HARD-CAPPED at WAL_CONTACT_CAP. Measured on this tree with
// the default soft_wall_cutoff = 0.40 and k_wal = 50: the cap binds at exactly
// 1.0 A of overlap, beyond which dE/do == 0 EXACTLY. Two atoms may then
// interpenetrate arbitrarily far at zero marginal cost, while cf.com keeps
// growing linearly with buried Voronoi area (vcfunction.cpp: contribution =
// yval * area, unbounded). Repulsion therefore does NOT dominate attraction at
// short range -- it is constant where attraction is increasing -- so the
// gradient points at maximum burial everywhere in the overlap regime. That is
// non-physical, and no reweighting can repair it: a weight cannot outrun a zero
// derivative. Measured consequence (1gpk): the engine's own docked pose carries
// 21.8x the wall term of the crystal pose and still scores 14.7x better.
//
// THE SHAPE, AND WHY IT IS THIS SHAPE. The softening exists for a real reason,
// recorded at vcfunction.cpp:1081 -- a bare r^-12 spikes just inside cr, so a
// near-native pose carrying 0.2-0.4 A of crystallographic coordinate error can
// score worse than a decoy. That protection is KEPT EXACTLY: for overlap up to
// o_soft this function is the same Hermite smoothstep, bit for bit. Beyond
// o_soft the quadratic-then-capped continuation is replaced by the TRUE r^-12
// increment, so the Pauli wall is recovered where it belongs -- at real
// interpenetration -- and diverges as d -> 0.
//
//   o <= o_soft : k_wal * o_soft^2 * t^2 * (3 - 2t),  t = o/o_soft   (unchanged)
//   o >  o_soft : k_wal * o_soft^2 + [ raw_r12(d) - raw_r12(cr - o_soft) ]
//
// C0-continuous by construction: the bracket vanishes at o == o_soft.
//
// NO FITTED CONSTANT IS INTRODUCED. Every quantity is pre-existing: the contact
// radius cr (= r_min, where the wall is zero), KWALL_D, and o_soft. The shape
// follows from matching the physical form at the transition, not from anything
// tuned to make a test pass.
//
// NUMERICAL GUARD, NOT A PHYSICAL CEILING. d is floored at WALL_D_FLOOR purely
// so d^-12 stays representable; it is deliberately far outside the physically
// meaningful range and must never bind on a real pose. At the floor the energy
// is ~1e18 while the largest per-contact value measured on real Astex poses is
// ~1e3 -- fifteen orders of magnitude of headroom. A test asserts it does not
// bind; if it ever does, that is a bug to investigate, not a value to clamp.
inline double wall_energy_fitness_physical(double d, double cr,
                                           float soft_wall_cutoff,
                                           double k_wal_override = 0.0,
                                           bool c1_match = false)
{
	// Below this separation d^-12 stops being representable in double with the
	// KWALL_D prefactor. NUMERICAL ONLY -- see the note above.
	constexpr double WALL_D_FLOOR = 0.10;   // A

	const double o_soft = (soft_wall_cutoff > 0.0f)
	                          ? static_cast<double>(soft_wall_cutoff) : 0.0;
	const double k_wal  = (k_wal_override > 0.0) ? k_wal_override : WAL_CONTACT_CAP;

	if (d >= cr) return 0.0;               // no overlap, no wall
	const double o = cr - d;

	const double d_eff = (d < WALL_D_FLOOR) ? WALL_D_FLOOR : d;
	const double d_s   = cr - o_soft;      // transition separation
	const double base  = k_wal * o_soft * o_soft;   // V, the join value

	// ── c1_match == false: BYTE-IDENTICAL to the shipped form ───────────
	// Two MEASURED defects live in this branch; both are fixed only under
	// c1_match. Recorded here so the default is understood as legacy, not
	// as correct:
	//
	// (1) DERIVATIVE BREAK AT THE JOIN. The smoothstep has dE/do == 0 at
	//     t == 1 by construction while the r^-12 increment leaves at
	//     12*KWALL_D/d_s^13. Measured cr=3.50, o_soft=0.40, k_wal=50:
	//     0.029994 -> 4.916560, a 164x jump. C0 holds, C1 does not.
	//
	// (2) SOFTER THAN THE FORM IT REPLACED, over a third of the range.
	//     Because the increment is offset to vanish at the join, it starts
	//     from zero and climbs more slowly than the capped quadratic it
	//     superseded. Measured over 2500 samples of o in (0, 2.5]:
	//     softer than legacy_capped on 789, worst deficit -26.4924 against
	//     a per-contact cap of 50. Crossover near o = 1.2 A, so the band
	//     0.41-1.15 A -- where real steric clashes sit -- is CONCEDED.
	if (!c1_match) {
		if (o_soft > 0.0 && o <= o_soft) {     // unchanged near-contact smoothstep
			const double t = o / o_soft;
			return k_wal * o_soft * o_soft * t * t * (3.0 - 2.0 * t);
		}
		return base + (wall_energy_raw_r12(d_eff, cr) - wall_energy_raw_r12(d_s, cr));
	}

	// ── c1_match == true: C1-matched soft branch + quadratic continuation ──
	//
	//   o <= o_soft : c2*t^2 + c3*t^3,  c2 = 3V - S_t,  c3 = S_t - 2V
	//                 S_t = s_eff*o_soft,  s_eff = min(s, 3*k_wal*o_soft)
	//   o >  o_soft : k_wal*o^2 + [ raw_r12(d) - raw_r12(d_s) ]
	//
	// The upper branch CONTINUES the quadratic rather than freezing it at its
	// join value. Since the legacy uncapped form IS k_wal*o^2 beyond the join
	// (verified: identical to 0.000e+00 at o = 0.50/0.75/1.00/1.50), adding a
	// non-negative increment makes this form NEVER SOFTER than what it replaced,
	// above the join, as a theorem rather than a tabulation. Measured 1/2500
	// residual at -2.4e-07, which is the float o_soft (see below), not structure.
	//
	// s_eff is CLAMPED at the admissible boundary -- see
	// wall_c1_slope_is_admissible above for why, and for the measured
	// consequence of not clamping (an attractive well on O-O and N-O).
	//
	// FLOAT JOIN. soft_wall_cutoff is a float, so o_soft widens to
	// 0.40000000596046448, not 0.40. Any test comparing against a double 0.40
	// literal inherits ~1.5e-08 relative, which is exactly the 2.98e-08
	// "violation" a probe reported on 2026-09-12 before the cause was found.
	//
	// SOFT-BAND DEFICIT, in closed form. Subtracting the smoothstep gives
	// h - smoothstep = -S_t*t^2*(1-t), which is <= 0 for every S_t >= 0 -- so
	// "never exceeds the shipped smoothstep" is ALGEBRAIC, not sampled. It peaks
	// at t = 2/3 with depth (4/27)*S_t, and since S_t <= 3V the deficit can
	// never exceed 4V/9, i.e. 44.4% of the join value, for ANY pair at ANY cr.
	// t = 2/3 is also where the two slopes cross (o = 0.2667 at o_soft = 0.40),
	// because the deficit is stationary exactly where the slopes are equal.
	//
	// COST. 12*KWALL_D/d_s^13 depends only on cr and o_soft, both fixed per
	// atom-type pair at setup, so it hoists out of the inner loop. Recomputing
	// it per contact instead measures 1.90x (1.364 -> 2.585 ns/contact) and is
	// the wrong implementation.
	if (o_soft > 0.0 && o <= o_soft) {
		const double s_raw = 2.0 * k_wal * o_soft + wall_r12_overlap_slope(d_s);
		const double s_max = 3.0 * k_wal * o_soft;
		const double s_eff = (s_raw > s_max) ? s_max : s_raw;
		const double S_t   = s_eff * o_soft;
		const double c2    = 3.0 * base - S_t;
		const double c3    = S_t - 2.0 * base;
		const double t     = o / o_soft;
		return c2 * t * t + c3 * t * t * t;
	}

	return k_wal * o * o
	     + (wall_energy_raw_r12(d_eff, cr) - wall_energy_raw_r12(d_s, cr));
}

// Fitness wall energy for clash tally / CF.wal accumulation.
// soft_wall_cutoff = 0.0 recovers legacy capped r^-12 (per-contact ceiling).
// soft_wall_cutoff > 0 applies the v43 overlap Hermite cubic ramp:
//   o <= o_soft: E = k_wal * o_soft^2 * t^2 * (3 - 2t),  t = o/o_soft
//   o >  o_soft: E = k_wal * o_soft^2 + k_wal * (2*o_soft*delta + delta^2)
//
// coercive=true  (FLEXAIDDS_WAL_COERCIVE): removes the WAL_CONTACT_CAP ceiling
//   so the quadratic penalty can overcome unbounded CF.com overpacking.
// k_wal_override > 0 (FLEXAIDDS_WAL_STIFF): replaces the welded k_wal=WAL_CONTACT_CAP
//   with the supplied value, allowing stiffness sweeps at benchmark time.
inline double soft_wall_fitness_energy(double d, double cr, float soft_wall_cutoff,
                                       bool coercive = false,
                                       double k_wal_override = 0.0)
{
	if (soft_wall_cutoff > 0.0f) {
		const double o      = cr - d;
		const double o_soft = static_cast<double>(soft_wall_cutoff);
		const double k_wal  = (k_wal_override > 0.0) ? k_wal_override : WAL_CONTACT_CAP;
		double Ewall_sc;
		if (o <= o_soft) {
			const double t = o / o_soft;
			Ewall_sc = k_wal * o_soft * o_soft * t * t * (3.0 - 2.0 * t);
		} else {
			const double base  = k_wal * o_soft * o_soft;
			const double delta = o - o_soft;
			Ewall_sc = base + k_wal * (2.0 * o_soft * delta + delta * delta);
		}
		// FLEXAIDDS_WAL_COERCIVE: skip the cap so deep clashes overwhelm CF.com.
		// Default (coercive=false): preserve existing WAL_CONTACT_CAP ceiling.
		if (!coercive && Ewall_sc > WAL_CONTACT_CAP) return WAL_CONTACT_CAP;
		return Ewall_sc;
	}

	const double Ewall_raw = wall_energy_raw_r12(d, cr);
	return (Ewall_raw > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall_raw;
}
