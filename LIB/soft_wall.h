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
                                           double k_wal_override = 0.0)
{
	// Below this separation d^-12 stops being representable in double with the
	// KWALL_D prefactor. NUMERICAL ONLY -- see the note above.
	constexpr double WALL_D_FLOOR = 0.10;   // A

	const double o_soft = (soft_wall_cutoff > 0.0f)
	                          ? static_cast<double>(soft_wall_cutoff) : 0.0;
	const double k_wal  = (k_wal_override > 0.0) ? k_wal_override : WAL_CONTACT_CAP;

	if (d >= cr) return 0.0;               // no overlap, no wall
	const double o = cr - d;

	if (o_soft > 0.0 && o <= o_soft) {     // unchanged near-contact smoothstep
		const double t = o / o_soft;
		return k_wal * o_soft * o_soft * t * t * (3.0 - 2.0 * t);
	}

	const double d_eff = (d < WALL_D_FLOOR) ? WALL_D_FLOOR : d;
	const double d_s   = cr - o_soft;      // transition separation
	const double base  = k_wal * o_soft * o_soft;
	return base + (wall_energy_raw_r12(d_eff, cr) - wall_energy_raw_r12(d_s, cr));
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
