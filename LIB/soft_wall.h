// soft_wall.h — Shared overlap-based soft-core clash potential (v43)
//
// Used by Vcontacts pre-filter (get_contlist4 clash_value) and vcfunction
// fitness WAL accumulation so both paths agree on clash tallies.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <cstdlib>
#include <string_view>

// WAL_CONTACT_CAP applies ONLY to the legacy capped r^-12 branch
// (soft_wall_cutoff == 0.0).  It is NOT applied to the soft-core overlap
// branch below — capping that branch flattens the wall past ~1 A of
// overlap, zeroing the fitness gradient and letting the GA bury the ligand
// for unbounded attractive (CF.com) gain.  See tests/test_soft_wall.cpp
// for the regression this guards against.
constexpr double WAL_CONTACT_CAP = 50.0;

// Soft-core wall stiffness (k_wal).  Decoupled from WAL_CONTACT_CAP: this is
// a curvature parameter (energy grows as k_wal * overlap^2), not a ceiling.
// Overridable via FLEXAID_KWAL for calibration sweeps; defaults to 50.0
// (numerically the same magnitude as the legacy cap, but semantically
// unrelated — the soft-core branch is uncapped).
constexpr double K_WAL_STIFF_DEFAULT = 50.0;

inline double k_wal_stiff()
{
	static const double k = []() -> double {
		const char* env = std::getenv("FLEXAID_KWAL");
		if (env) {
			char* end = nullptr;
			const double v = std::strtod(env, &end);
			if (end != env && v > 0.0) return v;
		}
		return K_WAL_STIFF_DEFAULT;
	}();
	return k;
}

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

// Fitness wall energy for clash tally / CF.wal accumulation.
// soft_wall_cutoff = 0.0 recovers legacy capped r^-12 (per-contact ceiling).
// soft_wall_cutoff > 0 applies the v43 overlap Hermite cubic ramp:
//   o <= o_soft: E = k_wal * o_soft^2 * t^2 * (3 - 2t),  t = o/o_soft
//   o >  o_soft: E = k_wal * o_soft^2 + k_wal * (2*o_soft*delta + delta^2)
// Soft-core branch is intentionally UNCAPPED (see k_wal_stiff / tests).
inline double soft_wall_fitness_energy(double d, double cr, float soft_wall_cutoff)
{
	if (soft_wall_cutoff > 0.0f) {
		const double o = cr - d;
		if (o <= 0.0) return 0.0;
		const double o_soft = static_cast<double>(soft_wall_cutoff);
		const double k_wal  = k_wal_stiff();
		if (o <= o_soft) {
			const double t = o / o_soft;
			return k_wal * o_soft * o_soft * t * t * (3.0 - 2.0 * t);
		}
		// Deep region: base + k_wal*(2*o_soft*delta + delta^2) is the
		// expanded form of k_wal*o^2 (delta = o - o_soft), kept in this
		// shifted form for numerical continuity with the ramp at o_soft.
		// UNCAPPED — the prior min(E, WAL_CONTACT_CAP) here flattened the
		// wall past ~1 A overlap and killed the GA's gradient away from
		// buried poses.
		const double base  = k_wal * o_soft * o_soft;
		const double delta = o - o_soft;
		return base + k_wal * (2.0 * o_soft * delta + delta * delta);
	}

	const double Ewall_raw = wall_energy_raw_r12(d, cr);
	return (Ewall_raw > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall_raw;
}
