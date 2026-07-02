// soft_wall.h — Shared overlap-based soft-core clash potential (v43)
//
// Used by Vcontacts pre-filter (get_contlist4 clash_value) and vcfunction
// fitness WAL accumulation so both paths agree on clash tallies.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>

constexpr double WAL_CONTACT_CAP = 50.0;

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
inline double soft_wall_fitness_energy(double d, double cr, float soft_wall_cutoff)
{
	if (soft_wall_cutoff > 0.0f) {
		const double o      = cr - d;
		const double o_soft = static_cast<double>(soft_wall_cutoff);
		constexpr double k_wal = WAL_CONTACT_CAP;
		double Ewall_sc;
		if (o <= o_soft) {
			const double t = o / o_soft;
			Ewall_sc = k_wal * o_soft * o_soft * t * t * (3.0 - 2.0 * t);
		} else {
			const double base  = k_wal * o_soft * o_soft;
			const double delta = o - o_soft;
			Ewall_sc = base + k_wal * (2.0 * o_soft * delta + delta * delta);
		}
		return (Ewall_sc > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall_sc;
	}

	const double Ewall_raw = wall_energy_raw_r12(d, cr);
	return (Ewall_raw > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall_raw;
}