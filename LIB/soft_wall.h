// soft_wall.h — Shared overlap-based soft-core clash potential
//
// Used by Vcontacts pre-filter (get_contlist4 clash_value), vcfunction
// fitness WAL accumulation, and CPU/CUDA/Metal batch eval paths so clash
// tallies agree across backends.
//
// Soft-core (soft_wall_cutoff > 0):
//   E = k_wal · max(o, 0)²   with o = cr − d
//   Uncapped pure quadratic — C¹ (actually C∞) for o > 0, E'(0+)=0.
//   Restores GA gradient for deep burials (o ≳ 1 Å) that the old
//   min(E, WAL_CONTACT_CAP) path flattened to a dead zone.
//   soft_wall_cutoff only selects soft-core vs legacy (value > 0 enables).
//
// Legacy (soft_wall_cutoff == 0):
//   min(KWALL · (d⁻¹² − cr⁻¹²), WAL_CONTACT_CAP)
//
// Protocol / claim note: default soft_wall_cutoff=0.40 means production
// CF/scoring-proxy paths use soft-core. Deep-clash CF ranking differs from
// pre-uncap soft-core campaigns (o ≳ 1 Å). Stamp k_wal + cutoff in receipts.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <cstdlib>
#include <string_view>

// WAL_CONTACT_CAP applies ONLY to the legacy capped r^-12 branch
// (soft_wall_cutoff == 0.0). Soft-core is intentionally uncapped.
constexpr double WAL_CONTACT_CAP = 50.0;

// Soft-core wall stiffness (k_wal). Curvature parameter: energy grows as
// k_wal * overlap^2. Overridable via FLEXAID_KWAL / FLEXAIDDS_K_WAL env
// (read each call — no process-static cache) or FA->k_wal_stiff.
constexpr double K_WAL_STIFF_DEFAULT = 50.0;

/// Resolve k_wal from explicit value, else FLEXAID_KWAL / FLEXAIDDS_K_WAL, else default.
inline double resolve_k_wal(double k_explicit = 0.0) noexcept
{
	if (k_explicit > 0.0 && std::isfinite(k_explicit))
		return k_explicit;
	for (const char* name : {"FLEXAIDDS_K_WAL", "FLEXAID_KWAL"}) {
		const char* env = std::getenv(name);
		if (!env || !env[0])
			continue;
		char* end = nullptr;
		const double v = std::strtod(env, &end);
		if (end != env && v > 0.0 && std::isfinite(v))
			return v;
	}
	return K_WAL_STIFF_DEFAULT;
}

/// Deprecated alias: prefer resolve_k_wal / explicit FA->k_wal_stiff.
inline double k_wal_stiff() noexcept
{
	return resolve_k_wal(0.0);
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
//
// soft_wall_cutoff > 0: uncapped soft-core E = k_wal * max(cr-d, 0)^2
//   (C¹ at o=0; pure quadratic — matches deep k·o² value and slope everywhere)
// soft_wall_cutoff == 0: legacy capped r^-12
//
// k_wal_explicit <= 0 → resolve from env / default (see resolve_k_wal).
inline double soft_wall_fitness_energy(double d, double cr, float soft_wall_cutoff,
                                       double k_wal_explicit = 0.0)
{
	if (soft_wall_cutoff > 0.0f) {
		const double o = cr - d;
		if (o <= 0.0) return 0.0;
		const double k_wal = resolve_k_wal(k_wal_explicit);
		return k_wal * o * o;
	}

	const double Ewall_raw = wall_energy_raw_r12(d, cr);
	return (Ewall_raw > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall_raw;
}

/// Float variant for CPU/GPU batch kernels (same physics).
inline float soft_wall_fitness_energy_f(float d, float cr, float soft_wall_cutoff,
                                        float k_wal)
{
	if (soft_wall_cutoff > 0.0f) {
		const float o = cr - d;
		if (o <= 0.0f) return 0.0f;
		const float k = (k_wal > 0.0f) ? k_wal : static_cast<float>(K_WAL_STIFF_DEFAULT);
		return k * o * o;
	}
	// Legacy capped r^-12
	constexpr float KWALL_F = 1.0e6f;
	constexpr float CAP = static_cast<float>(WAL_CONTACT_CAP);
	if (!(d > 0.0f) || !(cr > 0.0f)) return CAP;
	const float d2 = d * d;
	const float d4 = d2 * d2;
	const float d6 = d4 * d2;
	const float inv_d12 = 1.0f / (d6 * d6);
	const float cr2 = cr * cr;
	const float cr4 = cr2 * cr2;
	const float cr6 = cr4 * cr2;
	const float inv_cr12 = 1.0f / (cr6 * cr6);
	const float raw = KWALL_F * (inv_d12 - inv_cr12);
	return (raw > CAP) ? CAP : raw;
}
