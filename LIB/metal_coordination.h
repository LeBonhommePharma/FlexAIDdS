// metal_coordination.h — Geometry-aware metal ion coordination potential
//
// Implements a Gaussian well potential for metal-ligand coordination bonds,
// with per-metal/per-donor ideal distances and well depths derived from
// Cambridge Structural Database metalloprotein surveys (Harding 2001/2006,
// Acta Cryst D57:401-411, D62:678-682) and Li-Merz 12-6-4 LJ-type parameters
// (JCTC 2014 10:289-297).
//
// Gaussian well (NOT Morse):
//   E(r) = D * exp[ -((r - r0) / sigma)^2 ]
//
// Rationale for Gaussian over Morse:
//   1. Metal-ligand bonds are predominantly electrostatic/ion-dipole, not
//      covalent — the Morse repulsive wall double-counts with the existing
//      Lennard-Jones wall term in vcfunction.cpp.
//   2. Consistent with hbond_potential.h which also uses a Gaussian bell.
//   3. The sigma parameter is more intuitive than Morse alpha (sigma in A).
//
// Well depths represent the NON-ELECTROSTATIC component of coordination
// (charge transfer, polarization, orbital covalency). When Coulomb is
// active, these supplement rather than replace the electrostatic term.
// Calibrated so that total interaction energies match Li-Merz 12-6-4
// ion-water binding benchmarks and AutoDock4Zn validation sets.
//
// Also provides coordination-number tracking and a positive quadratic
// CN-deviation penalty applied after the pairwise loop completes.
//
// Integration: called from vcfunction.cpp during the pairwise contact loop.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <array>
#include <optional>

// Need the full atom_struct definition for accessing .type member
#include "flexaid.h"

namespace metal_coord {

// ─── Coordination geometry types ────────────────────────────────────────────
enum class Geometry : uint8_t {
    TETRAHEDRAL,           // CN=4, ideal angle=109.5
    SQUARE_PLANAR,         // CN=4, ideal angle=90.0
    TRIGONAL_BIPYRAMIDAL,  // CN=5, ideal angle=90/120
    OCTAHEDRAL,            // CN=6, ideal angle=90.0
    PENTAGONAL_BIPYRAMIDAL,// CN=7, ideal angle=72.0/90.0
    SQUARE_ANTIPRISMATIC   // CN=8, ideal angle=70.5/99.6
};

// ─── Per-metal parameters ───────────────────────────────────────────────────
struct MetalParams {
    int      sybyl_type;    // SYBYL type (28=MG, 35=ZN, 36=CA, 37=FE)
    int      ideal_cn;      // preferred coordination number
    Geometry geometry;      // preferred geometry
    double   angle_primary; // primary ideal L-M-L angle (degrees)
};

// Lookup table: indexed by SYBYL type.  Only entries 28-38 are metals.
// Returns nullptr if the SYBYL type is not a recognized metal.
inline const MetalParams* get_metal_params(int sybyl_type) noexcept {
    // SYBYL: 28=MG, 29=SR, 30=CU, 31=MN, 32=HG, 33=CD, 34=NI, 35=ZN, 36=CA, 37=FE, 38=CO
    static constexpr MetalParams table[] = {
        {28, 6, Geometry::OCTAHEDRAL,            90.0},  // Mg2+
        {29, 0, Geometry::OCTAHEDRAL,             0.0},  // Sr -- no params (placeholder)
        {30, 4, Geometry::TETRAHEDRAL,           109.5}, // Cu2+
        {31, 6, Geometry::OCTAHEDRAL,            90.0},  // Mn2+
        {32, 0, Geometry::TETRAHEDRAL,            0.0},  // Hg -- placeholder
        {33, 0, Geometry::OCTAHEDRAL,             0.0},  // Cd -- placeholder
        {34, 6, Geometry::OCTAHEDRAL,            90.0},  // Ni2+
        {35, 4, Geometry::TETRAHEDRAL,           109.5}, // Zn2+
        {36, 7, Geometry::PENTAGONAL_BIPYRAMIDAL, 72.0}, // Ca2+
        {37, 6, Geometry::OCTAHEDRAL,            90.0},  // Fe2+/3+
        {38, 6, Geometry::OCTAHEDRAL,            90.0},  // Co2+
    };
    if (sybyl_type < 28 || sybyl_type > 38) return nullptr;
    const MetalParams& p = table[sybyl_type - 28];
    if (p.ideal_cn == 0) return nullptr;  // placeholder entry
    return &p;
}

// ─── Per metal-donor pair: ideal distance and well depth ────────────────────
struct DonorAffinity {
    double ideal_dist;   // ideal M-L distance (Angstroms)
    double well_depth;   // energy well depth (kcal/mol, negative = favorable)
};

// Donor SYBYL types:
//  6: N.1, 7: N.2, 8: N.3, 9: N.4 (EXCLUDED), 10: N.AR, 11: N.AM, 12: N.PL3
// 13: O.2 (carbonyl)  14: O.3 (hydroxyl/ether)  15: O.CO2 (carboxylate)
// 16: O.AR  17: S.2  18: S.3  22: P.3
//
// Well depth rationale:
// These represent the NON-electrostatic component of coordination energy:
// charge transfer, polarization, and orbital covalency. The electrostatic
// component is handled by the Coulomb term when charges are available.
//
// Calibration:
// - AutoDock4Zn uses 0.35 kcal/mol VDW well (Coulomb does most work)
// - Li-Merz 12-6-4: C4 term provides ~2-5 kcal/mol charge-induced-dipole
// - We use -3 to -6 kcal/mol to capture polarization + charge transfer
//   that Coulomb alone cannot model (especially for soft acid/base pairs
//   like Zn-S and Fe-N)
//
// Distances from:
// - Harding (2001) Acta Cryst D57:401-411 (Table 1-4)
// - Harding (2006) Acta Cryst D62:678-682

// Returns the DonorAffinity for a (metal, donor) pair, or nullopt if unknown.
inline std::optional<DonorAffinity> get_donor_affinity(int metal_sybyl,
                                                        int donor_sybyl) noexcept {
    // ── N.4 exclusion: quaternary N has NO lone pair, cannot coordinate ──
    if (donor_sybyl == 9) return std::nullopt;

    // ── N donor strength multipliers ──────────────────────────────────────
    // Not all N types coordinate equally. Applied to generic N well_depth.
    // N.AR (10, His imidazole): 1.0 (strong, sigma-donor)
    // N.PL3 (12): 0.6 (moderate)
    // N.1/N.2/N.3 (6-8): 0.3 (rare in proteins)
    // N.AM (11, backbone amide): 0.2 (very weak, lone pair delocalized into C=O)
    auto n_strength = [](int donor) -> double {
        switch (donor) {
            case 10: return 1.0;   // N.AR (His)
            case 12: return 0.6;   // N.PL3
            case 11: return 0.2;   // N.AM
            default: return 0.3;   // N.1, N.2, N.3
        }
    };

    // ── S donor filtering ──────────────────────────────────────────────────
    // Only S.3 (thiolate/thioether) and S.2 (thioether) coordinate metals.
    // S.O (sulfoxide) coordinates via O, not S; S.O2 (sulfone) is inert;
    // S.AR (thiophene) rarely coordinates in proteins.
    if (donor_sybyl >= 19 && donor_sybyl <= 21) return std::nullopt;  // S.O, S.O2, S.AR

    // ── Ca2+ (36) — strongly prefers O donors ──────────────────────────────
    // Harding (2001): Ca-O(carboxylate mono) 2.36, Ca-O(carbonyl) 2.36-2.40,
    //                 Ca-O(hydroxyl) 2.40-2.48, Ca-O(water) 2.39
    // Ca-N coordination is weak (rare, only in EDTA-like chelates)
    if (metal_sybyl == 36) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.36, -5.0};   // O.CO2 (carboxylate)
            case 13: return DonorAffinity{2.38, -4.0};   // O.2 (carbonyl)
            case 14: return DonorAffinity{2.43, -3.5};   // O.3 (hydroxyl/water)
            case 16: return DonorAffinity{2.40, -3.0};   // O.AR
            case 22: return DonorAffinity{2.55, -3.5};   // P.3 (phosphate)
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.50, -1.5 * n_strength(donor_sybyl)};
                if (donor_sybyl == 17 || donor_sybyl == 18)  // S.2, S.3
                    return DonorAffinity{2.80, -1.5};  // Ca-S very rare
                return std::nullopt;
        }
    }

    // ── Zn2+ (35) — coordinates N, O, S ────────────────────────────────────
    // Harding (2001): Zn-O(carboxylate) 2.00-2.05, Zn-N(His) 2.05,
    //                 Zn-S(Cys) 2.30, CN=4 tetrahedral
    // Zn-S has significant covalent character (soft-soft HSAB)
    if (metal_sybyl == 35) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.00, -4.0};   // O.CO2
            case 13: return DonorAffinity{2.05, -3.5};   // O.2
            case 14: return DonorAffinity{2.10, -3.0};   // O.3
            case 18: return DonorAffinity{2.30, -6.0};   // S.3 (Cys thiolate)
            case 17: return DonorAffinity{2.30, -5.0};   // S.2
            case 22: return DonorAffinity{2.15, -3.0};   // P.3
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.05, -5.0 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    // ── Mg2+ (28) — strict octahedral, strong O preference ─────────────────
    // Harding (2001): Mg-O(carboxylate) 2.06, Mg-O(water) 2.08-2.09
    // Mg2+ is a hard acid — negligible covalent character
    if (metal_sybyl == 28) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.06, -4.5};   // O.CO2
            case 13: return DonorAffinity{2.07, -4.0};   // O.2
            case 14: return DonorAffinity{2.09, -3.5};   // O.3
            case 22: return DonorAffinity{2.10, -4.0};   // P.3
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.20, -2.0 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    // ── Fe2+/3+ (37) — coordinates N, O, S ─────────────────────────────────
    // Harding: Fe-O(carboxylate) 2.08, Fe-S(Cys) 2.25-2.30
    // Fe-S/N have significant covalent character (d-orbital overlap)
    if (metal_sybyl == 37) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.08, -4.5};   // O.CO2
            case 13: return DonorAffinity{2.10, -4.0};   // O.2
            case 14: return DonorAffinity{2.12, -3.5};   // O.3
            case 18: return DonorAffinity{2.30, -6.0};   // S.3
            case 17: return DonorAffinity{2.30, -5.0};   // S.2
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.15, -5.0 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    // ── Cu2+ (30) — strong N/S coordinator ─────────────────────────────────
    if (metal_sybyl == 30) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{1.97, -4.0};   // O.CO2
            case 13: return DonorAffinity{1.97, -3.5};   // O.2
            case 14: return DonorAffinity{2.00, -3.0};   // O.3
            case 18: return DonorAffinity{2.15, -5.5};   // S.3
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.02, -5.0 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    // ── Mn2+ (31) — octahedral, moderate O/N ───────────────────────────────
    if (metal_sybyl == 31) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.14, -4.0};   // O.CO2
            case 13: return DonorAffinity{2.15, -3.5};   // O.2
            case 14: return DonorAffinity{2.18, -3.0};   // O.3
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.22, -2.0 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    // ── Ni2+ (34) — octahedral, strong N/S ─────────────────────────────────
    if (metal_sybyl == 34) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.04, -4.0};   // O.CO2
            case 13: return DonorAffinity{2.06, -3.5};   // O.2
            case 14: return DonorAffinity{2.08, -3.0};   // O.3
            case 18: return DonorAffinity{2.20, -5.5};   // S.3
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.10, -4.5 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    // ── Co2+ (38) — octahedral, moderate-strong N/S ────────────────────────
    if (metal_sybyl == 38) {
        switch (donor_sybyl) {
            case 15: return DonorAffinity{2.08, -4.0};   // O.CO2
            case 13: return DonorAffinity{2.10, -3.5};   // O.2
            case 14: return DonorAffinity{2.12, -3.0};   // O.3
            case 18: return DonorAffinity{2.25, -5.5};   // S.3
            default:
                if (donor_sybyl >= 6 && donor_sybyl <= 12)
                    return DonorAffinity{2.15, -3.5 * n_strength(donor_sybyl)};
                return std::nullopt;
        }
    }

    return std::nullopt;
}

// ─── Query: is this SYBYL type a metal ion? ────────────────────────────────
inline bool is_metal_type(int sybyl_type) noexcept {
    return get_metal_params(sybyl_type) != nullptr;
}

// ─── Query: is this SYBYL type a potential coordination donor? ──────────────
inline bool is_coord_donor_type(int sybyl_type) noexcept {
    // N types (6-12, but NOT N.4=9), O types (13-16), S.2/S.3 (17-18), P.3 (22)
    if (sybyl_type == 9) return false;  // N.4: no lone pair
    if (sybyl_type >= 6 && sybyl_type <= 18) return true;
    if (sybyl_type == 22) return true;  // P.3
    return false;
}

// ─── Gaussian well potential for metal-ligand coordination ─────────────────
//
// E(r) = D * exp[ -((r - r0) / sigma)^2 ]
//
// where D = well_depth (negative = favorable), r0 = ideal distance,
// sigma = Gaussian width (controls coordination shell tolerance).
//
// Properties:
//   - Minimum at r = r0 with E = D (the well depth)
//   - E -> 0 as r -> infinity (smooth decay, no hard cutoff)
//   - E -> 0 as r -> 0 (NO repulsive wall — handled by existing wall term)
//   - sigma ~0.3-0.4 A for tight coordination (Zn2+), ~0.5-0.6 for loose (Ca2+)
//
inline double gaussian_well(double r, double r0, double well_depth,
                            double sigma) noexcept {
    double dx = (r - r0) / sigma;
    return well_depth * std::exp(-(dx * dx));
}

// ─── Distance cutoff for performance ───────────────────────────────────────
// Beyond r0 + 2.5*sigma, the Gaussian is < 0.2% of well depth — skip.
inline constexpr double GAUSSIAN_CUTOFF_NSIGMA = 2.5;

// ─── Main scoring function ──────────────────────────────────────────────────
//
// Computes metal coordination energy for a pairwise contact between atoms
// idx_a and idx_b.  Automatically determines which (if any) is the metal
// and which is the donor.
//
// Returns 0.0 if neither atom is a metal, if the metal-donor pair has
// no parameterized affinity, or if the distance exceeds the cutoff.
//
// Parameters:
//   atoms:  full atom array
//   idx_a, idx_b: indices of the contacting atoms
//   dist:   distance between atoms (Angstroms)
//   weight: global weight multiplier (FA->metal_coord_weight)
//   sigma:  Gaussian width (FA->metal_coord_sigma, default 0.45 A)
//
inline double compute_metal_coord_energy(
    const atom_struct* atoms,
    int idx_a, int idx_b,
    double dist,
    double weight,
    double sigma) noexcept
{
    int type_a = atoms[idx_a].type;
    int type_b = atoms[idx_b].type;

    // Identify metal and donor
    int metal_type = 0, donor_type = 0;
    if (is_metal_type(type_a) && is_coord_donor_type(type_b)) {
        metal_type = type_a;
        donor_type = type_b;
    } else if (is_metal_type(type_b) && is_coord_donor_type(type_a)) {
        metal_type = type_b;
        donor_type = type_a;
    } else {
        return 0.0;  // not a metal-donor pair
    }

    // Look up affinity parameters
    auto aff = get_donor_affinity(metal_type, donor_type);
    if (!aff) return 0.0;

    // Distance cutoff: skip if far beyond coordination shell
    if (dist > aff->ideal_dist + GAUSSIAN_CUTOFF_NSIGMA * sigma)
        return 0.0;

    // Compute Gaussian well potential
    return weight * gaussian_well(dist, aff->ideal_dist, aff->well_depth, sigma);
}

// ─── Coordination number penalty ────────────────────────────────────────────
//
// Quadratic penalty for deviating from the ideal coordination number.
// Applied per metal atom after the pairwise loop, NOT inside it.
//
// Returns POSITIVE energy (unfavorable) for any deviation from ideal CN.
// In FlexAIDdS, lower CF = better, so this penalty makes under- or over-
// coordinated metals score worse.
//
// E_cn = cn_weight * (actual_cn - ideal_cn)^2
//
inline double cn_penalty(int actual_cn, int ideal_cn, double cn_weight) noexcept {
    int delta = actual_cn - ideal_cn;
    return cn_weight * static_cast<double>(delta * delta);
}

// ─── Coordination cutoff distance ───────────────────────────────────────────
// A contact is "coordinating" if dist < ideal_dist + COORD_CUTOFF_MARGIN
inline constexpr double COORD_CUTOFF_MARGIN = 0.8;  // Angstroms

inline bool is_coordinating(double dist, double ideal_dist) noexcept {
    return dist < (ideal_dist + COORD_CUTOFF_MARGIN);
}

} // namespace metal_coord
