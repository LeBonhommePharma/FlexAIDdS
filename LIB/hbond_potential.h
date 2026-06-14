// hbond_potential.h — Angular-dependent hydrogen bond potential
//
// Implements a Gaussian bell potential for H-bond scoring that accounts for
// both donor-acceptor distance and D-H...A angle. Also differentiates
// standard H-bonds from salt bridges based on charge bin classification.
//
// Integration: called from vcfunction.cpp during the pairwise contact loop.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <cstdint>
#include "atom_typing_256.h"

// Forward-declare atom_struct to avoid circular include with flexaid.h
struct atom_struct;

namespace hbond {

// Compute the angle (in degrees) between vectors (a->b) and (a->c)
// where a, b, c are 3D coordinate arrays.
inline double angle_deg(const float* a, const float* b, const float* c) {
    double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
    double ac[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
    double dot = ab[0]*ac[0] + ab[1]*ac[1] + ab[2]*ac[2];
    double mag_ab = std::sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    double mag_ac = std::sqrt(ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2]);
    if (mag_ab < 1e-8 || mag_ac < 1e-8) return 0.0;
    double cos_theta = dot / (mag_ab * mag_ac);
    // Clamp to [-1, 1] for numerical safety
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    return std::acos(cos_theta) * 180.0 / 3.14159265358979323846;
}

// Find the index of a bonded hydrogen atom for a given donor atom.
// Returns the internal atom index of the H, or -1 if none found.
// The atom.bond[] array: bond[0] = count of bonded atoms, bond[1..6] = indices.
// H atoms are identified by element name "H" or single-letter match.
inline int find_bonded_hydrogen(const atom_struct* atoms, const atom_struct& donor) {
    int nbonds = donor.bond[0];
    for (int b = 1; b <= nbonds && b <= 6; ++b) {
        int idx = donor.bond[b];
        if (idx < 0) continue;
        // Check element name for hydrogen
        if (atoms[idx].element[0] == 'H' ||
            (atoms[idx].element[0] == ' ' && atoms[idx].element[1] == 'H') ||
            (atoms[idx].name[0] == 'H')) {
            return idx;
        }
    }
    return -1;
}

// Compute the angular-dependent H-bond energy between two contacting atoms.
//
// The Gaussian bell potential:
//   E_hb = weight * exp(-0.5 * ((d - d0) / sigma_d)^2)
//                 * exp(-0.5 * ((theta - theta0) / sigma_theta)^2)
//
// Salt bridge detection: if one atom is anionic (Q_ANIONIC) and the other
// cationic (Q_CATIONIC), use the salt_bridge_weight instead of hbond_weight.
//
// Parameters from FA_Global: use_hbond, hbond_optimal_dist, hbond_optimal_angle,
// hbond_sigma_dist, hbond_sigma_angle, hbond_weight, hbond_salt_bridge_weight.
//
// Returns 0.0 if neither atom is H-bond capable or if disabled.
inline double compute_hbond_energy(
    const atom_struct* atoms,
    int idx_a, int idx_b,
    double dist,
    double optimal_dist,
    double optimal_angle,
    double sigma_dist,
    double sigma_angle,
    double weight,
    double salt_bridge_weight)
{
    const atom_struct& a = atoms[idx_a];
    const atom_struct& b = atoms[idx_b];

    // Both atoms must be H-bond capable (donor on one side, acceptor on other).
    // FIX (v51): tightened from (!hb_a && !hb_b) — the old guard allowed the
    // Gaussian bonus to fire when only ONE atom was H-bond capable (e.g. N
    // touching a carbon), inflating scores for non-polar contacts.
    bool hb_a = atom256::get_hbond(a.type256);
    bool hb_b = atom256::get_hbond(b.type256);
    if (!hb_a || !hb_b) return 0.0;

    // Distance Gaussian component
    double dd = (dist - optimal_dist) / sigma_dist;
    double E_dist = std::exp(-0.5 * dd * dd);

    // Angular component: find the hydrogen to compute D-H...A angle
    // Try atom A as donor first, then atom B
    double best_angle_term = 0.0;

    if (hb_a) {
        int h_idx = find_bonded_hydrogen(atoms, a);
        if (h_idx >= 0) {
            // D-H...A angle: angle at H between D and A
            double theta = angle_deg(atoms[h_idx].coor, a.coor, b.coor);
            double da = (theta - optimal_angle) / sigma_angle;
            double term = std::exp(-0.5 * da * da);
            if (term > best_angle_term) best_angle_term = term;
        }
    }
    if (hb_b) {
        int h_idx = find_bonded_hydrogen(atoms, b);
        if (h_idx >= 0) {
            double theta = angle_deg(atoms[h_idx].coor, b.coor, a.coor);
            double da = (theta - optimal_angle) / sigma_angle;
            double term = std::exp(-0.5 * da * da);
            if (term > best_angle_term) best_angle_term = term;
        }
    }

    // If no hydrogen found on either side, use a reduced distance-only term
    // (acceptor-acceptor contacts can still form water-mediated bridges)
    if (best_angle_term == 0.0) {
        best_angle_term = 0.3; // reduced penalty for geometry-unknown contacts
    }

    // Determine weight: salt bridge vs standard H-bond
    double w = weight;
    uint8_t qbin_a = atom256::get_charge_bin(a.type256);
    uint8_t qbin_b = atom256::get_charge_bin(b.type256);
    if ((qbin_a == atom256::Q_NEGATIVE && qbin_b == atom256::Q_POSITIVE) ||
        (qbin_a == atom256::Q_POSITIVE && qbin_b == atom256::Q_NEGATIVE)) {
        w = salt_bridge_weight;
    }

    return w * E_dist * best_angle_term;
}

} // namespace hbond
