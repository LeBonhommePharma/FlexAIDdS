// hbond_potential.h — Angular-dependent hydrogen bond potential
//
// Implements a Gaussian bell potential for H-bond scoring that accounts for
// both donor-acceptor distance and D-H...A angle. Also differentiates
// standard H-bonds from salt bridges based on atom charges.
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

inline bool is_hydrogen_atom(const atom_struct& atom) {
    return atom.element[0] == 'H' ||
           (atom.element[0] == ' ' && atom.element[1] == 'H') ||
           atom.name[0] == 'H';
}

// Find the index of a bonded hydrogen atom for a given donor atom. Returns the
// internal atom index of the H, or -1 if none found. The atom.bond[] array uses
// bond[0] as count and bond[1..6] as internal atom indices.
inline int find_bonded_hydrogen(const atom_struct* atoms, const atom_struct& donor) {
    int nbonds = donor.bond[0];
    for (int b = 1; b <= nbonds && b <= 6; ++b) {
        int idx = donor.bond[b];
        if (idx < 0) continue;
        if (is_hydrogen_atom(atoms[idx])) {
            return idx;
        }
    }
    return -1;
}

// Heavy-atom fallback for input structures without explicit hydrogens. Place a
// virtual donor H opposite the normalized sum of donor-heavy-neighbour vectors.
// This is conservative: if the donor topology is missing or symmetric enough
// that the direction is undefined, the caller falls back to the reduced
// distance-only term.
inline bool virtual_hydrogen_coord(const atom_struct* atoms,
                                   const atom_struct& donor,
                                   float out[3]) {
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    int n_heavy = 0;

    int nbonds = donor.bond[0];
    for (int b = 1; b <= nbonds && b <= 6; ++b) {
        int idx = donor.bond[b];
        if (idx < 0) continue;
        const atom_struct& nb = atoms[idx];
        if (is_hydrogen_atom(nb)) continue;

        double dx = nb.coor[0] - donor.coor[0];
        double dy = nb.coor[1] - donor.coor[1];
        double dz = nb.coor[2] - donor.coor[2];
        double len = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (len < 1e-8) continue;
        vx += dx / len;
        vy += dy / len;
        vz += dz / len;
        ++n_heavy;
    }

    if (n_heavy == 0) return false;
    double mag = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (mag < 1e-8) return false;

    constexpr float VIRTUAL_DH_LEN = 1.0f;
    out[0] = donor.coor[0] - VIRTUAL_DH_LEN * static_cast<float>(vx / mag);
    out[1] = donor.coor[1] - VIRTUAL_DH_LEN * static_cast<float>(vy / mag);
    out[2] = donor.coor[2] - VIRTUAL_DH_LEN * static_cast<float>(vz / mag);
    return true;
}

inline double donor_angle_term(const atom_struct* atoms,
                               const atom_struct& donor,
                               const atom_struct& acceptor,
                               double optimal_angle,
                               double sigma_angle) {
    int h_idx = find_bonded_hydrogen(atoms, donor);
    const float* h_coord = nullptr;
    float virtual_h[3] = {0.0f, 0.0f, 0.0f};

    if (h_idx >= 0) {
        h_coord = atoms[h_idx].coor;
    } else if (virtual_hydrogen_coord(atoms, donor, virtual_h)) {
        h_coord = virtual_h;
    }

    if (!h_coord) return 0.0;

    double theta = angle_deg(h_coord, donor.coor, acceptor.coor);
    double da = (theta - optimal_angle) / sigma_angle;
    return std::exp(-0.5 * da * da);
}

inline bool is_charge_supported_salt_bridge(const atom_struct& a,
                                            const atom_struct& b) {
    constexpr float Q_SALT = 0.30f;
    return (a.charge <= -Q_SALT && b.charge >= Q_SALT) ||
           (b.charge <= -Q_SALT && a.charge >= Q_SALT);
}

// Compute the angular-dependent H-bond energy between two contacting atoms.
//
// The Gaussian bell potential:
//   E_hb = weight * exp(-0.5 * ((d - d0) / sigma_d)^2)
//                 * exp(-0.5 * ((theta - theta0) / sigma_theta)^2)
//
// Salt bridge detection: if one atom is anionic and the other cationic, use
// the salt_bridge_weight instead of hbond_weight.
//
// Parameters from FA_Global: use_hbond, hbond_optimal_dist, hbond_optimal_angle,
// hbond_sigma_dist, hbond_sigma_angle, hbond_weight, hbond_salt_bridge_weight.
//
// Returns 0.0 unless the encoded roles form a donor-acceptor pair.
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

    const bool donor_a = atom256::get_hbond_donor(a.type256);
    const bool acceptor_a = atom256::get_hbond_acceptor(a.type256);
    const bool donor_b = atom256::get_hbond_donor(b.type256);
    const bool acceptor_b = atom256::get_hbond_acceptor(b.type256);

    const bool a_to_b = donor_a && acceptor_b;
    const bool b_to_a = donor_b && acceptor_a;
    if (!a_to_b && !b_to_a) return 0.0;

    // Distance Gaussian component
    double dd = (dist - optimal_dist) / sigma_dist;
    double E_dist = std::exp(-0.5 * dd * dd);

    // Angular component: use explicit H when available, otherwise a
    // conservative virtual-H direction from donor heavy-atom topology.
    double best_angle_term = 0.0;

    if (a_to_b) {
        double term = donor_angle_term(atoms, a, b, optimal_angle, sigma_angle);
        if (term > best_angle_term) best_angle_term = term;
    }
    if (b_to_a) {
        double term = donor_angle_term(atoms, b, a, optimal_angle, sigma_angle);
        if (term > best_angle_term) best_angle_term = term;
    }

    // If no explicit or virtual hydrogen direction is available, use a reduced
    // distance-only term instead of pretending angular information exists.
    if (best_angle_term == 0.0) {
        best_angle_term = 0.3; // reduced penalty for geometry-unknown contacts
    }

    // Determine weight: salt bridge vs standard H-bond. Use actual charges,
    // not type256: v56 uses those top bits for donor/acceptor roles.
    double w = weight;
    if (is_charge_supported_salt_bridge(a, b)) {
        w = salt_bridge_weight;
    }

    double energy = w * E_dist * best_angle_term;
    constexpr double HBOND_PAIR_MIN = -2.0;
    if (energy < HBOND_PAIR_MIN) return HBOND_PAIR_MIN;
    return energy;
}

} // namespace hbond
