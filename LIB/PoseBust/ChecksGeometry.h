// ChecksGeometry.h — PoseBusters dock-suite geometry plausibility (clean-room)
//
// Binary check keys align with the public PoseBusters dock geometry suite:
//   bond_lengths_within_bounds, bond_angles_within_bounds, no_internal_clash,
//   flatness_passes (aromatic), flat_double_bonds
//
// Implementation uses independent covalent/vdW and hybridization heuristics;
// it does not copy RDKit distance-geometry or PoseBusters source.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Types.h"

#include <vector>

namespace flexaids::posebust {

// Relative widen for bond-length window around covalent-radius sum.
// Pass when |d − (r1+r2)| ≤ threshold_bad_bond_length * (r1+r2).
inline constexpr float kThresholdBadBondLength = 0.25f;

// Relative angle tolerance: pass when |Δθ| / θ_ideal ≤ this value.
inline constexpr float kThresholdBadAngle = 0.25f;

// Steric clash: heavy–heavy pairs with graph distance ≥ 4 must satisfy
// d ≥ kClashVdwScale * (vdW_i + vdW_j).
inline constexpr float kClashVdwScale = 0.70f;
inline constexpr int   kClashMinGraphDistance = 4;

// Planarity: max out-of-plane distance (Å) for aromatic rings.
inline constexpr float kThresholdFlatnessAngstrom = 0.25f;

// Double-bond planarity: |sin(φ)| of first-substituent torsion must be ≤ this.
// sin(14.5°) ≈ 0.25 — same numeric scale as the ring flatness threshold.
inline constexpr float kThresholdDoubleBondSinPhi = 0.25f;

/// Bondi-style van der Waals radius (Å) for atomic number Z.
/// Unknown / zero Z returns a conservative default (2.0 Å).
[[nodiscard]] float vdw_radius(int Z) noexcept;

/// Append binary geometry checks for intramolecular distance geometry:
///   - bond_lengths_within_bounds
///   - bond_angles_within_bounds
///   - no_internal_clash
/// Hydrogens are ignored for clashes; angles that only involve hydrogens
/// (or have H as the central atom) are skipped.
void check_distance_geometry(const Molecule& pred, std::vector<CheckItem>& out);

/// Append planarity checks:
///   - flatness_passes (aromatic / sp2 ring flatness)
///   - flat_double_bonds
void check_flatness(const Molecule& pred, std::vector<CheckItem>& out);

} // namespace flexaids::posebust
