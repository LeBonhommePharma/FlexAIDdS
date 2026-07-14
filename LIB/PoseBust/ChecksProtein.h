// ChecksProtein.h — Protein–ligand geometry checks for FlexAIDdS PoseBust
//
// Clean-room C++26 implementation (Apache-2.0). Independent of any third-party
// PoseBusters / RDKit / OpenBabel source. Binary check *keys* follow the public
// dock-suite vocabulary used for benchmark reporting.
//
// Depends on Types.h:
//   Atom:     { float x,y,z; int atomic_num; bool is_h; ... }
//   Molecule: { std::vector<Atom> atoms; ... }
//   CheckItem:{ std::string key, label; bool passed; std::string detail; }
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Types.h"

#include <vector>

namespace flexaids::posebust {

/// Bondi van der Waals radius in Å for atomic number Z.
/// Unknown / out-of-range Z returns a conservative default (1.70 Å).
[[nodiscard]] float vdw_radius(int Z) noexcept;

/// Intermolecular steric check between ligand and protein heavy atoms.
///
/// Appends to \p out:
///   - "minimum_distance_to_protein": min heavy–heavy lig–prot distance
///     (absolute floor + relative vdW scale)
///     (soft vdW scale 0.75·(ri+rj) is reported in detail only).
///     Pairs use a 5 Å cell list, or brute-force O(n·m) when n·m < 5e6.
void check_intermolecular_distance(const Molecule& ligand,
                                   const Molecule& protein,
                                   std::vector<CheckItem>& out);

/// Approximate volume-overlap, pocket-presence, and cofactor/water keys.
///
/// Appends to \p out:
///   - "volume_overlap_with_protein"
///   - "protein-ligand_maximum_distance" (pocket presence within 5 Å)
///   - cofactor/water min-distance + volume keys (vacuous pass on apo crop)
void check_volume_overlap(const Molecule& ligand,
                          const Molecule& protein,
                          std::vector<CheckItem>& out);

}  // namespace flexaids::posebust
