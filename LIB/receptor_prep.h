// =============================================================================
// receptor_prep.h — Binding-site rotamer pre-relaxation for apo receptors
//
// Addresses the apo-strain false-minimum problem: the crystal binding pose
// requires small sidechain adjustments that the apo PDB has not made.  CF.wal
// hard-wall penalties inflate for near-native poses that overlap unmoved apo
// sidechains, causing peripheral decoys (that avoid those atoms) to score better.
//
// Algorithm: greedy single-residue Dunbrack 2010 backbone-independent rotamer
// search over all protein residues whose Cα falls within a given radius of the
// oracle binding site centroid.  Zero external dependencies — pure C++ stdlib.
//
// Usage from DatasetRunner::run():
//   #include "receptor_prep.h"
//   int n = receptor_prep::prep_receptor_rotamers(
//               entry.receptor_path, entry.binding_site_path,
//               out_dir + "/" + pdb_id + "_prepped.pdb");
//   if (n >= 0) effective_receptor = prepped;
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#pragma once
#include <string>

namespace receptor_prep {

/// Pre-relax binding-site sidechain rotamers.
///
/// @param receptor_pdb    Apo receptor PDB file (input)
/// @param oracle_site_pdb Oracle binding-site PDB (used only for centroid)
/// @param out_pdb         Output PDB with optimised sidechain coords
/// @param radius_ang      Cα-to-centroid radius for pocket residue selection (Å)
/// @param top_n           Number of Dunbrack rotamers to evaluate per residue
/// @param vdw_tol         VDW overlap tolerance in Å (overlap = rA+rB-tol-d)
///
/// @return Number of residues whose rotamer was changed (≥ 0 on success, -1 on error)
int prep_receptor_rotamers(
    const std::string& receptor_pdb,
    const std::string& oracle_site_pdb,
    const std::string& out_pdb,
    float radius_ang = 5.5f,
    int   top_n      = 5,
    float vdw_tol    = 0.40f);

} // namespace receptor_prep
