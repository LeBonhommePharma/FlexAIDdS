// =============================================================================
// native_score.h — IdealPopulation native-pose CF scoring diagnostic
//
// Evaluates the contact function (CF) at the crystal / reference ligand pose
// WITHOUT running the genetic algorithm.  Answers the question:
//   "Is the scorer broken for this ligand?"   (cf << 0 → scorer works)
//   "Did we find the right pose?"             (that's the GA's job, not this)
//
// Trigger:  set env var  FLEXAIDDS_SCORE_NATIVE=1  before launching FlexAIDdS.
//
// Coordinates: always uses atoms[].coor[] — the crystal ligand coordinates in
//   the PROCESSED internal frame (same frame as the cleft grid).  FLEXAIDDS_RMSDST
//   is intentionally ignored for coordinate input: that file holds raw
//   crystallographic coordinates which are in the wrong frame and cause +500M
//   steric-clash explosions.
//
// Output:  one line on stderr  (DatasetRunner parses the cf= field):
//   [NATIVE_CF] cf=<total> breakdown=com:<v>,wal:<v>,sas:<v>,con:<v>
//
// The GA continues normally after this diagnostic — no extra subprocess needed.
//
// Copyright 2026 Le Bonhomme Pharma.  Licensed under Apache-2.0.
// =============================================================================

#pragma once

#include "flexaid.h"
#include "Vcontacts.h"

/// Score the native (crystal / reference) ligand pose using the contact
/// function, without running the GA.
///
/// Called from top.cpp when env var FLEXAIDDS_SCORE_NATIVE=1 is set.
/// Prints one [NATIVE_CF] line to stderr and returns.  The GA then runs
/// normally; the subsequent ic2cf(FA->opt_par) call in top.cpp resets coor[].
///
/// @param FA         FlexAID global state (npar, map_par, num_grd, res_cnt …)
/// @param VC         Voronoi contacts global state
/// @param atoms      atom array  (ligand coor[] and dis/ang/dih are written
///                   temporarily, then restored before returning)
/// @param residue    residue array
/// @param cleftgrid  cleft grid-point array  (length FA->num_grd)
void score_native_pose(FA_Global* FA, VC_Global* VC, atom* atoms,
                       resid* residue, gridpoint* cleftgrid);

/// Populate atoms[].coor_ref from FLEXAIDDS_RMSDST (crystal SDF) and set
/// FA->refstructure=1 when FLEXAIDDS_DUMP_POP is truthy.
///
/// Audit-only: does NOT inject crystal coordinates into the GA search
/// (atoms[].coor, seeds, reference_ligand.file remain unchanged). Required
/// for .rrd / .pop.tsv gates in cluster.cpp that need calc_rmsd vs native.
///
/// Returns true if refstructure is active after the call (already set, or newly
/// loaded). Returns false if DUMP_POP is off, RMSDST missing, or load failed.
bool load_dump_pop_refstructure(FA_Global* FA, atom* atoms, resid* residue);

