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
// Optional: set  FLEXAIDDS_RMSDST=<path>  to load reference coordinates from a
//   separate file (SDF or PDB format).  If unset, FlexAIDdS uses the
//   coordinates already loaded from the input ligand file (ideal when the input
//   IS the crystal SDF, as in DatasetRunner benchmark runs).
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
