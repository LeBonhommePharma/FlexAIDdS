// PoseQualityContext.ts — Cross-platform pose quality context types
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import type { ScientificProvenance } from './BindingPopulation.js';

/** Summary of a single pose for quality analysis. */
export interface PoseProfile {
  /** Rank within the binding mode (0 = best) */
  rank: number;
  /** Complementarity-function diagnostic in arbitrary CF units. */
  cfScore: number;
  /** Boltzmann weight (fraction, 0.0-1.0) */
  boltzmannWeight: number;
  /** RMSD to the mode centroid (Angstroms) */
  rmsdToCentroid: number;
}

/** Full pose quality context for a binding mode. */
export interface PoseQualityContext {
  /** Binding mode index */
  modeIndex: number;
  /** Top poses (up to 5), sorted by Boltzmann weight descending */
  topPoses: PoseProfile[];
  /** Total poses in this mode */
  totalPoses: number;
  /** Weight of the dominant pose (highest Boltzmann weight) */
  dominantPoseWeight: number;
  /** Whether the top pose dominates (weight > 0.5) */
  hasDominantPose: boolean;
  /** RMSD spread of top 5 poses (max - min RMSD to centroid) */
  rmsdSpread: number;
  /** Whether CF score ranking matches Boltzmann weight ranking */
  scoreWeightAligned: boolean;
  /** Spearman rank correlation between CF score and Boltzmann weight (top 5) */
  scoreWeightCorrelation: number;
  /** Mode score; physical free energy only with supporting provenance. */
  modeFreeEnergy: number;
  /** False when the source did not supply thermodynamic moments. */
  thermodynamicsAvailable?: boolean;
  /** Scientific provenance for the mode thermodynamic record. */
  scientificProvenance?: ScientificProvenance;
}
