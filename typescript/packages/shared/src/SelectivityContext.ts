// SelectivityContext.ts — Cross-platform selectivity context types
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import type { ScientificProvenance } from './BindingPopulation.js';
import type { DeltaDeltaG } from './SelectivityAnalysis.js';

/** Summary of docking results for a single target. */
export interface TargetDockingSummary {
  /** Target name (e.g., "5HT2A", "D2R") */
  targetName: string;
  /** Best ensemble score; physical free energy only with supporting provenance. */
  bestFreeEnergy: number;
  /** Number of binding modes */
  modeCount: number;
  /** Configurational entropy (nats) */
  sConf: number;
  /** Vibrational entropy slot; physical units only with supporting provenance. */
  sVib: number;
  /** Whether entropy converged */
  isConverged: boolean;
  /** Cavity volume (cubic Angstroms, if known) */
  cavityVolume: number | null;
  /** Population size */
  populationSize: number;
  /** False when the source did not supply thermodynamic moments. */
  thermodynamicsAvailable?: boolean;
  /** Scientific provenance for this target's thermodynamic record. */
  scientificProvenance?: ScientificProvenance;
}

/** Context for multi-target selectivity analysis. */
export interface SelectivityContext {
  /** Ligand identifier */
  ligandName: string;
  /** Target summaries (typically 2-4 targets) */
  targets: TargetDockingSummary[];
  /** Pre-computed pairwise score differences. */
  deltaDeltas: DeltaDeltaG[];
}
