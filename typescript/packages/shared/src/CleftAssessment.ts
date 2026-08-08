// CleftAssessment.ts — Cross-platform cleft assessment types
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import type { ClaimValidity, ThermodynamicClaimSource } from './scientificProvenance.js';

/** Druggability tier for a binding cleft. */
export type DruggabilityTier = 'high' | 'moderate' | 'low' | 'undruggable';

/** Assessment of a binding cleft's druggability and properties. */
export interface CleftAssessment {
  /**
   * Geometric druggability score bucket.
   *
   * This is derived from pocket geometry alone. It is a potency-flavoured
   * label, so consumers must check `claimValidity` before presenting it as a
   * druggability verdict rather than a shape diagnostic.
   */
  druggability: DruggabilityTier;
  /** Summary of the cleft assessment */
  summary: string;
  /** Suggested ligand properties for this cleft */
  suggestedLigandProperties: string;
  /** Warnings about the cleft */
  warnings: string[];
  /** Derived claim level of the thermodynamic evidence behind this pocket. */
  claimValidity?: ClaimValidity;
  /** Evidence-bearing source; consumers must re-derive validity from it. */
  claimSource?: ThermodynamicClaimSource;
}
