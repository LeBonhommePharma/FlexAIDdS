// SelectivityAnalysis.ts — Cross-platform selectivity analysis types
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import type {
  ClaimValidity,
  ThermodynamicClaimSource,
} from './BindingPopulation.js';

/** Driver of selectivity between targets. */
export type SelectivityDriver = 'enthalpic' | 'entropic' | 'mixed' | 'inconclusive';

/** Numeric difference between two target ensemble scores. */
export interface DeltaDeltaG {
  /** First target identifier */
  targetA: string;
  /** Second target identifier */
  targetB: string;
  /** Score difference; physical delta-delta-G only for binding-physical inputs. */
  ddg: number;
}

/** Selectivity analysis across targets. */
export interface SelectivityAnalysis {
  /** Preferred target identifier */
  preferredTarget: string;
  /** Numeric difference between the top two target scores. */
  deltaG: number;
  /** Thermodynamic driver of selectivity */
  driver: SelectivityDriver;
  /** Explanation of selectivity */
  explanation: string;
  /** Suggestion for improving selectivity by design */
  designSuggestion: string;
  /** Derived scientific claim level for the target comparison. */
  claimValidity?: ClaimValidity;
  /** Evidence-bearing target sources; UI consumers must re-derive validity. */
  claimSources?: ThermodynamicClaimSource[];
}
