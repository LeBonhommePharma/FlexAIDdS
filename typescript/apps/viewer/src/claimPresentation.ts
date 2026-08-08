// claimPresentation.ts — evidence-derived viewer claim gates
// SPDX-License-Identifier: Apache-2.0

import { allowsBindingClaims } from '@bonhomme/shared';
import type {
  BindingPopulation,
  CleftAssessment,
  RefereeVerdict,
  SelectivityAnalysis,
} from '@bonhomme/shared';

/** Re-derive the population claim level instead of trusting display metadata. */
export function populationAllowsPhysicalAffinity(population: BindingPopulation): boolean {
  return allowsBindingClaims(population.globalThermodynamics);
}

/** A pairwise affinity claim requires two independently valid target sources. */
export function selectivityAllowsPhysicalAffinity(analysis: SelectivityAnalysis): boolean {
  const sources = analysis.claimSources ?? [];
  return sources.length >= 2
    && sources.every((source) => allowsBindingClaims(source));
}

/** Fleet/referee payloads without their evidence-bearing source fail closed. */
export function refereeAllowsPhysicalAffinity(verdict: RefereeVerdict): boolean {
  return verdict.claimSource !== undefined && allowsBindingClaims(verdict.claimSource);
}

/**
 * A druggability verdict is a potency claim about an unmeasured molecule.
 * Without binding-physical evidence the tier is a geometry diagnostic only.
 */
export function cleftAllowsDruggabilityClaim(assessment: CleftAssessment): boolean {
  return assessment.claimSource !== undefined
    && allowsBindingClaims(assessment.claimSource);
}
