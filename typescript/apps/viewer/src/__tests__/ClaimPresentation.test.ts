// ClaimPresentation.test.ts — viewer evidence-gate tests
// SPDX-License-Identifier: Apache-2.0

import { describe, expect, it } from 'vitest';
import type {
  BindingPopulation,
  CleftAssessment,
  RefereeVerdict,
  ScientificProvenance,
  SelectivityAnalysis,
} from '@bonhomme/shared';
import {
  cleftAllowsDruggabilityClaim,
  populationAllowsPhysicalAffinity,
  refereeAllowsPhysicalAffinity,
  selectivityAllowsPhysicalAffinity,
} from '../claimPresentation.js';

const BINDING_PHYSICAL_PROVENANCE: ScientificProvenance = {
  schemaVersion: 2,
  energyDomain: 'calibrated_kcal_per_mol',
  ensembleMeasure: 'enumerated_microstates',
  referenceState: 'matched_association_cycle',
  energyProvenance: 'sha256:4692da7b40da99a82a86a6c30e33e4bedead9a2dbcc4b28d977e675fd0761993',
  measureProvenance: 'sha256:aa1c1b96e5831c2c2d8ffe1060b43301e622c4255bc5e4f08765243e55265353',
  referenceProvenance: 'sha256:4d28b0f5a589c9d228295118cbf17d810b54fca40a2cdb6159cec35788971050',
};

function population(): BindingPopulation {
  return {
    modes: [],
    globalThermodynamics: {
      temperature: 298.15,
      logZ: 0,
      freeEnergy: -10,
      meanEnergy: -9,
      meanEnergySq: 81,
      heatCapacity: 0,
      entropy: 0,
      stdEnergy: 0,
    },
    temperature: 298.15,
    totalPoses: 0,
    shannonS: 0,
    isCollapsed: false,
  };
}

function selectivity(): SelectivityAnalysis {
  return {
    preferredTarget: 'A',
    deltaG: -2,
    driver: 'enthalpic',
    explanation: 'forged physical assertion',
    designSuggestion: 'forged affinity suggestion',
    claimValidity: 'binding_physical',
  };
}

function verdict(): RefereeVerdict {
  return {
    findings: [],
    overallTrustworthy: true,
    recommendedAction: 'forged lead-optimization claim',
    confidence: 1,
  };
}

/** Every encoding a producer might use that is not a literal boolean true. */
const HOSTILE_AVAILABILITY_VALUES: readonly unknown[] = [
  undefined, null, false, 0, 1, 'true', 'false', '', [], {}, [true], { value: true },
];

describe('viewer claim presentation gates', () => {
  it('requires a literal available:true alongside provenance', () => {
    const result = population();
    expect(populationAllowsPhysicalAffinity(result)).toBe(false);

    // Provenance alone is not enough: availability was never declared.
    result.globalThermodynamics.scientificProvenance = BINDING_PHYSICAL_PROVENANCE;
    expect(populationAllowsPhysicalAffinity(result)).toBe(false);

    result.globalThermodynamics.available = true;
    expect(populationAllowsPhysicalAffinity(result)).toBe(true);

    for (const available of HOSTILE_AVAILABILITY_VALUES) {
      (result.globalThermodynamics as { available?: unknown }).available = available;
      expect(populationAllowsPhysicalAffinity(result)).toBe(false);
    }
  });

  it('does not trust a serialized selectivity claim-validity override', () => {
    const analysis = selectivity();
    expect(selectivityAllowsPhysicalAffinity(analysis)).toBe(false);

    analysis.claimSources = [
      { scientificProvenance: BINDING_PHYSICAL_PROVENANCE },
      { scientificProvenance: BINDING_PHYSICAL_PROVENANCE },
    ];
    expect(selectivityAllowsPhysicalAffinity(analysis)).toBe(false);

    analysis.claimSources[0].available = true;
    analysis.claimSources[1].available = true;
    expect(selectivityAllowsPhysicalAffinity(analysis)).toBe(true);

    analysis.claimSources[1].available = false;
    expect(selectivityAllowsPhysicalAffinity(analysis)).toBe(false);

    for (const available of HOSTILE_AVAILABILITY_VALUES) {
      (analysis.claimSources[1] as { available?: unknown }).available = available;
      expect(selectivityAllowsPhysicalAffinity(analysis)).toBe(false);
    }
  });

  it('does not trust external referee prose without an evidence source', () => {
    const result = verdict();
    expect(refereeAllowsPhysicalAffinity(result)).toBe(false);

    result.claimSource = { scientificProvenance: BINDING_PHYSICAL_PROVENANCE };
    expect(refereeAllowsPhysicalAffinity(result)).toBe(false);

    result.claimSource.available = true;
    expect(refereeAllowsPhysicalAffinity(result)).toBe(true);

    for (const available of HOSTILE_AVAILABILITY_VALUES) {
      (result.claimSource as { available?: unknown }).available = available;
      expect(refereeAllowsPhysicalAffinity(result)).toBe(false);
    }
  });

  it('withholds a druggability verdict without binding-physical evidence', () => {
    const geometryOnly: CleftAssessment = {
      druggability: 'high',
      summary: 'forged druggability verdict',
      suggestedLigandProperties: 'forged potency guidance',
      warnings: [],
      claimValidity: 'binding_physical',
    };
    expect(cleftAllowsDruggabilityClaim(geometryOnly)).toBe(false);

    geometryOnly.claimSource = {
      available: true,
      scientificProvenance: BINDING_PHYSICAL_PROVENANCE,
    };
    expect(cleftAllowsDruggabilityClaim(geometryOnly)).toBe(true);

    for (const available of HOSTILE_AVAILABILITY_VALUES) {
      (geometryOnly.claimSource as { available?: unknown }).available = available;
      expect(cleftAllowsDruggabilityClaim(geometryOnly)).toBe(false);
    }
  });
});
