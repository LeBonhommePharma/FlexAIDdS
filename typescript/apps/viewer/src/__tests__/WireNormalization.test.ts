// WireNormalization.test.ts — producer -> normalizer -> viewer golden tests
//
// These pin the whole path a real payload takes: a Python-shaped docking
// result (snake_case, unit-suffixed keys, nested provenance) is pushed through
// the single shared normalizer and then rendered by the viewer analyzers. The
// point is that no dialect and no hostile availability encoding can reach the
// UI as a physical affinity, potency or Kd claim.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import { describe, expect, it } from 'vitest';
import {
  claimValidityForRecord,
  normalizeDockingRecord,
  type BindingPopulation,
  type NormalizedDockingRecord,
  type ScientificProvenance,
  type ThermodynamicClaimSource,
} from '@bonhomme/shared';
import { IntelligenceEngine, RuleBasedReferee } from '../IntelligenceEngine.js';
import { BindingModeAnalyzer } from '../intelligence/BindingModeAnalyzer.js';
import { SelectivityAnalyzer } from '../intelligence/SelectivityAnalyzer.js';
import {
  populationAllowsPhysicalAffinity,
  refereeAllowsPhysicalAffinity,
} from '../claimPresentation.js';

const ENERGY_SHA = 'sha256:4692da7b40da99a82a86a6c30e33e4bedead9a2dbcc4b28d977e675fd0761993';
const MEASURE_SHA = 'sha256:aa1c1b96e5831c2c2d8ffe1060b43301e622c4255bc5e4f08765243e55265353';
const REFERENCE_SHA = 'sha256:4d28b0f5a589c9d228295118cbf17d810b54fca40a2cdb6159cec35788971050';

const BINDING_PHYSICAL_PROVENANCE: ScientificProvenance = {
  schemaVersion: 2,
  energyDomain: 'calibrated_kcal_per_mol',
  ensembleMeasure: 'enumerated_microstates',
  referenceState: 'matched_association_cycle',
  energyProvenance: ENERGY_SHA,
  measureProvenance: MEASURE_SHA,
  referenceProvenance: REFERENCE_SHA,
};

/** Every non-`true` availability encoding a producer might emit. */
const HOSTILE_AVAILABILITY_VALUES: readonly unknown[] = [
  undefined, null, false, 0, 1, 'true', 'false', '', [], {}, [true], { value: true },
];

/**
 * Shape emitted by python/flexaidds DockingResult.to_json(): snake_case mode
 * list, flat legacy fields, and a nested ThermodynamicBreakdown.to_dict().
 * Note there is no `available` key anywhere — Python never emits one.
 */
function pythonPayload(options: {
  provenance?: Record<string, unknown>;
  available?: unknown;
} = {}): Record<string, unknown> {
  const provenance = options.provenance ?? {
    schema_version: 2,
    energy_domain: 'cf_arbitrary_units',
    ensemble_measure: 'optimizer_samples',
    reference_state: 'none',
    energy_provenance: 'FlexAID CF contact-function score',
    measure_provenance: 'GA optimizer samples',
    reference_provenance: '',
    claim_validity: 'binding_physical', // forged label; must be ignored
  };
  const breakdown: Record<string, unknown> = {
    temperature_K: 298.15,
    logZ_config: 3.25,
    G_config_kcal_mol: -12.4,
    H_eff_kcal_mol: -13.2,
    S_config_kcal_mol_K: 0.02,
    minus_T_S_config_kcal_mol: 0.8,
    Cv_kcal_mol_K: 0.021,
    sigma_E_kcal_mol: 1.1,
    G_vib_kcal_mol: 0.0,
    G_natural_kcal_mol: 0.0,
    G_other_kcal_mol: 0.0,
    G_total_kcal_mol: -12.4,
    scientific_provenance: provenance,
  };
  if ('available' in options) breakdown.available = options.available;

  return {
    source_dir: '/tmp/results',
    temperature: 298.15,
    n_modes: 1,
    metadata: {},
    binding_modes: [{
      mode_id: 0,
      rank: 1,
      n_poses: 12,
      free_energy: -12.4,
      proxy_free_energy: -12.4,
      soft_beta_G: -12.1,
      enthalpy: -13.2,
      entropy: 0.02,
      heat_capacity: 0.021,
      std_energy: 1.1,
      best_cf: -142.7,
      temperature: 298.15,
      scientific_provenance: provenance,
      best_pose_path: '/tmp/results/mode0_pose1.pdb',
      thermodynamics: breakdown,
    }],
    global_thermodynamics: breakdown,
  };
}

/** Adapt the canonical record onto the viewer's BindingPopulation shape. */
function toPopulation(record: NormalizedDockingRecord): BindingPopulation {
  return {
    modes: record.modes.map((mode) => ({
      index: mode.index,
      size: mode.size,
      freeEnergy: mode.freeEnergy,
      entropy: mode.entropy,
      enthalpy: mode.enthalpy,
      heatCapacity: mode.heatCapacity,
      probability: mode.probability,
      thermodynamics: mode.thermodynamics,
    })),
    globalThermodynamics: record.globalThermodynamics,
    temperature: record.temperature,
    totalPoses: record.totalPoses,
    shannonS: 0.4,
    isCollapsed: false,
  };
}

describe('python producer -> shared normalizer -> viewer', () => {
  it('carries numbers through unchanged', () => {
    const record = normalizeDockingRecord(pythonPayload());
    expect(record.temperature).toBe(298.15);
    expect(record.temperatureKnown).toBe(true);
    expect(record.modes[0].size).toBe(12);
    expect(record.modes[0].freeEnergy).toBe(-12.4);
    expect(record.modes[0].enthalpy).toBe(-13.2);
    expect(record.globalThermodynamics.freeEnergy).toBe(-12.4);
    expect(record.globalThermodynamics.entropy).toBe(0.02);
    expect(record.globalThermodynamics.heatCapacity).toBe(0.021);
  });

  it('renders default Python output as a proxy diagnostic in every surface', async () => {
    const population = toPopulation(normalizeDockingRecord(pythonPayload()));

    expect(populationAllowsPhysicalAffinity(population)).toBe(false);

    const analysis = await IntelligenceEngine.analyze(population);
    expect(analysis.claimValidity).toBe('proxy_only');
    const oracleText = `${analysis.bullets.join('\n')}\n${analysis.inputSummary}`;
    expect(oracleText).not.toContain('kcal/mol');
    expect(oracleText).not.toMatch(/(?:Strong|Moderate|Weak) binding affinity/);
    expect(oracleText).toMatch(/physical affinity unavailable/i);

    const verdict = RuleBasedReferee.referee(population);
    expect(refereeAllowsPhysicalAffinity(verdict)).toBe(false);
    expect(verdict.findings.find((f) => f.category === 'affinity')?.title)
      .toBe('Physical affinity unavailable');
    const verdictText = verdict.findings.map((f) => f.detail).join('\n');
    expect(verdictText).not.toContain('kcal/mol');

    const narrative = BindingModeAnalyzer.analyze(population);
    expect(narrative.selectivityInsight).toContain('Physical affinity unavailable');
    expect(narrative.modeDescriptions[0].characterization).not.toContain('kcal/mol');
  });

  it('stays proxy-only when Python forges receipts but declares no availability', async () => {
    const record = normalizeDockingRecord(pythonPayload({
      provenance: {
        schema_version: 2,
        energy_domain: 'calibrated_kcal_per_mol',
        ensemble_measure: 'enumerated_microstates',
        reference_state: 'matched_association_cycle',
        energy_provenance: ENERGY_SHA,
        measure_provenance: MEASURE_SHA,
        reference_provenance: REFERENCE_SHA,
      },
    }));
    expect(record.globalThermodynamics.available).toBe(false);
    expect(record.globalThermodynamics.unavailableReason)
      .toBe('availability_not_declared_by_producer');

    const analysis = await IntelligenceEngine.analyze(toPopulation(record));
    expect(analysis.claimValidity).toBe('proxy_only');
    expect(analysis.bullets.join('\n')).not.toContain('kcal/mol');
  });

  it('speaks physically only when receipts and available:true both arrive', async () => {
    const record = normalizeDockingRecord(pythonPayload({
      available: true,
      provenance: {
        schema_version: 2,
        energy_domain: 'calibrated_kcal_per_mol',
        ensemble_measure: 'enumerated_microstates',
        reference_state: 'matched_association_cycle',
        energy_provenance: ENERGY_SHA,
        measure_provenance: MEASURE_SHA,
        reference_provenance: REFERENCE_SHA,
      },
    }));
    expect(record.globalThermodynamics.available).toBe(true);

    const population = toPopulation(record);
    expect(populationAllowsPhysicalAffinity(population)).toBe(true);

    const analysis = await IntelligenceEngine.analyze(population);
    expect(analysis.claimValidity).toBe('binding_physical');
    expect(analysis.bullets[0]).toContain('Strong binding affinity');
    expect(analysis.bullets[0]).toContain('kcal/mol');
  });

  it('fails closed for every hostile availability encoding', async () => {
    for (const available of HOSTILE_AVAILABILITY_VALUES) {
      const record = normalizeDockingRecord(pythonPayload({
        available,
        provenance: {
          schema_version: 2,
          energy_domain: 'calibrated_kcal_per_mol',
          ensemble_measure: 'enumerated_microstates',
          reference_state: 'matched_association_cycle',
          energy_provenance: ENERGY_SHA,
          measure_provenance: MEASURE_SHA,
          reference_provenance: REFERENCE_SHA,
        },
      }));
      expect(record.globalThermodynamics.available).toBe(false);
      expect(claimValidityForRecord(record.globalThermodynamics)).toBe('proxy_only');

      const population = toPopulation(record);
      expect(populationAllowsPhysicalAffinity(population)).toBe(false);
      const analysis = await IntelligenceEngine.analyze(population);
      expect(analysis.claimValidity).toBe('proxy_only');
      expect(analysis.bullets.join('\n')).not.toContain('kcal/mol');
    }
  });

  it('gives SDK camelCase, viewer `modes` and Python snake_case the same verdict', () => {
    const physical: ThermodynamicClaimSource & Record<string, unknown> = {
      available: true,
      temperature: 298.15,
      logZ: 3.25,
      freeEnergy: -12.4,
      meanEnergy: -13.2,
      meanEnergySq: 180,
      heatCapacity: 0.021,
      entropy: 0.02,
      stdEnergy: 1.1,
      scientificProvenance: BINDING_PHYSICAL_PROVENANCE,
    };
    const mode = {
      index: 0, size: 12, freeEnergy: -12.4, entropy: 0.02,
      enthalpy: -13.2, heatCapacity: 0.021, probability: 1,
      thermodynamics: physical,
    };

    const verdicts = [
      normalizeDockingRecord({ bindingModes: [mode], globalThermodynamics: physical, temperature: 298.15 }),
      normalizeDockingRecord({ modes: [mode], globalThermodynamics: physical, temperature: 298.15 }),
      normalizeDockingRecord({ binding_modes: [mode], global_thermodynamics: physical, temperature_K: 298.15 }),
    ].map((record) => claimValidityForRecord(record.globalThermodynamics));

    expect(verdicts).toEqual(['binding_physical', 'binding_physical', 'binding_physical']);
  });

  it('keeps a proxy selectivity driver inconclusive end to end', () => {
    const analysis = SelectivityAnalyzer.analyze({
      ligandName: 'LSD',
      targets: [
        {
          targetName: '5HT2A', bestFreeEnergy: -9.5, modeCount: 3, sConf: 1.8,
          sVib: 0.002, isConverged: true, cavityVolume: 450, populationSize: 500,
          scientificProvenance: BINDING_PHYSICAL_PROVENANCE,
        },
        {
          targetName: 'D2R', bestFreeEnergy: -6.2, modeCount: 3, sConf: 1.2,
          sVib: 0.002, isConverged: true, cavityVolume: 450, populationSize: 500,
          scientificProvenance: BINDING_PHYSICAL_PROVENANCE,
        },
      ],
      deltaDeltas: [{ targetA: '5HT2A', targetB: 'D2R', ddg: -3.3 }],
    });

    // Receipts are valid, but neither target declared availability.
    expect(analysis.claimValidity).toBe('proxy_only');
    expect(analysis.driver).toBe('inconclusive');
    expect(analysis.explanation).not.toContain('kcal/mol');
    expect(analysis.deltaG).toBeCloseTo(-3.3, 5);
  });
});
