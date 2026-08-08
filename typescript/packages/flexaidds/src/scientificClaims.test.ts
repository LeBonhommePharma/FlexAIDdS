import { describe, expect, it } from 'vitest';

import { parseDockingResultJSON, parseResultFile } from './resultLoader.js';
import { StatMechEngine } from './StatMechEngine.js';
import {
  allowsBindingClaims,
  allowsCanonicalClaims,
  deriveClaimValidity,
  normalizeScientificProvenance,
  type ScientificProvenance,
} from './types.js';

const calibratedCanonical: ScientificProvenance = {
  schemaVersion: 2,
  energyDomain: 'calibrated_kcal_per_mol',
  ensembleMeasure: 'enumerated_microstates',
  referenceState: 'bound_only',
  energyProvenance: 'sha256:6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b',
  measureProvenance: 'sha256:d4735e3a265e16eee03f59718b9b5d03019c07d8b6c51f90da3a666eec13ab35',
  referenceProvenance: '',
};

describe('scientific claim firewall', () => {
  it('fails closed when provenance is absent', () => {
    expect(deriveClaimValidity()).toBe('proxy_only');
    expect(allowsCanonicalClaims({})).toBe(false);
    expect(allowsBindingClaims({})).toBe(false);
  });

  it('fails closed on malformed, trivial, or non-string evidence', () => {
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: '\u00a0\u2003',
    })).toBe('proxy_only');
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: 2,
    } as unknown as ScientificProvenance)).toBe('proxy_only');
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: `sha256:${'0'.repeat(64)}`,
    })).toBe('proxy_only');
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: `sha256:${'ab'.repeat(32)}`,
    })).toBe('proxy_only');
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: 'sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855',
    })).toBe('proxy_only');
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: 'sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a',
    })).toBe('proxy_only');
    expect(deriveClaimValidity({
      ...calibratedCanonical,
      energyProvenance: 'SHA256:6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b',
    })).toBe('proxy_only');
  });

  it('normalizes snake-case metadata and ignores serialized validity', () => {
    const provenance = normalizeScientificProvenance({
      schema_version: 2,
      energy_domain: 'cf_arbitrary_units',
      ensemble_measure: 'optimizer_samples',
      reference_state: 'bound_only',
      energy_provenance: 'CF score',
      measure_provenance: 'GA records',
      reference_provenance: 'bound only',
      claim_validity: 'binding_physical',
    });
    expect(deriveClaimValidity(provenance)).toBe('proxy_only');

    const result = parseDockingResultJSON(JSON.stringify({
      bindingModes: [],
      globalThermodynamics: {
        temperature: 300, logZ: 0, freeEnergy: -10,
        meanEnergy: -9, meanEnergySq: 81, heatCapacity: 0,
        entropy: 0, stdEnergy: 0,
        scientific_provenance: {
          schema_version: 2,
          energy_domain: 'calibrated_kcal_per_mol',
          ensemble_measure: 'enumerated_microstates',
          reference_state: 'bound_only',
          energy_provenance: 'sha256:6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b',
          measure_provenance: 'sha256:d4735e3a265e16eee03f59718b9b5d03019c07d8b6c51f90da3a666eec13ab35',
          reference_provenance: '',
        },
      },
      temperature: 300, populationSize: 0, timestamp: 'fixture',
    }));
    expect(deriveClaimValidity(result.globalThermodynamics.scientificProvenance))
      .toBe('canonical_physical');
  });

  it('permits canonical but not binding claims without a matched cycle', () => {
    expect(deriveClaimValidity(calibratedCanonical)).toBe('canonical_physical');
  });

  it('requires a matched association reference for binding claims', () => {
    const provenance: ScientificProvenance = {
      ...calibratedCanonical,
      referenceState: 'matched_association_cycle',
      referenceProvenance: 'sha256:4e07408562bedb8b60ce05c1decfe3ad16b72230967de01f640b7e4729b49fce',
    };
    expect(deriveClaimValidity(provenance)).toBe('binding_physical');
  });

  it('propagates provenance without changing engine numerics', () => {
    const proxy = new StatMechEngine(300);
    const physical = new StatMechEngine(300, calibratedCanonical);
    for (const energy of [-10, -8]) {
      proxy.addSample(energy);
      physical.addSample(energy);
    }

    const proxyResult = proxy.compute();
    const physicalResult = physical.compute();
    expect(physicalResult.freeEnergy).toBe(proxyResult.freeEnergy);
    expect(physicalResult.meanEnergy).toBe(proxyResult.meanEnergy);
    expect(deriveClaimValidity(proxyResult.scientificProvenance))
      .toBe('proxy_only');
    expect(deriveClaimValidity(physicalResult.scientificProvenance))
      .toBe('canonical_physical');

    const emptyPhysical = new StatMechEngine(300, calibratedCanonical).compute();
    expect(emptyPhysical.available).toBe(false);
    expect(allowsCanonicalClaims(emptyPhysical)).toBe(false);
  });

  it('marks raw CF result files explicitly unavailable and proxy-only', () => {
    const result = parseResultFile('CLUSTER|1|CF|-10.50|POSES|2\n');
    expect(result.globalThermodynamics.available).toBe(false);
    expect(result.globalThermodynamics.unavailableReason)
      .toBe('raw_result_contains_cf_scores_only');
    expect(deriveClaimValidity(result.globalThermodynamics.scientificProvenance))
      .toBe('proxy_only');
    expect(result.bindingModes[0].scoreDomain).toBe('cf_arbitrary_units');
    expect(result.bindingModes[0].thermodynamicsAvailable).toBe(false);
    // 300 K is the historical hard-coded default, deliberately unchanged by
    // the claim firewall. temperatureKnown is the metadata that flags it.
    expect(result.temperature).toBe(300);
    expect(result.temperatureKnown).toBe(false);
    expect(result.globalThermodynamics.temperature).toBe(300);
  });

  it('keeps .rrd pose counting identical to the pre-firewall parser', () => {
    // Regression pin: Chunk 0 is metadata-only. `size` must still come from
    // counted POSE lines with the historical `|| 1` floor, and cluster header
    // fields must still be read from the historical offset.
    const counted = parseResultFile([
      'CLUSTER|1|CF|-10.50|POSES|2',
      'POSE|1|CF|-10.50',
      'POSE|2|CF|-10.00',
    ].join('\n'));
    expect(counted.bindingModes[0].size).toBe(2);
    expect(counted.bindingModes[0].freeEnergy).toBe(0);

    const headerOnly = parseResultFile('CLUSTER|1|CF|-10.50|POSES|0\n');
    expect(headerOnly.bindingModes[0].size).toBe(1);
    expect(headerOnly.globalThermodynamics.available).toBe(false);
    expect(headerOnly.globalThermodynamics.unavailableReason)
      .toBe('raw_result_contains_cf_scores_only');
  });
});

// ─── Wire dialects ──────────────────────────────────────────────────────────

/** Shape produced by python/flexaidds DockingResult.to_json(). */
function pythonPayload(overrides: {
  provenance?: Record<string, unknown>;
  available?: unknown;
} = {}): Record<string, unknown> {
  const provenance = overrides.provenance ?? {
    schema_version: 2,
    energy_domain: 'cf_arbitrary_units',
    ensemble_measure: 'optimizer_samples',
    reference_state: 'none',
    energy_provenance: 'FlexAID CF contact-function score',
    measure_provenance: 'GA optimizer samples',
    reference_provenance: '',
    claim_validity: 'binding_physical', // forged label; must be ignored
  };
  const thermodynamics: Record<string, unknown> = {
    temperature_K: 298.15,
    logZ_config: 3.25,
    G_config_kcal_mol: -12.4,
    H_eff_kcal_mol: -13.2,
    S_config_kcal_mol_K: 0.02,
    Cv_kcal_mol_K: 0.021,
    sigma_E_kcal_mol: 1.1,
    scientific_provenance: provenance,
  };
  if ('available' in overrides) thermodynamics.available = overrides.available;

  return {
    source_dir: '/tmp/results',
    temperature: 298.15,
    n_modes: 1,
    binding_modes: [{
      mode_id: 0,
      rank: 1,
      n_poses: 12,
      free_energy: -12.4,
      enthalpy: -13.2,
      entropy: 0.02,
      heat_capacity: 0.021,
      std_energy: 1.1,
      temperature: 298.15,
      scientific_provenance: provenance,
      thermodynamics,
    }],
    global_thermodynamics: thermodynamics,
  };
}

describe('result loader wire dialects', () => {
  it('adopts a Python snake_case payload without changing numbers', () => {
    const result = parseDockingResultJSON(JSON.stringify(pythonPayload()));

    expect(result.bindingModes).toHaveLength(1);
    expect(result.bindingModes[0].size).toBe(12);
    expect(result.bindingModes[0].freeEnergy).toBe(-12.4);
    expect(result.bindingModes[0].enthalpy).toBe(-13.2);
    expect(result.bindingModes[0].entropy).toBe(0.02);
    expect(result.bindingModes[0].heatCapacity).toBe(0.021);
    expect(result.temperature).toBe(298.15);
    expect(result.temperatureKnown).toBe(true);
    expect(result.globalThermodynamics.freeEnergy).toBe(-12.4);
    expect(result.globalThermodynamics.meanEnergy).toBe(-13.2);
    expect(result.globalThermodynamics.stdEnergy).toBe(1.1);

    // Python declares no availability, so the record fails closed.
    expect(result.globalThermodynamics.available).toBe(false);
    expect(result.bindingModes[0].thermodynamicsAvailable).toBe(false);
    expect(allowsCanonicalClaims(result.globalThermodynamics)).toBe(false);
  });

  it('never lets a forged claim_validity or availability open the gate', () => {
    const receipts = {
      schema_version: 2,
      energy_domain: 'calibrated_kcal_per_mol',
      ensemble_measure: 'enumerated_microstates',
      reference_state: 'matched_association_cycle',
      energy_provenance: 'sha256:4692da7b40da99a82a86a6c30e33e4bedead9a2dbcc4b28d977e675fd0761993',
      measure_provenance: 'sha256:aa1c1b96e5831c2c2d8ffe1060b43301e622c4255bc5e4f08765243e55265353',
      reference_provenance: 'sha256:4d28b0f5a589c9d228295118cbf17d810b54fca40a2cdb6159cec35788971050',
    };
    const hostile: readonly unknown[] = [
      undefined, null, false, 0, 1, 'true', 'false', '', [], {}, [true], { value: true },
    ];
    for (const available of hostile) {
      const result = parseDockingResultJSON(JSON.stringify(
        pythonPayload({ provenance: receipts, available }),
      ));
      expect(allowsBindingClaims(result.globalThermodynamics)).toBe(false);
      expect(allowsCanonicalClaims(result.globalThermodynamics)).toBe(false);
    }

    const honest = parseDockingResultJSON(JSON.stringify(
      pythonPayload({ provenance: receipts, available: true }),
    ));
    expect(allowsBindingClaims(honest.globalThermodynamics)).toBe(true);
  });

  it('reads SDK camelCase and viewer `modes` identically', () => {
    const thermo = {
      available: true,
      temperature: 300, logZ: 1, freeEnergy: -8, meanEnergy: -7,
      meanEnergySq: 50, heatCapacity: 0.02, entropy: 0.1, stdEnergy: 1,
      scientificProvenance: calibratedCanonical,
    };
    const mode = {
      index: 0, size: 3, freeEnergy: -8, entropy: 0.1,
      enthalpy: -7, heatCapacity: 0.02, thermodynamics: thermo,
    };

    const sdk = parseDockingResultJSON(JSON.stringify({
      bindingModes: [mode], globalThermodynamics: thermo, temperature: 300,
      populationSize: 3, timestamp: 'fixture',
    }));
    const viewer = parseDockingResultJSON(JSON.stringify({
      modes: [mode], globalThermodynamics: thermo, temperature: 300,
      populationSize: 3, timestamp: 'fixture',
    }));

    expect(viewer.bindingModes).toEqual(sdk.bindingModes);
    expect(viewer.globalThermodynamics).toEqual(sdk.globalThermodynamics);
    expect(allowsCanonicalClaims(sdk.globalThermodynamics)).toBe(true);
    expect(allowsBindingClaims(sdk.globalThermodynamics)).toBe(false);
    expect(sdk.bindingModes[0].freeEnergy).toBe(-8);
  });
});
