// Scientific provenance contract tests for the shared viewer model.
// SPDX-License-Identifier: Apache-2.0

import assert from 'node:assert/strict';
import test from 'node:test';

import {
  PROXY_ONLY_PROVENANCE,
  allowsBindingClaims,
  allowsCanonicalClaims,
  claimValidityForRecord,
  claimValidityForWireRecord,
  deriveClaimValidity,
  deserializePopulation,
  hasStrictAvailability,
  normalizeAvailability,
  normalizeDockingRecord,
  normalizeScientificProvenance,
  normalizeThermodynamicRecord,
  serializePopulation,
} from './dist/index.js';

const ENERGY_SHA = `sha256:${'4692da7b40da99a82a86a6c30e33e4bedead9a2dbcc4b28d977e675fd0761993'}`;
const MEASURE_SHA = `sha256:${'aa1c1b96e5831c2c2d8ffe1060b43301e622c4255bc5e4f08765243e55265353'}`;
const REFERENCE_SHA = `sha256:${'4d28b0f5a589c9d228295118cbf17d810b54fca40a2cdb6159cec35788971050'}`;

const bindingPhysical = {
  schemaVersion: 2,
  energyDomain: 'calibrated_kcal_per_mol',
  ensembleMeasure: 'enumerated_microstates',
  referenceState: 'matched_association_cycle',
  energyProvenance: ENERGY_SHA,
  measureProvenance: MEASURE_SHA,
  referenceProvenance: REFERENCE_SHA,
};

function populationPayload(scientificProvenance) {
  return {
    modes: [
      {
        index: 0,
        size: 3,
        freeEnergy: -8,
        entropy: 0.1,
        enthalpy: -7,
        heatCapacity: 0.02,
        probability: 0.7,
      },
      {
        index: 1,
        size: 2,
        freeEnergy: -5,
        entropy: 0.2,
        enthalpy: -4,
        heatCapacity: 0.03,
        probability: 0.3,
      },
    ],
    globalThermodynamics: {
      temperature: 298.15,
      logZ: 1,
      freeEnergy: -8,
      meanEnergy: -7,
      meanEnergySq: 50,
      heatCapacity: 0.02,
      entropy: 0.1,
      stdEnergy: 1,
      ...(scientificProvenance === undefined ? {} : { scientificProvenance }),
    },
    temperature: 298.15,
    totalPoses: 5,
    shannonS: 0.4,
    isCollapsed: false,
  };
}

test('derives only evidence-backed schema-v2 claim levels', () => {
  assert.equal(deriveClaimValidity(), 'proxy_only');
  assert.equal(deriveClaimValidity(bindingPhysical), 'binding_physical');
  assert.equal(
    deriveClaimValidity({ ...bindingPhysical, referenceState: 'bound_only' }),
    'canonical_physical',
  );
  assert.equal(
    deriveClaimValidity({
      ...PROXY_ONLY_PROVENANCE,
      claimValidity: 'binding_physical',
    }),
    'proxy_only',
  );
});

test('rejects malformed schemas, prose receipts, filler, and low-diversity digests', () => {
  for (const schemaVersion of [1, 2.9, '2', true]) {
    assert.equal(
      deriveClaimValidity({ ...bindingPhysical, schemaVersion }),
      'proxy_only',
    );
  }

  const invalidEvidence = [
    'receipt',
    'sha256:abc123',
    `SHA256:${ENERGY_SHA.slice('sha256:'.length)}`,
    `sha256:${'0'.repeat(64)}`,
    `sha256:${'ab'.repeat(32)}`,
    'sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855',
    'sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a',
    '\u00a0',
    '\u2003',
    '参照',
    123,
    ['sha256', ENERGY_SHA],
    true,
  ];
  for (const evidence of invalidEvidence) {
    assert.equal(
      deriveClaimValidity({ ...bindingPhysical, energyProvenance: evidence }),
      'proxy_only',
    );
  }
});

test('accepts upper-case hex but fails closed when moments are unavailable', () => {
  const provenance = {
    ...bindingPhysical,
    energyProvenance: ENERGY_SHA.toUpperCase().replace('SHA256:', 'sha256:'),
  };
  assert.equal(deriveClaimValidity(provenance), 'binding_physical');

  // Availability is a separate gate: provenance alone never authorizes a
  // physical presentation.
  assert.equal(allowsCanonicalClaims({ scientificProvenance: provenance }), false);
  assert.equal(allowsBindingClaims({ scientificProvenance: provenance }), false);
  assert.equal(
    allowsCanonicalClaims({ scientificProvenance: provenance, available: true }),
    true,
  );
  assert.equal(
    allowsBindingClaims({ scientificProvenance: provenance, available: true }),
    true,
  );
});

test('only a literal boolean true satisfies the availability gate', () => {
  const provenance = { ...bindingPhysical };
  const hostile = [
    undefined, null, false, 0, 1, 'true', 'false', '', [], {}, [true], { value: true },
    Number.NaN, 'TRUE', 1n,
  ];
  for (const available of hostile) {
    assert.equal(hasStrictAvailability(available), false);
    assert.equal(
      allowsCanonicalClaims({ scientificProvenance: provenance, available }),
      false,
      `availability ${String(available)} must fail closed`,
    );
    assert.equal(
      allowsBindingClaims({ scientificProvenance: provenance, available }),
      false,
    );
    assert.equal(
      claimValidityForRecord({ scientificProvenance: provenance, available }),
      'proxy_only',
    );
  }
  assert.equal(hasStrictAvailability(true), true);
  assert.equal(
    claimValidityForRecord({ scientificProvenance: provenance, available: true }),
    'binding_physical',
  );
  assert.equal(claimValidityForRecord(undefined), 'proxy_only');
  assert.equal(claimValidityForRecord(null), 'proxy_only');
});

test('normalizes snake_case metadata and ignores serialized validity', () => {
  const normalized = normalizeScientificProvenance({
    schema_version: 2,
    energy_domain: 'calibrated_kcal_per_mol',
    ensemble_measure: 'weighted_quadrature',
    reference_state: 'matched_association_cycle',
    energy_provenance: ENERGY_SHA,
    measure_provenance: MEASURE_SHA,
    reference_provenance: REFERENCE_SHA,
    claim_validity: 'proxy_only',
  });
  assert.equal(deriveClaimValidity(normalized), 'binding_physical');
});

test('legacy deserialization inserts proxy provenance without changing values or order', () => {
  const payload = populationPayload(undefined);
  const restored = deserializePopulation(JSON.stringify(payload));

  assert.deepEqual(restored.modes.map((mode) => mode.index), [0, 1]);
  assert.deepEqual(restored.modes.map((mode) => mode.freeEnergy), [-8, -5]);
  assert.equal(restored.globalThermodynamics.freeEnergy, -8);
  assert.deepEqual(
    restored.globalThermodynamics.scientificProvenance,
    PROXY_ONLY_PROVENANCE,
  );
  assert.equal(allowsBindingClaims(restored.globalThermodynamics), false);

  const serialized = JSON.parse(serializePopulation(payload));
  assert.deepEqual(
    serialized.globalThermodynamics.scientificProvenance,
    PROXY_ONLY_PROVENANCE,
  );
});

// ─── Wire normalization goldens ─────────────────────────────────────────────
//
// Key names below were read out of python/flexaidds/models.py
// (DockingResult.to_json / _binding_mode_json_record) and
// python/flexaidds/thermodynamics.py (Thermodynamics.to_dict,
// ThermodynamicBreakdown.to_dict, ScientificProvenance.to_dict). If Python
// renames a field these goldens must fail.

function pythonProvenance(overrides = {}) {
  return {
    schema_version: 2,
    energy_domain: 'cf_arbitrary_units',
    ensemble_measure: 'optimizer_samples',
    reference_state: 'none',
    energy_provenance: 'FlexAID CF contact-function score',
    measure_provenance: 'GA optimizer samples',
    reference_provenance: '',
    // Python serializes a derived label; it must never be trusted.
    claim_validity: 'binding_physical',
    ...overrides,
  };
}

function pythonDockingPayload(provenance = pythonProvenance()) {
  return {
    source_dir: '/tmp/results',
    temperature: 298.15,
    n_modes: 2,
    metadata: {},
    binding_modes: [
      {
        mode_id: 0,
        rank: 1,
        n_poses: 12,
        free_energy: -8.4,
        proxy_free_energy: -8.4,
        soft_beta_G: -8.1,
        enthalpy: -9.2,
        entropy: 0.0031,
        heat_capacity: 0.021,
        std_energy: 1.1,
        best_cf: -142.7,
        temperature: 298.15,
        scientific_provenance: provenance,
        best_pose_path: '/tmp/results/mode0_pose1.pdb',
        thermodynamics: {
          temperature_K: 298.15,
          logZ_config: 3.25,
          G_config_kcal_mol: -8.4,
          H_eff_kcal_mol: -9.2,
          S_config_kcal_mol_K: 0.0031,
          minus_T_S_config_kcal_mol: 0.8,
          Cv_kcal_mol_K: 0.021,
          sigma_E_kcal_mol: 1.1,
          G_vib_kcal_mol: 0.0,
          G_natural_kcal_mol: 0.0,
          G_other_kcal_mol: 0.0,
          G_total_kcal_mol: -8.4,
          component_sum_kcal_mol: 0.0,
          components_complete: false,
          component_means: {},
          affinity: null,
          scientific_provenance: provenance,
        },
      },
      {
        mode_id: 1,
        rank: 2,
        n_poses: 5,
        free_energy: -6.1,
        enthalpy: -6.8,
        entropy: 0.0022,
        heat_capacity: 0.018,
        std_energy: 0.9,
        temperature: 298.15,
        scientific_provenance: provenance,
      },
    ],
  };
}

test('normalizes a real Python docking payload without touching numbers', () => {
  const normalized = normalizeDockingRecord(pythonDockingPayload());

  assert.equal(normalized.temperature, 298.15);
  assert.equal(normalized.temperatureKnown, true);
  assert.deepEqual(normalized.modes.map((m) => m.index), [0, 1]);
  assert.deepEqual(normalized.modes.map((m) => m.size), [12, 5]);
  assert.deepEqual(normalized.modes.map((m) => m.freeEnergy), [-8.4, -6.1]);
  assert.deepEqual(normalized.modes.map((m) => m.enthalpy), [-9.2, -6.8]);
  assert.deepEqual(normalized.modes.map((m) => m.entropy), [0.0031, 0.0022]);
  assert.deepEqual(normalized.modes.map((m) => m.heatCapacity), [0.021, 0.018]);
  assert.equal(normalized.totalPoses, 17);

  // Unit-suffixed ledger keys map onto the canonical slots.
  const first = normalized.modes[0].thermodynamics;
  assert.equal(first.temperature, 298.15);
  assert.equal(first.logZ, 3.25);
  assert.equal(first.freeEnergy, -8.4);
  assert.equal(first.meanEnergy, -9.2);
  assert.equal(first.entropy, 0.0031);
  assert.equal(first.heatCapacity, 0.021);
  assert.equal(first.stdEnergy, 1.1);

  // Python declares no availability at all -> deliberate false, never true.
  for (const mode of normalized.modes) {
    assert.equal(mode.thermodynamics.available, false);
    assert.equal(
      mode.thermodynamics.unavailableReason,
      'availability_not_declared_by_producer',
    );
    assert.equal(claimValidityForRecord(mode.thermodynamics), 'proxy_only');
  }
  // No global thermodynamics in the Python payload at all.
  assert.equal(normalized.globalThermodynamics.available, false);
});

test('ignores a serialized Python claim_validity even with real receipts', () => {
  const forged = pythonProvenance({
    energy_domain: 'calibrated_kcal_per_mol',
    ensemble_measure: 'enumerated_microstates',
    reference_state: 'matched_association_cycle',
    energy_provenance: ENERGY_SHA,
    measure_provenance: MEASURE_SHA,
    reference_provenance: REFERENCE_SHA,
  });
  const normalized = normalizeDockingRecord(pythonDockingPayload(forged));

  // Receipts are real, so provenance derives binding_physical — but the
  // producer never declared availability, so the record stays proxy-only.
  assert.equal(
    deriveClaimValidity(normalized.modes[0].thermodynamics.scientificProvenance),
    'binding_physical',
  );
  assert.equal(normalized.modes[0].thermodynamics.available, false);
  assert.equal(claimValidityForRecord(normalized.modes[0].thermodynamics), 'proxy_only');
});

test('accepts all three mode-list dialects with identical results', () => {
  const thermo = {
    available: true,
    temperature: 300,
    logZ: 1,
    freeEnergy: -8,
    meanEnergy: -7,
    meanEnergySq: 50,
    heatCapacity: 0.02,
    entropy: 0.1,
    stdEnergy: 1,
    scientificProvenance: bindingPhysical,
  };
  const mode = { index: 0, size: 3, freeEnergy: -8, entropy: 0.1, enthalpy: -7, heatCapacity: 0.02, probability: 1, thermodynamics: thermo };

  const sdk = normalizeDockingRecord({ bindingModes: [mode], globalThermodynamics: thermo, temperature: 300 });
  const viewer = normalizeDockingRecord({ modes: [mode], globalThermodynamics: thermo, temperature: 300 });
  const python = normalizeDockingRecord({ binding_modes: [mode], global_thermodynamics: thermo, temperature_K: 300 });

  assert.deepEqual(sdk, viewer);
  assert.deepEqual(sdk, python);
  assert.equal(sdk.globalThermodynamics.available, true);
  assert.equal(claimValidityForRecord(sdk.globalThermodynamics), 'binding_physical');
});

test('a nested false availability cannot be promoted by the outer mode', () => {
  const normalized = normalizeDockingRecord({
    bindingModes: [{
      index: 0,
      size: 1,
      thermodynamicsAvailable: true,
      thermodynamics: {
        available: false,
        temperature: 300,
        freeEnergy: -8,
        meanEnergy: -7,
        entropy: 0.1,
        scientificProvenance: bindingPhysical,
      },
    }],
  });
  assert.equal(normalized.modes[0].thermodynamics.available, false);
  assert.equal(claimValidityForRecord(normalized.modes[0].thermodynamics), 'proxy_only');
});

test('available:true without complete moments still fails closed', () => {
  const record = normalizeThermodynamicRecord({
    available: true,
    scientificProvenance: bindingPhysical,
    freeEnergy: -8,
  });
  assert.equal(record.available, false);
  assert.equal(record.unavailableReason, 'thermodynamic_moments_missing_or_non_numeric');
  assert.equal(record.temperature, 0);
});

test('numeric-looking strings are not coerced into moments', () => {
  const record = normalizeThermodynamicRecord({
    available: true,
    temperature: '300',
    freeEnergy: '-8',
    meanEnergy: -7,
    entropy: 0.1,
    scientificProvenance: bindingPhysical,
  });
  assert.equal(record.available, false);
  assert.equal(record.temperature, 0);
  assert.equal(record.freeEnergy, 0);
});

test('every non-true availability encoding is reported, never assumed', () => {
  assert.deepEqual(normalizeAvailability(undefined), {
    available: false,
    unavailableReason: 'record_is_not_an_object',
  });
  assert.deepEqual(normalizeAvailability({}), {
    available: false,
    unavailableReason: 'availability_not_declared_by_producer',
  });
  for (const available of [null, false, 0, 1, 'true', 'false', '', [], {}]) {
    assert.deepEqual(
      normalizeAvailability({ available }),
      { available: false, unavailableReason: 'availability_declared_non_true' },
      `availability ${String(available)} must fail closed`,
    );
  }
  assert.deepEqual(normalizeAvailability({ available: true }), { available: true });
});

test('claimValidityForWireRecord gates snake_case exactly like camelCase', () => {
  const camel = {
    available: true,
    temperature: 300, logZ: 0, freeEnergy: -8, meanEnergy: -7,
    meanEnergySq: 50, heatCapacity: 0, entropy: 0.1, stdEnergy: 0,
    scientificProvenance: bindingPhysical,
  };
  const snake = {
    available: true,
    temperature_K: 300, log_Z: 0, free_energy_kcal_mol: -8,
    enthalpy_kcal_mol: -7, mean_energy_sq: 50,
    heat_capacity_kcal_mol_K: 0, entropy_kcal_mol_K: 0.1,
    std_energy_kcal_mol: 0,
    scientific_provenance: {
      schema_version: 2,
      energy_domain: 'calibrated_kcal_per_mol',
      ensemble_measure: 'enumerated_microstates',
      reference_state: 'matched_association_cycle',
      energy_provenance: ENERGY_SHA,
      measure_provenance: MEASURE_SHA,
      reference_provenance: REFERENCE_SHA,
    },
  };
  assert.equal(claimValidityForWireRecord(camel), 'binding_physical');
  assert.equal(claimValidityForWireRecord(snake), 'binding_physical');
  assert.equal(claimValidityForWireRecord({ ...snake, available: 'true' }), 'proxy_only');
});

test('deserializePopulation sets availability deliberately', () => {
  const restored = deserializePopulation(JSON.stringify(populationPayload(bindingPhysical)));
  assert.equal(restored.globalThermodynamics.available, false);
  assert.equal(
    restored.globalThermodynamics.unavailableReason,
    'availability_not_declared_by_producer',
  );
  assert.equal(allowsBindingClaims(restored.globalThermodynamics), false);
  // Numbers and order are untouched.
  assert.equal(restored.globalThermodynamics.freeEnergy, -8);
  assert.deepEqual(restored.modes.map((mode) => mode.index), [0, 1]);
});

// Verbatim output of python/flexaidds DockingResult.to_json() (captured from
// the real dataclasses, not hand-written), including the `null` optional
// fields Python emits. If these keys move, this golden fails.
const PYTHON_TO_JSON_GOLDEN = `{
 "source_dir": "/tmp/results",
 "temperature": 298.15,
 "n_modes": 1,
 "metadata": {},
 "binding_modes": [
  {
   "mode_id": 0,
   "rank": 1,
   "n_poses": 0,
   "free_energy": -12.4,
   "proxy_free_energy": null,
   "soft_beta_G": null,
   "enthalpy": -13.2,
   "entropy": 0.02,
   "heat_capacity": 0.021,
   "std_energy": 1.1,
   "best_cf": null,
   "temperature": 298.15,
   "scientific_provenance": {
    "schema_version": 2,
    "energy_domain": "unclassified",
    "ensemble_measure": "unclassified",
    "reference_state": "none",
    "energy_provenance": "",
    "measure_provenance": "",
    "reference_provenance": "",
    "claim_validity": "proxy_only"
   },
   "best_pose_path": null,
   "thermodynamics": {
    "temperature_K": 298.15,
    "logZ_config": 3.25,
    "G_config_kcal_mol": -12.4,
    "H_eff_kcal_mol": -13.2,
    "S_config_kcal_mol_K": 0.02,
    "minus_T_S_config_kcal_mol": 0.7999999999999989,
    "Cv_kcal_mol_K": 0.021,
    "sigma_E_kcal_mol": 1.1,
    "G_vib_kcal_mol": 0.0,
    "G_natural_kcal_mol": 0.0,
    "G_other_kcal_mol": 0.0,
    "G_total_kcal_mol": -12.4,
    "component_sum_kcal_mol": 0.0,
    "components_complete": false,
    "component_means": {},
    "affinity": null,
    "scientific_provenance": {
     "schema_version": 2,
     "energy_domain": "unclassified",
     "ensemble_measure": "unclassified",
     "reference_state": "none",
     "energy_provenance": "",
     "measure_provenance": "",
     "reference_provenance": "",
     "claim_validity": "proxy_only"
    }
   }
  }
 ],
 "grand_log_xi": null,
 "ligand_occupancies": {},
 "selectivities": {},
 "per_ligand_results": {},
 "empty_probability": null,
 "mean_occupancy": null
}`;

test('handles verbatim Python to_json() output including its nulls', () => {
  const normalized = normalizeDockingRecord(JSON.parse(PYTHON_TO_JSON_GOLDEN));

  assert.equal(normalized.temperature, 298.15);
  assert.equal(normalized.temperatureKnown, true);
  assert.equal(normalized.modes.length, 1);
  assert.equal(normalized.modes[0].index, 0);
  assert.equal(normalized.modes[0].size, 0);
  assert.equal(normalized.modes[0].freeEnergy, -12.4);
  assert.equal(normalized.modes[0].enthalpy, -13.2);
  assert.equal(normalized.modes[0].entropy, 0.02);
  assert.equal(normalized.modes[0].heatCapacity, 0.021);

  const ledger = normalized.modes[0].thermodynamics;
  assert.equal(ledger.temperature, 298.15);
  assert.equal(ledger.logZ, 3.25);
  assert.equal(ledger.freeEnergy, -12.4);
  assert.equal(ledger.meanEnergy, -13.2);
  assert.equal(ledger.entropy, 0.02);
  assert.equal(ledger.heatCapacity, 0.021);
  assert.equal(ledger.stdEnergy, 1.1);

  // `null` optional fields must never be read as numbers.
  assert.equal(typeof normalized.modes[0].probability, 'number');
  assert.equal(normalized.modes[0].probability, 0);

  // Python emits no availability and unclassified provenance: proxy-only.
  assert.equal(ledger.available, false);
  assert.equal(ledger.unavailableReason, 'availability_not_declared_by_producer');
  assert.equal(claimValidityForRecord(ledger), 'proxy_only');
  assert.deepEqual(ledger.scientificProvenance, PROXY_ONLY_PROVENANCE);
});
