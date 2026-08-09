// wireNormalization.ts — the single runtime trust boundary for docking payloads
//
// Producers disagree on spelling: the TypeScript SDK emits camelCase
// (`bindingModes`), the shared/viewer model emits `modes`, and the Python
// package emits snake_case with unit suffixes (`binding_modes`,
// `heat_capacity_kcal_mol_K`, `G_config_kcal_mol`, ...). Before this module
// there was no single place where those shapes met, so each consumer invented
// its own defaults — including defaulting availability to "present".
//
// Everything that crosses into TypeScript from a file, a socket, or another
// language goes through here. Two invariants:
//   1. Numbers are re-keyed, never recomputed. A missing moment becomes 0 and
//      forces the record unavailable; it is never synthesised.
//   2. Availability is always set deliberately. It is `true` only when the
//      producer wrote a literal boolean `true`; otherwise it is `false` with a
//      machine-readable reason.
//
// Verified against python/flexaidds/thermodynamics.py (Thermodynamics.to_dict,
// ThermodynamicBreakdown.to_dict) and python/flexaidds/models.py
// (DockingResult.to_json / _binding_mode_json_record).
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import {
  PROXY_ONLY_PROVENANCE,
  asRecord,
  hasStrictAvailability,
  normalizeScientificProvenance,
  type ClaimValidity,
  type ScientificProvenance,
} from './scientificProvenance.js';
import { deriveClaimValidity } from './scientificProvenance.js';

type UnknownRecord = Record<string, unknown>;

/** Reason codes attached to a record that cannot be presented physically. */
export type UnavailableReason =
  | 'record_is_not_an_object'
  | 'availability_not_declared_by_producer'
  | 'availability_declared_non_true'
  | 'thermodynamic_moments_missing_or_non_numeric';

/** Canonical thermodynamic record. Availability is always explicit. */
export interface NormalizedThermodynamicRecord {
  available: boolean;
  unavailableReason?: UnavailableReason;
  scientificProvenance: ScientificProvenance;
  temperature: number;
  logZ: number;
  freeEnergy: number;
  meanEnergy: number;
  meanEnergySq: number;
  heatCapacity: number;
  entropy: number;
  stdEnergy: number;
}

/** Canonical binding-mode record. */
export interface NormalizedBindingMode {
  index: number;
  size: number;
  freeEnergy: number;
  entropy: number;
  enthalpy: number;
  heatCapacity: number;
  probability: number;
  thermodynamics: NormalizedThermodynamicRecord;
}

/** Canonical docking payload, whatever dialect the producer used. */
export interface NormalizedDockingRecord {
  modes: NormalizedBindingMode[];
  globalThermodynamics: NormalizedThermodynamicRecord;
  temperature: number;
  /** False when no producer key declared a temperature. */
  temperatureKnown: boolean;
  totalPoses: number;
}

/** Only finite JSON numbers count; strings are not silently coerced. */
function finiteNumber(value: unknown): number | undefined {
  return typeof value === 'number' && Number.isFinite(value) ? value : undefined;
}

/** First key present with a finite numeric value, in declaration order. */
function pickNumber(raw: UnknownRecord, keys: readonly string[]): number | undefined {
  for (const key of keys) {
    if (!(key in raw)) continue;
    const value = finiteNumber(raw[key]);
    if (value !== undefined) return value;
  }
  return undefined;
}

/** Availability keys understood on the wire, most specific first. */
const AVAILABILITY_KEYS = [
  'available',
  'thermodynamicsAvailable',
  'thermodynamics_available',
] as const;

/**
 * Decide availability for one untrusted record.
 *
 * Absent, null, false, 0, 1, "true", "false", [] and {} all fail closed. The
 * distinction between "never declared" and "declared not-true" is preserved so
 * a UI can say which one happened.
 */
export function normalizeAvailability(
  value: unknown,
): { available: boolean; unavailableReason?: UnavailableReason } {
  const raw = asRecord(value);
  if (!raw) return { available: false, unavailableReason: 'record_is_not_an_object' };

  for (const key of AVAILABILITY_KEYS) {
    if (!(key in raw)) continue;
    if (hasStrictAvailability(raw[key])) return { available: true };
    return { available: false, unavailableReason: 'availability_declared_non_true' };
  }
  return { available: false, unavailableReason: 'availability_not_declared_by_producer' };
}

// Accepted spellings per moment. camelCase (SDK/shared) first, then the
// Python unit-suffixed JSON keys, then Python raw attribute names, then the
// ThermodynamicBreakdown ledger names.
const TEMPERATURE_KEYS = ['temperature', 'temperature_K', 'temperatureK'] as const;
const LOG_Z_KEYS = ['logZ', 'log_Z', 'logZ_config'] as const;
const FREE_ENERGY_KEYS = [
  'freeEnergy', 'free_energy_kcal_mol', 'free_energy', 'G_config_kcal_mol',
] as const;
const MEAN_ENERGY_KEYS = [
  'meanEnergy', 'enthalpy_kcal_mol', 'mean_energy', 'H_eff_kcal_mol',
] as const;
const MEAN_ENERGY_SQ_KEYS = ['meanEnergySq', 'mean_energy_sq'] as const;
const HEAT_CAPACITY_KEYS = [
  'heatCapacity', 'heat_capacity_kcal_mol_K', 'heat_capacity_kcal_mol_K2',
  'heat_capacity', 'Cv_kcal_mol_K',
] as const;
const ENTROPY_KEYS = [
  'entropy', 'entropy_kcal_mol_K', 'S_config_kcal_mol_K',
] as const;
const STD_ENERGY_KEYS = [
  'stdEnergy', 'std_energy_kcal_mol', 'std_energy', 'sigma_E_kcal_mol',
] as const;

// Moments a record must carry before "available: true" can mean anything.
const REQUIRED_MOMENTS: ReadonlyArray<readonly string[]> = [
  TEMPERATURE_KEYS, FREE_ENERGY_KEYS, MEAN_ENERGY_KEYS, ENTROPY_KEYS,
];

/** Combine availability signals: every declaration must be strictly true. */
function intersectAvailability(
  ...signals: Array<{ available: boolean; unavailableReason?: UnavailableReason } | undefined>
): { available: boolean; unavailableReason?: UnavailableReason } {
  const declared = signals.filter((s): s is NonNullable<typeof s> => s !== undefined);
  if (declared.length === 0) {
    return { available: false, unavailableReason: 'availability_not_declared_by_producer' };
  }
  const blocked = declared.find((s) => !s.available);
  return blocked ?? { available: true };
}

/** Canonicalize one thermodynamic record from any supported dialect. */
export function normalizeThermodynamicRecord(
  value: unknown,
  availabilityOverride?: { available: boolean; unavailableReason?: UnavailableReason },
): NormalizedThermodynamicRecord {
  const raw = asRecord(value);
  if (!raw) {
    return {
      available: false,
      unavailableReason: availabilityOverride?.available === false
        ? availabilityOverride.unavailableReason
        : 'record_is_not_an_object',
      scientificProvenance: { ...PROXY_ONLY_PROVENANCE },
      temperature: 0, logZ: 0, freeEnergy: 0,
      meanEnergy: 0, meanEnergySq: 0,
      heatCapacity: 0, entropy: 0, stdEnergy: 0,
    };
  }

  const availability = availabilityOverride ?? normalizeAvailability(raw);
  const momentsComplete = REQUIRED_MOMENTS.every(
    (keys) => pickNumber(raw, keys) !== undefined,
  );

  return {
    available: availability.available && momentsComplete,
    unavailableReason: availability.available && !momentsComplete
      ? 'thermodynamic_moments_missing_or_non_numeric'
      : availability.unavailableReason,
    scientificProvenance: normalizeScientificProvenance(
      raw.scientificProvenance ?? raw.scientific_provenance,
    ),
    temperature: pickNumber(raw, TEMPERATURE_KEYS) ?? 0,
    logZ: pickNumber(raw, LOG_Z_KEYS) ?? 0,
    freeEnergy: pickNumber(raw, FREE_ENERGY_KEYS) ?? 0,
    meanEnergy: pickNumber(raw, MEAN_ENERGY_KEYS) ?? 0,
    meanEnergySq: pickNumber(raw, MEAN_ENERGY_SQ_KEYS) ?? 0,
    heatCapacity: pickNumber(raw, HEAT_CAPACITY_KEYS) ?? 0,
    entropy: pickNumber(raw, ENTROPY_KEYS) ?? 0,
    stdEnergy: pickNumber(raw, STD_ENERGY_KEYS) ?? 0,
  };
}

const MODE_INDEX_KEYS = ['index', 'modeId', 'mode_id', 'rank'] as const;
const MODE_SIZE_KEYS = ['size', 'nPoses', 'n_poses'] as const;
const MODE_FREE_ENERGY_KEYS = [
  'freeEnergy', 'free_energy', 'proxyFreeEnergy', 'proxy_free_energy',
] as const;
const MODE_ENTHALPY_KEYS = ['enthalpy', 'meanEnergy', 'mean_energy'] as const;

/** Canonicalize one binding mode from any supported dialect. */
function normalizeBindingMode(value: unknown, fallbackIndex: number): NormalizedBindingMode {
  const raw = asRecord(value) ?? {};

  // A mode may declare availability on itself (SDK `thermodynamicsAvailable`)
  // or only inside its nested record (shared/Python `thermodynamics`). Both
  // declarations must agree on "true"; neither spelling can promote the other.
  const nested = asRecord(raw.thermodynamics);
  const outerDeclares = AVAILABILITY_KEYS.some((key) => key in raw);
  const nestedDeclares = nested !== undefined
    && AVAILABILITY_KEYS.some((key) => key in nested);
  const availability = intersectAvailability(
    outerDeclares ? normalizeAvailability(raw) : undefined,
    nestedDeclares ? normalizeAvailability(nested) : undefined,
  );
  const thermodynamics = normalizeThermodynamicRecord(
    nested ?? (outerDeclares ? raw : undefined),
    availability,
  );

  return {
    index: pickNumber(raw, MODE_INDEX_KEYS) ?? fallbackIndex,
    size: pickNumber(raw, MODE_SIZE_KEYS) ?? 0,
    freeEnergy: pickNumber(raw, MODE_FREE_ENERGY_KEYS) ?? 0,
    entropy: pickNumber(raw, ENTROPY_KEYS) ?? 0,
    enthalpy: pickNumber(raw, MODE_ENTHALPY_KEYS) ?? 0,
    heatCapacity: pickNumber(raw, HEAT_CAPACITY_KEYS) ?? 0,
    probability: pickNumber(raw, ['probability']) ?? 0,
    thermodynamics,
  };
}

const MODE_LIST_KEYS = ['bindingModes', 'modes', 'binding_modes'] as const;
const GLOBAL_THERMO_KEYS = [
  'globalThermodynamics', 'global_thermodynamics', 'thermodynamics',
] as const;
const TOTAL_POSES_KEYS = ['totalPoses', 'total_poses', 'populationSize', 'population_size'] as const;

/**
 * Canonicalize a whole docking payload: SDK camelCase, shared/viewer `modes`,
 * or Python snake_case. Mode order is preserved exactly as produced.
 */
export function normalizeDockingRecord(value: unknown): NormalizedDockingRecord {
  const raw = asRecord(value) ?? {};

  let rawModes: unknown[] = [];
  for (const key of MODE_LIST_KEYS) {
    if (Array.isArray(raw[key])) {
      rawModes = raw[key] as unknown[];
      break;
    }
  }
  const modes = rawModes.map(normalizeBindingMode);

  let rawGlobal: unknown;
  for (const key of GLOBAL_THERMO_KEYS) {
    if (key in raw) {
      rawGlobal = raw[key];
      break;
    }
  }
  const globalThermodynamics = normalizeThermodynamicRecord(rawGlobal);

  const temperature = pickNumber(raw, TEMPERATURE_KEYS);
  const totalPoses = pickNumber(raw, TOTAL_POSES_KEYS)
    ?? modes.reduce((sum, mode) => sum + mode.size, 0);

  return {
    modes,
    globalThermodynamics,
    temperature: temperature ?? globalThermodynamics.temperature,
    temperatureKnown: temperature !== undefined,
    totalPoses,
  };
}

/**
 * Claim level for an untrusted wire record, normalizing first.
 *
 * Both the shared model and the SDK route their record-level gating through
 * this function so a snake_case payload can never out-rank a camelCase one.
 */
export function claimValidityForWireRecord(value: unknown): ClaimValidity {
  const normalized = normalizeThermodynamicRecord(value);
  if (!normalized.available) return 'proxy_only';
  return deriveClaimValidity(normalized.scientificProvenance);
}
