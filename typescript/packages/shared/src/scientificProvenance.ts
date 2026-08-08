// scientificProvenance.ts — canonical TypeScript mirror of LIB/statmech.h
//
// This module is the ONLY place in the TypeScript workspace that decides what
// a numeric result is allowed to be called. C++ (LIB/statmech.h and
// LIB/statmech.cpp) is the source of truth; this file mirrors it exactly and
// must not weaken it. Every other package/app re-exports from here so the
// predicates cannot drift apart.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

/** Scientific schema version for thermodynamic provenance. */
export const SCIENTIFIC_PROVENANCE_SCHEMA_VERSION = 2 as const;

/** Domain of the values supplied to a statistical-mechanics calculation. */
export type EnergyDomain =
  | 'unclassified'
  | 'cf_arbitrary_units'
  | 'calibrated_kcal_per_mol'
  | 'model_scale';

/** Measure represented by the records in an ensemble. */
export type EnsembleMeasure =
  | 'unclassified'
  | 'optimizer_samples'
  | 'enumerated_microstates'
  | 'weighted_quadrature';

/** Reference-state completeness for association claims. */
export type ReferenceState =
  | 'none'
  | 'bound_only'
  | 'matched_association_cycle';

/** Strongest scientific claim supported by the declared provenance. */
export type ClaimValidity =
  | 'proxy_only'
  | 'canonical_physical'
  | 'binding_physical';

/**
 * Evidence required before ensemble diagnostics may be interpreted physically.
 *
 * `claimValidity` is deliberately absent. Consumers derive it from the
 * evidence fields rather than trusting a serialized assertion.
 */
export interface ScientificProvenance {
  schemaVersion: typeof SCIENTIFIC_PROVENANCE_SCHEMA_VERSION;
  energyDomain: EnergyDomain;
  ensembleMeasure: EnsembleMeasure;
  referenceState: ReferenceState;
  energyProvenance: string;
  measureProvenance: string;
  referenceProvenance: string;
}

/**
 * Minimal evidence-bearing source used by derived analyses.
 *
 * `available` is a transport / record-integrity signal, NOT a substitute for
 * provenance. Both must hold before anything is presented physically.
 */
export interface ThermodynamicClaimSource {
  available?: boolean;
  scientificProvenance?: ScientificProvenance;
}

/** Explicit fail-closed provenance for legacy or unclassified results. */
export const PROXY_ONLY_PROVENANCE: ScientificProvenance = Object.freeze({
  schemaVersion: SCIENTIFIC_PROVENANCE_SCHEMA_VERSION,
  energyDomain: 'unclassified',
  ensembleMeasure: 'unclassified',
  referenceState: 'none',
  energyProvenance: '',
  measureProvenance: '',
  referenceProvenance: '',
});

type UnknownRecord = Record<string, unknown>;

/** Narrow an untrusted wire value to a plain record (arrays included). */
export function asRecord(value: unknown): UnknownRecord | undefined {
  return typeof value === 'object' && value !== null
    ? value as UnknownRecord
    : undefined;
}

const KNOWN_EMPTY_OR_FILLER_SHA256 = new Set([
  // SHA-256 of the empty byte string — never identifies a real artifact.
  'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855',
  // Historical entropy.help placeholder digest.
  '3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a',
]);

/**
 * Require a non-trivial SHA-256 artifact identity, not a prose assertion.
 *
 * Mirrors has_artifact_sha256() in LIB/statmech.cpp: exact 71-char length,
 * lowercase literal `sha256:` prefix, 64 hex digits in either case, at least
 * three distinct nibble values, and no known filler digest. No trimming.
 */
export function hasArtifactSha256(value: unknown): value is string {
  if (typeof value !== 'string') return false;
  if (value.length !== 71) return false;
  const match = /^sha256:([0-9a-fA-F]{64})$/.exec(value);
  if (!match) return false;

  const digest = match[1].toLowerCase();
  return new Set(digest).size >= 3
    && !KNOWN_EMPTY_OR_FILLER_SHA256.has(digest);
}

/**
 * Strict transport gate: only a literal boolean `true` counts as available.
 *
 * `undefined`, `null`, `false`, `0`, `1`, `"true"`, `"false"`, `[]` and `{}`
 * must all fail closed — a record that never declared availability has not
 * earned a physical presentation.
 */
export function hasStrictAvailability(value: unknown): boolean {
  return value === true;
}

/**
 * Normalize untrusted camelCase or snake_case scientific metadata.
 *
 * Unknown schemas and values are downgraded to the explicit proxy default.
 * Serialized `claimValidity` / `claim_validity` fields are ignored.
 */
export function normalizeScientificProvenance(value: unknown): ScientificProvenance {
  const raw = asRecord(value);
  if (!raw) return { ...PROXY_ONLY_PROVENANCE };

  const schemaVersion = raw.schemaVersion ?? raw.schema_version;
  if (schemaVersion !== SCIENTIFIC_PROVENANCE_SCHEMA_VERSION) {
    return { ...PROXY_ONLY_PROVENANCE };
  }

  const energyDomain = raw.energyDomain ?? raw.energy_domain;
  const ensembleMeasure = raw.ensembleMeasure ?? raw.ensemble_measure;
  const referenceState = raw.referenceState ?? raw.reference_state;
  const energyProvenance = raw.energyProvenance ?? raw.energy_provenance;
  const measureProvenance = raw.measureProvenance ?? raw.measure_provenance;
  const referenceProvenance = raw.referenceProvenance ?? raw.reference_provenance;

  return {
    schemaVersion: SCIENTIFIC_PROVENANCE_SCHEMA_VERSION,
    energyDomain: (
      energyDomain === 'cf_arbitrary_units'
      || energyDomain === 'calibrated_kcal_per_mol'
      || energyDomain === 'model_scale'
    ) ? energyDomain : 'unclassified',
    ensembleMeasure: (
      ensembleMeasure === 'optimizer_samples'
      || ensembleMeasure === 'enumerated_microstates'
      || ensembleMeasure === 'weighted_quadrature'
    ) ? ensembleMeasure : 'unclassified',
    referenceState: (
      referenceState === 'bound_only'
      || referenceState === 'matched_association_cycle'
    ) ? referenceState : 'none',
    energyProvenance: typeof energyProvenance === 'string' ? energyProvenance : '',
    measureProvenance: typeof measureProvenance === 'string' ? measureProvenance : '',
    referenceProvenance: typeof referenceProvenance === 'string' ? referenceProvenance : '',
  };
}

/** Derive claim validity; serialized validity labels are never authoritative. */
export function deriveClaimValidity(provenance?: unknown): ClaimValidity {
  const normalized = normalizeScientificProvenance(provenance);
  const canonical =
    normalized.schemaVersion === SCIENTIFIC_PROVENANCE_SCHEMA_VERSION
    && normalized.energyDomain === 'calibrated_kcal_per_mol'
    && (normalized.ensembleMeasure === 'enumerated_microstates'
      || normalized.ensembleMeasure === 'weighted_quadrature')
    && hasArtifactSha256(normalized.energyProvenance)
    && hasArtifactSha256(normalized.measureProvenance);

  if (!canonical) return 'proxy_only';

  const binding =
    normalized.referenceState === 'matched_association_cycle'
    && hasArtifactSha256(normalized.referenceProvenance);
  return binding ? 'binding_physical' : 'canonical_physical';
}

/**
 * Strongest claim a whole record supports: availability AND provenance.
 *
 * This is the single predicate every viewer/SDK consumer must use; there are
 * deliberately no per-call-site variants of the availability check.
 */
export function claimValidityForRecord(
  source?: ThermodynamicClaimSource | null,
): ClaimValidity {
  if (!source) return 'proxy_only';
  const raw = source as UnknownRecord;
  if (!hasStrictAvailability(raw.available)) return 'proxy_only';
  return deriveClaimValidity(raw.scientificProvenance);
}

/** Whether a thermodynamic record supports canonical physical claims. */
export function allowsCanonicalClaims(
  source?: ThermodynamicClaimSource | null,
): boolean {
  return claimValidityForRecord(source) !== 'proxy_only';
}

/** Whether a thermodynamic record supports physical association claims. */
export function allowsBindingClaims(
  source?: ThermodynamicClaimSource | null,
): boolean {
  return claimValidityForRecord(source) === 'binding_physical';
}
