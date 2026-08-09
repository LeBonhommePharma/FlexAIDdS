// index.ts — FlexAIDdS TypeScript SDK entry point
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

export type {
  EnergyDomain,
  EnsembleMeasure,
  ReferenceState,
  ClaimValidity,
  ScientificProvenance,
  ThermodynamicClaimSource,
  UnavailableReason,
  NormalizedThermodynamicRecord,
  NormalizedBindingMode,
  NormalizedDockingRecord,
  ThermodynamicResult,
  VibrationalEntropyResult,
  WHAMBinResult,
  TIPoint,
  PoseResult,
  BindingModeResult,
  DockingResult,
  BindingEntropyScore,
  WorkChunk,
  DeviceCapability,
} from './types.js';

export {
  SCIENTIFIC_PROVENANCE_SCHEMA_VERSION,
  PROXY_ONLY_PROVENANCE,
  hasArtifactSha256,
  hasStrictAvailability,
  normalizeScientificProvenance,
  deriveClaimValidity,
  claimValidityForRecord,
  allowsCanonicalClaims,
  allowsBindingClaims,
  normalizeAvailability,
  normalizeThermodynamicRecord,
  normalizeDockingRecord,
  claimValidityForWireRecord,
} from './types.js';

export { StatMechEngine } from './StatMechEngine.js';
export {
  parseResultFile,
  parseDockingResultJSON,
  adoptDockingPayload,
} from './resultLoader.js';
export { kB_kcal, kB_SI, hbar_SI, NA } from './constants.js';
