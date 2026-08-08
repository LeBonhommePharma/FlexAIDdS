// index.ts — Shared types entry point
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

export type {
  BindingPopulation,
  BindingMode,
  Pose,
  Thermodynamics,
  TargetModification,
  HealthCorrelation,
  ShannonEntropyDecomposition,
  EnergyDomain,
  EnsembleMeasure,
  ReferenceState,
  ClaimValidity,
  ScientificProvenance,
  ThermodynamicClaimSource,
} from './BindingPopulation.js';

export type {
  RefereeSeverity,
  RefereeCategory,
  RefereeFinding,
  RefereeVerdict,
  TemperatureSensitivity,
  ComparativeVerdict,
} from './RefereeVerdict.js';

export type {
  ModeDescription,
  BindingModeNarrative,
} from './BindingModeNarrative.js';

export type {
  DruggabilityTier,
  CleftAssessment,
} from './CleftAssessment.js';

export type {
  GAAdvice,
  ConvergenceCoaching,
} from './ConvergenceCoaching.js';

export type {
  FleetExplanation,
} from './FleetExplanation.js';

export type {
  HealthEntropyInsight,
} from './HealthEntropyInsight.js';

export type {
  VibrationalInsight,
} from './VibrationalInsight.js';

export type {
  SelectivityDriver,
  DeltaDeltaG,
  SelectivityAnalysis,
} from './SelectivityAnalysis.js';

export type {
  CampaignSummary,
} from './CampaignSummary.js';

export type {
  PoseQualityReport,
} from './PoseQualityReport.js';

export type {
  CleftFeatures,
} from './CleftFeatures.js';

export type {
  GAProgressSnapshot,
} from './GAProgressSnapshot.js';

export type {
  PoseProfile,
  PoseQualityContext,
} from './PoseQualityContext.js';

export type {
  TargetDockingSummary,
  SelectivityContext,
} from './SelectivityContext.js';

export {
  serializePopulation,
  deserializePopulation,
} from './BindingPopulation.js';

// Single source of truth for the claim firewall. Do not re-implement any of
// these predicates at a call site — import them from here.
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
} from './scientificProvenance.js';

// The one runtime boundary for untrusted producer payloads.
export {
  normalizeAvailability,
  normalizeThermodynamicRecord,
  normalizeDockingRecord,
  claimValidityForWireRecord,
} from './wireNormalization.js';

export type {
  UnavailableReason,
  NormalizedThermodynamicRecord,
  NormalizedBindingMode,
  NormalizedDockingRecord,
} from './wireNormalization.js';
