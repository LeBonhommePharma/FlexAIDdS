// RefereeVerdict.ts — Cross-platform referee verdict types
//
// Mirrors the Swift @Generable ThermoRefereeTypes and
// CrossPlatformRefereeVerdict for TypeScript consumers.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import type { ThermodynamicClaimSource } from './BindingPopulation.js';

/** Severity level for a referee finding. */
export type RefereeSeverity = 'pass' | 'advisory' | 'warning' | 'critical';

/** Category of thermodynamic referee finding. */
export type RefereeCategory =
  | 'convergence'
  | 'entropyBalance'
  | 'compensation'
  | 'histogram'
  | 'modeBalance'
  | 'affinity'
  | 'recommendation';

/** A single referee finding with typed severity and category. */
export interface RefereeFinding {
  /** Short title (3-8 words) */
  title: string;
  /** Detailed explanation with quantitative values */
  detail: string;
  /** Severity of this finding */
  severity: RefereeSeverity;
  /** Category of thermodynamic concern */
  category: RefereeCategory;
}

/** Complete referee verdict. */
export interface RefereeVerdict {
  /** List of findings (typically 3-5) */
  findings: RefereeFinding[];
  /** Can we trust these thermodynamic values? */
  overallTrustworthy: boolean;
  /** Recommended next action */
  recommendedAction: string;
  /** Confidence in the verdict (0.0-1.0) */
  confidence: number;
  /** Evidence-bearing source; consumers must derive validity from it. */
  claimSource?: ThermodynamicClaimSource;
}

/** Temperature sensitivity analysis result. */
export interface TemperatureSensitivity {
  /** Does free energy change significantly with temperature? */
  isSensitive: boolean;
  /** Assessment of the sensitivity */
  assessment: string;
  /** Recommended temperature range */
  recommendedTempRange: string;
}

/** Comparative verdict between two campaigns. */
export interface ComparativeVerdict {
  /** Did binding improve, worsen, or stay similar? */
  bindingTrend: string;
  /** Did entropy converge better or worse? */
  convergenceTrend: string;
  /** Overall comparison summary */
  summary: string;
  /** Is the new run an improvement? */
  improved: boolean;
}
