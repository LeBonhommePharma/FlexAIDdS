// IntelligenceEngine.ts — Client-side intelligence for BonhommeViewer
//
// On Apple devices with FoundationModels, delegates to the native bridge.
// On web, provides rule-based analysis with structured confidence metadata,
// trend tracking, and enthalpy-entropy compensation detection.
// Matches Swift IntelligenceOracle and RuleBasedOracle APIs.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import {
  allowsBindingClaims,
  claimValidityForRecord,
} from '@bonhomme/shared';
import type {
  BindingPopulation,
  ClaimValidity,
  HealthCorrelation,
  ShannonEntropyDecomposition,
  ThermodynamicClaimSource,
} from '@bonhomme/shared';
import type { RefereeFinding, RefereeVerdict, RefereeSeverity, RefereeCategory } from '@bonhomme/shared';

/** Confidence level for each analysis bullet. */
export type AnalysisConfidence = 'high' | 'moderate' | 'low';

/** Category of analysis bullet. */
export type BulletCategory =
  | 'binding'
  | 'entropy'
  | 'health'
  | 'modification'
  | 'fleet'
  | 'trend';

/** A single analysis bullet with metadata. */
export interface AnalysisBullet {
  text: string;
  confidence: AnalysisConfidence;
  category: BulletCategory;
}

/** Structured oracle analysis result. */
export interface OracleAnalysis {
  bullets: string[];
  structuredBullets: AnalysisBullet[];
  overallConfidence: AnalysisConfidence;
  inputSummary: string;
  timestamp: string;
  /** Derived from evidence, never from a serialized claim label. */
  claimValidity: ClaimValidity;
  /** Evidence-bearing source used to re-derive claim validity. */
  claimSource: ThermodynamicClaimSource;
}

/** Stored analysis entry for trend tracking. */
interface AnalysisEntry {
  key: string;
  analysis: OracleAnalysis;
}

export class IntelligenceEngine {
  /** Analysis history for trend detection. */
  private static history: AnalysisEntry[] = [];
  private static readonly MAX_HISTORY = 50;

  /**
   * Referee a binding population and return structured analysis.
   *
   * Uses rule-based referee analysis (same thresholds as the Swift RuleBasedOracle)
   * with Shannon entropy decomposition, convergence diagnostics, and
   * vibrational/configurational entropy split when available.
   */
  static async analyze(
    population: BindingPopulation,
    health?: HealthCorrelation,
    campaignKey?: string,
  ): Promise<OracleAnalysis> {
    const structured: AnalysisBullet[] = [];
    const thermo = population.globalThermodynamics;
    // Availability and provenance are both required; the strict gate lives in
    // @bonhomme/shared so this call site cannot soften it.
    const claimValidity = claimValidityForRecord(thermo);
    const bindingPhysical = allowsBindingClaims(thermo);
    const decomp: ShannonEntropyDecomposition | undefined = health?.shannonDecomposition;
    const isConverged = decomp?.isConverged ?? false;

    // Bullet 1: Free energy assessment — gate confidence on convergence
    const energyStable =
      thermo.stdEnergy !== undefined &&
      thermo.stdEnergy < Math.abs(thermo.freeEnergy) * 0.5;
    const fConfidence: AnalysisConfidence =
      energyStable && isConverged ? 'high' : energyStable ? 'moderate' : 'low';

    if (!bindingPhysical) {
      structured.push({
        text: `Ensemble/CF diagnostic = ${thermo.freeEnergy.toFixed(1)} (source units); physical affinity unavailable because binding_physical provenance is absent.`,
        confidence: fConfidence,
        category: 'binding',
      });
    } else if (thermo.freeEnergy < -10) {
      structured.push({
        text: `Strong binding affinity (F = ${thermo.freeEnergy.toFixed(1)} kcal/mol)${isConverged ? '' : ' — but entropy has NOT converged, F may shift with more sampling'}.`,
        confidence: fConfidence,
        category: 'binding',
      });
    } else if (thermo.freeEnergy < -5) {
      structured.push({
        text: `Moderate binding affinity (F = ${thermo.freeEnergy.toFixed(1)} kcal/mol)${isConverged ? ' — converged ensemble' : ' — entropy not converged, consider more GA generations'}.`,
        confidence: fConfidence,
        category: 'binding',
      });
    } else {
      structured.push({
        text: `Weak binding affinity (F = ${thermo.freeEnergy.toFixed(1)} kcal/mol) — consider structural optimization${isConverged ? '' : ' after convergence is reached'}.`,
        confidence: fConfidence,
        category: 'binding',
      });
    }

    // Bullet 2: Shannon entropy decomposition referee
    if (decomp) {
      // Convergence check — primary referee concern
      if (!decomp.isConverged) {
        structured.push({
          text: bindingPhysical
            ? `Entropy NOT converged (rate = ${decomp.convergenceRate.toFixed(4)}). Increase GA generations or population size before trusting thermodynamic values.`
            : `Entropy NOT converged (rate = ${decomp.convergenceRate.toFixed(4)}). Increase GA generations or population size before interpreting ensemble/CF diagnostics; physical affinity unavailable.`,
          confidence: 'high',
          category: 'entropy',
        });
      }

      // Vibrational dominance check
      const kB = 0.001987206;
      const sConfPhysical = decomp.configurational * kB;
      if (decomp.vibrational > sConfPhysical * 3.0 && decomp.vibrational > 0.001) {
        structured.push({
          text: bindingPhysical
            ? `Vibrational entropy dominates (S_vib = ${decomp.vibrational.toFixed(6)} >> S_conf = ${sConfPhysical.toFixed(6)} kcal/mol/K). Protein backbone flexibility drives entropy — ligand conformational space may be under-explored.`
            : `Ensemble/CF entropy diagnostic: S_vib = ${decomp.vibrational.toFixed(6)} >> S_conf = ${sConfPhysical.toFixed(6)} (source units). Physical affinity unavailable; physical energy units are not authorized.`,
          confidence: 'high',
          category: 'entropy',
        });
      }

      // Sparse histogram check
      if (decomp.totalBins > 0) {
        const occupancyRatio = decomp.occupiedBins / decomp.totalBins;
        if (occupancyRatio < 0.5) {
          structured.push({
            text: `Sparse energy histogram (${decomp.occupiedBins}/${decomp.totalBins} bins, ${(occupancyRatio * 100).toFixed(0)}% occupied). Energy landscape poorly sampled — increase ensemble size.`,
            confidence: 'high',
            category: 'entropy',
          });
        }
      }

      // Per-mode entropy imbalance
      if (decomp.perModeEntropy.length >= 2) {
        const maxS = Math.max(...decomp.perModeEntropy);
        const positiveEntropies = decomp.perModeEntropy.filter((s) => s > 0);
        const minS = positiveEntropies.length > 0 ? Math.min(...positiveEntropies) : 0;
        if (minS > 0 && maxS > minS * 10) {
          const ratio = maxS / minS;
          structured.push({
            text: `Per-mode entropy imbalance (${ratio.toFixed(1)}x ratio). One binding mode absorbs most conformational diversity — check for kinetic trapping.`,
            confidence: 'moderate',
            category: 'entropy',
          });
        }
      }

      // Converged and well-sampled: report summary
      if (decomp.isConverged && structured.filter((b) => b.category === 'entropy').length === 0) {
        const modeCount = population.modes.length;
        structured.push({
          text: bindingPhysical
            ? `Shannon entropy converged: S_conf = ${decomp.configurational.toFixed(4)} nats, S_vib = ${decomp.vibrational.toFixed(6)} kcal/mol/K across ${modeCount} modes (${decomp.hardwareBackend} backend).`
            : `Ensemble/CF entropy diagnostic converged: S_conf = ${decomp.configurational.toFixed(4)} nats, S_vib = ${decomp.vibrational.toFixed(6)} (source units) across ${modeCount} modes (${decomp.hardwareBackend} backend); physical affinity unavailable.`,
          confidence: 'high',
          category: 'entropy',
        });
      }
    } else {
      // Fallback: scalar-only entropy (no decomposition)
      const modeCount = population.modes.length;
      const sConfidence: AnalysisConfidence = modeCount >= 3 ? 'moderate' : 'low';

      if (population.isCollapsed) {
        structured.push({
          text: bindingPhysical
            ? `Entropy collapsed to ${modeCount} mode(s) — high specificity but check for enthalpy-entropy compensation.`
            : `Ensemble/CF diagnostic collapsed to ${modeCount} mode(s). Physical affinity unavailable; do not infer specificity or enthalpy-entropy compensation.`,
          confidence: sConfidence,
          category: 'entropy',
        });
      } else if (population.shannonS > 0.5) {
        structured.push({
          text: bindingPhysical
            ? `High conformational entropy (S = ${population.shannonS.toFixed(4)}) — population still exploring. More sampling may refine.`
            : `High ensemble/CF configurational diagnostic (S slot = ${population.shannonS.toFixed(4)} source units) — population still exploring; physical affinity unavailable.`,
          confidence: 'moderate',
          category: 'entropy',
        });
      } else {
        structured.push({
          text: bindingPhysical
            ? `Moderate entropy with ${modeCount} binding modes — population converging. Enable ShannonThermoStack for decomposed referee analysis.`
            : `Moderate ensemble/CF entropy diagnostic across ${modeCount} modes. Physical affinity unavailable; enable ShannonThermoStack for decomposed diagnostic analysis.`,
          confidence: sConfidence,
          category: 'entropy',
        });
      }
    }

    // Bullet 3: Enthalpy-entropy compensation detection
    if (thermo.freeEnergy < -5 && thermo.entropy > 0.01) {
      structured.push({
        text: bindingPhysical
          ? `Enthalpy-entropy compensation: strong binding (F = ${thermo.freeEnergy.toFixed(1)}) offset by conformational flexibility (S = ${thermo.entropy.toFixed(4)}). Net \u0394G may be less favorable than F alone suggests.`
          : `Ensemble/CF diagnostic co-variation: score = ${thermo.freeEnergy.toFixed(1)}, entropy slot = ${thermo.entropy.toFixed(4)} (source units). Physical affinity unavailable; this does not establish enthalpy-entropy compensation.`,
        confidence: 'moderate',
        category: 'binding',
      });
    }

    // Bullet 4: Target modifications with population impact
    if (population.targetModifications?.length) {
      const mods = population.targetModifications;
      const modSummary = mods.map((m) => `${m.type}@${m.residueName}${m.residueNumber}`).join(', ');
      const deltaF = mods
        .filter((m) => m.effect?.deltaFreeEnergy !== undefined)
        .reduce((sum, m) => sum + (m.effect?.deltaFreeEnergy ?? 0), 0);
      const deltaS = mods
        .filter((m) => m.effect?.deltaEntropy !== undefined)
        .reduce((sum, m) => sum + (m.effect?.deltaEntropy ?? 0), 0);

      let modText = `${mods.length} target modification(s) (${modSummary})`;
      if (deltaF !== 0 || deltaS !== 0) {
        modText += bindingPhysical
          ? ` — net population shift: \u0394F=${deltaF.toFixed(2)} kcal/mol, \u0394S=${deltaS.toFixed(4)} kcal/mol/K`
          : ` — net ensemble/CF diagnostic shift: score=${deltaF.toFixed(2)}, entropy=${deltaS.toFixed(4)} (source units); physical affinity unavailable`;
      }
      modText += bindingPhysical
        ? '. Population recalculated with PTM/glycan effects on the binding landscape.'
        : '. Population diagnostics recalculated with PTM/glycan inputs.';

      structured.push({
        text: modText,
        confidence: 'moderate',
        category: 'modification',
      });
    }

    // Bullet 5: Health correlation
    if (health?.hrvSDNN != null) {
      if (population.isCollapsed && health.hrvSDNN > 60) {
        structured.push({
          text: `Entropy collapse correlates with good HRV (${health.hrvSDNN.toFixed(0)} ms) — system recovering. Gentle activity recommended.`,
          confidence: 'moderate',
          category: 'health',
        });
      } else if (health.hrvSDNN < 40) {
        structured.push({
          text: `Low HRV (${health.hrvSDNN.toFixed(0)} ms) — prioritize rest before interpreting docking results.`,
          confidence: 'high',
          category: 'health',
        });
      } else {
        structured.push({
          text: `HRV at ${health.hrvSDNN.toFixed(0)} ms — stable physiological state for analysis.`,
          confidence: 'moderate',
          category: 'health',
        });
      }
    } else if (!population.targetModifications?.length) {
      structured.push({
        text: 'Connect HealthKit for entropy-health correlation. Enable fleet mode for distributed compute.',
        confidence: 'low',
        category: 'fleet',
      });
    }

    const overallConfidence = IntelligenceEngine.computeOverallConfidence(structured);
    let inputSummary = bindingPhysical
      ? `T=${thermo.temperature}K, F=${thermo.freeEnergy.toFixed(2)} kcal/mol`
      : `T=${thermo.temperature}K, diagnostic=${thermo.freeEnergy.toFixed(2)} source-units, physical-affinity=unavailable`;
    if (decomp) {
      inputSummary += bindingPhysical
        ? `, S_conf=${decomp.configurational.toFixed(4)}nats, S_vib=${decomp.vibrational.toFixed(6)}kcal/mol/K`
        : `, S_conf=${decomp.configurational.toFixed(4)}nats, S_vib-diagnostic=${decomp.vibrational.toFixed(6)} source-units`;
      inputSummary += decomp.isConverged ? ' [converged]' : ' [not converged]';
    }

    const analysis: OracleAnalysis = {
      bullets: structured.map((b) => b.text),
      structuredBullets: structured,
      overallConfidence,
      inputSummary,
      timestamp: new Date().toISOString(),
      claimValidity,
      claimSource: {
        available: thermo.available,
        scientificProvenance: thermo.scientificProvenance,
      },
    };

    // Record for trend tracking
    if (campaignKey) {
      IntelligenceEngine.recordAnalysis(campaignKey, analysis);
    }

    return analysis;
  }

  /**
   * Compare current results with a previous analysis for the same campaign.
   */
  static compareTrend(
    campaignKey: string,
    current: OracleAnalysis,
  ): AnalysisBullet | null {
    const previous = IntelligenceEngine.lastAnalysis(campaignKey);
    if (!previous) return null;

    // Parse the unchanged numeric slot from physical or diagnostic summaries.
    const extractScore = (summary: string): number | null => {
      const match = summary.match(/(?:F|diagnostic)=(-?[\d.]+)/);
      return match ? parseFloat(match[1]) : null;
    };

    const prevF = extractScore(previous.inputSummary);
    const currF = extractScore(current.inputSummary);

    if (prevF !== null && currF !== null) {
      const delta = currF - prevF;
      const improved = delta < 0;
      if (!allowsBindingClaims(previous.claimSource)
        || !allowsBindingClaims(current.claimSource)) {
        return {
          text: `Compared to previous run: ensemble/CF diagnostic shift = ${delta.toFixed(2)} source units (${improved ? 'decreased' : 'increased'}); physical affinity unavailable.`,
          confidence: Math.abs(delta) > 1.0 ? 'high' : 'moderate',
          category: 'trend',
        };
      }
      return {
        text: `Compared to previous run: \u0394F = ${delta.toFixed(2)} kcal/mol (${improved ? 'improved' : 'worsened'} binding).`,
        confidence: Math.abs(delta) > 1.0 ? 'high' : 'moderate',
        category: 'trend',
      };
    }

    return null;
  }

  /**
   * Get analysis history for a campaign key.
   */
  static getHistory(campaignKey: string): OracleAnalysis[] {
    return IntelligenceEngine.history
      .filter((e) => e.key === campaignKey)
      .map((e) => e.analysis);
  }

  // ─── Private helpers ───────────────────────────────────────────────────────

  private static computeOverallConfidence(bullets: AnalysisBullet[]): AnalysisConfidence {
    const highCount = bullets.filter((b) => b.confidence === 'high').length;
    const lowCount = bullets.filter((b) => b.confidence === 'low').length;
    if (highCount >= 2) return 'high';
    if (lowCount >= 2) return 'low';
    return 'moderate';
  }

  private static recordAnalysis(key: string, analysis: OracleAnalysis): void {
    IntelligenceEngine.history.push({ key, analysis });
    if (IntelligenceEngine.history.length > IntelligenceEngine.MAX_HISTORY) {
      IntelligenceEngine.history.shift();
    }
  }

  private static lastAnalysis(key: string): OracleAnalysis | null {
    for (let i = IntelligenceEngine.history.length - 1; i >= 0; i--) {
      if (IntelligenceEngine.history[i].key === key) {
        return IntelligenceEngine.history[i].analysis;
      }
    }
    return null;
  }
}

// ─── Rule-Based Referee (Cross-Platform) ──────────────────────────────────────

/**
 * Deterministic thermodynamic referee for web/non-Apple platforms.
 * Mirrors the Swift RuleBasedReferee with identical threshold logic.
 * Produces typed RefereeVerdict with findings, severity, and confidence.
 */
export class RuleBasedReferee {
  /**
   * Run deterministic referee analysis on a binding population.
   */
  static referee(
    population: BindingPopulation,
    health?: HealthCorrelation,
  ): RefereeVerdict {
    const findings: RefereeFinding[] = [];
    const thermo = population.globalThermodynamics;
    const bindingPhysical = allowsBindingClaims(thermo);
    let trustworthy = bindingPhysical;
    const decomp: ShannonEntropyDecomposition | undefined = health?.shannonDecomposition;
    const kB = 0.001987206;

    // 1. Convergence check (highest priority)
    if (decomp) {
      if (!decomp.isConverged) {
        trustworthy = false;
        findings.push({
          title: 'Entropy not converged',
          detail: bindingPhysical
            ? `Shannon entropy has not reached a plateau (rate = ${decomp.convergenceRate.toFixed(4)}). Increase GA generations before trusting F or S values.`
            : `Shannon entropy has not reached a plateau (rate = ${decomp.convergenceRate.toFixed(4)}). Increase GA generations before interpreting ensemble/CF diagnostics; physical affinity unavailable.`,
          severity: 'critical',
          category: 'convergence',
        });
      } else {
        findings.push({
          title: 'Entropy converged',
          detail: bindingPhysical
            ? `Shannon entropy plateau reached on ${decomp.hardwareBackend} backend. Thermodynamic values are reliable.`
            : `Shannon entropy plateau reached on ${decomp.hardwareBackend} backend for the ensemble/CF diagnostic. Physical affinity unavailable without binding_physical provenance.`,
          severity: 'pass',
          category: 'convergence',
        });
      }
    } else {
      findings.push({
        title: 'No Shannon decomposition',
        detail: 'ShannonThermoStack data unavailable. Cannot assess convergence — enable ShannonThermoStack for reliable referee analysis.',
        severity: 'warning',
        category: 'convergence',
      });
    }

    // 2. Histogram quality
    if (decomp && decomp.totalBins > 0) {
      const occupancy = decomp.occupiedBins / decomp.totalBins;
      if (occupancy < 0.3) {
        trustworthy = false;
        findings.push({
          title: 'Critically sparse histogram',
          detail: `Only ${decomp.occupiedBins}/${decomp.totalBins} bins occupied (${(occupancy * 100).toFixed(0)}%). Energy landscape severely under-sampled. Double the population size.`,
          severity: 'critical',
          category: 'histogram',
        });
      } else if (occupancy < 0.5) {
        findings.push({
          title: 'Sparse histogram',
          detail: `${decomp.occupiedBins}/${decomp.totalBins} bins occupied (${(occupancy * 100).toFixed(0)}%). Consider increasing ensemble size for better coverage.`,
          severity: 'warning',
          category: 'histogram',
        });
      }
    }

    // 3. Entropy balance (S_vib vs S_conf)
    if (decomp) {
      const sConfPhysical = decomp.configurational * kB;
      if (decomp.vibrational > sConfPhysical * 3.0 && decomp.vibrational > 0.001) {
        findings.push({
          title: 'Vibrational entropy dominates',
          detail: bindingPhysical
            ? `S_vib (${decomp.vibrational.toFixed(6)}) >> S_conf (${sConfPhysical.toFixed(6)}) kcal/mol/K. Protein backbone flexibility drives entropy; ligand conformational space may be under-explored.`
            : `Ensemble/CF entropy diagnostic: S_vib (${decomp.vibrational.toFixed(6)}) >> S_conf (${sConfPhysical.toFixed(6)}) in source units. Physical affinity unavailable; physical energy units are not authorized.`,
          severity: 'warning',
          category: 'entropyBalance',
        });
      }
    }

    // 4. Mode entropy imbalance
    if (decomp && decomp.perModeEntropy.length >= 2) {
      const maxS = Math.max(...decomp.perModeEntropy);
      const positiveEntropies = decomp.perModeEntropy.filter((s) => s > 0);
      const minS = positiveEntropies.length > 0 ? Math.min(...positiveEntropies) : 0;
      if (minS > 0 && maxS > minS * 10) {
        const ratio = maxS / minS;
        findings.push({
          title: 'Mode entropy imbalance',
          detail: `Entropy ratio across modes is ${ratio.toFixed(1)}x. One binding mode absorbs most conformational diversity — check for kinetic trapping.`,
          severity: 'warning',
          category: 'modeBalance',
        });
      }
    }

    // 5. Enthalpy-entropy compensation
    if (thermo.freeEnergy < -5 && thermo.entropy > 0.01) {
      findings.push({
        title: bindingPhysical ? 'Enthalpy-entropy compensation' : 'Diagnostic co-variation only',
        detail: bindingPhysical
          ? `Strong binding (F = ${thermo.freeEnergy.toFixed(1)} kcal/mol) offset by conformational flexibility (S = ${thermo.entropy.toFixed(4)}). Net \u0394G may be less favorable than F alone suggests.`
          : `Ensemble/CF diagnostic score ${thermo.freeEnergy.toFixed(1)} co-varies with entropy slot ${thermo.entropy.toFixed(4)} (source units). Physical affinity unavailable; compensation is not established.`,
        severity: 'advisory',
        category: 'compensation',
      });
    }

    // 6. Binding affinity
    const converged = decomp?.isConverged ?? false;
    if (!bindingPhysical) {
      findings.push({
        title: 'Physical affinity unavailable',
        detail: `Ensemble/CF diagnostic = ${thermo.freeEnergy.toFixed(1)} (source units). Missing binding_physical provenance; do not interpret this value as a physical association free energy.`,
        severity: 'warning',
        category: 'affinity',
      });
    } else if (thermo.freeEnergy < -10) {
      findings.push({
        title: 'Strong binding affinity',
        detail: `F = ${thermo.freeEnergy.toFixed(1)} kcal/mol${converged ? ' (converged)' : ' — convergence not confirmed'}.`,
        severity: converged ? 'pass' : 'advisory',
        category: 'affinity',
      });
    } else if (thermo.freeEnergy < -5) {
      findings.push({
        title: 'Moderate binding affinity',
        detail: `F = ${thermo.freeEnergy.toFixed(1)} kcal/mol — reasonable drug candidate${converged ? '' : ', pending convergence'}.`,
        severity: 'advisory',
        category: 'affinity',
      });
    } else {
      findings.push({
        title: 'Weak binding affinity',
        detail: `F = ${thermo.freeEnergy.toFixed(1)} kcal/mol — consider structural optimization.`,
        severity: 'warning',
        category: 'affinity',
      });
    }

    // Recommendation
    let recommendedAction: string;
    if (!bindingPhysical) {
      recommendedAction = 'Physical affinity unavailable. Treat values as ensemble/CF diagnostics and supply calibrated energy, ensemble-measure, and matched association-cycle SHA-256 artifact identities.';
    } else if (!trustworthy) {
      recommendedAction = 'Do not trust current thermodynamic values. Increase GA generations and population size, then re-run.';
    } else if (!decomp) {
      recommendedAction = 'Enable ShannonThermoStack for decomposed entropy analysis.';
    } else if (thermo.freeEnergy > -5) {
      recommendedAction = 'Consider structural optimization of the ligand or alternative binding sites.';
    } else {
      recommendedAction = 'Results are reliable. Proceed with lead optimization.';
    }

    let confidence: number;
    if (!trustworthy) confidence = 0.3;
    else if (decomp?.isConverged) confidence = 0.9;
    else confidence = 0.6;

    return {
      findings,
      overallTrustworthy: trustworthy,
      recommendedAction,
      confidence,
      claimSource: {
        available: thermo.available,
        scientificProvenance: thermo.scientificProvenance,
      },
    };
  }
}
