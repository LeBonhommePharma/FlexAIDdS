// SelectivityAnalyzer.ts — Rule-based multi-target selectivity analysis
//
// Ports the Swift RuleBasedSelectivityAnalyst to TypeScript for web parity.
// Compares docking results across protein targets to determine
// selectivity drivers (enthalpic vs entropic).
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import { claimValidityForRecord } from '@bonhomme/shared';
import type {
  ClaimValidity,
  SelectivityContext,
  TargetDockingSummary,
  ThermodynamicClaimSource,
} from '@bonhomme/shared';
import type { SelectivityDriver, SelectivityAnalysis } from '@bonhomme/shared';

/** A target row carries its evidence in two fields; keep them together. */
function claimSourceFor(target: TargetDockingSummary): ThermodynamicClaimSource {
  return {
    available: target.thermodynamicsAvailable,
    scientificProvenance: target.scientificProvenance,
  };
}

// Named thresholds (matching Swift RuleBasedSelectivityAnalyst)
const INCONCLUSIVE_ENTHALPIC_THRESHOLD = 0.5; // kcal/mol
const INCONCLUSIVE_ENTROPIC_THRESHOLD = 0.0001; // kcal/mol/K
const ENTROPIC_DOMINANCE_RATIO = 1.3;
const ENTHALPIC_SIGNIFICANCE_THRESHOLD = 2.0; // kcal/mol
const kB = 0.001987206;

export class SelectivityAnalyzer {
  /**
   * Analyze selectivity across targets using threshold logic.
   */
  static analyze(context: SelectivityContext): SelectivityAnalysis {
    if (context.targets.length < 2) {
      return {
        preferredTarget: context.targets[0]?.targetName ?? 'unknown',
        deltaG: 0,
        driver: 'inconclusive',
        explanation: 'Only one target ensemble/CF diagnostic is available. Physical affinity unavailable; selectivity analysis requires at least two binding_physical targets.',
        designSuggestion: 'Dock additional targets and supply matched binding_physical provenance before affinity or selectivity claims.',
        claimValidity: 'proxy_only',
        claimSources: context.targets.map(claimSourceFor),
      };
    }

    const sorted = [...context.targets].sort((a, b) => a.bestFreeEnergy - b.bestFreeEnergy);
    const preferred = sorted[0];
    const second = sorted[1];
    const ddg = preferred.bestFreeEnergy - second.bestFreeEnergy;
    // Strict gate from @bonhomme/shared: a target with no declared
    // availability is proxy-only, exactly like one that declared false.
    const preferredValidity = claimValidityForRecord(claimSourceFor(preferred));
    const secondValidity = claimValidityForRecord(claimSourceFor(second));
    const claimValidity: ClaimValidity =
      preferredValidity === 'binding_physical' && secondValidity === 'binding_physical'
        ? 'binding_physical'
        : preferredValidity !== 'proxy_only' && secondValidity !== 'proxy_only'
          ? 'canonical_physical'
          : 'proxy_only';
    const bindingPhysical = claimValidity === 'binding_physical';
    const claimSources = [preferred, second].map(claimSourceFor);

    // Determine driver
    const sConfPhysA = preferred.sConf * kB;
    const sConfPhysB = second.sConf * kB;
    const entropicDiff = Math.abs(sConfPhysA - sConfPhysB);
    const enthalpicDiff = Math.abs(ddg);

    let driver: SelectivityDriver;
    if (enthalpicDiff < INCONCLUSIVE_ENTHALPIC_THRESHOLD && entropicDiff < INCONCLUSIVE_ENTROPIC_THRESHOLD) {
      driver = 'inconclusive';
    } else if (entropicDiff > enthalpicDiff * 0.001 && preferred.sConf > second.sConf * ENTROPIC_DOMINANCE_RATIO) {
      driver = 'entropic';
    } else if (enthalpicDiff > ENTHALPIC_SIGNIFICANCE_THRESHOLD) {
      driver = 'enthalpic';
    } else {
      driver = 'mixed';
    }

    // "enthalpic" / "entropic" / "mixed" are thermodynamic-driver claims about
    // an association process. Proxy CF diagnostics cannot support one, so the
    // published driver is inconclusive regardless of the internal threshold
    // outcome. The threshold arithmetic above is unchanged.
    const publishedDriver: SelectivityDriver = bindingPhysical ? driver : 'inconclusive';

    // Explanation
    let explanation: string;
    if (!bindingPhysical) {
      explanation = `Ensemble/CF diagnostic ranks ${preferred.targetName} (${preferred.bestFreeEnergy.toFixed(1)}) ahead of ${second.targetName} (${second.bestFreeEnergy.toFixed(1)}), with diagnostic difference ${ddg.toFixed(2)} source units. Physical affinity unavailable; no binding selectivity or thermodynamic driver is established.`;
    } else {
      explanation = `${context.ligandName} prefers ${preferred.targetName} (F = ${preferred.bestFreeEnergy.toFixed(1)} kcal/mol) over ${second.targetName} (F = ${second.bestFreeEnergy.toFixed(1)} kcal/mol). `;
      explanation += `DDG = ${ddg.toFixed(2)} kcal/mol. `;

      switch (driver) {
        case 'entropic':
          explanation += `Selectivity is entropy-driven: ${preferred.targetName} has broader conformational ensemble (S_conf = ${preferred.sConf.toFixed(4)} vs ${second.sConf.toFixed(4)} nats).`;
          break;
        case 'enthalpic':
          explanation += `Selectivity is enthalpy-driven: stronger direct interactions at ${preferred.targetName}.`;
          break;
        case 'mixed':
          explanation += 'Both enthalpy and entropy contribute to selectivity.';
          break;
        case 'inconclusive':
          explanation += 'Selectivity is marginal — results may not be significant.';
          break;
      }
    }

    // Design suggestion
    let designSuggestion: string;
    if (!bindingPhysical) {
      designSuggestion = 'Treat the ordering as an ensemble/CF diagnostic only. Supply matched binding_physical provenance before affinity-driven design.';
    } else {
      switch (driver) {
        case 'entropic':
          designSuggestion = `Rigidify ligand to reduce entropy-driven selectivity for ${preferred.targetName}, or add flexible groups to enhance selectivity.`;
          break;
        case 'enthalpic':
          designSuggestion = `Optimize interaction geometry at ${second.targetName} binding site to improve affinity and shift selectivity.`;
          break;
        case 'mixed':
          designSuggestion = 'Balanced selectivity — consider fragment-based approach targeting unique pocket features of each receptor.';
          break;
        case 'inconclusive':
          designSuggestion = 'Increase sampling (more GA generations) and verify convergence before drawing selectivity conclusions.';
          break;
      }
    }

    // Convergence gating
    const convergedBoth = preferred.isConverged && second.isConverged;
    if (!convergedBoth) {
      const notConverged = [preferred, second].filter((t) => !t.isConverged).map((t) => t.targetName);
      return {
        preferredTarget: preferred.targetName,
        deltaG: ddg,
        driver: 'inconclusive',
        explanation: explanation + ` Note: ${notConverged.join(', ')} not converged — selectivity assessment unreliable.`,
        designSuggestion: bindingPhysical
          ? 'Achieve convergence at all targets before selectivity-driven optimization.'
          : 'Achieve convergence, then supply binding_physical provenance; until then this remains an ensemble/CF diagnostic and physical affinity is unavailable.',
        claimValidity,
        claimSources,
      };
    }

    return {
      preferredTarget: preferred.targetName,
      deltaG: ddg,
      driver: publishedDriver,
      explanation,
      designSuggestion,
      claimValidity,
      claimSources,
    };
  }
}
