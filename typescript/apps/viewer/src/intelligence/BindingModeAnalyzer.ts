// BindingModeAnalyzer.ts — Rule-based binding mode narrative generator
//
// Ports the Swift RuleBasedModeNarrator to TypeScript for web parity.
// Analyzes binding modes from a BindingPopulation and generates
// plain-English characterizations with optimization hints.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import { allowsBindingClaims } from '@bonhomme/shared';
import type { BindingPopulation } from '@bonhomme/shared';
import type { ModeDescription, BindingModeNarrative } from '@bonhomme/shared';

/** Internal mode profile for analysis. */
interface ModeProfile {
  index: number;
  poseCount: number;
  freeEnergy: number;
  entropy: number;
  enthalpy: number;
  boltzmannWeight: number;
  isEntropyDriven: boolean;
  isEnthalpyDriven: boolean;
}

const kBkcal = 0.001987206;

export class BindingModeAnalyzer {
  /**
   * Analyze binding modes from a population and generate narrative.
   */
  static analyze(population: BindingPopulation): BindingModeNarrative {
    const profiles = BindingModeAnalyzer.buildProfiles(population);
    const bindingPhysical = allowsBindingClaims(population.globalThermodynamics);
    const descriptions = profiles.map((mode) => BindingModeAnalyzer.describeMode(mode, bindingPhysical));

    // Find dominant mode (highest Boltzmann weight)
    const dominantIdx = profiles.reduce(
      (best, mode, i) => (mode.boltzmannWeight > profiles[best].boltzmannWeight ? i : best),
      0,
    );

    const dominant = profiles[dominantIdx];
    let selectivityInsight: string;
    if (dominant) {
      selectivityInsight = bindingPhysical
        ? `Mode ${dominant.index + 1} dominates (${(dominant.boltzmannWeight * 100).toFixed(0)}% weight). Focus SAR optimization on this binding geometry.`
        : `Mode ${dominant.index + 1} dominates the ensemble/CF diagnostic (${(dominant.boltzmannWeight * 100).toFixed(0)}% diagnostic weight). Physical affinity unavailable; retain this ordering only as a pose-ensemble diagnostic.`;
    } else {
      selectivityInsight = bindingPhysical
        ? 'No dominant mode — population is diverse across geometries.'
        : 'No dominant mode in the ensemble/CF diagnostic; physical affinity unavailable.';
    }

    const confidence = profiles.length >= 2 ? 0.8 : 0.5;
    return { modeDescriptions: descriptions, selectivityInsight, confidence };
  }

  private static buildProfiles(population: BindingPopulation): ModeProfile[] {
    const T = population.temperature;
    // Compute Boltzmann weights from free energies
    const modes = population.modes.slice(0, 3); // top 3 by free energy
    const totalWeight = modes.reduce((sum, m) => sum + Math.exp(-m.freeEnergy / (kBkcal * T)), 0);

    return modes.map((mode, i) => {
      const weight = totalWeight > 0 ? Math.exp(-mode.freeEnergy / (kBkcal * T)) / totalWeight : 0;
      const tsContribution = Math.abs(T * mode.entropy);
      const hContribution = Math.abs(mode.enthalpy);
      return {
        index: i,
        poseCount: mode.size,
        freeEnergy: mode.freeEnergy,
        entropy: mode.entropy,
        enthalpy: mode.enthalpy,
        boltzmannWeight: weight,
        isEntropyDriven: tsContribution > hContribution * 1.5,
        isEnthalpyDriven: hContribution > tsContribution * 1.5,
      };
    });
  }

  private static describeMode(mode: ModeProfile, bindingPhysical: boolean): ModeDescription {
    // Driving force
    let driving: string;
    if (mode.isEntropyDriven) {
      driving = 'entropy-driven (conformational flexibility stabilizes this geometry)';
    } else if (mode.isEnthalpyDriven) {
      driving = 'enthalpy-driven (direct interactions dominate stability)';
    } else {
      driving = 'balanced enthalpy-entropy contributions';
    }

    // Tightness
    const tightness = mode.poseCount < 5 ? 'tight cluster' : mode.poseCount < 20 ? 'moderate cluster' : 'broad ensemble';

    if (!bindingPhysical) {
      const diagnosticDriver = mode.isEntropyDriven
        ? 'high entropy-slot influence'
        : mode.isEnthalpyDriven
          ? 'high direct-score influence'
          : 'balanced score/entropy slots';
      return {
        characterization: `Mode ${mode.index + 1}: ${tightness} (${mode.poseCount} poses), ensemble/CF diagnostic = ${mode.freeEnergy.toFixed(1)} (source units), ${diagnosticDriver}. Diagnostic weight ${(mode.boltzmannWeight * 100).toFixed(0)}%.`,
        optimizationHint: 'Physical affinity unavailable; validate calibrated energy, ensemble measure, and matched association-cycle artifacts before potency or selectivity claims.',
      };
    }

    const characterization = `Mode ${mode.index + 1}: ${tightness} (${mode.poseCount} poses), F = ${mode.freeEnergy.toFixed(1)} kcal/mol, ${driving}. Boltzmann weight ${(mode.boltzmannWeight * 100).toFixed(0)}%.`;

    // Optimization hint
    let optimizationHint: string;
    if (mode.isEntropyDriven) {
      optimizationHint = 'Rigidify ligand to shift from entropy- to enthalpy-driven binding for selectivity.';
    } else if (mode.isEnthalpyDriven) {
      optimizationHint = 'Add flexible substituents to improve entropy contribution if affinity needs boosting.';
    } else {
      optimizationHint = 'Well-balanced mode — optimize interaction geometry for potency.';
    }

    return { characterization, optimizationHint };
  }
}
