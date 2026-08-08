// CleftAnalyzer.ts — Rule-based binding cleft druggability assessor
//
// Ports the Swift RuleBasedCleftAssessor to TypeScript for web parity.
// Evaluates pocket geometry to determine druggability tier
// and suggest ligand properties.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import { allowsBindingClaims, claimValidityForRecord } from '@bonhomme/shared';
import type { CleftFeatures, ThermodynamicClaimSource } from '@bonhomme/shared';
import type { DruggabilityTier, CleftAssessment } from '@bonhomme/shared';

export class CleftAnalyzer {
  /**
   * Score a binding cleft from pocket geometry.
   *
   * "Druggable" and "suggested ligand" are potency claims about a molecule
   * that has not been measured, so they are only spoken when `claimSource`
   * carries binding-physical evidence. Without it the same geometric score is
   * reported as a shape diagnostic. The threshold arithmetic is identical in
   * both cases — only the language changes.
   */
  static assess(
    cleft: CleftFeatures,
    claimSource?: ThermodynamicClaimSource,
  ): CleftAssessment {
    const warnings: string[] = [];
    let score = 0.0;

    // Volume scoring (200-1000 A^3 is ideal for drug-like molecules)
    if (cleft.volume >= 200 && cleft.volume <= 1000) {
      score += 0.3;
    } else if (cleft.volume < 200) {
      warnings.push(`Pocket too small (${cleft.volume.toFixed(0)} A^3) for drug-like molecules.`);
    } else if (cleft.volume > 1200) {
      warnings.push(`Pocket very large (${cleft.volume.toFixed(0)} A^3) — may lack specificity.`);
    } else {
      score += 0.15;
    }

    // Hydrophobicity scoring (40-80% is drug-like)
    if (cleft.hydrophobicFraction >= 0.4 && cleft.hydrophobicFraction <= 0.8) {
      score += 0.25;
    } else if (cleft.hydrophobicFraction > 0.8) {
      score += 0.15;
      warnings.push('Heavily hydrophobic pocket — limited polar interaction sites.');
    } else {
      score += 0.1;
    }

    // Depth scoring (deeper pockets bind better)
    if (cleft.depth > 6.0) {
      score += 0.2;
    } else if (cleft.depth > 3.0) {
      score += 0.1;
    } else {
      warnings.push(`Shallow pocket (${cleft.depth.toFixed(1)} A) — geometry may not constrain docked poses.`);
    }

    // Solvent exposure (lower generally gives a more enclosed docking cleft)
    if (cleft.solventExposure < 0.3) {
      score += 0.15;
    } else if (cleft.solventExposure > 0.6) {
      warnings.push(`High solvent exposure (${(cleft.solventExposure * 100).toFixed(0)}%) — desolvation penalty.`);
    } else {
      score += 0.05;
    }

    // Anchor residues
    if (cleft.anchorResidueCount >= 4) {
      score += 0.1;
    }

    // Tier determination
    let druggability: DruggabilityTier;
    if (score >= 0.7) {
      druggability = 'high';
    } else if (score >= 0.45) {
      druggability = 'moderate';
    } else if (score >= 0.2) {
      druggability = 'low';
    } else {
      druggability = 'undruggable';
    }

    const bindingPhysical = allowsBindingClaims(claimSource);
    const claimValidity = claimValidityForRecord(claimSource);

    // Summary
    const volumeDesc = cleft.volume < 400 ? 'Small' : cleft.volume < 800 ? 'Medium' : 'Large';
    const shapeDesc = cleft.elongation < 0.3 ? 'spherical' : cleft.elongation < 0.7 ? 'oval' : 'elongated';
    const geometry = `${volumeDesc} ${shapeDesc} pocket (${cleft.volume.toFixed(0)} A^3, ${(cleft.hydrophobicFraction * 100).toFixed(0)}% hydrophobic, ${cleft.anchorResidueCount} anchor residues).`;
    const summary = bindingPhysical
      ? `${geometry} Druggability: ${druggability}.`
      : `${geometry} Geometric enclosure score: ${druggability} (shape diagnostic, source units). Druggability and physical affinity are unavailable without binding-physical provenance.`;

    // Suggested ligand properties
    let suggestedLigandProperties: string;
    if (!bindingPhysical) {
      suggestedLigandProperties =
        'Unavailable: ligand-property guidance is a potency claim and requires '
        + 'binding-physical provenance. Pocket shape and polarity above are '
        + 'geometric diagnostics only.';
    } else if (cleft.hydrophobicFraction > 0.6) {
      suggestedLigandProperties = 'Lipophilic compounds, MW 300-500 Da, LogP > 2';
    } else if (cleft.hydrophobicFraction < 0.4) {
      suggestedLigandProperties = 'Polar compounds with H-bond donors/acceptors, MW 200-400 Da';
    } else {
      suggestedLigandProperties = 'Balanced compounds, MW 300-500 Da, 2-3 H-bond acceptors';
    }

    return {
      druggability,
      summary,
      suggestedLigandProperties,
      warnings,
      claimValidity,
      claimSource,
    };
  }
}
