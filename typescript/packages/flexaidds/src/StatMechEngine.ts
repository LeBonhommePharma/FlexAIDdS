// StatMechEngine.ts — pure-TypeScript statistical-mechanics diagnostics
//
// Provides a TypeScript-native API matching the Swift FlexAIDRunner.
//
// WASM IS NOT SUPPORTED. There is no compiled-C++ path behind this class and
// nothing here is provenance-equivalent to LIB/statmech.cpp: results carry
// whatever provenance the caller supplies, which for current FlexAID output is
// proxy-only. Do not advertise a WASM/native mode until one exists and can
// emit sha256 artifact receipts.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import {
  PROXY_ONLY_PROVENANCE,
  type ScientificProvenance,
  type ThermodynamicResult,
  type TIPoint,
  type WHAMBinResult,
} from './types.js';
import { kB_kcal } from './constants.js';

interface State {
  energy: number;
  count: number;
}

/**
 * Pure TypeScript implementation of the StatMechEngine.
 *
 * Computes full thermodynamic properties from a conformational ensemble
 * using Boltzmann statistics with log-sum-exp numerical stability.
 *
 * Usage:
 * ```ts
 * const engine = new StatMechEngine(300.0);
 * engine.addSample(-10.0);
 * engine.addSample(-8.0, 2);
 * const result = engine.compute();
 * console.log(`Ensemble F-like value: ${result.freeEnergy}`);
 * ```
 */
export class StatMechEngine {
  private readonly T: number;
  private readonly beta: number;
  private ensemble: State[] = [];
  private provenanceValue: ScientificProvenance;

  constructor(
    temperatureK: number = 300.0,
    provenance: ScientificProvenance = PROXY_ONLY_PROVENANCE,
  ) {
    this.T = temperatureK;
    this.beta = 1.0 / (kB_kcal * temperatureK);
    this.provenanceValue = { ...provenance };
  }

  /** Scientific interpretation attached to computed results. */
  get provenance(): ScientificProvenance {
    return { ...this.provenanceValue };
  }

  /** Replace interpretation metadata without changing numerical state. */
  setProvenance(provenance: ScientificProvenance): void {
    this.provenanceValue = { ...provenance };
  }

  /** Add a sampled conformation to the ensemble. */
  addSample(energy: number, multiplicity: number = 1): void {
    this.ensemble.push({ energy, count: multiplicity });
  }

  /** Remove all samples. */
  clear(): void {
    this.ensemble = [];
  }

  /** Number of states in the ensemble. */
  get size(): number {
    return this.ensemble.length;
  }

  /** Temperature in Kelvin. */
  get temperature(): number {
    return this.T;
  }

  /** Compute full thermodynamic properties. */
  compute(): ThermodynamicResult {
    if (this.ensemble.length === 0) {
      return {
        available: false,
        unavailableReason: 'empty_ensemble',
        scientificProvenance: { ...this.provenanceValue },
        temperature: this.T, logZ: 0, freeEnergy: 0,
        meanEnergy: 0, meanEnergySq: 0,
        heatCapacity: 0, entropy: 0, stdEnergy: 0,
      };
    }

    // Log-sum-exp for numerical stability
    const logTerms: number[] = [];
    for (const state of this.ensemble) {
      const logTerm = -this.beta * state.energy + Math.log(state.count);
      logTerms.push(logTerm);
    }

    const logZ = this.logSumExp(logTerms);
    const freeEnergy = -this.T * kB_kcal * logZ;

    // Boltzmann weights and expectation values
    let meanE = 0;
    let meanE2 = 0;
    for (let i = 0; i < this.ensemble.length; i++) {
      const w = Math.exp(logTerms[i] - logZ);
      meanE += w * this.ensemble[i].energy;
      meanE2 += w * this.ensemble[i].energy * this.ensemble[i].energy;
    }

    const kT2 = kB_kcal * this.T * this.T;
    const variance = meanE2 - meanE * meanE;
    const heatCapacity = variance / kT2;
    const entropy = (meanE - freeEnergy) / this.T;
    const stdEnergy = Math.sqrt(Math.max(0, variance));

    return {
      available: true,
      scientificProvenance: { ...this.provenanceValue },
      temperature: this.T,
      logZ,
      freeEnergy,
      meanEnergy: meanE,
      meanEnergySq: meanE2,
      heatCapacity,
      entropy,
      stdEnergy,
    };
  }

  /** Get normalized Boltzmann weights for all samples. */
  boltzmannWeights(): number[] {
    if (this.ensemble.length === 0) return [];

    const logTerms = this.ensemble.map(
      (s) => -this.beta * s.energy + Math.log(s.count),
    );
    const logZ = this.logSumExp(logTerms);
    return logTerms.map((lt) => Math.exp(lt - logZ));
  }

  /** Helmholtz free energy from a raw energy array (static). */
  static helmholtz(energies: number[], temperatureK: number): number {
    const beta = 1.0 / (kB_kcal * temperatureK);
    const logTerms = energies.map((e) => -beta * e);
    const maxTerm = Math.max(...logTerms);
    const logZ = maxTerm + Math.log(logTerms.reduce((s, t) => s + Math.exp(t - maxTerm), 0));
    return -temperatureK * kB_kcal * logZ;
  }

  /** Thermodynamic integration via trapezoidal rule (static). */
  static thermodynamicIntegration(points: TIPoint[]): number {
    if (points.length < 2) return 0;
    let integral = 0;
    for (let i = 1; i < points.length; i++) {
      const dLambda = points[i].lambda - points[i - 1].lambda;
      integral += 0.5 * (points[i].dVdLambda + points[i - 1].dVdLambda) * dLambda;
    }
    return integral;
  }

  /** Numerically stable log(sum(exp(x))). */
  private logSumExp(x: number[]): number {
    const maxVal = Math.max(...x);
    if (!isFinite(maxVal)) return -Infinity;
    const sumExp = x.reduce((s, xi) => s + Math.exp(xi - maxVal), 0);
    return maxVal + Math.log(sumExp);
  }
}
