// resultLoader.ts — Parse FlexAID output files into TypeScript types
//
// Port of python/flexaidds/results_io.py for browser/Node.js usage.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import {
  PROXY_ONLY_PROVENANCE,
  normalizeDockingRecord,
  type DockingResult,
  type BindingModeResult,
  type NormalizedThermodynamicRecord,
  type ThermodynamicResult,
} from './types.js';

/**
 * Parse a FlexAID .rrd (result) file into a DockingResult.
 *
 * File format (pipe-delimited):
 * ```
 * CLUSTER|1|CF|-10.50|POSES|25|RMSD|1.23
 * POSE|1|CF|-10.50|RANK|1|RMSD|0.00
 * ...
 * ```
 */
export function parseResultFile(content: string): DockingResult {
  const lines = content.split('\n').filter((l) => l.trim().length > 0);
  const modes: BindingModeResult[] = [];
  let currentMode: Partial<BindingModeResult> | null = null;
  let poseCount = 0;

  for (const line of lines) {
    const parts = line.split('|').map((p) => p.trim());

    if (parts[0] === 'CLUSTER') {
      // Save previous mode
      if (currentMode) {
        modes.push(finishMode(currentMode, poseCount));
      }

      currentMode = {};
      poseCount = 0;

      // Parse cluster fields
      for (let i = 1; i < parts.length - 1; i += 2) {
        const key = parts[i];
        const val = parts[i + 1];
        if (key === 'CF') currentMode.freeEnergy = parseFloat(val);
        if (key === 'POSES') poseCount = parseInt(val, 10);
      }
    } else if (parts[0] === 'POSE') {
      poseCount++;
    }
  }

  // Save last mode
  if (currentMode) {
    modes.push(finishMode(currentMode, poseCount));
  }

  // Compute global thermodynamics from all modes
  const globalThermo = computeGlobalThermo(modes);

  return {
    bindingModes: modes,
    globalThermodynamics: globalThermo,
    // 300 K is a hard-coded default, not a value read from the file. The
    // numeric default is preserved for wire compatibility; temperatureKnown
    // is the metadata that stops a consumer treating it as measured.
    temperature: 300.0,
    temperatureKnown: false,
    populationSize: modes.reduce((s, m) => s + m.size, 0),
    timestamp: new Date().toISOString(),
  };
}

type UnknownRecord = Record<string, unknown>;

/**
 * Fill canonical camelCase numeric slots from the shared normalizer without
 * ever overwriting a value the producer already supplied in camelCase.
 *
 * Producer-supplied numbers therefore pass through byte-identical; only a
 * foreign dialect (Python snake_case / unit-suffixed keys) gets re-keyed.
 */
function adoptNumericSlots<T extends UnknownRecord>(
  raw: UnknownRecord,
  canonical: UnknownRecord,
  keys: readonly string[],
): T {
  const merged: UnknownRecord = { ...raw };
  for (const key of keys) {
    if (typeof merged[key] === 'number' && Number.isFinite(merged[key])) continue;
    merged[key] = canonical[key];
  }
  return merged as T;
}

const THERMO_SLOTS = [
  'temperature', 'logZ', 'freeEnergy', 'meanEnergy',
  'meanEnergySq', 'heatCapacity', 'entropy', 'stdEnergy',
] as const;

const MODE_SLOTS = ['size', 'freeEnergy', 'entropy', 'enthalpy', 'heatCapacity'] as const;

function adoptThermodynamicRecord(
  raw: unknown,
  canonical: NormalizedThermodynamicRecord,
): ThermodynamicResult {
  const source = (raw && typeof raw === 'object' ? raw : {}) as UnknownRecord;
  const adopted = adoptNumericSlots<ThermodynamicResult & UnknownRecord>(
    source,
    canonical as unknown as UnknownRecord,
    THERMO_SLOTS,
  );
  // Availability and provenance are decided by the shared firewall, never by
  // whatever the producer asserted about itself.
  adopted.available = canonical.available;
  if (!canonical.available && adopted.unavailableReason === undefined) {
    adopted.unavailableReason = canonical.unavailableReason;
  }
  adopted.scientificProvenance = canonical.scientificProvenance;
  delete adopted.scientific_provenance;
  return adopted;
}

/**
 * Parse a JSON docking result from any producer dialect.
 *
 * Accepts SDK camelCase `bindingModes`, shared/viewer `modes`, and Python
 * snake_case `binding_modes` with unit-suffixed thermodynamic keys. All of
 * them go through the single shared wire normalizer, so availability is
 * always decided deliberately and provenance is never taken on trust.
 */
export function parseDockingResultJSON(json: string): DockingResult {
  return adoptDockingPayload(JSON.parse(json));
}

/** Runtime-boundary adoption of an already-parsed payload. */
export function adoptDockingPayload(payload: unknown): DockingResult {
  const raw = (payload && typeof payload === 'object' ? payload : {}) as UnknownRecord;
  const canonical = normalizeDockingRecord(raw);

  const rawModes: unknown[] = [raw.bindingModes, raw.modes, raw.binding_modes]
    .find((candidate): candidate is unknown[] => Array.isArray(candidate)) ?? [];

  const bindingModes = canonical.modes.map((mode, i) => {
    const source = (rawModes[i] && typeof rawModes[i] === 'object'
      ? rawModes[i]
      : {}) as UnknownRecord;
    const adopted = adoptNumericSlots<BindingModeResult & UnknownRecord>(
      source,
      mode as unknown as UnknownRecord,
      MODE_SLOTS,
    );
    adopted.thermodynamicsAvailable = mode.thermodynamics.available;
    adopted.thermodynamics = adoptThermodynamicRecord(
      source.thermodynamics,
      mode.thermodynamics,
    );
    return adopted;
  });

  const result: DockingResult & UnknownRecord = {
    ...raw,
    bindingModes,
    globalThermodynamics: adoptThermodynamicRecord(
      raw.globalThermodynamics ?? raw.global_thermodynamics,
      canonical.globalThermodynamics,
    ),
    temperature: canonical.temperature,
    temperatureKnown: canonical.temperatureKnown,
    populationSize: typeof raw.populationSize === 'number'
      ? raw.populationSize
      : canonical.totalPoses,
    timestamp: typeof raw.timestamp === 'string' ? raw.timestamp : '',
  };
  delete result.modes;
  delete result.binding_modes;
  delete result.global_thermodynamics;
  return result;
}

function finishMode(partial: Partial<BindingModeResult>, poseCount: number): BindingModeResult {
  return {
    size: poseCount || 1,
    freeEnergy: partial.freeEnergy ?? 0,
    thermodynamicsAvailable: false,
    scoreDomain: 'cf_arbitrary_units',
    entropy: partial.entropy ?? 0,
    enthalpy: partial.enthalpy ?? partial.freeEnergy ?? 0,
    heatCapacity: partial.heatCapacity ?? 0,
    thermodynamics: partial.thermodynamics,
  };
}

function computeGlobalThermo(modes: BindingModeResult[]): ThermodynamicResult {
  if (modes.length === 0) {
    return {
      available: false,
      unavailableReason: 'raw_result_contains_no_thermodynamic_record',
      scientificProvenance: PROXY_ONLY_PROVENANCE,
      temperature: 300, logZ: 0, freeEnergy: 0,
      meanEnergy: 0, meanEnergySq: 0,
      heatCapacity: 0, entropy: 0, stdEnergy: 0,
    };
  }

  // Preserve legacy numeric placeholders for wire compatibility, but mark them
  // explicitly unavailable. A list of cluster CF scores is not a partition
  // function and cannot be promoted to physical thermodynamics.
  const totalSize = modes.reduce((s, m) => s + m.size, 0);
  const weightedEnergy = modes.reduce((s, m) => s + m.freeEnergy * m.size, 0) / totalSize;

  return {
    available: false,
    unavailableReason: 'raw_result_contains_cf_scores_only',
    scientificProvenance: {
      ...PROXY_ONLY_PROVENANCE,
      energyDomain: 'cf_arbitrary_units',
      ensembleMeasure: 'optimizer_samples',
      energyProvenance: 'FlexAID result-file CF field',
      measureProvenance: 'clustered GA result records',
    },
    temperature: 300,
    logZ: 0,
    freeEnergy: Math.min(...modes.map((m) => m.freeEnergy)),
    meanEnergy: weightedEnergy,
    meanEnergySq: 0,
    heatCapacity: modes.reduce((s, m) => s + m.heatCapacity, 0) / modes.length,
    entropy: modes.reduce((s, m) => s + m.entropy, 0) / modes.length,
    stdEnergy: 0,
  };
}
