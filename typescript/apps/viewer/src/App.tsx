// App.tsx — BonhommeViewer main application component
//
// PWA for visualizing FlexAIDdS binding population results.
// Uses Mol* for 3D molecular rendering.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import React, { useState } from 'react';
import type { BindingPopulation, HealthCorrelation } from '@bonhomme/shared';
import { deserializePopulation } from '@bonhomme/shared';
import { IntelligenceEngine } from './IntelligenceEngine';
import { FleetDashboard } from './FleetDashboard';
import { MolstarViewer as MolstarViewerComponent, populationColorHex } from './MolstarViewer';
import { populationAllowsPhysicalAffinity } from './claimPresentation';

type View = 'population' | 'fleet' | 'health';

export function App() {
  const [population, setPopulation] = useState<BindingPopulation | null>(null);
  const [health, setHealth] = useState<HealthCorrelation | null>(null);
  const [view, setView] = useState<View>('population');
  const [oracleAnalysis, setOracleAnalysis] = useState<string[]>([]);

  // Load population from file or iCloud watcher
  const handleFileLoad = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;
    const text = await file.text();
    try {
      const pop = deserializePopulation(text);
      setPopulation(pop);
      // Auto-analyze with oracle
      const analysis = await IntelligenceEngine.analyze(pop, health ?? undefined);
      setOracleAnalysis(analysis.bullets);
    } catch {
      console.error('Failed to parse population file');
    }
  };

  return (
    <div style={{ fontFamily: 'system-ui', padding: '1rem' }}>
      <header style={{ display: 'flex', gap: '1rem', alignItems: 'center', marginBottom: '1rem' }}>
        <h1 style={{ margin: 0 }}>BonhommeViewer</h1>
        <nav style={{ display: 'flex', gap: '0.5rem' }}>
          <button onClick={() => setView('population')} disabled={view === 'population'}>
            Population
          </button>
          <button onClick={() => setView('fleet')} disabled={view === 'fleet'}>
            Fleet
          </button>
          <button onClick={() => setView('health')} disabled={view === 'health'}>
            Health
          </button>
        </nav>
        {/* JSON only: the viewer has no .rrd parser, and advertising one
            would invite users to load a file it cannot interpret. */}
        <input type="file" accept=".json" onChange={handleFileLoad} />
      </header>

      {view === 'population' && population && (
        <PopulationView population={population} oracleAnalysis={oracleAnalysis} />
      )}

      {view === 'fleet' && <FleetDashboard />}

      {view === 'health' && health && (
        <HealthView health={health} />
      )}

      {!population && view === 'population' && (
        <div style={{ textAlign: 'center', padding: '4rem', color: '#666' }}>
          <p>Load a population JSON file or connect to iCloud Drive watcher.</p>
        </div>
      )}
    </div>
  );
}

function PopulationView({ population, oracleAnalysis }: {
  population: BindingPopulation;
  oracleAnalysis: string[];
}) {
  const bindingPhysical = populationAllowsPhysicalAffinity(population);

  return (
    <div>
      <section>
        <h2>{bindingPhysical ? 'Global Thermodynamics' : 'Ensemble / CF Diagnostics'}</h2>
        <table>
          <tbody>
            <tr><td>Temperature</td><td>{population.temperature} K</td></tr>
            <tr>
              <td>{bindingPhysical ? 'Free Energy' : 'Ensemble/CF diagnostic'}</td>
              <td>
                {population.globalThermodynamics.freeEnergy.toFixed(3)}{' '}
                {bindingPhysical ? 'kcal/mol' : 'source units'}
              </td>
            </tr>
            <tr>
              <td>{bindingPhysical ? 'Shannon S' : 'Shannon diagnostic'}</td>
              <td>
                {population.shannonS.toFixed(6)}{' '}
                {bindingPhysical ? 'kcal/mol/K' : 'source units'}
              </td>
            </tr>
            {!bindingPhysical && (
              <tr><td>Claim status</td><td>Physical affinity unavailable</td></tr>
            )}
            <tr><td>Entropy Collapsed</td><td>{population.isCollapsed ? 'YES' : 'No'}</td></tr>
            <tr><td>Binding Modes</td><td>{population.modes.length}</td></tr>
            <tr><td>Total Poses</td><td>{population.totalPoses}</td></tr>
          </tbody>
        </table>
      </section>

      {oracleAnalysis.length > 0 && (
        <section>
          <h2>Intelligence Oracle</h2>
          <ul>
            {oracleAnalysis.map((bullet, i) => (
              <li key={i}>{bullet}</li>
            ))}
          </ul>
        </section>
      )}

      <section>
        <h2>Binding Modes</h2>
        <table>
          <thead>
            <tr>
              <th>#</th>
              <th>Poses</th>
              <th>{bindingPhysical ? 'F (kcal/mol)' : 'Ensemble/CF diagnostic (source units)'}</th>
              <th>{bindingPhysical ? 'S (kcal/mol/K)' : 'Entropy diagnostic (source units)'}</th>
              <th>Cv</th>
              <th>Probability</th>
            </tr>
          </thead>
          <tbody>
            {population.modes.map((mode, i) => (
              <tr key={i}>
                <td>{i + 1}</td>
                <td>{mode.size}</td>
                <td>{mode.freeEnergy.toFixed(3)}</td>
                <td>{mode.entropy.toFixed(6)}</td>
                <td>{mode.heatCapacity.toFixed(6)}</td>
                <td style={{ color: populationColorHex(mode.probability) }}>
                  {(mode.probability * 100).toFixed(1)}%
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </section>

      {population.targetModifications && population.targetModifications.length > 0 && (
        <section>
          <h2>Target PTM / Glycan Modifications</h2>
          <ul>
            {population.targetModifications.map((mod, i) => (
              <li key={i}>
                {mod.type} at {mod.residueName}{mod.residueNumber} (chain {mod.chainID})
                {mod.composition && ` — ${mod.composition}`}
                {mod.effect && (bindingPhysical
                  ? ` — ΔF=${mod.effect.deltaFreeEnergy.toFixed(2)} kcal/mol, ΔS=${mod.effect.deltaEntropy.toFixed(4)} kcal/mol/K`
                  : ` — ensemble/CF diagnostic shift=${mod.effect.deltaFreeEnergy.toFixed(2)}, entropy diagnostic=${mod.effect.deltaEntropy.toFixed(4)} (source units); physical affinity unavailable`)}
              </li>
            ))}
          </ul>
        </section>
      )}

      <section>
        <h2>3D Viewer</h2>
        <MolstarViewerComponent population={population} />
      </section>
    </div>
  );
}

function HealthView({ health }: { health: HealthCorrelation }) {
  return (
    <div>
      <h2>Health Correlation</h2>
      <table>
        <tbody>
          <tr><td>Shannon S</td><td>{health.shannonS.toFixed(6)}</td></tr>
          {health.hrvSDNN && <tr><td>HRV SDNN</td><td>{health.hrvSDNN.toFixed(0)} ms</td></tr>}
          {health.sleepHours && <tr><td>Sleep</td><td>{health.sleepHours.toFixed(1)} hours</td></tr>}
          {health.restingHeartRate && <tr><td>Resting HR</td><td>{health.restingHeartRate.toFixed(0)} bpm</td></tr>}
        </tbody>
      </table>

      {health.oracleAnalysis && (
        <section>
          <h3>Oracle Analysis</h3>
          <ul>
            {health.oracleAnalysis.map((b, i) => <li key={i}>{b}</li>)}
          </ul>
        </section>
      )}

      {health.fitnessRecommendation && (
        <section>
          <h3>Fitness Recommendation: {health.fitnessRecommendation.intensity}</h3>
          <ul>
            {health.fitnessRecommendation.activities.map((a, i) => <li key={i}>{a}</li>)}
          </ul>
        </section>
      )}
    </div>
  );
}
