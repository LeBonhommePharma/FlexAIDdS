// FleetDashboard.tsx — Real-time fleet status dashboard
//
// Displays device status, active jobs, chunk completion, orphan recovery,
// battery-aware scheduling indicators, and per-device telemetry.
// Polls fleet status from iCloud watcher JSON or WebSocket relay.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import React, { useState, useEffect, useCallback, useRef } from 'react';
import type { DeviceCapability, WorkChunk } from '@bonhomme/flexaidds';
import type { RefereeVerdict } from '@bonhomme/shared';
import { IntelligencePanel } from './IntelligencePanel.js';

interface FleetMetrics {
  jobID: string;
  totalChunks: number;
  completedChunks: number;
  failedChunks: number;
  orphanedChunks: number;
  activeDevices: number;
  totalTFLOPS: number;
  meanChunkTimeSeconds: number | null;
  estimatedRemainingSeconds: number | null;
  timestamp: string;
}

interface FleetState {
  devices: DeviceCapability[];
  activeChunks: WorkChunk[];
  completedChunks: number;
  totalChunks: number;
  failedChunks: number;
  orphanedChunks: number;
  metrics: FleetMetrics | null;
  refereeVerdict: RefereeVerdict | null;
  connectionStatus: 'connected' | 'polling' | 'disconnected';
  lastUpdate: Date | null;
}

interface FleetDashboardProps {
  /** URL to poll for fleet status JSON (e.g., iCloud shared status file) */
  statusURL?: string;
  /** WebSocket URL for real-time updates (preferred over polling) */
  wsURL?: string;
  /** Polling interval in milliseconds (default: 5000) */
  pollIntervalMs?: number;
}

export function FleetDashboard({
  statusURL,
  wsURL,
  pollIntervalMs = 5000,
}: FleetDashboardProps) {
  const [fleet, setFleet] = useState<FleetState>({
    devices: [],
    activeChunks: [],
    completedChunks: 0,
    totalChunks: 0,
    failedChunks: 0,
    orphanedChunks: 0,
    metrics: null,
    refereeVerdict: null,
    connectionStatus: 'disconnected',
    lastUpdate: null,
  });

  const wsRef = useRef<WebSocket | null>(null);

  // Parse fleet status from JSON (shared format between polling and WebSocket)
  const applyFleetUpdate = useCallback((data: Record<string, unknown>) => {
    const devices = (data.devices ?? []) as DeviceCapability[];
    const chunks = (data.activeChunks ?? []) as WorkChunk[];
    const metrics = (data.metrics ?? null) as FleetMetrics | null;
    const refereeVerdict = (data.refereeVerdict ?? null) as RefereeVerdict | null;

    setFleet((prev) => ({
      ...prev,
      devices,
      activeChunks: chunks,
      completedChunks: metrics?.completedChunks ?? (data.completedChunks as number) ?? prev.completedChunks,
      totalChunks: metrics?.totalChunks ?? (data.totalChunks as number) ?? prev.totalChunks,
      failedChunks: metrics?.failedChunks ?? 0,
      orphanedChunks: metrics?.orphanedChunks ?? 0,
      metrics,
      refereeVerdict,
      lastUpdate: new Date(),
    }));
  }, []);

  // WebSocket connection (preferred: real-time updates)
  useEffect(() => {
    if (!wsURL) return;

    const ws = new WebSocket(wsURL);
    wsRef.current = ws;

    ws.onopen = () => {
      setFleet((prev) => ({ ...prev, connectionStatus: 'connected' }));
    };

    ws.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data);
        applyFleetUpdate(data);
      } catch {
        // Ignore malformed messages
      }
    };

    ws.onclose = () => {
      setFleet((prev) => ({ ...prev, connectionStatus: 'disconnected' }));
    };

    return () => {
      ws.close();
      wsRef.current = null;
    };
  }, [wsURL, applyFleetUpdate]);

  // HTTP polling fallback (when no WebSocket)
  useEffect(() => {
    if (wsURL || !statusURL) return;

    setFleet((prev) => ({ ...prev, connectionStatus: 'polling' }));

    const poll = async () => {
      try {
        const res = await fetch(statusURL);
        if (res.ok) {
          const data = await res.json();
          applyFleetUpdate(data);
        }
      } catch {
        // Network error — fleet may be offline
      }
    };

    poll(); // Initial fetch
    const interval = setInterval(poll, pollIntervalMs);
    return () => clearInterval(interval);
  }, [statusURL, wsURL, pollIntervalMs, applyFleetUpdate]);

  const thermalColor = (state: string) => {
    switch (state) {
      case 'nominal': return '#4B0082';   // purple blue (cool)
      case 'fair': return '#6B0060';      // mid-purple
      case 'serious': return '#700030';   // mid-burgundy
      case 'critical': return '#800020';  // burgundy red (hot)
      default: return '#9ca3af';
    }
  };

  const batteryIcon = (level: number | undefined, charging: boolean | undefined) => {
    if (level === undefined) return null; // Desktop Mac
    const pct = Math.round(level * 100);
    const color = pct < 20 ? '#ef4444' : pct < 50 ? '#fbbf24' : '#4ade80';
    return (
      <span style={{ color }}>
        {pct}%{charging ? ' (charging)' : ''}
      </span>
    );
  };

  const formatTime = (seconds: number | null | undefined) => {
    if (!seconds) return '--';
    if (seconds < 60) return `${Math.round(seconds)}s`;
    if (seconds < 3600) return `${Math.round(seconds / 60)}m`;
    return `${(seconds / 3600).toFixed(1)}h`;
  };

  const progressPct = fleet.totalChunks > 0
    ? (fleet.completedChunks / fleet.totalChunks) * 100
    : 0;

  return (
    <div>
      <h2>Bonhomme Fleet</h2>

      {/* Connection status */}
      <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginBottom: '1rem' }}>
        <span style={{
          width: '10px', height: '10px', borderRadius: '50%',
          backgroundColor: fleet.connectionStatus === 'connected' ? '#4ade80'
            : fleet.connectionStatus === 'polling' ? '#fbbf24' : '#ef4444',
          display: 'inline-block',
        }} />
        <span style={{ fontSize: '0.85rem', color: '#999' }}>
          {fleet.connectionStatus === 'connected' ? 'Live (WebSocket)'
            : fleet.connectionStatus === 'polling' ? 'Polling'
            : 'Disconnected'}
          {fleet.lastUpdate && ` — updated ${fleet.lastUpdate.toLocaleTimeString()}`}
        </span>
      </div>

      {/* Fleet summary metrics */}
      {fleet.metrics && (
        <section style={{
          display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(140px, 1fr))',
          gap: '0.75rem', marginBottom: '1.5rem',
        }}>
          {[
            { label: 'Devices', value: fleet.metrics.activeDevices.toString() },
            { label: 'TFLOPS', value: fleet.metrics.totalTFLOPS.toFixed(1) },
            { label: 'Completed', value: `${fleet.metrics.completedChunks}/${fleet.metrics.totalChunks}` },
            { label: 'Failed', value: fleet.metrics.failedChunks.toString(), color: fleet.metrics.failedChunks > 0 ? '#ef4444' : undefined },
            { label: 'Orphaned', value: fleet.metrics.orphanedChunks.toString(), color: fleet.metrics.orphanedChunks > 0 ? '#f97316' : undefined },
            { label: 'Avg Time', value: formatTime(fleet.metrics.meanChunkTimeSeconds) },
            { label: 'ETA', value: formatTime(fleet.metrics.estimatedRemainingSeconds) },
          ].map((item) => (
            <div key={item.label} style={{
              background: '#1a1a2e', borderRadius: '8px', padding: '0.75rem', textAlign: 'center',
            }}>
              <div style={{ fontSize: '0.75rem', color: '#666' }}>{item.label}</div>
              <div style={{ fontSize: '1.25rem', fontWeight: 'bold', color: item.color ?? '#fff' }}>{item.value}</div>
            </div>
          ))}
        </section>
      )}

      {/* Devices grid */}
      <section>
        <h3>Devices ({fleet.devices.length})</h3>
        {fleet.devices.length === 0 ? (
          <p style={{ color: '#666' }}>
            No fleet devices connected. Enable iCloud Drive and fleet mode on each device.
          </p>
        ) : (
          <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fill, minmax(250px, 1fr))', gap: '1rem' }}>
            {fleet.devices.map((device) => (
              <div key={device.deviceID} style={{
                border: '1px solid #333',
                borderRadius: '8px',
                padding: '1rem',
                background: '#1a1a2e',
                opacity: device.thermalState === 'critical' ? 0.5 : 1,
              }}>
                <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                  <strong>{device.model}</strong>
                  <span style={{
                    width: '12px', height: '12px', borderRadius: '50%',
                    backgroundColor: thermalColor(device.thermalState),
                    display: 'inline-block',
                  }} />
                </div>
                <div style={{ fontSize: '0.85rem', color: '#999', marginTop: '0.5rem' }}>
                  <div>{device.estimatedTFLOPS.toFixed(1)} TFLOPS</div>
                  <div>{device.availableMemoryGB.toFixed(0)} GB RAM</div>
                  <div>Weight: {(device.computeWeight * 100).toFixed(0)}%</div>
                  <div>Thermal: {device.thermalState}</div>
                  {(device as unknown as Record<string, unknown>).batteryLevel !== undefined && (
                    <div>Battery: {batteryIcon(
                      (device as unknown as Record<string, unknown>).batteryLevel as number,
                      (device as unknown as Record<string, unknown>).isCharging as boolean,
                    )}</div>
                  )}
                </div>
              </div>
            ))}
          </div>
        )}
      </section>

      {/* Job progress bar */}
      {fleet.totalChunks > 0 && (
        <section style={{ marginTop: '1.5rem' }}>
          <h3>Job Progress</h3>
          <div style={{ background: '#333', borderRadius: '4px', height: '24px', overflow: 'hidden', position: 'relative' }}>
            {/* Completed (green) */}
            <div style={{
              background: 'linear-gradient(90deg, #800020, #4B0082)',
              height: '100%',
              width: `${progressPct}%`,
              transition: 'width 0.3s',
              position: 'absolute',
              left: 0,
            }} />
            {/* Failed (red, stacked after completed) */}
            {fleet.failedChunks > 0 && (
              <div style={{
                background: '#ef4444',
                height: '100%',
                width: `${(fleet.failedChunks / fleet.totalChunks) * 100}%`,
                transition: 'width 0.3s',
                position: 'absolute',
                left: `${progressPct}%`,
              }} />
            )}
          </div>
          <p>
            {fleet.completedChunks} / {fleet.totalChunks} chunks completed
            {fleet.failedChunks > 0 && <span style={{ color: '#ef4444' }}> ({fleet.failedChunks} failed)</span>}
            {fleet.orphanedChunks > 0 && <span style={{ color: '#f97316' }}> ({fleet.orphanedChunks} orphaned, retrying)</span>}
          </p>
        </section>
      )}

      {/* Intelligence analysis panel */}
      {fleet.refereeVerdict && (
        <section style={{ marginTop: '1.5rem' }}>
          <IntelligencePanel verdict={fleet.refereeVerdict} title="Fleet Intelligence" />
        </section>
      )}

      {/* Active chunks table */}
      {fleet.activeChunks.length > 0 && (
        <section style={{ marginTop: '1.5rem' }}>
          <h3>Active Chunks</h3>
          <table style={{ width: '100%', borderCollapse: 'collapse' }}>
            <thead>
              <tr style={{ borderBottom: '1px solid #444' }}>
                <th style={{ textAlign: 'left', padding: '0.5rem' }}>Chunk</th>
                <th style={{ textAlign: 'left', padding: '0.5rem' }}>Status</th>
                <th style={{ textAlign: 'left', padding: '0.5rem' }}>Device</th>
                <th style={{ textAlign: 'right', padding: '0.5rem' }}>Chromosomes</th>
                <th style={{ textAlign: 'right', padding: '0.5rem' }}>Retries</th>
                <th style={{ textAlign: 'left', padding: '0.5rem' }}>Priority</th>
              </tr>
            </thead>
            <tbody>
              {fleet.activeChunks.map((chunk) => {
                const chunkAny = chunk as unknown as Record<string, unknown>;
                const status = chunkAny.status as string;
                const statusColor = status === 'completed' ? '#4ade80'
                  : status === 'running' ? '#60a5fa'
                  : status === 'failed' || status === 'permanentlyFailed' ? '#ef4444'
                  : status === 'orphaned' ? '#f97316'
                  : '#999';

                return (
                  <tr key={chunk.id} style={{ borderBottom: '1px solid #222' }}>
                    <td style={{ padding: '0.5rem' }}>{chunk.index + 1}/{chunk.totalChunks}</td>
                    <td style={{ padding: '0.5rem', color: statusColor }}>{status}</td>
                    <td style={{ padding: '0.5rem' }}>{chunk.claimedBy ?? '\u2014'}</td>
                    <td style={{ padding: '0.5rem', textAlign: 'right' }}>{chunk.gaParameters.numChromosomes}</td>
                    <td style={{ padding: '0.5rem', textAlign: 'right' }}>{(chunkAny.retryCount as number) ?? 0}</td>
                    <td style={{ padding: '0.5rem' }}>{(chunkAny.priority as string) ?? 'normal'}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </section>
      )}
    </div>
  );
}
