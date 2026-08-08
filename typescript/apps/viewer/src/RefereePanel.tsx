// RefereePanel.tsx — Thermodynamic referee verdict display component
//
// Renders RefereeVerdict findings as severity-colored cards with
// trust badge, confidence bar, and recommended action banner.
// Integrates into FleetDashboard for fleet-level quality monitoring.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import React from 'react';
import type { RefereeVerdict, RefereeFinding, RefereeSeverity } from '@bonhomme/shared';
import { refereeAllowsPhysicalAffinity } from './claimPresentation.js';

// ─── Severity color mapping ────────────────────────────────────────────────

const SEVERITY_COLORS: Record<RefereeSeverity, { bg: string; border: string; text: string; badge: string }> = {
  pass: { bg: '#e8f5e9', border: '#4caf50', text: '#1b5e20', badge: '#4caf50' },
  advisory: { bg: '#e3f2fd', border: '#2196f3', text: '#0d47a1', badge: '#2196f3' },
  warning: { bg: '#fff8e1', border: '#ffc107', text: '#e65100', badge: '#ff9800' },
  critical: { bg: '#ffebee', border: '#f44336', text: '#b71c1c', badge: '#f44336' },
};

const SEVERITY_ICONS: Record<RefereeSeverity, string> = {
  pass: '\u2713',      // ✓
  advisory: '\u2139',  // ℹ
  warning: '\u26A0',   // ⚠
  critical: '\u2716',  // ✖
};

// ─── Sub-components ─────────────────────────────────────────────────────────

function FindingCard({ finding }: { finding: RefereeFinding }) {
  const colors = SEVERITY_COLORS[finding.severity] ?? SEVERITY_COLORS.advisory;

  return (
    <div
      style={{
        border: `2px solid ${colors.border}`,
        borderRadius: '8px',
        padding: '12px 16px',
        marginBottom: '8px',
        backgroundColor: colors.bg,
      }}
    >
      <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '4px' }}>
        <span
          style={{
            display: 'inline-flex',
            alignItems: 'center',
            justifyContent: 'center',
            width: '24px',
            height: '24px',
            borderRadius: '50%',
            backgroundColor: colors.badge,
            color: '#fff',
            fontSize: '14px',
            fontWeight: 'bold',
          }}
        >
          {SEVERITY_ICONS[finding.severity]}
        </span>
        <span style={{ fontWeight: 600, color: colors.text, fontSize: '14px' }}>
          {finding.title}
        </span>
        <span
          style={{
            marginLeft: 'auto',
            fontSize: '11px',
            padding: '2px 8px',
            borderRadius: '12px',
            backgroundColor: colors.badge,
            color: '#fff',
            textTransform: 'uppercase',
            fontWeight: 600,
          }}
        >
          {finding.severity}
        </span>
      </div>
      <p style={{ margin: 0, fontSize: '13px', color: '#333', lineHeight: 1.4 }}>
        {finding.detail}
      </p>
      <span
        style={{
          display: 'inline-block',
          marginTop: '4px',
          fontSize: '11px',
          color: '#666',
          fontStyle: 'italic',
        }}
      >
        {finding.category}
      </span>
    </div>
  );
}

function ConfidenceBar({ confidence }: { confidence: number }) {
  const pct = Math.round(confidence * 100);
  const color = pct >= 80 ? '#4caf50' : pct >= 50 ? '#ff9800' : '#f44336';

  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
      <span style={{ fontSize: '12px', color: '#666', minWidth: '80px' }}>Confidence</span>
      <div
        style={{
          flex: 1,
          height: '8px',
          backgroundColor: '#e0e0e0',
          borderRadius: '4px',
          overflow: 'hidden',
        }}
      >
        <div
          style={{
            width: `${pct}%`,
            height: '100%',
            backgroundColor: color,
            borderRadius: '4px',
            transition: 'width 0.3s ease',
          }}
        />
      </div>
      <span style={{ fontSize: '12px', fontWeight: 600, color, minWidth: '36px' }}>
        {pct}%
      </span>
    </div>
  );
}

function TrustBadge({ trustworthy }: { trustworthy: boolean }) {
  return (
    <span
      style={{
        display: 'inline-flex',
        alignItems: 'center',
        gap: '6px',
        padding: '4px 12px',
        borderRadius: '16px',
        fontSize: '13px',
        fontWeight: 600,
        backgroundColor: trustworthy ? '#e8f5e9' : '#ffebee',
        color: trustworthy ? '#1b5e20' : '#b71c1c',
        border: `1px solid ${trustworthy ? '#4caf50' : '#f44336'}`,
      }}
    >
      {trustworthy ? '\u2713 Trustworthy' : '\u2716 Not Trustworthy'}
    </span>
  );
}

// ─── Main Component ─────────────────────────────────────────────────────────

export interface RefereePanelProps {
  /** The referee verdict to display */
  verdict: RefereeVerdict | null;
  /** Optional title override */
  title?: string;
  /** Whether to show in compact mode (fewer details) */
  compact?: boolean;
}

export function RefereePanel({ verdict, title = 'Thermodynamic Referee', compact = false }: RefereePanelProps) {
  if (!verdict) {
    return (
      <div
        style={{
          padding: '16px',
          borderRadius: '8px',
          backgroundColor: '#f5f5f5',
          color: '#999',
          textAlign: 'center',
          fontSize: '13px',
        }}
      >
        No referee verdict available. Run docking analysis to generate a verdict.
      </div>
    );
  }

  const bindingPhysical = refereeAllowsPhysicalAffinity(verdict);
  const displayedVerdict: RefereeVerdict = bindingPhysical
    ? verdict
    : {
        findings: [{
          title: 'Physical affinity unavailable',
          detail: 'This result is an ensemble/CF diagnostic only. Binding-physical provenance is absent or unavailable, so physical affinity and energy-unit claims are suppressed.',
          severity: 'warning',
          category: 'affinity',
        }],
        overallTrustworthy: false,
        recommendedAction: 'Treat values as ensemble/CF diagnostics and supply valid SHA-256 artifact identities for calibrated energy, ensemble measure, and the matched association cycle.',
        confidence: Math.min(verdict.confidence, 0.3),
        claimSource: verdict.claimSource,
      };

  return (
    <div
      style={{
        border: '1px solid #e0e0e0',
        borderRadius: '12px',
        padding: '16px',
        backgroundColor: '#fafafa',
      }}
    >
      {/* Header */}
      <div
        style={{
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'space-between',
          marginBottom: '12px',
        }}
      >
        <h3 style={{ margin: 0, fontSize: '16px', fontWeight: 600 }}>{title}</h3>
        <TrustBadge trustworthy={displayedVerdict.overallTrustworthy} />
      </div>

      {/* Confidence bar */}
      <div style={{ marginBottom: '12px' }}>
        <ConfidenceBar confidence={displayedVerdict.confidence} />
      </div>

      {/* Findings */}
      <div style={{ marginBottom: '12px' }}>
        {displayedVerdict.findings.map((finding, i) => (
          <FindingCard key={i} finding={finding} />
        ))}
      </div>

      {/* Recommended Action Banner */}
      {!compact && (
        <div
          style={{
            padding: '10px 14px',
            borderRadius: '8px',
            backgroundColor: displayedVerdict.overallTrustworthy ? '#e8f5e9' : '#fff3e0',
            border: `1px solid ${displayedVerdict.overallTrustworthy ? '#c8e6c9' : '#ffe0b2'}`,
            fontSize: '13px',
            lineHeight: 1.4,
          }}
        >
          <strong>Recommended Action:</strong> {displayedVerdict.recommendedAction}
        </div>
      )}
    </div>
  );
}
