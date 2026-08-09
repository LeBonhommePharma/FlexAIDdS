// IntelligencePanel.tsx — Unified intelligence feature display component
//
// Renders outputs from all Intelligence module features in a single panel.
// Each section only appears when the corresponding data is provided.
// Composes with existing RefereePanel for verdict display.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import React from 'react';
import type { RefereeVerdict } from '@bonhomme/shared';
import type { BindingModeNarrative } from '@bonhomme/shared';
import type { CleftAssessment, DruggabilityTier } from '@bonhomme/shared';
import type { ConvergenceCoaching, GAAdvice } from '@bonhomme/shared';
import type { SelectivityAnalysis, SelectivityDriver } from '@bonhomme/shared';
import type { PoseQualityReport } from '@bonhomme/shared';
import { RefereePanel } from './RefereePanel.js';
import {
  cleftAllowsDruggabilityClaim,
  selectivityAllowsPhysicalAffinity,
} from './claimPresentation.js';

// ─── Props ──────────────────────────────────────────────────────────────────

export interface IntelligencePanelProps {
  /** Thermodynamic referee verdict */
  verdict?: RefereeVerdict | null;
  /** Binding mode narrative */
  modeNarrative?: BindingModeNarrative | null;
  /** Cleft druggability assessment */
  cleftAssessment?: CleftAssessment | null;
  /** GA convergence coaching */
  convergenceCoaching?: ConvergenceCoaching | null;
  /** Multi-target selectivity analysis */
  selectivityAnalysis?: SelectivityAnalysis | null;
  /** Pose quality report */
  poseQuality?: PoseQualityReport | null;
  /** Panel title */
  title?: string;
}

// ─── Color mappings ─────────────────────────────────────────────────────────

const DRUGGABILITY_COLORS: Record<DruggabilityTier, { bg: string; text: string }> = {
  high: { bg: '#e8f5e9', text: '#1b5e20' },
  moderate: { bg: '#e3f2fd', text: '#0d47a1' },
  low: { bg: '#fff8e1', text: '#e65100' },
  undruggable: { bg: '#ffebee', text: '#b71c1c' },
};

const ADVICE_COLORS: Record<GAAdvice, { bg: string; text: string }> = {
  continueRun: { bg: '#e8f5e9', text: '#1b5e20' },
  stopEarly: { bg: '#e3f2fd', text: '#0d47a1' },
  increasePopulation: { bg: '#fff8e1', text: '#e65100' },
  increaseMutationRate: { bg: '#fff8e1', text: '#e65100' },
  restart: { bg: '#ffebee', text: '#b71c1c' },
};

const ADVICE_LABELS: Record<GAAdvice, string> = {
  continueRun: 'Continue Run',
  stopEarly: 'Stop Early',
  increasePopulation: 'Increase Population',
  increaseMutationRate: 'Increase Mutation Rate',
  restart: 'Restart',
};

const DRIVER_COLORS: Record<SelectivityDriver, { bg: string; text: string }> = {
  enthalpic: { bg: '#e3f2fd', text: '#0d47a1' },
  entropic: { bg: '#f3e5f5', text: '#4a148c' },
  mixed: { bg: '#fff8e1', text: '#e65100' },
  inconclusive: { bg: '#f5f5f5', text: '#616161' },
};

// ─── Shared sub-components ──────────────────────────────────────────────────

function SectionHeader({ title }: { title: string }) {
  return (
    <h4 style={{ margin: '0 0 8px', fontSize: '14px', fontWeight: 600, color: '#333' }}>
      {title}
    </h4>
  );
}

function Badge({ label, bg, text }: { label: string; bg: string; text: string }) {
  return (
    <span
      style={{
        display: 'inline-block',
        padding: '3px 10px',
        borderRadius: '12px',
        fontSize: '12px',
        fontWeight: 600,
        backgroundColor: bg,
        color: text,
        textTransform: 'uppercase',
      }}
    >
      {label}
    </span>
  );
}

function ConfidenceBar({ confidence }: { confidence: number }) {
  const pct = Math.round(confidence * 100);
  const color = pct >= 80 ? '#4caf50' : pct >= 50 ? '#ff9800' : '#f44336';
  return (
    <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginTop: '6px' }}>
      <span style={{ fontSize: '11px', color: '#666', minWidth: '70px' }}>Confidence</span>
      <div style={{ flex: 1, height: '6px', backgroundColor: '#e0e0e0', borderRadius: '3px', overflow: 'hidden' }}>
        <div style={{ width: `${pct}%`, height: '100%', backgroundColor: color, borderRadius: '3px' }} />
      </div>
      <span style={{ fontSize: '11px', fontWeight: 600, color, minWidth: '32px' }}>{pct}%</span>
    </div>
  );
}

function DetailText({ text }: { text: string }) {
  return <p style={{ margin: '4px 0', fontSize: '13px', color: '#333', lineHeight: 1.4 }}>{text}</p>;
}

// ─── Section: Binding Modes ─────────────────────────────────────────────────

function BindingModesSection({ narrative }: { narrative: BindingModeNarrative }) {
  return (
    <div style={{ marginBottom: '16px' }}>
      <SectionHeader title="Binding Modes" />
      {narrative.modeDescriptions.map((desc, i) => (
        <div
          key={i}
          style={{
            border: '1px solid #e0e0e0',
            borderRadius: '6px',
            padding: '10px 12px',
            marginBottom: '6px',
            backgroundColor: '#fff',
          }}
        >
          <DetailText text={desc.characterization} />
          <p style={{ margin: '2px 0 0', fontSize: '12px', color: '#666', fontStyle: 'italic' }}>
            {desc.optimizationHint}
          </p>
        </div>
      ))}
      <DetailText text={narrative.selectivityInsight} />
      <ConfidenceBar confidence={narrative.confidence} />
    </div>
  );
}

// ─── Section: Cleft Assessment ──────────────────────────────────────────────

const DIAGNOSTIC_BADGE = { bg: '#f5f5f5', text: '#616161' };

function CleftSection({ assessment }: { assessment: CleftAssessment }) {
  // The tier itself is the druggability claim; without binding-physical
  // evidence it is presented as a shape diagnostic, not a verdict.
  const druggabilityClaimed = cleftAllowsDruggabilityClaim(assessment);
  const colors = druggabilityClaimed
    ? DRUGGABILITY_COLORS[assessment.druggability]
    : DIAGNOSTIC_BADGE;
  return (
    <div style={{ marginBottom: '16px' }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '8px' }}>
        <SectionHeader title={druggabilityClaimed ? 'Pocket Quality' : 'Pocket Geometry'} />
        <Badge
          label={druggabilityClaimed ? assessment.druggability : 'geometry diagnostic'}
          bg={colors.bg}
          text={colors.text}
        />
      </div>
      <DetailText text={assessment.summary} />
      <p style={{ margin: '4px 0', fontSize: '12px', color: '#555' }}>
        <strong>Suggested ligand:</strong> {assessment.suggestedLigandProperties}
      </p>
      {assessment.warnings.length > 0 && (
        <div style={{ marginTop: '6px' }}>
          {assessment.warnings.map((w, i) => (
            <p key={i} style={{ margin: '2px 0', fontSize: '12px', color: '#e65100' }}>
              {'\u26A0'} {w}
            </p>
          ))}
        </div>
      )}
    </div>
  );
}

// ─── Section: Convergence Coaching ──────────────────────────────────────────

function ConvergenceSection({ coaching }: { coaching: ConvergenceCoaching }) {
  const colors = ADVICE_COLORS[coaching.advice];
  return (
    <div style={{ marginBottom: '16px' }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '8px' }}>
        <SectionHeader title="Convergence" />
        <Badge label={ADVICE_LABELS[coaching.advice]} bg={colors.bg} text={colors.text} />
      </div>
      <DetailText text={coaching.reasoning} />
      {coaching.estimatedGenerationsRemaining != null && (
        <p style={{ margin: '4px 0', fontSize: '12px', color: '#555' }}>
          <strong>Estimated remaining:</strong> {coaching.estimatedGenerationsRemaining} generations
        </p>
      )}
      <ConfidenceBar confidence={coaching.confidence} />
    </div>
  );
}

// ─── Section: Selectivity ───────────────────────────────────────────────────

/** Fixed text shown instead of producer prose that we cannot vouch for. */
const PROXY_SELECTIVITY_EXPLANATION =
  'Target ordering is an ensemble/CF diagnostic. Physical affinity, potency '
  + 'and selectivity are unavailable: the supplied records do not carry '
  + 'binding-physical provenance.';
const PROXY_SELECTIVITY_SUGGESTION =
  'Supply calibrated energy, ensemble-measure and matched association-cycle '
  + 'SHA-256 artifact identities for both targets before any affinity-driven '
  + 'design decision.';

function SelectivitySection({ analysis }: { analysis: SelectivityAnalysis }) {
  const bindingPhysical = selectivityAllowsPhysicalAffinity(analysis);
  const displayedDriver = bindingPhysical ? analysis.driver : 'inconclusive';
  const colors = DRIVER_COLORS[displayedDriver];
  // `explanation` / `designSuggestion` are free text from whoever produced the
  // payload. For proxy data they are the most likely place for an unbacked
  // affinity or potency claim to leak into the UI, so they are replaced
  // wholesale rather than annotated.
  const explanation = bindingPhysical
    ? analysis.explanation
    : PROXY_SELECTIVITY_EXPLANATION;
  const designSuggestion = bindingPhysical
    ? analysis.designSuggestion
    : PROXY_SELECTIVITY_SUGGESTION;
  return (
    <div style={{ marginBottom: '16px' }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '8px' }}>
        <SectionHeader title={bindingPhysical ? 'Selectivity' : 'Target Ensemble / CF Diagnostic'} />
        <Badge label={bindingPhysical ? analysis.driver : 'diagnostic'} bg={colors.bg} text={colors.text} />
      </div>
      <p style={{ margin: '4px 0', fontSize: '13px', color: '#333' }}>
        <strong>{bindingPhysical ? 'Preferred:' : 'Top diagnostic rank:'}</strong>{' '}
        {analysis.preferredTarget}{' '}
        {bindingPhysical
          ? `(DDG = ${analysis.deltaG.toFixed(2)} kcal/mol)`
          : `(ensemble/CF diagnostic difference = ${analysis.deltaG.toFixed(2)} source units; physical affinity unavailable)`}
      </p>
      <DetailText text={explanation} />
      <p style={{ margin: '4px 0', fontSize: '12px', color: '#555', fontStyle: 'italic' }}>
        {designSuggestion}
      </p>
    </div>
  );
}

// ─── Section: Pose Quality ──────────────────────────────────────────────────

function PoseQualitySection({ report }: { report: PoseQualityReport }) {
  const consensusColor =
    report.poseConsensus === 'strong' ? '#4caf50' :
    report.poseConsensus === 'moderate' ? '#ff9800' :
    report.poseConsensus === 'ambiguous' ? '#f44336' : '#9e9e9e';
  return (
    <div style={{ marginBottom: '16px' }}>
      <div style={{ display: 'flex', alignItems: 'center', gap: '8px', marginBottom: '8px' }}>
        <SectionHeader title="Pose Quality" />
        <Badge
          label={report.poseConsensus}
          bg={report.poseConsensus === 'strong' ? '#e8f5e9' : report.poseConsensus === 'moderate' ? '#fff8e1' : '#ffebee'}
          text={consensusColor}
        />
      </div>
      <DetailText text={report.topPoseSummary} />
      <DetailText text={report.scoreWeightAlignment} />
      <ConfidenceBar confidence={report.confidenceInTopPose} />
      <p style={{ margin: '6px 0 0', fontSize: '12px', color: '#555', fontStyle: 'italic' }}>
        {report.medicinalChemistryNote}
      </p>
    </div>
  );
}

// ─── Main Component ─────────────────────────────────────────────────────────

export function IntelligencePanel({
  verdict,
  modeNarrative,
  cleftAssessment,
  convergenceCoaching,
  selectivityAnalysis,
  poseQuality,
  title = 'Intelligence Analysis',
}: IntelligencePanelProps) {
  const hasAny = verdict || modeNarrative || cleftAssessment || convergenceCoaching || selectivityAnalysis || poseQuality;

  if (!hasAny) {
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
        No intelligence analysis available. Run docking analysis to generate insights.
      </div>
    );
  }

  return (
    <div
      style={{
        border: '1px solid #e0e0e0',
        borderRadius: '12px',
        padding: '16px',
        backgroundColor: '#fafafa',
      }}
    >
      <h3 style={{ margin: '0 0 16px', fontSize: '16px', fontWeight: 600 }}>{title}</h3>

      {verdict && (
        <div style={{ marginBottom: '16px' }}>
          <RefereePanel verdict={verdict} title="Thermodynamic Referee" compact />
        </div>
      )}

      {modeNarrative && <BindingModesSection narrative={modeNarrative} />}
      {cleftAssessment && <CleftSection assessment={cleftAssessment} />}
      {convergenceCoaching && <ConvergenceSection coaching={convergenceCoaching} />}
      {selectivityAnalysis && <SelectivitySection analysis={selectivityAnalysis} />}
      {poseQuality && <PoseQualitySection report={poseQuality} />}
    </div>
  );
}
