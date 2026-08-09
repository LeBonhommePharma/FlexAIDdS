// CampaignJournalist.swift — Multi-run docking campaign progress summarizer
//
// Summarizes a series of docking campaigns into narrative progress reports
// for lab notebooks and team communication. Reads from AnalysisHistory
// and produces structured summaries using Apple FoundationModels.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Campaign Context

/// Summary of a single run within a campaign.
public struct RunSnapshot: Sendable, Codable {
    /// Evidence governing physical interpretation of this run.
    public let scientificProvenance: ScientificProvenance?
    public let runIndex: Int
    /// Free-energy-like value in the declared energy domain.
    public let freeEnergy: Double
    /// Entropy-like value in the declared energy domain per kelvin.
    public let entropy: Double
    public let isConverged: Bool
    public let modeCount: Int
    public let confidence: String          // "high", "moderate", "low"
    public let timestamp: Date

    public init(runIndex: Int, freeEnergy: Double, entropy: Double,
                isConverged: Bool, modeCount: Int, confidence: String,
                timestamp: Date, scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.runIndex = runIndex
        self.freeEnergy = freeEnergy
        self.entropy = entropy
        self.isConverged = isConverged
        self.modeCount = modeCount
        self.confidence = confidence
        self.timestamp = timestamp
    }
}

/// Full campaign context for summarization.
public struct CampaignContext: Sendable, Codable {
    /// Campaign identifier (e.g., "5HT2A-psilocin")
    public let campaignKey: String
    /// All run snapshots, oldest first (capped to last 8)
    public let runs: [RunSnapshot]
    /// Best free energy across all runs
    public let bestFreeEnergy: Double
    /// Index of the best run
    public let bestRunIndex: Int
    /// Delta-F from first to last run
    public let deltaFFirstToLast: Double
    /// Whether recent runs are within noise (|ΔF| < 0.5 in the declared
    /// energy domain over the last 3 runs)
    public let isStagnating: Bool
    /// Number of converged runs out of total
    public let convergedRunCount: Int

    /// Strongest domain shared by every run. Mixed or missing provenance fails closed.
    public var claimValidity: ClaimValidity {
        guard !runs.isEmpty else { return .proxyOnly }
        if runs.allSatisfy({ $0.scientificProvenance?.claimValidity == .bindingPhysical }) {
            return .bindingPhysical
        }
        if runs.allSatisfy({
            let validity = $0.scientificProvenance?.claimValidity ?? .proxyOnly
            return validity == .canonicalPhysical || validity == .bindingPhysical
        }) {
            return .canonicalPhysical
        }
        return .proxyOnly
    }

    public init(campaignKey: String, runs: [RunSnapshot]) {
        self.campaignKey = campaignKey
        self.runs = Array(runs.suffix(8))

        let best = self.runs.min(by: { $0.freeEnergy < $1.freeEnergy })
        self.bestFreeEnergy = best?.freeEnergy ?? 0
        self.bestRunIndex = best?.runIndex ?? 0

        if let first = self.runs.first, let last = self.runs.last {
            self.deltaFFirstToLast = last.freeEnergy - first.freeEnergy
        } else {
            self.deltaFFirstToLast = 0
        }

        // Stagnation: last 3 runs within 0.5 kcal/mol of each other
        let recent = self.runs.suffix(3)
        if recent.count >= 3 {
            let energies = recent.map(\.freeEnergy)
            let range = (energies.max() ?? 0) - (energies.min() ?? 0)
            self.isStagnating = range < 0.5
        } else {
            self.isStagnating = false
        }

        self.convergedRunCount = self.runs.filter(\.isConverged).count
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Campaign progress summary.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct CampaignSummary: Sendable, Codable {
    /// Campaign identifier
    public var campaignKey: String
    /// Number of runs summarized
    public var runCount: Int
    /// 2-3 sentence progress narrative
    public var progressNarrative: String
    /// Best result description
    public var bestResult: String
    /// Overall trend
    public var trend: String
    /// What to do next
    public var nextStepRecommendation: String
    /// Whether results are publication-quality
    public var readyForPublication: Bool
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor CampaignJournalistActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a scientific writing assistant summarizing molecular docking \
        campaign progress. All thermodynamic values are pre-computed. \
        Write concise, lab-notebook-quality summaries. Report binding-affinity \
        and molecular-optimization conclusions only for binding_physical inputs. \
        Be quantitative and cite specific F values and run numbers. \
        Recommend next steps based on the trajectory.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Summarize a docking campaign.
    public func summarize(context: CampaignContext) async throws -> CampaignSummary {
        if context.claimValidity != .bindingPhysical {
            let canonical = context.claimValidity == .canonicalPhysical
            return CampaignSummary(
                campaignKey: context.campaignKey,
                runCount: context.runs.count,
                progressNarrative: canonical
                    ? "Canonical-ensemble campaign diagnostic; binding affinity is unavailable without a matched association cycle."
                    : "Proxy-only optimizer campaign diagnostic; physical binding affinity is unavailable.",
                bestResult: canonical
                    ? "Lowest canonical F recorded; not a binding-affinity ranking."
                    : "Lowest input-scale proxy score recorded; physical units unavailable.",
                trend: "diagnostic only",
                nextStepRecommendation: "Use convergence/sampling diagnostics only; molecular and lead optimization are unavailable.",
                readyForPublication: false
            )
        }
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            // Truncate to last 4 runs
            let truncated = CampaignContext(
                campaignKey: context.campaignKey,
                runs: Array(context.runs.suffix(4))
            )
            prompt = buildPrompt(context: truncated)
        }
        return try await session.respond(to: prompt, generating: CampaignSummary.self).content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(context: CampaignContext) -> String {
        var p = "Summarize this docking campaign. Produce a CampaignSummary.\n"
        p += "Campaign: \(context.campaignKey), \(context.runs.count) runs\n"
        p += "Best: F=\(String(format: "%.2f", context.bestFreeEnergy))kcal/mol (run \(context.bestRunIndex + 1))\n"
        p += "ΔF first→last: \(String(format: "%.2f", context.deltaFFirstToLast))kcal/mol\n"
        p += "Converged: \(context.convergedRunCount)/\(context.runs.count) runs\n"
        if context.isStagnating { p += "FLAG: Last 3 runs stagnating (<0.5 kcal/mol variation)\n" }

        p += "\nRun history:"
        for run in context.runs {
            p += "\n  Run \(run.runIndex + 1): F=\(String(format: "%.2f", run.freeEnergy))kcal/mol, "
            p += "S=\(String(format: "%.4f", run.entropy)), \(run.modeCount) modes, "
            p += "converged=\(run.isConverged), conf=\(run.confidence)"
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent campaign summary.
public struct CrossPlatformCampaignSummary: Sendable, Codable {
    public let campaignKey: String
    public let runCount: Int
    public let progressNarrative: String
    public let bestResult: String
    public let trend: String
    public let nextStepRecommendation: String
    public let readyForPublication: Bool

    public init(campaignKey: String, runCount: Int, progressNarrative: String,
                bestResult: String, trend: String, nextStepRecommendation: String,
                readyForPublication: Bool) {
        self.campaignKey = campaignKey
        self.runCount = runCount
        self.progressNarrative = progressNarrative
        self.bestResult = bestResult
        self.trend = trend
        self.nextStepRecommendation = nextStepRecommendation
        self.readyForPublication = readyForPublication
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic campaign journalist for non-Apple platforms.
public struct RuleBasedCampaignJournalist: Sendable {

    public init() {}

    /// Build a CampaignContext from AnalysisHistory entries.
    /// Snapshots with unparseable free energy are skipped.
    ///
    /// `OracleAnalysis` is a rendered narrative record: its `inputSummary`
    /// text is not evidence, so provenance can never be recovered by parsing
    /// it. The caller must hand in the run's actual `ScientificProvenance`;
    /// omitting it leaves every snapshot proxy-only, which is the fail-closed
    /// outcome rather than a silent upgrade.
    public func buildContext(
        campaignKey: String,
        analyses: [OracleAnalysis],
        scientificProvenance: ScientificProvenance? = nil
    ) -> CampaignContext {
        let snapshots: [RunSnapshot] = analyses.enumerated().compactMap { i, analysis in
            // Parse F and S from inputSummary (format: "T=300K, F=-10.20 kcal/mol, ...")
            guard let freeEnergy = parseValue(from: analysis.inputSummary, key: "F=") else {
                // Skip snapshots with unparseable free energy rather than defaulting to 0
                return nil
            }
            let entropy = parseValue(from: analysis.inputSummary, key: "S_conf=") ?? parseValue(from: analysis.inputSummary, key: "S=") ?? 0
            let isConverged = analysis.inputSummary.contains("[converged]")
            let modeCount = analysis.structuredBullets.filter { $0.category == .entropy }.count > 0 ? 3 : 1

            return RunSnapshot(
                runIndex: i,
                freeEnergy: freeEnergy,
                entropy: entropy,
                isConverged: isConverged,
                modeCount: modeCount,
                confidence: analysis.overallConfidence.rawValue,
                timestamp: analysis.timestamp,
                scientificProvenance: scientificProvenance
            )
        }

        return CampaignContext(campaignKey: campaignKey, runs: snapshots)
    }

    /// Summarize a campaign using threshold logic.
    public func summarize(context: CampaignContext) -> CrossPlatformCampaignSummary {
        guard !context.runs.isEmpty else {
            return CrossPlatformCampaignSummary(
                campaignKey: context.campaignKey,
                runCount: 0,
                progressNarrative: "No runs available for campaign '\(context.campaignKey)'.",
                bestResult: "No data",
                trend: "insufficient data",
                nextStepRecommendation: "Run at least one docking campaign to begin analysis.",
                readyForPublication: false
            )
        }
        let runs = context.runs
        let bindingPhysical = context.claimValidity == .bindingPhysical
        let canonicalPhysical = context.claimValidity == .canonicalPhysical

        // Trend detection
        let trend: String
        if runs.count < 2 {
            trend = "insufficient data"
        } else if context.deltaFFirstToLast < -2.0 {
            trend = "improving"
        } else if context.deltaFFirstToLast > 1.0 {
            trend = "regressing"
        } else if context.isStagnating {
            trend = "stagnating"
        } else {
            trend = "stable"
        }

        // Progress narrative
        var narrative = "Campaign '\(context.campaignKey)' across \(runs.count) runs: "
        if trend == "improving" {
            if bindingPhysical {
                narrative += "binding affinity improved from F = \(String(format: "%.1f", runs.first?.freeEnergy ?? 0)) to \(String(format: "%.1f", runs.last?.freeEnergy ?? 0)) kcal/mol (ΔF = \(String(format: "%.1f", context.deltaFFirstToLast))). "
            } else if canonicalPhysical {
                narrative += "canonical F decreased from \(String(format: "%.1f", runs.first?.freeEnergy ?? 0)) to \(String(format: "%.1f", runs.last?.freeEnergy ?? 0)) kcal/mol; binding affinity remains unavailable. "
            } else {
                narrative += "optimizer proxy score decreased from \(String(format: "%.1f", runs.first?.freeEnergy ?? 0)) to \(String(format: "%.1f", runs.last?.freeEnergy ?? 0)) input-scale units; binding affinity remains unavailable. "
            }
        } else if trend == "stagnating" {
            // The 0.5 threshold is applied in whatever domain the scores live
            // in; only a canonical/binding record may call that domain kcal/mol.
            let units = (bindingPhysical || canonicalPhysical) ? "kcal/mol" : "input-scale"
            narrative += "last 3 runs are within noise (<0.5 \(units) variation), suggesting the sampled conformational space is exhausted. "
        } else if trend == "regressing" {
            if bindingPhysical {
                narrative += "binding affinity has worsened from F = \(String(format: "%.1f", runs.first?.freeEnergy ?? 0)) to \(String(format: "%.1f", runs.last?.freeEnergy ?? 0)) kcal/mol — review parameter changes. "
            } else {
                narrative += "the diagnostic score increased across runs; this is not a binding-affinity trend. "
            }
        } else {
            narrative += "results are stable across runs. "
        }
        narrative += "Convergence achieved in \(context.convergedRunCount)/\(runs.count) runs."

        // Best result
        let bestResult = (bindingPhysical || canonicalPhysical
            ? "Run \(context.bestRunIndex + 1): F = \(String(format: "%.2f", context.bestFreeEnergy)) kcal/mol"
            : "Run \(context.bestRunIndex + 1): proxy score = \(String(format: "%.2f", context.bestFreeEnergy)) input-scale")
            + (runs[safe: context.bestRunIndex]?.isConverged == true ? " (converged)" : " (not converged)")

        // Next step
        let recommendation: String
        if !bindingPhysical {
            recommendation = context.convergedRunCount < runs.count
                ? "Improve convergence and sampling. Binding-affinity and molecular/lead-optimization guidance are unavailable."
                : "Diagnostic campaign complete. Supply a calibrated matched association cycle before binding-affinity or molecular-optimization guidance."
        } else if context.isStagnating && context.convergedRunCount == runs.count {
            recommendation = "Conformational space exhausted. Consider modifying the ligand scaffold or targeting alternative binding sites."
        } else if context.convergedRunCount < runs.count / 2 {
            recommendation = "Most runs did not converge. Increase GA generations and population size."
        } else if trend == "improving" {
            recommendation = "Continue optimization — trajectory shows consistent improvement."
        } else if trend == "regressing" {
            recommendation = "Revert recent parameter changes and restart from the best-performing configuration."
        } else {
            recommendation = "Results are stable. Proceed with lead optimization of the best binding mode."
        }

        // Publication readiness
        let ready = bindingPhysical
            && context.convergedRunCount >= max(runs.count - 1, 1)
            && context.bestFreeEnergy < -5.0
            && (context.isStagnating || runs.count >= 5)

        return CrossPlatformCampaignSummary(
            campaignKey: context.campaignKey,
            runCount: runs.count,
            progressNarrative: narrative,
            bestResult: bestResult,
            trend: trend,
            nextStepRecommendation: recommendation,
            readyForPublication: ready
        )
    }

    // MARK: - Helpers

    private func parseValue(from summary: String, key: String) -> Double? {
        guard let range = summary.range(of: key) else { return nil }
        let after = summary[range.upperBound...]
        // Extract the number (possibly negative, with decimals)
        var numStr = ""
        for ch in after {
            if ch == "-" || ch == "." || ch.isNumber {
                numStr.append(ch)
            } else if !numStr.isEmpty {
                break
            }
        }
        return Double(numStr)
    }
}

// MARK: - Array Safe Index

private extension Array {
    subscript(safe index: Index) -> Element? {
        indices.contains(index) ? self[index] : nil
    }
}
