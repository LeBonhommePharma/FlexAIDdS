// HealthEntropyInsight.swift — Health-entropy correlation narratives
//
// Explains correlations between researcher health metrics (HRV, sleep)
// and docking analysis patterns using Apple FoundationModels.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Health-Entropy Context

/// Health-entropy snapshot for a single analysis session.
public struct HealthEntropySnapshot: Sendable, Codable {
    /// Evidence governing the docking-side values below. Health metrics are
    /// instrument readings and are unaffected; `shannonS` and `freeEnergy`
    /// are only physical when this record authorizes a canonical claim.
    public let scientificProvenance: ScientificProvenance?

    /// Strongest claim supported by the declared provenance.
    public var claimValidity: ClaimValidity {
        scientificProvenance?.claimValidity ?? .proxyOnly
    }

    /// Unit label for `freeEnergy`.
    public var energyUnitLabel: String {
        claimValidity == .proxyOnly ? "input-scale units" : "kcal/mol"
    }

    /// Unit label for `shannonS`.
    public var entropyUnitLabel: String {
        claimValidity == .proxyOnly ? "input-scale units/K" : "kcal/mol/K"
    }

    public let hrvSDNN: Double?        // ms
    public let restingHR: Double?      // bpm
    public let sleepHours: Double?
    /// Entropy-like diagnostic in the declared energy domain per kelvin.
    public let shannonS: Double
    public let isConverged: Bool
    public let modeCount: Int
    /// Free-energy-like diagnostic in the declared energy domain.
    public let freeEnergy: Double
    public let timestamp: Date

    public init(hrvSDNN: Double?, restingHR: Double?, sleepHours: Double?,
                shannonS: Double, isConverged: Bool, modeCount: Int,
                freeEnergy: Double, timestamp: Date = Date(),
                scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.hrvSDNN = hrvSDNN
        self.restingHR = restingHR
        self.sleepHours = sleepHours
        self.shannonS = shannonS
        self.isConverged = isConverged
        self.modeCount = modeCount
        self.freeEnergy = freeEnergy
        self.timestamp = timestamp
    }
}

/// Health-entropy context with trend information.
public struct HealthEntropyContext: Sendable, Codable {
    /// Current snapshot
    public let current: HealthEntropySnapshot
    /// Recent history (last 5 snapshots, oldest first)
    public let history: [HealthEntropySnapshot]
    /// HRV trend: "improving", "declining", "stable", "unknown"
    public let hrvTrend: String
    /// Convergence rate trend: "improving", "declining", "stable", "unknown"
    public let convergenceTrend: String

    public init(current: HealthEntropySnapshot, history: [HealthEntropySnapshot],
                hrvTrend: String, convergenceTrend: String) {
        self.current = current
        self.history = history
        self.hrvTrend = hrvTrend
        self.convergenceTrend = convergenceTrend
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Health-entropy insight from the on-device model.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct HealthEntropyInsight: Sendable, Codable {
    /// Summary of health-entropy correlations
    public var correlationSummary: String
    /// Wellness recommendation based on trends
    public var wellnessRecommendation: String
    /// Note about data quality relative to health state
    public var dataQualityNote: String
    /// Confidence in the insight (0.0-1.0)
    public var confidence: Double
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor HealthEntropyInsightActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a wellness-science advisor correlating molecular docking quality \
        with researcher health metrics. All values are pre-computed. \
        Explain how HRV, sleep, and resting heart rate relate to docking \
        result quality (convergence, entropy). Provide practical wellness \
        recommendations. Be empathetic and evidence-based.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Generate health-entropy insight.
    public func analyze(context: HealthEntropyContext) async throws -> HealthEntropyInsight {
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            // Truncate history to last 2 snapshots
            let truncated = HealthEntropyContext(
                current: context.current,
                history: Array(context.history.suffix(2)),
                hrvTrend: context.hrvTrend,
                convergenceTrend: context.convergenceTrend
            )
            prompt = buildPrompt(context: truncated)
        }
        return try await session.respond(to: prompt, generating: HealthEntropyInsight.self).content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(context: HealthEntropyContext) -> String {
        let c = context.current
        var p = "Analyze this health-entropy correlation. Produce a HealthEntropyInsight.\n"
        p += "Claim validity: \(c.claimValidity.rawValue) (docking values are diagnostics, not affinities)\n"
        p += "Current analysis: F=\(String(format: "%.2f", c.freeEnergy)) \(c.energyUnitLabel), "
        p += "S=\(String(format: "%.4f", c.shannonS)) \(c.entropyUnitLabel), \(c.modeCount) modes, "
        p += "converged=\(c.isConverged)\n"

        var metrics: [String] = []
        if let hrv = c.hrvSDNN { metrics.append("HRV: \(String(format: "%.0f", hrv))ms") }
        if let hr = c.restingHR { metrics.append("HR: \(String(format: "%.0f", hr))bpm") }
        if let sleep = c.sleepHours { metrics.append("Sleep: \(String(format: "%.1f", sleep))h") }
        if !metrics.isEmpty { p += "\n" + metrics.joined(separator: ", ") }

        p += "\nHRV trend: \(context.hrvTrend), Convergence trend: \(context.convergenceTrend)"

        if context.history.count >= 2 {
            p += "\nRecent history (\(context.history.count) sessions):"
            for snap in context.history.suffix(3) {
                p += "\n  F=\(String(format: "%.1f", snap.freeEnergy)), "
                p += "converged=\(snap.isConverged)"
                if let h = snap.hrvSDNN { p += ", HRV=\(String(format: "%.0f", h))" }
                if let s = snap.sleepHours { p += ", sleep=\(String(format: "%.1f", s))h" }
            }
        }

        if let hrv = c.hrvSDNN, hrv < 30 { p += "\nFLAG: Very low HRV — acute stress" }
        if let sleep = c.sleepHours, sleep < 5 { p += "\nFLAG: Poor sleep — cognitive load elevated" }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent health-entropy insight.
public struct CrossPlatformHealthEntropyInsight: Sendable, Codable {
    public let correlationSummary: String
    public let wellnessRecommendation: String
    public let dataQualityNote: String
    public let confidence: Double

    public init(correlationSummary: String, wellnessRecommendation: String,
                dataQualityNote: String, confidence: Double) {
        self.correlationSummary = correlationSummary
        self.wellnessRecommendation = wellnessRecommendation
        self.dataQualityNote = dataQualityNote
        self.confidence = confidence
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic health-entropy insight for non-Apple platforms.
public struct RuleBasedHealthEntropyInsight: Sendable {

    public init() {}

    /// Generate health-entropy insight using threshold logic.
    public func analyze(context: HealthEntropyContext) -> CrossPlatformHealthEntropyInsight {
        let c = context.current

        // Correlation summary
        var summary = ""
        if let hrv = c.hrvSDNN {
            if hrv > 60 && c.isConverged {
                summary = "Good HRV (\(String(format: "%.0f", hrv))ms) correlates with converged docking results."
            } else if hrv < 30 {
                summary = "Very low HRV (\(String(format: "%.0f", hrv))ms) during analysis — stress may affect result interpretation."
            } else {
                summary = "HRV at \(String(format: "%.0f", hrv))ms — physiological state is moderate."
            }
        } else {
            summary = "No HRV data available. Connect HealthKit for health-entropy correlation."
        }

        if context.hrvTrend == "improving" && context.convergenceTrend == "improving" {
            summary += " Both health and convergence are improving — positive trend."
        }

        // Wellness recommendation
        var recommendation: String
        if let sleep = c.sleepHours, sleep < 6 {
            recommendation = "Sleep was \(String(format: "%.1f", sleep)) hours. Analyses after 7+ hours sleep tend to produce better-converged results."
        } else if let hrv = c.hrvSDNN, hrv < 40 {
            recommendation = "Consider recovery before critical analysis decisions. Low HRV suggests autonomic stress."
        } else {
            recommendation = "Health metrics are within normal range. Good conditions for analytical work."
        }

        // Data quality note
        var quality: String
        if !c.isConverged {
            quality = "Results not converged — data quality is the limiting factor, not health state."
        } else if let hrv = c.hrvSDNN, hrv < 30 {
            quality = "Results are valid but review with fresh eyes — acute stress may affect subjective interpretation."
        } else {
            quality = "Results are reliable and health state supports good analytical judgment."
        }
        if c.claimValidity == .proxyOnly {
            // The docking side of this correlation is an optimizer diagnostic.
            quality += " Docking values here are proxy diagnostics in \(c.energyUnitLabel), not physical free energies or affinities."
        }

        let confidence: Double
        if c.hrvSDNN != nil && c.sleepHours != nil && context.history.count >= 3 {
            confidence = 0.8
        } else if c.hrvSDNN != nil || c.sleepHours != nil {
            confidence = 0.5
        } else {
            confidence = 0.3
        }

        return CrossPlatformHealthEntropyInsight(
            correlationSummary: summary,
            wellnessRecommendation: recommendation,
            dataQualityNote: quality,
            confidence: confidence
        )
    }
}
