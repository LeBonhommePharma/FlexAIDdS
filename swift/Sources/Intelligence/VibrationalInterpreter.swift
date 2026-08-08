// VibrationalInterpreter.swift — Vibrational mode interpretation for binding
//
// Explains which protein motions matter for binding and why, using
// Apple FoundationModels. Pre-computes all normal mode data on CPU;
// the LLM interprets structural implications.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Vibrational Context

/// Summary of a single normal mode for LLM interpretation.
public struct NormalModeSummary: Sendable, Codable {
    /// Mode index (0-based)
    public let index: Int
    /// Frequency in cm⁻¹
    public let frequencyCm: Double
    /// Entropy contribution (kcal/mol/K)
    public let entropyContribution: Double
    /// Fraction of total vibrational entropy
    public let entropyFraction: Double
    /// Whether this mode shifts significantly upon ligand binding
    public let shiftsOnBinding: Bool
    /// Direction of shift: "restricted" (lower entropy) or "enhanced" (higher)
    public let shiftDirection: String?

    public init(index: Int, frequencyCm: Double, entropyContribution: Double,
                entropyFraction: Double, shiftsOnBinding: Bool, shiftDirection: String?) {
        self.index = index
        self.frequencyCm = frequencyCm
        self.entropyContribution = entropyContribution
        self.entropyFraction = entropyFraction
        self.shiftsOnBinding = shiftsOnBinding
        self.shiftDirection = shiftDirection
    }
}

/// Full vibrational analysis context.
public struct VibrationalContext: Sendable, Codable {
    /// Evidence governing physical binding and design interpretation.
    public let scientificProvenance: ScientificProvenance?
    /// Total vibrational entropy S_vib (kcal/mol/K)
    public let totalSVib: Double
    /// Configurational entropy S_conf (nats)
    public let totalSConf: Double
    /// Vibrational dominance ratio (S_vib / S_conf_physical)
    public let vibrationalDominance: Double
    /// Temperature (K)
    public let temperature: Double
    /// Top 5 modes by entropy contribution
    public let topModes: [NormalModeSummary]
    /// Total number of non-trivial modes
    public let totalModeCount: Int
    /// Whether binding restricts dominant motions (net S_vib decrease)
    public let bindingRestrictsMotion: Bool

    public init(totalSVib: Double, totalSConf: Double, vibrationalDominance: Double,
                temperature: Double, topModes: [NormalModeSummary], totalModeCount: Int,
                bindingRestrictsMotion: Bool,
                scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.totalSVib = totalSVib
        self.totalSConf = totalSConf
        self.vibrationalDominance = vibrationalDominance
        self.temperature = temperature
        self.topModes = topModes
        self.totalModeCount = totalModeCount
        self.bindingRestrictsMotion = bindingRestrictsMotion
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Vibrational mode insight.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct VibrationalInsight: Sendable, Codable {
    /// Description of the dominant protein motions
    public var dominantMotionDescription: String
    /// How ligand binding affects these motions
    public var bindingImpact: String
    /// Implications for drug design
    public var designImplication: String
    /// Whether binding is entropy-driven from the vibrational perspective
    public var isEntropicallyDriven: Bool
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor VibrationalInterpreterActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a structural biophysics advisor interpreting protein normal mode \
        vibrations and their impact on ligand binding. All frequencies and entropy \
        values are pre-computed — DO NOT calculate. Explain which motions matter \
        for binding, whether ligand binding restricts important dynamics, and \
        what this means for drug design only for binding_physical provenance. \
        Otherwise return diagnostic-only language and state that binding/design claims are unavailable.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Interpret vibrational mode data.
    public func interpret(context: VibrationalContext) async throws -> VibrationalInsight {
        let validity = context.scientificProvenance?.claimValidity ?? .proxyOnly
        if validity != .bindingPhysical {
            return VibrationalInsight(
                dominantMotionDescription: validity == .canonicalPhysical
                    ? "Canonical vibrational-mode diagnostic; no matched association cycle is available."
                    : "Proxy-only vibrational-mode diagnostic; physical units and binding interpretation are unavailable.",
                bindingImpact: "Binding impact is unavailable without binding_physical provenance.",
                designImplication: "Molecular and lead-optimization guidance is unavailable.",
                isEntropicallyDriven: false
            )
        }
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            let truncated = VibrationalContext(
                totalSVib: context.totalSVib, totalSConf: context.totalSConf,
                vibrationalDominance: context.vibrationalDominance,
                temperature: context.temperature,
                topModes: Array(context.topModes.prefix(3)),
                totalModeCount: context.totalModeCount,
                bindingRestrictsMotion: context.bindingRestrictsMotion,
                scientificProvenance: context.scientificProvenance
            )
            prompt = buildPrompt(context: truncated)
        }
        return try await session.respond(to: prompt, generating: VibrationalInsight.self).content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(context: VibrationalContext) -> String {
        var p = "Interpret these vibrational modes. Produce a VibrationalInsight.\n"
        p += "T=\(context.temperature)K, \(context.totalModeCount) modes total\n"
        p += "S_vib=\(String(format: "%.6f", context.totalSVib))kcal/mol/K, "
        p += "S_conf=\(String(format: "%.4f", context.totalSConf))nats\n"
        p += "Vibrational dominance: \(String(format: "%.1f", context.vibrationalDominance))x\n"
        p += "Binding restricts motion: \(context.bindingRestrictsMotion ? "YES" : "NO")\n"

        p += "\nTop modes:"
        for mode in context.topModes {
            p += "\n  Mode \(mode.index): \(String(format: "%.1f", mode.frequencyCm))cm⁻¹, "
            p += "ΔS=\(String(format: "%.6f", mode.entropyContribution))kcal/mol/K "
            p += "(\(String(format: "%.0f", mode.entropyFraction * 100))%)"
            if mode.shiftsOnBinding {
                p += " [shifts: \(mode.shiftDirection ?? "unknown")]"
            }
        }

        if context.vibrationalDominance > 3.0 {
            p += "\nFLAG: Vibrational entropy dominates — backbone flexibility exceeds ligand conformational space"
        }
        if context.bindingRestrictsMotion {
            p += "\nFLAG: Binding restricts dominant motions — entropic cost of binding"
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent vibrational insight.
public struct CrossPlatformVibrationalInsight: Sendable, Codable {
    public let dominantMotionDescription: String
    public let bindingImpact: String
    public let designImplication: String
    public let isEntropicallyDriven: Bool

    public init(dominantMotionDescription: String, bindingImpact: String,
                designImplication: String, isEntropicallyDriven: Bool) {
        self.dominantMotionDescription = dominantMotionDescription
        self.bindingImpact = bindingImpact
        self.designImplication = designImplication
        self.isEntropicallyDriven = isEntropicallyDriven
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic vibrational interpreter for non-Apple platforms.
public struct RuleBasedVibrationalInterpreter: Sendable {

    public init() {}

    /// Interpret vibrational mode data using threshold logic.
    public func interpret(context: VibrationalContext) -> CrossPlatformVibrationalInsight {
        let validity = context.scientificProvenance?.claimValidity ?? .proxyOnly
        if validity != .bindingPhysical {
            let diagnostic: String
            if let top = context.topModes.first {
                diagnostic = validity == .canonicalPhysical
                    ? "Canonical vibrational diagnostic: mode \(top.index) is dominant; binding impact is unavailable without a matched association cycle."
                    : "Proxy-only mode diagnostic: mode \(top.index) is dominant at reported frequency \(String(format: "%.1f", top.frequencyCm)); physical binding interpretation is unavailable."
            } else {
                diagnostic = "No significant diagnostic vibrational modes detected."
            }
            return CrossPlatformVibrationalInsight(
                dominantMotionDescription: diagnostic,
                bindingImpact: "Binding impact is unavailable without binding_physical provenance.",
                designImplication: "Molecular and lead-optimization guidance is unavailable.",
                isEntropicallyDriven: false
            )
        }

        guard !context.topModes.isEmpty else {
            return CrossPlatformVibrationalInsight(
                dominantMotionDescription: "No significant vibrational modes detected.",
                bindingImpact: "Insufficient data to assess binding impact on protein dynamics.",
                designImplication: "Enable ENCoM vibrational analysis for entropy-informed design.",
                isEntropicallyDriven: false
            )
        }
        // Dominant motion description
        let dominantMotion: String
        if let top = context.topModes.first {
            let freqType: String
            if top.frequencyCm < 50 {
                freqType = "very low-frequency collective breathing"
            } else if top.frequencyCm < 200 {
                freqType = "low-frequency domain motion"
            } else if top.frequencyCm < 500 {
                freqType = "moderate-frequency loop flexibility"
            } else {
                freqType = "high-frequency side-chain vibration"
            }
            dominantMotion = "Dominant motion is \(freqType) (mode \(top.index), \(String(format: "%.1f", top.frequencyCm)) cm⁻¹) contributing \(String(format: "%.0f", top.entropyFraction * 100))% of vibrational entropy."
        } else {
            dominantMotion = "No significant vibrational modes detected."
        }

        // Binding impact
        let impact: String
        if context.bindingRestrictsMotion {
            let restrictedModes = context.topModes.filter { $0.shiftsOnBinding && $0.shiftDirection == "restricted" }
            if restrictedModes.isEmpty {
                impact = "Ligand binding restricts protein motion — entropic cost estimated from overall S_vib decrease."
            } else {
                impact = "Ligand binding restricts \(restrictedModes.count) dominant mode(s) — costs \(String(format: "%.4f", context.totalSVib)) kcal/mol/K in vibrational entropy."
            }
        } else {
            impact = "Ligand binding does not significantly restrict protein dynamics. Vibrational entropy preserved."
        }

        // Design implication
        let design: String
        if context.vibrationalDominance > 3.0 {
            design = "Backbone flexibility dominates the entropy landscape. Consider flexible-receptor docking or ensemble-averaged scoring."
        } else if context.bindingRestrictsMotion {
            design = "Binding incurs entropic cost from restricted motion. A smaller fragment may preserve loop flexibility while maintaining key interactions."
        } else {
            design = "Vibrational entropy is balanced with conformational entropy. Standard lead optimization appropriate."
        }

        let isEntropicallyDriven = context.vibrationalDominance > 2.0

        return CrossPlatformVibrationalInsight(
            dominantMotionDescription: dominantMotion,
            bindingImpact: impact,
            designImplication: design,
            isEntropicallyDriven: isEntropicallyDriven
        )
    }
}
