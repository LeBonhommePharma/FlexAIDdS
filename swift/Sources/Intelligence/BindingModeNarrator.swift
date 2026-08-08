// BindingModeNarrator.swift — Natural-language binding mode interpretation
//
// Generates plain-English explanations of what distinguishes each binding mode
// using Apple FoundationModels. Pre-computes all thermodynamic data on CPU;
// the LLM only interprets pre-digested values.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Mode Context

/// Summary of a single binding mode for LLM interpretation.
public struct ModeProfile: Sendable, Codable {
    public let index: Int
    public let poseCount: Int
    public let freeEnergy: Double        // kcal/mol
    public let entropy: Double           // kcal/mol/K
    public let enthalpy: Double          // kcal/mol
    public let heatCapacity: Double
    public let boltzmannWeight: Double   // fraction of partition function
    public let isEntropyDriven: Bool     // |−TS| > |ΔH| contribution
    public let isEnthalpyDriven: Bool    // |ΔH| > |−TS| contribution

    public init(index: Int, poseCount: Int, freeEnergy: Double, entropy: Double,
                enthalpy: Double, heatCapacity: Double, boltzmannWeight: Double,
                temperature: Double) {
        self.index = index
        self.poseCount = poseCount
        self.freeEnergy = freeEnergy
        self.entropy = entropy
        self.enthalpy = enthalpy
        self.heatCapacity = heatCapacity
        self.boltzmannWeight = boltzmannWeight
        let tsContribution = abs(temperature * entropy)
        let hContribution = abs(enthalpy)
        self.isEntropyDriven = tsContribution > hContribution * 1.5
        self.isEnthalpyDriven = hContribution > tsContribution * 1.5
    }
}

/// All mode data pre-computed on CPU before prompting.
public struct PreComputedModeContext: Sendable, Codable {
    /// Evidence governing whether mode statistics support binding/SAR claims.
    public let scientificProvenance: ScientificProvenance?
    public let temperature: Double
    public let totalModes: Int
    public let totalPoses: Int
    public let globalFreeEnergy: Double
    public let modes: [ModeProfile]  // top 3 by free energy
    public let entropyImbalance: Double  // max/min ratio
    public let dominantModeIndex: Int    // mode absorbing most Boltzmann weight

    public var claimValidity: ClaimValidity {
        scientificProvenance?.claimValidity ?? .proxyOnly
    }

    public init(
        temperature: Double,
        totalModes: Int,
        totalPoses: Int,
        globalFreeEnergy: Double,
        modes: [ModeProfile],
        entropyImbalance: Double,
        dominantModeIndex: Int,
        scientificProvenance: ScientificProvenance? = nil
    ) {
        self.scientificProvenance = scientificProvenance
        self.temperature = temperature
        self.totalModes = totalModes
        self.totalPoses = totalPoses
        self.globalFreeEnergy = globalFreeEnergy
        self.modes = modes
        self.entropyImbalance = entropyImbalance
        self.dominantModeIndex = dominantModeIndex
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Description of a single binding mode.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct ModeDescription: Sendable, Codable {
    /// Plain-English characterization (e.g., "tight enthalpy-driven lock")
    public var characterization: String
    /// SAR optimization hint
    public var optimizationHint: String
}

/// Narrative covering all top binding modes.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct BindingModeNarrative: Sendable, Codable {
    /// Per-mode descriptions (indexed same as input modes)
    public var modeDescriptions: [ModeDescription]
    /// Which mode to prioritize for lead optimization and why
    public var selectivityInsight: String
    /// Confidence in the narrative (0.0-1.0)
    public var confidence: Double
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor BindingModeNarrator {
    private let session: LanguageModelSession
    private var latestClaimValidity: ClaimValidity = .proxyOnly

    private static let instructions = """
        You are a medicinal chemistry advisor interpreting molecular docking binding modes. \
        All thermodynamic values are pre-computed — DO NOT perform arithmetic. \
        Explain each binding mode in plain English: what drives its stability \
        (enthalpy vs entropy), how tight/loose it is, and what a chemist should \
        optimize. Binding affinity and molecular/SAR optimization are permitted \
        only for binding_physical provenance. Be concise and actionable.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Generate narratives for the top binding modes.
    public func narrate(context: PreComputedModeContext) async throws -> BindingModeNarrative {
        latestClaimValidity = context.claimValidity
        if context.claimValidity != .bindingPhysical {
            let domain = context.claimValidity == .canonicalPhysical
                ? "canonical-ensemble"
                : "optimizer/input-scale"
            return BindingModeNarrative(
                modeDescriptions: context.modes.map { mode in
                    ModeDescription(
                        characterization: "Mode \(mode.index + 1): \(domain) cluster diagnostic with \(mode.poseCount) poses.",
                        optimizationHint: "Molecular and SAR optimization unavailable without a calibrated matched association cycle."
                    )
                },
                selectivityInsight: "Binding-mode prioritization for lead optimization is unavailable; use these clusters only as \(domain) diagnostics.",
                confidence: 1.0
            )
        }
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            prompt = buildTruncatedPrompt(context: context)
        }
        return try await session.respond(to: prompt, generating: BindingModeNarrative.self).content
    }

    /// Follow-up question about the binding modes.
    public func followUp(_ question: String) async throws -> String {
        if latestClaimValidity != .bindingPhysical {
            return "Binding affinity, selectivity, and molecular/SAR optimization are unavailable without binding_physical provenance."
        }
        let response = try await session.respond(to: question)
        return response.content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildTruncatedPrompt(context: PreComputedModeContext) -> String {
        // Truncate to top 2 modes
        let truncatedContext = PreComputedModeContext(
            temperature: context.temperature, totalModes: context.totalModes,
            totalPoses: context.totalPoses, globalFreeEnergy: context.globalFreeEnergy,
            modes: Array(context.modes.prefix(2)), entropyImbalance: context.entropyImbalance,
            dominantModeIndex: context.dominantModeIndex,
            scientificProvenance: context.scientificProvenance
        )
        return buildPrompt(context: truncatedContext)
    }

    private func buildPrompt(context: PreComputedModeContext) -> String {
        var p = """
            Describe these binding modes. Produce a BindingModeNarrative.

            Global: T=\(context.temperature)K, F=\(String(format: "%.2f", context.globalFreeEnergy))kcal/mol, \
            \(context.totalModes) modes, \(context.totalPoses) poses
            """

        for mode in context.modes {
            p += "\nMode \(mode.index): \(mode.poseCount) poses, "
            p += "F=\(String(format: "%.2f", mode.freeEnergy))kcal/mol, "
            p += "S=\(String(format: "%.4f", mode.entropy))kcal/mol/K, "
            p += "H=\(String(format: "%.2f", mode.enthalpy))kcal/mol, "
            p += "w=\(String(format: "%.1f", mode.boltzmannWeight * 100))%"
            if mode.isEntropyDriven { p += " [entropy-driven]" }
            if mode.isEnthalpyDriven { p += " [enthalpy-driven]" }
        }

        if context.entropyImbalance > 5.0 {
            p += "\nFLAG: Entropy imbalance \(String(format: "%.1f", context.entropyImbalance))x — one mode dominates diversity"
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent mode description.
public struct CrossPlatformModeDescription: Sendable, Codable {
    public let characterization: String
    public let optimizationHint: String

    public init(characterization: String, optimizationHint: String) {
        self.characterization = characterization
        self.optimizationHint = optimizationHint
    }
}

/// Platform-independent binding mode narrative.
public struct CrossPlatformModeNarrative: Sendable, Codable {
    public let modeDescriptions: [CrossPlatformModeDescription]
    public let selectivityInsight: String
    public let confidence: Double

    public init(modeDescriptions: [CrossPlatformModeDescription], selectivityInsight: String, confidence: Double) {
        self.modeDescriptions = modeDescriptions
        self.selectivityInsight = selectivityInsight
        self.confidence = confidence
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic mode narrator for non-Apple platforms.
public struct RuleBasedModeNarrator: Sendable {

    public init() {}

    /// Build a mode context from a DockingResult.
    public func buildContext(from result: DockingResult) -> PreComputedModeContext {
        let totalPoses = result.bindingModes.reduce(0) { $0 + $1.size }
        let totalWeight = result.bindingModes.reduce(0.0) { $0 + exp(-$1.freeEnergy / (kBkcal * result.temperature)) }

        let profiles: [ModeProfile] = result.bindingModes.prefix(3).enumerated().map { i, mode in
            let weight = totalWeight > 0 ? exp(-mode.freeEnergy / (kBkcal * result.temperature)) / totalWeight : 0
            return ModeProfile(
                index: i, poseCount: mode.size,
                freeEnergy: mode.freeEnergy, entropy: mode.entropy,
                enthalpy: mode.enthalpy, heatCapacity: mode.heatCapacity,
                boltzmannWeight: weight, temperature: result.temperature
            )
        }

        let entropies = result.bindingModes.map(\.entropy)
        let imbalance: Double = {
            guard entropies.count >= 2,
                  let maxS = entropies.max(),
                  let minS = entropies.filter({ $0 > 0 }).min(),
                  minS > 0 else { return 1.0 }
            return maxS / minS
        }()

        let dominantIdx = profiles.enumerated().max(by: { $0.element.boltzmannWeight < $1.element.boltzmannWeight })?.offset ?? 0

        return PreComputedModeContext(
            temperature: result.temperature,
            totalModes: result.bindingModes.count,
            totalPoses: totalPoses,
            globalFreeEnergy: result.globalThermodynamics.freeEnergy,
            modes: profiles,
            entropyImbalance: imbalance,
            dominantModeIndex: dominantIdx,
            scientificProvenance: result.globalThermodynamics.scientificProvenance
        )
    }

    /// Generate deterministic mode narrative.
    public func narrate(context: PreComputedModeContext) -> CrossPlatformModeNarrative {
        let bindingPhysical = context.claimValidity == .bindingPhysical
        let canonicalPhysical = context.claimValidity == .canonicalPhysical
        let descriptions = context.modes.map { mode -> CrossPlatformModeDescription in
            if !bindingPhysical {
                let characterization = canonicalPhysical
                    ? "Mode \(mode.index + 1): canonical-ensemble cluster with \(mode.poseCount) poses, F = \(String(format: "%.1f", mode.freeEnergy)) kcal/mol. This is not a binding-affinity claim."
                    : "Mode \(mode.index + 1): optimizer cluster diagnostic with \(mode.poseCount) poses and input score \(String(format: "%.1f", mode.freeEnergy))."
                return CrossPlatformModeDescription(
                    characterization: characterization,
                    optimizationHint: "Molecular and SAR optimization unavailable without binding_physical provenance."
                )
            }

            let driving: String
            if mode.isEntropyDriven {
                driving = "entropy-driven (conformational flexibility stabilizes this geometry)"
            } else if mode.isEnthalpyDriven {
                driving = "enthalpy-driven (direct interactions dominate stability)"
            } else {
                driving = "balanced enthalpy-entropy contributions"
            }

            let tightness = mode.poseCount < 5 ? "tight cluster" : mode.poseCount < 20 ? "moderate cluster" : "broad ensemble"

            let characterization = "Mode \(mode.index + 1): \(tightness) (\(mode.poseCount) poses), F = \(String(format: "%.1f", mode.freeEnergy)) kcal/mol, \(driving). Boltzmann weight \(String(format: "%.0f", mode.boltzmannWeight * 100))%."

            let hint: String
            if mode.isEntropyDriven {
                hint = "Rigidify ligand to shift from entropy- to enthalpy-driven binding for selectivity."
            } else if mode.isEnthalpyDriven {
                hint = "Add flexible substituents to improve entropy contribution if affinity needs boosting."
            } else {
                hint = "Well-balanced mode — optimize interaction geometry for potency."
            }

            return CrossPlatformModeDescription(characterization: characterization, optimizationHint: hint)
        }

        let dominant = context.modes[safe: context.dominantModeIndex]
        let insight: String
        if let dom = dominant {
            insight = bindingPhysical
                ? "Mode \(dom.index + 1) dominates (\(String(format: "%.0f", dom.boltzmannWeight * 100))% weight). Focus SAR optimization on this binding geometry."
                : "Mode \(dom.index + 1) dominates the diagnostic weight (\(String(format: "%.0f", dom.boltzmannWeight * 100))%). Binding-mode and SAR prioritization are unavailable."
        } else {
            insight = "No dominant mode — population is diverse across geometries."
        }

        let confidence: Double = context.modes.count >= 2 ? 0.8 : 0.5
        return CrossPlatformModeNarrative(modeDescriptions: descriptions, selectivityInsight: insight, confidence: confidence)
    }
}

// MARK: - Array Safe Index

private extension Array {
    subscript(safe index: Index) -> Element? {
        indices.contains(index) ? self[index] : nil
    }
}
