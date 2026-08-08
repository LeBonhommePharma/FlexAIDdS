// SelectivityAnalyst.swift — Multi-target selectivity analysis
//
// Compares docking results across protein targets to explain selectivity.
// Uses Apple FoundationModels for nuanced interpretation of ΔΔG,
// entropy differences, and structural drivers.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Target Summary

/// Summary of docking results for a single target.
public struct TargetDockingSummary: Sendable, Codable {
    /// Evidence governing whether this target supports binding/selectivity claims.
    public let scientificProvenance: ScientificProvenance?
    /// Target name (e.g., "5HT2A", "D2R")
    public let targetName: String
    /// Best free energy across modes (kcal/mol)
    public let bestFreeEnergy: Double
    /// Number of binding modes
    public let modeCount: Int
    /// Configurational entropy (nats)
    public let sConf: Double
    /// Vibrational entropy (kcal/mol/K)
    public let sVib: Double
    /// Whether entropy converged
    public let isConverged: Bool
    /// Cavity volume (ų, if known)
    public let cavityVolume: Double?
    /// Population size
    public let populationSize: Int

    public init(targetName: String, bestFreeEnergy: Double, modeCount: Int,
                sConf: Double, sVib: Double, isConverged: Bool,
                cavityVolume: Double?, populationSize: Int,
                scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.targetName = targetName
        self.bestFreeEnergy = bestFreeEnergy
        self.modeCount = modeCount
        self.sConf = sConf
        self.sVib = sVib
        self.isConverged = isConverged
        self.cavityVolume = cavityVolume
        self.populationSize = populationSize
    }
}

/// Pre-computed pairwise ΔΔG between two targets.
public struct DeltaDeltaG: Sendable, Codable, Equatable {
    public let targetA: String
    public let targetB: String
    public let ddg: Double

    public init(targetA: String, targetB: String, ddg: Double) {
        self.targetA = targetA
        self.targetB = targetB
        self.ddg = ddg
    }
}

/// Context for multi-target selectivity analysis.
public struct SelectivityContext: Sendable, Codable {
    /// Ligand identifier
    public let ligandName: String
    /// Target summaries (typically 2-4 targets)
    public let targets: [TargetDockingSummary]
    /// Pre-computed pairwise ΔΔG values
    public let deltaDeltas: [DeltaDeltaG]

    public var allowsBindingClaims: Bool {
        targets.count >= 2 && targets.allSatisfy {
            $0.scientificProvenance?.claimValidity == .bindingPhysical
        }
    }

    /// Weakest claim across the compared targets. A selectivity statement is
    /// only as strong as the least-evidenced target in the comparison.
    public var claimValidity: ClaimValidity {
        guard targets.count >= 2 else { return .proxyOnly }
        if allowsBindingClaims { return .bindingPhysical }
        if targets.allSatisfy({
            let validity = $0.scientificProvenance?.claimValidity ?? .proxyOnly
            return validity == .canonicalPhysical || validity == .bindingPhysical
        }) {
            return .canonicalPhysical
        }
        return .proxyOnly
    }

    public init(ligandName: String, targets: [TargetDockingSummary]) {
        self.ligandName = ligandName
        self.targets = targets
        // Pre-compute pairwise ΔΔG
        var deltas: [DeltaDeltaG] = []
        for i in 0..<targets.count {
            for j in (i+1)..<targets.count {
                let ddg = targets[i].bestFreeEnergy - targets[j].bestFreeEnergy
                deltas.append(DeltaDeltaG(
                    targetA: targets[i].targetName,
                    targetB: targets[j].targetName,
                    ddg: ddg
                ))
            }
        }
        self.deltaDeltas = deltas
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// What drives selectivity between targets.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public enum SelectivityDriver: String, Sendable, Codable {
    case enthalpic
    case entropic
    case mixed
    case inconclusive
}

/// Selectivity analysis result.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct SelectivityAnalysis: Sendable, Codable {
    /// Preferred target name
    public var preferredTarget: String
    /// Retained numeric compatibility field: the raw score difference between
    /// the two best targets. It is a physical ΔΔG in kcal/mol ONLY when
    /// `deltaGClaimValidity` is `binding_physical`; otherwise it is an
    /// input-scale proxy diagnostic and must never be presented as affinity,
    /// potency, Kd or Ki.
    public var deltaG: Double
    /// Fail-closed marker for `deltaG`, written deterministically by the
    /// analyst after generation. The model never authors this value.
    public var deltaGClaimValidity: String
    /// What drives the selectivity
    public var driver: SelectivityDriver
    /// Plain-English explanation
    public var explanation: String
    /// Design suggestion for modulating selectivity
    public var designSuggestion: String
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor SelectivityAnalystActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a pharmacology advisor analyzing drug-target selectivity from \
        molecular docking results. All energies and entropies are pre-computed. \
        Explain which target is preferred, what drives selectivity (enthalpy vs entropy), \
        and how to modify the ligand only when every target is binding_physical. \
        Otherwise state that physical selectivity and optimization are unavailable. Be concise and cite \
        specific ΔΔG values. This is for psychopharmacology (serotonin/dopamine receptors).
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Analyze selectivity across targets.
    public func analyze(context: SelectivityContext) async throws -> SelectivityAnalysis {
        if !context.allowsBindingClaims {
            return SelectivityAnalysis(
                preferredTarget: "unavailable",
                deltaG: context.deltaDeltas.first?.ddg ?? 0,
                deltaGClaimValidity: context.claimValidity.rawValue,
                driver: .inconclusive,
                explanation: "Physical target selectivity and ΔΔG are unavailable; retained score differences are proxy diagnostics only.",
                designSuggestion: "Molecular and selectivity optimization are unavailable until every target has binding_physical provenance."
            )
        }
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            // Truncate to first 2 targets
            let truncated = SelectivityContext(
                ligandName: context.ligandName,
                targets: Array(context.targets.prefix(2))
            )
            prompt = buildPrompt(context: truncated)
        }
        var analysis = try await session.respond(
            to: prompt, generating: SelectivityAnalysis.self).content
        // Overwrite the generated marker: claim validity is decided by the
        // evidence gate above, never by the language model.
        analysis.deltaGClaimValidity = context.claimValidity.rawValue
        return analysis
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(context: SelectivityContext) -> String {
        var p = "Analyze selectivity for \(context.ligandName). Produce a SelectivityAnalysis.\n"

        for target in context.targets {
            p += "\n\(target.targetName): F=\(String(format: "%.2f", target.bestFreeEnergy))kcal/mol, "
            p += "\(target.modeCount) modes, S_conf=\(String(format: "%.4f", target.sConf))nats, "
            p += "S_vib=\(String(format: "%.6f", target.sVib))kcal/mol/K, "
            p += "converged=\(target.isConverged)"
            if let vol = target.cavityVolume {
                p += ", pocket=\(String(format: "%.0f", vol))ų"
            }
        }

        for delta in context.deltaDeltas {
            p += "\nΔΔG(\(delta.targetA) vs \(delta.targetB)) = \(String(format: "%.2f", delta.ddg)) kcal/mol"
        }

        // Entropy-driven selectivity check
        if context.targets.count >= 2 {
            let sorted = context.targets.sorted { $0.bestFreeEnergy < $1.bestFreeEnergy }
            let best = sorted[0], second = sorted[1]
            if best.sConf > second.sConf * 1.5 {
                p += "\nFLAG: Preferred target has higher S_conf — selectivity may be entropy-driven"
            }
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent selectivity driver.
public enum CrossPlatformSelectivityDriver: String, Sendable, Codable {
    case enthalpic, entropic, mixed, inconclusive
}

/// Platform-independent selectivity analysis.
public struct CrossPlatformSelectivityAnalysis: Sendable, Codable {
    public let preferredTarget: String

    /// Retained numeric compatibility field: the raw score difference between
    /// the two best targets. Physical ΔΔG in kcal/mol ONLY when
    /// `deltaGClaimValidity == .bindingPhysical`; otherwise an input-scale
    /// proxy diagnostic that must never be read as affinity or potency.
    public let deltaG: Double

    /// Fail-closed claim marker travelling with `deltaG`. Defaults to
    /// proxy-only, including for records decoded from older payloads that
    /// predate this field.
    public let deltaGClaimValidity: ClaimValidity

    /// Whether `deltaG` may be presented as a physical binding quantity.
    public var isDeltaGPhysical: Bool { deltaGClaimValidity == .bindingPhysical }

    /// Presentation-safe unit label for `deltaG`.
    public var deltaGUnitLabel: String {
        isDeltaGPhysical ? "kcal/mol" : "input-scale units (proxy diagnostic)"
    }

    public let driver: CrossPlatformSelectivityDriver
    public let explanation: String
    public let designSuggestion: String

    public init(preferredTarget: String, deltaG: Double,
                driver: CrossPlatformSelectivityDriver,
                explanation: String, designSuggestion: String,
                deltaGClaimValidity: ClaimValidity = .proxyOnly) {
        self.preferredTarget = preferredTarget
        self.deltaG = deltaG
        self.deltaGClaimValidity = deltaGClaimValidity
        self.driver = driver
        self.explanation = explanation
        self.designSuggestion = designSuggestion
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        preferredTarget = try container.decode(String.self, forKey: .preferredTarget)
        deltaG = try container.decode(Double.self, forKey: .deltaG)
        driver = try container.decode(
            CrossPlatformSelectivityDriver.self, forKey: .driver)
        explanation = try container.decode(String.self, forKey: .explanation)
        designSuggestion = try container.decode(String.self, forKey: .designSuggestion)
        // A missing or hostile marker must never upgrade the retained number.
        deltaGClaimValidity = (try? container.decodeIfPresent(
            ClaimValidity.self, forKey: .deltaGClaimValidity)) ?? .proxyOnly
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic selectivity analyst for non-Apple platforms.
public struct RuleBasedSelectivityAnalyst: Sendable {

    // Named thresholds for selectivity determination
    private static let inconclusiveEnthalpicThreshold = 0.5  // kcal/mol
    private static let inconclusiveEntropicThreshold = 0.0001  // kcal/mol/K
    private static let entropicDominanceRatio = 1.3
    private static let enthalpicSignificanceThreshold = 2.0  // kcal/mol

    public init() {}

    /// Analyze selectivity using threshold logic.
    public func analyze(context: SelectivityContext) -> CrossPlatformSelectivityAnalysis {
        guard context.targets.count >= 2 else {
            return CrossPlatformSelectivityAnalysis(
                preferredTarget: context.targets.first?.targetName ?? "unknown",
                deltaG: 0,
                driver: .inconclusive,
                explanation: "Only one target available — selectivity analysis requires at least two targets.",
                designSuggestion: "Dock against additional targets for selectivity analysis.",
                deltaGClaimValidity: .proxyOnly
            )
        }

        let sorted = context.targets.sorted { $0.bestFreeEnergy < $1.bestFreeEnergy }
        let preferred = sorted[0]
        let second = sorted[1]
        let ddg = preferred.bestFreeEnergy - second.bestFreeEnergy

        guard context.allowsBindingClaims else {
            return CrossPlatformSelectivityAnalysis(
                preferredTarget: "unavailable",
                deltaG: ddg,
                driver: .inconclusive,
                explanation: "Physical target selectivity and ΔΔG are unavailable. Raw score difference \(String(format: "%.2f", ddg)) is retained as a proxy input-scale diagnostic only.",
                designSuggestion: "Molecular and selectivity optimization are unavailable until every target has binding_physical provenance.",
                deltaGClaimValidity: context.claimValidity
            )
        }

        // Determine driver
        let kB = 0.001987206
        let sConfPhysA = preferred.sConf * kB
        let sConfPhysB = second.sConf * kB
        let entropicDiff = abs(sConfPhysA - sConfPhysB)
        let enthalpicDiff = abs(ddg)

        let driver: CrossPlatformSelectivityDriver
        if enthalpicDiff < 0.5 && entropicDiff < 0.0001 {
            driver = .inconclusive
        } else if entropicDiff > enthalpicDiff * 0.001 && preferred.sConf > second.sConf * 1.3 {
            driver = .entropic
        } else if enthalpicDiff > 2.0 {
            driver = .enthalpic
        } else {
            driver = .mixed
        }

        // Explanation
        var explanation = "\(context.ligandName) prefers \(preferred.targetName) (F = \(String(format: "%.1f", preferred.bestFreeEnergy)) kcal/mol) over \(second.targetName) (F = \(String(format: "%.1f", second.bestFreeEnergy)) kcal/mol). "
        explanation += "ΔΔG = \(String(format: "%.2f", ddg)) kcal/mol. "

        switch driver {
        case .entropic:
            explanation += "Selectivity is entropy-driven: \(preferred.targetName) has broader conformational ensemble (S_conf = \(String(format: "%.4f", preferred.sConf)) vs \(String(format: "%.4f", second.sConf)) nats)."
        case .enthalpic:
            explanation += "Selectivity is enthalpy-driven: stronger direct interactions at \(preferred.targetName)."
        case .mixed:
            explanation += "Both enthalpy and entropy contribute to selectivity."
        case .inconclusive:
            explanation += "Selectivity is marginal — results may not be significant."
        }

        // Design suggestion
        let suggestion: String
        switch driver {
        case .entropic:
            suggestion = "Rigidify ligand to reduce entropy-driven selectivity for \(preferred.targetName), or add flexible groups to enhance selectivity."
        case .enthalpic:
            suggestion = "Optimize interaction geometry at \(second.targetName) binding site to improve affinity and shift selectivity."
        case .mixed:
            suggestion = "Balanced selectivity — consider fragment-based approach targeting unique pocket features of each receptor."
        case .inconclusive:
            suggestion = "Increase sampling (more GA generations) and verify convergence before drawing selectivity conclusions."
        }

        let convergedBoth = preferred.isConverged && second.isConverged
        if !convergedBoth {
            let notConverged = [preferred, second].filter { !$0.isConverged }.map(\.targetName)
            return CrossPlatformSelectivityAnalysis(
                preferredTarget: preferred.targetName,
                deltaG: ddg,
                driver: .inconclusive,
                explanation: explanation + " Note: \(notConverged.joined(separator: ", ")) not converged — selectivity assessment unreliable.",
                designSuggestion: "Achieve convergence at all targets before selectivity-driven optimization.",
                deltaGClaimValidity: context.claimValidity
            )
        }

        return CrossPlatformSelectivityAnalysis(
            preferredTarget: preferred.targetName,
            deltaG: ddg,
            driver: driver,
            explanation: explanation,
            designSuggestion: suggestion,
            deltaGClaimValidity: context.claimValidity
        )
    }
}
