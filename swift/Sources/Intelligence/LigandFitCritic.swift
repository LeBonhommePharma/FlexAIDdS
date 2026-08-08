// LigandFitCritic.swift — Pose quality annotator for top-ranked docking poses
//
// Evaluates the quality of top-ranked poses using Boltzmann weight distribution,
// CF score alignment, and pose diversity. Produces medicinal-chemistry-friendly
// annotations using Apple FoundationModels.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Pose Context

/// Summary of a single pose for LLM interpretation.
public struct PoseProfile: Sendable, Codable {
    /// Rank within the binding mode (0 = best)
    public let rank: Int
    /// Complementarity-function score in the declared energy domain
    /// (arbitrary CF units unless the context provenance says otherwise).
    public let cfScore: Double
    /// Boltzmann weight (fraction, 0.0-1.0)
    public let boltzmannWeight: Double
    /// RMSD to the mode centroid (Angstroms)
    public let rmsdToCentroid: Double

    public init(rank: Int, cfScore: Double, boltzmannWeight: Double, rmsdToCentroid: Double) {
        self.rank = rank
        self.cfScore = cfScore
        self.boltzmannWeight = boltzmannWeight
        self.rmsdToCentroid = rmsdToCentroid
    }
}

/// Full pose quality context for a binding mode.
public struct PoseQualityContext: Sendable, Codable {
    /// Evidence governing whether these pose statistics may be presented as
    /// physical energies or used for medicinal-chemistry decisions.
    public let scientificProvenance: ScientificProvenance?

    /// Strongest claim supported by the declared provenance. Missing evidence
    /// fails closed to proxy-only.
    public var claimValidity: ClaimValidity {
        scientificProvenance?.claimValidity ?? .proxyOnly
    }

    /// Unit label for `cfScore` / `modeFreeEnergy`. kcal/mol is spoken only
    /// when a calibrated, measure-receipted ensemble authorizes it.
    public var energyUnitLabel: String {
        claimValidity == .proxyOnly
            ? "CF arbitrary units"
            : "kcal/mol"
    }

    /// Binding mode index
    public let modeIndex: Int
    /// Top poses (up to 5), sorted by Boltzmann weight descending
    public let topPoses: [PoseProfile]
    /// Total poses in this mode
    public let totalPoses: Int
    /// Weight of the dominant pose (highest Boltzmann weight)
    public let dominantPoseWeight: Double
    /// Whether the top pose dominates (weight > 0.5)
    public let hasDominantPose: Bool
    /// RMSD spread of top 5 poses (max - min RMSD to centroid)
    public let rmsdSpread: Double
    /// Whether CF score ranking matches Boltzmann weight ranking
    public let scoreWeightAligned: Bool
    /// Spearman rank correlation between CF score and Boltzmann weight (top 5)
    public let scoreWeightCorrelation: Double
    /// Mode ensemble free-energy-like value in the declared energy domain.
    public let modeFreeEnergy: Double

    public init(modeIndex: Int, topPoses: [PoseProfile], totalPoses: Int,
                modeFreeEnergy: Double,
                scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.modeIndex = modeIndex
        self.topPoses = topPoses
        self.totalPoses = totalPoses
        self.modeFreeEnergy = modeFreeEnergy

        self.dominantPoseWeight = topPoses.first?.boltzmannWeight ?? 0
        self.hasDominantPose = self.dominantPoseWeight > 0.5

        let rmsds = topPoses.map(\.rmsdToCentroid)
        self.rmsdSpread = (rmsds.max() ?? 0) - (rmsds.min() ?? 0)

        // Check if CF rank matches Boltzmann rank
        let cfRanked = topPoses.sorted { $0.cfScore < $1.cfScore }
        let weightRanked = topPoses.sorted { $0.boltzmannWeight > $1.boltzmannWeight }
        if topPoses.count >= 2 {
            self.scoreWeightAligned = cfRanked.first?.rank == weightRanked.first?.rank
        } else {
            self.scoreWeightAligned = true
        }

        // Simplified rank correlation (Spearman-like)
        // Uses index-based ranking to avoid dictionary key collisions
        if topPoses.count >= 2 {
            let cfSorted = topPoses.sorted { $0.cfScore < $1.cfScore }
            let wSorted = topPoses.sorted { $0.boltzmannWeight > $1.boltzmannWeight }
            var d2Sum = 0.0
            let n = Double(topPoses.count)
            for pose in topPoses {
                let cfR = cfSorted.firstIndex(where: { $0.rank == pose.rank }) ?? 0
                let wR = wSorted.firstIndex(where: { $0.rank == pose.rank }) ?? 0
                let d = Double(cfR - wR)
                d2Sum += d * d
            }
            self.scoreWeightCorrelation = 1.0 - (6.0 * d2Sum) / (n * (n * n - 1))
        } else {
            self.scoreWeightCorrelation = 1.0
        }
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Pose quality report.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct PoseQualityReport: Sendable, Codable {
    /// Summary of the top pose
    public var topPoseSummary: String
    /// Consensus assessment
    public var poseConsensus: String
    /// Score-weight alignment assessment
    public var scoreWeightAlignment: String
    /// Confidence in the top pose prediction (0.0-1.0)
    public var confidenceInTopPose: Double
    /// Practical advice for the medicinal chemist
    public var medicinalChemistryNote: String
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor LigandFitCriticActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a computational chemistry advisor evaluating docking pose quality. \
        All scores and weights are pre-computed — DO NOT calculate. \
        Assess whether the top-ranked pose is reliable enough for medicinal chemistry \
        decisions. Consider Boltzmann weight dominance, CF score agreement, \
        RMSD spread, and pose consensus. Provide practical advice for synthesis \
        prioritization only for binding_physical inputs; otherwise report \
        sampling diagnostics and state that affinity guidance is unavailable. \
        Be concise and honest about uncertainty.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Evaluate pose quality for a binding mode.
    ///
    /// Synthesis prioritization and structure-based design are binding claims,
    /// so the model is only consulted for `binding_physical` provenance. Every
    /// other case returns a deterministic sampling diagnostic.
    public func evaluate(context: PoseQualityContext) async throws -> PoseQualityReport {
        guard context.claimValidity == .bindingPhysical else {
            let domain = context.claimValidity == .canonicalPhysical
                ? "canonical-ensemble"
                : "optimizer/input-scale"
            return PoseQualityReport(
                topPoseSummary: "Mode \(context.modeIndex + 1): \(domain) pose-clustering diagnostic over \(context.totalPoses) poses; pose scores are not binding affinities.",
                poseConsensus: context.hasDominantPose ? "dominant-sample" : "distributed-samples",
                scoreWeightAlignment: "Score/weight agreement is reported as a sampling diagnostic only.",
                confidenceInTopPose: 0.0,
                medicinalChemistryNote: "Medicinal-chemistry, synthesis-prioritization and structure-based-design guidance are unavailable without a calibrated matched association cycle."
            )
        }
        var prompt = buildPrompt(context: context)
        if estimateTokenCount(prompt) > 3800 {
            prompt = buildPrompt(context: PoseQualityContext(
                modeIndex: context.modeIndex,
                topPoses: Array(context.topPoses.prefix(3)),
                totalPoses: context.totalPoses,
                modeFreeEnergy: context.modeFreeEnergy,
                scientificProvenance: context.scientificProvenance
            ))
        }
        return try await session.respond(to: prompt, generating: PoseQualityReport.self).content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(context: PoseQualityContext) -> String {
        let units = context.energyUnitLabel
        var p = "Evaluate pose quality for mode \(context.modeIndex + 1). Produce a PoseQualityReport.\n"
        p += "Claim validity: \(context.claimValidity.rawValue)\n"
        p += "Mode F=\(String(format: "%.2f", context.modeFreeEnergy)) \(units), \(context.totalPoses) poses total\n"
        p += "RMSD spread: \(String(format: "%.1f", context.rmsdSpread))Å\n"
        p += "Score-weight correlation: \(String(format: "%.2f", context.scoreWeightCorrelation))\n"
        p += "Score-weight aligned: \(context.scoreWeightAligned ? "YES" : "NO")\n"
        if context.hasDominantPose {
            p += "Dominant pose: \(String(format: "%.0f", context.dominantPoseWeight * 100))% Boltzmann weight\n"
        } else {
            p += "No dominant pose — weight distributed across ensemble\n"
        }

        p += "\nTop poses:"
        for pose in context.topPoses {
            p += "\n  Rank \(pose.rank + 1): CF=\(String(format: "%.1f", pose.cfScore)) \(units), "
            p += "w=\(String(format: "%.1f", pose.boltzmannWeight * 100))%, "
            p += "RMSD=\(String(format: "%.1f", pose.rmsdToCentroid))Å"
        }

        if !context.scoreWeightAligned {
            p += "\nFLAG: Best CF score and highest Boltzmann weight disagree — entropy fights enthalpy"
        }
        if context.rmsdSpread > 3.0 {
            p += "\nFLAG: Large RMSD spread (\(String(format: "%.1f", context.rmsdSpread))Å) — ligand orientation not fully resolved"
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent pose quality report.
public struct CrossPlatformPoseQualityReport: Sendable, Codable {
    public let topPoseSummary: String
    public let poseConsensus: String
    public let scoreWeightAlignment: String
    public let confidenceInTopPose: Double
    public let medicinalChemistryNote: String

    public init(topPoseSummary: String, poseConsensus: String,
                scoreWeightAlignment: String, confidenceInTopPose: Double,
                medicinalChemistryNote: String) {
        self.topPoseSummary = topPoseSummary
        self.poseConsensus = poseConsensus
        self.scoreWeightAlignment = scoreWeightAlignment
        self.confidenceInTopPose = confidenceInTopPose
        self.medicinalChemistryNote = medicinalChemistryNote
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic pose quality critic for non-Apple platforms.
public struct RuleBasedLigandFitCritic: Sendable {

    // Named thresholds for deterministic logic
    private static let strongConsensusRMSD = 2.0
    private static let moderateConsensusRMSD = 3.0
    private static let highCorrelation = 0.8
    private static let moderateCorrelation = 0.4
    private static let largeRMSDSpread = 3.0
    private static let maxConfidence = 0.95

    public init() {}

    /// Build a PoseQualityContext from a DockingResult binding mode.
    public func buildContext(
        modeIndex: Int,
        mode: BindingModeResult,
        poses: [PoseResult],
        temperature: Double
    ) -> PoseQualityContext {
        // Guard against invalid temperature (division by zero)
        let safeTemp = temperature > 0 ? temperature : 298.15
        // Compute Boltzmann weights from CF scores
        let beta = 1.0 / (kBkcal * safeTemp)
        let energies = poses.map(\.cf)
        let minE = energies.min() ?? 0
        let expWeights = energies.map { exp(-beta * ($0 - minE)) }
        let totalWeight = expWeights.reduce(0, +)
        let weights = totalWeight > 0 ? expWeights.map { $0 / totalWeight } : Array(repeating: 0.0, count: poses.count)

        // Build profiles sorted by weight
        var profiles: [(PoseProfile, Double)] = zip(poses.indices, zip(poses, weights)).map { i, pw in
            let (pose, weight) = pw
            return (PoseProfile(
                rank: i,
                cfScore: pose.cf,
                boltzmannWeight: weight,
                rmsdToCentroid: Double(pose.reachDist)
            ), weight)
        }
        profiles.sort { $0.1 > $1.1 }

        let topProfiles = Array(profiles.prefix(5).map(\.0))

        return PoseQualityContext(
            modeIndex: modeIndex,
            topPoses: topProfiles,
            totalPoses: poses.count,
            modeFreeEnergy: mode.freeEnergy,
            // Carry the mode's own evidence; a mode without thermodynamics
            // fails closed to proxy-only rather than borrowing a default.
            scientificProvenance: mode.thermodynamics?.scientificProvenance
        )
    }

    /// Evaluate pose quality using threshold logic.
    public func evaluate(context: PoseQualityContext) -> CrossPlatformPoseQualityReport {
        guard !context.topPoses.isEmpty else {
            return CrossPlatformPoseQualityReport(
                topPoseSummary: "No poses available for this binding mode.",
                poseConsensus: "weak",
                scoreWeightAlignment: "No data available.",
                confidenceInTopPose: 0.0,
                medicinalChemistryNote: "No poses to evaluate. Check docking parameters and re-run."
            )
        }
        let top = context.topPoses.first
        // Medicinal-chemistry and structure-based-design advice are binding
        // claims; kcal/mol requires at least a canonical-physical record.
        let bindingPhysical = context.claimValidity == .bindingPhysical
        let units = context.energyUnitLabel

        // Top pose summary
        let topSummary: String
        if let pose = top {
            topSummary = "Pose \(pose.rank + 1) (CF = \(String(format: "%.1f", pose.cfScore)) \(units)) carries \(String(format: "%.0f", pose.boltzmannWeight * 100))% Boltzmann weight at RMSD \(String(format: "%.1f", pose.rmsdToCentroid)) Å from mode centroid."
        } else {
            topSummary = "No poses available for this binding mode."
        }

        // Consensus
        let consensus: String
        if context.hasDominantPose && context.rmsdSpread < 2.0 {
            consensus = "strong"
        } else if context.hasDominantPose || context.rmsdSpread < 3.0 {
            consensus = "moderate"
        } else if context.topPoses.count >= 3 && !context.hasDominantPose {
            consensus = "ambiguous"
        } else {
            consensus = "weak"
        }

        // Score-weight alignment
        let alignment: String
        if context.scoreWeightCorrelation > 0.8 {
            alignment = "Well aligned — best CF scores correspond to highest Boltzmann weights. Entropy and enthalpy agree."
        } else if context.scoreWeightCorrelation > 0.4 {
            alignment = "Partially aligned — some entropy-enthalpy tension in pose ranking."
        } else {
            alignment = "Misaligned — entropy fights enthalpy. The energetically best pose is not the most thermodynamically populated. Consider ensemble-averaged analysis."
        }

        // Confidence
        var confidence = 0.5
        if context.hasDominantPose { confidence += 0.2 }
        if context.scoreWeightAligned { confidence += 0.15 }
        if context.rmsdSpread < 2.0 { confidence += 0.1 }
        if context.totalPoses >= 10 { confidence += 0.05 }
        confidence = min(max(confidence, 0.0), 0.95)

        // Medicinal chemistry note
        let note: String
        if !bindingPhysical {
            // No calibrated matched association cycle: the ensemble supports
            // sampling diagnostics, never synthesis or SAR prioritization.
            note = "Medicinal-chemistry guidance is unavailable (claim_validity=\(context.claimValidity.rawValue)). These clusters are sampling diagnostics in \(units); they are not binding affinities and must not drive synthesis or structure-based design."
        } else if consensus == "strong" && context.scoreWeightAligned {
            note = "High-confidence binding pose. The predicted geometry is suitable for structure-based design and interaction analysis."
        } else if consensus == "ambiguous" {
            note = "Multiple competing poses with similar weights. Consider whether the ligand has internal symmetry or if the binding site accommodates multiple orientations. Validate with experimental data before committing to a single geometry."
        } else if !context.scoreWeightAligned {
            note = "The energetically best pose differs from the thermodynamically most populated one. This suggests conformational entropy plays a significant role. Consider rigidifying the ligand scaffold to resolve the ambiguity."
        } else if context.rmsdSpread > 3.0 {
            note = "Large positional spread (\(String(format: "%.1f", context.rmsdSpread)) Å) indicates the ligand samples multiple orientations within this mode. The binding orientation is not fully resolved — increase sampling or constrain the docking."
        } else {
            note = "Moderate pose quality. Results are usable for preliminary SAR but verify key interactions with molecular dynamics or experimental data."
        }

        return CrossPlatformPoseQualityReport(
            topPoseSummary: topSummary,
            poseConsensus: consensus,
            scoreWeightAlignment: alignment,
            confidenceInTopPose: confidence,
            medicinalChemistryNote: note
        )
    }
}
