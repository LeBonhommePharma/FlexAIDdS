// ConvergenceCoach.swift — GA convergence monitoring and advice
//
// Monitors GA generational progress and advises when to stop, adjust
// parameters, or restart. Uses Apple FoundationModels for nuanced
// interpretation of convergence patterns.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed GA Progress

/// Snapshot of GA generational progress for LLM interpretation.
public struct GAProgressSnapshot: Sendable, Codable {
    /// Evidence governing the energy domain of the fitness values below.
    /// GA fitness is a raw optimizer objective, so this fails closed to
    /// proxy-only unless a calibrated ensemble record is supplied.
    public let scientificProvenance: ScientificProvenance?

    /// Strongest claim supported by the declared provenance.
    public var claimValidity: ClaimValidity {
        scientificProvenance?.claimValidity ?? .proxyOnly
    }

    /// Unit label for the fitness values. GA fitness is a contact-function
    /// score; kcal/mol is spoken only for calibrated, receipted ensembles.
    public var fitnessUnitLabel: String {
        claimValidity == .proxyOnly ? "CF arbitrary units" : "kcal/mol"
    }

    /// Current generation number
    public let currentGeneration: Int
    /// Maximum generations configured
    public let maxGenerations: Int
    /// Best fitness (CF score, negative = better) at current generation,
    /// expressed in the declared energy domain.
    public let bestFitness: Double
    /// Mean fitness of current population in the declared energy domain.
    public let meanFitness: Double
    /// Population diversity (Shannon entropy of fitness distribution)
    public let populationDiversity: Double
    /// Generations since last improvement in best fitness
    public let generationsSinceImprovement: Int
    /// Best fitness trajectory (last 10 generations, oldest first)
    public let fitnessTrajectory: [Double]
    /// Diversity trajectory (last 10 generations, oldest first)
    public let diversityTrajectory: [Double]
    /// Whether fitness is still improving (slope of recent trajectory)
    public let isImproving: Bool
    /// Whether diversity has collapsed (below threshold)
    public let isDiversityCollapsed: Bool
    /// Population size (number of chromosomes)
    public let populationSize: Int

    public init(currentGeneration: Int, maxGenerations: Int, bestFitness: Double,
                meanFitness: Double, populationDiversity: Double,
                generationsSinceImprovement: Int, fitnessTrajectory: [Double],
                diversityTrajectory: [Double], isImproving: Bool,
                isDiversityCollapsed: Bool, populationSize: Int,
                scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.currentGeneration = currentGeneration
        self.maxGenerations = maxGenerations
        self.bestFitness = bestFitness
        self.meanFitness = meanFitness
        self.populationDiversity = populationDiversity
        self.generationsSinceImprovement = generationsSinceImprovement
        self.fitnessTrajectory = fitnessTrajectory
        self.diversityTrajectory = diversityTrajectory
        self.isImproving = isImproving
        self.isDiversityCollapsed = isDiversityCollapsed
        self.populationSize = populationSize
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Advice from the convergence coach.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public enum GAAdvice: String, Sendable, Codable {
    case continueRun
    case stopEarly
    case increasePopulation
    case increaseMutationRate
    case restart
}

/// Full convergence coaching result.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct ConvergenceCoaching: Sendable, Codable {
    /// Recommended action
    public var advice: GAAdvice
    /// Explanation for the recommendation
    public var reasoning: String
    /// Estimated generations remaining to convergence (nil if unknown)
    public var estimatedGenerationsRemaining: Int?
    /// Confidence in the recommendation (0.0-1.0)
    public var confidence: Double
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor ConvergenceCoachActor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a genetic algorithm convergence advisor for molecular docking. \
        All fitness and diversity values are pre-computed — DO NOT perform arithmetic. \
        Assess whether the GA should continue, stop early, or adjust parameters. \
        Key signals: stagnation (no improvement for many generations), \
        diversity collapse (premature convergence), and fitness plateau. \
        Be concise and actionable.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Get convergence advice for the current GA state.
    public func coach(snapshot: GAProgressSnapshot) async throws -> ConvergenceCoaching {
        var prompt = buildPrompt(snapshot: snapshot)
        if estimateTokenCount(prompt) > 3800 {
            // Truncate trajectory to last 5 entries
            let truncated = GAProgressSnapshot(
                currentGeneration: snapshot.currentGeneration,
                maxGenerations: snapshot.maxGenerations,
                bestFitness: snapshot.bestFitness,
                meanFitness: snapshot.meanFitness,
                populationDiversity: snapshot.populationDiversity,
                generationsSinceImprovement: snapshot.generationsSinceImprovement,
                fitnessTrajectory: Array(snapshot.fitnessTrajectory.suffix(5)),
                diversityTrajectory: Array(snapshot.diversityTrajectory.suffix(5)),
                isImproving: snapshot.isImproving,
                isDiversityCollapsed: snapshot.isDiversityCollapsed,
                populationSize: snapshot.populationSize,
                scientificProvenance: snapshot.scientificProvenance
            )
            prompt = buildPrompt(snapshot: truncated)
        }
        return try await session.respond(to: prompt, generating: ConvergenceCoaching.self).content
    }

    private func estimateTokenCount(_ text: String) -> Int {
        Int(ceil(Double(text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count) * 1.3))
    }

    private func buildPrompt(snapshot: GAProgressSnapshot) -> String {
        // GA fitness is the raw optimizer objective. Never label it kcal/mol
        // unless the supplied provenance authorizes a calibrated energy scale.
        let units = snapshot.fitnessUnitLabel
        var p = "Assess this GA run. Produce ConvergenceCoaching.\n"
        p += "Claim validity: \(snapshot.claimValidity.rawValue) (fitness is an optimizer objective, not an affinity)\n"
        p += "Generation: \(snapshot.currentGeneration)/\(snapshot.maxGenerations)\n"
        p += "Best fitness: \(String(format: "%.2f", snapshot.bestFitness)) \(units)\n"
        p += "Mean fitness: \(String(format: "%.2f", snapshot.meanFitness)) \(units)\n"
        p += "Population: \(snapshot.populationSize) chromosomes\n"
        p += "Diversity: \(String(format: "%.3f", snapshot.populationDiversity))"
        p += snapshot.isDiversityCollapsed ? " [COLLAPSED]" : ""
        p += "\nStagnation: \(snapshot.generationsSinceImprovement) generations without improvement"
        p += "\nImproving: \(snapshot.isImproving ? "YES" : "NO")"

        if !snapshot.fitnessTrajectory.isEmpty {
            p += "\nFitness trajectory (last \(snapshot.fitnessTrajectory.count)): "
            p += snapshot.fitnessTrajectory.map { String(format: "%.1f", $0) }.joined(separator: " → ")
        }

        if snapshot.isDiversityCollapsed {
            p += "\nFLAG: Diversity collapsed — premature convergence risk"
        }
        if snapshot.generationsSinceImprovement > snapshot.maxGenerations / 4 {
            p += "\nFLAG: Stagnated for >\(snapshot.maxGenerations / 4) generations"
        }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent GA advice.
public enum CrossPlatformGAAdvice: String, Sendable, Codable {
    case continueRun, stopEarly, increasePopulation, increaseMutationRate, restart
}

/// Platform-independent convergence coaching.
public struct CrossPlatformConvergenceCoaching: Sendable, Codable {
    public let advice: CrossPlatformGAAdvice
    public let reasoning: String
    public let estimatedGenerationsRemaining: Int?
    public let confidence: Double

    public init(advice: CrossPlatformGAAdvice, reasoning: String,
                estimatedGenerationsRemaining: Int?, confidence: Double) {
        self.advice = advice
        self.reasoning = reasoning
        self.estimatedGenerationsRemaining = estimatedGenerationsRemaining
        self.confidence = confidence
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic convergence coach for non-Apple platforms.
public struct RuleBasedConvergenceCoach: Sendable {

    public init() {}

    /// Assess GA convergence using threshold logic.
    public func coach(snapshot: GAProgressSnapshot) -> CrossPlatformConvergenceCoaching {
        // Same gate as the model path: the fitness number keeps its value, but
        // its unit label is only physical when the provenance says so.
        let units = snapshot.fitnessUnitLabel
        let progress = Double(snapshot.currentGeneration) / Double(max(snapshot.maxGenerations, 1))
        let stagnationRatio = Double(snapshot.generationsSinceImprovement) / Double(max(snapshot.maxGenerations, 1))

        // Early run: always continue
        if progress < 0.2 {
            return CrossPlatformConvergenceCoaching(
                advice: .continueRun,
                reasoning: "Run is early (\(String(format: "%.0f", progress * 100))% complete). Allow more exploration time.",
                estimatedGenerationsRemaining: snapshot.maxGenerations - snapshot.currentGeneration,
                confidence: 0.7
            )
        }

        // Diversity collapsed + stagnated → restart or increase mutation
        if snapshot.isDiversityCollapsed && stagnationRatio > 0.15 {
            if snapshot.populationSize < 200 {
                return CrossPlatformConvergenceCoaching(
                    advice: .increasePopulation,
                    reasoning: "Diversity collapsed with only \(snapshot.populationSize) chromosomes. Increase population to maintain exploration.",
                    estimatedGenerationsRemaining: nil,
                    confidence: 0.8
                )
            }
            return CrossPlatformConvergenceCoaching(
                advice: .increaseMutationRate,
                reasoning: "Diversity collapsed at generation \(snapshot.currentGeneration) despite \(snapshot.populationSize) chromosomes. Increase mutation rate to escape local minimum.",
                estimatedGenerationsRemaining: nil,
                confidence: 0.75
            )
        }

        // Long stagnation → stop early
        if stagnationRatio > 0.25 && !snapshot.isImproving {
            return CrossPlatformConvergenceCoaching(
                advice: .stopEarly,
                reasoning: "No improvement for \(snapshot.generationsSinceImprovement) generations (\(String(format: "%.0f", stagnationRatio * 100))% of run). Best fitness \(String(format: "%.2f", snapshot.bestFitness)) \(units) is likely the search optimum.",
                estimatedGenerationsRemaining: 0,
                confidence: 0.85
            )
        }

        // Still improving → continue
        if snapshot.isImproving {
            let remaining = snapshot.maxGenerations - snapshot.currentGeneration
            return CrossPlatformConvergenceCoaching(
                advice: .continueRun,
                reasoning: "Fitness still improving. Best: \(String(format: "%.2f", snapshot.bestFitness)) \(units). \(remaining) generations remaining.",
                estimatedGenerationsRemaining: remaining,
                confidence: 0.8
            )
        }

        // Default: continue but with low confidence
        return CrossPlatformConvergenceCoaching(
            advice: .continueRun,
            reasoning: "Run at \(String(format: "%.0f", progress * 100))% progress. Fitness plateaued but diversity maintained — may still improve.",
            estimatedGenerationsRemaining: snapshot.maxGenerations - snapshot.currentGeneration,
            confidence: 0.5
        )
    }
}
