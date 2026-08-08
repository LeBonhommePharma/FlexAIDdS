// IntelligenceOracle.swift — Apple Intelligence integration for FlexAIDdS
//
// Uses the on-device FoundationModels framework (macOS 26+ / iOS 26+)
// to analyze BindingPopulation thermodynamics + health correlation.
// Returns structured analysis with confidence scoring, trend tracking,
// and multi-turn conversational context.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS
import FlexAIDCore
import HealthIntegration

// MARK: - Oracle Result

/// Confidence level for oracle analysis bullets.
public enum AnalysisConfidence: String, Sendable, Codable, CaseIterable {
    case high     // Strong signal, well-supported by data
    case moderate // Reasonable inference, some uncertainty
    case low      // Speculative, needs more data or sampling
}

/// A single analysis bullet with metadata.
public struct AnalysisBullet: Sendable, Codable {
    public let text: String
    public let confidence: AnalysisConfidence
    public let category: BulletCategory

    public enum BulletCategory: String, Sendable, Codable {
        case binding       // Free energy / affinity assessment
        case entropy       // Conformational or vibrational entropy
        case health        // HealthKit correlation
        case modification  // PTM / glycan effects
        case fleet         // Fleet compute recommendation
        case trend         // Historical comparison
        case provenance    // Scientific-domain / reference-state validity
    }

    public init(text: String, confidence: AnalysisConfidence, category: BulletCategory) {
        self.text = text
        self.confidence = confidence
        self.category = category
    }
}

/// Structured analysis from the Intelligence oracle.
public struct OracleAnalysis: Sendable, Codable {
    /// Analysis bullets (typically 3-5)
    public let bullets: [String]

    /// Structured bullets with confidence and category metadata
    public let structuredBullets: [AnalysisBullet]

    /// Timestamp of analysis
    public let timestamp: Date

    /// Input summary used for the analysis
    public let inputSummary: String

    /// Overall confidence of the analysis
    public let overallConfidence: AnalysisConfidence

    public init(bullets: [String], inputSummary: String) {
        self.bullets = bullets
        self.structuredBullets = bullets.map {
            AnalysisBullet(text: $0, confidence: .moderate, category: .binding)
        }
        self.timestamp = Date()
        self.inputSummary = inputSummary
        self.overallConfidence = .moderate
    }

    public init(structuredBullets: [AnalysisBullet], inputSummary: String, overallConfidence: AnalysisConfidence) {
        self.structuredBullets = structuredBullets
        self.bullets = structuredBullets.map(\.text)
        self.timestamp = Date()
        self.inputSummary = inputSummary
        self.overallConfidence = overallConfidence
    }
}

// MARK: - Analysis History (trend tracking)

/// Stores past analyses for trend detection across docking campaigns.
/// Uses per-campaign-key FIFO (max 10 per key, 50 total) to prevent unbounded growth.
/// Internally keyed by campaign key for O(1) per-key lookups and eviction.
public actor AnalysisHistory {
    public static let shared = AnalysisHistory()

    private struct Entry {
        let key: String
        let source: String
        let analysis: OracleAnalysis
    }

    /// Per-key FIFO queues for O(1) key-scoped operations
    private var keyedEntries: [String: [(source: String, analysis: OracleAnalysis)]] = [:]
    /// Global insertion-order FIFO for total-count eviction
    private var insertionOrder: [String] = []
    private let maxEntriesPerKey = 10
    private let maxTotalEntries = 50

    private init() {}

    /// Record an analysis keyed by receptor-ligand pair.
    /// - Parameters:
    ///   - key: Campaign identifier (e.g., "5HT2A-psilocin")
    ///   - analysis: The analysis result
    ///   - source: Origin of the analysis ("text" for IntelligenceOracle, "referee" for ThermoReferee)
    public func record(key: String, analysis: OracleAnalysis, source: String = "text") {
        keyedEntries[key, default: []].append((source: source, analysis: analysis))
        insertionOrder.append(key)

        // Per-key FIFO: keep at most maxEntriesPerKey per campaign key
        if let count = keyedEntries[key]?.count, count > maxEntriesPerKey {
            keyedEntries[key]?.removeFirst(count - maxEntriesPerKey)
        }

        // Global FIFO: cap total entries
        while totalCount > maxTotalEntries, let oldestKey = insertionOrder.first {
            insertionOrder.removeFirst()
            if keyedEntries[oldestKey]?.isEmpty == false {
                keyedEntries[oldestKey]?.removeFirst()
                if keyedEntries[oldestKey]?.isEmpty == true {
                    keyedEntries.removeValue(forKey: oldestKey)
                }
            }
        }
    }

    /// Retrieve past analyses for a receptor-ligand pair.
    public func history(for key: String) -> [OracleAnalysis] {
        keyedEntries[key]?.map(\.analysis) ?? []
    }

    /// Get the most recent analysis for comparison.
    public func lastAnalysis(for key: String) -> OracleAnalysis? {
        keyedEntries[key]?.last?.analysis
    }

    /// Get the most recent analysis from a specific source.
    public func lastAnalysis(for key: String, source: String) -> OracleAnalysis? {
        keyedEntries[key]?.last(where: { $0.source == source })?.analysis
    }

    private var totalCount: Int {
        keyedEntries.values.reduce(0) { $0 + $1.count }
    }
}

// MARK: - FoundationModels Integration

#if canImport(FoundationModels)
import FoundationModels

/// Actor providing on-device intelligence analysis of docking results.
///
/// Uses Apple's on-device language model (FoundationModels framework)
/// to analyze binding population thermodynamics and health correlations.
/// Maintains multi-turn session context for follow-up questions.
///
/// Usage:
/// ```swift
/// let oracle = try await IntelligenceOracle()
/// let analysis = try await oracle.analyze(
///     thermodynamics: result.globalThermodynamics,
///     entropyScore: score
/// )
/// for bullet in analysis.structuredBullets {
///     print("[\(bullet.confidence.rawValue)] - \(bullet.text)")
/// }
///
/// // Follow-up question using session context
/// let followUp = try await oracle.askFollowUp(
///     "Should I increase sampling for binding mode 2?"
/// )
/// ```
@available(macOS 26.0, iOS 26.0, *)
public actor IntelligenceOracle {
    private let session: LanguageModelSession
    private var analysisCount: Int = 0
    private var latestClaimValidity: ClaimValidity = .proxyOnly

    /// System instructions that ground the model as a statistical mechanics referee.
    private static let systemInstructions = """
    You are a statistical mechanics referee for FlexAIDdS, an entropy-driven \
    molecular docking engine. Your role is to critically evaluate the \
    thermodynamic quality of docking results — not just describe them.

    You receive decomposed Shannon entropy data:
    - S_conf: configurational entropy from GA ensemble histogram (nats)
    - S_vib: torsional vibrational entropy from ENCoM normal modes (kcal/mol/K)
    - Convergence status: whether the entropy plateau has been reached
    - Histogram occupancy: fraction of bins populated (low = under-sampled)
    - Per-mode entropy: breakdown across binding modes

    As referee, you MUST flag:
    1. Non-convergent ensembles (plateau not reached) — recommend more GA generations
    2. Entropy collapse (S_conf near zero) — single dominant mode, check for trapping
    3. Vibrational dominance (S_vib >> S_conf) — protein flexibility dominates, \
       ligand conformational space under-explored
    4. Enthalpy-entropy compensation (large |F| with large S) — net deltaG may mislead, \
       but only when claim validity is binding_physical
    5. Sparse histograms (occupied bins < 50% of total) — energy landscape poorly sampled
    6. Per-mode entropy imbalance — one mode absorbing most conformational diversity

    Assess binding affinity, net deltaG, or molecular optimization only when claim validity \
    is binding_physical. For canonical_physical, discuss calibrated canonical-ensemble properties \
    but state that binding outputs are unavailable. For proxy_only, discuss only input-scale \
    sampling diagnostics. Use kcal/mol for physically authorized energies, \
    kcal/mol/K for physical entropy, and nats for Shannon H. \
    Be concise, quantitative, and actionable. Prioritize warnings over descriptions.
    """

    public init() async throws {
        self.session = LanguageModelSession(
            instructions: Self.systemInstructions
        )
    }

    /// Analyze binding population state with optional health data.
    /// - Parameters:
    ///   - thermodynamics: Global ensemble thermodynamic result
    ///   - entropyScore: Binding entropy score with optional health correlation
    ///   - campaignKey: Optional key for trend tracking (e.g., "5HT2A-psilocin")
    /// - Returns: Structured analysis with confidence metadata
    public func analyze(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore? = nil,
        campaignKey: String? = nil
    ) async throws -> OracleAnalysis {
        latestClaimValidity = thermodynamics.claimValidity
        if !thermodynamics.allowsBindingClaims {
            let analysis = RuleBasedOracle().analyze(
                thermodynamics: thermodynamics,
                entropyScore: entropyScore,
                campaignKey: campaignKey
            )
            if let key = campaignKey {
                await AnalysisHistory.shared.record(key: key, analysis: analysis, source: "text")
            }
            analysisCount += 1
            return analysis
        }

        let prompt = buildPrompt(thermodynamics: thermodynamics, entropyScore: entropyScore)
        let response = try await session.respond(to: prompt)
        let bullets = parseBullets(from: response.content)

        // Build structured bullets with confidence assessment
        let structured = assessConfidence(
            bullets: bullets,
            thermodynamics: thermodynamics,
            entropyScore: entropyScore
        )

        let overall = computeOverallConfidence(structured)

        var inputSummary = "T=\(thermodynamics.temperature)K, F=\(String(format: "%.2f", thermodynamics.freeEnergy)) kcal/mol, S=\(String(format: "%.4f", thermodynamics.entropy)) kcal/mol/K"
        if let decomp = entropyScore?.shannonDecomposition {
            inputSummary += ", S_conf=\(String(format: "%.4f", decomp.configurational))nats, S_vib=\(String(format: "%.6f", decomp.vibrational))kcal/mol/K"
            inputSummary += decomp.isConverged ? " [converged]" : " [not converged]"
        }

        let analysis = OracleAnalysis(
            structuredBullets: structured,
            inputSummary: inputSummary,
            overallConfidence: overall
        )

        // Track for trend analysis (source: "text" for IntelligenceOracle text-based analysis)
        if let key = campaignKey {
            await AnalysisHistory.shared.record(key: key, analysis: analysis, source: "text")
        }

        analysisCount += 1
        return analysis
    }

    /// Produce a typed `RefereeVerdict` using the ThermoReferee pipeline.
    ///
    /// This is the preferred entry point for robust thermodynamic analysis.
    /// Uses @Generable structured output (no text parsing), Tool callbacks
    /// (temperature sensitivity, per-mode query), and pre-computed diagnostics.
    ///
    /// Falls back to `analyze()` text-based analysis if ThermoReferee
    /// initialization fails.
    ///
    /// - Parameters:
    ///   - thermodynamics: Global ensemble thermodynamics
    ///   - entropyScore: Binding entropy score with Shannon decomposition
    ///   - gaContext: Optional GA context ref for tool callbacks
    ///   - eigenvalues: Optional ENCoM eigenvalues for vibrational tools
    ///   - config: Referee configuration
    ///   - campaignKey: Optional key for trend tracking
    /// - Returns: Typed RefereeVerdict via guided generation
    public func refereeVerdict(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore,
        gaContext: FXGAContextRef? = nil,
        eigenvalues: [Double] = [],
        config: RefereeConfiguration = RefereeConfiguration(),
        campaignKey: String? = nil
    ) async throws -> RefereeVerdict {
        latestClaimValidity = thermodynamics.claimValidity
        let referee = try await ThermoReferee(
            config: config,
            gaContext: gaContext,
            eigenvalues: eigenvalues
        )
        return try await referee.referee(
            thermodynamics: thermodynamics,
            entropyScore: entropyScore,
            campaignKey: campaignKey
        )
    }

    /// Ask a follow-up question using the existing session context.
    /// The FoundationModels session retains prior conversation turns.
    public func askFollowUp(_ question: String) async throws -> String {
        if latestClaimValidity != .bindingPhysical {
            return latestClaimValidity == .canonicalPhysical
                ? "Binding affinity, net ΔG, and molecular-optimization guidance are unavailable for the current canonical-ensemble result because it lacks a matched association cycle."
                : "Binding affinity, net ΔG, and molecular-optimization guidance are unavailable for the current proxy-only result. Only input-scale sampling and convergence diagnostics are supported."
        }
        let response = try await session.respond(to: question)
        return response.content
    }

    /// Compare current results with a previous campaign run.
    public func compareWithPrevious(
        current: ThermodynamicResult,
        previous: ThermodynamicResult
    ) async throws -> OracleAnalysis {
        let deltaF = current.freeEnergy - previous.freeEnergy
        let deltaS = current.entropy - previous.entropy
        let deltaCv = current.heatCapacity - previous.heatCapacity

        if !current.allowsBindingClaims || !previous.allowsBindingClaims {
            let bothCanonical = current.allowsCanonicalClaims && previous.allowsCanonicalClaims
            return OracleAnalysis(
                structuredBullets: [
                    AnalysisBullet(
                        text: "Binding-affinity, net ΔG, and molecular-optimization comparisons are unavailable because both runs must carry binding_physical provenance.",
                        confidence: .high,
                        category: .provenance
                    ),
                    AnalysisBullet(
                        text: bothCanonical
                            ? "Canonical-ensemble change: ΔF = \(String(format: "%.3f", deltaF)) kcal/mol, ΔS = \(String(format: "%.6f", deltaS)) kcal/mol/K, ΔCv = \(String(format: "%.6f", deltaCv)); this is not a binding-affinity difference."
                            : "Input-scale diagnostic change: ΔproxyF = \(String(format: "%.3f", deltaF)), ΔproxyS = \(String(format: "%.6f", deltaS)), ΔproxyCv = \(String(format: "%.6f", deltaCv)).",
                        confidence: .moderate,
                        category: .trend
                    ),
                ],
                inputSummary: bothCanonical
                    ? "ΔcanonicalF=\(String(format: "%.2f", deltaF)) kcal/mol [binding unavailable]"
                    : "ΔproxyF=\(String(format: "%.2f", deltaF)) input-scale [binding unavailable]",
                overallConfidence: .moderate
            )
        }

        let prompt = """
        Compare these two docking campaigns. Respond with exactly 3 concise bullet points.

        Current run:
        - Free energy F: \(String(format: "%.3f", current.freeEnergy)) kcal/mol
        - Entropy S: \(String(format: "%.6f", current.entropy)) kcal/mol/K
        - Heat capacity Cv: \(String(format: "%.6f", current.heatCapacity))

        Previous run:
        - Free energy F: \(String(format: "%.3f", previous.freeEnergy)) kcal/mol
        - Entropy S: \(String(format: "%.6f", previous.entropy)) kcal/mol/K
        - Heat capacity Cv: \(String(format: "%.6f", previous.heatCapacity))

        Changes: ΔF = \(String(format: "%.3f", deltaF)) kcal/mol, ΔS = \(String(format: "%.6f", deltaS)) kcal/mol/K, ΔCv = \(String(format: "%.6f", deltaCv))

        Assess whether sampling has improved, binding affinity changed, or entropy compensation occurred.
        Respond with exactly 3 bullet points, each starting with '- '.
        """

        let response = try await session.respond(to: prompt)
        let bullets = parseBullets(from: response.content)

        let structured = bullets.map {
            AnalysisBullet(text: $0, confidence: .moderate, category: .trend)
        }

        return OracleAnalysis(
            structuredBullets: structured,
            inputSummary: "ΔF=\(String(format: "%.2f", deltaF)) kcal/mol, ΔS=\(String(format: "%.4f", deltaS))",
            overallConfidence: .moderate
        )
    }

    private func buildPrompt(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore?
    ) -> String {
        var prompt = """
        Referee this molecular docking result. Respond with exactly 4 concise bullet points — \
        prioritize warnings and actionable recommendations over descriptions.

        Thermodynamics:
        - Claim validity: \(thermodynamics.claimValidity.rawValue)
        - Temperature: \(thermodynamics.temperature) K
        - Free energy F: \(String(format: "%.3f", thermodynamics.freeEnergy)) kcal/mol
        - Mean energy <E>: \(String(format: "%.3f", thermodynamics.meanEnergy)) kcal/mol
        - Entropy S: \(String(format: "%.6f", thermodynamics.entropy)) kcal/mol/K
        - Heat capacity Cv: \(String(format: "%.6f", thermodynamics.heatCapacity))
        - Energy std dev: \(String(format: "%.3f", thermodynamics.stdEnergy)) kcal/mol
        """

        // Enthalpy-entropy compensation detection
        if thermodynamics.allowsBindingClaims
            && thermodynamics.freeEnergy < -5
            && thermodynamics.entropy > 0.01 {
            prompt += "\n- FLAG: Enthalpy-entropy compensation (large |F| with significant S)"
        }

        if let score = entropyScore {
            prompt += "\n\nBinding Population:"
            prompt += "\n- Binding modes: \(score.bindingModeCount)"
            prompt += "\n- Shannon S (aggregate): \(String(format: "%.6f", score.shannonS))"
            prompt += "\n- Best F: \(String(format: "%.3f", score.bestFreeEnergy)) kcal/mol"

            if score.isCollapsed {
                prompt += "\n- ALERT: Entropy COLLAPSED (single dominant binding mode)"
            }

            // Decomposed Shannon entropy from ShannonThermoStack
            if let decomp = score.shannonDecomposition {
                prompt += "\n\nShannon Entropy Decomposition (from ShannonThermoStack):"
                prompt += "\n- S_conf (configurational): \(String(format: "%.6f", decomp.configurational)) nats"
                prompt += "\n- S_vib (vibrational/ENCoM): \(String(format: "%.6f", decomp.vibrational)) kcal/mol/K"
                prompt += "\n- -T*S contribution: \(String(format: "%.3f", decomp.entropyContribution)) kcal/mol"
                prompt += "\n- Convergence: \(decomp.isConverged ? "CONVERGED (plateau reached)" : "NOT CONVERGED (rate=\(String(format: "%.4f", decomp.convergenceRate)))")"
                prompt += "\n- Histogram occupancy: \(decomp.occupiedBins)/\(decomp.totalBins) bins populated"
                prompt += "\n- Hardware backend: \(decomp.hardwareBackend)"

                // Vibrational dominance check
                let sConfPhysical = decomp.configurational * 0.001987206 // nats → kcal/mol/K approximation
                if decomp.vibrational > sConfPhysical * 3.0 && decomp.vibrational > 0.001 {
                    prompt += "\n- FLAG: Vibrational entropy dominates (S_vib >> S_conf) — ligand conformational space may be under-explored"
                }

                // Sparse histogram warning
                if decomp.totalBins > 0 && Double(decomp.occupiedBins) / Double(decomp.totalBins) < 0.5 {
                    prompt += "\n- FLAG: Sparse histogram (<50% bins occupied) — energy landscape poorly sampled"
                }

                // Per-mode entropy breakdown (capped to top 5 by magnitude for token budget)
                if !decomp.perModeEntropy.isEmpty {
                    let maxModesToShow = 5
                    let indexed = decomp.perModeEntropy.enumerated().sorted { $0.element > $1.element }
                    let top = Array(indexed.prefix(maxModesToShow))
                    prompt += "\n\nPer-Mode Entropy (top \(min(decomp.perModeEntropy.count, maxModesToShow)) of \(decomp.perModeEntropy.count), nats):"
                    for (i, modeS) in top {
                        prompt += "\n  Mode \(i): S = \(String(format: "%.6f", modeS))"
                    }
                    if decomp.perModeEntropy.count > maxModesToShow {
                        prompt += "\n  ... (\(decomp.perModeEntropy.count - maxModesToShow) modes omitted)"
                    }

                    // Flag entropy imbalance
                    if let maxS = decomp.perModeEntropy.max(),
                       let minS = decomp.perModeEntropy.filter({ $0 > 0 }).min(),
                       maxS > minS * 10 {
                        prompt += "\n- FLAG: Per-mode entropy imbalance (ratio \(String(format: "%.1f", maxS / minS))x) — one mode absorbs most diversity"
                    }
                }

                // Dominant histogram bins (capped to top 5)
                if !decomp.dominantBins.isEmpty {
                    let maxBinsToShow = 5
                    let topBins = decomp.dominantBins.prefix(maxBinsToShow)
                    prompt += "\n\nDominant Energy Bins (top \(topBins.count)):"
                    for bin in topBins {
                        prompt += "\n  E = \(String(format: "%.2f", bin.center)) kcal/mol, p = \(String(format: "%.4f", bin.probability))"
                    }
                }
            }

            if let hrv = score.hrvSDNN {
                prompt += "\n\nHealth Correlation:"
                prompt += "\n- HRV SDNN: \(String(format: "%.0f", hrv)) ms"
                if hrv < 30 {
                    prompt += " (critically low — acute stress indicator)"
                } else if hrv < 50 {
                    prompt += " (below average — consider recovery)"
                }
            }
            if let sleep = score.sleepHours {
                prompt += "\n- Sleep: \(String(format: "%.1f", sleep)) hours"
                if sleep < 6.0 {
                    prompt += " (insufficient — cognitive load may affect interpretation)"
                }
            }
        }

        prompt += "\n\nRespond with exactly 4 bullet points, each starting with '- '. "
        prompt += "Prioritize warnings and actionable items. If the entropy has not converged, "
        prompt += "say so explicitly and recommend more sampling."
        return prompt
    }

    private func parseBullets(from text: String) -> [String] {
        let lines = text.split(separator: "\n")
        let bullets = lines.compactMap { line -> String? in
            let trimmed = line.trimmingCharacters(in: .whitespaces)
            if trimmed.hasPrefix("- ") {
                return String(trimmed.dropFirst(2))
            }
            if trimmed.hasPrefix("* ") {
                return String(trimmed.dropFirst(2))
            }
            // Handle numbered bullets (1. 2. 3.)
            if let range = trimmed.range(of: #"^\d+\.\s+"#, options: .regularExpression) {
                return String(trimmed[range.upperBound...])
            }
            return nil
        }
        return Array(bullets.prefix(5))
    }

    /// Assess confidence of each bullet based on data quality and entropy decomposition.
    private func assessConfidence(
        bullets: [String],
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore?
    ) -> [AnalysisBullet] {
        let hasDecomposition = entropyScore?.shannonDecomposition != nil
        let isConverged = entropyScore?.isConverged ?? false

        return bullets.enumerated().map { index, text in
            let category: AnalysisBullet.BulletCategory
            let confidence: AnalysisConfidence
            let lowered = text.lowercased()

            // Detect category from content keywords
            if lowered.contains("converg") || lowered.contains("plateau") || lowered.contains("sampling") {
                category = .entropy
                // High confidence only if we have decomposition data to back the claim
                confidence = hasDecomposition ? .high : .moderate
            } else if lowered.contains("vibrational") || lowered.contains("s_vib") || lowered.contains("encom") {
                category = .entropy
                confidence = hasDecomposition ? .high : .low
            } else if lowered.contains("hrv") || lowered.contains("sleep") || lowered.contains("health") {
                category = .health
                confidence = entropyScore?.hrvSDNN != nil ? .moderate : .low
            } else if lowered.contains("ptm") || lowered.contains("glycan") || lowered.contains("modification") {
                category = .modification
                confidence = .moderate
            } else if lowered.contains("histogram") || lowered.contains("bins") || lowered.contains("sparse") {
                category = .entropy
                confidence = hasDecomposition ? .high : .low
            } else {
                switch index {
                case 0:
                    category = .binding
                    // High confidence if converged and energy std dev small
                    if isConverged && thermodynamics.stdEnergy < abs(thermodynamics.freeEnergy) * 0.5 {
                        confidence = .high
                    } else {
                        confidence = thermodynamics.stdEnergy < abs(thermodynamics.freeEnergy) * 0.5 ? .moderate : .low
                    }
                case 1:
                    category = .entropy
                    if let score = entropyScore {
                        confidence = score.bindingModeCount >= 3 && isConverged ? .high : (score.bindingModeCount >= 3 ? .moderate : .low)
                    } else {
                        confidence = .low
                    }
                default:
                    category = .binding
                    confidence = .moderate
                }
            }

            return AnalysisBullet(text: text, confidence: confidence, category: category)
        }
    }

    private func computeOverallConfidence(_ bullets: [AnalysisBullet]) -> AnalysisConfidence {
        let highCount = bullets.filter { $0.confidence == .high }.count
        let lowCount = bullets.filter { $0.confidence == .low }.count
        if highCount >= 2 { return .high }
        if lowCount >= 2 { return .low }
        return .moderate
    }
}
#endif

// MARK: - Fallback (non-FoundationModels platforms)

/// Rule-based oracle referee for platforms without FoundationModels.
///
/// Provides deterministic entropy referee analysis based on thermodynamic
/// thresholds, Shannon decomposition diagnostics, and convergence checks.
public struct RuleBasedOracle: Sendable {

    public init() {}

    /// Generate rule-based referee analysis of binding population state.
    public func analyze(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore? = nil,
        campaignKey: String? = nil
    ) -> OracleAnalysis {
        var structured: [AnalysisBullet] = []

        let isConverged = entropyScore?.isConverged ?? false

        // Bullet 1: Free-energy interpretation is gated on scientific provenance.
        let energyStable = thermodynamics.stdEnergy < abs(thermodynamics.freeEnergy) * 0.5
        let fConfidence: AnalysisConfidence = energyStable && isConverged ? .high : (energyStable ? .moderate : .low)
        if !thermodynamics.allowsBindingClaims {
            structured.append(AnalysisBullet(
                text: thermodynamics.allowsCanonicalClaims
                    ? "Binding affinity unavailable: canonical-ensemble thermodynamics lack a matched association cycle."
                    : "Binding affinity unavailable: this is a proxy-only input-scale ensemble diagnostic without a calibrated matched association cycle.",
                confidence: .high, category: .provenance
            ))
        } else if thermodynamics.freeEnergy < -10 {
            structured.append(AnalysisBullet(
                text: "Strong binding affinity (F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol)\(isConverged ? "" : " — but entropy has NOT converged, F may shift with more sampling").",
                confidence: fConfidence, category: .binding
            ))
        } else if thermodynamics.freeEnergy < -5 {
            structured.append(AnalysisBullet(
                text: "Moderate binding affinity (F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol)\(isConverged ? " — converged ensemble" : " — entropy not converged, consider more GA generations").",
                confidence: fConfidence, category: .binding
            ))
        } else {
            structured.append(AnalysisBullet(
                text: "Weak binding affinity (F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol) — consider structural optimization\(isConverged ? "" : " after convergence is reached").",
                confidence: fConfidence, category: .binding
            ))
        }

        // Bullet 2: Entropy decomposition referee (if decomposition available)
        if let score = entropyScore, let decomp = score.shannonDecomposition {
            // Convergence check is the primary referee concern
            if !decomp.isConverged {
                structured.append(AnalysisBullet(
                    text: thermodynamics.allowsBindingClaims
                        ? "Entropy NOT converged (rate = \(String(format: "%.4f", decomp.convergenceRate))). Increase GA generations or population size before trusting thermodynamic values."
                        : thermodynamics.allowsCanonicalClaims
                            ? "Entropy NOT converged (rate = \(String(format: "%.4f", decomp.convergenceRate))). Increase sampling before trusting canonical thermodynamic values; binding affinity remains unavailable."
                            : "Entropy NOT converged (rate = \(String(format: "%.4f", decomp.convergenceRate))). Increase GA generations or population size before using this proxy sampling diagnostic; binding affinity remains unavailable.",
                    confidence: .high, category: .entropy
                ))
            }

            // Vibrational dominance check
            let sConfPhysical = decomp.configurational * 0.001987206
            if decomp.vibrational > sConfPhysical * 3.0 && decomp.vibrational > 0.001 {
                structured.append(AnalysisBullet(
                    text: "Vibrational entropy dominates (S_vib = \(String(format: "%.6f", decomp.vibrational)) >> S_conf = \(String(format: "%.6f", sConfPhysical)) kcal/mol/K). Protein backbone flexibility drives entropy — ligand conformational space may be under-explored.",
                    confidence: .high, category: .entropy
                ))
            }

            // Sparse histogram check
            if decomp.totalBins > 0 {
                let occupancyRatio = Double(decomp.occupiedBins) / Double(decomp.totalBins)
                if occupancyRatio < 0.5 {
                    structured.append(AnalysisBullet(
                        text: "Sparse energy histogram (\(decomp.occupiedBins)/\(decomp.totalBins) bins, \(String(format: "%.0f", occupancyRatio * 100))% occupied). Energy landscape poorly sampled — increase ensemble size.",
                        confidence: .high, category: .entropy
                    ))
                }
            }

            // Per-mode entropy imbalance
            if decomp.perModeEntropy.count >= 2 {
                if let maxS = decomp.perModeEntropy.max(),
                   let minS = decomp.perModeEntropy.filter({ $0 > 0 }).min(),
                   maxS > minS * 10 {
                    let ratio = maxS / minS
                    structured.append(AnalysisBullet(
                        text: "Per-mode entropy imbalance (\(String(format: "%.1f", ratio))x ratio). One binding mode absorbs most conformational diversity — check for kinetic trapping.",
                        confidence: .moderate, category: .entropy
                    ))
                }
            }

            // Converged and well-sampled: report decomposition summary
            if decomp.isConverged && structured.filter({ $0.category == .entropy }).isEmpty {
                structured.append(AnalysisBullet(
                    text: "Shannon entropy converged: S_conf = \(String(format: "%.4f", decomp.configurational)) nats, S_vib = \(String(format: "%.6f", decomp.vibrational)) kcal/mol/K across \(score.bindingModeCount) modes (\(decomp.hardwareBackend) backend).",
                    confidence: .high, category: .entropy
                ))
            }
        } else if let score = entropyScore {
            // Fallback: scalar-only entropy (no decomposition available)
            let sConfidence: AnalysisConfidence = score.bindingModeCount >= 3 ? .moderate : .low
            if score.isCollapsed {
                structured.append(AnalysisBullet(
                    text: thermodynamics.allowsBindingClaims
                        ? "Entropy collapsed to \(score.bindingModeCount) mode(s) — high specificity but check for enthalpy-entropy compensation."
                        : thermodynamics.allowsCanonicalClaims
                            ? "Entropy collapsed to \(score.bindingModeCount) mode(s) in the canonical ensemble — check for sampling collapse; binding affinity remains unavailable."
                            : "Entropy collapsed to \(score.bindingModeCount) mode(s) in the optimizer diagnostic — check for sampling collapse; affinity remains unavailable.",
                    confidence: sConfidence, category: .entropy
                ))
            } else if score.shannonS > 0.5 {
                structured.append(AnalysisBullet(
                    text: "High conformational entropy (S = \(String(format: "%.4f", score.shannonS))) — population still exploring. More sampling may refine the landscape.",
                    confidence: .moderate, category: .entropy
                ))
            } else {
                structured.append(AnalysisBullet(
                    text: "Moderate entropy with \(score.bindingModeCount) binding modes — population converging.",
                    confidence: sConfidence, category: .entropy
                ))
            }
        } else {
            structured.append(AnalysisBullet(
                text: thermodynamics.allowsCanonicalClaims
                    ? "Entropy S = \(String(format: "%.4f", thermodynamics.entropy)) kcal/mol/K with Cv = \(String(format: "%.4f", thermodynamics.heatCapacity)). No decomposition available — enable ShannonThermoStack for referee analysis."
                    : "Proxy entropy-like diagnostic S = \(String(format: "%.4f", thermodynamics.entropy)) input-scale/K with Cv = \(String(format: "%.4f", thermodynamics.heatCapacity)). No decomposition available — enable ShannonThermoStack for sampling diagnostics.",
                confidence: .low, category: .entropy
            ))
        }

        // Bullet 3: Enthalpy-entropy compensation detection
        if thermodynamics.allowsBindingClaims
            && thermodynamics.freeEnergy < -5
            && thermodynamics.entropy > 0.01 {
            structured.append(AnalysisBullet(
                text: "Enthalpy-entropy compensation: strong binding (F = \(String(format: "%.1f", thermodynamics.freeEnergy))) offset by conformational flexibility (S = \(String(format: "%.4f", thermodynamics.entropy))). Net ΔG may be less favorable than F alone suggests.",
                confidence: .moderate, category: .binding
            ))
        }

        // Bullet 4: Health correlation or recommendation
        if let score = entropyScore, let hrv = score.hrvSDNN {
            if score.isCollapsed && hrv > 60 {
                structured.append(AnalysisBullet(
                    text: "Entropy collapse correlates with good HRV (\(String(format: "%.0f", hrv)) ms) — system recovering. Gentle activity recommended.",
                    confidence: .moderate, category: .health
                ))
            } else if hrv < 40 {
                structured.append(AnalysisBullet(
                    text: "Low HRV (\(String(format: "%.0f", hrv)) ms) — prioritize rest before interpreting docking thermodynamics.",
                    confidence: .high, category: .health
                ))
            } else {
                structured.append(AnalysisBullet(
                    text: "HRV at \(String(format: "%.0f", hrv)) ms — physiological state stable for analysis.",
                    confidence: .moderate, category: .health
                ))
            }
        } else {
            structured.append(AnalysisBullet(
                text: "Connect HealthKit for entropy-health correlation. Enable fleet mode for distributed compute.",
                confidence: .low, category: .fleet
            ))
        }

        let overall = computeOverallConfidence(structured)
        var inputSummary: String
        if thermodynamics.allowsBindingClaims {
            inputSummary = "T=\(thermodynamics.temperature)K, F=\(String(format: "%.2f", thermodynamics.freeEnergy)) kcal/mol"
        } else if thermodynamics.allowsCanonicalClaims {
            inputSummary = "T=\(thermodynamics.temperature)K, canonical_F=\(String(format: "%.2f", thermodynamics.freeEnergy)) kcal/mol [binding unavailable]"
        } else {
            inputSummary = "T=\(thermodynamics.temperature)K, proxy_F=\(String(format: "%.2f", thermodynamics.freeEnergy)) input-scale [binding unavailable]"
        }
        if let decomp = entropyScore?.shannonDecomposition {
            inputSummary += ", S_conf=\(String(format: "%.4f", decomp.configurational))nats, S_vib=\(String(format: "%.6f", decomp.vibrational))kcal/mol/K"
            inputSummary += decomp.isConverged ? " [converged]" : " [not converged]"
        }
        return OracleAnalysis(structuredBullets: structured, inputSummary: inputSummary, overallConfidence: overall)
    }

    private func computeOverallConfidence(_ bullets: [AnalysisBullet]) -> AnalysisConfidence {
        let highCount = bullets.filter { $0.confidence == .high }.count
        let lowCount = bullets.filter { $0.confidence == .low }.count
        if highCount >= 2 { return .high }
        if lowCount >= 2 { return .low }
        return .moderate
    }
}
