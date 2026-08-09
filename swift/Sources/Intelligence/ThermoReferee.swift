// ThermoReferee.swift — Orchestrator for on-device thermodynamic referee analysis
//
// ThermoReferee is the main entry point for robust thermodynamic analysis
// using Apple FoundationModels. It:
//   1. Pre-computes all thermodynamic data on CPU/GPU (never asks the LLM to do math)
//   2. Builds a concise token-budget-aware prompt (<2K tokens)
//   3. Registers Tool callbacks so the model can request additional computations
//   4. Uses @Generable for typed RefereeVerdict output (no text parsing)
//   5. Falls back to RuleBasedReferee on non-Apple platforms
//
// The on-device ~3B model is used for LANGUAGE UNDERSTANDING (interpreting
// the thermodynamic landscape), not for COMPUTATION (all math stays in C++/Metal).
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS
import FlexAIDCore
import HealthIntegration

// MARK: - Referee Configuration

/// Configuration for the ThermoReferee analysis pipeline.
public struct RefereeConfiguration: Sendable {
    /// Maximum number of findings to request (3-6)
    public var maxFindings: Int = 5

    /// Whether to enable tool callbacks (temperature sensitivity, per-mode query)
    public var enableTools: Bool = true

    /// Temperature range for sensitivity analysis (K)
    public var sensitivityTempLow: Double = 280.0
    public var sensitivityTempHigh: Double = 320.0

    /// Campaign key for trend tracking
    public var campaignKey: String?

    public init() {}
}

// MARK: - Pre-Computed Thermodynamic Context

/// All thermodynamic data pre-computed on CPU/GPU before prompting the model.
/// The LLM never does arithmetic — it only interprets pre-digested values.
public struct PreComputedThermoContext: Sendable, Codable {
    /// Derived scientific validity; never supplied or overridden by an LLM.
    public let claimValidity: ClaimValidity

    // Global ensemble
    public let temperature: Double
    public let freeEnergy: Double
    public let entropy: Double
    public let heatCapacity: Double
    public let meanEnergy: Double
    public let stdEnergy: Double

    // Shannon decomposition
    public let sConf: Double            // nats
    public let sVib: Double             // kcal/mol/K
    public let entropyContribution: Double  // -T*S kcal/mol
    public let isConverged: Bool
    public let convergenceRate: Double
    public let occupiedBins: Int
    public let totalBins: Int
    public let hardwareBackend: String

    // Pre-computed diagnostics (CPU-side, not LLM-computed)
    public let sConfPhysical: Double    // S_conf in kcal/mol/K (converted from nats)
    public let vibrationalDominance: Double  // S_vib / S_conf ratio (inf if S_conf ~ 0)
    public let histogramOccupancy: Double    // occupiedBins / totalBins
    public let modeCount: Int
    public let perModeEntropy: [Double]
    public let modeEntropyImbalance: Double  // max/min ratio (inf if only one mode)

    // Binding assessment
    public let bestModeFreeEnergy: Double
    public let hasEnthalpyEntropyCompensation: Bool

    // Health correlation (optional)
    public let hrvSDNN: Double?
    public let sleepHours: Double?
}

// MARK: - ThermoReferee Actor (FoundationModels)

#if canImport(FoundationModels)
import FoundationModels

@available(macOS 26.0, iOS 26.0, *)
public actor ThermoReferee {
    private var session: LanguageModelSession?
    private let config: RefereeConfiguration
    private var gaContext: FXGAContextRef?
    private var eigenvalues: [Double]
    private var latestClaimValidity: ClaimValidity = .proxyOnly

    /// Referee system instructions — concise to fit within 4K token budget.
    /// All computation results are pre-digested; the model only interprets.
    private static let refereeInstructions = """
        You are a statistical mechanics referee for FlexAIDdS molecular docking. \
        All thermodynamic values are pre-computed — DO NOT perform arithmetic. \
        Evaluate the quality of the docking result and produce a structured verdict.

        Key criteria (in priority order):
        1. CONVERGENCE: Has the Shannon entropy reached a plateau? If not, results are unreliable.
        2. HISTOGRAM QUALITY: Are >50% of energy bins occupied? Sparse = under-sampled.
        3. ENTROPY BALANCE: Is S_vib >> S_conf? If so, backbone flexibility dominates and \
           ligand conformational space is under-explored.
        4. MODE BALANCE: Does one binding mode absorb most entropy? May indicate kinetic trapping.
        5. COMPENSATION: Only for binding_physical results, large |F| with large S may suggest \
           enthalpy-entropy compensation and a less favorable net ΔG.
        6. CLAIM DOMAIN: Assess physical binding affinity only when claimValidity is binding_physical. \
           Otherwise describe F as a proxy ensemble diagnostic and state that affinity is unavailable.

        Be concise. Each finding should have a clear severity (pass/advisory/warning/critical) \
        and actionable recommendation. Set overallTrustworthy=false if convergence or \
        histogram quality fails.
        """

    public init(config: RefereeConfiguration = RefereeConfiguration(),
                gaContext: FXGAContextRef? = nil,
                eigenvalues: [Double] = []) async throws {
        self.config = config
        self.gaContext = gaContext
        self.eigenvalues = eigenvalues

        // Create session with tools if GA context available
        if config.enableTools, let ctx = gaContext {
            var tools: [any Tool] = [
                HelmholtzCalculatorTool()
            ]
            tools.append(RecomputeAtTemperatureTool(gaContext: ctx))
            tools.append(PerModeEntropyTool(gaContext: ctx))
            if !eigenvalues.isEmpty {
                tools.append(VibrationalEntropyTool(eigenvalues: eigenvalues))
            }
            self.session = LanguageModelSession(
                tools: tools,
                instructions: Self.refereeInstructions
            )
        } else {
            self.session = LanguageModelSession(
                instructions: Self.refereeInstructions
            )
        }
    }

    /// Run the full referee pipeline: pre-compute → prompt → structured verdict.
    ///
    /// - Parameters:
    ///   - thermodynamics: Global ensemble thermodynamics
    ///   - entropyScore: Binding entropy score with Shannon decomposition
    ///   - campaignKey: Optional key for trend tracking
    /// - Returns: Typed `RefereeVerdict` via @Generable guided generation
    public func referee(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore,
        campaignKey: String? = nil
    ) async throws -> RefereeVerdict {
        latestClaimValidity = thermodynamics.claimValidity
        if !thermodynamics.allowsBindingClaims {
            let isCanonical = thermodynamics.allowsCanonicalClaims
            let verdict = RefereeVerdict(
                findings: [RefereeFinding(
                    title: "Binding affinity unavailable",
                    detail: isCanonical
                        ? "This canonical_physical result supports calibrated canonical-ensemble interpretation, but no matched association cycle is available."
                        : "This proxy_only result is an input-scale sampling diagnostic; no calibrated matched association cycle is available.",
                    severity: .advisory,
                    category: .recommendation
                )],
                overallTrustworthy: false,
                recommendedAction: isCanonical
                    ? "Use canonical-ensemble and convergence diagnostics only. Binding-affinity and molecular-optimization guidance are unavailable."
                    : "Use input-scale convergence and sampling diagnostics only. Binding-affinity and molecular-optimization guidance are unavailable.",
                confidence: 1.0
            )
            if let key = campaignKey ?? config.campaignKey {
                let analysis = OracleAnalysis(
                    structuredBullets: [AnalysisBullet(
                        text: isCanonical
                            ? "Binding affinity unavailable: canonical-ensemble result without a matched association cycle."
                            : "Binding affinity unavailable: proxy-only sampling diagnostic.",
                        confidence: .high,
                        category: .provenance
                    )],
                    inputSummary: "claim_validity=\(thermodynamics.claimValidity.rawValue) [binding unavailable]",
                    overallConfidence: .high
                )
                await AnalysisHistory.shared.record(key: key, analysis: analysis, source: "referee")
            }
            return verdict
        }

        guard let session else {
            throw RefereeError.sessionUnavailable
        }

        // Step 1: Pre-compute all diagnostics on CPU (no LLM math)
        let context = preCompute(thermodynamics: thermodynamics, entropyScore: entropyScore)

        // Step 2: Build token-efficient prompt (may throw if budget exceeded)
        let prompt = try buildRefereePrompt(context: context)

        // Step 3: Request structured verdict via @Generable guided generation
        let verdict = try await session.respond(
            to: prompt,
            generating: RefereeVerdict.self
        ).content

        // Step 4: Record for trend tracking
        if let key = campaignKey ?? config.campaignKey {
            let oracleBullets = verdict.findings.map { finding in
                AnalysisBullet(
                    text: "[\(finding.severity.rawValue.uppercased())] \(finding.title): \(finding.detail)",
                    confidence: severityToConfidence(finding.severity),
                    category: categoryToBulletCategory(finding.category)
                )
            }
            let analysis = OracleAnalysis(
                structuredBullets: oracleBullets,
                inputSummary: context.summaryString,
                overallConfidence: verdict.overallTrustworthy ? .high : .low
            )
            await AnalysisHistory.shared.record(key: key, analysis: analysis, source: "referee")
        }

        return verdict
    }

    /// Ask a follow-up question about the last verdict.
    /// The FoundationModels session retains full context including tool call history.
    public func followUp(_ question: String) async throws -> String {
        if latestClaimValidity != .bindingPhysical {
            return latestClaimValidity == .canonicalPhysical
                ? "Binding affinity, net ΔG, and molecular-optimization guidance are unavailable for the current canonical-ensemble result because it lacks a matched association cycle."
                : "Binding affinity, net ΔG, and molecular-optimization guidance are unavailable for the current proxy-only result. Only input-scale sampling and convergence diagnostics are supported."
        }
        guard let session else { throw RefereeError.sessionUnavailable }
        let response = try await session.respond(to: question)
        return response.content
    }

    /// Run temperature sensitivity analysis.
    /// Recomputes at low and high temps, then asks the model to assess sensitivity.
    public func temperatureSensitivity(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore
    ) async throws -> TemperatureSensitivity {
        latestClaimValidity = thermodynamics.claimValidity
        if !thermodynamics.allowsBindingClaims {
            return TemperatureSensitivity(
                isSensitive: false,
                assessment: thermodynamics.allowsCanonicalClaims
                    ? "Binding temperature sensitivity is unavailable: canonical thermodynamics lack a matched association cycle."
                    : "Binding temperature sensitivity is unavailable for a proxy-only sampling diagnostic.",
                recommendedTempRange: "Unavailable without calibrated matched-association provenance."
            )
        }

        guard let session else { throw RefereeError.sessionUnavailable }

        let context = preCompute(thermodynamics: thermodynamics, entropyScore: entropyScore)
        let prompt = """
            Assess temperature sensitivity of these docking results. \
            Current: T=\(context.temperature)K, F=\(String(format: "%.2f", context.freeEnergy)), \
            S_conf=\(String(format: "%.4f", context.sConf))nats, converged=\(context.isConverged). \
            Use the recompute_at_temperature tool at \(config.sensitivityTempLow)K and \(config.sensitivityTempHigh)K \
            to check how results change.
            """

        return try await session.respond(
            to: prompt,
            generating: TemperatureSensitivity.self
        ).content
    }

    /// Compare two docking campaigns and produce a structured verdict.
    public func compare(
        current: ThermodynamicResult, currentScore: BindingEntropyScore,
        previous: ThermodynamicResult, previousScore: BindingEntropyScore
    ) async throws -> ComparativeVerdict {
        latestClaimValidity = current.allowsBindingClaims && previous.allowsBindingClaims
            ? .bindingPhysical
            : .proxyOnly
        if !current.allowsBindingClaims || !previous.allowsBindingClaims {
            return ComparativeVerdict(
                bindingTrend: "unavailable",
                convergenceTrend: "Sampling diagnostics may be compared separately.",
                summary: "Binding-affinity and molecular-optimization comparisons require binding_physical provenance for both runs.",
                improved: false
            )
        }

        guard let session else { throw RefereeError.sessionUnavailable }

        let curCtx = preCompute(thermodynamics: current, entropyScore: currentScore)
        let prevCtx = preCompute(thermodynamics: previous, entropyScore: previousScore)

        let prompt = """
            Compare two docking campaigns:
            CURRENT: F=\(String(format: "%.2f", curCtx.freeEnergy)), S_conf=\(String(format: "%.4f", curCtx.sConf))nats, \
            converged=\(curCtx.isConverged), \(curCtx.modeCount) modes, occupancy=\(String(format: "%.0f", curCtx.histogramOccupancy * 100))%
            PREVIOUS: F=\(String(format: "%.2f", prevCtx.freeEnergy)), S_conf=\(String(format: "%.4f", prevCtx.sConf))nats, \
            converged=\(prevCtx.isConverged), \(prevCtx.modeCount) modes, occupancy=\(String(format: "%.0f", prevCtx.histogramOccupancy * 100))%
            ΔF = \(String(format: "%.3f", curCtx.freeEnergy - prevCtx.freeEnergy)) kcal/mol
            """

        return try await session.respond(
            to: prompt,
            generating: ComparativeVerdict.self
        ).content
    }

    // MARK: - Pre-Computation (CPU/GPU, no LLM)

    private func preCompute(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore
    ) -> PreComputedThermoContext {
        let decomp = entropyScore.shannonDecomposition
        let sConf = decomp?.configurational ?? 0
        let sVib = decomp?.vibrational ?? 0
        let kB = 0.001987206
        let sConfPhysical = sConf * kB

        let vibrationalDominance: Double = sConfPhysical > 1e-10 ? sVib / sConfPhysical : 0.0
        let occupiedBins = decomp?.occupiedBins ?? 0
        let totalBins = decomp?.totalBins ?? 1
        let histogramOccupancy = Double(occupiedBins) / Double(max(totalBins, 1))

        let perModeEntropy = decomp?.perModeEntropy ?? []
        let modeEntropyImbalance: Double = {
            guard perModeEntropy.count >= 2,
                  let maxS = perModeEntropy.max(),
                  let minS = perModeEntropy.filter({ $0 > 0 }).min(),
                  minS > 0 else { return 1.0 }
            return maxS / minS
        }()

        return PreComputedThermoContext(
            claimValidity: thermodynamics.claimValidity,
            temperature: thermodynamics.temperature,
            freeEnergy: thermodynamics.freeEnergy,
            entropy: thermodynamics.entropy,
            heatCapacity: thermodynamics.heatCapacity,
            meanEnergy: thermodynamics.meanEnergy,
            stdEnergy: thermodynamics.stdEnergy,
            sConf: sConf,
            sVib: sVib,
            entropyContribution: decomp?.entropyContribution ?? 0,
            isConverged: decomp?.isConverged ?? false,
            convergenceRate: decomp?.convergenceRate ?? 1.0,
            occupiedBins: occupiedBins,
            totalBins: totalBins,
            hardwareBackend: decomp?.hardwareBackend ?? "unknown",
            sConfPhysical: sConfPhysical,
            vibrationalDominance: vibrationalDominance,
            histogramOccupancy: histogramOccupancy,
            modeCount: entropyScore.bindingModeCount,
            perModeEntropy: perModeEntropy,
            modeEntropyImbalance: modeEntropyImbalance,
            bestModeFreeEnergy: entropyScore.bestFreeEnergy,
            hasEnthalpyEntropyCompensation: thermodynamics.allowsBindingClaims
                && thermodynamics.freeEnergy < -5
                && thermodynamics.entropy > 0.01,
            hrvSDNN: entropyScore.hrvSDNN,
            sleepHours: entropyScore.sleepHours
        )
    }

    /// Estimate token count for a string (word count × 1.3 heuristic for English + numbers).
    private func estimateTokenCount(_ text: String) -> Int {
        let wordCount = text.split(whereSeparator: { $0.isWhitespace || $0.isNewline }).count
        return Int(ceil(Double(wordCount) * 1.3))
    }

    /// Build a compact prompt that fits within the 4K token budget.
    /// Pre-computed flags are stated as facts; the model interprets, not computes.
    /// Per-mode entropy is capped to top 5 modes by magnitude to stay within budget.
    private func buildRefereePrompt(context: PreComputedThermoContext) throws -> String {
        var p = """
            Referee this docking result. Produce a RefereeVerdict with up to \(config.maxFindings) findings.

            Claim validity: \(context.claimValidity.rawValue)
            T=\(context.temperature)K, F=\(String(format: "%.2f", context.freeEnergy))kcal/mol, \
            <E>=\(String(format: "%.2f", context.meanEnergy))kcal/mol, σ=\(String(format: "%.2f", context.stdEnergy))kcal/mol
            S_conf=\(String(format: "%.4f", context.sConf))nats (\(String(format: "%.6f", context.sConfPhysical))kcal/mol/K), \
            S_vib=\(String(format: "%.6f", context.sVib))kcal/mol/K, -TS=\(String(format: "%.2f", context.entropyContribution))kcal/mol
            """

        // State pre-computed flags as facts
        p += "\nConverged: \(context.isConverged ? "YES" : "NO (rate=\(String(format: "%.4f", context.convergenceRate)))")"
        p += "\nHistogram: \(context.occupiedBins)/\(context.totalBins) bins (\(String(format: "%.0f", context.histogramOccupancy * 100))%)"
        p += "\nModes: \(context.modeCount), best F=\(String(format: "%.2f", context.bestModeFreeEnergy))kcal/mol"

        if context.vibrationalDominance.isFinite && context.vibrationalDominance > 3.0 {
            p += "\nFLAG: S_vib >> S_conf (ratio=\(String(format: "%.1f", context.vibrationalDominance))x)"
        }
        if context.histogramOccupancy < 0.5 {
            p += "\nFLAG: Sparse histogram (<50% occupied)"
        }
        if context.modeEntropyImbalance > 10.0 {
            p += "\nFLAG: Mode entropy imbalance (\(String(format: "%.1f", context.modeEntropyImbalance))x)"
        }
        if context.hasEnthalpyEntropyCompensation {
            p += "\nFLAG: Enthalpy-entropy compensation (F=\(String(format: "%.1f", context.freeEnergy)), S=\(String(format: "%.4f", context.entropy)))"
        }

        // Per-mode entropy: cap to top 5 modes sorted by entropy magnitude
        if !context.perModeEntropy.isEmpty {
            let maxModesToShow = 5
            let indexed = context.perModeEntropy.enumerated().sorted { $0.element > $1.element }
            let top = indexed.prefix(maxModesToShow)
            p += "\nTop modes by S:"
            for (i, s) in top {
                p += " M\(i)=\(String(format: "%.4f", s))"
            }
            if context.perModeEntropy.count > maxModesToShow {
                p += " ...(\(context.perModeEntropy.count - maxModesToShow) omitted)"
            }
        }

        if let hrv = context.hrvSDNN {
            p += "\nHRV: \(String(format: "%.0f", hrv))ms"
        }

        p += "\nBackend: \(context.hardwareBackend)"

        // Token budget guard: system instructions ~150 tokens + prompt must stay < 3800
        let estimated = estimateTokenCount(p)
        if estimated > 3000 {
            // Emergency truncation: remove per-mode line
            if let modeRange = p.range(of: "\nTop modes by S:.*", options: .regularExpression) {
                p.removeSubrange(modeRange)
            }
        }
        if estimateTokenCount(p) > 3800 {
            throw RefereeError.contextWindowExceeded
        }

        return p
    }

    // MARK: - Helpers

    private func severityToConfidence(_ severity: RefereeSeverity) -> AnalysisConfidence {
        switch severity {
        case .pass: return .high
        case .advisory: return .moderate
        case .warning: return .moderate
        case .critical: return .high
        }
    }

    private func categoryToBulletCategory(_ category: RefereeCategory) -> AnalysisBullet.BulletCategory {
        switch category {
        case .convergence, .entropyBalance, .histogram, .modeBalance:
            return .entropy
        case .compensation, .affinity:
            return .binding
        case .recommendation:
            return .fleet
        }
    }
}

// MARK: - Errors

public enum RefereeError: Error, LocalizedError {
    case sessionUnavailable
    case contextWindowExceeded

    public var errorDescription: String? {
        switch self {
        case .sessionUnavailable:
            return "FoundationModels session unavailable — use RuleBasedReferee instead"
        case .contextWindowExceeded:
            return "Prompt exceeded 4K token budget — reduce input data"
        }
    }
}

#endif

// MARK: - Rule-Based Referee (Cross-Platform Fallback)

/// Deterministic referee for platforms without FoundationModels.
/// Produces the same `CrossPlatformRefereeVerdict` structure using threshold logic.
public struct RuleBasedReferee: Sendable {

    public init() {}

    /// Run deterministic referee analysis.
    public func referee(
        thermodynamics: ThermodynamicResult,
        entropyScore: BindingEntropyScore
    ) -> CrossPlatformRefereeVerdict {
        var findings: [CrossPlatformRefereeFinding] = []
        var trustworthy = true

        let decomp = entropyScore.shannonDecomposition
        let kB = 0.001987206

        // 1. Convergence check (highest priority)
        if let d = decomp {
            if !d.isConverged {
                trustworthy = false
                let detail = thermodynamics.allowsBindingClaims
                    ? "Shannon entropy has not reached a plateau (rate = \(String(format: "%.4f", d.convergenceRate))). Increase GA generations before trusting F or S values."
                    : thermodynamics.allowsCanonicalClaims
                        ? "Shannon entropy has not reached a plateau (rate = \(String(format: "%.4f", d.convergenceRate))). Increase sampling before trusting canonical F or S; binding affinity remains unavailable."
                        : "Shannon entropy has not reached a plateau (rate = \(String(format: "%.4f", d.convergenceRate))). Increase GA generations before using this input-scale sampling diagnostic; binding affinity remains unavailable."
                findings.append(CrossPlatformRefereeFinding(
                    title: "Entropy not converged",
                    detail: detail,
                    severity: "critical",
                    category: "convergence"
                ))
            } else {
                let detail = thermodynamics.allowsBindingClaims
                    ? "Shannon entropy plateau reached on \(d.hardwareBackend) backend. Physical thermodynamic interpretation is supported by the declared binding provenance."
                    : thermodynamics.allowsCanonicalClaims
                        ? "Shannon entropy plateau reached on \(d.hardwareBackend) backend. Canonical thermodynamic interpretation is supported, but binding affinity remains unavailable."
                        : "Shannon entropy plateau reached on \(d.hardwareBackend) backend. The optimizer-sampling diagnostic is converged, but binding affinity remains unavailable."
                findings.append(CrossPlatformRefereeFinding(
                    title: "Entropy converged",
                    detail: detail,
                    severity: "pass",
                    category: "convergence"
                ))
            }
        } else {
            findings.append(CrossPlatformRefereeFinding(
                title: "No Shannon decomposition",
                detail: "ShannonThermoStack data unavailable. Cannot assess convergence — enable ShannonThermoStack for reliable referee analysis.",
                severity: "warning",
                category: "convergence"
            ))
        }

        // 2. Histogram quality
        if let d = decomp, d.totalBins > 0 {
            let occupancy = Double(d.occupiedBins) / Double(d.totalBins)
            if occupancy < 0.3 {
                trustworthy = false
                findings.append(CrossPlatformRefereeFinding(
                    title: "Critically sparse histogram",
                    detail: "Only \(d.occupiedBins)/\(d.totalBins) bins occupied (\(String(format: "%.0f", occupancy * 100))%). Energy landscape severely under-sampled. Double the population size.",
                    severity: "critical",
                    category: "histogram"
                ))
            } else if occupancy < 0.5 {
                findings.append(CrossPlatformRefereeFinding(
                    title: "Sparse histogram",
                    detail: "\(d.occupiedBins)/\(d.totalBins) bins occupied (\(String(format: "%.0f", occupancy * 100))%). Consider increasing ensemble size for better coverage.",
                    severity: "warning",
                    category: "histogram"
                ))
            }
        }

        // 3. Entropy balance (S_vib vs S_conf)
        if let d = decomp {
            let sConfPhysical = d.configurational * kB
            if d.vibrational > sConfPhysical * 3.0 && d.vibrational > 0.001 {
                findings.append(CrossPlatformRefereeFinding(
                    title: "Vibrational entropy dominates",
                    detail: "S_vib (\(String(format: "%.6f", d.vibrational))) >> S_conf (\(String(format: "%.6f", sConfPhysical))) kcal/mol/K. Protein backbone flexibility drives entropy; ligand conformational space may be under-explored.",
                    severity: "warning",
                    category: "entropyBalance"
                ))
            }
        }

        // 4. Mode entropy imbalance
        if let d = decomp, d.perModeEntropy.count >= 2 {
            if let maxS = d.perModeEntropy.max(),
               let minS = d.perModeEntropy.filter({ $0 > 0 }).min(),
               maxS > minS * 10 {
                let ratio = maxS / minS
                findings.append(CrossPlatformRefereeFinding(
                    title: "Mode entropy imbalance",
                    detail: "Entropy ratio across modes is \(String(format: "%.1f", ratio))x. One binding mode absorbs most conformational diversity — check for kinetic trapping.",
                    severity: "warning",
                    category: "modeBalance"
                ))
            }
        }

        // 5. Enthalpy-entropy compensation
        if thermodynamics.allowsBindingClaims
            && thermodynamics.freeEnergy < -5
            && thermodynamics.entropy > 0.01 {
            findings.append(CrossPlatformRefereeFinding(
                title: "Enthalpy-entropy compensation",
                detail: "Strong binding (F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol) offset by conformational flexibility (S = \(String(format: "%.4f", thermodynamics.entropy))). Net ΔG may be less favorable than F alone suggests.",
                severity: "advisory",
                category: "compensation"
            ))
        }

        // 6. Binding claims are forbidden unless a matched, calibrated
        // association cycle is present. Legacy/CF results fail closed.
        let converged = decomp?.isConverged ?? false
        if !thermodynamics.allowsBindingClaims {
            findings.append(CrossPlatformRefereeFinding(
                title: "Binding affinity unavailable",
                detail: thermodynamics.allowsCanonicalClaims
                    ? "The result supports canonical-ensemble thermodynamics, but no matched association cycle was supplied."
                    : "The result is a proxy-only input-scale ensemble diagnostic; no calibrated matched association cycle was supplied.",
                severity: "advisory",
                category: "provenance"
            ))
        } else if thermodynamics.freeEnergy < -10 {
            findings.append(CrossPlatformRefereeFinding(
                title: "Strong binding affinity",
                detail: "F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol\(converged ? " (converged)" : " — convergence not confirmed").",
                severity: converged ? "pass" : "advisory",
                category: "affinity"
            ))
        } else if thermodynamics.freeEnergy < -5 {
            findings.append(CrossPlatformRefereeFinding(
                title: "Moderate binding affinity",
                detail: "F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol — reasonable drug candidate\(converged ? "" : ", pending convergence").",
                severity: "advisory",
                category: "affinity"
            ))
        } else {
            findings.append(CrossPlatformRefereeFinding(
                title: "Weak binding affinity",
                detail: "F = \(String(format: "%.1f", thermodynamics.freeEnergy)) kcal/mol — consider structural optimization.",
                severity: "warning",
                category: "affinity"
            ))
        }

        // Recommendation
        let action: String
        if !thermodynamics.allowsBindingClaims {
            if !trustworthy {
                action = thermodynamics.allowsCanonicalClaims
                    ? "Improve sampling before interpreting canonical thermodynamics. Binding affinity and molecular optimization are unavailable."
                    : "Proxy diagnostic only: improve sampling before interpreting the optimizer landscape. Binding affinity and molecular optimization are unavailable."
            } else {
                action = thermodynamics.allowsCanonicalClaims
                    ? "Canonical thermodynamics only. Binding affinity and molecular optimization are unavailable without a matched association cycle."
                    : "Proxy diagnostic only. Binding affinity and molecular optimization are unavailable without a calibrated matched association cycle."
            }
        } else if !trustworthy {
            action = "Do not trust current thermodynamic values. Increase GA generations and population size, then re-run."
        } else if decomp == nil {
            action = "Enable ShannonThermoStack for decomposed entropy analysis."
        } else if thermodynamics.freeEnergy > -5 {
            action = "Consider structural optimization of the ligand or alternative binding sites."
        } else {
            action = "Results are reliable. Proceed with lead optimization."
        }

        let confidence: Double
        if !trustworthy { confidence = 0.3 }
        else if decomp?.isConverged == true { confidence = 0.9 }
        else { confidence = 0.6 }

        return CrossPlatformRefereeVerdict(
            findings: findings,
            overallTrustworthy: trustworthy,
            recommendedAction: action,
            confidence: confidence
        )
    }
}

// MARK: - DockingRunner + Referee Convenience

extension DockingRunner {
    /// Run the GA, extract Shannon decomposition, and produce a referee verdict in one call.
    ///
    /// This is the simplest entry point for docking + quality assessment:
    /// ```swift
    /// let runner = DockingRunner(configPath: "config.inp", gaPath: "ga.inp")
    /// let (result, verdict) = try await runner.runWithReferee()
    /// if !verdict.overallTrustworthy {
    ///     print("Warning: \(verdict.recommendedAction)")
    /// }
    /// ```
    public func runWithReferee(
        config: RefereeConfiguration = RefereeConfiguration()
    ) async throws -> (DockingResult, CrossPlatformRefereeVerdict) {
        let (result, decomposition) = try await runWithDecomposition()

        var score = BindingEntropyScore(
            shannonS: result.globalThermodynamics.entropy,
            temperature: result.temperature,
            bindingModeCount: result.bindingModes.count,
            bestFreeEnergy: result.bindingModes.first?.freeEnergy ?? 0,
            heatCapacity: result.globalThermodynamics.heatCapacity
        )
        score.shannonDecomposition = decomposition

        let referee = RuleBasedReferee()
        let verdict = referee.referee(
            thermodynamics: result.globalThermodynamics,
            entropyScore: score
        )

        return (result, verdict)
    }
}

// MARK: - Convenience Extension on PreComputedThermoContext

extension PreComputedThermoContext {
    /// Compact summary string for OracleAnalysis.inputSummary.
    var summaryString: String {
        var s: String
        if claimValidity == .bindingPhysical {
            s = "T=\(temperature)K, F=\(String(format: "%.2f", freeEnergy))kcal/mol"
        } else {
            s = "T=\(temperature)K, proxy_F=\(String(format: "%.2f", freeEnergy)) input-scale"
        }
        s += ", S_conf=\(String(format: "%.4f", sConf))nats, S_vib=\(String(format: "%.6f", sVib))kcal/mol/K"
        s += isConverged ? " [converged]" : " [not converged]"
        return s
    }
}
