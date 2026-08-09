// ThermoRefereeTests.swift — Unit tests for RuleBasedReferee deterministic referee
//
// Tests the rule-based (non-LLM) referee logic with identical thresholds
// to the FoundationModels ThermoReferee and TypeScript RuleBasedReferee.
// Verifies cross-platform parity of findings, severity, and trustworthiness.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import XCTest
@testable import FlexAIDdS
@testable import HealthIntegration
@testable import Intelligence

final class ThermoRefereeTests: XCTestCase {

    // MARK: - Helpers

    private func makeThermo(
        freeEnergy: Double = -8.0, entropy: Double = 0.005,
        temperature: Double = 298.15, stdEnergy: Double = 1.0,
        heatCapacity: Double = 0.1, meanEnergy: Double = -7.0,
        scientificProvenance: ScientificProvenance? = nil
    ) -> ThermodynamicResult {
        ThermodynamicResult(
            temperature: temperature, logZ: 10.0,
            freeEnergy: freeEnergy, meanEnergy: meanEnergy,
            meanEnergySq: meanEnergy * meanEnergy + stdEnergy * stdEnergy,
            heatCapacity: heatCapacity, entropy: entropy,
            stdEnergy: stdEnergy,
            scientificProvenance: scientificProvenance
        )
    }

    private var bindingPhysicalProvenance: ScientificProvenance {
        ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .matchedAssociationCycle,
            energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
            measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4",
            referenceProvenance: "sha256:01d26f66e709d388d7b971de6204680694e7cca10b1570506a5738c6909a6442"
        )
    }

    private var canonicalPhysicalProvenance: ScientificProvenance {
        ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .boundOnly,
            energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
            measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4"
        )
    }

    private func makeScore(
        decomposition: ShannonEntropyDecomposition? = nil,
        bindingModeCount: Int = 5, shannonS: Double = 1.0,
        bestFreeEnergy: Double = -10.0
    ) -> BindingEntropyScore {
        var score = BindingEntropyScore(
            shannonS: shannonS, temperature: 298.15,
            bindingModeCount: bindingModeCount,
            bestFreeEnergy: bestFreeEnergy,
            heatCapacity: 0.1
        )
        score.shannonDecomposition = decomposition
        return score
    }

    private func makeDecomp(
        configurational: Double = 2.5, vibrational: Double = 0.0005,
        entropyContribution: Double = -1.5, isConverged: Bool = true,
        convergenceRate: Double = 0.001, occupiedBins: Int = 18,
        totalBins: Int = 20, perModeEntropy: [Double] = [0.5, 0.4, 0.3, 0.2, 0.1]
    ) -> ShannonEntropyDecomposition {
        ShannonEntropyDecomposition(
            configurational: configurational, vibrational: vibrational,
            entropyContribution: entropyContribution, isConverged: isConverged,
            convergenceRate: convergenceRate, hardwareBackend: "scalar",
            occupiedBins: occupiedBins, totalBins: totalBins,
            perModeEntropy: perModeEntropy
        )
    }

    // MARK: - Convergence Tests

    func testRuleBasedReferee_nonConverged_flagsCritical() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(isConverged: false, convergenceRate: 0.05)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo()

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        // Must flag critical convergence issue
        let convergenceFinding = verdict.findings.first { $0.category == "convergence" }
        XCTAssertNotNil(convergenceFinding)
        XCTAssertEqual(convergenceFinding?.severity, "critical")
        XCTAssertFalse(verdict.overallTrustworthy, "Non-converged should be untrustworthy")
        XCTAssertLessThanOrEqual(verdict.confidence, 0.5)
    }

    func testRuleBasedReferee_convergedGoodHistogram_passesTrust() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(
            isConverged: true, occupiedBins: 18, totalBins: 20,
            perModeEntropy: [0.5, 0.4, 0.3, 0.2, 0.1]
        )
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo(freeEnergy: -12.0)

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        XCTAssertTrue(verdict.overallTrustworthy, "Converged with good histogram should be trustworthy")
        let convergenceFinding = verdict.findings.first { $0.category == "convergence" }
        XCTAssertEqual(convergenceFinding?.severity, "pass")
        XCTAssertGreaterThanOrEqual(verdict.confidence, 0.8)
    }

    // MARK: - Histogram Tests

    func testRuleBasedReferee_sparseHistogram_flagsWarning() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(occupiedBins: 8, totalBins: 20)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo()

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let histFinding = verdict.findings.first { $0.category == "histogram" }
        XCTAssertNotNil(histFinding, "Should flag sparse histogram")
        XCTAssertEqual(histFinding?.severity, "warning")
    }

    func testRuleBasedReferee_criticallySparseHistogram_flagsCritical() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(occupiedBins: 4, totalBins: 20)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo()

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let histFinding = verdict.findings.first { $0.category == "histogram" }
        XCTAssertNotNil(histFinding, "Should flag critically sparse histogram")
        XCTAssertEqual(histFinding?.severity, "critical")
        XCTAssertFalse(verdict.overallTrustworthy)
    }

    // MARK: - Entropy Balance Tests

    func testRuleBasedReferee_vibrationalDominance_flagsWarning() {
        let referee = RuleBasedReferee()
        // S_vib = 0.01 >> S_conf_physical = 0.5 * 0.001987 ≈ 0.001 → ratio ~10x
        let decomp = makeDecomp(configurational: 0.5, vibrational: 0.01)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo()

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let balanceFinding = verdict.findings.first { $0.category == "entropyBalance" }
        XCTAssertNotNil(balanceFinding, "Should flag vibrational dominance")
        XCTAssertEqual(balanceFinding?.severity, "warning")
    }

    // MARK: - Mode Imbalance Tests

    func testRuleBasedReferee_modeImbalance_detects10xRatio() {
        let referee = RuleBasedReferee()
        // Mode entropies with >10x ratio
        let decomp = makeDecomp(perModeEntropy: [5.0, 0.3, 0.2])
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo()

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let modeFinding = verdict.findings.first { $0.category == "modeBalance" }
        XCTAssertNotNil(modeFinding, "Should detect mode imbalance >10x")
        XCTAssertEqual(modeFinding?.severity, "warning")
        XCTAssertTrue(modeFinding?.detail.contains("25.0x") ?? false, "Should report 25x ratio")
    }

    // MARK: - Compensation Tests

    func testRuleBasedReferee_enthalpyEntropyCompensation() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp()
        let score = makeScore(decomposition: decomp)
        // F < -5 and S > 0.01 triggers compensation
        let thermo = makeThermo(
            freeEnergy: -12.0,
            entropy: 0.02,
            scientificProvenance: bindingPhysicalProvenance
        )

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let compFinding = verdict.findings.first { $0.category == "compensation" }
        XCTAssertNotNil(compFinding, "Should detect enthalpy-entropy compensation")
        XCTAssertEqual(compFinding?.severity, "advisory")
    }

    func testRuleBasedReferee_proxySuppressesCompensationAndNetDeltaG() {
        let referee = RuleBasedReferee()
        let score = makeScore(decomposition: makeDecomp())
        let thermo = makeThermo(freeEnergy: -12.0, entropy: 0.02)

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)
        let output = (verdict.findings.map(\.detail) + [verdict.recommendedAction])
            .joined(separator: " ")
            .lowercased()

        XCTAssertNil(verdict.findings.first { $0.category == "compensation" })
        XCTAssertFalse(output.contains("net δg"))
        XCTAssertFalse(output.contains("strong binding"))
        XCTAssertTrue(output.contains("proxy diagnostic"))
        XCTAssertTrue(output.contains("unavailable"))
    }

    // MARK: - No Decomposition Tests

    func testRuleBasedReferee_noDecomposition_advisoryWarning() {
        let referee = RuleBasedReferee()
        let score = makeScore(decomposition: nil)
        let thermo = makeThermo()

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let convergenceFinding = verdict.findings.first { $0.category == "convergence" }
        XCTAssertNotNil(convergenceFinding)
        XCTAssertEqual(convergenceFinding?.severity, "warning")
        XCTAssertTrue(convergenceFinding?.detail.contains("ShannonThermoStack") ?? false)
    }

    // MARK: - Affinity Tests

    func testRuleBasedReferee_proxySuppressesAffinityClaim() {
        let referee = RuleBasedReferee()
        let score = makeScore(decomposition: makeDecomp(isConverged: true))
        let thermo = makeThermo(freeEnergy: -15.0)

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        XCTAssertNil(verdict.findings.first { $0.category == "affinity" })
        let provenanceFinding = verdict.findings.first { $0.category == "provenance" }
        XCTAssertNotNil(provenanceFinding)
        XCTAssertTrue(provenanceFinding?.title.contains("unavailable") ?? false)
    }

    func testRuleBasedReferee_canonicalIsPhysicalButBindingUnavailable() {
        let referee = RuleBasedReferee()
        let score = makeScore(decomposition: makeDecomp(isConverged: true))
        let thermo = makeThermo(
            freeEnergy: -15.0,
            scientificProvenance: canonicalPhysicalProvenance
        )

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)
        let output = (verdict.findings.map(\.detail) + [verdict.recommendedAction])
            .joined(separator: " ")
            .lowercased()

        XCTAssertEqual(thermo.claimValidity, .canonicalPhysical)
        XCTAssertNil(verdict.findings.first { $0.category == "affinity" })
        XCTAssertTrue(output.contains("canonical"))
        XCTAssertTrue(output.contains("binding affinity") && output.contains("unavailable"))
        XCTAssertFalse(output.contains("proxy diagnostic"))
        XCTAssertFalse(output.contains("input-scale"))
    }

    func testRuleBasedReferee_strongAffinityRequiresMatchedCycle() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(isConverged: true)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo(
            freeEnergy: -15.0,
            scientificProvenance: bindingPhysicalProvenance
        )

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let affinityFinding = verdict.findings.first { $0.category == "affinity" }
        XCTAssertNotNil(affinityFinding)
        XCTAssertEqual(affinityFinding?.severity, "pass")
        XCTAssertTrue(affinityFinding?.title.contains("Strong") ?? false)
    }

    func testRuleBasedReferee_weakAffinity() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(isConverged: true)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo(
            freeEnergy: -2.0,
            scientificProvenance: bindingPhysicalProvenance
        )

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        let affinityFinding = verdict.findings.first { $0.category == "affinity" }
        XCTAssertNotNil(affinityFinding)
        XCTAssertEqual(affinityFinding?.severity, "warning")
        XCTAssertTrue(affinityFinding?.title.contains("Weak") ?? false)
    }

    // MARK: - Cross-Platform Verdict Codable

    func testCrossPlatformVerdictCodable() throws {
        let original = CrossPlatformRefereeVerdict(
            findings: [
                CrossPlatformRefereeFinding(
                    title: "Entropy converged",
                    detail: "Plateau reached on Metal backend.",
                    severity: "pass",
                    category: "convergence"
                ),
                CrossPlatformRefereeFinding(
                    title: "Sparse histogram",
                    detail: "8/20 bins occupied (40%).",
                    severity: "warning",
                    category: "histogram"
                )
            ],
            overallTrustworthy: true,
            recommendedAction: "Proceed with lead optimization.",
            confidence: 0.85
        )

        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformRefereeVerdict.self, from: data)

        XCTAssertEqual(original, decoded)
        XCTAssertEqual(decoded.findings.count, 2)
        XCTAssertEqual(decoded.findings[0].severity, "pass")
        XCTAssertEqual(decoded.findings[1].category, "histogram")
        XCTAssertEqual(decoded.overallTrustworthy, true)
        XCTAssertEqual(decoded.confidence, 0.85, accuracy: 1e-10)
    }

    // MARK: - PreComputedThermoContext Validation

    func testPreComputedContext_vibrationalDominanceRatio() {
        // Verify the pre-computation math for vibrational dominance
        let kB = 0.001987206
        let sConf: Double = 2.0  // nats
        let sVib: Double = 0.01  // kcal/mol/K
        let sConfPhysical = sConf * kB  // ≈ 0.003974

        let ratio = sVib / sConfPhysical
        XCTAssertGreaterThan(ratio, 2.5, "S_vib should be ~2.5x S_conf for these values")
        XCTAssertLessThan(ratio, 2.6)
    }

    func testPreComputedContext_histogramOccupancy() {
        let occupiedBins = 12
        let totalBins = 20
        let occupancy = Double(occupiedBins) / Double(totalBins)
        XCTAssertEqual(occupancy, 0.6, accuracy: 1e-10)
        XCTAssertGreaterThan(occupancy, 0.5, "12/20 should pass the 50% threshold")
    }

    // MARK: - Recommendation Logic

    func testRuleBasedReferee_untrustedRecommendation() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(isConverged: false)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo(scientificProvenance: bindingPhysicalProvenance)

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        XCTAssertTrue(verdict.recommendedAction.contains("Do not trust"))
    }

    func testRuleBasedReferee_proxyRecommendationIsDiagnosticOnly() {
        let referee = RuleBasedReferee()
        let decomp = makeDecomp(isConverged: true)
        let score = makeScore(decomposition: decomp)
        let thermo = makeThermo(freeEnergy: -8.0)

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        XCTAssertTrue(verdict.recommendedAction.contains("Proxy diagnostic only"))
        XCTAssertTrue(verdict.recommendedAction.contains("unavailable"))
        XCTAssertFalse(verdict.recommendedAction.contains("optimization of the ligand"))
        XCTAssertFalse(verdict.recommendedAction.contains("lead optimization"))
    }

    func testRuleBasedReferee_bindingPhysicalPermitsLeadOptimization() {
        let referee = RuleBasedReferee()
        let score = makeScore(decomposition: makeDecomp(isConverged: true))
        let thermo = makeThermo(
            freeEnergy: -8.0,
            scientificProvenance: bindingPhysicalProvenance
        )

        let verdict = referee.referee(thermodynamics: thermo, entropyScore: score)

        XCTAssertTrue(verdict.recommendedAction.contains("lead optimization"))
    }

    func testRuleBasedOracle_proxyNeverEmitsAffinityNetDeltaGOrOptimization() {
        let oracle = RuleBasedOracle()
        // shannonS must be below BindingEntropyScore.isCollapsed's 0.1
        // threshold for the collapse branch under test to be reachable.
        let score = makeScore(decomposition: nil, bindingModeCount: 1, shannonS: 0.05)
        let thermo = makeThermo(freeEnergy: -15.0, entropy: 0.02)

        let analysis = oracle.analyze(thermodynamics: thermo, entropyScore: score)
        let output = analysis.bullets.joined(separator: " ").lowercased()

        XCTAssertTrue(output.contains("binding affinity unavailable"))
        XCTAssertTrue(output.contains("sampling collapse"))
        XCTAssertFalse(output.contains("strong binding"))
        XCTAssertFalse(output.contains("net δg"))
        XCTAssertFalse(output.contains("structural optimization"))
        XCTAssertTrue(analysis.inputSummary.contains("proxy_F"))
        XCTAssertTrue(analysis.inputSummary.contains("input-scale"))
        XCTAssertFalse(analysis.inputSummary.contains("kcal/mol"))
    }

    func testRuleBasedOracle_bindingPhysicalPermitsAffinityAndCompensation() {
        let oracle = RuleBasedOracle()
        let thermo = makeThermo(
            freeEnergy: -15.0,
            entropy: 0.02,
            scientificProvenance: bindingPhysicalProvenance
        )

        let analysis = oracle.analyze(thermodynamics: thermo)
        let output = analysis.bullets.joined(separator: " ").lowercased()

        XCTAssertTrue(output.contains("strong binding affinity"))
        XCTAssertTrue(output.contains("net δg"))
        XCTAssertTrue(analysis.inputSummary.contains("kcal/mol"))
    }

    func testRuleBasedOracle_canonicalUsesPhysicalUnitsWithoutBindingClaim() {
        let oracle = RuleBasedOracle()
        let thermo = makeThermo(
            freeEnergy: -15.0,
            entropy: 0.02,
            scientificProvenance: canonicalPhysicalProvenance
        )

        let analysis = oracle.analyze(thermodynamics: thermo)
        let output = analysis.bullets.joined(separator: " ").lowercased()

        XCTAssertTrue(output.contains("canonical"))
        XCTAssertTrue(output.contains("binding affinity unavailable"))
        XCTAssertFalse(output.contains("strong binding"))
        XCTAssertFalse(output.contains("net δg"))
        XCTAssertTrue(analysis.inputSummary.contains("canonical_F"))
        XCTAssertTrue(analysis.inputSummary.contains("kcal/mol"))
        XCTAssertFalse(analysis.inputSummary.contains("proxy_F"))
    }
}
