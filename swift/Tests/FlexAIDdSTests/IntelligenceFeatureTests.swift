// IntelligenceFeatureTests.swift — Unit tests for all Intelligence RuleBased* fallbacks
//
// Tests deterministic rule-based logic for BindingModeNarrator, CleftAssessor,
// ConvergenceCoach, FleetExplainer, HealthEntropyInsight, VibrationalInterpreter,
// SelectivityAnalyst, CampaignJournalist, and LigandFitCritic.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import XCTest
@testable import FlexAIDdS
@testable import Intelligence

private let intelligenceBindingProvenance = ScientificProvenance(
    energyDomain: .calibratedKcalPerMol,
    ensembleMeasure: .enumeratedMicrostates,
    referenceState: .matchedAssociationCycle,
    energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
    measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4",
    referenceProvenance: "sha256:01d26f66e709d388d7b971de6204680694e7cca10b1570506a5738c6909a6442"
)

// MARK: - BindingModeNarrator Tests

final class BindingModeNarratorTests: XCTestCase {

    private func makeMode(
        index: Int = 0, poseCount: Int = 10, freeEnergy: Double = -8.0,
        entropy: Double = 0.005, enthalpy: Double = -6.0, heatCapacity: Double = 0.1,
        boltzmannWeight: Double = 0.6, temperature: Double = 298.15
    ) -> ModeProfile {
        ModeProfile(
            index: index, poseCount: poseCount, freeEnergy: freeEnergy,
            entropy: entropy, enthalpy: enthalpy, heatCapacity: heatCapacity,
            boltzmannWeight: boltzmannWeight, temperature: temperature
        )
    }

    private func makeContext(
        modes: [ModeProfile],
        temperature: Double = 298.15,
        scientificProvenance: ScientificProvenance? = intelligenceBindingProvenance
    ) -> PreComputedModeContext {
        let entropies = modes.map(\.entropy)
        let imbalance: Double = {
            guard entropies.count >= 2,
                  let maxS = entropies.max(),
                  let minS = entropies.filter({ $0 > 0 }).min(),
                  minS > 0 else { return 1.0 }
            return maxS / minS
        }()
        return PreComputedModeContext(
            temperature: temperature,
            totalModes: modes.count,
            totalPoses: modes.reduce(0) { $0 + $1.poseCount },
            globalFreeEnergy: modes.first?.freeEnergy ?? 0,
            modes: modes,
            entropyImbalance: imbalance,
            dominantModeIndex: 0,
            scientificProvenance: scientificProvenance
        )
    }

    func testSingleMode() {
        let narrator = RuleBasedModeNarrator()
        let context = makeContext(modes: [makeMode()])
        let result = narrator.narrate(context: context)

        XCTAssertEqual(result.modeDescriptions.count, 1)
        XCTAssertEqual(result.confidence, 0.5, "Single mode should have 0.5 confidence")
    }

    func testEntropyDrivenMode() {
        let narrator = RuleBasedModeNarrator()
        // |−TS| = 298.15 * 0.1 = 29.8 >> |ΔH| = 2.0 → entropy-driven
        let mode = makeMode(entropy: 0.1, enthalpy: -2.0)
        let context = makeContext(modes: [mode])
        let result = narrator.narrate(context: context)

        XCTAssertTrue(result.modeDescriptions[0].characterization.contains("entropy-driven"))
        XCTAssertTrue(result.modeDescriptions[0].optimizationHint.contains("Rigidify"))
    }

    func testEnthalpyDrivenMode() {
        let narrator = RuleBasedModeNarrator()
        // |ΔH| = 20.0 >> |−TS| = 298.15 * 0.001 = 0.3 → enthalpy-driven
        let mode = makeMode(entropy: 0.001, enthalpy: -20.0)
        let context = makeContext(modes: [mode])
        let result = narrator.narrate(context: context)

        XCTAssertTrue(result.modeDescriptions[0].characterization.contains("enthalpy-driven"))
        XCTAssertTrue(result.modeDescriptions[0].optimizationHint.contains("flexible"))
    }

    func testThreeModesRanking() {
        let narrator = RuleBasedModeNarrator()
        let modes = [
            makeMode(index: 0, boltzmannWeight: 0.6),
            makeMode(index: 1, boltzmannWeight: 0.3),
            makeMode(index: 2, boltzmannWeight: 0.1)
        ]
        let context = makeContext(modes: modes)
        let result = narrator.narrate(context: context)

        XCTAssertEqual(result.modeDescriptions.count, 3)
        XCTAssertEqual(result.confidence, 0.8, "Multiple modes should have 0.8 confidence")
        XCTAssertTrue(result.selectivityInsight.contains("Mode 1"))
    }

    func testEmptyModes() {
        let narrator = RuleBasedModeNarrator()
        let context = makeContext(modes: [])
        let result = narrator.narrate(context: context)

        XCTAssertTrue(result.modeDescriptions.isEmpty)
        XCTAssertTrue(result.selectivityInsight.contains("No dominant mode"))
    }

    func testProxyModeNarrativeSuppressesPhysicalAndOptimizationClaims() {
        let narrator = RuleBasedModeNarrator()
        let context = makeContext(
            modes: [makeMode(entropy: 0.1, enthalpy: -20)],
            scientificProvenance: nil
        )

        let result = narrator.narrate(context: context)
        let output = (result.modeDescriptions.flatMap { [$0.characterization, $0.optimizationHint] }
            + [result.selectivityInsight]).joined(separator: " ").lowercased()

        XCTAssertTrue(output.contains("optimizer cluster diagnostic"))
        XCTAssertTrue(output.contains("unavailable"))
        XCTAssertFalse(output.contains("kcal/mol"))
        XCTAssertFalse(output.contains("affinity needs boosting"))
        XCTAssertFalse(output.contains("focus sar optimization"))
    }
}

// MARK: - CleftAssessor Tests

final class CleftAssessorTests: XCTestCase {

    // Druggability tiers and ligand-property prescriptions are drug-design
    // claims, so the tiering logic is only reachable with binding_physical
    // evidence. Tests that exercise the tiers supply it explicitly.
    private func makeCleft(
        volume: Double = 500, depth: Double = 8.0, sphereCount: Int = 50,
        maxSphereRadius: Double = 3.0, hydrophobicFraction: Double = 0.6,
        anchorResidueCount: Int = 6, elongation: Double = 0.4, solventExposure: Double = 0.2,
        provenance: ScientificProvenance? = intelligenceBindingProvenance
    ) -> CleftFeatures {
        CleftFeatures(
            volume: volume, depth: depth, sphereCount: sphereCount,
            maxSphereRadius: maxSphereRadius, hydrophobicFraction: hydrophobicFraction,
            anchorResidueCount: anchorResidueCount, elongation: elongation,
            solventExposure: solventExposure,
            scientificProvenance: provenance
        )
    }

    func testDrugDesignClaimsUnavailableWithoutProvenance() {
        let assessor = RuleBasedCleftAssessor()
        // Default (no provenance) must fail closed: no tier, no prescription.
        let cleft = makeCleft(volume: 600, depth: 10.0, hydrophobicFraction: 0.6,
                              anchorResidueCount: 8, solventExposure: 0.15,
                              provenance: nil)

        XCTAssertEqual(cleft.claimValidity, .proxyOnly)
        XCTAssertFalse(cleft.allowsDrugDesignClaims)

        let result = assessor.assess(cleft: cleft)
        XCTAssertEqual(result.druggability, .unavailable)
        XCTAssertTrue(result.suggestedLigandProperties.contains("unavailable"))
        XCTAssertFalse(result.summary.contains("Druggability: high"))
    }

    func testProxyProvenanceStillFailsClosed() {
        let assessor = RuleBasedCleftAssessor()
        // CF-optimizer provenance is structurally proxy-only.
        let cleft = makeCleft(provenance: .proxyContactFunction)

        XCTAssertEqual(cleft.claimValidity, .proxyOnly)
        XCTAssertEqual(assessor.assess(cleft: cleft).druggability, .unavailable)
    }

    func testHighDruggability() {
        let assessor = RuleBasedCleftAssessor()
        let cleft = makeCleft(volume: 600, depth: 10.0, hydrophobicFraction: 0.6,
                              anchorResidueCount: 8, solventExposure: 0.15)
        let result = assessor.assess(cleft: cleft)

        XCTAssertEqual(result.druggability, .high)
        XCTAssertTrue(result.warnings.isEmpty)
    }

    func testUndruggableSmallVolume() {
        let assessor = RuleBasedCleftAssessor()
        let cleft = makeCleft(volume: 80, depth: 2.0, hydrophobicFraction: 0.2,
                              anchorResidueCount: 1, solventExposure: 0.7)
        let result = assessor.assess(cleft: cleft)

        XCTAssertTrue(result.druggability == .low || result.druggability == .undruggable)
        XCTAssertFalse(result.warnings.isEmpty, "Should have warnings about small volume")
    }

    func testHighSolventExposure() {
        let assessor = RuleBasedCleftAssessor()
        let cleft = makeCleft(solventExposure: 0.75)
        let result = assessor.assess(cleft: cleft)

        let hasExposureWarning = result.warnings.contains { $0.contains("solvent") || $0.contains("exposure") }
        XCTAssertTrue(hasExposureWarning)
    }

    func testHeavilyHydrophobic() {
        let assessor = RuleBasedCleftAssessor()
        let cleft = makeCleft(hydrophobicFraction: 0.9)
        let result = assessor.assess(cleft: cleft)

        let hasHydrophobicWarning = result.warnings.contains { $0.lowercased().contains("hydrophobic") }
        XCTAssertTrue(hasHydrophobicWarning)
        XCTAssertTrue(result.suggestedLigandProperties.contains("Lipophilic"))
    }

    func testEdgeCaseVolume() {
        let assessor = RuleBasedCleftAssessor()
        // Volume between 1000 and 1200 — borderline
        let cleft = makeCleft(volume: 1100)
        let result = assessor.assess(cleft: cleft)

        // Should still be assessed (moderate or high depending on other features)
        XCTAssertNotEqual(result.druggability, .undruggable)
    }
}

// MARK: - ConvergenceCoach Tests

final class ConvergenceCoachTests: XCTestCase {

    private func makeSnapshot(
        currentGeneration: Int = 500, maxGenerations: Int = 1000,
        bestFitness: Double = -8.0, meanFitness: Double = -5.0,
        populationDiversity: Double = 2.5, generationsSinceImprovement: Int = 50,
        fitnessTrajectory: [Double] = [-7.0, -7.5, -8.0],
        diversityTrajectory: [Double] = [3.0, 2.8, 2.5],
        isImproving: Bool = true, isDiversityCollapsed: Bool = false,
        populationSize: Int = 300
    ) -> GAProgressSnapshot {
        GAProgressSnapshot(
            currentGeneration: currentGeneration, maxGenerations: maxGenerations,
            bestFitness: bestFitness, meanFitness: meanFitness,
            populationDiversity: populationDiversity,
            generationsSinceImprovement: generationsSinceImprovement,
            fitnessTrajectory: fitnessTrajectory, diversityTrajectory: diversityTrajectory,
            isImproving: isImproving, isDiversityCollapsed: isDiversityCollapsed,
            populationSize: populationSize
        )
    }

    func testContinueRunEarly() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = makeSnapshot(currentGeneration: 100) // 10% progress
        let result = coach.coach(snapshot: snapshot)

        XCTAssertEqual(result.advice, .continueRun)
        XCTAssertTrue(result.reasoning.contains("early"))
    }

    func testGAFitnessIsNeverLabelledKcalPerMolWithoutProvenance() {
        let coach = RuleBasedConvergenceCoach()
        // GA fitness is a raw contact-function objective. Without calibrated,
        // receipted provenance it must not be presented in kcal/mol.
        let snapshot = makeSnapshot(currentGeneration: 800,
                                    generationsSinceImprovement: 300,
                                    isImproving: false)

        XCTAssertEqual(snapshot.claimValidity, .proxyOnly)
        XCTAssertEqual(snapshot.fitnessUnitLabel, "CF arbitrary units")

        let result = coach.coach(snapshot: snapshot)
        XCTAssertFalse(result.reasoning.contains("kcal/mol"))
        XCTAssertTrue(result.reasoning.contains("CF arbitrary units"))
    }

    func testGAFitnessKcalPerMolRequiresCalibratedProvenance() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = GAProgressSnapshot(
            currentGeneration: 800, maxGenerations: 1000,
            bestFitness: -8.0, meanFitness: -5.0, populationDiversity: 2.5,
            generationsSinceImprovement: 300,
            fitnessTrajectory: [-7.0, -7.5, -8.0],
            diversityTrajectory: [3.0, 2.8, 2.5],
            isImproving: false, isDiversityCollapsed: false, populationSize: 300,
            scientificProvenance: intelligenceBindingProvenance
        )

        XCTAssertEqual(snapshot.claimValidity, .bindingPhysical)
        XCTAssertTrue(coach.coach(snapshot: snapshot).reasoning.contains("kcal/mol"))
    }

    func testStopEarly() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = makeSnapshot(
            currentGeneration: 800,
            generationsSinceImprovement: 300, // 30% stagnation
            isImproving: false
        )
        let result = coach.coach(snapshot: snapshot)

        XCTAssertEqual(result.advice, .stopEarly)
        XCTAssertEqual(result.estimatedGenerationsRemaining, 0)
    }

    func testDiversityCollapse() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = makeSnapshot(
            currentGeneration: 500,
            generationsSinceImprovement: 200,
            isDiversityCollapsed: true,
            populationSize: 100 // small population
        )
        let result = coach.coach(snapshot: snapshot)

        XCTAssertEqual(result.advice, .increasePopulation)
    }

    func testIncreaseMutationRate() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = makeSnapshot(
            currentGeneration: 500,
            generationsSinceImprovement: 200,
            isDiversityCollapsed: true,
            populationSize: 500 // large population, still collapsed
        )
        let result = coach.coach(snapshot: snapshot)

        XCTAssertEqual(result.advice, .increaseMutationRate)
    }

    func testContinueImproving() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = makeSnapshot(
            currentGeneration: 500,
            generationsSinceImprovement: 10,
            isImproving: true
        )
        let result = coach.coach(snapshot: snapshot)

        XCTAssertEqual(result.advice, .continueRun)
        XCTAssertNotNil(result.estimatedGenerationsRemaining)
    }

    func testSingleGeneration() {
        let coach = RuleBasedConvergenceCoach()
        let snapshot = makeSnapshot(currentGeneration: 1, maxGenerations: 1000)
        let result = coach.coach(snapshot: snapshot)

        XCTAssertEqual(result.advice, .continueRun)
    }
}

// MARK: - FleetExplainer Tests

final class FleetExplainerTests: XCTestCase {

    private func makeDevice(
        model: String = "MacBook Pro", tflops: Double = 5.0,
        thermalState: String = "nominal", batteryPercent: Int? = 85,
        isCharging: Bool = false, computeWeight: Double = 0.4,
        chunksAssigned: Int = 40, chunksCompleted: Int = 30, chunksFailed: Int = 0
    ) -> DeviceSummary {
        DeviceSummary(
            model: model, tflops: tflops, thermalState: thermalState,
            batteryPercent: batteryPercent, isCharging: isCharging,
            computeWeight: computeWeight, chunksAssigned: chunksAssigned,
            chunksCompleted: chunksCompleted, chunksFailed: chunksFailed
        )
    }

    private func makeContext(devices: [DeviceSummary]) -> FleetStatusContext {
        FleetStatusContext(
            devices: devices, totalChunks: 100, completedChunks: 60,
            failedChunks: 0, orphanedChunks: 0, etaSeconds: 1800,
            totalTFLOPS: devices.reduce(0) { $0 + $1.tflops },
            refereeRecommendation: nil
        )
    }

    func testAllHealthy() {
        let explainer = RuleBasedFleetExplainer()
        let devices = [makeDevice(model: "Mac Pro"), makeDevice(model: "MacBook Air", tflops: 2.0)]
        let context = makeContext(devices: devices)
        let result = explainer.explain(context: context)

        XCTAssertTrue(result.bottleneckAnalysis.contains("No bottlenecks"))
        XCTAssertTrue(result.actionItems.contains { $0.contains("No action") })
    }

    func testThrottledDevice() {
        let explainer = RuleBasedFleetExplainer()
        let devices = [
            makeDevice(model: "MacBook Air", thermalState: "serious"),
            makeDevice(model: "Mac Mini", thermalState: "nominal")
        ]
        let context = makeContext(devices: devices)
        let result = explainer.explain(context: context)

        XCTAssertTrue(result.bottleneckAnalysis.contains("MacBook Air"))
        XCTAssertTrue(result.bottleneckAnalysis.contains("throttled"))
    }

    func testLowBattery() {
        let explainer = RuleBasedFleetExplainer()
        let devices = [makeDevice(model: "iPad Pro", batteryPercent: 10, isCharging: false)]
        let context = makeContext(devices: devices)
        let result = explainer.explain(context: context)

        XCTAssertTrue(result.actionItems.contains { $0.contains("Plug in") })
    }

    func testEmptyDevices() {
        let explainer = RuleBasedFleetExplainer()
        let context = FleetStatusContext(
            devices: [], totalChunks: 100, completedChunks: 0,
            failedChunks: 0, orphanedChunks: 0, etaSeconds: nil,
            totalTFLOPS: 0, refereeRecommendation: nil
        )
        let result = explainer.explain(context: context)

        XCTAssertTrue(result.allocationRationale.contains("No devices"))
    }

    func testETAFormatting() {
        let explainer = RuleBasedFleetExplainer()
        let context = FleetStatusContext(
            devices: [makeDevice()], totalChunks: 100, completedChunks: 50,
            failedChunks: 0, orphanedChunks: 0, etaSeconds: 5400,
            totalTFLOPS: 5.0, refereeRecommendation: nil
        )
        let result = explainer.explain(context: context)

        XCTAssertTrue(result.estimatedCompletion.contains("hour"))
    }
}

// MARK: - HealthEntropyInsight Tests

final class HealthEntropyInsightTests: XCTestCase {

    private func makeSnapshot(
        hrvSDNN: Double? = 55, restingHR: Double? = 65, sleepHours: Double? = 7.5,
        shannonS: Double = 0.005, isConverged: Bool = true, modeCount: Int = 3,
        freeEnergy: Double = -8.0
    ) -> HealthEntropySnapshot {
        HealthEntropySnapshot(
            hrvSDNN: hrvSDNN, restingHR: restingHR, sleepHours: sleepHours,
            shannonS: shannonS, isConverged: isConverged, modeCount: modeCount,
            freeEnergy: freeEnergy
        )
    }

    private func makeContext(
        current: HealthEntropySnapshot? = nil,
        history: [HealthEntropySnapshot] = [],
        hrvTrend: String = "stable",
        convergenceTrend: String = "stable"
    ) -> HealthEntropyContext {
        HealthEntropyContext(
            current: current ?? makeSnapshot(),
            history: history,
            hrvTrend: hrvTrend,
            convergenceTrend: convergenceTrend
        )
    }

    func testGoodHRVConverged() {
        let insight = RuleBasedHealthEntropyInsight()
        let context = makeContext(current: makeSnapshot(hrvSDNN: 70, isConverged: true))
        let result = insight.analyze(context: context)

        XCTAssertTrue(result.correlationSummary.contains("Good HRV"))
        XCTAssertTrue(result.dataQualityNote.contains("reliable"))
    }

    func testLowHRVStress() {
        let insight = RuleBasedHealthEntropyInsight()
        let context = makeContext(current: makeSnapshot(hrvSDNN: 25))
        let result = insight.analyze(context: context)

        XCTAssertTrue(result.correlationSummary.contains("Very low HRV"))
    }

    func testPoorSleep() {
        let insight = RuleBasedHealthEntropyInsight()
        let context = makeContext(current: makeSnapshot(sleepHours: 4.5))
        let result = insight.analyze(context: context)

        XCTAssertTrue(result.wellnessRecommendation.contains("Sleep was"))
        XCTAssertTrue(result.wellnessRecommendation.contains("4.5"))
    }

    func testNoHealthData() {
        let insight = RuleBasedHealthEntropyInsight()
        let context = makeContext(current: makeSnapshot(hrvSDNN: nil, restingHR: nil, sleepHours: nil))
        let result = insight.analyze(context: context)

        XCTAssertTrue(result.correlationSummary.contains("No HRV data"))
        XCTAssertEqual(result.confidence, 0.3, accuracy: 0.01)
    }

    func testImprovingTrend() {
        let insight = RuleBasedHealthEntropyInsight()
        let context = makeContext(
            hrvTrend: "improving",
            convergenceTrend: "improving"
        )
        let result = insight.analyze(context: context)

        XCTAssertTrue(result.correlationSummary.contains("positive trend"))
    }

    func testConfidenceWithFullData() {
        let insight = RuleBasedHealthEntropyInsight()
        let history = (0..<3).map { _ in makeSnapshot() }
        let context = makeContext(history: history)
        let result = insight.analyze(context: context)

        XCTAssertEqual(result.confidence, 0.8, accuracy: 0.01)
    }
}

// MARK: - VibrationalInterpreter Tests

final class VibrationalInterpreterTests: XCTestCase {

    private func makeMode(
        index: Int = 0, frequencyCm: Double = 30.0, entropyContribution: Double = 0.0001,
        entropyFraction: Double = 0.4, shiftsOnBinding: Bool = false, shiftDirection: String? = nil
    ) -> NormalModeSummary {
        NormalModeSummary(
            index: index, frequencyCm: frequencyCm, entropyContribution: entropyContribution,
            entropyFraction: entropyFraction, shiftsOnBinding: shiftsOnBinding,
            shiftDirection: shiftDirection
        )
    }

    private func makeContext(
        topModes: [NormalModeSummary]? = nil, totalSVib: Double = 0.001,
        totalSConf: Double = 2.0, vibrationalDominance: Double = 1.5,
        temperature: Double = 298.15, bindingRestrictsMotion: Bool = false,
        scientificProvenance: ScientificProvenance? = intelligenceBindingProvenance
    ) -> VibrationalContext {
        VibrationalContext(
            totalSVib: totalSVib, totalSConf: totalSConf,
            vibrationalDominance: vibrationalDominance,
            temperature: temperature,
            topModes: topModes ?? [makeMode()],
            totalModeCount: 20,
            bindingRestrictsMotion: bindingRestrictsMotion,
            scientificProvenance: scientificProvenance
        )
    }

    func testLowFrequencyDominant() {
        let interpreter = RuleBasedVibrationalInterpreter()
        let context = makeContext(topModes: [makeMode(frequencyCm: 20.0, entropyFraction: 0.5)])
        let result = interpreter.interpret(context: context)

        XCTAssertTrue(result.dominantMotionDescription.contains("very low-frequency"))
    }

    func testHighFrequency() {
        let interpreter = RuleBasedVibrationalInterpreter()
        let context = makeContext(topModes: [makeMode(frequencyCm: 600.0)])
        let result = interpreter.interpret(context: context)

        XCTAssertTrue(result.dominantMotionDescription.contains("high-frequency"))
    }

    func testBindingRestricts() {
        let interpreter = RuleBasedVibrationalInterpreter()
        let modes = [makeMode(shiftsOnBinding: true, shiftDirection: "restricted")]
        let context = makeContext(topModes: modes, bindingRestrictsMotion: true)
        let result = interpreter.interpret(context: context)

        XCTAssertTrue(result.bindingImpact.contains("restricts"))
    }

    func testVibrationalDominance() {
        let interpreter = RuleBasedVibrationalInterpreter()
        let context = makeContext(vibrationalDominance: 4.0)
        let result = interpreter.interpret(context: context)

        XCTAssertTrue(result.designImplication.contains("Backbone flexibility"))
        XCTAssertTrue(result.isEntropicallyDriven)
    }

    func testEmptyModes() {
        let interpreter = RuleBasedVibrationalInterpreter()
        let context = VibrationalContext(
            totalSVib: 0.001, totalSConf: 2.0, vibrationalDominance: 1.0,
            temperature: 298.15, topModes: [], totalModeCount: 0,
            bindingRestrictsMotion: false
        )
        let result = interpreter.interpret(context: context)

        XCTAssertTrue(result.dominantMotionDescription.contains("No significant"))
    }

    func testProxyVibrationalOutputIsDiagnosticAndUnavailable() {
        let interpreter = RuleBasedVibrationalInterpreter()
        let context = makeContext(
            vibrationalDominance: 10,
            bindingRestrictsMotion: true,
            scientificProvenance: nil
        )

        let result = interpreter.interpret(context: context)
        let output = [result.dominantMotionDescription, result.bindingImpact, result.designImplication]
            .joined(separator: " ").lowercased()

        XCTAssertTrue(output.contains("proxy-only") && output.contains("diagnostic"))
        XCTAssertTrue(output.contains("unavailable"))
        XCTAssertFalse(output.contains("kcal/mol"))
        XCTAssertFalse(output.contains("standard lead optimization appropriate"))
        XCTAssertFalse(result.isEntropicallyDriven)
    }
}

// MARK: - SelectivityAnalyst Tests

final class SelectivityAnalystTests: XCTestCase {

    private func makeTarget(
        name: String = "5HT2A", bestFreeEnergy: Double = -10.0,
        modeCount: Int = 5, sConf: Double = 2.5, sVib: Double = 0.001,
        isConverged: Bool = true, cavityVolume: Double? = 500, populationSize: Int = 300,
        scientificProvenance: ScientificProvenance? = intelligenceBindingProvenance
    ) -> TargetDockingSummary {
        TargetDockingSummary(
            targetName: name, bestFreeEnergy: bestFreeEnergy, modeCount: modeCount,
            sConf: sConf, sVib: sVib, isConverged: isConverged,
            cavityVolume: cavityVolume, populationSize: populationSize,
            scientificProvenance: scientificProvenance
        )
    }

    func testEnthalpicDriver() {
        let analyst = RuleBasedSelectivityAnalyst()
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [
                makeTarget(name: "5HT2A", bestFreeEnergy: -15.0, sConf: 2.5),
                makeTarget(name: "D2R", bestFreeEnergy: -8.0, sConf: 2.5)
            ]
        )
        let result = analyst.analyze(context: context)

        XCTAssertEqual(result.preferredTarget, "5HT2A")
        XCTAssertEqual(result.driver, .enthalpic)
        XCTAssertTrue(result.explanation.contains("enthalpy-driven"))
        XCTAssertEqual(result.deltaGClaimValidity, .bindingPhysical)
        XCTAssertTrue(result.isDeltaGPhysical)
        XCTAssertEqual(result.deltaGUnitLabel, "kcal/mol")
    }

    func testRetainedDeltaGIsMarkedUnavailableWithoutProvenance() {
        let analyst = RuleBasedSelectivityAnalyst()
        // One un-evidenced target is enough to disqualify the comparison.
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [
                makeTarget(name: "5HT2A", bestFreeEnergy: -15.0),
                makeTarget(name: "D2R", bestFreeEnergy: -8.0,
                           scientificProvenance: nil)
            ]
        )
        let result = analyst.analyze(context: context)

        XCTAssertEqual(result.preferredTarget, "unavailable")
        XCTAssertEqual(result.deltaGClaimValidity, .proxyOnly)
        XCTAssertFalse(result.isDeltaGPhysical)
        XCTAssertFalse(result.deltaGUnitLabel.contains("kcal/mol"))
        // The number itself is retained for compatibility, unchanged.
        XCTAssertEqual(result.deltaG, -7.0, accuracy: 1e-9)
    }

    func testDecodedDeltaGMarkerFailsClosedWhenAbsent() throws {
        // A payload predating the marker must not be read as physical.
        let json = """
        {"preferredTarget":"5HT2A","deltaG":-2.5,"driver":"enthalpic",
         "explanation":"x","designSuggestion":"y"}
        """.data(using: .utf8)!
        let decoded = try JSONDecoder().decode(
            CrossPlatformSelectivityAnalysis.self, from: json)

        XCTAssertEqual(decoded.deltaG, -2.5, accuracy: 1e-10)
        XCTAssertEqual(decoded.deltaGClaimValidity, .proxyOnly)
        XCTAssertFalse(decoded.isDeltaGPhysical)
    }

    func testEntropicDriver() {
        let analyst = RuleBasedSelectivityAnalyst()
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [
                makeTarget(name: "5HT2A", bestFreeEnergy: -10.0, sConf: 5.0),
                makeTarget(name: "D2R", bestFreeEnergy: -9.5, sConf: 1.0)
            ]
        )
        let result = analyst.analyze(context: context)

        XCTAssertEqual(result.preferredTarget, "5HT2A")
        XCTAssertEqual(result.driver, .entropic)
    }

    func testInconclusive() {
        let analyst = RuleBasedSelectivityAnalyst()
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [
                makeTarget(name: "5HT2A", bestFreeEnergy: -10.0, sConf: 2.5),
                makeTarget(name: "D2R", bestFreeEnergy: -10.1, sConf: 2.5)
            ]
        )
        let result = analyst.analyze(context: context)

        XCTAssertTrue(result.driver == .inconclusive || result.driver == .mixed)
    }

    func testSingleTarget() {
        let analyst = RuleBasedSelectivityAnalyst()
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [makeTarget(name: "5HT2A")]
        )
        let result = analyst.analyze(context: context)

        XCTAssertEqual(result.driver, .inconclusive)
        XCTAssertTrue(result.explanation.contains("Only one target"))
    }

    func testNotConverged() {
        let analyst = RuleBasedSelectivityAnalyst()
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [
                makeTarget(name: "5HT2A", bestFreeEnergy: -12.0, isConverged: false),
                makeTarget(name: "D2R", bestFreeEnergy: -8.0, isConverged: true)
            ]
        )
        let result = analyst.analyze(context: context)

        XCTAssertEqual(result.driver, .inconclusive)
        XCTAssertTrue(result.explanation.contains("not converged"))
    }

    func testCodableRoundTrip() throws {
        let context = SelectivityContext(
            ligandName: "psilocin",
            targets: [
                makeTarget(name: "5HT2A", bestFreeEnergy: -10.0),
                makeTarget(name: "D2R", bestFreeEnergy: -8.0)
            ]
        )

        let data = try JSONEncoder().encode(context)
        let decoded = try JSONDecoder().decode(SelectivityContext.self, from: data)

        XCTAssertEqual(decoded.ligandName, context.ligandName)
        XCTAssertEqual(decoded.targets.count, 2)
        XCTAssertEqual(decoded.deltaDeltas.count, 1)
        XCTAssertEqual(decoded.deltaDeltas[0].targetA, "5HT2A")
        XCTAssertEqual(decoded.deltaDeltas[0].targetB, "D2R")
        XCTAssertEqual(decoded.deltaDeltas[0].ddg, -2.0, accuracy: 1e-10)
    }

    func testProxySelectivitySuppressesPreferenceDeltaGAndOptimizationClaims() {
        let analyst = RuleBasedSelectivityAnalyst()
        let context = SelectivityContext(
            ligandName: "hostile",
            targets: [
                makeTarget(name: "A", bestFreeEnergy: -50, scientificProvenance: nil),
                makeTarget(name: "B", bestFreeEnergy: 50, scientificProvenance: nil),
            ]
        )

        let result = analyst.analyze(context: context)
        let output = [result.explanation, result.designSuggestion].joined(separator: " ").lowercased()

        XCTAssertEqual(result.preferredTarget, "unavailable")
        XCTAssertEqual(result.driver, .inconclusive)
        XCTAssertTrue(output.contains("proxy") && output.contains("unavailable"))
        XCTAssertFalse(output.contains("enthalpy-driven selectivity"))
        XCTAssertFalse(output.contains("improve affinity"))
    }
}

// MARK: - CampaignJournalist Tests

final class CampaignJournalistTests: XCTestCase {

    private func makeRun(
        index: Int = 0, freeEnergy: Double = -8.0, entropy: Double = 0.005,
        isConverged: Bool = true, modeCount: Int = 3, confidence: String = "high",
        scientificProvenance: ScientificProvenance? = intelligenceBindingProvenance
    ) -> RunSnapshot {
        RunSnapshot(
            runIndex: index, freeEnergy: freeEnergy, entropy: entropy,
            isConverged: isConverged, modeCount: modeCount,
            confidence: confidence, timestamp: Date(),
            scientificProvenance: scientificProvenance
        )
    }

    func testImproving() {
        let journalist = RuleBasedCampaignJournalist()
        let runs = [
            makeRun(index: 0, freeEnergy: -5.0),
            makeRun(index: 1, freeEnergy: -6.0),
            makeRun(index: 2, freeEnergy: -8.0)
        ]
        let context = CampaignContext(campaignKey: "5HT2A-psilocin", runs: runs)
        let result = journalist.summarize(context: context)

        XCTAssertEqual(result.trend, "improving")
        XCTAssertTrue(result.nextStepRecommendation.contains("Continue"))
    }

    func testStagnating() {
        let journalist = RuleBasedCampaignJournalist()
        let runs = [
            makeRun(index: 0, freeEnergy: -8.0),
            makeRun(index: 1, freeEnergy: -8.1),
            makeRun(index: 2, freeEnergy: -8.2),
            makeRun(index: 3, freeEnergy: -8.15)
        ]
        let context = CampaignContext(campaignKey: "5HT2A-psilocin", runs: runs)
        let result = journalist.summarize(context: context)

        XCTAssertEqual(result.trend, "stagnating")
    }

    func testRegressing() {
        let journalist = RuleBasedCampaignJournalist()
        let runs = [
            makeRun(index: 0, freeEnergy: -10.0),
            makeRun(index: 1, freeEnergy: -8.0)
        ]
        let context = CampaignContext(campaignKey: "5HT2A-psilocin", runs: runs)
        let result = journalist.summarize(context: context)

        XCTAssertEqual(result.trend, "regressing")
        XCTAssertTrue(result.nextStepRecommendation.contains("Revert"))
    }

    func testSingleRun() {
        let journalist = RuleBasedCampaignJournalist()
        let runs = [makeRun(index: 0)]
        let context = CampaignContext(campaignKey: "5HT2A-psilocin", runs: runs)
        let result = journalist.summarize(context: context)

        XCTAssertEqual(result.trend, "insufficient data")
    }

    func testPublicationReady() {
        let journalist = RuleBasedCampaignJournalist()
        let runs = (0..<5).map { i in makeRun(index: i, freeEnergy: -8.0 - Double(i) * 0.1, isConverged: true) }
        let context = CampaignContext(campaignKey: "5HT2A-psilocin", runs: runs)
        let result = journalist.summarize(context: context)

        XCTAssertTrue(result.readyForPublication)
    }

    func testEmptyRuns() {
        let journalist = RuleBasedCampaignJournalist()
        let context = CampaignContext(campaignKey: "test", runs: [])
        let result = journalist.summarize(context: context)

        XCTAssertEqual(result.runCount, 0)
        XCTAssertTrue(result.progressNarrative.contains("No runs"))
    }

    func testProxyCampaignSuppressesAffinityUnitsAndLeadOptimization() {
        let journalist = RuleBasedCampaignJournalist()
        let runs = [
            makeRun(index: 0, freeEnergy: -5, scientificProvenance: nil),
            makeRun(index: 1, freeEnergy: -10, scientificProvenance: nil),
        ]
        let result = journalist.summarize(
            context: CampaignContext(campaignKey: "hostile", runs: runs)
        )
        let output = [result.progressNarrative, result.bestResult, result.nextStepRecommendation]
            .joined(separator: " ").lowercased()

        XCTAssertTrue(output.contains("proxy") && output.contains("unavailable"))
        XCTAssertFalse(output.contains("binding affinity improved"))
        XCTAssertFalse(output.contains("kcal/mol"))
        XCTAssertFalse(output.contains("proceed with lead optimization"))
        XCTAssertFalse(result.readyForPublication)
    }
}

// MARK: - LigandFitCritic Tests

final class LigandFitCriticTests: XCTestCase {

    private func makePose(rank: Int, cfScore: Double, boltzmannWeight: Double,
                          rmsdToCentroid: Double = 1.0) -> PoseProfile {
        PoseProfile(rank: rank, cfScore: cfScore, boltzmannWeight: boltzmannWeight,
                    rmsdToCentroid: rmsdToCentroid)
    }

    private func makeContext(
        poses: [PoseProfile], modeIndex: Int = 0, modeFreeEnergy: Double = -8.0
    ) -> PoseQualityContext {
        PoseQualityContext(
            modeIndex: modeIndex,
            topPoses: poses,
            totalPoses: poses.count,
            modeFreeEnergy: modeFreeEnergy
        )
    }

    func testDominantPose() {
        let critic = RuleBasedLigandFitCritic()
        let poses = [
            makePose(rank: 0, cfScore: -10.0, boltzmannWeight: 0.7, rmsdToCentroid: 0.5),
            makePose(rank: 1, cfScore: -8.0, boltzmannWeight: 0.2, rmsdToCentroid: 1.0),
            makePose(rank: 2, cfScore: -6.0, boltzmannWeight: 0.1, rmsdToCentroid: 1.5)
        ]
        let context = makeContext(poses: poses)
        let result = critic.evaluate(context: context)

        XCTAssertTrue(result.poseConsensus == "strong" || result.poseConsensus == "moderate")
        XCTAssertGreaterThan(result.confidenceInTopPose, 0.7)
    }

    func testAmbiguousPoses() {
        let critic = RuleBasedLigandFitCritic()
        let poses = [
            makePose(rank: 0, cfScore: -8.0, boltzmannWeight: 0.25, rmsdToCentroid: 1.0),
            makePose(rank: 1, cfScore: -7.8, boltzmannWeight: 0.25, rmsdToCentroid: 3.0),
            makePose(rank: 2, cfScore: -7.5, boltzmannWeight: 0.25, rmsdToCentroid: 5.0),
            makePose(rank: 3, cfScore: -7.2, boltzmannWeight: 0.25, rmsdToCentroid: 2.0)
        ]
        let context = makeContext(poses: poses)
        let result = critic.evaluate(context: context)

        XCTAssertTrue(result.poseConsensus == "ambiguous" || result.poseConsensus == "weak")
    }

    func testScoreWeightAligned() {
        let critic = RuleBasedLigandFitCritic()
        // Best CF score matches highest Boltzmann weight
        let poses = [
            makePose(rank: 0, cfScore: -12.0, boltzmannWeight: 0.8),
            makePose(rank: 1, cfScore: -8.0, boltzmannWeight: 0.2)
        ]
        let context = makeContext(poses: poses)

        XCTAssertTrue(context.scoreWeightAligned)
        XCTAssertGreaterThan(context.scoreWeightCorrelation, 0.9)
    }

    func testScoreWeightMisaligned() {
        let critic = RuleBasedLigandFitCritic()
        // Best CF score does NOT match highest Boltzmann weight
        let poses = [
            makePose(rank: 0, cfScore: -12.0, boltzmannWeight: 0.2),
            makePose(rank: 1, cfScore: -6.0, boltzmannWeight: 0.8)
        ]
        let context = makeContext(poses: poses)

        XCTAssertFalse(context.scoreWeightAligned)
        XCTAssertLessThan(context.scoreWeightCorrelation, 0.0)
    }

    func testEmptyPoses() {
        let critic = RuleBasedLigandFitCritic()
        let context = makeContext(poses: [])
        let result = critic.evaluate(context: context)

        XCTAssertEqual(result.confidenceInTopPose, 0.0)
        XCTAssertTrue(result.topPoseSummary.contains("No poses"))
    }

    func testPoseScoresAreNotKcalPerMolWithoutProvenance() {
        let critic = RuleBasedLigandFitCritic()
        let poses = [
            makePose(rank: 0, cfScore: -10.0, boltzmannWeight: 0.7, rmsdToCentroid: 0.5),
            makePose(rank: 1, cfScore: -8.0, boltzmannWeight: 0.2, rmsdToCentroid: 1.0)
        ]
        let context = makeContext(poses: poses)

        XCTAssertEqual(context.claimValidity, .proxyOnly)
        XCTAssertEqual(context.energyUnitLabel, "CF arbitrary units")

        let result = critic.evaluate(context: context)
        XCTAssertFalse(result.topPoseSummary.contains("kcal/mol"))
        XCTAssertTrue(result.medicinalChemistryNote.contains("unavailable"))
    }

    func testMedicinalChemistryAdviceRequiresBindingProvenance() {
        let critic = RuleBasedLigandFitCritic()
        let poses = [
            makePose(rank: 0, cfScore: -10.0, boltzmannWeight: 0.7, rmsdToCentroid: 0.5),
            makePose(rank: 1, cfScore: -8.0, boltzmannWeight: 0.2, rmsdToCentroid: 0.6)
        ]
        let context = PoseQualityContext(
            modeIndex: 0, topPoses: poses, totalPoses: poses.count,
            modeFreeEnergy: -8.0,
            scientificProvenance: intelligenceBindingProvenance
        )

        let result = critic.evaluate(context: context)
        XCTAssertTrue(result.topPoseSummary.contains("kcal/mol"))
        XCTAssertFalse(result.medicinalChemistryNote.contains("unavailable"))
    }

    func testSpearmanCorrelationPerfect() {
        // Perfect correlation: CF rank matches weight rank exactly
        let poses = [
            makePose(rank: 0, cfScore: -10.0, boltzmannWeight: 0.5),
            makePose(rank: 1, cfScore: -8.0, boltzmannWeight: 0.3),
            makePose(rank: 2, cfScore: -6.0, boltzmannWeight: 0.2)
        ]
        let context = makeContext(poses: poses)

        XCTAssertEqual(context.scoreWeightCorrelation, 1.0, accuracy: 0.01)
    }

    func testConfidenceClamped() {
        let critic = RuleBasedLigandFitCritic()
        let poses = [
            makePose(rank: 0, cfScore: -10.0, boltzmannWeight: 0.7, rmsdToCentroid: 0.5),
            makePose(rank: 1, cfScore: -8.0, boltzmannWeight: 0.2, rmsdToCentroid: 0.6)
        ]
        let context = makeContext(poses: poses)
        let result = critic.evaluate(context: context)

        XCTAssertGreaterThanOrEqual(result.confidenceInTopPose, 0.0)
        XCTAssertLessThanOrEqual(result.confidenceInTopPose, 0.95)
    }
}

// MARK: - Codable Round-Trip Tests

final class IntelligenceCodableTests: XCTestCase {

    func testCrossPlatformModeNarrativeCodable() throws {
        let original = CrossPlatformModeNarrative(
            modeDescriptions: [
                CrossPlatformModeDescription(characterization: "tight lock", optimizationHint: "rigidify")
            ],
            selectivityInsight: "Mode 1 dominates.",
            confidence: 0.8
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformModeNarrative.self, from: data)

        XCTAssertEqual(decoded.modeDescriptions.count, 1)
        XCTAssertEqual(decoded.confidence, 0.8, accuracy: 1e-10)
    }

    func testCrossPlatformCleftAssessmentCodable() throws {
        let original = CrossPlatformCleftAssessment(
            druggability: .high, summary: "Drug-like pocket",
            suggestedLigandProperties: "MW 300-500", warnings: ["Shallow"]
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformCleftAssessment.self, from: data)

        XCTAssertEqual(decoded.druggability, .high)
        XCTAssertEqual(decoded.warnings.count, 1)
    }

    func testCrossPlatformConvergenceCoachingCodable() throws {
        let original = CrossPlatformConvergenceCoaching(
            advice: .stopEarly, reasoning: "Stagnated",
            estimatedGenerationsRemaining: 0, confidence: 0.9
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformConvergenceCoaching.self, from: data)

        XCTAssertEqual(decoded.advice, .stopEarly)
        XCTAssertEqual(decoded.estimatedGenerationsRemaining, 0)
    }

    func testCrossPlatformFleetExplanationCodable() throws {
        let original = CrossPlatformFleetExplanation(
            allocationRationale: "By TFLOPS", bottleneckAnalysis: "None",
            actionItems: ["No action"], estimatedCompletion: "~30 min"
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformFleetExplanation.self, from: data)

        XCTAssertEqual(decoded.actionItems.count, 1)
    }

    func testCrossPlatformHealthEntropyInsightCodable() throws {
        let original = CrossPlatformHealthEntropyInsight(
            correlationSummary: "Good HRV", wellnessRecommendation: "Rest",
            dataQualityNote: "Reliable", confidence: 0.8
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformHealthEntropyInsight.self, from: data)

        XCTAssertEqual(decoded.confidence, 0.8, accuracy: 1e-10)
    }

    func testCrossPlatformVibrationalInsightCodable() throws {
        let original = CrossPlatformVibrationalInsight(
            dominantMotionDescription: "Loop breathing",
            bindingImpact: "Restricts motion",
            designImplication: "Use smaller fragment",
            isEntropicallyDriven: true
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformVibrationalInsight.self, from: data)

        XCTAssertTrue(decoded.isEntropicallyDriven)
    }

    func testCrossPlatformSelectivityAnalysisCodable() throws {
        let original = CrossPlatformSelectivityAnalysis(
            preferredTarget: "5HT2A", deltaG: -2.5,
            driver: .enthalpic, explanation: "Stronger interactions",
            designSuggestion: "Optimize geometry"
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformSelectivityAnalysis.self, from: data)

        XCTAssertEqual(decoded.driver, .enthalpic)
        XCTAssertEqual(decoded.deltaG, -2.5, accuracy: 1e-10)
    }

    func testCrossPlatformCampaignSummaryCodable() throws {
        let original = CrossPlatformCampaignSummary(
            campaignKey: "5HT2A-psilocin", runCount: 5,
            progressNarrative: "Improving", bestResult: "Run 3",
            trend: "improving", nextStepRecommendation: "Continue",
            readyForPublication: true
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformCampaignSummary.self, from: data)

        XCTAssertEqual(decoded.campaignKey, "5HT2A-psilocin")
        XCTAssertTrue(decoded.readyForPublication)
    }

    func testCrossPlatformPoseQualityReportCodable() throws {
        let original = CrossPlatformPoseQualityReport(
            topPoseSummary: "Pose 1 dominates",
            poseConsensus: "strong",
            scoreWeightAlignment: "Well aligned",
            confidenceInTopPose: 0.85,
            medicinalChemistryNote: "High confidence"
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(CrossPlatformPoseQualityReport.self, from: data)

        XCTAssertEqual(decoded.poseConsensus, "strong")
        XCTAssertEqual(decoded.confidenceInTopPose, 0.85, accuracy: 1e-10)
    }

    func testDeltaDeltaGCodable() throws {
        let original = DeltaDeltaG(targetA: "5HT2A", targetB: "D2R", ddg: -3.5)
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(DeltaDeltaG.self, from: data)

        XCTAssertEqual(decoded, original)
    }

    func testHealthEntropySnapshotCodableWithNils() throws {
        let original = HealthEntropySnapshot(
            hrvSDNN: nil, restingHR: nil, sleepHours: nil,
            shannonS: 0.005, isConverged: true, modeCount: 3, freeEnergy: -8.0
        )
        let data = try JSONEncoder().encode(original)
        let decoded = try JSONDecoder().decode(HealthEntropySnapshot.self, from: data)

        XCTAssertNil(decoded.hrvSDNN)
        XCTAssertNil(decoded.restingHR)
        XCTAssertNil(decoded.sleepHours)
        XCTAssertEqual(decoded.freeEnergy, -8.0, accuracy: 1e-10)
    }
}
