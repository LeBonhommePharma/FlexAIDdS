// IntelligencePipelineTests.swift — Tests for unified IntelligencePipeline
//
// Tests that IntelligencePipeline correctly orchestrates all rule-based
// features, isolates failures, and produces correct aggregate results.
//
// NOTE: The FlexAIDdSTests test target in Package.swift currently depends only
// on "FlexAIDdS", not "Intelligence". The @testable import Intelligence below
// will fail to build until the test target dependency is updated in Package.swift:
//
//   .testTarget(
//       name: "FlexAIDdSTests",
//       dependencies: ["FlexAIDdS", "Intelligence"],
//       path: "Tests/FlexAIDdSTests"
//   )
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import XCTest
import FlexAIDdS
@testable import Intelligence

final class IntelligencePipelineTests: XCTestCase {

    // MARK: - Factories

    private func makeBindingMode(
        freeEnergy: Double = -8.5,
        entropy: Double = 0.005,
        enthalpy: Double = -10.0,
        heatCapacity: Double = 0.02,
        size: Int = 15
    ) -> BindingModeResult {
        BindingModeResult(
            freeEnergy: freeEnergy,
            entropy: entropy,
            enthalpy: enthalpy,
            heatCapacity: heatCapacity,
            size: size,
            poses: (0..<size).map { i in
                PoseResult(cf: freeEnergy + Double(i) * 0.5, reachDist: Float(i) * 0.3)
            }
        )
    }

    private func makeDockingResult(modeCount: Int = 2) -> DockingResult {
        let modes = (0..<modeCount).map { i in
            makeBindingMode(freeEnergy: -8.5 + Double(i) * 2.0)
        }
        return DockingResult(
            bindingModes: modes,
            globalThermodynamics: ThermodynamicResult(
                temperature: 298.15, logZ: 10.0, freeEnergy: -8.0,
                meanEnergy: -9.5, meanEnergySq: 92.0,
                heatCapacity: 0.03, entropy: 0.006, stdEnergy: 1.2
            ),
            temperature: 298.15,
            populationSize: 500
        )
    }

    private func makeGASnapshot() -> GAProgressSnapshot {
        GAProgressSnapshot(
            currentGeneration: 150,
            maxGenerations: 200,
            bestFitness: -12.5,
            meanFitness: -8.3,
            populationDiversity: 0.65,
            generationsSinceImprovement: 60,
            fitnessTrajectory: [-10.0, -11.0, -11.5, -12.0, -12.5],
            diversityTrajectory: [0.8, 0.75, 0.7, 0.68, 0.65],
            isImproving: false,
            isDiversityCollapsed: false,
            populationSize: 300
        )
    }

    private func makeCleftFeatures(
        scientificProvenance: ScientificProvenance? = nil
    ) -> CleftFeatures {
        CleftFeatures(
            volume: 500, depth: 7.0, sphereCount: 25,
            maxSphereRadius: 3.5, hydrophobicFraction: 0.55,
            anchorResidueCount: 6, elongation: 0.4, solventExposure: 0.25,
            scientificProvenance: scientificProvenance
        )
    }

    private static let pipelineBindingProvenance = ScientificProvenance(
        energyDomain: .calibratedKcalPerMol,
        ensembleMeasure: .enumeratedMicrostates,
        referenceState: .matchedAssociationCycle,
        energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
        measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4",
        referenceProvenance: "sha256:01d26f66e709d388d7b971de6204680694e7cca10b1570506a5738c6909a6442"
    )

    // MARK: - Pipeline Tests

    func testPipelineGeneratesModeNarrative() {
        let pipeline = IntelligencePipeline()
        let result = pipeline.analyze(dockingResult: makeDockingResult(modeCount: 2))

        XCTAssertNotNil(result.modeNarrative, "Pipeline should produce a mode narrative for non-empty modes")
        XCTAssertEqual(result.modeNarrative?.modeDescriptions.count, 2,
                       "Mode narrative should contain one description per binding mode")
        XCTAssertTrue(result.errors.isEmpty, "No errors expected: \(result.errors)")
    }

    func testPipelineGeneratesPoseQuality() {
        let pipeline = IntelligencePipeline()
        let docking = makeDockingResult(modeCount: 2)
        let result = pipeline.analyze(dockingResult: docking)

        XCTAssertFalse(result.poseQualityReports.isEmpty,
                       "Should produce pose quality reports when modes have poses")
        XCTAssertLessThanOrEqual(result.poseQualityReports.count, 3,
                                 "Pipeline evaluates at most 3 modes for pose quality")
    }

    func testPipelineWithConvergenceCoaching() {
        let pipeline = IntelligencePipeline()
        let result = pipeline.analyze(
            dockingResult: makeDockingResult(),
            gaSnapshot: makeGASnapshot()
        )

        XCTAssertNotNil(result.convergenceCoaching,
                        "Should produce convergence coaching when GAProgressSnapshot is provided")
        // Snapshot: gen 150/200, 60 gens since improvement, isImproving=false
        // 60/200 = 30% stagnation → stopEarly
        XCTAssertEqual(result.convergenceCoaching?.advice, .stopEarly)
    }

    func testPipelineWithCleftAssessment() {
        let pipeline = IntelligencePipeline()
        let result = pipeline.analyze(
            dockingResult: makeDockingResult(),
            cleftFeatures: makeCleftFeatures(
                scientificProvenance: Self.pipelineBindingProvenance)
        )

        XCTAssertNotNil(result.cleftAssessment,
                        "Should produce cleft assessment when CleftFeatures is provided")
        // Volume 500, depth 7.0, hydrophobic 0.55, anchors 6, low solvent exposure
        XCTAssertEqual(result.cleftAssessment?.druggability, .high)
    }

    func testPipelineCleftAssessmentFailsClosedWithoutProvenance() {
        let pipeline = IntelligencePipeline()
        let result = pipeline.analyze(
            dockingResult: makeDockingResult(),
            cleftFeatures: makeCleftFeatures()
        )

        // Same geometry, no evidence: no druggability tier may be issued.
        XCTAssertEqual(result.cleftAssessment?.druggability, .unavailable)
    }

    func testPipelineWithEmptyBindingModes() {
        let pipeline = IntelligencePipeline()
        let empty = DockingResult(
            bindingModes: [],
            globalThermodynamics: ThermodynamicResult(
                temperature: 298.15, logZ: 0, freeEnergy: 0,
                meanEnergy: 0, meanEnergySq: 0,
                heatCapacity: 0, entropy: 0, stdEnergy: 0
            ),
            temperature: 298.15,
            populationSize: 0
        )
        let result = pipeline.analyze(dockingResult: empty)

        XCTAssertNil(result.modeNarrative, "Mode narrative should be nil for empty binding modes")
        XCTAssertTrue(result.errors.contains { $0.contains("BindingModeNarrator") },
                      "Errors should mention BindingModeNarrator when no modes available")
    }

    func testPipelineIsolatesFailures() {
        let pipeline = IntelligencePipeline()
        let empty = DockingResult(
            bindingModes: [],
            globalThermodynamics: ThermodynamicResult(
                temperature: 298.15, logZ: 0, freeEnergy: 0,
                meanEnergy: 0, meanEnergySq: 0,
                heatCapacity: 0, entropy: 0, stdEnergy: 0
            ),
            temperature: 298.15,
            populationSize: 0
        )
        // Even with empty modes, convergence coaching should still work if snapshot is provided
        let result = pipeline.analyze(
            dockingResult: empty,
            gaSnapshot: makeGASnapshot()
        )

        XCTAssertNil(result.modeNarrative, "Mode narrative should fail for empty modes")
        XCTAssertNotNil(result.convergenceCoaching,
                        "Convergence coaching should succeed independently of mode narrative failure")
    }

    // MARK: - Codable Round-Trip

    func testIntelligenceResultCodable() throws {
        let pipeline = IntelligencePipeline()
        let original = pipeline.analyze(
            dockingResult: makeDockingResult(),
            gaSnapshot: makeGASnapshot(),
            cleftFeatures: makeCleftFeatures()
        )

        let encoder = JSONEncoder()
        let data = try encoder.encode(original)

        let decoder = JSONDecoder()
        let decoded = try decoder.decode(IntelligenceResult.self, from: data)

        // Verify mode narrative round-trips
        XCTAssertEqual(decoded.modeNarrative?.modeDescriptions.count,
                       original.modeNarrative?.modeDescriptions.count)
        XCTAssertEqual(decoded.modeNarrative?.confidence,
                       original.modeNarrative?.confidence)

        // Verify convergence coaching round-trips
        XCTAssertEqual(decoded.convergenceCoaching?.advice,
                       original.convergenceCoaching?.advice)
        XCTAssertEqual(decoded.convergenceCoaching?.estimatedGenerationsRemaining,
                       original.convergenceCoaching?.estimatedGenerationsRemaining)

        // Verify cleft assessment round-trips
        XCTAssertEqual(decoded.cleftAssessment?.druggability,
                       original.cleftAssessment?.druggability)
        XCTAssertEqual(decoded.cleftAssessment?.warnings.count,
                       original.cleftAssessment?.warnings.count)

        // Verify pose quality reports round-trip
        XCTAssertEqual(decoded.poseQualityReports.count,
                       original.poseQualityReports.count)

        // Verify errors round-trip
        XCTAssertEqual(decoded.errors.count, original.errors.count)
        XCTAssertEqual(decoded.errors, original.errors)
    }
}
