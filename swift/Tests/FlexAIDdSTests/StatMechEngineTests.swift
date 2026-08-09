// StatMechEngineTests.swift — Unit tests for the Swift StatMechEngine wrapper
//
// Mirrors test cases from tests/test_statmech.cpp to verify the C shim layer.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import XCTest
@testable import FlexAIDdS

final class StatMechEngineTests: XCTestCase {

    // MARK: - Basic Thermodynamics

    func testTwoStateSystemAt300K() async {
        let runner = FlexAIDRunner(temperature: 300.0)

        // Two-state system: E1 = -10.0, E2 = -8.0 kcal/mol
        await runner.addSample(energy: -10.0)
        await runner.addSample(energy: -8.0)

        let result = await runner.compute()

        XCTAssertEqual(result.temperature, 300.0, accuracy: 1e-10)
        XCTAssertEqual(result.scientificProvenance?.energyDomain, .unclassified)
        XCTAssertEqual(result.scientificProvenance?.ensembleMeasure, .unclassified)
        XCTAssertEqual(result.claimValidity, .proxyOnly)
        // Free energy should be lower than both energies (entropy contribution)
        XCTAssertLessThan(result.freeEnergy, -8.0)
        // Mean energy should be between the two
        XCTAssertGreaterThanOrEqual(result.meanEnergy, -10.0)
        XCTAssertLessThanOrEqual(result.meanEnergy, -8.0)
        // Entropy should be positive
        XCTAssertGreaterThan(result.entropy, 0)
        // Heat capacity should be non-negative
        XCTAssertGreaterThanOrEqual(result.heatCapacity, 0)
    }

    func testSingleStateEnsemble() async {
        let runner = FlexAIDRunner(temperature: 300.0)
        await runner.addSample(energy: -5.0)

        let result = await runner.compute()

        // Single state: F = E, S = 0, Cv = 0
        XCTAssertEqual(result.freeEnergy, -5.0, accuracy: 1e-10)
        XCTAssertEqual(result.meanEnergy, -5.0, accuracy: 1e-10)
        XCTAssertEqual(result.entropy, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.heatCapacity, 0.0, accuracy: 1e-10)
    }

    func testDegenerateStates() async {
        let runner = FlexAIDRunner(temperature: 300.0)
        await runner.addSample(energy: -5.0, multiplicity: 10)

        let result = await runner.compute()

        // All degenerate: <E> = E, S = kB * ln(10)
        XCTAssertEqual(result.meanEnergy, -5.0, accuracy: 1e-10)
        let expectedS = kBkcal * log(10.0)
        XCTAssertEqual(result.entropy, expectedS, accuracy: 1e-8)
    }

    // MARK: - Boltzmann Weights

    func testBoltzmannWeightsNormalize() async {
        let runner = FlexAIDRunner(temperature: 300.0)
        await runner.addSample(energy: -10.0)
        await runner.addSample(energy: -8.0)
        await runner.addSample(energy: -6.0)

        let weights = await runner.boltzmannWeights()
        XCTAssertEqual(weights.count, 3)

        let total = weights.reduce(0, +)
        XCTAssertEqual(total, 1.0, accuracy: 1e-10)

        // Lower energy should have higher weight
        XCTAssertGreaterThan(weights[0], weights[1])
        XCTAssertGreaterThan(weights[1], weights[2])
    }

    // MARK: - Sample Management

    func testClearResetsSamples() async {
        let runner = FlexAIDRunner(temperature: 300.0)
        await runner.addSample(energy: -10.0)
        let countBeforeClear = await runner.sampleCount
        XCTAssertEqual(countBeforeClear, 1)

        await runner.clear()
        let countAfterClear = await runner.sampleCount
        XCTAssertEqual(countAfterClear, 0)
    }

    // MARK: - Static Functions

    func testHelmholtzFromArray() {
        let energies = [-10.0, -8.0, -6.0]
        let F = FlexAIDRunner.helmholtz(energies: energies, temperature: 300.0)

        // F should be lower than the minimum energy
        XCTAssertLessThan(F, -6.0)
        // F should be lower than or close to the lowest energy
        XCTAssertLessThanOrEqual(F, -10.0 + 1.0) // within ~kT of lowest
    }

    func testThermodynamicIntegration() {
        // Linear integrand: dV/dlambda = 2*lambda
        // Integral from 0 to 1 = 1.0
        let points = [
            TIPoint(lambda: 0.0, dVdLambda: 0.0),
            TIPoint(lambda: 0.5, dVdLambda: 1.0),
            TIPoint(lambda: 1.0, dVdLambda: 2.0),
        ]

        let deltaG = FlexAIDRunner.thermodynamicIntegration(points: points)
        XCTAssertEqual(deltaG, 1.0, accuracy: 1e-10)
    }

    // MARK: - Accessors

    func testTemperatureAndBeta() async {
        let runner = FlexAIDRunner(temperature: 300.0)
        let temperature = await runner.temperature
        XCTAssertEqual(temperature, 300.0, accuracy: 1e-10)

        let expectedBeta = 1.0 / (kBkcal * 300.0)
        let beta = await runner.beta
        XCTAssertEqual(beta, expectedBeta, accuracy: 1e-6)
    }

    // MARK: - Boltzmann LUT

    func testBoltzmannLookupAccuracy() {
        let beta = 1.0 / (kBkcal * 300.0)
        let lut = BoltzmannLookup(beta: beta, eMin: -20.0, eMax: 0.0)

        // Test several energies
        for energy in stride(from: -20.0, through: 0.0, by: 1.0) {
            let exact = exp(-beta * energy)
            let lookup = lut.lookup(energy: energy)
            XCTAssertEqual(lookup, exact, accuracy: exact * 0.01) // 1% tolerance
        }
    }

    // MARK: - ENCoM

    func testENCoMTotalEntropy() {
        let sConf = 0.005  // kcal/mol/K
        let sVib = 0.003   // kcal/mol/K
        let sTotal = ENCoMRunner.totalEntropy(configurational: sConf, vibrational: sVib)
        XCTAssertEqual(sTotal, sConf + sVib, accuracy: 1e-15)
    }

    func testENCoMFreeEnergyCorrection() {
        let fElec = -10.0  // kcal/mol
        let sVib = 0.005   // kcal/mol/K
        let T = 300.0

        let fTotal = ENCoMRunner.freeEnergyWithVibrations(
            electronic: fElec, vibrationalEntropy: sVib, temperature: T)

        // F_total = F_elec - T*S_vib = -10.0 - 300.0*0.005 = -11.5
        XCTAssertEqual(fTotal, fElec - T * sVib, accuracy: 1e-10)
    }

    // MARK: - Model Codable

    func testThermodynamicResultCodable() throws {
        let result = ThermodynamicResult(
            temperature: 300.0, logZ: 10.0, freeEnergy: -15.5,
            meanEnergy: -12.0, meanEnergySq: 150.0,
            heatCapacity: 0.002, entropy: 0.005, stdEnergy: 1.2)

        let encoded = try JSONEncoder().encode(result)
        let decoded = try JSONDecoder().decode(ThermodynamicResult.self, from: encoded)

        XCTAssertEqual(result, decoded)
    }

    func testExplicitProxyProvenanceRoundTripsAcrossCBridge() async {
        let runner = FlexAIDRunner(
            temperature: 300.0,
            scientificProvenance: .proxyContactFunction
        )
        await runner.addSample(energy: -5.0)

        let result = await runner.compute()

        XCTAssertEqual(result.scientificProvenance, .proxyContactFunction)
        XCTAssertEqual(result.claimValidity, .proxyOnly)
    }

    func testExplicitPhysicalProvenanceRoundTripsAcrossCBridge() async {
        let provenance = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .matchedAssociationCycle,
            energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
            measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4",
            referenceProvenance: "sha256:01d26f66e709d388d7b971de6204680694e7cca10b1570506a5738c6909a6442"
        )
        let runner = FlexAIDRunner(
            temperature: 300.0,
            scientificProvenance: provenance
        )
        await runner.addSample(energy: -5.0)

        let result = await runner.compute()

        XCTAssertEqual(result.scientificProvenance, provenance)
        XCTAssertEqual(result.claimValidity, .bindingPhysical)
    }

    func testOutOfRangeSchemaFailsClosedAcrossCBridge() async {
        let provenance = ScientificProvenance(
            schemaVersion: Int.max,
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
            measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4"
        )
        let runner = FlexAIDRunner(scientificProvenance: provenance)
        await runner.addSample(energy: -5.0)

        let result = await runner.compute()

        XCTAssertEqual(result.scientificProvenance?.schemaVersion, 0)
        XCTAssertEqual(result.claimValidity, .proxyOnly)
    }

    func testEmbeddedNULCannotAuthorizeAcrossCBridge() async {
        let provenance = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .matchedAssociationCycle,
            energyProvenance: "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c\0junk",
            measureProvenance: "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4",
            referenceProvenance: "sha256:01d26f66e709d388d7b971de6204680694e7cca10b1570506a5738c6909a6442"
        )
        let runner = FlexAIDRunner(scientificProvenance: provenance)
        await runner.addSample(energy: -5.0)

        let result = await runner.compute()

        XCTAssertEqual(result.scientificProvenance?.energyProvenance, "")
        XCTAssertEqual(result.claimValidity, .proxyOnly)
    }
}
