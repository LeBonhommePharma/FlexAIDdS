// FleetSchedulerTests.swift — Unit tests for the fleet scheduler
//
// Tests work splitting, encryption round-trips, and device capability detection.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import XCTest
@testable import FleetScheduler
@testable import FlexAIDdS
#if canImport(CryptoKit)
import CryptoKit
#endif

final class FleetSchedulerTests: XCTestCase {

    // MARK: - Work Splitting

    func testSplitWorkProportional() async {
        let scheduler = FleetScheduler()

        let devices = [
            DeviceCapability(
                deviceID: "mac", model: "MacBookPro18,3",
                estimatedTFLOPS: 5.0, hasGPU: true,
                availableMemoryGB: 16.0, thermalState: .nominal,
                computeWeight: 0.75),
            DeviceCapability(
                deviceID: "iphone", model: "iPhone15,3",
                estimatedTFLOPS: 1.5, hasGPU: true,
                availableMemoryGB: 6.0, thermalState: .nominal,
                computeWeight: 0.15),
            DeviceCapability(
                deviceID: "ipad", model: "iPad13,18",
                estimatedTFLOPS: 0.8, hasGPU: true,
                availableMemoryGB: 4.0, thermalState: .fair,
                computeWeight: 0.10),
        ]

        let chunks = await scheduler.splitWork(
            totalChromosomes: 1000,
            maxGenerations: 100,
            devices: devices,
            configData: Data()
        )

        XCTAssertEqual(chunks.count, 3)

        // Total chromosomes should sum to 1000
        let totalChrom = chunks.reduce(0) { $0 + $1.gaParameters.numChromosomes }
        XCTAssertEqual(totalChrom, 1000)

        // Mac should get the most chromosomes
        XCTAssertGreaterThan(chunks[0].gaParameters.numChromosomes, chunks[1].gaParameters.numChromosomes)
        XCTAssertGreaterThan(chunks[1].gaParameters.numChromosomes, chunks[2].gaParameters.numChromosomes)
    }

    func testCriticalThermalDevicesExcluded() async {
        let scheduler = FleetScheduler()

        let devices = [
            DeviceCapability(
                deviceID: "mac", model: "MacBookPro18,3",
                estimatedTFLOPS: 5.0, hasGPU: true,
                availableMemoryGB: 16.0, thermalState: .nominal,
                computeWeight: 0.75),
            DeviceCapability(
                deviceID: "hot-iphone", model: "iPhone15,3",
                estimatedTFLOPS: 1.5, hasGPU: true,
                availableMemoryGB: 6.0, thermalState: .critical,
                computeWeight: 0.0),
        ]

        let chunks = await scheduler.splitWork(
            totalChromosomes: 500,
            maxGenerations: 50,
            devices: devices,
            configData: Data()
        )

        // Only 1 chunk — critical device excluded
        XCTAssertEqual(chunks.count, 1)
        XCTAssertEqual(chunks[0].gaParameters.numChromosomes, 500)
    }

    // MARK: - Encryption

    #if canImport(CryptoKit)
    func testEncryptionRoundTrip() throws {
        let key = SymmetricKey(size: .bits256)

        let chunk = WorkChunk(
            jobID: UUID(),
            index: 0,
            totalChunks: 1,
            configData: "test config".data(using: .utf8)!,
            gaParameters: GAChunkParameters(numChromosomes: 100, maxGenerations: 50, seed: 42)
        )

        let encrypted = try chunk.encrypt(using: key)
        let decrypted = try WorkChunk.decrypt(from: encrypted, using: key)

        XCTAssertEqual(chunk.id, decrypted.id)
        XCTAssertEqual(chunk.jobID, decrypted.jobID)
        XCTAssertEqual(chunk.gaParameters.numChromosomes, decrypted.gaParameters.numChromosomes)
        XCTAssertEqual(chunk.gaParameters.seed, decrypted.gaParameters.seed)
    }
    #endif

    // MARK: - Device Capability

    func testCurrentDeviceDetection() {
        let device = DeviceCapability.current()

        XCTAssertFalse(device.deviceID.isEmpty)
        XCTAssertGreaterThan(device.availableMemoryGB, 0)
        XCTAssertGreaterThan(device.estimatedTFLOPS, 0)
    }

    // MARK: - Work Chunk Models

    func testWorkChunkCodable() throws {
        let chunk = WorkChunk(
            jobID: UUID(),
            index: 3,
            totalChunks: 10,
            configData: Data([1, 2, 3]),
            gaParameters: GAChunkParameters(numChromosomes: 200, maxGenerations: 100, seed: 99)
        )

        let data = try JSONEncoder().encode(chunk)
        let decoded = try JSONDecoder().decode(WorkChunk.self, from: data)

        XCTAssertEqual(chunk.id, decoded.id)
        XCTAssertEqual(chunk.index, decoded.index)
        XCTAssertEqual(chunk.totalChunks, decoded.totalChunks)
        XCTAssertEqual(chunk.gaParameters.seed, decoded.gaParameters.seed)
    }
}
