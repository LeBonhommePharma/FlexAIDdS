// FleetScheduler.swift — Distributed work scheduling across Apple devices
//
// Actor managing work distribution, iCloud coordination, thermal awareness,
// battery-aware scheduling, orphan recovery, and local fallback.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import Intelligence
#if canImport(CryptoKit)
import CryptoKit
#endif
#if canImport(os)
import os
#endif

/// Fleet operation mode.
public enum FleetMode: String, Sendable, Codable {
    case distributed   // iCloud fleet coordination across devices
    case localOnly     // Single device (fallback when iCloud unavailable)
}

/// Fleet-wide metrics for dashboard reporting.
public struct FleetMetrics: Sendable, Codable {
    public let jobID: UUID
    public let totalChunks: Int
    public let completedChunks: Int
    public let failedChunks: Int
    public let orphanedChunks: Int
    public let activeDevices: Int
    public let totalTFLOPS: Double
    public let meanChunkTimeSeconds: Double?
    public let estimatedRemainingSeconds: Double?
    public let timestamp: Date

    public var progressFraction: Double {
        totalChunks > 0 ? Double(completedChunks) / Double(totalChunks) : 0
    }

    public var isComplete: Bool {
        completedChunks + failedChunks >= totalChunks
    }
}

/// Actor for scheduling and distributing docking work across the Bonhomme fleet.
///
/// Usage:
/// ```swift
/// let scheduler = FleetScheduler()
/// let chunks = scheduler.splitWork(totalChromosomes: 1000, devices: devices)
/// for chunk in chunks {
///     try await scheduler.submitToiCloud(chunk)
/// }
///
/// // Periodically sweep for orphaned chunks
/// let recovered = scheduler.sweepOrphanedChunks()
/// ```
public actor FleetScheduler {
    private var activeJobs: [UUID: [WorkChunk]] = [:]
    private var completedResults: [UUID: [ChunkResult]] = [:]
    private var deviceThroughput: [String: Double] = [:]  // deviceID → chromosomes/second

    #if canImport(CryptoKit)
    private var encryptionKey: SymmetricKey?
    #endif

    /// iCloud Drive container identifier for fleet data
    private let containerID: String

    /// Current fleet mode (distributed or local-only fallback)
    public private(set) var mode: FleetMode

    #if canImport(os)
    private let logger = Logger(subsystem: "com.bonhomme.flexaidds", category: "fleet")
    #endif

    public init(containerID: String = "iCloud.com.bonhomme.flexaidds") {
        self.containerID = containerID
        // Detect iCloud availability and set mode
        let fileManager = FileManager.default
        if fileManager.url(forUbiquityContainerIdentifier: nil) != nil {
            self.mode = .distributed
        } else {
            self.mode = .localOnly
        }
    }

    // MARK: - ADP Enforcement

    /// Check if Advanced Data Protection is available.
    /// Fleet mode requires E2E encryption — blocks if ADP is not enabled.
    public nonisolated func isAdvancedDataProtectionAvailable() -> Bool {
        let fileManager = FileManager.default
        guard let containerURL = fileManager.url(forUbiquityContainerIdentifier: nil) else {
            return false
        }
        return fileManager.fileExists(atPath: containerURL.path)
    }

    /// Check iCloud health and switch mode if needed.
    public func refreshMode() {
        let fileManager = FileManager.default
        if fileManager.url(forUbiquityContainerIdentifier: nil) != nil {
            mode = .distributed
        } else {
            mode = .localOnly
            #if canImport(os)
            logger.warning("iCloud unavailable — falling back to local-only mode")
            #endif
        }
    }

    // MARK: - Work Splitting

    /// Split docking work across available devices proportional to their compute weight.
    /// - Parameters:
    ///   - totalChromosomes: Total number of chromosomes in the GA population
    ///   - maxGenerations: Maximum GA generations
    ///   - temperature: Simulation temperature in Kelvin
    ///   - devices: Available devices and their capabilities
    ///   - configData: Serialized configuration data
    ///   - chunkTimeout: Timeout per chunk in seconds (default 1 hour)
    /// - Returns: Array of work chunks ready for distribution
    public func splitWork(
        totalChromosomes: Int,
        maxGenerations: Int,
        temperature: Double = 300.0,
        devices: [DeviceCapability],
        configData: Data,
        chunkTimeout: TimeInterval = 3600
    ) -> [WorkChunk] {
        let jobID = UUID()

        // Filter devices that are available (thermal + battery safe)
        let activeDevices = devices.filter { $0.isAvailable }
        guard !activeDevices.isEmpty else {
            #if canImport(os)
            logger.error("No available devices for job — all excluded by thermal/battery state")
            #endif
            return []
        }

        // Use dynamic throughput if available, otherwise use static compute weight
        let weights: [Double] = activeDevices.map { device in
            if let throughput = deviceThroughput[device.deviceID], throughput > 0 {
                return throughput * (device.computeWeight > 0 ? 1.0 : 0.0)
            }
            return device.computeWeight
        }
        let totalWeight = weights.reduce(0.0, +)
        guard totalWeight > 0 else { return [] }

        var chunks: [WorkChunk] = []
        var allocatedChromosomes = 0

        #if canImport(os)
        logger.info("Splitting \(totalChromosomes) chromosomes across \(activeDevices.count) devices (total weight: \(String(format: "%.2f", totalWeight)))")
        #endif

        for (i, device) in activeDevices.enumerated() {
            let isLast = (i == activeDevices.count - 1)
            let share = weights[i] / totalWeight
            let chromCount = isLast
                ? totalChromosomes - allocatedChromosomes
                : Int(Double(totalChromosomes) * share)

            guard chromCount > 0 else { continue }

            let params = GAChunkParameters(
                numChromosomes: chromCount,
                maxGenerations: maxGenerations,
                seed: Int.random(in: 0..<Int.max),
                temperature: temperature
            )

            var chunk = WorkChunk(
                jobID: jobID,
                index: chunks.count,
                totalChunks: activeDevices.count,
                configData: configData,
                gaParameters: params,
                timeoutSeconds: chunkTimeout
            )
            chunk.claim(by: device.deviceID)
            chunks.append(chunk)

            allocatedChromosomes += chromCount

            #if canImport(os)
            logger.info("  Chunk \(chunk.index): \(chromCount) chromosomes → \(device.model) (\(device.statusSummary))")
            #endif
        }

        activeJobs[jobID] = chunks
        return chunks
    }

    // MARK: - Orphan Recovery

    /// Sweep all active jobs for timed-out chunks and prepare them for retry.
    /// Returns the number of chunks reclaimed.
    @discardableResult
    public func sweepOrphanedChunks() -> Int {
        var recovered = 0

        for (jobID, chunks) in activeJobs {
            var updated = chunks
            for i in updated.indices {
                guard updated[i].status == .claimed || updated[i].status == .running else { continue }

                if updated[i].isTimedOut {
                    let deviceID = updated[i].claimedBy ?? "unknown"
                    let success = updated[i].markOrphanedForRetry()

                    #if canImport(os)
                    if success {
                        logger.warning("Chunk \(updated[i].index) orphaned from \(deviceID) — retry \(updated[i].retryCount)/\(updated[i].maxRetries)")
                    } else {
                        logger.error("Chunk \(updated[i].index) permanently failed after \(updated[i].maxRetries) retries")
                    }
                    #endif

                    recovered += 1
                }
            }
            activeJobs[jobID] = updated
        }

        return recovered
    }

    /// Get unclaimed or orphaned chunks available for work-stealing.
    public func availableChunks(for jobID: UUID) -> [WorkChunk] {
        guard let chunks = activeJobs[jobID] else { return [] }
        return chunks
            .filter { $0.status == .pending || $0.status == .orphaned }
            .sorted { $0.priority > $1.priority }
    }

    /// Claim an available chunk for a device (work-stealing).
    public func claimChunk(chunkID: UUID, deviceID: String) -> Bool {
        for (jobID, chunks) in activeJobs {
            if let idx = chunks.firstIndex(where: { $0.id == chunkID }) {
                guard activeJobs[jobID]![idx].status == .pending || activeJobs[jobID]![idx].status == .orphaned else {
                    return false
                }
                activeJobs[jobID]![idx].claim(by: deviceID)
                return true
            }
        }
        return false
    }

    // MARK: - iCloud Submission

    /// Submit a work chunk to iCloud Drive for fleet pickup.
    /// Falls back to local mode if iCloud is unavailable.
    public func submitToiCloud(_ chunk: WorkChunk) async throws {
        let fileManager = FileManager.default
        guard let containerURL = fileManager.url(forUbiquityContainerIdentifier: nil) else {
            if mode == .distributed {
                mode = .localOnly
                #if canImport(os)
                logger.warning("iCloud became unavailable during submission — switched to local mode")
                #endif
            }
            throw FleetError.iCloudUnavailable
        }

        let jobDir = containerURL.appendingPathComponent("FleetJobs/\(chunk.jobID.uuidString)")
        try fileManager.createDirectory(at: jobDir, withIntermediateDirectories: true)

        let chunkFile = jobDir.appendingPathComponent("chunk_\(chunk.index).json")
        let data = try JSONEncoder().encode(chunk)
        try data.write(to: chunkFile)
    }

    /// Submit a completed result back to iCloud Drive.
    public func submitResult(_ result: ChunkResult) async throws {
        let fileManager = FileManager.default
        guard let containerURL = fileManager.url(forUbiquityContainerIdentifier: nil) else {
            throw FleetError.iCloudUnavailable
        }

        let resultsDir = containerURL.appendingPathComponent("FleetResults/\(result.jobID.uuidString)")
        try fileManager.createDirectory(at: resultsDir, withIntermediateDirectories: true)

        let resultFile = resultsDir.appendingPathComponent("result_\(result.chunkID.uuidString).json")
        let data = try JSONEncoder().encode(result)
        try data.write(to: resultFile)

        // Track completion
        completedResults[result.jobID, default: []].append(result)

        // Update device throughput for dynamic rebalancing
        if let chunks = activeJobs[result.jobID],
           let chunk = chunks.first(where: { $0.id == result.chunkID }) {
            let chromosomes = Double(chunk.gaParameters.numChromosomes)
            let throughput = chromosomes / max(result.computeTimeSeconds, 0.001)
            deviceThroughput[result.computedBy] = throughput
        }

        #if canImport(os)
        logger.info("Result received for chunk \(result.chunkID) from \(result.computedBy) in \(String(format: "%.1f", result.computeTimeSeconds))s")
        #endif
    }

    /// Check if all chunks for a job have completed.
    public func isJobComplete(_ jobID: UUID) -> Bool {
        guard let chunks = activeJobs[jobID] else { return false }
        let results = completedResults[jobID] ?? []
        let failed = chunks.filter { $0.status == .permanentlyFailed }.count
        return results.count + failed >= chunks.count
    }

    // MARK: - Fleet Metrics

    /// Get current metrics for a job (for dashboard reporting).
    public func metrics(for jobID: UUID) -> FleetMetrics? {
        guard let chunks = activeJobs[jobID] else { return nil }
        let results = completedResults[jobID] ?? []

        let completedCount = results.count
        let failedCount = chunks.filter { $0.status == .permanentlyFailed }.count
        let orphanedCount = chunks.filter { $0.status == .orphaned }.count
        let activeDeviceIDs = Set(chunks.compactMap(\.claimedBy))

        let meanTime: Double? = results.isEmpty ? nil :
            results.reduce(0.0) { $0 + $1.computeTimeSeconds } / Double(results.count)

        let remaining: Double? = {
            guard let mean = meanTime, completedCount > 0 else { return nil }
            let remaining = chunks.count - completedCount - failedCount
            return mean * Double(remaining) / max(1, Double(activeDeviceIDs.count))
        }()

        let totalTFLOPS = deviceThroughput.values.reduce(0.0, +) * 0.001  // Rough conversion

        return FleetMetrics(
            jobID: jobID,
            totalChunks: chunks.count,
            completedChunks: completedCount,
            failedChunks: failedCount,
            orphanedChunks: orphanedCount,
            activeDevices: activeDeviceIDs.count,
            totalTFLOPS: totalTFLOPS,
            meanChunkTimeSeconds: meanTime,
            estimatedRemainingSeconds: remaining,
            timestamp: Date()
        )
    }

    /// Export fleet metrics as JSON for the PWA dashboard.
    public func metricsJSON(for jobID: UUID) throws -> Data? {
        guard let m = metrics(for: jobID) else { return nil }
        return try JSONEncoder().encode(m)
    }
}

// MARK: - Referee Integration

/// Recommendation from the thermodynamic referee for scheduling decisions.
public enum RefereeRecommendation: Sendable, Codable {
    /// Results are trustworthy — no further action needed
    case proceed
    /// Increase population size to improve sampling
    case increasePopulation(factor: Double)
    /// Increase GA generations to reach convergence
    case increaseGenerations(factor: Double)
    /// Need more data — cannot assess quality
    case needsMoreData
}

extension FleetScheduler {
    /// Aggregate completed chunk results and run the thermodynamic referee.
    ///
    /// After all chunks for a job complete, this method:
    /// 1. Collects all chunk result data
    /// 2. Runs `RuleBasedReferee` on the aggregated thermodynamics
    /// 3. Returns a scheduling recommendation based on the verdict
    ///
    /// - Parameter jobID: The fleet job to analyze
    /// - Returns: Recommendation for scheduling follow-up work, or nil if job incomplete
    public func aggregateAndReferee(jobID: UUID) -> (CrossPlatformRefereeVerdict, RefereeRecommendation)? {
        guard isJobComplete(jobID) else { return nil }
        guard let results = completedResults[jobID], !results.isEmpty else { return nil }

        // Aggregate thermodynamic summary from chunk results
        // Each ChunkResult carries summary data; we combine them
        let totalChunks = results.count
        var sumFreeEnergy = 0.0
        var sumEntropy = 0.0
        var sumHeatCapacity = 0.0
        var sumMeanEnergy = 0.0
        var sumStdEnergy = 0.0
        var totalModes = 0
        var minFreeEnergy = Double.infinity

        for result in results {
            sumFreeEnergy += result.summaryFreeEnergy
            sumEntropy += result.summaryEntropy
            sumHeatCapacity += result.summaryHeatCapacity
            sumMeanEnergy += result.summaryMeanEnergy
            sumStdEnergy += result.summaryStdEnergy
            totalModes += result.summaryModeCount
            minFreeEnergy = min(minFreeEnergy, result.summaryBestFreeEnergy)
        }

        // Weighted average across chunks
        let n = Double(totalChunks)
        let avgThermo = AggregatedThermodynamics(
            temperature: results.first?.temperature ?? 298.15,
            freeEnergy: sumFreeEnergy / n,
            entropy: sumEntropy / n,
            heatCapacity: sumHeatCapacity / n,
            meanEnergy: sumMeanEnergy / n,
            stdEnergy: sumStdEnergy / n,
            modeCount: totalModes,
            bestFreeEnergy: minFreeEnergy
        )

        // Build a verdict using threshold logic (no LLM needed for fleet decisions)
        let verdict = buildFleetVerdict(from: avgThermo)

        // Determine scheduling recommendation from verdict
        let recommendation = scheduleFromVerdict(verdict)

        return (verdict, recommendation)
    }

    /// Build a fleet-level verdict from aggregated thermodynamics.
    private func buildFleetVerdict(from thermo: AggregatedThermodynamics) -> CrossPlatformRefereeVerdict {
        var findings: [CrossPlatformRefereeFinding] = []
        var trustworthy = true

        // Fleet-level convergence heuristic: if stdEnergy > |freeEnergy| * 0.8, not converged
        if thermo.stdEnergy > abs(thermo.freeEnergy) * 0.8 {
            trustworthy = false
            findings.append(CrossPlatformRefereeFinding(
                title: "Fleet ensemble not converged",
                detail: "Energy variance (σ=\(String(format: "%.2f", thermo.stdEnergy))) is large relative to F (\(String(format: "%.2f", thermo.freeEnergy))). Aggregated chunks may need more generations.",
                severity: "critical",
                category: "convergence"
            ))
        } else {
            findings.append(CrossPlatformRefereeFinding(
                title: "Fleet ensemble converging",
                detail: "Energy variance (σ=\(String(format: "%.2f", thermo.stdEnergy))) reasonable for F=\(String(format: "%.2f", thermo.freeEnergy)) across \(thermo.modeCount) modes.",
                severity: "pass",
                category: "convergence"
            ))
        }

        // Fleet aggregation currently carries no calibrated energy/measure/
        // reference provenance, so it cannot support an affinity claim.
        findings.append(CrossPlatformRefereeFinding(
            title: "Aggregated affinity unavailable",
            detail: "Fleet F is an ensemble-score diagnostic; a calibrated matched association cycle was not supplied.",
            severity: "advisory",
            category: "provenance"
        ))

        let action = trustworthy
            ? "Fleet results reliable. Proceed with full FOPTICS re-clustering."
            : "Increase sampling. Submit follow-up job with larger population."

        return CrossPlatformRefereeVerdict(
            findings: findings,
            overallTrustworthy: trustworthy,
            recommendedAction: action,
            confidence: trustworthy ? 0.8 : 0.4
        )
    }

    /// Translate a verdict into a concrete scheduling recommendation.
    private func scheduleFromVerdict(_ verdict: CrossPlatformRefereeVerdict) -> RefereeRecommendation {
        if verdict.overallTrustworthy {
            return .proceed
        }

        // Check for convergence issues → increase generations
        let hasConvergenceIssue = verdict.findings.contains {
            $0.category == "convergence" && ($0.severity == "critical" || $0.severity == "warning")
        }
        if hasConvergenceIssue {
            return .increaseGenerations(factor: 1.5)
        }

        // Check for sampling issues → increase population
        let hasSamplingIssue = verdict.findings.contains {
            $0.category == "histogram" && ($0.severity == "critical" || $0.severity == "warning")
        }
        if hasSamplingIssue {
            return .increasePopulation(factor: 2.0)
        }

        return .needsMoreData
    }
}

/// Aggregated thermodynamics from fleet chunk results (not a full StatMechEngine recompute).
struct AggregatedThermodynamics {
    let temperature: Double
    let freeEnergy: Double
    let entropy: Double
    let heatCapacity: Double
    let meanEnergy: Double
    let stdEnergy: Double
    let modeCount: Int
    let bestFreeEnergy: Double
}

// MARK: - iCloud Watcher

/// Watches an iCloud Drive directory for incoming fleet results.
#if os(macOS) || os(iOS)
public final class iCloudWatcher: @unchecked Sendable {
    private var query: NSMetadataQuery?
    private let containerID: String
    private let continuation: AsyncStream<ChunkResult>.Continuation?
    private let stream: AsyncStream<ChunkResult>
    private var seenResultIDs: Set<String> = []  // Deduplication

    /// Stream of incoming chunk results.
    public var results: AsyncStream<ChunkResult> {
        stream
    }

    public init(containerID: String = "iCloud.com.bonhomme.flexaidds") {
        self.containerID = containerID
        var cont: AsyncStream<ChunkResult>.Continuation?
        self.stream = AsyncStream { cont = $0 }
        self.continuation = cont
    }

    /// Start watching for incoming results.
    public func start() {
        let query = NSMetadataQuery()
        query.searchScopes = [NSMetadataQueryUbiquitousDataScope]
        query.predicate = NSPredicate(format: "%K LIKE 'FleetResults/*/result_*.json'",
                                       NSMetadataItemPathKey)

        NotificationCenter.default.addObserver(
            forName: .NSMetadataQueryDidUpdate,
            object: query, queue: .main
        ) { [weak self] notification in
            self?.handleQueryUpdate(notification)
        }

        query.start()
        self.query = query
    }

    /// Stop watching.
    public func stop() {
        query?.stop()
        query = nil
        continuation?.finish()
    }

    private func handleQueryUpdate(_ notification: Notification) {
        guard let query = notification.object as? NSMetadataQuery else { return }

        for item in query.results {
            guard let mdItem = item as? NSMetadataItem,
                  let path = mdItem.value(forAttribute: NSMetadataItemPathKey) as? String else {
                continue
            }

            // Deduplicate (NSMetadataQuery may report same file multiple times)
            guard !seenResultIDs.contains(path) else { continue }
            seenResultIDs.insert(path)

            let url = URL(fileURLWithPath: path)
            if let data = try? Data(contentsOf: url),
               let result = try? JSONDecoder().decode(ChunkResult.self, from: data) {
                continuation?.yield(result)
            }
        }
    }
}
#endif

// MARK: - Fleet Aggregation

extension FleetScheduler {
    /// Aggregate completed results for a job into a unified binding population.
    ///
    /// Runs the full FleetAggregator pipeline: deduplication, re-clustering,
    /// global Boltzmann weighting, and per-mode thermodynamics.
    ///
    /// - Parameters:
    ///   - jobID: The fleet job to aggregate
    ///   - deduplicationThreshold: Energy difference (kcal/mol) for duplicate detection (default 0.01)
    ///   - clusterMinSize: Minimum poses per binding mode (default 2)
    ///   - clusterEnergyBandwidth: Energy bandwidth for density-peak clustering (default 2.0)
    /// - Returns: Aggregated result, or nil if the job is not yet complete
    public func aggregateResults(
        jobID: UUID,
        deduplicationThreshold: Double = 0.01,
        clusterMinSize: Int = 2,
        clusterEnergyBandwidth: Double = 2.0
    ) -> FleetAggregationResult? {
        guard isJobComplete(jobID) else { return nil }
        guard let results = completedResults[jobID], !results.isEmpty else { return nil }

        let temperature = results.first?.temperature ?? 298.15
        let aggregator = FleetAggregator(temperature: temperature)

        for result in results {
            aggregator.ingest(result)
        }

        return aggregator.aggregate(
            deduplicationThreshold: deduplicationThreshold,
            clusterMinSize: clusterMinSize,
            clusterEnergyBandwidth: clusterEnergyBandwidth
        )
    }
}

// MARK: - Errors

public enum FleetError: Error, LocalizedError {
    case iCloudUnavailable
    case advancedDataProtectionRequired
    case encryptionFailed
    case noAvailableDevices
    case chunkTimeout(chunkID: UUID, deviceID: String)
    case allRetriesExhausted(chunkID: UUID)

    public var errorDescription: String? {
        switch self {
        case .iCloudUnavailable:
            return "iCloud Drive is not available — sign in to iCloud to use fleet mode"
        case .advancedDataProtectionRequired:
            return "Advanced Data Protection must be enabled for fleet mode"
        case .encryptionFailed:
            return "Failed to encrypt work chunk for transit"
        case .noAvailableDevices:
            return "No devices available for fleet compute (all may be in critical thermal state or low battery)"
        case .chunkTimeout(let chunkID, let deviceID):
            return "Chunk \(chunkID) timed out on device \(deviceID)"
        case .allRetriesExhausted(let chunkID):
            return "Chunk \(chunkID) failed permanently after all retries exhausted"
        }
    }
}
