// FleetAggregator.swift — Merge distributed ChunkResult populations into a unified ensemble
//
// Collects poses from all fleet devices, deduplicates by RMSD, re-clusters
// via density-peak heuristic, and recomputes global partition function /
// Boltzmann weights / thermodynamics across the merged population.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Aggregated Pose

/// A pose extracted from a ChunkResult for cross-device merging.
public struct AggregatedPose: Sendable, Codable, Hashable, Identifiable {
    public let id: UUID
    /// Energy-like scalar in the declared energy domain of the aggregation
    /// (see `FleetAggregator.scientificProvenance`). Fleet chunks carry raw
    /// optimizer scores, so this is not kcal/mol unless provenance says so.
    public let energy: Double
    /// Source chunk ID
    public let sourceChunkID: UUID
    /// Source device ID
    public let sourceDeviceID: String
    /// Cluster assignment after re-clustering (-1 = noise)
    public var clusterID: Int
    /// Boltzmann weight in the merged ensemble
    public var boltzmannWeight: Double

    public init(
        energy: Double, sourceChunkID: UUID,
        sourceDeviceID: String, clusterID: Int = -1
    ) {
        self.id = UUID()
        self.energy = energy
        self.sourceChunkID = sourceChunkID
        self.sourceDeviceID = sourceDeviceID
        self.clusterID = clusterID
        self.boltzmannWeight = 0.0
    }
}

// MARK: - Aggregated Binding Mode

/// A binding mode (cluster of poses) after fleet-wide re-clustering.
public struct AggregatedBindingMode: Sendable, Codable, Identifiable {
    public let id: Int
    /// Poses in this binding mode
    public let poses: [AggregatedPose]
    /// Number of contributing devices
    public let deviceCount: Int
    /// Mode-level thermodynamics
    public let thermodynamics: ModeThermodynamics

    /// Per-mode thermodynamic summary.
    ///
    /// These are the same identities the C++ engine evaluates, so they are
    /// bound by the same claim firewall: the arithmetic is valid for any
    /// finite scalar list, but physical units require provenance.
    public struct ModeThermodynamics: Sendable, Codable {
        /// Evidence for the values below. Missing metadata fails closed.
        public let scientificProvenance: ScientificProvenance?
        /// F = -kT ln Z_mode in the declared energy domain.
        public let freeEnergy: Double
        /// S = (<E> - F) / T in the declared energy domain per kelvin.
        public let entropy: Double
        /// Boltzmann-weighted mean energy <E> in the declared energy domain.
        public let meanEnergy: Double
        /// Boltzmann weight of this mode in the global ensemble
        public let globalWeight: Double

        /// Strongest claim supported by the declared provenance.
        public var claimValidity: ClaimValidity {
            scientificProvenance?.claimValidity ?? .proxyOnly
        }
    }
}

// MARK: - Aggregation Result

/// Complete result of fleet-wide aggregation.
public struct FleetAggregationResult: Sendable, Codable {
    /// All merged poses (with Boltzmann weights and cluster assignments)
    public let poses: [AggregatedPose]
    /// Re-clustered binding modes sorted by free energy
    public let bindingModes: [AggregatedBindingMode]
    /// Global thermodynamics across the merged ensemble
    public let globalThermodynamics: GlobalThermodynamics
    /// Number of poses removed as duplicates
    public let deduplicatedCount: Int
    /// Number of contributing chunks
    public let chunkCount: Int
    /// Source job ID
    public let jobID: UUID

    /// Global thermodynamic summary.
    ///
    /// Mirrors `FlexAIDdS.ThermodynamicResult`: numeric fields are always
    /// computed, physical interpretation is gated on provenance.
    public struct GlobalThermodynamics: Sendable, Codable {
        /// Evidence for the values below. Missing metadata fails closed.
        public let scientificProvenance: ScientificProvenance?
        public let temperature: Double
        public let logZ: Double
        /// Declared-domain values; kcal/mol only for calibrated provenance.
        public let freeEnergy: Double
        public let meanEnergy: Double
        public let entropy: Double
        public let heatCapacity: Double
        public let stdEnergy: Double

        /// Strongest claim supported by the declared provenance.
        public var claimValidity: ClaimValidity {
            scientificProvenance?.claimValidity ?? .proxyOnly
        }

        public var allowsCanonicalClaims: Bool {
            claimValidity == .canonicalPhysical || claimValidity == .bindingPhysical
        }

        public var allowsBindingClaims: Bool {
            claimValidity == .bindingPhysical
        }
    }
}

// MARK: - Fleet Aggregator

/// Merges distributed ChunkResult populations into a unified ensemble.
///
/// Usage:
/// ```swift
/// let aggregator = FleetAggregator(temperature: 298.15)
/// for result in chunkResults {
///     aggregator.ingest(result)
/// }
/// let merged = aggregator.aggregate(
///     deduplicationThreshold: 0.5,
///     clusterMinSize: 3
/// )
/// ```
public final class FleetAggregator: Sendable {
    /// Temperature (K) for Boltzmann weighting
    public let temperature: Double

    /// Evidence describing the energy domain and ensemble measure of the
    /// ingested chunk energies. Fleet chunks are optimizer output, so the
    /// default is `nil` — proxy-only — and a caller must supply the real
    /// record from the producing engine to unlock physical presentation.
    /// This aggregator never manufactures provenance for its own output.
    public let scientificProvenance: ScientificProvenance?

    /// Boltzmann constant (kcal mol^-1 K^-1)
    private let kB: Double = 0.001987204

    private let _poses: LockedBox<[AggregatedPose]> = .init([])
    private let _chunkIDs: LockedBox<Set<UUID>> = .init([])

    public init(
        temperature: Double = 298.15,
        scientificProvenance: ScientificProvenance? = nil
    ) {
        self.temperature = temperature
        self.scientificProvenance = scientificProvenance
    }

    /// Number of ingested poses.
    public var poseCount: Int { _poses.value.count }

    /// Number of ingested chunks.
    public var chunkCount: Int { _chunkIDs.value.count }

    // MARK: - Ingestion

    /// Extract poses from a ChunkResult and add them to the pool.
    ///
    /// The `resultData` field is decoded as a JSON array of pose energies.
    /// For production, this should decode full PoseResult structures.
    public func ingest(_ result: ChunkResult) {
        guard !_chunkIDs.value.contains(result.chunkID) else { return } // skip duplicate chunks
        _chunkIDs.mutate { $0.insert(result.chunkID) }

        // Decode pose energies from resultData
        let energies: [Double]
        if let decoded = try? JSONDecoder().decode([Double].self, from: result.resultData) {
            energies = decoded
        } else {
            // Fallback: use summary data to synthesize a representative pose
            energies = [result.summaryMeanEnergy]
        }

        let newPoses = energies.map { energy in
            AggregatedPose(
                energy: energy,
                sourceChunkID: result.chunkID,
                sourceDeviceID: result.computedBy
            )
        }

        _poses.mutate { $0.append(contentsOf: newPoses) }
    }

    /// Directly add poses (for testing or when poses are already extracted).
    public func addPoses(_ poses: [AggregatedPose]) {
        _poses.mutate { $0.append(contentsOf: poses) }
    }

    // MARK: - Aggregation Pipeline

    /// Run the full aggregation pipeline:
    /// 1. Deduplicate poses by energy proximity
    /// 2. Cluster via density-peak heuristic
    /// 3. Compute global and per-mode thermodynamics
    ///
    /// - Parameters:
    ///   - deduplicationThreshold: Energy difference (kcal/mol) below which two poses
    ///     from the same device are considered duplicates (default 0.01)
    ///   - clusterMinSize: Minimum poses to form a binding mode (default 2)
    ///   - clusterEnergyBandwidth: Energy bandwidth for density-peak clustering (default 2.0 kcal/mol)
    /// - Returns: Aggregated result with merged binding modes and thermodynamics
    public func aggregate(
        deduplicationThreshold: Double = 0.01,
        clusterMinSize: Int = 2,
        clusterEnergyBandwidth: Double = 2.0
    ) -> FleetAggregationResult {
        let jobID = _chunkIDs.value.first ?? UUID()

        // Step 1: Deduplicate
        let (deduplicated, removedCount) = deduplicatePoses(
            _poses.value,
            threshold: deduplicationThreshold
        )

        // Step 2: Cluster
        var clustered = clusterByEnergyDensity(
            deduplicated,
            bandwidth: clusterEnergyBandwidth,
            minSize: clusterMinSize
        )

        // Step 3: Compute Boltzmann weights across merged ensemble
        let beta = 1.0 / (kB * temperature)
        let logZ = logSumExp(clustered.map { -beta * $0.energy })

        for i in clustered.indices {
            clustered[i].boltzmannWeight = exp(-beta * clustered[i].energy - logZ)
        }

        // Step 4: Build binding modes with per-mode thermodynamics
        let modeIDs = Set(clustered.map(\.clusterID)).sorted()
        var bindingModes: [AggregatedBindingMode] = []

        for modeID in modeIDs where modeID >= 0 {
            let modePoses = clustered.filter { $0.clusterID == modeID }
            guard modePoses.count >= clusterMinSize else { continue }

            let modeThermo = computeModeThermodynamics(
                modePoses, beta: beta, globalLogZ: logZ
            )
            let deviceCount = Set(modePoses.map(\.sourceDeviceID)).count

            bindingModes.append(AggregatedBindingMode(
                id: modeID,
                poses: modePoses,
                deviceCount: deviceCount,
                thermodynamics: modeThermo
            ))
        }

        // Sort by free energy (most favorable first)
        bindingModes.sort { $0.thermodynamics.freeEnergy < $1.thermodynamics.freeEnergy }

        // Step 5: Global thermodynamics
        let globalThermo = computeGlobalThermodynamics(clustered, beta: beta, logZ: logZ)

        return FleetAggregationResult(
            poses: clustered,
            bindingModes: bindingModes,
            globalThermodynamics: globalThermo,
            deduplicatedCount: removedCount,
            chunkCount: _chunkIDs.value.count,
            jobID: jobID
        )
    }

    // MARK: - Deduplication

    /// Remove near-duplicate poses from the same device (retry/resubmission artifacts).
    ///
    /// Two poses from the same device with |E_i - E_j| < threshold are duplicates;
    /// only the first is kept.
    func deduplicatePoses(
        _ poses: [AggregatedPose],
        threshold: Double
    ) -> (poses: [AggregatedPose], removedCount: Int) {
        guard !poses.isEmpty else { return ([], 0) }

        // Group by source device
        var byDevice: [String: [AggregatedPose]] = [:]
        for pose in poses {
            byDevice[pose.sourceDeviceID, default: []].append(pose)
        }

        var deduplicated: [AggregatedPose] = []
        var removedCount = 0

        for (_, devicePoses) in byDevice {
            // Sort by energy for efficient dedup
            let sorted = devicePoses.sorted { $0.energy < $1.energy }
            var kept: [AggregatedPose] = []

            for pose in sorted {
                let isDuplicate = kept.contains { abs($0.energy - pose.energy) < threshold }
                if isDuplicate {
                    removedCount += 1
                } else {
                    kept.append(pose)
                }
            }
            deduplicated.append(contentsOf: kept)
        }

        return (deduplicated, removedCount)
    }

    // MARK: - Density-Peak Clustering

    /// Simple energy-based density-peak clustering.
    ///
    /// Finds local density maxima in the energy histogram, then assigns
    /// each pose to the nearest density peak within the bandwidth.
    func clusterByEnergyDensity(
        _ poses: [AggregatedPose],
        bandwidth: Double,
        minSize: Int
    ) -> [AggregatedPose] {
        guard poses.count >= minSize else {
            return poses.map { var p = $0; p.clusterID = -1; return p }
        }

        let sorted = poses.sorted { $0.energy < $1.energy }
        let energies = sorted.map(\.energy)

        // Compute local density for each pose (Gaussian kernel)
        var densities = [Double](repeating: 0, count: energies.count)
        for i in energies.indices {
            for j in energies.indices where i != j {
                let d = abs(energies[i] - energies[j]) / bandwidth
                densities[i] += exp(-0.5 * d * d)
            }
        }

        // Find density peaks (local maxima with density above median)
        let medianDensity = densities.sorted()[densities.count / 2]
        var peaks: [Int] = []
        for i in densities.indices {
            guard densities[i] > medianDensity else { continue }
            // A peak: higher density than all neighbors within bandwidth
            let isPeak = !densities.indices.contains { j in
                j != i
                && abs(energies[i] - energies[j]) < bandwidth
                && densities[j] > densities[i]
            }
            if isPeak { peaks.append(i) }
        }

        // Assign each pose to nearest peak within bandwidth, or noise (-1)
        var result = sorted
        for i in result.indices {
            var bestPeak = -1
            var bestDist = Double.infinity
            for peakIdx in peaks {
                let dist = abs(energies[i] - energies[peakIdx])
                if dist < bandwidth * 2.0 && dist < bestDist {
                    bestDist = dist
                    bestPeak = peaks.firstIndex(of: peakIdx)!
                }
            }
            result[i].clusterID = bestPeak
        }

        // Prune clusters below minSize → noise
        let clusterSizes = Dictionary(result.map { ($0.clusterID, 1) }, uniquingKeysWith: +)
        for i in result.indices {
            if let size = clusterSizes[result[i].clusterID], size < minSize {
                result[i].clusterID = -1
            }
        }

        return result
    }

    // MARK: - Thermodynamics

    /// Compute per-mode thermodynamics.
    private func computeModeThermodynamics(
        _ poses: [AggregatedPose], beta: Double, globalLogZ: Double
    ) -> AggregatedBindingMode.ModeThermodynamics {
        let modeLogZ = logSumExp(poses.map { -beta * $0.energy })
        let freeEnergy = -modeLogZ / beta

        // <E>_mode
        var meanE = 0.0
        for pose in poses {
            let w = exp(-beta * pose.energy - modeLogZ)
            meanE += w * pose.energy
        }

        let entropy = (meanE - freeEnergy) / temperature

        // Mode weight in global ensemble: Z_mode / Z_global
        let globalWeight = exp(modeLogZ - globalLogZ)

        return AggregatedBindingMode.ModeThermodynamics(
            scientificProvenance: scientificProvenance,
            freeEnergy: freeEnergy,
            entropy: entropy,
            meanEnergy: meanE,
            globalWeight: globalWeight
        )
    }

    /// Compute global ensemble thermodynamics.
    private func computeGlobalThermodynamics(
        _ poses: [AggregatedPose], beta: Double, logZ: Double
    ) -> FleetAggregationResult.GlobalThermodynamics {
        let freeEnergy = -logZ / beta

        // <E> and <E^2>
        var meanE = 0.0
        var meanESq = 0.0
        for pose in poses {
            let w = exp(-beta * pose.energy - logZ)
            meanE += w * pose.energy
            meanESq += w * pose.energy * pose.energy
        }

        let entropy = (meanE - freeEnergy) / temperature
        let variance = meanESq - meanE * meanE
        let heatCapacity = variance * beta * beta
        let stdEnergy = sqrt(max(0, variance))

        return FleetAggregationResult.GlobalThermodynamics(
            scientificProvenance: scientificProvenance,
            temperature: temperature,
            logZ: logZ,
            freeEnergy: freeEnergy,
            meanEnergy: meanE,
            entropy: entropy,
            heatCapacity: heatCapacity,
            stdEnergy: stdEnergy
        )
    }

    // MARK: - Numerics

    /// Log-sum-exp for numerical stability: log(sum(exp(x_i)))
    private func logSumExp(_ values: [Double]) -> Double {
        guard let maxVal = values.max() else { return -.infinity }
        let sumExp = values.reduce(0.0) { $0 + exp($1 - maxVal) }
        return maxVal + log(sumExp)
    }
}

// MARK: - Thread-Safe Box

/// Minimal thread-safe container for Sendable conformance.
final class LockedBox<T: Sendable>: @unchecked Sendable {
    private var _value: T
    private let lock = NSLock()

    init(_ value: T) { _value = value }

    var value: T {
        lock.lock()
        defer { lock.unlock() }
        return _value
    }

    func mutate(_ transform: (inout T) -> Void) {
        lock.lock()
        defer { lock.unlock() }
        transform(&_value)
    }
}
