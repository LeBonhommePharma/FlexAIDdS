// ShannonEntropyDecomposition.swift — shared ShannonThermoStack value model
// SPDX-License-Identifier: Apache-2.0

import Foundation

/// Decomposed entropy report from ShannonThermoStack.
///
/// This lives in the core Swift module because DockingRunner produces it;
/// HealthIntegration and Intelligence consume the same value without a
/// dependency cycle.
public struct ShannonEntropyDecomposition: Sendable, Codable, Hashable {
    /// Evidence carried across the C bridge from the ensemble that produced
    /// these values. Missing metadata fails closed to proxy-only.
    public let scientificProvenance: ScientificProvenance?

    /// Strongest claim supported by the declared provenance.
    public var claimValidity: ClaimValidity {
        scientificProvenance?.claimValidity ?? .proxyOnly
    }

    /// Configurational entropy from GA ensemble histogram (nats)
    public let configurational: Double

    /// Torsional vibrational entropy in the declared energy domain per kelvin;
    /// physical kcal mol^-1 K^-1 requires calibrated provenance.
    public let vibrational: Double

    /// Combined -T*S contribution in the declared energy domain.
    public let entropyContribution: Double

    /// Whether the Shannon entropy has reached a convergence plateau
    public let isConverged: Bool

    /// Relative change in entropy over the last convergence window
    public let convergenceRate: Double

    /// Hardware backend used for computation
    public let hardwareBackend: String

    /// Number of non-zero histogram bins
    public let occupiedBins: Int

    /// Total histogram bins used
    public let totalBins: Int

    /// Per-binding-mode Shannon entropy breakdown (nats)
    public let perModeEntropy: [Double]

    /// Top populated bins represented as (center, probability) pairs
    public let dominantBins: [(center: Double, probability: Double)]

    private enum CodingKeys: String, CodingKey {
        case configurational, vibrational, entropyContribution, isConverged
        case convergenceRate, hardwareBackend, occupiedBins, totalBins
        case perModeEntropy, dominantBinCenters, dominantBinProbabilities
        case scientificProvenance
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        // Absent or hostile provenance fails closed to proxy-only.
        scientificProvenance = try? container.decodeIfPresent(
            ScientificProvenance.self, forKey: .scientificProvenance)
        configurational = try container.decode(Double.self, forKey: .configurational)
        vibrational = try container.decode(Double.self, forKey: .vibrational)
        entropyContribution = try container.decode(Double.self, forKey: .entropyContribution)
        isConverged = try container.decode(Bool.self, forKey: .isConverged)
        convergenceRate = try container.decode(Double.self, forKey: .convergenceRate)
        hardwareBackend = try container.decode(String.self, forKey: .hardwareBackend)
        occupiedBins = try container.decode(Int.self, forKey: .occupiedBins)
        totalBins = try container.decode(Int.self, forKey: .totalBins)
        perModeEntropy = try container.decode([Double].self, forKey: .perModeEntropy)
        let centers = try container.decode([Double].self, forKey: .dominantBinCenters)
        let probabilities = try container.decode([Double].self, forKey: .dominantBinProbabilities)
        dominantBins = zip(centers, probabilities).map { ($0, $1) }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encodeIfPresent(scientificProvenance, forKey: .scientificProvenance)
        try container.encode(configurational, forKey: .configurational)
        try container.encode(vibrational, forKey: .vibrational)
        try container.encode(entropyContribution, forKey: .entropyContribution)
        try container.encode(isConverged, forKey: .isConverged)
        try container.encode(convergenceRate, forKey: .convergenceRate)
        try container.encode(hardwareBackend, forKey: .hardwareBackend)
        try container.encode(occupiedBins, forKey: .occupiedBins)
        try container.encode(totalBins, forKey: .totalBins)
        try container.encode(perModeEntropy, forKey: .perModeEntropy)
        try container.encode(dominantBins.map(\.center), forKey: .dominantBinCenters)
        try container.encode(dominantBins.map(\.probability), forKey: .dominantBinProbabilities)
    }

    public static func == (
        lhs: ShannonEntropyDecomposition,
        rhs: ShannonEntropyDecomposition
    ) -> Bool {
        lhs.configurational == rhs.configurational
            && lhs.vibrational == rhs.vibrational
            && lhs.entropyContribution == rhs.entropyContribution
            && lhs.isConverged == rhs.isConverged
            && lhs.hardwareBackend == rhs.hardwareBackend
            && lhs.occupiedBins == rhs.occupiedBins
            && lhs.perModeEntropy == rhs.perModeEntropy
            && lhs.scientificProvenance == rhs.scientificProvenance
    }

    public func hash(into hasher: inout Hasher) {
        hasher.combine(configurational)
        hasher.combine(vibrational)
        hasher.combine(entropyContribution)
        hasher.combine(isConverged)
        hasher.combine(hardwareBackend)
        hasher.combine(scientificProvenance)
    }

    public init(
        configurational: Double,
        vibrational: Double,
        entropyContribution: Double,
        isConverged: Bool,
        convergenceRate: Double,
        hardwareBackend: String,
        occupiedBins: Int,
        totalBins: Int,
        perModeEntropy: [Double] = [],
        dominantBins: [(center: Double, probability: Double)] = [],
        scientificProvenance: ScientificProvenance? = nil
    ) {
        self.scientificProvenance = scientificProvenance
        self.configurational = configurational
        self.vibrational = vibrational
        self.entropyContribution = entropyContribution
        self.isConverged = isConverged
        self.convergenceRate = convergenceRate
        self.hardwareBackend = hardwareBackend
        self.occupiedBins = occupiedBins
        self.totalBins = totalBins
        self.perModeEntropy = perModeEntropy
        self.dominantBins = dominantBins
    }
}
