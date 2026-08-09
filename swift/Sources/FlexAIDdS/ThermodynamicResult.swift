// ThermodynamicResult.swift — Swift model for statistical mechanics thermodynamics
//
// Maps from the C FXThermodynamics struct to a Swift-native Sendable, Codable type.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import FlexAIDCore

/// Full thermodynamic analysis of a conformational ensemble.
///
/// Computed from a partition function Z(T) via log-sum-exp for numerical stability:
/// - Free energy: F = -kT ln Z
/// - Entropy: S = (<E> - F) / T
/// - Heat capacity: C_v = (<E^2> - <E>^2) / (kT^2)
public struct ThermodynamicResult: Sendable, Codable, Hashable {
    /// Scientific provenance. Missing legacy metadata fails closed to proxy-only.
    public let scientificProvenance: ScientificProvenance?

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

    /// Temperature in Kelvin
    public let temperature: Double

    /// Natural log of the partition function ln(Z)
    public let logZ: Double

    /// Helmholtz-like F = -kT ln Z in the declared energy domain.
    /// Units are kcal/mol only for calibrated provenance.
    public let freeEnergy: Double

    /// Boltzmann-weighted mean energy <E> in the declared energy domain.
    public let meanEnergy: Double

    /// Mean squared energy <E^2>
    public let meanEnergySq: Double

    /// Heat capacity C_v = (<E^2> - <E>^2) / (kT^2)
    public let heatCapacity: Double

    /// Conformational entropy-like value S = (<E> - F) / T in the declared
    /// energy domain per kelvin; physical kcal mol^-1 K^-1 requires calibration.
    public let entropy: Double

    /// Standard deviation in the declared energy domain.
    public let stdEnergy: Double

    /// Initialize from a C FXThermodynamics struct
    init(from c: FXThermodynamics) {
        self.scientificProvenance = ScientificProvenance(from: c.scientific_provenance)
        self.temperature = c.temperature
        self.logZ = c.log_Z
        self.freeEnergy = c.free_energy
        self.meanEnergy = c.mean_energy
        self.meanEnergySq = c.mean_energy_sq
        self.heatCapacity = c.heat_capacity
        self.entropy = c.entropy
        self.stdEnergy = c.std_energy
    }

    /// Initialize with explicit values
    public init(
        temperature: Double, logZ: Double, freeEnergy: Double,
        meanEnergy: Double, meanEnergySq: Double,
        heatCapacity: Double, entropy: Double, stdEnergy: Double,
        scientificProvenance: ScientificProvenance? = nil
    ) {
        self.scientificProvenance = scientificProvenance
        self.temperature = temperature
        self.logZ = logZ
        self.freeEnergy = freeEnergy
        self.meanEnergy = meanEnergy
        self.meanEnergySq = meanEnergySq
        self.heatCapacity = heatCapacity
        self.entropy = entropy
        self.stdEnergy = stdEnergy
    }
}
