// BindingEntropyScore.swift — Entropy-health correlation data
//
// Links docking entropy metrics with HealthKit biometrics.
// Provides decomposed Shannon entropy (configurational + vibrational)
// for the IntelligenceOracle referee pipeline.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

/// Compatibility name for the decomposition produced by DockingRunner.
public typealias ShannonEntropyDecomposition = FlexAIDdS.ShannonEntropyDecomposition

// MARK: - Binding Entropy Score

/// Correlates binding population entropy with health biometrics.
///
/// This is the core data object that connects molecular docking results
/// with HealthKit measurements (HRV, sleep, etc.) and the decomposed
/// Shannon entropy for the IntelligenceOracle referee pipeline.
public struct BindingEntropyScore: Sendable, Codable, Hashable, Identifiable {
    public let id: UUID

    // MARK: - Docking Metrics

    /// Shannon configurational entropy of the binding population
    public let shannonS: Double

    /// Temperature of the simulation (K)
    public let temperature: Double

    /// Number of binding modes in the population
    public let bindingModeCount: Int

    /// Helmholtz free energy of the lowest-energy binding mode (kcal/mol)
    public let bestFreeEnergy: Double

    /// Heat capacity of the global ensemble
    public let heatCapacity: Double

    // MARK: - Shannon Decomposition (from ShannonThermoStack)

    /// Decomposed entropy with configurational/vibrational split and convergence data.
    /// Available when `DockingRunner` has access to the ShannonThermoStack output.
    public var shannonDecomposition: ShannonEntropyDecomposition?

    // MARK: - Health Metrics (optional, filled when HealthKit data available)

    /// Heart rate variability SDNN (ms), if available
    public var hrvSDNN: Double?

    /// Resting heart rate (bpm), if available
    public var restingHeartRate: Double?

    /// Sleep duration (hours), if available
    public var sleepHours: Double?

    /// Timestamp of measurement
    public let timestamp: Date

    /// Whether the entropy has "collapsed" (shannonS below threshold)
    public var isCollapsed: Bool {
        shannonS < 0.1  // Threshold for entropy collapse
    }

    /// Whether the Shannon entropy has converged (plateau reached)
    public var isConverged: Bool {
        shannonDecomposition?.isConverged ?? false
    }

    /// Vibrational entropy contribution, if decomposition available (kcal/mol/K)
    public var vibrationalEntropy: Double? {
        shannonDecomposition?.vibrational
    }

    /// Configurational entropy contribution, if decomposition available (nats)
    public var configurationalEntropy: Double? {
        shannonDecomposition?.configurational
    }

    public init(
        shannonS: Double, temperature: Double,
        bindingModeCount: Int, bestFreeEnergy: Double,
        heatCapacity: Double, timestamp: Date = Date()
    ) {
        self.id = UUID()
        self.shannonS = shannonS
        self.temperature = temperature
        self.bindingModeCount = bindingModeCount
        self.bestFreeEnergy = bestFreeEnergy
        self.heatCapacity = heatCapacity
        self.timestamp = timestamp
    }

    /// Create from a DockingResult with optional Shannon decomposition.
    public init(from result: DockingResult, decomposition: ShannonEntropyDecomposition? = nil) {
        self.id = UUID()
        self.shannonS = result.globalThermodynamics.entropy
        self.temperature = result.temperature
        self.bindingModeCount = result.bindingModes.count
        self.bestFreeEnergy = result.bindingModes.first?.freeEnergy ?? 0
        self.heatCapacity = result.globalThermodynamics.heatCapacity
        self.timestamp = result.timestamp
        self.shannonDecomposition = decomposition
    }
}
