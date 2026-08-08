// Models.swift — Swift data models for FlexAIDdS
//
// Value types for vibrational entropy, poses, binding modes, and physical constants.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import FlexAIDCore
import Foundation

// MARK: - Physical Constants

/// Boltzmann constant in kcal mol^-1 K^-1
public let kBkcal: Double = FX_KB_KCAL

/// Boltzmann constant in J K^-1
public let kBSI: Double = FX_KB_SI

// MARK: - Vibrational Entropy

/// Result of ENCoM vibrational entropy calculation via Schlitter formula.
public struct VibrationalEntropyResult: Sendable, Codable, Hashable {
    /// Vibrational entropy (kcal mol^-1 K^-1)
    public let entropy: Double

    /// Vibrational entropy in SI units (J mol^-1 K^-1)
    public let entropySI: Double

    /// Effective frequency omega_eff (rad/s)
    public let effectiveFrequency: Double

    /// Number of non-zero normal modes (3N - 6)
    public let modeCount: Int

    /// Temperature (K)
    public let temperature: Double

    init(from c: FXVibrationalEntropy) {
        self.entropy = c.S_vib_kcal_mol_K
        self.entropySI = c.S_vib_J_mol_K
        self.effectiveFrequency = c.omega_eff
        self.modeCount = Int(c.n_modes)
        self.temperature = c.temperature
    }
}

// MARK: - WHAM Free Energy Profile

/// A bin from weighted histogram analysis (WHAM).
public struct WHAMBinResult: Sendable, Codable, Hashable {
    /// Center of the coordinate bin
    public let coordCenter: Double

    /// Count of samples in this bin
    public let count: Double

    /// Free energy at this coordinate (kcal/mol)
    public let freeEnergy: Double

    init(from c: FXWHAMBin) {
        self.coordCenter = c.coord_center
        self.count = c.count
        self.freeEnergy = c.free_energy
    }
}

// MARK: - Thermodynamic Integration

/// A point for thermodynamic integration (TI).
public struct TIPoint: Sendable, Codable, Hashable {
    /// Coupling parameter lambda in [0, 1]
    public let lambda: Double

    /// Ensemble average <dV/dlambda> at this lambda
    public let dVdLambda: Double

    public init(lambda: Double, dVdLambda: Double) {
        self.lambda = lambda
        self.dVdLambda = dVdLambda
    }
}

// MARK: - Pose

/// Lightweight view of a molecular pose from the GA ensemble.
public struct PoseResult: Sendable, Codable, Hashable {
    /// Index in the chromosome array
    public let chromIndex: Int

    /// OPTICS clustering order
    public let order: Int

    /// Reachability distance for density-based clustering
    public let reachDist: Float

    /// Complementarity-function score in arbitrary CF units (negative = favorable)
    public let cf: Double

    public init(
        cf: Double,
        reachDist: Float,
        chromIndex: Int = 0,
        order: Int = 0
    ) {
        self.chromIndex = chromIndex
        self.order = order
        self.reachDist = reachDist
        self.cf = cf
    }

    init(from c: FXPoseInfo) {
        self.chromIndex = Int(c.chrom_index)
        self.order = Int(c.order)
        self.reachDist = c.reach_dist
        self.cf = c.cf
    }
}

// MARK: - Binding Mode

/// Summary of a binding mode (cluster of poses) with legacy ensemble diagnostics.
public struct BindingModeResult: Sendable, Codable, Hashable {
    /// Number of poses in this binding mode
    public let size: Int

    /// Legacy ensemble score; not physical Helmholtz F without provenance.
    public let freeEnergy: Double

    /// Legacy entropy-like value in the input score domain per kelvin.
    public let entropy: Double

    /// Boltzmann-weighted mean energy <E> (kcal/mol)
    public let enthalpy: Double

    /// Heat capacity C_v
    public let heatCapacity: Double

    /// Full thermodynamic result (when available)
    public var thermodynamics: ThermodynamicResult?

    /// Poses belonging to this binding mode.
    public internal(set) var poses: [PoseResult]

    public init(
        freeEnergy: Double,
        entropy: Double,
        enthalpy: Double,
        heatCapacity: Double,
        size: Int,
        poses: [PoseResult] = [],
        thermodynamics: ThermodynamicResult? = nil
    ) {
        self.size = size
        self.freeEnergy = freeEnergy
        self.entropy = entropy
        self.enthalpy = enthalpy
        self.heatCapacity = heatCapacity
        self.thermodynamics = thermodynamics
        self.poses = poses
    }

    init(from c: FXBindingModeInfo) {
        self.size = Int(c.size)
        self.freeEnergy = c.free_energy
        self.entropy = c.entropy
        self.enthalpy = c.enthalpy
        self.heatCapacity = c.heat_capacity
        self.thermodynamics = nil
        self.poses = []
    }
}

// MARK: - Docking Result

/// Complete docking result with binding-population ensemble diagnostics.
public struct DockingResult: Sendable, Codable {
    /// All binding modes found, sorted by free energy
    public let bindingModes: [BindingModeResult]

    /// Global ensemble thermodynamics across all binding modes
    public let globalThermodynamics: ThermodynamicResult

    /// Temperature used for the simulation (K)
    public let temperature: Double

    /// Number of chromosomes in the GA population
    public let populationSize: Int

    /// Timestamp of the docking run
    public let timestamp: Date

    public init(
        bindingModes: [BindingModeResult],
        globalThermodynamics: ThermodynamicResult,
        temperature: Double,
        populationSize: Int,
        timestamp: Date = Date()
    ) {
        self.bindingModes = bindingModes
        self.globalThermodynamics = globalThermodynamics
        self.temperature = temperature
        self.populationSize = populationSize
        self.timestamp = timestamp
    }
}
