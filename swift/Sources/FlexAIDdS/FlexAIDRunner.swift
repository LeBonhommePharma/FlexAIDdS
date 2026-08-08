// FlexAIDRunner.swift — Swift actor wrapping the StatMechEngine
//
// Thread-safe access to the C++ StatMechEngine via actor isolation.
// The underlying engine is NOT thread-safe, so all access is serialized.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import FlexAIDCore

/// Actor providing thread-safe access to the FlexAID statistical mechanics engine.
///
/// Usage:
/// ```swift
/// let runner = FlexAIDRunner(temperature: 300.0)
/// runner.addSample(energy: -10.5)
/// runner.addSample(energy: -8.3, multiplicity: 2)
/// let result = runner.compute()
/// print("Free energy: \(result.freeEnergy) kcal/mol")
/// print("Entropy: \(result.entropy) kcal/mol/K")
/// ```
public actor FlexAIDRunner {
    private let engineRef: FXStatMechEngineRef

    /// Create a new statistical mechanics engine at the given temperature.
    /// Numeric inputs are unclassified unless the producer explicitly supplies
    /// provenance; no energy domain or ensemble measure is inferred.
    public init(
        temperature: Double = 300.0,
        scientificProvenance: ScientificProvenance? = nil
    ) {
        guard let engine = fx_statmech_create(temperature) else {
            preconditionFailure("Failed to create StatMechEngine")
        }
        self.engineRef = engine
        if let provenance = scientificProvenance {
            let schemaVersion = Int32(exactly: provenance.schemaVersion) ?? 0
            let energyProvenance = Self.bridgeEvidence(provenance.energyProvenance)
            let measureProvenance = Self.bridgeEvidence(provenance.measureProvenance)
            let referenceProvenance = Self.bridgeEvidence(provenance.referenceProvenance)
            energyProvenance.withCString { energyEvidence in
                measureProvenance.withCString { measureEvidence in
                    referenceProvenance.withCString { referenceEvidence in
                        fx_statmech_set_scientific_provenance(
                            engine,
                            schemaVersion,
                            Self.cEnergyDomain(provenance.energyDomain),
                            Self.cEnsembleMeasure(provenance.ensembleMeasure),
                            Self.cReferenceState(provenance.referenceState),
                            energyEvidence,
                            measureEvidence,
                            referenceEvidence
                        )
                    }
                }
            }
        }
    }

    /// C strings cannot preserve bytes after an embedded NUL. Empty the field
    /// instead of allowing a hostile prefix to become an authoritative receipt.
    private nonisolated static func bridgeEvidence(_ value: String) -> String {
        value.utf8.contains(0) ? "" : value
    }

    private nonisolated static func cEnergyDomain(_ value: EnergyDomain) -> Int32 {
        switch value {
        case .unclassified: return Int32(FX_ENERGY_DOMAIN_UNCLASSIFIED)
        case .cfArbitraryUnits: return Int32(FX_ENERGY_DOMAIN_CF_ARBITRARY_UNITS)
        case .calibratedKcalPerMol: return Int32(FX_ENERGY_DOMAIN_CALIBRATED_KCAL_PER_MOL)
        case .modelScale: return Int32(FX_ENERGY_DOMAIN_MODEL_SCALE)
        }
    }

    private nonisolated static func cEnsembleMeasure(_ value: EnsembleMeasure) -> Int32 {
        switch value {
        case .unclassified: return Int32(FX_ENSEMBLE_MEASURE_UNCLASSIFIED)
        case .optimizerSamples: return Int32(FX_ENSEMBLE_MEASURE_OPTIMIZER_SAMPLES)
        case .enumeratedMicrostates: return Int32(FX_ENSEMBLE_MEASURE_ENUMERATED_MICROSTATES)
        case .weightedQuadrature: return Int32(FX_ENSEMBLE_MEASURE_WEIGHTED_QUADRATURE)
        }
    }

    private nonisolated static func cReferenceState(_ value: ReferenceState) -> Int32 {
        switch value {
        case .none: return Int32(FX_REFERENCE_STATE_NONE)
        case .boundOnly: return Int32(FX_REFERENCE_STATE_BOUND_ONLY)
        case .matchedAssociationCycle: return Int32(FX_REFERENCE_STATE_MATCHED_ASSOCIATION_CYCLE)
        }
    }

    deinit {
        fx_statmech_destroy(engineRef)
    }

    // MARK: - Sample Management

    /// Add a sampled conformation to the ensemble.
    /// - Parameters:
    ///   - energy: Value in the explicitly declared energy domain (negative = favorable)
    ///   - multiplicity: Degeneracy / sampling count (default: 1)
    public func addSample(energy: Double, multiplicity: Int = 1) {
        fx_statmech_add_sample(engineRef, energy, Int32(multiplicity))
    }

    /// Remove all samples from the ensemble.
    public func clear() {
        fx_statmech_clear(engineRef)
    }

    /// Number of samples in the ensemble.
    public var sampleCount: Int {
        Int(fx_statmech_size(engineRef))
    }

    // MARK: - Thermodynamic Computation

    /// Compute full thermodynamic properties of the current ensemble.
    /// - Returns: Complete thermodynamic analysis (F, S, <E>, C_v, etc.)
    public func compute() -> ThermodynamicResult {
        ThermodynamicResult(from: fx_statmech_compute(engineRef))
    }

    /// Get Boltzmann weights for all samples (same order as insertion).
    /// - Returns: Array of normalized Boltzmann probabilities
    public func boltzmannWeights() -> [Double] {
        var count: Int32 = 0
        guard let ptr = fx_statmech_boltzmann_weights(engineRef, &count) else { return [] }
        defer { fx_free_doubles(ptr) }
        return Array(UnsafeBufferPointer(start: ptr, count: Int(count)))
    }

    /// Compute an ensemble free-energy difference vs a reference ensemble.
    /// This is not a binding ΔG claim unless both results carry a validated,
    /// matched association-cycle provenance.
    /// - Parameter reference: Another FlexAIDRunner with a reference ensemble
    /// - Returns: Difference in the declared energy domain
    public func deltaG(relativeTo reference: FlexAIDRunner) async -> Double {
        let refEngine = await reference.engineRef
        return fx_statmech_delta_G(engineRef, refEngine)
    }

    // MARK: - Accessors

    /// Temperature in Kelvin
    public var temperature: Double {
        fx_statmech_temperature(engineRef)
    }

    /// Inverse temperature beta = 1/(kT) in (kcal/mol)^-1
    public var beta: Double {
        fx_statmech_beta(engineRef)
    }

    // MARK: - Static / Pure Functions (nonisolated)

    /// Compute Helmholtz free energy from a raw energy array.
    /// - Parameters:
    ///   - energies: Array of energy values (kcal/mol)
    ///   - temperature: Temperature in Kelvin
    /// - Returns: Helmholtz free energy F (kcal/mol)
    public nonisolated static func helmholtz(energies: [Double], temperature: Double) -> Double {
        energies.withUnsafeBufferPointer { buf in
            fx_statmech_helmholtz(buf.baseAddress, Int32(buf.count), temperature)
        }
    }

    /// Thermodynamic integration via trapezoidal rule.
    /// - Parameter points: Array of (lambda, <dV/dlambda>) points
    /// - Returns: Delta-G in kcal/mol
    public nonisolated static func thermodynamicIntegration(points: [TIPoint]) -> Double {
        let cPoints = points.map { FXTIPoint(lambda: $0.lambda, dV_dlambda: $0.dVdLambda) }
        return cPoints.withUnsafeBufferPointer { buf in
            fx_statmech_thermodynamic_integration(buf.baseAddress, Int32(buf.count))
        }
    }

    /// Boltzmann-reweighted PMF for free energy profiles along a reaction coordinate.
    /// - Parameters:
    ///   - energies: Energy values for each sample
    ///   - coordinates: Coordinate values for each sample
    ///   - temperature: Temperature in Kelvin
    ///   - nBins: Number of histogram bins
    ///   - maxIter: Maximum iterations (default: 1000)
    ///   - tolerance: Convergence tolerance (default: 1e-6)
    /// - Returns: Array of bins with free energy profile
    public nonisolated static func boltzmannPMF(
        energies: [Double], coordinates: [Double],
        temperature: Double, nBins: Int,
        maxIter: Int = 1000, tolerance: Double = 1e-6
    ) -> [WHAMBinResult] {
        precondition(energies.count == coordinates.count, "Energies and coordinates must have equal length")

        var outCount: Int32 = 0
        let ptr = energies.withUnsafeBufferPointer { eBuf in
            coordinates.withUnsafeBufferPointer { cBuf in
                fx_statmech_boltzmann_pmf(eBuf.baseAddress, cBuf.baseAddress,
                                 Int32(eBuf.count), temperature, Int32(nBins),
                                 Int32(maxIter), tolerance, &outCount)
            }
        }

        guard let ptr = ptr else { return [] }
        defer { fx_free_wham_bins(ptr) }

        return (0..<Int(outCount)).map { i in
            WHAMBinResult(from: ptr[i])
        }
    }
}

// MARK: - ENCoM Vibrational Entropy

/// Static methods for ENCoM vibrational entropy calculations.
/// No instance state needed — all methods are pure functions.
public enum ENCoMRunner {
    /// Compute vibrational entropy from normal mode eigenvalues.
    /// - Parameters:
    ///   - eigenvalues: Array of normal mode eigenvalues
    ///   - temperature: Temperature in Kelvin (default: 300.0)
    ///   - cutoff: Skip modes with eigenvalue below this threshold (default: 1e-6)
    /// - Returns: Vibrational entropy result
    public static func computeVibrationalEntropy(
        eigenvalues: [Double],
        temperature: Double = 300.0,
        cutoff: Double = 1e-6
    ) -> VibrationalEntropyResult {
        eigenvalues.withUnsafeBufferPointer { buf in
            VibrationalEntropyResult(from:
                fx_encom_compute_vibrational_entropy(
                    buf.baseAddress, Int32(buf.count),
                    temperature, cutoff))
        }
    }

    /// Combine configurational entropy (from StatMechEngine) with vibrational entropy.
    /// - Parameters:
    ///   - configurational: S_conf in kcal mol^-1 K^-1
    ///   - vibrational: S_vib in kcal mol^-1 K^-1
    /// - Returns: Total entropy S_total = S_conf + S_vib
    public static func totalEntropy(configurational: Double, vibrational: Double) -> Double {
        fx_encom_total_entropy(configurational, vibrational)
    }

    /// Free energy with vibrational correction.
    /// - Parameters:
    ///   - electronic: F_electronic from BindingMode (kcal/mol)
    ///   - vibrationalEntropy: S_vib (kcal mol^-1 K^-1)
    ///   - temperature: Temperature in Kelvin
    /// - Returns: F_total = F_elec - T * S_vib
    public static func freeEnergyWithVibrations(
        electronic: Double,
        vibrationalEntropy: Double,
        temperature: Double
    ) -> Double {
        fx_encom_free_energy_with_vibrations(electronic, vibrationalEntropy, temperature)
    }
}

// MARK: - Boltzmann Lookup Table

/// Pre-tabulated exp(-beta*E) for O(1) inner-loop evaluation.
public final class BoltzmannLookup: @unchecked Sendable {
    private let lutRef: FXBoltzmannLUTRef

    /// Create a Boltzmann lookup table.
    /// - Parameters:
    ///   - beta: Inverse temperature 1/(kT)
    ///   - eMin: Minimum energy in the table
    ///   - eMax: Maximum energy in the table
    ///   - nBins: Number of bins (default: 10000)
    public init(beta: Double, eMin: Double, eMax: Double, nBins: Int = 10000) {
        self.lutRef = fx_lut_create(beta, eMin, eMax, Int32(nBins))
    }

    deinit {
        fx_lut_destroy(lutRef)
    }

    /// Look up exp(-beta * energy) in O(1) time.
    /// - Parameter energy: Energy value in kcal/mol
    /// - Returns: Boltzmann factor exp(-beta * energy)
    public func lookup(energy: Double) -> Double {
        fx_lut_lookup(lutRef, energy)
    }
}
