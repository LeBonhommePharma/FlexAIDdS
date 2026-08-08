// ScientificProvenance.swift — fail-closed scientific claim domains
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDCore

/// Domain of values supplied to statistical-mechanics calculations.
public enum EnergyDomain: String, Sendable, Codable, Hashable {
    case unclassified
    case cfArbitraryUnits = "cf_arbitrary_units"
    case calibratedKcalPerMol = "calibrated_kcal_per_mol"
    case modelScale = "model_scale"
}

/// Measure represented by records in an ensemble.
public enum EnsembleMeasure: String, Sendable, Codable, Hashable {
    case unclassified
    case optimizerSamples = "optimizer_samples"
    case enumeratedMicrostates = "enumerated_microstates"
    case weightedQuadrature = "weighted_quadrature"
}

/// Completeness of the reference state used for association claims.
public enum ReferenceState: String, Sendable, Codable, Hashable {
    case none
    case boundOnly = "bound_only"
    case matchedAssociationCycle = "matched_association_cycle"
}

/// Strongest scientific claim supported by declared provenance.
public enum ClaimValidity: String, Sendable, Codable, Hashable {
    case proxyOnly = "proxy_only"
    case canonicalPhysical = "canonical_physical"
    case bindingPhysical = "binding_physical"
}

/// Metadata required before thermodynamic values may be interpreted physically.
/// Claim validity is derived and is never accepted as a serialized override.
public struct ScientificProvenance: Sendable, Codable, Hashable {
    public static let currentSchemaVersion = 2

    public let schemaVersion: Int
    public let energyDomain: EnergyDomain
    public let ensembleMeasure: EnsembleMeasure
    public let referenceState: ReferenceState
    public let energyProvenance: String
    public let measureProvenance: String
    public let referenceProvenance: String

    public init(
        schemaVersion: Int = ScientificProvenance.currentSchemaVersion,
        energyDomain: EnergyDomain = .unclassified,
        ensembleMeasure: EnsembleMeasure = .unclassified,
        referenceState: ReferenceState = .none,
        energyProvenance: String = "",
        measureProvenance: String = "",
        referenceProvenance: String = ""
    ) {
        self.schemaVersion = schemaVersion
        self.energyDomain = energyDomain
        self.ensembleMeasure = ensembleMeasure
        self.referenceState = referenceState
        self.energyProvenance = energyProvenance
        self.measureProvenance = measureProvenance
        self.referenceProvenance = referenceProvenance
    }

    /// Mirror of `statmech::make_contact_function_optimizer_provenance`
    /// (LIB/statmech.cpp). C++ is the single source of truth for this
    /// vocabulary, so the descriptive strings must match byte-for-byte: a
    /// Swift-side paraphrase would make the same ensemble serialize two
    /// different evidence records. The strings are deliberately prose rather
    /// than `sha256:` receipts, so this factory is structurally incapable of
    /// authorizing a canonical or binding physical claim.
    public static func contactFunctionOptimizer(
        referenceState: ReferenceState = .boundOnly
    ) -> ScientificProvenance {
        ScientificProvenance(
            energyDomain: .cfArbitraryUnits,
            ensembleMeasure: .optimizerSamples,
            referenceState: referenceState,
            energyProvenance:
                "FlexAID Voronoi/contact-function score; no kcal/mol calibration",
            measureProvenance:
                "optimizer-selected, deduplicated and/or clustered GA pose records",
            referenceProvenance: referenceState == .boundOnly
                ? "bound pose ensemble only"
                : "no matched association-cycle artifact"
        )
    }

    public static let proxyContactFunction = ScientificProvenance
        .contactFunctionOptimizer(referenceState: .boundOnly)

    private static let rejectedArtifactDigests: Set<String> = [
        // Empty artifacts cannot substantiate a scientific calibration receipt.
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        // Historical entropy.help example filler, never derived from an artifact.
        "3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a",
    ]

    /// Physical claims require an exact artifact identity, not descriptive
    /// prose. The stable receipt form is `sha256:` followed by 64 hexadecimal
    /// digits. Known placeholders and low-diversity digests are rejected.
    private static func hasArtifactIdentity(_ value: String) -> Bool {
        let bytes = Array(value.utf8)
        let prefix = Array("sha256:".utf8)
        guard bytes.count == prefix.count + 64,
              Array(bytes.prefix(prefix.count)) == prefix
        else { return false }

        let digestBytes = bytes.dropFirst(prefix.count)
        guard digestBytes.allSatisfy({ byte in
            (byte >= 0x30 && byte <= 0x39) ||
            (byte >= 0x41 && byte <= 0x46) ||
            (byte >= 0x61 && byte <= 0x66)
        }) else { return false }

        let digest = String(decoding: digestBytes, as: UTF8.self).lowercased()
        return !Self.rejectedArtifactDigests.contains(digest)
            && Set(digest.utf8).count > 2
    }

    public var claimValidity: ClaimValidity {
        guard schemaVersion == Self.currentSchemaVersion,
              energyDomain == .calibratedKcalPerMol,
              ensembleMeasure == .enumeratedMicrostates || ensembleMeasure == .weightedQuadrature,
              Self.hasArtifactIdentity(energyProvenance),
              Self.hasArtifactIdentity(measureProvenance)
        else { return .proxyOnly }

        guard referenceState == .matchedAssociationCycle,
              Self.hasArtifactIdentity(referenceProvenance)
        else { return .canonicalPhysical }

        return .bindingPhysical
    }

    public var allowsCanonicalPhysicalClaim: Bool {
        claimValidity == .canonicalPhysical || claimValidity == .bindingPhysical
    }

    public var allowsBindingPhysicalClaim: Bool {
        claimValidity == .bindingPhysical
    }

    public var isProxyOnly: Bool {
        claimValidity == .proxyOnly
    }

    // MARK: - Stable schema-v2 wire format

    private enum CodingKeys: String, CodingKey {
        case schemaVersion = "schema_version"
        case energyDomain = "energy_domain"
        case ensembleMeasure = "ensemble_measure"
        case referenceState = "reference_state"
        case energyProvenance = "energy_provenance"
        case measureProvenance = "measure_provenance"
        case referenceProvenance = "reference_provenance"
        case claimValidity = "claim_validity"
    }

    private enum LegacyCodingKeys: String, CodingKey {
        case schemaVersion
        case energyDomain
        case ensembleMeasure
        case referenceState
        case energyProvenance
        case measureProvenance
        case referenceProvenance
        case claimValidity
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(schemaVersion, forKey: .schemaVersion)
        try container.encode(energyDomain.rawValue, forKey: .energyDomain)
        try container.encode(ensembleMeasure.rawValue, forKey: .ensembleMeasure)
        try container.encode(referenceState.rawValue, forKey: .referenceState)
        try container.encode(energyProvenance, forKey: .energyProvenance)
        try container.encode(measureProvenance, forKey: .measureProvenance)
        try container.encode(referenceProvenance, forKey: .referenceProvenance)
        // Informational only. Decoding always derives validity from evidence.
        try container.encode(claimValidity.rawValue, forKey: .claimValidity)
    }

    public init(from decoder: Decoder) throws {
        let canonical = try decoder.container(keyedBy: CodingKeys.self)
        let legacy = try decoder.container(keyedBy: LegacyCodingKeys.self)

        func decode<T: Decodable>(
            _ type: T.Type,
            canonicalKey: CodingKeys,
            legacyKey: LegacyCodingKeys
        ) -> T? {
            // A present canonical key is authoritative. If it has a hostile
            // type, discard it and fail closed instead of consulting legacy.
            if canonical.contains(canonicalKey) {
                return try? canonical.decode(type, forKey: canonicalKey)
            }
            return try? legacy.decodeIfPresent(type, forKey: legacyKey)
        }

        let schemaVersion = decode(
            Int.self,
            canonicalKey: .schemaVersion,
            legacyKey: .schemaVersion
        ) ?? 0
        let energyDomain = decode(
            String.self,
            canonicalKey: .energyDomain,
            legacyKey: .energyDomain
        ).flatMap(EnergyDomain.init(rawValue:)) ?? .unclassified
        let ensembleMeasure = decode(
            String.self,
            canonicalKey: .ensembleMeasure,
            legacyKey: .ensembleMeasure
        ).flatMap(EnsembleMeasure.init(rawValue:)) ?? .unclassified
        let referenceState = decode(
            String.self,
            canonicalKey: .referenceState,
            legacyKey: .referenceState
        ).flatMap(ReferenceState.init(rawValue:)) ?? .none

        self.init(
            schemaVersion: schemaVersion,
            energyDomain: energyDomain,
            ensembleMeasure: ensembleMeasure,
            referenceState: referenceState,
            energyProvenance: decode(
                String.self,
                canonicalKey: .energyProvenance,
                legacyKey: .energyProvenance
            ) ?? "",
            measureProvenance: decode(
                String.self,
                canonicalKey: .measureProvenance,
                legacyKey: .measureProvenance
            ) ?? "",
            referenceProvenance: decode(
                String.self,
                canonicalKey: .referenceProvenance,
                legacyKey: .referenceProvenance
            ) ?? ""
        )
        // Deliberately do not decode either spelling of claimValidity.
    }

    // MARK: - C bridge

    /// Adopt evidence produced by the C++ engine. Public so feature modules
    /// can gate presentation on the record the bridge actually carried,
    /// instead of constructing a stand-in provenance locally.
    public init(from c: FXScientificProvenance) {
        self.init(
            schemaVersion: Int(c.schema_version),
            energyDomain: Self.energyDomain(from: c.energy_domain),
            ensembleMeasure: Self.ensembleMeasure(from: c.ensemble_measure),
            referenceState: Self.referenceState(from: c.reference_state),
            energyProvenance: Self.string(fromCStringStorage: c.energy_provenance),
            measureProvenance: Self.string(fromCStringStorage: c.measure_provenance),
            referenceProvenance: Self.string(fromCStringStorage: c.reference_provenance)
        )
    }

    private static func energyDomain(from rawValue: Int32) -> EnergyDomain {
        switch rawValue {
        case Int32(FX_ENERGY_DOMAIN_CF_ARBITRARY_UNITS): return .cfArbitraryUnits
        case Int32(FX_ENERGY_DOMAIN_CALIBRATED_KCAL_PER_MOL): return .calibratedKcalPerMol
        case Int32(FX_ENERGY_DOMAIN_MODEL_SCALE): return .modelScale
        default: return .unclassified
        }
    }

    private static func ensembleMeasure(from rawValue: Int32) -> EnsembleMeasure {
        switch rawValue {
        case Int32(FX_ENSEMBLE_MEASURE_OPTIMIZER_SAMPLES): return .optimizerSamples
        case Int32(FX_ENSEMBLE_MEASURE_ENUMERATED_MICROSTATES): return .enumeratedMicrostates
        case Int32(FX_ENSEMBLE_MEASURE_WEIGHTED_QUADRATURE): return .weightedQuadrature
        default: return .unclassified
        }
    }

    private static func referenceState(from rawValue: Int32) -> ReferenceState {
        switch rawValue {
        case Int32(FX_REFERENCE_STATE_BOUND_ONLY): return .boundOnly
        case Int32(FX_REFERENCE_STATE_MATCHED_ASSOCIATION_CYCLE): return .matchedAssociationCycle
        default: return .none
        }
    }

    private static func string<T>(fromCStringStorage storage: T) -> String {
        var storage = storage
        return withUnsafeBytes(of: &storage) { rawBuffer in
            let bytes = rawBuffer.prefix { $0 != 0 }
            return String(decoding: bytes, as: UTF8.self)
        }
    }
}
