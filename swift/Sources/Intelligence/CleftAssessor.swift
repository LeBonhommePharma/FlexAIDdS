// CleftAssessor.swift — Binding site druggability assessment
//
// Assesses pocket quality and druggability from cavity geometry using
// Apple FoundationModels. Pre-computes all geometric features on CPU;
// the LLM only provides natural-language interpretation.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import Foundation
import FlexAIDdS

// MARK: - Pre-Computed Cleft Features

/// Summary of a detected binding cleft for LLM interpretation.
public struct CleftFeatures: Sendable, Codable {
    /// Evidence governing whether these geometric features may be turned into
    /// druggability / ligand-design claims. Geometry alone never authorizes a
    /// potency statement, so this fails closed to proxy-only.
    public let scientificProvenance: ScientificProvenance?

    /// Strongest claim supported by the declared provenance.
    public var claimValidity: ClaimValidity {
        scientificProvenance?.claimValidity ?? .proxyOnly
    }

    /// Druggability and ligand-property prescriptions are binding claims.
    public var allowsDrugDesignClaims: Bool {
        claimValidity == .bindingPhysical
    }

    /// Volume in cubic Angstroms
    public let volume: Double
    /// Effective pocket depth (Angstroms)
    public let depth: Double
    /// Number of probe spheres defining the pocket
    public let sphereCount: Int
    /// Largest probe sphere radius (Angstroms)
    public let maxSphereRadius: Double
    /// Fraction of pocket surface that is hydrophobic (0.0-1.0)
    public let hydrophobicFraction: Double
    /// Number of anchor residues lining the pocket
    public let anchorResidueCount: Int
    /// Pocket shape: 0 = spherical, 1 = elongated
    public let elongation: Double
    /// Fraction of pocket exposed to solvent (0.0-1.0)
    public let solventExposure: Double

    public init(volume: Double, depth: Double, sphereCount: Int,
                maxSphereRadius: Double, hydrophobicFraction: Double,
                anchorResidueCount: Int, elongation: Double, solventExposure: Double,
                scientificProvenance: ScientificProvenance? = nil) {
        self.scientificProvenance = scientificProvenance
        self.volume = volume
        self.depth = depth
        self.sphereCount = sphereCount
        self.maxSphereRadius = maxSphereRadius
        self.hydrophobicFraction = hydrophobicFraction
        self.anchorResidueCount = anchorResidueCount
        self.elongation = elongation
        self.solventExposure = solventExposure
    }
}

// MARK: - Output Types

#if canImport(FoundationModels)
import FoundationModels

/// Druggability tier.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public enum DruggabilityTier: String, Sendable, Codable {
    case high
    case moderate
    case low
    case undruggable
    /// No druggability verdict may be issued for this record. Distinct from
    /// `undruggable`, which is itself a drug-design claim.
    case unavailable
}

/// Assessment of a binding site's druggability.
@available(macOS 26.0, iOS 26.0, *)
@Generable
public struct CleftAssessment: Sendable, Codable {
    /// Druggability tier
    public var druggability: DruggabilityTier
    /// Plain-English pocket description
    public var summary: String
    /// Suggested ligand properties (e.g., "lipophilic, <400 Da")
    public var suggestedLigandProperties: String
    /// Any warnings about the binding site
    public var warnings: [String]
}

// MARK: - FoundationModels Actor

@available(macOS 26.0, iOS 26.0, *)
public actor CleftDruggabilityAssessor {
    private let session: LanguageModelSession

    private static let instructions = """
        You are a structural biology advisor assessing binding site druggability. \
        All geometric values are pre-computed — DO NOT perform calculations. \
        Assess whether the pocket can accommodate a drug-like small molecule. \
        Consider volume (200-1000 ų is drug-like), hydrophobicity, depth, \
        anchor residues, and solvent exposure. Emit druggability tiers and \
        ligand-property prescriptions only for binding_physical inputs. \
        Be concise and actionable.
        """

    public init() {
        self.session = LanguageModelSession(instructions: Self.instructions)
    }

    /// Assess druggability of a binding cleft.
    ///
    /// Druggability tiers and ligand-property prescriptions are drug-design
    /// claims. Without a calibrated matched association cycle the pocket can
    /// only be described geometrically, so the model is not consulted.
    public func assess(cleft: CleftFeatures) async throws -> CleftAssessment {
        guard cleft.allowsDrugDesignClaims else {
            let shape = cleft.elongation < 0.3 ? "spherical"
                : cleft.elongation < 0.7 ? "oval" : "elongated"
            return CleftAssessment(
                druggability: .unavailable,
                summary: "Geometry-only cavity descriptor: \(String(format: "%.0f", cleft.volume)) ų \(shape) pocket, \(String(format: "%.1f", cleft.depth)) Å deep, \(String(format: "%.0f", cleft.hydrophobicFraction * 100))% hydrophobic surface. No druggability tier is assigned (claim_validity=\(cleft.claimValidity.rawValue)).",
                suggestedLigandProperties: "Ligand-property prescriptions are unavailable without binding_physical provenance.",
                warnings: ["Druggability, affinity and potency claims are unavailable; this record is a cavity-geometry diagnostic only."]
            )
        }
        let prompt = buildPrompt(cleft: cleft)
        return try await session.respond(to: prompt, generating: CleftAssessment.self).content
    }

    private func buildPrompt(cleft: CleftFeatures) -> String {
        var p = "Assess this binding pocket. Produce a CleftAssessment.\n"
        p += "Claim validity: \(cleft.claimValidity.rawValue)\n"
        p += "Volume: \(String(format: "%.0f", cleft.volume)) ų\n"
        p += "Depth: \(String(format: "%.1f", cleft.depth)) Å\n"
        p += "Shape: \(cleft.elongation < 0.3 ? "spherical" : cleft.elongation < 0.7 ? "oval" : "elongated") (index \(String(format: "%.2f", cleft.elongation)))\n"
        p += "Hydrophobic: \(String(format: "%.0f", cleft.hydrophobicFraction * 100))%\n"
        p += "Anchor residues: \(cleft.anchorResidueCount)\n"
        p += "Solvent exposure: \(String(format: "%.0f", cleft.solventExposure * 100))%\n"
        p += "Max sphere: \(String(format: "%.1f", cleft.maxSphereRadius)) Å, \(cleft.sphereCount) probes\n"

        if cleft.volume < 200 { p += "FLAG: Very small pocket (<200 ų)\n" }
        if cleft.volume > 1200 { p += "FLAG: Very large pocket (>1200 ų) — may lack specificity\n" }
        if cleft.solventExposure > 0.6 { p += "FLAG: Highly solvent-exposed — binding may be weak\n" }
        if cleft.hydrophobicFraction > 0.8 { p += "FLAG: Heavily hydrophobic — select lipophilic compounds\n" }

        return p
    }
}
#endif

// MARK: - Cross-Platform Output

/// Platform-independent druggability tier.
public enum CrossPlatformDruggabilityTier: String, Sendable, Codable {
    case high, moderate, low, undruggable
    /// No druggability verdict may be issued for this record.
    case unavailable
}

/// Platform-independent cleft assessment.
public struct CrossPlatformCleftAssessment: Sendable, Codable {
    public let druggability: CrossPlatformDruggabilityTier
    public let summary: String
    public let suggestedLigandProperties: String
    public let warnings: [String]

    public init(druggability: CrossPlatformDruggabilityTier, summary: String,
                suggestedLigandProperties: String, warnings: [String]) {
        self.druggability = druggability
        self.summary = summary
        self.suggestedLigandProperties = suggestedLigandProperties
        self.warnings = warnings
    }
}

// MARK: - Rule-Based Fallback

/// Deterministic druggability assessor for non-Apple platforms.
public struct RuleBasedCleftAssessor: Sendable {

    public init() {}

    /// Assess druggability using geometric thresholds.
    ///
    /// Cavity geometry can be reported at any provenance level, but a
    /// druggability tier and a ligand-property prescription are drug-design
    /// claims and require `binding_physical` evidence.
    public func assess(cleft: CleftFeatures) -> CrossPlatformCleftAssessment {
        guard cleft.allowsDrugDesignClaims else {
            let shape = cleft.elongation < 0.3 ? "spherical"
                : cleft.elongation < 0.7 ? "oval" : "elongated"
            return CrossPlatformCleftAssessment(
                druggability: .unavailable,
                summary: "Geometry-only cavity descriptor: \(String(format: "%.0f", cleft.volume)) ų \(shape) pocket, \(String(format: "%.1f", cleft.depth)) Å deep, \(String(format: "%.0f", cleft.hydrophobicFraction * 100))% hydrophobic surface, \(cleft.anchorResidueCount) lining residues. No druggability tier is assigned (claim_validity=\(cleft.claimValidity.rawValue)).",
                suggestedLigandProperties: "Ligand-property prescriptions are unavailable without binding_physical provenance.",
                warnings: ["Druggability, affinity and potency claims are unavailable; this record is a cavity-geometry diagnostic only."]
            )
        }

        var warnings: [String] = []
        var score = 0.0

        // Volume scoring (200-1000 ų is ideal for drug-like molecules)
        if cleft.volume >= 200 && cleft.volume <= 1000 {
            score += 0.3
        } else if cleft.volume < 200 {
            warnings.append("Pocket too small (\(String(format: "%.0f", cleft.volume)) ų) for drug-like molecules.")
        } else if cleft.volume > 1200 {
            warnings.append("Pocket very large (\(String(format: "%.0f", cleft.volume)) ų) — may lack specificity.")
        } else {
            score += 0.15
        }

        // Hydrophobicity scoring (40-80% is drug-like)
        if cleft.hydrophobicFraction >= 0.4 && cleft.hydrophobicFraction <= 0.8 {
            score += 0.25
        } else if cleft.hydrophobicFraction > 0.8 {
            score += 0.15
            warnings.append("Heavily hydrophobic pocket — limited polar interaction sites.")
        } else {
            score += 0.1
        }

        // Depth scoring (deeper pockets bind better)
        if cleft.depth > 6.0 { score += 0.2 }
        else if cleft.depth > 3.0 { score += 0.1 }
        else { warnings.append("Shallow pocket (\(String(format: "%.1f", cleft.depth)) Å) — binding may be transient.") }

        // Solvent exposure (lower is better for binding)
        if cleft.solventExposure < 0.3 { score += 0.15 }
        else if cleft.solventExposure > 0.6 {
            warnings.append("High solvent exposure (\(String(format: "%.0f", cleft.solventExposure * 100))%) — desolvation penalty.")
        } else { score += 0.05 }

        // Anchor residues
        if cleft.anchorResidueCount >= 4 { score += 0.1 }

        let tier: CrossPlatformDruggabilityTier
        if score >= 0.7 { tier = .high }
        else if score >= 0.45 { tier = .moderate }
        else if score >= 0.2 { tier = .low }
        else { tier = .undruggable }

        let volumeDesc = cleft.volume < 400 ? "small" : cleft.volume < 800 ? "medium" : "large"
        let shapeDesc = cleft.elongation < 0.3 ? "spherical" : cleft.elongation < 0.7 ? "oval" : "elongated"
        let summary = "\(volumeDesc.capitalized) \(shapeDesc) pocket (\(String(format: "%.0f", cleft.volume)) ų, \(String(format: "%.0f", cleft.hydrophobicFraction * 100))% hydrophobic, \(cleft.anchorResidueCount) anchor residues). Druggability: \(tier.rawValue)."

        let ligandProps: String
        if cleft.hydrophobicFraction > 0.6 {
            ligandProps = "Lipophilic compounds, MW 300-500 Da, LogP > 2"
        } else if cleft.hydrophobicFraction < 0.4 {
            ligandProps = "Polar compounds with H-bond donors/acceptors, MW 200-400 Da"
        } else {
            ligandProps = "Balanced compounds, MW 300-500 Da, 2-3 H-bond acceptors"
        }

        return CrossPlatformCleftAssessment(
            druggability: tier, summary: summary,
            suggestedLigandProperties: ligandProps, warnings: warnings
        )
    }
}
