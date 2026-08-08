import XCTest
@testable import FlexAIDdS

final class ScientificProvenanceTests: XCTestCase {
    private let energyReceipt = "sha256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c"
    private let measureReceipt = "sha256:7b27545430c950e8f5b4ba83ae3e2ad5e9fe32b83b625d8efea0e668af2782f4"
    private let referenceReceipt = "sha256:01d26f66e709d388d7b971de6204680694e7cca10b1570506a5738c6909a6442"

    func testMissingProvenanceFailsClosed() {
        let result = ThermodynamicResult(
            temperature: 300, logZ: 1, freeEnergy: -1,
            meanEnergy: -0.5, meanEnergySq: 0.25,
            heatCapacity: 0, entropy: 0, stdEnergy: 0
        )
        XCTAssertEqual(result.claimValidity, .proxyOnly)
        XCTAssertFalse(result.allowsCanonicalClaims)
        XCTAssertFalse(result.allowsBindingClaims)
    }

    func testCalibratedEnumerationAllowsCanonicalOnly() {
        let provenance = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .boundOnly,
            energyProvenance: energyReceipt,
            measureProvenance: measureReceipt
        )
        XCTAssertEqual(provenance.claimValidity, .canonicalPhysical)
    }

    func testMatchedCycleAllowsBindingClaim() {
        let provenance = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .weightedQuadrature,
            referenceState: .matchedAssociationCycle,
            energyProvenance: energyReceipt,
            measureProvenance: measureReceipt,
            referenceProvenance: referenceReceipt
        )
        XCTAssertEqual(provenance.claimValidity, .bindingPhysical)
    }

    func testLowercasePrefixAcceptsUppercaseHexDigest() {
        let uppercaseEnergyReceipt = "sha256:" + energyReceipt.dropFirst(7).uppercased()
        let provenance = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .boundOnly,
            energyProvenance: uppercaseEnergyReceipt,
            measureProvenance: measureReceipt
        )

        XCTAssertEqual(provenance.claimValidity, .canonicalPhysical)
    }

    func testWrongSchemaVersionFailsClosed() {
        let provenance = ScientificProvenance(
            schemaVersion: 1,
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            referenceState: .matchedAssociationCycle,
            energyProvenance: energyReceipt,
            measureProvenance: measureReceipt,
            referenceProvenance: referenceReceipt
        )
        XCTAssertEqual(provenance.claimValidity, .proxyOnly)
    }

    func testDescriptiveAndMalformedEvidenceFailsClosed() {
        let descriptiveOnly = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            energyProvenance: "calibrated Hamiltonian receipt",
            measureProvenance: measureReceipt
        )
        XCTAssertEqual(descriptiveOnly.claimValidity, .proxyOnly)

        let shortDigest = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .enumeratedMicrostates,
            energyProvenance: "sha256:abc123",
            measureProvenance: measureReceipt
        )
        XCTAssertEqual(shortDigest.claimValidity, .proxyOnly)
    }

    func testPlaceholderAndLowDiversityDigestsFailClosed() {
        let invalidDigests = [
            "sha256:" + String(repeating: "0", count: 64),
            "sha256:" + String(repeating: "ab", count: 32),
            "SHA256:e638aaee2a68410cdc827397b2aa095cf227090f494f535748c756ac49e6da3c",
            "sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a",
            "sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
        ]

        for digest in invalidDigests {
            let provenance = ScientificProvenance(
                energyDomain: .calibratedKcalPerMol,
                ensembleMeasure: .enumeratedMicrostates,
                energyProvenance: digest,
                measureProvenance: measureReceipt
            )
            XCTAssertEqual(provenance.claimValidity, .proxyOnly, digest)
        }
    }

    func testSchemaV2EncodingUsesCanonicalSnakeCaseAndDerivedValidity() throws {
        let provenance = ScientificProvenance(
            energyDomain: .calibratedKcalPerMol,
            ensembleMeasure: .weightedQuadrature,
            referenceState: .matchedAssociationCycle,
            energyProvenance: energyReceipt,
            measureProvenance: measureReceipt,
            referenceProvenance: referenceReceipt
        )
        let result = ThermodynamicResult(
            temperature: 300, logZ: 1, freeEnergy: -1,
            meanEnergy: -0.5, meanEnergySq: 0.25,
            heatCapacity: 0, entropy: 0, stdEnergy: 0,
            scientificProvenance: provenance
        )

        let data = try JSONEncoder().encode(result)
        let root = try XCTUnwrap(JSONSerialization.jsonObject(with: data) as? [String: Any])
        let encoded = try XCTUnwrap(root["scientificProvenance"] as? [String: Any])

        XCTAssertEqual(encoded["schema_version"] as? Int, 2)
        XCTAssertEqual(encoded["energy_domain"] as? String, "calibrated_kcal_per_mol")
        XCTAssertEqual(encoded["ensemble_measure"] as? String, "weighted_quadrature")
        XCTAssertEqual(encoded["reference_state"] as? String, "matched_association_cycle")
        XCTAssertEqual(encoded["energy_provenance"] as? String, energyReceipt)
        XCTAssertEqual(encoded["measure_provenance"] as? String, measureReceipt)
        XCTAssertEqual(encoded["reference_provenance"] as? String, referenceReceipt)
        XCTAssertEqual(encoded["claim_validity"] as? String, "binding_physical")
        XCTAssertNil(encoded["schemaVersion"])
        XCTAssertNil(encoded["claimValidity"])
    }

    func testSerializedClaimValidityCannotOverrideDerivedValidity() throws {
        let hostile: [String: Any] = [
            "schema_version": 2,
            "energy_domain": "unclassified",
            "ensemble_measure": "unclassified",
            "reference_state": "none",
            "energy_provenance": energyReceipt,
            "measure_provenance": measureReceipt,
            "reference_provenance": referenceReceipt,
            "claim_validity": "binding_physical",
        ]
        let decoded = try JSONDecoder().decode(
            ScientificProvenance.self,
            from: JSONSerialization.data(withJSONObject: hostile)
        )
        XCTAssertEqual(decoded.claimValidity, .proxyOnly)

        var valid = hostile
        valid["energy_domain"] = "calibrated_kcal_per_mol"
        valid["ensemble_measure"] = "enumerated_microstates"
        valid["reference_state"] = "matched_association_cycle"
        valid["claim_validity"] = "proxy_only"
        let validDecoded = try JSONDecoder().decode(
            ScientificProvenance.self,
            from: JSONSerialization.data(withJSONObject: valid)
        )
        XCTAssertEqual(validDecoded.claimValidity, .bindingPhysical)
    }

    func testMissingAndHostileCanonicalFieldsFailClosed() throws {
        let missingSchema: [String: Any] = [
            "energy_domain": "calibrated_kcal_per_mol",
            "ensemble_measure": "enumerated_microstates",
            "reference_state": "matched_association_cycle",
            "energy_provenance": energyReceipt,
            "measure_provenance": measureReceipt,
            "reference_provenance": referenceReceipt,
            "claim_validity": "binding_physical",
        ]
        let missingDecoded = try JSONDecoder().decode(
            ScientificProvenance.self,
            from: JSONSerialization.data(withJSONObject: missingSchema)
        )
        XCTAssertEqual(missingDecoded.claimValidity, .proxyOnly)

        var hostile = missingSchema
        hostile["schema_version"] = 2.5
        hostile["energy_domain"] = ["calibrated_kcal_per_mol"]
        hostile["energy_provenance"] = ["receipt": energyReceipt]
        let hostileDecoded = try JSONDecoder().decode(
            ScientificProvenance.self,
            from: JSONSerialization.data(withJSONObject: hostile)
        )
        XCTAssertEqual(hostileDecoded.schemaVersion, 0)
        XCTAssertEqual(hostileDecoded.energyDomain, .unclassified)
        XCTAssertEqual(hostileDecoded.energyProvenance, "")
        XCTAssertEqual(hostileDecoded.claimValidity, .proxyOnly)
    }

    func testLegacyCamelCaseDecodesButReencodesCanonically() throws {
        let legacy: [String: Any] = [
            "schemaVersion": 2,
            "energyDomain": "calibrated_kcal_per_mol",
            "ensembleMeasure": "weighted_quadrature",
            "referenceState": "matched_association_cycle",
            "energyProvenance": energyReceipt,
            "measureProvenance": measureReceipt,
            "referenceProvenance": referenceReceipt,
            "claimValidity": "proxy_only",
        ]
        let decoded = try JSONDecoder().decode(
            ScientificProvenance.self,
            from: JSONSerialization.data(withJSONObject: legacy)
        )
        XCTAssertEqual(decoded.claimValidity, .bindingPhysical)

        let encodedData = try JSONEncoder().encode(decoded)
        let encoded = try XCTUnwrap(
            JSONSerialization.jsonObject(with: encodedData) as? [String: Any]
        )
        XCTAssertEqual(encoded["claim_validity"] as? String, "binding_physical")
        XCTAssertNil(encoded["claimValidity"])
        XCTAssertNil(encoded["energyProvenance"])
    }

    func testHostileCanonicalKeyDoesNotFallBackToLegacy() throws {
        let mixed: [String: Any] = [
            "schema_version": "2",
            "schemaVersion": 2,
            "energy_domain": ["calibrated_kcal_per_mol"],
            "energyDomain": "calibrated_kcal_per_mol",
            "ensembleMeasure": "enumerated_microstates",
            "referenceState": "matched_association_cycle",
            "energyProvenance": energyReceipt,
            "measureProvenance": measureReceipt,
            "referenceProvenance": referenceReceipt,
            "claimValidity": "binding_physical",
        ]
        let decoded = try JSONDecoder().decode(
            ScientificProvenance.self,
            from: JSONSerialization.data(withJSONObject: mixed)
        )
        XCTAssertEqual(decoded.schemaVersion, 0)
        XCTAssertEqual(decoded.energyDomain, .unclassified)
        XCTAssertEqual(decoded.claimValidity, .proxyOnly)
    }
}
