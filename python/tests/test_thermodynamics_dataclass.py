"""Tests for the pure-Python Thermodynamics dataclass.

No C++ extension required – these run in every environment.
"""

from __future__ import annotations

import pytest

from flexaidds.thermodynamics import (
    ClaimValidity,
    EnergyDomain,
    EnsembleMeasure,
    ReferenceState,
    ScientificProvenance,
    StatMechEngine,
    ThermodynamicBreakdown,
    Thermodynamics,
)

ENERGY_RECEIPT = "sha256:6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b"
MEASURE_RECEIPT = "sha256:d4735e3a265e16eee03f59718b9b5d03019c07d8b6c51f90da3a666eec13ab35"
REFERENCE_RECEIPT = "sha256:4e07408562bedb8b60ce05c1decfe3ad16b72230967de01f640b7e4729b49fce"


def _make(**overrides) -> Thermodynamics:
    defaults = dict(
        temperature=300.0,
        log_Z=-16.8,
        free_energy=-10.0,
        mean_energy=-9.5,
        mean_energy_sq=90.5,
        heat_capacity=0.05,
        entropy=0.001667,
        std_energy=0.3,
    )
    defaults.update(overrides)
    return Thermodynamics(**defaults)


def _physical_provenance(**overrides) -> ScientificProvenance:
    defaults = dict(
        schema_version=2,
        energy_domain=EnergyDomain.CALIBRATED_KCAL_PER_MOL,
        ensemble_measure=EnsembleMeasure.ENUMERATED_MICROSTATES,
        reference_state=ReferenceState.MATCHED_ASSOCIATION_CYCLE,
        energy_provenance=ENERGY_RECEIPT,
        measure_provenance=MEASURE_RECEIPT,
        reference_provenance=REFERENCE_RECEIPT,
    )
    defaults.update(overrides)
    return ScientificProvenance(**defaults)


class TestThermodynamicsDataclass:
    def test_binding_free_energy_is_alias(self):
        t = _make(free_energy=-10.0)
        assert t.binding_free_energy == t.free_energy

    def test_entropy_term_is_T_times_S(self):
        T, S = 300.0, 0.002
        t = _make(temperature=T, entropy=S)
        assert t.entropy_term == pytest.approx(T * S)

    def test_to_dict_has_expected_keys(self):
        d = _make().to_dict()
        expected = {
            "temperature_K", "log_Z", "free_energy_kcal_mol",
            "enthalpy_kcal_mol", "mean_energy_sq",
            "entropy_kcal_mol_K",
            "heat_capacity_kcal_mol_K", "std_energy_kcal_mol",
            "scientific_provenance",
        }
        assert expected == set(d.keys())

    def test_to_dict_values_match_fields(self):
        t = _make()
        d = t.to_dict()
        assert d["temperature_K"] == t.temperature
        assert d["free_energy_kcal_mol"] == t.free_energy
        assert d["enthalpy_kcal_mol"] == t.mean_energy
        assert d["mean_energy_sq"] == t.mean_energy_sq
        assert d["entropy_kcal_mol_K"] == t.entropy
        assert d["heat_capacity_kcal_mol_K"] == t.heat_capacity
        assert "heat_capacity_kcal_mol_K2" not in d
        assert d["std_energy_kcal_mol"] == t.std_energy

    def test_default_provenance_is_explicitly_proxy_only(self):
        t = _make()
        metadata = t.to_dict()["scientific_provenance"]

        assert metadata == {
            "schema_version": 2,
            "energy_domain": "unclassified",
            "ensemble_measure": "unclassified",
            "reference_state": "none",
            "energy_provenance": "",
            "measure_provenance": "",
            "reference_provenance": "",
            "claim_validity": "proxy_only",
        }
        assert t.claim_validity is ClaimValidity.PROXY_ONLY
        assert t.is_proxy_only()
        assert not t.allows_canonical_claims()
        assert not t.allows_binding_claims()

    def test_entropy_term_zero_when_entropy_zero(self):
        t = _make(entropy=0.0)
        assert t.entropy_term == 0.0

    def test_fields_are_readable(self):
        """All declared fields are readable (basic sanity)."""
        t = _make()
        for attr in ("temperature", "log_Z", "free_energy", "mean_energy",
                     "mean_energy_sq", "heat_capacity", "entropy", "std_energy"):
            assert hasattr(t, attr)


class TestThermodynamicsFromDict:
    """Tests for Thermodynamics.from_dict() round-trip deserialization."""

    def test_round_trip_via_to_dict(self):
        """from_dict(to_dict()) reproduces all scalar fields."""
        original = _make(provenance=_physical_provenance())
        d = original.to_dict()
        restored = Thermodynamics.from_dict(d)
        assert restored.temperature == pytest.approx(original.temperature)
        assert restored.log_Z == pytest.approx(original.log_Z)
        assert restored.free_energy == pytest.approx(original.free_energy)
        assert restored.mean_energy == pytest.approx(original.mean_energy)
        assert restored.mean_energy_sq == pytest.approx(original.mean_energy_sq)
        assert restored.heat_capacity == pytest.approx(original.heat_capacity)
        assert restored.entropy == pytest.approx(original.entropy)
        assert restored.std_energy == pytest.approx(original.std_energy)
        assert restored.provenance == original.provenance
        assert restored.claim_validity is ClaimValidity.BINDING_PHYSICAL

    def test_accepts_raw_attribute_names(self):
        """from_dict works with raw attribute names (no unit suffixes)."""
        data = dict(
            temperature=310.0,
            log_Z=-15.0,
            free_energy=-8.0,
            mean_energy=-7.5,
            mean_energy_sq=57.0,
            heat_capacity=0.04,
            entropy=0.0016,
            std_energy=0.25,
        )
        t = Thermodynamics.from_dict(data)
        assert t.temperature == pytest.approx(310.0)
        assert t.free_energy == pytest.approx(-8.0)
        assert t.mean_energy == pytest.approx(-7.5)

    def test_suffixed_keys_take_priority(self):
        """When both suffixed and raw keys exist, suffixed wins."""
        data = dict(
            temperature_K=300.0,
            temperature=999.0,  # should be ignored
            log_Z=-16.0,
            free_energy_kcal_mol=-10.0,
            free_energy=-999.0,
            enthalpy_kcal_mol=-9.5,
            mean_energy_sq=90.0,
            heat_capacity_kcal_mol_K=0.05,
            heat_capacity_kcal_mol_K2=999.0,
            heat_capacity=888.0,
            entropy_kcal_mol_K=0.002,
            std_energy_kcal_mol=0.3,
        )
        t = Thermodynamics.from_dict(data)
        assert t.temperature == pytest.approx(300.0)
        assert t.free_energy == pytest.approx(-10.0)
        assert t.heat_capacity == pytest.approx(0.05)

    def test_accepts_legacy_heat_capacity_K2_key(self):
        """The historical, dimensionally incorrect K2 key remains readable."""
        data = _make().to_dict()
        value = data.pop("heat_capacity_kcal_mol_K")
        data["heat_capacity_kcal_mol_K2"] = value

        restored = Thermodynamics.from_dict(data)

        assert restored.heat_capacity == pytest.approx(value)

    def test_missing_key_raises_key_error(self):
        """from_dict raises KeyError when a required field is absent."""
        data = dict(temperature_K=300.0)  # missing everything else
        with pytest.raises(KeyError):
            Thermodynamics.from_dict(data)

    def test_derived_properties_after_from_dict(self):
        """Derived properties (entropy_term, binding_free_energy) work."""
        data = dict(
            temperature=300.0,
            log_Z=-16.0,
            free_energy=-10.0,
            mean_energy=-9.5,
            mean_energy_sq=90.0,
            heat_capacity=0.05,
            entropy=0.002,
            std_energy=0.3,
        )
        t = Thermodynamics.from_dict(data)
        assert t.binding_free_energy == pytest.approx(-10.0)
        assert t.entropy_term == pytest.approx(300.0 * 0.002)


class TestScientificClaimGates:
    @pytest.mark.parametrize(
        "ensemble_measure",
        [
            EnsembleMeasure.ENUMERATED_MICROSTATES,
            EnsembleMeasure.WEIGHTED_QUADRATURE,
        ],
    )
    def test_canonical_claim_requires_calibrated_domain_measure_and_provenance(
        self, ensemble_measure
    ):
        provenance = _physical_provenance(
            ensemble_measure=ensemble_measure,
            reference_state=ReferenceState.NONE,
            reference_provenance="",
        )

        assert provenance.allows_canonical_claims()
        assert provenance.allows_canonical_physical_claim()
        assert not provenance.allows_binding_claims()
        assert provenance.claim_validity is ClaimValidity.CANONICAL_PHYSICAL

    @pytest.mark.parametrize(
        "override",
        [
            {"schema_version": 1},
            {"energy_domain": EnergyDomain.UNCLASSIFIED},
            {"energy_domain": EnergyDomain.CF_ARBITRARY_UNITS},
            {"energy_domain": EnergyDomain.MODEL_SCALE},
            {"ensemble_measure": EnsembleMeasure.UNCLASSIFIED},
            {"ensemble_measure": EnsembleMeasure.OPTIMIZER_SAMPLES},
            {"energy_provenance": ""},
            {"energy_provenance": "   "},
            {"energy_provenance": "\u00a0"},
            {"energy_provenance": "校准"},
            {"energy_provenance": "receipt"},
            {"energy_provenance": "sha256:" + "0" * 64},
            {"energy_provenance": "sha256:" + "ab" * 32},
            {"energy_provenance": "sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
            {"energy_provenance": "sha256:3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a"},
            {"measure_provenance": ""},
            {"measure_provenance": "\u2003"},
        ],
    )
    def test_canonical_claim_fails_closed_when_evidence_is_insufficient(
        self, override
    ):
        provenance = _physical_provenance(**override)

        assert not provenance.allows_canonical_claims()
        assert not provenance.allows_binding_claims()
        assert provenance.claim_validity is ClaimValidity.PROXY_ONLY

    def test_binding_claim_requires_matched_cycle_and_reference_provenance(self):
        binding = _physical_provenance()
        bound_only = _physical_provenance(
            reference_state=ReferenceState.BOUND_ONLY
        )
        missing_reference = _physical_provenance(reference_provenance="")

        assert binding.allows_binding_claims()
        assert binding.allows_binding_physical_claim()
        assert binding.claim_validity is ClaimValidity.BINDING_PHYSICAL
        assert not bound_only.allows_binding_claims()
        assert bound_only.claim_validity is ClaimValidity.CANONICAL_PHYSICAL
        assert not missing_reference.allows_binding_claims()
        assert missing_reference.claim_validity is ClaimValidity.CANONICAL_PHYSICAL

    def test_binding_claim_rejects_unicode_only_reference_provenance(self):
        for reference_provenance in ("\u00a0", "\u2003", "参照"):
            provenance = _physical_provenance(
                reference_provenance=reference_provenance
            )
            assert not provenance.allows_binding_claims()
            assert provenance.claim_validity is ClaimValidity.CANONICAL_PHYSICAL

    def test_only_structured_sha256_receipts_authorize_physical_claims(self):
        provenance = _physical_provenance(
            energy_provenance=ENERGY_RECEIPT.upper().replace("SHA256:", "sha256:"),
            measure_provenance=MEASURE_RECEIPT,
            reference_provenance=REFERENCE_RECEIPT,
        )

        assert provenance.allows_canonical_claims()
        assert provenance.allows_binding_claims()
        assert provenance.claim_validity is ClaimValidity.BINDING_PHYSICAL

    def test_breakdown_round_trip_preserves_provenance_and_claim_gate(self):
        provenance = _physical_provenance()
        breakdown = ThermodynamicBreakdown.from_thermodynamics(
            _make(provenance=provenance)
        )

        restored = ThermodynamicBreakdown.from_dict(breakdown.to_dict())

        assert restored.provenance == provenance
        assert restored.allows_binding_claims()
        assert restored.claim_validity is ClaimValidity.BINDING_PHYSICAL

    def test_serialized_claim_validity_override_is_ignored(self):
        payload = _make().to_dict()
        payload["scientific_provenance"]["claim_validity"] = "binding_physical"

        restored = Thermodynamics.from_dict(payload)

        assert restored.claim_validity is ClaimValidity.PROXY_ONLY
        assert (
            restored.to_dict()["scientific_provenance"]["claim_validity"]
            == "proxy_only"
        )

    def test_missing_metadata_schema_version_fails_closed(self):
        payload = _make().to_dict()
        payload["scientific_provenance"] = _physical_provenance().to_dict()
        payload["scientific_provenance"].pop("schema_version")

        restored = Thermodynamics.from_dict(payload)

        assert restored.claim_validity is ClaimValidity.PROXY_ONLY

    def test_unknown_metadata_values_fail_closed(self):
        payload = _make().to_dict()
        payload["scientific_provenance"].update(
            schema_version=999,
            energy_domain="experimental_free_energy",
            ensemble_measure="magic_sampler",
            reference_state="assumed_standard_state",
            energy_provenance="asserted",
            measure_provenance="asserted",
            reference_provenance="asserted",
            claim_validity="binding_physical",
        )

        restored = Thermodynamics.from_dict(payload)

        assert restored.provenance.energy_domain is EnergyDomain.UNCLASSIFIED
        assert restored.provenance.ensemble_measure is EnsembleMeasure.UNCLASSIFIED
        assert restored.provenance.reference_state is ReferenceState.NONE
        assert restored.claim_validity is ClaimValidity.PROXY_ONLY

    @pytest.mark.parametrize("schema_version", [2.0, 2.9, "2", True, [2]])
    def test_schema_version_coercions_cannot_authorize_claims(self, schema_version):
        payload = _make().to_dict()
        payload["scientific_provenance"] = _physical_provenance().to_dict()
        payload["scientific_provenance"]["schema_version"] = schema_version

        restored = Thermodynamics.from_dict(payload)

        assert restored.provenance.schema_version == 0
        assert restored.claim_validity is ClaimValidity.PROXY_ONLY

    @pytest.mark.parametrize("bad_evidence", [1, True, ["receipt"], {"id": "x"}])
    def test_non_string_evidence_cannot_authorize_claims(self, bad_evidence):
        payload = _make().to_dict()
        payload["scientific_provenance"] = _physical_provenance().to_dict()
        payload["scientific_provenance"]["energy_provenance"] = bad_evidence

        restored = Thermodynamics.from_dict(payload)

        assert restored.provenance.energy_provenance == ""
        assert restored.claim_validity is ClaimValidity.PROXY_ONLY

    def test_flat_metadata_is_accepted_but_new_output_is_nested(self):
        payload = _make().to_dict()
        payload.pop("scientific_provenance")
        payload.update(_physical_provenance().to_dict())

        restored = Thermodynamics.from_dict(payload)
        serialized = restored.to_dict()

        assert restored.claim_validity is ClaimValidity.BINDING_PHYSICAL
        assert "scientific_provenance" in serialized
        assert "energy_domain" not in serialized

    def test_engine_provenance_changes_only_claim_metadata(self):
        default_engine = StatMechEngine(300.0)
        physical_engine = StatMechEngine(300.0, _physical_provenance())
        for energy in (-10.0, -9.0):
            default_engine.add_sample(energy)
            physical_engine.add_sample(energy)

        default_result = default_engine.compute()
        physical_result = physical_engine.compute()

        for field in (
            "temperature",
            "log_Z",
            "free_energy",
            "mean_energy",
            "mean_energy_sq",
            "heat_capacity",
            "entropy",
            "std_energy",
        ):
            assert getattr(physical_result, field) == pytest.approx(
                getattr(default_result, field)
            )
        assert default_result.claim_validity is ClaimValidity.PROXY_ONLY
        assert physical_result.claim_validity is ClaimValidity.BINDING_PHYSICAL
