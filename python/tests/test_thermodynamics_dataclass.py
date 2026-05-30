"""Tests for the pure-Python Thermodynamics dataclass.

No C++ extension required – these run in every environment.
"""

from __future__ import annotations

import pytest

from flexaidds.thermodynamics import Thermodynamics


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
            "heat_capacity_kcal_mol_K2", "std_energy_kcal_mol",
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
        assert d["heat_capacity_kcal_mol_K2"] == t.heat_capacity
        assert d["std_energy_kcal_mol"] == t.std_energy

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
        original = _make()
        d = original.to_dict()
        # to_dict omits mean_energy_sq, so supply it manually
        d["mean_energy_sq"] = original.mean_energy_sq
        restored = Thermodynamics.from_dict(d)
        assert restored.temperature == pytest.approx(original.temperature)
        assert restored.log_Z == pytest.approx(original.log_Z)
        assert restored.free_energy == pytest.approx(original.free_energy)
        assert restored.mean_energy == pytest.approx(original.mean_energy)
        assert restored.mean_energy_sq == pytest.approx(original.mean_energy_sq)
        assert restored.heat_capacity == pytest.approx(original.heat_capacity)
        assert restored.entropy == pytest.approx(original.entropy)
        assert restored.std_energy == pytest.approx(original.std_energy)

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
            heat_capacity_kcal_mol_K2=0.05,
            entropy_kcal_mol_K=0.002,
            std_energy_kcal_mol=0.3,
        )
        t = Thermodynamics.from_dict(data)
        assert t.temperature == pytest.approx(300.0)
        assert t.free_energy == pytest.approx(-10.0)

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
