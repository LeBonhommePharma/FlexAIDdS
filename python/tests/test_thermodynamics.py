"""Tests for flexaidds.thermodynamics — Python-only paths + C++ integration.

These tests cover the pure-Python numerical logic and the module's public API.
Tests that require the compiled C++ extension are marked with needs_core / @requires_core.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

# ── C++ availability guard ────────────────────────────────────────────────────

def _real_core_available():
    """Check if the real C++ _core (not the test stub) is usable."""
    try:
        import flexaidds._core as c
        c.StatMechEngine(300.0)
        return True
    except Exception:
        return False

_CORE_AVAILABLE = _real_core_available()

needs_core = pytest.mark.skipif(not _CORE_AVAILABLE,
                                reason="C++ _core module not built")

try:
    from conftest import requires_core
except ImportError:
    requires_core = needs_core


# ── helpers ───────────────────────────────────────────────────────────────────

kB = 0.001987206  # kcal mol^-1 K^-1


def _boltzmann_weights(energies, temperature=300.0):
    """Reference implementation of Boltzmann weighting."""
    beta = 1.0 / (kB * temperature)
    shifted = [e - min(energies) for e in energies]
    weights = [math.exp(-beta * e) for e in shifted]
    z = sum(weights)
    return [w / z for w in weights]


def _helmholtz(energies, temperature=300.0):
    """Reference Helmholtz free energy via log-sum-exp."""
    beta = 1.0 / (kB * temperature)
    e_min = min(energies)
    log_z = math.log(sum(math.exp(-beta * (e - e_min)) for e in energies)) - beta * e_min
    return -kB * temperature * log_z


# ─────────────────────────────────────────────────────────────────────────────
# Module-level imports (no C++ needed)
# ─────────────────────────────────────────────────────────────────────────────

def test_import_thermodynamics():
    """The thermodynamics module must import without the C++ extension."""
    import importlib
    mod = importlib.import_module("flexaidds.thermodynamics")
    assert hasattr(mod, "Thermodynamics")
    assert hasattr(mod, "StatMechEngine")
    assert hasattr(mod, "BoltzmannLUT")
    assert hasattr(mod, "helmholtz_from_energies")
    assert hasattr(mod, "kB_kcal")


def test_kB_kcal_value():
    """kB_kcal must be within 0.1% of the NIST value."""
    from flexaidds.thermodynamics import kB_kcal
    assert abs(kB_kcal - 0.001987206) < 1e-8


def test_package_exports_stable_thermodynamics_facade():
    """Top-level thermodynamics types must not depend on extension presence."""
    import flexaidds
    from flexaidds import thermodynamics as thermo_module

    for name in (
        "StatMechEngine",
        "Thermodynamics",
        "ThermodynamicBreakdown",
        "EnergyDomain",
        "EnsembleMeasure",
        "ReferenceState",
        "ClaimValidity",
        "ScientificProvenance",
    ):
        assert getattr(flexaidds, name) is getattr(thermo_module, name)


@needs_core
def test_public_facade_delegates_to_native_engine_with_provenance():
    """Compiled installs retain the Python schema while using native math."""
    import flexaidds
    from flexaidds import _core

    provenance = flexaidds.ScientificProvenance(
        energy_domain=flexaidds.EnergyDomain.CALIBRATED_KCAL_PER_MOL,
        ensemble_measure=flexaidds.EnsembleMeasure.ENUMERATED_MICROSTATES,
        energy_provenance="sha256:6b86b273ff34fce19d6b804eff5a3f5747ada4eaa22f1d49c01e52ddb7875b4b",
        measure_provenance="sha256:d4735e3a265e16eee03f59718b9b5d03019c07d8b6c51f90da3a666eec13ab35",
    )
    engine = flexaidds.StatMechEngine(300.0, provenance)
    engine.add_sample(-10.0)
    result = engine.compute()

    assert isinstance(engine._engine, _core.StatMechEngine)
    assert isinstance(result, flexaidds.Thermodynamics)
    assert result.allows_canonical_claims()
    assert result.to_dict()["scientific_provenance"]["schema_version"] == 2


# ─────────────────────────────────────────────────────────────────────────────
# Thermodynamics dataclass
# ─────────────────────────────────────────────────────────────────────────────

class TestThermodynamicsDataclass:
    def _make(self, **kw):
        from flexaidds.thermodynamics import Thermodynamics
        defaults = dict(
            temperature=300.0, log_Z=0.0, free_energy=-10.0,
            mean_energy=-9.0, mean_energy_sq=82.0,
            heat_capacity=0.05, entropy=0.003, std_energy=0.5,
        )
        defaults.update(kw)
        return Thermodynamics(**defaults)

    def test_binding_free_energy_alias(self):
        t = self._make(free_energy=-12.3)
        assert t.binding_free_energy == t.free_energy == -12.3

    def test_entropy_term(self):
        t = self._make(temperature=300.0, entropy=0.01)
        assert abs(t.entropy_term - 3.0) < 1e-9

    def test_to_dict_keys(self):
        t = self._make()
        d = t.to_dict()
        for key in ("temperature_K", "free_energy_kcal_mol",
                    "entropy_kcal_mol_K", "heat_capacity_kcal_mol_K",
                    "scientific_provenance"):
            assert key in d
        assert "heat_capacity_kcal_mol_K2" not in d

    def test_to_dict_values(self):
        t = self._make(temperature=310.0, free_energy=-11.5)
        d = t.to_dict()
        assert d["temperature_K"] == 310.0
        assert d["free_energy_kcal_mol"] == -11.5


def test_thermodynamics_dataclass_legacy():
    """Thermodynamics dataclass properties and serialisation (legacy test)."""
    from flexaidds.thermodynamics import Thermodynamics

    thermo = Thermodynamics(
        temperature=300.0,
        log_Z=10.0,
        free_energy=-5.0,
        mean_energy=-4.5,
        mean_energy_sq=20.5,
        heat_capacity=0.1,
        entropy=0.001667,
        std_energy=0.3,
    )

    assert thermo.binding_free_energy == thermo.free_energy
    assert math.isclose(thermo.entropy_term, 300.0 * 0.001667, rel_tol=1e-9)

    d = thermo.to_dict()
    assert d["temperature_K"] == 300.0
    assert d["free_energy_kcal_mol"] == -5.0
    assert d["entropy_kcal_mol_K"] == 0.001667


# ─────────────────────────────────────────────────────────────────────────────
# Tests that require the compiled C++ extension
# ─────────────────────────────────────────────────────────────────────────────

@needs_core
class TestStatMechEngine:
    def test_empty_engine(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(300.0)
        assert eng.n_samples == 0
        assert abs(eng.temperature - 300.0) < 1e-9

    def test_add_samples_length(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(300.0)
        energies = [-10.5, -9.8, -10.2, -11.0]
        eng.add_samples(np.array(energies))
        assert eng.n_samples == len(energies)

    def test_free_energy_single_sample(self):
        """With one sample, F = E (partition function = 1)."""
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(300.0)
        eng.add_sample(-10.0)
        thermo = eng.compute()
        assert abs(thermo.free_energy - (-10.0)) < 1e-6

    def test_boltzmann_weights_sum_to_one(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(300.0)
        eng.add_samples(np.array([-10.0, -9.0, -11.0]))
        weights = eng.boltzmann_weights()
        assert abs(sum(weights) - 1.0) < 1e-9

    def test_clear_resets_count(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(300.0)
        eng.add_sample(-8.0)
        eng.clear()
        assert eng.n_samples == 0

    def test_delta_g_same_ensemble_is_zero(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng1 = StatMechEngine(300.0)
        eng2 = StatMechEngine(300.0)
        for e in [-10.0, -9.5, -10.5]:
            eng1.add_sample(e)
            eng2.add_sample(e)
        assert abs(eng1.delta_G(eng2)) < 1e-9

    def test_thermodynamics_entropy_nonnegative(self):
        """Configurational entropy should be >= 0 for a spread ensemble."""
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(300.0)
        eng.add_samples(np.linspace(-12.0, -8.0, 20))
        thermo = eng.compute()
        assert thermo.entropy >= 0.0


@needs_core
def test_statmech_engine_basic(sample_energies):
    """StatMechEngine returns sensible thermodynamics for a small ensemble."""
    from flexaidds.thermodynamics import StatMechEngine

    engine = StatMechEngine(temperature_K=300.0)
    engine.add_samples(sample_energies)

    thermo = engine.compute()

    # Free energy must be <= min(energies) at any temperature
    assert thermo.free_energy <= min(sample_energies) + 1e-9

    # Entropy must be non-negative
    assert thermo.entropy >= 0.0

    # Heat capacity must be non-negative
    assert thermo.heat_capacity >= 0.0

    # Mean energy must be between min and max
    assert min(sample_energies) - 1e-9 <= thermo.mean_energy <= max(sample_energies) + 1e-9


@needs_core
def test_statmech_boltzmann_weights_sum_to_one(sample_energies):
    """Boltzmann weights must sum to 1.0."""
    from flexaidds.thermodynamics import StatMechEngine

    engine = StatMechEngine(temperature_K=300.0)
    engine.add_samples(sample_energies)

    weights = engine.boltzmann_weights()
    assert abs(sum(weights) - 1.0) < 1e-10


@needs_core
def test_statmech_single_sample():
    """With one sample, F = E and S = 0."""
    from flexaidds.thermodynamics import StatMechEngine

    engine = StatMechEngine(temperature_K=300.0)
    engine.add_sample(-10.0)

    thermo = engine.compute()
    assert math.isclose(thermo.free_energy, -10.0, abs_tol=1e-9)
    assert abs(thermo.entropy) < 1e-10


@needs_core
def test_statmech_clear():
    """clear() resets the ensemble."""
    from flexaidds.thermodynamics import StatMechEngine

    engine = StatMechEngine(300.0)
    engine.add_sample(-10.0)
    assert len(engine) == 1
    engine.clear()
    assert len(engine) == 0


@needs_core
def test_statmech_delta_G():
    """delta_G computes a reasonable relative free energy."""
    from flexaidds.thermodynamics import StatMechEngine

    eng1 = StatMechEngine(300.0)
    for e in [-12.0, -11.5]:
        eng1.add_sample(e)

    eng2 = StatMechEngine(300.0)
    for e in [-10.0, -9.5]:
        eng2.add_sample(e)

    dG = eng1.delta_G(eng2)
    # eng1 has lower energies -> more negative dG
    assert dG < 0.0


# ── BoltzmannLUT ──────────────────────────────────────────────────────────────

@needs_core
class TestBoltzmannLUT:
    def test_lookup_matches_exact(self):
        from flexaidds.thermodynamics import BoltzmannLUT, kB_kcal
        beta = 1.0 / (kB_kcal * 300.0)
        lut = BoltzmannLUT(beta=beta, e_min=-20.0, e_max=5.0, n_bins=100000)
        energy = -12.5
        exact = math.exp(-beta * energy)
        approx = lut(energy)
        # LUT resolution: within 0.05% of exact value
        # (nearest-bin lookup with 100K bins over 25 kcal/mol range)
        assert abs(approx - exact) / exact < 5e-4


@needs_core
def test_boltzmann_lut_legacy():
    """BoltzmannLUT result approximates direct exp(-beta*E) within 1%."""
    from flexaidds.thermodynamics import BoltzmannLUT, kB_kcal

    T = 300.0
    beta = 1.0 / (kB_kcal * T)
    lut = BoltzmannLUT(beta=beta, e_min=-20.0, e_max=5.0, n_bins=10000)

    for E in [-15.0, -10.0, -5.0, 0.0]:
        expected = math.exp(-beta * E)
        got = lut(E)
        assert abs(got - expected) / (expected + 1e-30) < 0.01


@needs_core
def test_helmholtz_from_energies():
    """helmholtz_from_energies matches StatMechEngine.compute()."""
    from flexaidds.thermodynamics import StatMechEngine, helmholtz_from_energies

    energies = [-12.5, -11.8, -12.1, -10.9]
    T = 300.0

    engine = StatMechEngine(T)
    engine.add_samples(energies)
    expected_F = engine.compute().free_energy

    got_F = helmholtz_from_energies(energies, T)
    assert math.isclose(expected_F, got_F, rel_tol=1e-9)
