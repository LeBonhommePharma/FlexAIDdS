"""Tests for StatMechEngine and related thermodynamic utilities.

These tests verify *correctness* of the computed values against analytical
expectations derived from canonical statistical mechanics.  They rely on the
compiled C++ ``_core`` extension; all tests are automatically skipped when the
extension is not built (use ``pip install -e python/`` to build it).

Physical definitions used throughout
--------------------------------------
  β = 1 / (kB T)          (kB = 0.001987206 kcal mol⁻¹ K⁻¹)
  Z = Σ exp(−β Eᵢ)        (partition function)
  F = −kT ln Z             (Helmholtz free energy)
  ⟨E⟩ = Σ wᵢ Eᵢ           (Boltzmann-weighted mean)  where wᵢ = exp(−β Eᵢ)/Z
  S = (⟨E⟩ − F) / T       (configurational entropy)
  Cv = (⟨E²⟩ − ⟨E⟩²)/(kT²)
  σ_E = sqrt(Cv (kT)²)    = sqrt(⟨E²⟩ − ⟨E⟩²)
"""

from __future__ import annotations

import math

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Skip the whole module when the C++ extension is not available.
# The conftest.py skip_without_core marker is re-used here as a module-level
# pytestmark so every test inherits it automatically.
# ---------------------------------------------------------------------------

_core_available: bool
try:
    import flexaidds._core as _core_mod
    _core_mod.StatMechEngine(300.0)   # raises if it's the stub
    _core_available = True
except Exception:
    _core_available = False

pytestmark = pytest.mark.skipif(
    not _core_available,
    reason="C++ _core extension not built",
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

KB = 0.001987206  # kcal mol⁻¹ K⁻¹


def _make_engine(temperature: float, energies):
    """Return a StatMechEngine loaded with the given energies."""
    from flexaidds.thermodynamics import StatMechEngine
    eng = StatMechEngine(temperature_K=temperature)
    for e in energies:
        eng.add_sample(float(e))
    return eng


def _analytical(temperature: float, energies: list[float]) -> dict:
    """Compute ground-truth thermodynamics from first principles."""
    beta = 1.0 / (KB * temperature)
    # Numerically stable log-sum-exp
    x = [-beta * e for e in energies]
    max_x = max(x)
    log_Z = max_x + math.log(sum(math.exp(xi - max_x) for xi in x))
    Z = math.exp(log_Z)

    weights = [math.exp(xi - log_Z) for xi in x]  # normalised
    mean_E = sum(w * e for w, e in zip(weights, energies))
    mean_E2 = sum(w * e * e for w, e in zip(weights, energies))
    free_energy = -KB * temperature * log_Z
    entropy = (mean_E - free_energy) / temperature
    heat_cap = (mean_E2 - mean_E ** 2) / (KB * temperature ** 2)
    std_energy = math.sqrt(max(0.0, mean_E2 - mean_E ** 2))

    return {
        "log_Z": log_Z,
        "free_energy": free_energy,
        "mean_energy": mean_E,
        "mean_energy_sq": mean_E2,
        "entropy": entropy,
        "heat_capacity": heat_cap,
        "std_energy": std_energy,
        "weights": weights,
    }


# ===========================================================================
# StatMechEngine – thermodynamic identity checks
# ===========================================================================

class TestStatMechEngineThermodynamics:
    T = 300.0  # Kelvin
    ENERGIES = [-10.5, -9.8, -10.2, -11.0, -9.5]  # kcal/mol
    ABS_TOL = 1e-9

    @pytest.fixture(autouse=True)
    def engine_and_reference(self):
        self.eng = _make_engine(self.T, self.ENERGIES)
        self.thermo = self.eng.compute()
        self.ref = _analytical(self.T, self.ENERGIES)

    # --- F = −kT ln Z -------------------------------------------------------

    def test_free_energy_equals_minus_kT_log_Z(self):
        expected = -KB * self.T * self.thermo.log_Z
        assert self.thermo.free_energy == pytest.approx(expected, abs=self.ABS_TOL)

    def test_log_Z_matches_analytical(self):
        assert self.thermo.log_Z == pytest.approx(self.ref["log_Z"], abs=self.ABS_TOL)

    def test_free_energy_matches_analytical(self):
        assert self.thermo.free_energy == pytest.approx(self.ref["free_energy"], abs=self.ABS_TOL)

    # --- S = (⟨E⟩ − F) / T -------------------------------------------------

    def test_entropy_satisfies_S_eq_H_minus_F_over_T(self):
        expected = (self.thermo.mean_energy - self.thermo.free_energy) / self.T
        assert self.thermo.entropy == pytest.approx(expected, abs=self.ABS_TOL)

    def test_entropy_matches_analytical(self):
        assert self.thermo.entropy == pytest.approx(self.ref["entropy"], abs=self.ABS_TOL)

    # --- Cv and σ_E ---------------------------------------------------------

    def test_heat_capacity_matches_analytical(self):
        assert self.thermo.heat_capacity == pytest.approx(self.ref["heat_capacity"], rel=1e-6)

    def test_std_energy_equals_sqrt_variance(self):
        variance = self.thermo.mean_energy_sq - self.thermo.mean_energy ** 2
        assert self.thermo.std_energy == pytest.approx(math.sqrt(variance), abs=self.ABS_TOL)

    def test_std_energy_matches_analytical(self):
        assert self.thermo.std_energy == pytest.approx(self.ref["std_energy"], abs=self.ABS_TOL)

    # --- mean ⟨E⟩ and ⟨E²⟩ -------------------------------------------------

    def test_mean_energy_matches_analytical(self):
        assert self.thermo.mean_energy == pytest.approx(self.ref["mean_energy"], abs=self.ABS_TOL)

    def test_mean_energy_sq_matches_analytical(self):
        assert self.thermo.mean_energy_sq == pytest.approx(
            self.ref["mean_energy_sq"], abs=self.ABS_TOL
        )

    # --- temperature preserved ----------------------------------------------

    def test_temperature_preserved(self):
        assert self.thermo.temperature == pytest.approx(self.T)


# ===========================================================================
# StatMechEngine – Boltzmann weights
# ===========================================================================

class TestBoltzmannWeights:
    T = 300.0
    ENERGIES = [-10.5, -9.8, -10.2, -11.0, -9.5]

    @pytest.fixture(autouse=True)
    def engine(self):
        self.eng = _make_engine(self.T, self.ENERGIES)
        self.ref = _analytical(self.T, self.ENERGIES)

    def test_weights_sum_to_one(self):
        weights = self.eng.boltzmann_weights()
        assert float(np.sum(weights)) == pytest.approx(1.0, abs=1e-12)

    def test_weights_are_non_negative(self):
        weights = self.eng.boltzmann_weights()
        assert np.all(weights >= 0.0)

    def test_weights_match_analytical(self):
        weights = self.eng.boltzmann_weights()
        assert list(weights) == pytest.approx(self.ref["weights"], abs=1e-9)

    def test_lowest_energy_has_highest_weight(self):
        """The state with the most negative energy should have the largest weight."""
        weights = self.eng.boltzmann_weights()
        min_energy_idx = int(np.argmin(self.ENERGIES))
        assert int(np.argmax(weights)) == min_energy_idx

    def test_weight_count_matches_sample_count(self):
        weights = self.eng.boltzmann_weights()
        assert len(weights) == len(self.ENERGIES)


# ===========================================================================
# StatMechEngine – single-sample edge case
# ===========================================================================

class TestSingleSample:
    """For a single state, entropy = 0, F = H = E, Cv = 0."""

    T = 300.0
    E = -12.0

    @pytest.fixture(autouse=True)
    def engine(self):
        self.eng = _make_engine(self.T, [self.E])
        self.thermo = self.eng.compute()

    def test_free_energy_equals_energy(self):
        assert self.thermo.free_energy == pytest.approx(self.E, abs=1e-9)

    def test_mean_energy_equals_energy(self):
        assert self.thermo.mean_energy == pytest.approx(self.E, abs=1e-9)

    def test_entropy_is_zero(self):
        assert self.thermo.entropy == pytest.approx(0.0, abs=1e-9)

    def test_heat_capacity_is_zero(self):
        assert self.thermo.heat_capacity == pytest.approx(0.0, abs=1e-9)

    def test_std_energy_is_zero(self):
        assert self.thermo.std_energy == pytest.approx(0.0, abs=1e-9)

    def test_weight_is_one(self):
        weights = self.eng.boltzmann_weights()
        assert float(weights[0]) == pytest.approx(1.0, abs=1e-12)


# ===========================================================================
# StatMechEngine – delta_G antisymmetry
# ===========================================================================

class TestDeltaG:
    T = 300.0

    def test_delta_G_is_difference_of_free_energies(self):
        eng_a = _make_engine(self.T, [-10.0, -9.5, -10.5])
        eng_b = _make_engine(self.T, [-8.0, -7.5])
        delta = eng_a.delta_G(eng_b)
        F_a = eng_a.compute().free_energy
        F_b = eng_b.compute().free_energy
        assert delta == pytest.approx(F_a - F_b, abs=1e-9)

    def test_delta_G_is_antisymmetric(self):
        eng_a = _make_engine(self.T, [-10.0, -9.5])
        eng_b = _make_engine(self.T, [-8.0, -7.5])
        assert eng_a.delta_G(eng_b) == pytest.approx(-eng_b.delta_G(eng_a), abs=1e-9)

    def test_delta_G_to_self_is_zero(self):
        eng = _make_engine(self.T, [-10.0, -9.5, -10.5])
        assert eng.delta_G(eng) == pytest.approx(0.0, abs=1e-9)


# ===========================================================================
# StatMechEngine – clear() and re-use
# ===========================================================================

class TestClear:
    T = 300.0

    def test_clear_resets_sample_count(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(self.T)
        eng.add_sample(-10.0)
        eng.add_sample(-9.0)
        assert eng.n_samples == 2
        eng.clear()
        assert eng.n_samples == 0

    def test_after_clear_new_samples_give_correct_result(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(self.T)
        eng.add_sample(-5.0)    # first batch – to be discarded
        eng.add_sample(-4.0)
        eng.clear()

        # Load a fresh set of known energies
        fresh = [-10.0, -9.5]
        for e in fresh:
            eng.add_sample(e)

        thermo = eng.compute()
        ref = _analytical(self.T, fresh)
        assert thermo.free_energy == pytest.approx(ref["free_energy"], abs=1e-9)
        assert thermo.entropy == pytest.approx(ref["entropy"], abs=1e-9)


# ===========================================================================
# StatMechEngine – add_samples() via NumPy array
# ===========================================================================

class TestAddSamples:
    T = 300.0
    ENERGIES = [-10.5, -9.8, -10.2, -11.0]

    def test_add_samples_matches_sequential_add_sample(self):
        """add_samples(array) must give the same result as repeated add_sample()."""
        from flexaidds.thermodynamics import StatMechEngine

        eng_bulk = StatMechEngine(self.T)
        eng_bulk.add_samples(np.array(self.ENERGIES))

        eng_seq = StatMechEngine(self.T)
        for e in self.ENERGIES:
            eng_seq.add_sample(e)

        t_bulk = eng_bulk.compute()
        t_seq = eng_seq.compute()
        assert t_bulk.free_energy == pytest.approx(t_seq.free_energy, abs=1e-9)
        assert t_bulk.entropy == pytest.approx(t_seq.entropy, abs=1e-9)
        assert t_bulk.mean_energy == pytest.approx(t_seq.mean_energy, abs=1e-9)

    def test_add_samples_sample_count(self):
        from flexaidds.thermodynamics import StatMechEngine
        eng = StatMechEngine(self.T)
        eng.add_samples(np.array(self.ENERGIES))
        assert eng.n_samples == len(self.ENERGIES)


# ===========================================================================
# helmholtz_from_energies convenience function
# ===========================================================================

class TestHelmholtzFromEnergies:
    T = 300.0
    ENERGIES = [-10.5, -9.8, -10.2, -11.0]

    def test_matches_StatMechEngine_result(self):
        from flexaidds.thermodynamics import helmholtz_from_energies
        eng = _make_engine(self.T, self.ENERGIES)
        expected_F = eng.compute().free_energy
        result_F = helmholtz_from_energies(np.array(self.ENERGIES), self.T)
        assert result_F == pytest.approx(expected_F, abs=1e-9)

    def test_matches_analytical(self):
        from flexaidds.thermodynamics import helmholtz_from_energies
        ref = _analytical(self.T, self.ENERGIES)
        result = helmholtz_from_energies(np.array(self.ENERGIES), self.T)
        assert result == pytest.approx(ref["free_energy"], abs=1e-9)

    def test_single_energy_returns_that_energy(self):
        from flexaidds.thermodynamics import helmholtz_from_energies
        E = -7.5
        result = helmholtz_from_energies(np.array([E]), self.T)
        assert result == pytest.approx(E, abs=1e-9)


