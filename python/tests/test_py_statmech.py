"""Tests for _PyStatMechEngine – the pure-Python canonical ensemble fallback.

These tests verify the mathematical correctness of the pure-Python
implementation that runs when the C++ extension is not available.  They do not
require the compiled _core module and run in every environment.

Physical definitions used throughout
--------------------------------------
  β = 1 / (kB T)            kB = 0.001987206 kcal mol⁻¹ K⁻¹
  Z = Σ exp(−β Eᵢ)          partition function
  F = −kT ln Z               Helmholtz free energy
  ⟨E⟩ = Σ wᵢ Eᵢ             Boltzmann-weighted mean
  S = (⟨E⟩ − F) / T         configurational entropy
  Cv = (⟨E²⟩ − ⟨E⟩²)/(kT²) heat capacity at constant volume
  σ_E = sqrt(⟨E²⟩ − ⟨E⟩²)  standard deviation of energy
"""

from __future__ import annotations

import math

import pytest

from flexaidds.thermodynamics import _PyStatMechEngine, kB_kcal

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

KB = kB_kcal  # kcal mol⁻¹ K⁻¹


def _make(temperature: float, energies) -> _PyStatMechEngine:
    """Return a _PyStatMechEngine loaded with the given energies."""
    eng = _PyStatMechEngine(temperature)
    for e in energies:
        eng.add_sample(float(e))
    return eng


def _analytical(temperature: float, energies: list) -> dict:
    """Compute ground-truth thermodynamics from first principles."""
    beta = 1.0 / (KB * temperature)
    x = [-beta * e for e in energies]
    max_x = max(x)
    log_Z = max_x + math.log(sum(math.exp(xi - max_x) for xi in x))

    weights = [math.exp(xi - log_Z) for xi in x]
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
# _PyStatMechEngine – thermodynamic identity checks
# ===========================================================================

class TestPyStatMechEngineThermodynamics:
    T = 300.0
    ENERGIES = [-10.5, -9.8, -10.2, -11.0, -9.5]
    ABS_TOL = 1e-9

    @pytest.fixture(autouse=True)
    def engine_and_reference(self):
        self.eng = _make(self.T, self.ENERGIES)
        self.thermo = self.eng.compute()
        self.ref = _analytical(self.T, self.ENERGIES)

    def test_free_energy_equals_minus_kT_log_Z(self):
        expected = -KB * self.T * self.thermo.log_Z
        assert self.thermo.free_energy == pytest.approx(expected, abs=self.ABS_TOL)

    def test_log_Z_matches_analytical(self):
        assert self.thermo.log_Z == pytest.approx(self.ref["log_Z"], abs=self.ABS_TOL)

    def test_free_energy_matches_analytical(self):
        assert self.thermo.free_energy == pytest.approx(self.ref["free_energy"], abs=self.ABS_TOL)

    def test_entropy_satisfies_S_eq_H_minus_F_over_T(self):
        expected = (self.thermo.mean_energy - self.thermo.free_energy) / self.T
        assert self.thermo.entropy == pytest.approx(expected, abs=self.ABS_TOL)

    def test_entropy_matches_analytical(self):
        assert self.thermo.entropy == pytest.approx(self.ref["entropy"], abs=self.ABS_TOL)

    def test_heat_capacity_matches_analytical(self):
        assert self.thermo.heat_capacity == pytest.approx(self.ref["heat_capacity"], rel=1e-6)

    def test_std_energy_equals_sqrt_variance(self):
        variance = self.thermo.mean_energy_sq - self.thermo.mean_energy ** 2
        assert self.thermo.std_energy == pytest.approx(math.sqrt(variance), abs=self.ABS_TOL)

    def test_std_energy_matches_analytical(self):
        assert self.thermo.std_energy == pytest.approx(self.ref["std_energy"], abs=self.ABS_TOL)

    def test_mean_energy_matches_analytical(self):
        assert self.thermo.mean_energy == pytest.approx(self.ref["mean_energy"], abs=self.ABS_TOL)

    def test_mean_energy_sq_matches_analytical(self):
        assert self.thermo.mean_energy_sq == pytest.approx(
            self.ref["mean_energy_sq"], abs=self.ABS_TOL
        )

    def test_temperature_preserved(self):
        assert self.thermo.temperature == pytest.approx(self.T)

    def test_free_energy_leq_mean_energy(self):
        """F ≤ ⟨E⟩ always holds (entropy term is non-negative)."""
        assert self.thermo.free_energy <= self.thermo.mean_energy + 1e-12

    def test_entropy_non_negative(self):
        assert self.thermo.entropy >= -1e-12


# ===========================================================================
# _PyStatMechEngine – Boltzmann weights
# ===========================================================================

class TestPyBoltzmannWeights:
    T = 300.0
    ENERGIES = [-10.5, -9.8, -10.2, -11.0, -9.5]

    @pytest.fixture(autouse=True)
    def engine(self):
        self.eng = _make(self.T, self.ENERGIES)
        self.ref = _analytical(self.T, self.ENERGIES)

    def test_weights_sum_to_one(self):
        weights = self.eng.boltzmann_weights()
        assert abs(sum(weights) - 1.0) < 1e-12

    def test_weights_are_non_negative(self):
        weights = self.eng.boltzmann_weights()
        assert all(w >= 0.0 for w in weights)

    def test_weights_match_analytical(self):
        weights = self.eng.boltzmann_weights()
        assert weights == pytest.approx(self.ref["weights"], abs=1e-9)

    def test_lowest_energy_has_highest_weight(self):
        """The most negative energy must correspond to the largest Boltzmann weight."""
        weights = self.eng.boltzmann_weights()
        min_idx = self.ENERGIES.index(min(self.ENERGIES))
        assert weights.index(max(weights)) == min_idx

    def test_weight_count_matches_sample_count(self):
        weights = self.eng.boltzmann_weights()
        assert len(weights) == len(self.ENERGIES)

    def test_empty_engine_returns_empty_list(self):
        eng = _PyStatMechEngine(self.T)
        assert eng.boltzmann_weights() == []


# ===========================================================================
# _PyStatMechEngine – single-sample edge case
# ===========================================================================

class TestPySingleSample:
    """For a single state, entropy = 0, F = ⟨E⟩ = E, Cv = 0."""

    T = 300.0
    E = -12.0

    @pytest.fixture(autouse=True)
    def engine(self):
        self.eng = _make(self.T, [self.E])
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
        assert weights[0] == pytest.approx(1.0, abs=1e-12)


# ===========================================================================
# _PyStatMechEngine – delta_G
# ===========================================================================

class TestPyDeltaG:
    T = 300.0

    def test_delta_G_is_difference_of_free_energies(self):
        eng_a = _make(self.T, [-10.0, -9.5, -10.5])
        eng_b = _make(self.T, [-8.0, -7.5])
        delta = eng_a.delta_G(eng_b)
        F_a = eng_a.compute().free_energy
        F_b = eng_b.compute().free_energy
        assert delta == pytest.approx(F_a - F_b, abs=1e-9)

    def test_delta_G_is_antisymmetric(self):
        eng_a = _make(self.T, [-10.0, -9.5])
        eng_b = _make(self.T, [-8.0, -7.5])
        assert eng_a.delta_G(eng_b) == pytest.approx(-eng_b.delta_G(eng_a), abs=1e-9)

    def test_delta_G_to_self_is_zero(self):
        eng = _make(self.T, [-10.0, -9.5, -10.5])
        assert eng.delta_G(eng) == pytest.approx(0.0, abs=1e-9)

    def test_lower_energy_ensemble_has_negative_delta_G(self):
        """The ensemble with lower energies has more negative free energy."""
        eng_low = _make(self.T, [-12.0, -11.5])
        eng_high = _make(self.T, [-8.0, -7.5])
        assert eng_low.delta_G(eng_high) < 0.0


# ===========================================================================
# _PyStatMechEngine – clear() and re-use
# ===========================================================================

class TestPyClear:
    T = 300.0

    def test_clear_resets_sample_count(self):
        eng = _make(self.T, [-10.0, -9.0])
        assert eng.size == 2
        eng.clear()
        assert eng.size == 0

    def test_clear_raises_on_compute(self):
        eng = _make(self.T, [-10.0])
        eng.clear()
        with pytest.raises(RuntimeError, match="No samples"):
            eng.compute()

    def test_after_clear_new_samples_give_correct_result(self):
        eng = _make(self.T, [-5.0, -4.0])  # first batch – to be discarded
        eng.clear()
        fresh = [-10.0, -9.5]
        for e in fresh:
            eng.add_sample(e)
        thermo = eng.compute()
        ref = _analytical(self.T, fresh)
        assert thermo.free_energy == pytest.approx(ref["free_energy"], abs=1e-9)
        assert thermo.entropy == pytest.approx(ref["entropy"], abs=1e-9)


# ===========================================================================
# _PyStatMechEngine – multiplicity support
# ===========================================================================

class TestPyMultiplicity:
    T = 300.0

    def test_multiplicity_two_keeps_one_microstate(self):
        eng = _PyStatMechEngine(self.T)
        eng.add_sample(-10.0, multiplicity=2)
        assert eng.size == 1

    def test_multiplicity_equivalent_to_repeated_add(self):
        """add_sample(E, 3) must produce the same result as three add_sample(E, 1) calls."""
        eng_mult = _PyStatMechEngine(self.T)
        eng_mult.add_sample(-10.0, multiplicity=3)
        eng_mult.add_sample(-9.0, multiplicity=1)

        eng_seq = _PyStatMechEngine(self.T)
        for _ in range(3):
            eng_seq.add_sample(-10.0)
        eng_seq.add_sample(-9.0)

        t_m = eng_mult.compute()
        t_s = eng_seq.compute()
        assert t_m.free_energy == pytest.approx(t_s.free_energy, abs=1e-9)
        assert t_m.entropy == pytest.approx(t_s.entropy, abs=1e-9)
        assert t_m.mean_energy == pytest.approx(t_s.mean_energy, abs=1e-9)


# ===========================================================================
# _PyStatMechEngine – properties (temperature, beta, size)
# ===========================================================================

def test_py_temperature_property():
    eng = _PyStatMechEngine(310.0)
    assert eng.temperature == pytest.approx(310.0)


def test_py_beta_property():
    T = 300.0
    eng = _PyStatMechEngine(T)
    assert eng.beta == pytest.approx(1.0 / (KB * T), rel=1e-9)


def test_py_size_increments_with_each_add_sample():
    eng = _PyStatMechEngine(300.0)
    assert eng.size == 0
    eng.add_sample(-10.0)
    assert eng.size == 1
    eng.add_sample(-9.0)
    assert eng.size == 2


def test_py_compute_raises_when_empty():
    eng = _PyStatMechEngine(300.0)
    with pytest.raises(RuntimeError, match="No samples"):
        eng.compute()


# ===========================================================================
# _PyStatMechEngine – temperature sensitivity
# ===========================================================================

def test_py_higher_temperature_increases_entropy():
    """Configurational entropy grows with temperature for an identical energy set."""
    energies = [-10.5, -9.8, -10.2, -11.0, -9.5]
    eng_low = _make(200.0, energies)
    eng_high = _make(400.0, energies)
    assert eng_high.compute().entropy > eng_low.compute().entropy


def test_py_higher_temperature_increases_free_energy():
    """At higher T, F = ⟨E⟩ − TS becomes more negative (larger −TS term)."""
    energies = [-10.0, -9.0]
    eng_low = _make(200.0, energies)
    eng_high = _make(400.0, energies)
    # TS grows with T, so F = H - TS becomes more negative
    assert eng_high.compute().free_energy < eng_low.compute().free_energy
