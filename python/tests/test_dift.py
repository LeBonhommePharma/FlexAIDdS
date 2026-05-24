"""Tests for the flexaidds.dift module (pure-Python DiFT engine).

These tests verify the DiFT torsional parametrization engine — forward
transform, Shannon-collapse truncation, QM–MM refinement, per-bond
thermodynamics, Boltzmann inversion, and the GA scoring adapter — without
needing the compiled C++ extension.

SPDX-License-Identifier: Apache-2.0
"""

import math

import numpy as np
import pytest

from flexaidds.dift import (
    DiFTEngine,
    FourierTerm,
    TorsionalPotential,
    RotatableBondTorsion,
    score_torsional,
    make_bond_torsion,
    spectral_entropy,
    kB_kcal,
)

_TWO_PI = 2.0 * math.pi


# ── helpers ──────────────────────────────────────────────────────────────────

def make_profile(m, mean, mult, amp, phase):
    """M-point sample of V(φ) = mean + Σ Aₙ cos(nφ − ωₙ)."""
    k = np.arange(m, dtype=float)
    phi = k * (_TWO_PI / m)
    v = np.full(m, float(mean))
    for n, a, w in zip(mult, amp, phase):
        v += a * np.cos(n * phi - w)
    return v


def amp_of(spectrum, n):
    for t in spectrum:
        if t.multiplicity == n:
            return t.amplitude
    return 0.0


# ── forward transform ────────────────────────────────────────────────────────

def test_transform_recovers_amplitudes_pow2():
    engine = DiFTEngine(300.0)
    profile = make_profile(64, 1.5, [1, 3, 5], [2.0, 1.0, 0.5], [0.0, 0.7, -1.2])
    spectrum, mean = engine.transform(profile)
    assert mean == pytest.approx(1.5, abs=1e-9)
    assert amp_of(spectrum, 1) == pytest.approx(2.0, abs=1e-9)
    assert amp_of(spectrum, 3) == pytest.approx(1.0, abs=1e-9)
    assert amp_of(spectrum, 5) == pytest.approx(0.5, abs=1e-9)
    assert amp_of(spectrum, 2) == pytest.approx(0.0, abs=1e-9)


def test_transform_recovers_amplitudes_non_pow2():
    engine = DiFTEngine(300.0)
    profile = make_profile(36, -0.4, [2, 4], [1.3, 0.8], [0.3, 2.1])
    spectrum, mean = engine.transform(profile)
    assert mean == pytest.approx(-0.4, abs=1e-9)
    assert amp_of(spectrum, 2) == pytest.approx(1.3, abs=1e-8)
    assert amp_of(spectrum, 4) == pytest.approx(0.8, abs=1e-8)


def test_transform_round_trip():
    engine = DiFTEngine(300.0)
    m = 72
    profile = make_profile(m, 0.9, [1, 2, 6], [1.1, 0.6, 0.25], [-0.5, 1.0, 2.8])
    spectrum, mean = engine.transform(profile)
    pot = TorsionalPotential(terms=spectrum, mean=mean)
    model = pot.sample(m)
    np.testing.assert_allclose(model, profile, atol=1e-7)
    assert DiFTEngine.r_squared(profile, model) == pytest.approx(1.0, abs=1e-9)


def test_transform_rejects_short_profile():
    engine = DiFTEngine()
    with pytest.raises(ValueError):
        engine.transform([1.0])


# ── spectral Shannon entropy ─────────────────────────────────────────────────

def test_spectral_entropy_single_mode_collapses():
    engine = DiFTEngine()
    spectrum, _ = engine.transform(make_profile(64, 0.0, [2], [3.0], [0.4]))
    h = spectral_entropy(spectrum)
    assert h == pytest.approx(0.0, abs=1e-9)
    assert math.exp(h) == pytest.approx(1.0, abs=1e-9)


def test_spectral_entropy_equal_power_spreads():
    engine = DiFTEngine()
    spectrum, _ = engine.transform(
        make_profile(64, 0.0, [1, 2, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]))
    assert spectral_entropy(spectrum) == pytest.approx(math.log(3.0), abs=1e-9)


# ── parametrize: Shannon-collapse truncation ─────────────────────────────────

def test_parametrize_keeps_dominant_modes():
    engine = DiFTEngine()
    profile = make_profile(128, 0.2, [1, 4, 11],
                           [2.5, 1.8, 0.02], [0.1, -0.6, 1.4])
    pot = engine.parametrize(profile)
    assert 2 <= pot.n_terms <= 3
    assert 1.0 < pot.effective_modes < 3.0
    assert pot.r_squared > 0.99


def test_parametrize_max_multiplicity_guard():
    engine = DiFTEngine()
    profile = make_profile(64, 0.0, [2, 9], [1.0, 1.0], [0.0, 0.0])
    pot = engine.parametrize(profile, max_multiplicity=6)
    assert all(t.multiplicity <= 6 for t in pot.terms)


# ── iterative QM–MM refinement ───────────────────────────────────────────────

def test_refine_converges_to_qm():
    engine = DiFTEngine(300.0)
    m = 96
    qm = make_profile(m, 0.0, [1, 2, 3], [3.0, 1.5, 0.8], [0.2, -0.9, 1.7])
    mm = make_profile(m, 0.0, [1, 2], [2.0, 1.0], [0.2, -0.9])
    corrected = engine.refine(qm, mm, lambda_=0.5, r2_target=0.98)
    assert corrected.r_squared >= 0.98
    assert corrected.refinement_iters > 0
    total = mm + corrected.sample(m)
    assert DiFTEngine.r_squared(qm, total) >= 0.98


def test_refine_rejects_mismatched_grids():
    engine = DiFTEngine()
    with pytest.raises(ValueError):
        engine.refine([0.0] * 64, [0.0] * 32)


# ── per-bond torsional thermodynamics ────────────────────────────────────────

def test_thermo_free_rotor_zero_entropy():
    engine = DiFTEngine(300.0)
    pot = engine.parametrize([0.0] * 64)
    th = engine.thermodynamics(pot)
    assert th.partition_function == pytest.approx(1.0, abs=1e-6)
    assert th.entropy == pytest.approx(0.0, abs=1e-6)
    assert th.minus_TS == pytest.approx(0.0, abs=1e-6)


def test_thermo_confined_bond_pays_penalty():
    engine = DiFTEngine(300.0)
    pot = engine.parametrize(make_profile(128, 0.0, [1], [5.0], [0.0]))
    th = engine.thermodynamics(pot)
    assert th.partition_function < 1.0
    assert th.entropy < 0.0          # entropy loss vs free rotor
    assert th.minus_TS > 0.0         # → positive ΔG penalty
    assert th.mean_energy > 0.0


def test_thermo_deeper_well_pays_more():
    engine = DiFTEngine(300.0)
    shallow = engine.parametrize(make_profile(128, 0.0, [1], [1.0], [0.0]))
    deep = engine.parametrize(make_profile(128, 0.0, [1], [6.0], [0.0]))
    assert (engine.thermodynamics(deep).minus_TS
            > engine.thermodynamics(shallow).minus_TS)


# ── Boltzmann inversion ──────────────────────────────────────────────────────

def test_boltzmann_invert_uniform_is_flat():
    engine = DiFTEngine(300.0)
    energy = engine.boltzmann_invert([100.0] * 36)
    np.testing.assert_allclose(energy, 0.0, atol=1e-9)


def test_boltzmann_invert_peaked_becomes_well():
    engine = DiFTEngine(300.0)
    hist = [1.0] * 36
    hist[18] = 1000.0
    energy = engine.boltzmann_invert(hist)
    assert energy[18] == pytest.approx(0.0, abs=1e-9)
    assert np.all(energy >= -1e-12)
    assert np.all(np.delete(energy, 18) > 0.0)


def test_boltzmann_invert_rejects_empty():
    engine = DiFTEngine()
    with pytest.raises(ValueError):
        engine.boltzmann_invert([0.0] * 20)


# ── circular mean ────────────────────────────────────────────────────────────

def test_circular_mean_wraparound():
    angles = [math.radians(170.0), math.radians(-170.0)]
    assert abs(DiFTEngine.circular_mean(angles)) == pytest.approx(math.pi, abs=1e-6)


def test_circular_mean_simple():
    assert DiFTEngine.circular_mean([0.1, 0.2, 0.3]) == pytest.approx(0.2, abs=1e-6)


# ── GA scoring adapter ───────────────────────────────────────────────────────

def test_adapter_scores_energy_and_entropy():
    rbt0 = make_bond_torsion(make_profile(128, 0.0, [1], [3.0], [math.pi]), 0)
    rbt1 = make_bond_torsion(make_profile(128, 0.0, [2], [2.0], [math.pi]), 1)
    score = score_torsional([rbt0, rbt1], [0.0, 0.0], 300.0)
    assert score.n_bonds == 2
    assert score.energy == pytest.approx(0.0, abs=1e-3)
    assert score.minus_TS > 0.0
    assert score.total() > 0.0


def test_adapter_energy_rises_away_from_minimum():
    rbt = make_bond_torsion(make_profile(128, 0.0, [1], [4.0], [math.pi]), 0)
    e_min = score_torsional([rbt], [0.0]).energy
    e_top = score_torsional([rbt], [math.pi]).energy
    assert e_top > e_min + 1.0


def test_adapter_mismatched_input_returns_zero():
    rbt = make_bond_torsion([0.0] * 64, 0)
    score = score_torsional([rbt], [0.0, 0.0])
    assert score.n_bonds == 0
    assert score.total() == 0.0


# ── cross-check: C++ ⇔ Python parity (skipped if _core lacks DiFT) ───────────

def test_cpp_python_parity_if_available():
    """If the compiled core exposes DiFTEngine, results must match Python."""
    try:
        from flexaidds import _core
        cpp_engine = _core.DiFTEngine(300.0)
    except (ImportError, AttributeError):
        pytest.skip("_core.DiFTEngine not built")

    profile = make_profile(64, 0.3, [1, 3], [2.0, 1.0], [0.4, -0.8])
    py_pot = DiFTEngine(300.0).parametrize(list(profile))
    cpp_pot = cpp_engine.parametrize(list(profile))
    assert cpp_pot.r_squared == pytest.approx(py_pot.r_squared, abs=1e-6)
    assert cpp_pot.effective_modes == pytest.approx(py_pot.effective_modes, abs=1e-6)
