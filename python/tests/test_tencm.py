"""Tests for the flexaidds.tencm module (pure-Python fallback path).

These tests verify the TorsionalENM, ShannonThermoStack, and related
functions work correctly without the C++ extension.
"""

import math
import os
import tempfile

import pytest
import numpy as np

from flexaidds.tencm import (
    TorsionalENM,
    TorsionalNormalMode,
    Conformer,
    FullThermoResult,
    compute_boltzmann_shannon_entropy,
    compute_shannon_entropy,
    compute_torsional_vibrational_entropy,
    run_shannon_thermo_stack,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _write_synthetic_pdb(n_residues: int, perturb: float = 0.0) -> str:
    """Write a synthetic alpha-helix PDB and return the path."""
    fd, path = tempfile.mkstemp(suffix=".pdb")
    radius, rise = 2.3, 1.5
    turn = math.radians(100)
    with os.fdopen(fd, "w") as f:
        for r in range(n_residues):
            x = radius * math.cos(r * turn) + perturb * (0.5 if r % 3 == 0 else 0.0)
            y = radius * math.sin(r * turn) + perturb * (0.3 if r % 3 == 1 else 0.0)
            z = r * rise + perturb * (0.2 if r % 3 == 2 else 0.0)
            f.write(
                f"ATOM  {r+1:5d}  CA  ALA A{r+1:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        f.write("END\n")
    return path


# ── TorsionalENM Tests ───────────────────────────────────────────────────────

class TestTorsionalENM:
    def test_build_from_pdb(self):
        path = _write_synthetic_pdb(20)
        try:
            tenm = TorsionalENM()
            tenm.build_from_pdb(path)
            assert tenm.is_built
            assert tenm.n_residues == 20
            assert tenm.n_bonds == 19
            assert tenm.n_modes > 0
        finally:
            os.unlink(path)

    def test_not_built_initially(self):
        tenm = TorsionalENM()
        assert not tenm.is_built
        assert tenm.n_modes == 0

    def test_modes_have_positive_eigenvalues(self):
        path = _write_synthetic_pdb(30)
        try:
            tenm = TorsionalENM()
            tenm.build_from_pdb(path)
            for mode in tenm.modes:
                assert mode.eigenvalue > 0
        finally:
            os.unlink(path)

    def test_too_few_residues(self):
        path = _write_synthetic_pdb(3)
        try:
            tenm = TorsionalENM()
            with pytest.raises(ValueError):
                tenm.build_from_pdb(path)
        finally:
            os.unlink(path)

    def test_repr(self):
        tenm = TorsionalENM()
        assert "not built" in repr(tenm)


# ── Shannon Entropy Tests ────────────────────────────────────────────────────

class TestShannonEntropy:
    def test_uniform_distribution(self):
        # Uniform distribution over N bins → ln(N) nats
        values = list(range(100))
        H = compute_shannon_entropy(values, num_bins=10)
        # Should be close to ln(10) ≈ 2.30 for uniform
        assert H > 2.0
        assert H <= math.log(10) + 0.1

    def test_constant_distribution(self):
        # All same value → 0 entropy
        values = [5.0] * 100
        H = compute_shannon_entropy(values, num_bins=10)
        # All in one bin → H = 0
        assert H == 0.0

    def test_empty_returns_zero(self):
        assert compute_shannon_entropy([]) == 0.0

    def test_two_bins(self):
        # 50/50 split → ln(2) ≈ 0.693 nats
        values = [0.0] * 50 + [1.0] * 50
        H = compute_shannon_entropy(values, num_bins=2)
        assert abs(H - math.log(2)) < 0.01

    def test_boltzmann_entropy_equal_energies(self):
        H = compute_boltzmann_shannon_entropy([-10.0] * 5)
        assert H == pytest.approx(math.log(5.0))

    def test_boltzmann_entropy_dominant_basin(self):
        H = compute_boltzmann_shannon_entropy([-30.0, 0.0, 0.0, 0.0, 0.0])
        assert H < math.log(2.0)


# ── Torsional Vibrational Entropy Tests ──────────────────────────────────────

class TestTorsionalVibrationalEntropy:
    def test_positive_entropy(self):
        modes = [TorsionalNormalMode(eigenvalue=float(i + 1)) for i in range(10)]
        S = compute_torsional_vibrational_entropy(modes, temperature_K=300.0)
        assert S != 0.0
        assert math.isfinite(S)

    def test_empty_modes_zero(self):
        S = compute_torsional_vibrational_entropy([], temperature_K=300.0)
        assert S == 0.0

    def test_zero_eigenvalue_skipped(self):
        modes = [TorsionalNormalMode(eigenvalue=0.0)]
        S = compute_torsional_vibrational_entropy(modes, temperature_K=300.0)
        assert S == 0.0


# ── ShannonThermoStack Tests ─────────────────────────────────────────────────

class TestShannonThermoStack:
    def test_basic_result(self):
        energies = [-10.0, -12.0, -8.0, -11.0, -9.0]
        result = run_shannon_thermo_stack(energies, base_deltaG=-10.0)
        assert isinstance(result, FullThermoResult)
        assert math.isfinite(result.deltaG)
        assert result.shannonEntropy >= 0.0

    def test_with_tencm_model(self):
        path = _write_synthetic_pdb(20)
        try:
            tenm = TorsionalENM()
            tenm.build_from_pdb(path)

            energies = [-10.0, -12.0, -8.0]
            result = run_shannon_thermo_stack(
                energies, tencm_model=tenm, base_deltaG=-10.0)

            # Heuristic S_vib is reported but not folded into kcal/mol dG.
            assert result.torsionalVibEntropy != 0.0
            assert result.entropyContribution != 0.0
            assert math.isfinite(result.deltaG)
            expected_contrib = -298.15 * result.shannonEntropy * 0.001987206
            assert result.entropyContribution == pytest.approx(expected_contrib)
            assert result.deltaG == pytest.approx(-10.0 + expected_contrib)
            assert "excluded from dG" in result.report
        finally:
            os.unlink(path)

    def test_no_tencm_zero_vibrational(self):
        energies = [-10.0, -12.0]
        result = run_shannon_thermo_stack(energies, base_deltaG=-5.0)
        assert result.torsionalVibEntropy == 0.0
        # Shannon conf entropy contributes via S_conf = k_B * H_nats
        H = result.shannonEntropy
        kB = 0.001987206
        expected_S = H * kB
        expected_contrib = -298.15 * expected_S
        assert abs(result.entropyContribution - expected_contrib) < 1e-8
        assert abs(result.deltaG - (-5.0 + expected_contrib)) < 1e-8

    def test_equal_energies_are_not_entropy_collapse(self):
        result = run_shannon_thermo_stack([-10.0] * 5, base_deltaG=-5.0)
        assert result.shannonEntropy == pytest.approx(math.log(5.0))

    def test_report_string(self):
        energies = [-10.0, -12.0, -8.0]
        result = run_shannon_thermo_stack(energies, base_deltaG=-10.0)
        assert "ShannonThermoStack" in result.report
        assert "Shannon" in result.report
        assert "nats" in result.report
        assert "bits" not in result.report

    def test_repr(self):
        result = FullThermoResult(deltaG=-10.0, shannonEntropy=2.5)
        assert "ΔG=" in repr(result)


# ── TorsionalNormalMode Tests ────────────────────────────────────────────────

class TestTorsionalNormalMode:
    def test_repr(self):
        m = TorsionalNormalMode(eigenvalue=1.5)
        assert "1.5" in repr(m)

    def test_default_values(self):
        m = TorsionalNormalMode()
        assert m.eigenvalue == 0.0
        assert m.eigenvector == []
