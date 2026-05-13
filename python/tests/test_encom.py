"""Tests for ENCoM vibrational entropy — Python-level checks + C++ integration.

Pure-Python path tests run without the C++ extension.
Tests that specifically validate C++ agreement are marked with needs_core.
"""

from __future__ import annotations

import math
from pathlib import Path

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


# ─────────────────────────────────────────────────────────────────────────────
# NormalMode dataclass (pure Python)
# ─────────────────────────────────────────────────────────────────────────────

def test_normalmode_repr():
    from flexaidds.encom import NormalMode

    nm = NormalMode(index=7, eigenvalue=3.14)
    assert "NormalMode" in repr(nm)
    assert "7" in repr(nm)


def test_normalmode_defaults():
    from flexaidds.encom import NormalMode

    nm = NormalMode()
    assert nm.index == 0
    assert nm.eigenvalue == 0.0
    assert nm.eigenvector == []


# ─────────────────────────────────────────────────────────────────────────────
# VibrationalEntropy dataclass (pure Python)
# ─────────────────────────────────────────────────────────────────────────────

def test_vibrational_entropy_free_energy_correction():
    from flexaidds.encom import VibrationalEntropy

    vs = VibrationalEntropy(S_vib_kcal_mol_K=0.002, temperature=300.0)
    expected = -300.0 * 0.002  # = -0.6 kcal/mol
    assert math.isclose(vs.free_energy_correction, expected, abs_tol=1e-12)


def test_vibrational_entropy_repr():
    from flexaidds.encom import VibrationalEntropy

    vs = VibrationalEntropy(S_vib_kcal_mol_K=0.001, n_modes=12, temperature=298.0)
    r = repr(vs)
    assert "VibrationalEntropy" in r
    assert "12" in r


def test_frequency_calibration_metadata():
    from flexaidds.encom import FrequencyCalibration

    cal = FrequencyCalibration.calibrated_scale(
        1.0e12,
        "synthetic-test",
        "unit-test scale",
    )
    assert cal.calibrated is True
    assert cal.status == "calibrated"
    assert cal.eigenvalue_to_omega == pytest.approx(1.0e12)


# ─────────────────────────────────────────────────────────────────────────────
# ENCoMEngine – pure-Python path (no _core)
# ─────────────────────────────────────────────────────────────────────────────

def test_encom_total_entropy_additive():
    from flexaidds.encom import ENCoMEngine

    S_conf = 0.0030
    S_vib  = 0.0015
    total  = ENCoMEngine.total_entropy(S_conf, S_vib)
    assert math.isclose(total, S_conf + S_vib, abs_tol=1e-15)


def test_encom_free_energy_with_vibrations():
    from flexaidds.encom import ENCoMEngine

    F_elec = -10.0
    S_vib  = 0.002
    T      = 300.0
    F_total = ENCoMEngine.free_energy_with_vibrations(F_elec, S_vib, T)
    expected = F_elec - T * S_vib  # = -10.6
    assert math.isclose(F_total, expected, abs_tol=1e-10)


def test_encom_compute_vibrational_entropy_all_trivial():
    """All-zero eigenvalues (translations/rotations only) -> S_vib = 0."""
    from flexaidds.encom import ENCoMEngine, NormalMode

    modes = [NormalMode(index=i, eigenvalue=0.0, frequency=0.0) for i in range(6)]
    vs = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=300.0)
    assert vs.n_modes == 0
    assert vs.S_vib_kcal_mol_K == 0.0


def test_encom_compute_vibrational_entropy_single_mode():
    """One non-trivial mode -> n_modes == 1, S_vib > 0 at 300 K."""
    from flexaidds.encom import ENCoMEngine, NormalMode

    modes = [
        NormalMode(index=1, eigenvalue=0.0),
        NormalMode(index=2, eigenvalue=0.0),
        NormalMode(index=3, eigenvalue=4.5),
    ]
    vs = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=300.0)
    assert vs.n_modes == 1
    assert vs.S_vib_kcal_mol_K > 0.0


def test_encom_compute_vibrational_entropy_higher_T_more_entropy():
    """S_vib increases with temperature (more accessible states)."""
    from flexaidds.encom import ENCoMEngine, NormalMode

    modes = [NormalMode(index=i + 1, eigenvalue=float(i + 1)) for i in range(6)]

    vs_low  = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=200.0)
    vs_high = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=400.0)
    assert vs_high.S_vib_kcal_mol_K > vs_low.S_vib_kcal_mol_K


def test_encom_compute_vibrational_entropy_reports_model_scale_by_default():
    from flexaidds.encom import ENCoMEngine, NormalMode

    modes = [NormalMode(index=i + 1, eigenvalue=float(i + 1)) for i in range(6)]
    vs = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=300.0)
    assert vs.calibrated is False
    assert vs.absolute_claim_allowed is False
    assert vs.calibration_status == "model_scale_heuristic"


def test_encom_compute_vibrational_entropy_uses_calibrated_scale():
    from flexaidds.encom import ENCoMEngine, FrequencyCalibration, NormalMode

    modes = [NormalMode(index=i + 1, eigenvalue=float(i + 1)) for i in range(6)]
    base = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=300.0)
    cal = FrequencyCalibration.calibrated_scale(1.0e12, "synthetic-test")
    scaled = ENCoMEngine.compute_vibrational_entropy(
        modes,
        temperature_K=300.0,
        frequency_calibration=cal,
    )
    assert scaled.calibrated is True
    assert scaled.absolute_claim_allowed is True
    assert scaled.calibration_label == "synthetic-test"
    assert scaled.omega_eff == pytest.approx(base.omega_eff * 1.0e12)


def test_encom_load_modes_pure_python(encom_files):
    """load_modes parses eigenvalue/eigenvector files into NormalMode list."""
    from flexaidds.encom import ENCoMEngine

    ev_file, evc_file = encom_files
    modes = ENCoMEngine.load_modes(str(ev_file), str(evc_file))

    # 9 eigenvalues -> 9 modes; first 6 are zero (trivial)
    assert len(modes) == 9
    assert all(m.eigenvalue == 0.0 for m in modes[:6])
    assert modes[6].eigenvalue == pytest.approx(1.23, abs=1e-9)


def test_encom_load_modes_frequencies_match_eigenvalues(encom_files):
    """Frequency = sqrt(eigenvalue) for each mode."""
    from flexaidds.encom import ENCoMEngine

    ev_file, evc_file = encom_files
    modes = ENCoMEngine.load_modes(str(ev_file), str(evc_file))

    for m in modes:
        expected_freq = math.sqrt(max(m.eigenvalue, 0.0))
        assert math.isclose(m.frequency, expected_freq, abs_tol=1e-12)


# ── C++ ENCoM integration tests ──────────────────────────────────────────────

@needs_core
class TestNormalModeCpp:
    def _make(self, eigenvalue=1.0, frequency=100.0):
        from flexaidds import NormalMode
        return NormalMode(eigenvalue=eigenvalue, frequency=frequency)

    def test_eigenvalue_stored(self):
        nm = self._make(eigenvalue=2.5)
        assert abs(nm.eigenvalue - 2.5) < 1e-12

    def test_frequency_stored(self):
        nm = self._make(frequency=150.0)
        assert abs(nm.frequency - 150.0) < 1e-9

    def test_zero_eigenvalue_allowed(self):
        nm = self._make(eigenvalue=0.0, frequency=0.0)
        assert nm.eigenvalue == 0.0


@needs_core
class TestVibrationalEntropyCpp:
    def _make(self, s_vib=0.01, dG_vib=-3.0, temperature=300.0):
        import flexaidds._core as _core
        return _core.VibrationalEntropy(
            S_vib_kcal_mol_K=s_vib,
            dG_vib_kcal_mol=dG_vib,
            temperature=temperature,
        )

    def test_s_vib_stored(self):
        ve = self._make(s_vib=0.02)
        assert abs(ve.S_vib_kcal_mol_K - 0.02) < 1e-12

    def test_dG_vib_computed(self):
        """dG_vib = -T * S_vib (computed property)."""
        ve = self._make(s_vib=0.01, temperature=300.0)
        assert abs(ve.dG_vib_kcal_mol - (-300.0 * 0.01)) < 1e-12

    def test_temperature_stored(self):
        ve = self._make(temperature=310.0)
        assert abs(ve.temperature - 310.0) < 1e-9

    def test_ts_vib_equals_t_times_s(self):
        """−TΔS_vib should equal temperature x S_vib."""
        ve = self._make(s_vib=0.01, temperature=300.0)
        assert abs(ve.TS_vib_kcal_mol - 300.0 * 0.01) < 1e-9


@needs_core
class TestENCoMEngineCpp:
    def test_default_cutoff(self):
        import flexaidds._core as _core
        eng = _core.ENCoMEngine()
        # Default eigenvalue_cutoff should be small positive
        assert eng.eigenvalue_cutoff > 0.0

    def test_custom_cutoff(self):
        import flexaidds._core as _core
        eng = _core.ENCoMEngine(eigenvalue_cutoff=1e-5)
        assert abs(eng.eigenvalue_cutoff - 1e-5) < 1e-12

    def test_compute_with_empty_modes_zero(self):
        """No modes -> zero vibrational entropy."""
        import flexaidds._core as _core
        eng = _core.ENCoMEngine()
        result = eng.compute_vibrational_entropy([], temperature=300.0)
        assert result.S_vib_kcal_mol_K == 0.0

    def test_compute_with_positive_eigenvalues(self):
        """Positive eigenvalues -> positive vibrational entropy."""
        import flexaidds._core as _core
        eng = _core.ENCoMEngine()
        modes = [
            _core.NormalMode(eigenvalue=1.0, frequency=100.0),
            _core.NormalMode(eigenvalue=2.0, frequency=141.4),
            _core.NormalMode(eigenvalue=0.5, frequency=70.7),
        ]
        result = eng.compute_vibrational_entropy(modes, temperature=300.0)
        assert result.S_vib_kcal_mol_K >= 0.0

    def test_more_modes_higher_entropy(self):
        """Adding more vibrational modes should not decrease entropy."""
        import flexaidds._core as _core
        eng = _core.ENCoMEngine()

        def make_modes(n):
            return [_core.NormalMode(eigenvalue=float(i + 1), frequency=float((i + 1) * 50))
                    for i in range(n)]

        s3 = eng.compute_vibrational_entropy(make_modes(3), 300.0).S_vib_kcal_mol_K
        s6 = eng.compute_vibrational_entropy(make_modes(6), 300.0).S_vib_kcal_mol_K
        assert s6 >= s3

    def test_higher_temperature_higher_entropy(self):
        """Vibrational entropy increases with temperature."""
        from flexaidds import ENCoMEngine, NormalMode
        eng = ENCoMEngine()
        modes = [NormalMode(eigenvalue=1.0, frequency=100.0)]
        s_low = eng.compute_vibrational_entropy(modes, 200.0).S_vib_kcal_mol_K
        s_high = eng.compute_vibrational_entropy(modes, 400.0).S_vib_kcal_mol_K
        assert s_high >= s_low


# ── Python-level ENCoM smoke test (no C++ needed) ────────────────────────────

class TestENCoMPythonFallback:
    """Verify that ENCoM symbols are always available (pure-Python fallback)."""

    def test_encom_symbols_accessible(self):
        import flexaidds as fds
        # Should always be importable as usable classes (never None)
        assert hasattr(fds, "ENCoMEngine")
        assert hasattr(fds, "NormalMode")
        assert hasattr(fds, "VibrationalEntropy")

    def test_encom_symbols_fallback_when_no_core(self):
        if _CORE_AVAILABLE:
            pytest.skip("C++ core is built; checking fallback path not applicable")
        import flexaidds as fds
        from flexaidds.encom import ENCoMEngine, NormalMode, VibrationalEntropy
        assert fds.ENCoMEngine is ENCoMEngine
        assert fds.NormalMode is NormalMode
        assert fds.VibrationalEntropy is VibrationalEntropy
        # Verify they are actually usable
        mode = fds.NormalMode(index=1, eigenvalue=2.0)
        assert mode.eigenvalue == 2.0
        vs = fds.ENCoMEngine.total_entropy(0.003, 0.001)
        assert abs(vs - 0.004) < 1e-15
