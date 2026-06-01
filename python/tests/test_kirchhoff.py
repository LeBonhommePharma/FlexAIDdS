"""Kirchhoff / Robertson-Murphy dG(T) temperature extrapolation (Tm form)."""
import math
import pytest
from flexaidds.thermodynamics import (
    StabilityCurve, gibbs_helmholtz_dG, kirchhoff_dH, kirchhoff_dS,
)


def test_delta_g_zero_at_tm():
    s = StabilityCurve(Tm=320.0, dHm=100.0, dCp=2.0)
    assert abs(gibbs_helmholtz_dG(320.0, s)) < 1e-12


def test_zero_dcp_is_vant_hoff_line():
    s = StabilityCurve(Tm=320.0, dHm=100.0, dCp=0.0)
    for T in (298.15, 310.0, 330.0):
        assert math.isclose(gibbs_helmholtz_dG(T, s), 100.0 * (1 - T / 320.0), abs_tol=1e-12)


def test_delta_g_equals_dh_minus_tds():
    s = StabilityCurve(Tm=320.0, dHm=100.0, dCp=2.0)
    for T in (290.0, 298.15, 305.0, 320.0, 340.0):
        assert math.isclose(gibbs_helmholtz_dG(T, s),
                            kirchhoff_dH(T, s) - T * kirchhoff_dS(T, s), abs_tol=1e-10)


def test_numeric_reference():
    s = StabilityCurve(Tm=320.0, dHm=100.0, dCp=2.0)
    assert math.isclose(gibbs_helmholtz_dG(298.15, s), 5.301013297963166, abs_tol=1e-9)
    assert math.isclose(kirchhoff_dH(298.15, s), 56.3, abs_tol=1e-9)


def test_rejects_nonpositive_temperature():
    s = StabilityCurve(Tm=320.0, dHm=100.0, dCp=2.0)
    with pytest.raises(ValueError):
        gibbs_helmholtz_dG(-1.0, s)
    with pytest.raises(ValueError):
        kirchhoff_dS(0.0, s)


def test_top_level_export():
    import flexaidds
    assert flexaidds.StabilityCurve is not None
    assert flexaidds.gibbs_helmholtz_dG is not None
