"""Tests for pure-Python fallback stubs in flexaidds._fallback_types."""

from __future__ import annotations

import math

import pytest

from flexaidds._fallback_types import (
    BoltzmannLUT,
    DeltaCpFit,
    Replica,
    State,
    TemperatureScanPoint,
    TIPoint,
    WHAMBin,
)


def test_wham_bin_defaults():
    b = WHAMBin()
    assert b.coord_center == 0.0
    assert b.free_energy == 0.0
    assert b.count == 0


def test_ti_point_and_replica():
    ti = TIPoint(lambda_val=0.5, dV_dlambda=-1.2)
    assert ti.lambda_val == 0.5
    assert ti.dV_dlambda == pytest.approx(-1.2)

    r = Replica(id=3, temperature=310.0, beta=1.0 / 310.0, current_energy=-12.0)
    assert r.id == 3
    assert r.temperature == 310.0
    assert r.current_energy == -12.0


def test_state_defaults():
    s = State(energy=-8.5)
    assert s.energy == -8.5
    assert s.count == 1


def test_boltzmann_lut_lookup():
    lut = BoltzmannLUT(temperature=300.0)
    # At E=0, exp(0)=1
    assert lut.lookup(0.0) == pytest.approx(1.0)
    # Lower energy → higher Boltzmann weight
    assert lut.lookup(-5.0) > lut.lookup(0.0) > lut.lookup(5.0)
    # Positive energy factor matches manual computation
    kB = 0.0019872041
    beta = 1.0 / (kB * 300.0)
    assert lut.lookup(2.0) == pytest.approx(math.exp(-beta * 2.0))


def test_boltzmann_lut_zero_temperature():
    lut = BoltzmannLUT(temperature=0.0)
    # beta guard: returns exp(0) = 1 when T=0
    assert lut.lookup(10.0) == pytest.approx(1.0)


def test_temperature_scan_point_and_delta_cp_fit():
    pt = TemperatureScanPoint(
        T_K=298.15,
        logZ=10.0,
        G_kcal_mol=-5.9,
        H_kcal_mol=-8.0,
        S_kcal_mol_K=0.007,
        Cv_kcal_mol_K=0.1,
    )
    assert pt.T_K == 298.15
    assert pt.G_kcal_mol == pytest.approx(-5.9)

    fit = DeltaCpFit(
        delta_Cp_kcal_mol_K=0.25,
        T_ref_K=298.15,
        rmse_kcal_mol=0.05,
        model_derived=True,
        experimental=True,
    )
    assert fit.delta_Cp_kcal_mol_K == 0.25
    assert fit.model_derived is True
    assert fit.experimental is True
