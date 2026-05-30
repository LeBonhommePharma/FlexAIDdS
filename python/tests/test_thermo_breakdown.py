from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from flexaidds.models import BindingModeResult, DockingResult
from flexaidds.thermodynamics import StatMechEngine, ThermodynamicBreakdown, kB_kcal


def test_compute_breakdown_single_state_identities():
    engine = StatMechEngine(300.0)
    engine.add_sample(-12.0)

    b = engine.compute_breakdown()

    assert b.G_config_kcal_mol == pytest.approx(-kB_kcal * 300.0 * b.logZ_config)
    assert b.H_eff_kcal_mol == pytest.approx(-12.0)
    assert b.S_config_kcal_mol_K == pytest.approx(0.0)
    assert b.minus_T_S_config_kcal_mol == pytest.approx(b.G_config_kcal_mol - b.H_eff_kcal_mol)
    assert b.G_total_kcal_mol == pytest.approx(b.G_config_kcal_mol)


def test_compute_breakdown_two_equal_states():
    engine = StatMechEngine(300.0)
    engine.add_sample(-10.0)
    engine.add_sample(-10.0)

    b = engine.compute_breakdown()

    assert b.G_config_kcal_mol == pytest.approx(-10.0 - kB_kcal * 300.0 * math.log(2.0))
    assert b.H_eff_kcal_mol == pytest.approx(-10.0)
    assert b.S_config_kcal_mol_K == pytest.approx(kB_kcal * math.log(2.0))
    assert b.Cv_kcal_mol_K == pytest.approx(0.0)


def test_compute_breakdown_correction_sum_and_flags():
    engine = StatMechEngine(300.0)
    engine.add_sample(-10.0)
    thermo = engine.compute()

    b = engine.compute_breakdown(
        G_vib_kcal_mol=-0.25,
        G_natural_kcal_mol=0.1,
        G_other_kcal_mol=-0.05,
        has_vib=True,
        has_natural=True,
        has_other=True,
    )

    assert b.G_config_kcal_mol == pytest.approx(thermo.free_energy)
    assert b.G_total_kcal_mol == pytest.approx(thermo.free_energy - 0.25 + 0.1 - 0.05)
    assert b.has_vib is True
    assert b.has_natural is True
    assert b.has_other is True


def test_python_multiplicity_is_log_weighted_not_integer_duplicated():
    engine = StatMechEngine(300.0)
    engine.add_sample(-10.0, 0.5)
    engine.add_sample(-8.0, 1.5)

    b = engine.compute_breakdown()

    beta = 1.0 / (kB_kcal * 300.0)
    log_w = [math.log(0.5) - beta * -10.0, math.log(1.5) - beta * -8.0]
    m = max(log_w)
    log_z = m + math.log(sum(math.exp(x - m) for x in log_w))
    expected_g = -kB_kcal * 300.0 * log_z
    assert b.logZ_config == pytest.approx(log_z)
    assert b.G_config_kcal_mol == pytest.approx(expected_g)


def test_thermodynamic_breakdown_json_round_trip():
    breakdown = ThermodynamicBreakdown(
        temperature_K=300.0,
        logZ_config=2.0,
        G_config_kcal_mol=-1.2,
        H_eff_kcal_mol=-1.0,
        S_config_kcal_mol_K=0.001,
        minus_T_S_config_kcal_mol=-0.2,
        Cv_kcal_mol_K=0.01,
        sigma_E_kcal_mol=0.3,
        G_vib_kcal_mol=-0.1,
        G_total_kcal_mol=-1.3,
        has_vib=True,
    )

    mode = BindingModeResult(
        mode_id=1,
        rank=1,
        poses=[],
        free_energy=-1.3,
        thermodynamics=breakdown,
    )
    result = DockingResult(source_dir=Path("/tmp/flex"), binding_modes=[mode])

    payload = json.loads(result.to_json())
    nested = payload["binding_modes"][0]["thermodynamics"]
    assert nested["G_config_kcal_mol"] == pytest.approx(-1.2)
    assert nested["G_total_kcal_mol"] == pytest.approx(-1.3)

    restored = DockingResult.from_json(json.dumps(payload))
    assert restored.binding_modes[0].free_energy == pytest.approx(-1.3)
    assert restored.binding_modes[0].thermodynamics is not None
    assert restored.binding_modes[0].thermodynamics.G_config_kcal_mol == pytest.approx(-1.2)
