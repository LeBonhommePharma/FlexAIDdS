"""Smoke tests for the A1.1 entropy.help audit schema (thermo_audit.py).

These are intentionally lightweight so they run in pure-Python CI jobs
(no C++ bindings required).
"""

import pytest

from flexaidds.schemas.thermo_audit import (
    make_total_sampled_output,
    ThermodynamicOutputDC,
)


def test_make_and_validate_roundtrip():
    out = make_total_sampled_output(
        logZ=-42.0,
        mean_energy=-15.5,
        temperature_K=300.0,
        n_samples=5000,
        git_sha="deadbeef1234",
        timestamp="2026-05-27T12:00:00Z",
        gate_results={
            "gate5_convergence": {"passed": True, "delta_logZ": 1e-11},
            "gate6_crosscheck": {"passed": True},
        },
    )

    assert isinstance(out, ThermodynamicOutputDC)
    assert out.total_sampled.logZ_total_sampled == -42.0
    assert out.provenance.git_sha == "deadbeef1234"

    d = out.to_dict()
    out2 = ThermodynamicOutputDC.from_dict(d)
    out2.validate()

    assert out2.total_sampled.F_config_kcal_mol == pytest.approx(out.total_sampled.F_config_kcal_mol)


def test_validation_rejects_bad_consistency():
    with pytest.raises(ValueError, match="F_config inconsistency"):
        make_total_sampled_output(
            logZ=-42.0,
            mean_energy=-15.5,
            temperature_K=300.0,
            n_samples=100,
            git_sha="abc",
            timestamp="2026-01-01",
            gate_results={},
            # Deliberately wrong F by passing a bad value via direct construction
        )
    # The factory itself would have caught it; test direct bad object
    bad = ThermodynamicOutputDC.from_dict(
        {
            "total_sampled": {
                "logZ_total_sampled": -42.0,
                "F_config_kcal_mol": 999.0,  # wrong
                "H_eff_kcal_mol": -15.5,
                "S_config_kcal_mol_K": 0.01,
            },
            "temperature_K": 300.0,
            "n_samples_raw": 100,
            "provenance": {
                "temperature_K": 300.0,
                "n_samples": 100,
                "git_sha": "x",
                "timestamp": "t",
                "gate_results": {},
            },
        }
    )
    with pytest.raises(ValueError):
        bad.validate()
