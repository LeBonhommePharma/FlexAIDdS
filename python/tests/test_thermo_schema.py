"""Tests for the A1.1 entropy.help audit schema (thermo_audit.py).

Pure-Python CI jobs (no C++ bindings required).
"""

from __future__ import annotations

import math

import pytest

from flexaidds.schemas.thermo_audit import (
    ProvenanceDC,
    TotalSampledPartitionFunctionDC,
    ThermodynamicOutputDC,
    kB_kcal,
    make_total_sampled_output,
)


def _factory(**overrides):
    kwargs = dict(
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
    kwargs.update(overrides)
    return make_total_sampled_output(**kwargs)


def test_make_and_validate_roundtrip():
    out = _factory()

    assert isinstance(out, ThermodynamicOutputDC)
    assert out.total_sampled.logZ_total_sampled == -42.0
    assert out.provenance.git_sha == "deadbeef1234"
    assert out.n_samples_raw == 5000

    d = out.to_dict()
    out2 = ThermodynamicOutputDC.from_dict(d)
    out2.validate()

    assert out2.total_sampled.F_config_kcal_mol == pytest.approx(
        out.total_sampled.F_config_kcal_mol
    )
    assert out2.raw_ensemble_digest is None


def test_factory_identities_match_statmech():
    logZ = -12.5
    mean_E = -18.0
    T = 298.15
    out = _factory(logZ=logZ, mean_energy=mean_E, temperature_K=T, n_samples=10)
    F = -kB_kcal * T * logZ
    S = (mean_E - F) / T
    assert out.total_sampled.F_config_kcal_mol == pytest.approx(F)
    assert out.total_sampled.S_config_kcal_mol_K == pytest.approx(S)
    assert out.total_sampled.H_eff_kcal_mol == pytest.approx(mean_E)


def test_provenance_extras_and_digest():
    out = _factory(
        raw_ensemble_digest="sha256:abc",
        seed=42,
        runner_info="pytest",
        engine_version="2.0.0",
    )
    assert out.raw_ensemble_digest == "sha256:abc"
    assert out.provenance.seed == 42
    assert out.provenance.runner_info == "pytest"
    assert out.provenance.engine_version == "2.0.0"

    d = out.to_dict()
    assert d["provenance"]["seed"] == 42
    roundtrip = ThermodynamicOutputDC.from_dict(d)
    roundtrip.validate()


def test_validation_rejects_bad_consistency():
    good = _factory(n_samples=100, git_sha="abc", timestamp="2026-01-01", gate_results={})
    good.validate()

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
    with pytest.raises(ValueError, match="F_config inconsistency"):
        bad.validate()


def test_validation_rejects_bad_temperature_and_samples():
    base = _factory().to_dict()
    base["temperature_K"] = 0.0
    with pytest.raises(ValueError, match="temperature_K"):
        ThermodynamicOutputDC.from_dict(base).validate()

    base = _factory().to_dict()
    base["n_samples_raw"] = 0
    with pytest.raises(ValueError, match="n_samples_raw"):
        ThermodynamicOutputDC.from_dict(base).validate()


def test_validation_rejects_s_mismatch():
    d = _factory().to_dict()
    d["total_sampled"]["S_config_kcal_mol_K"] = 999.0
    with pytest.raises(ValueError, match="S_config"):
        ThermodynamicOutputDC.from_dict(d).validate()


def test_validation_rejects_nonfinite():
    d = _factory().to_dict()
    d["total_sampled"]["logZ_total_sampled"] = float("nan")
    # F identity will also fail; either error is fine
    with pytest.raises(ValueError):
        ThermodynamicOutputDC.from_dict(d).validate()

    d = _factory().to_dict()
    d["total_sampled"]["H_eff_kcal_mol"] = float("inf")
    # Fix F/S to match logZ so we hit the finite check
    T = d["temperature_K"]
    logZ = d["total_sampled"]["logZ_total_sampled"]
    F = -kB_kcal * T * logZ
    d["total_sampled"]["F_config_kcal_mol"] = F
    d["total_sampled"]["S_config_kcal_mol_K"] = (float("inf") - F) / T
    with pytest.raises(ValueError, match="finite"):
        ThermodynamicOutputDC.from_dict(d).validate()


def test_typed_dc_roundtrips():
    ts = TotalSampledPartitionFunctionDC(
        logZ_total_sampled=-1.0,
        F_config_kcal_mol=0.5,
        H_eff_kcal_mol=-2.0,
        S_config_kcal_mol_K=0.01,
    )
    assert TotalSampledPartitionFunctionDC.from_dict(ts.to_dict()) == ts

    prov = ProvenanceDC(
        temperature_K=300.0,
        n_samples=10,
        git_sha="abc",
        timestamp="t",
        gate_results={"g": True},
    )
    assert ProvenanceDC.from_dict(prov.to_dict()).git_sha == "abc"
    assert math.isclose(prov.temperature_K, 300.0)
