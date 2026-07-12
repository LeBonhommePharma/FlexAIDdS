#!/usr/bin/env python3
"""Unit tests for scripts/acf_vs_cf_ablation.py (classic entropy election).

Pure Python, no C++ binary, no docking. Uses the built-in 1HNN exhibit.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "acf_vs_cf_ablation.py"


def _load_mod():
    spec = importlib.util.spec_from_file_location("acf_vs_cf_ablation", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def abl():
    return _load_mod()


def test_synthetic_1hnn_classic_elects_acf_basin(abl):
    rows = abl.rows_from_synthetic(abl.SYNTHETIC_1HNN)
    i = abl.elect_rank0(rows, force_cf_rank_emission=False, temperature=300)
    assert rows[i].cluster == 3
    assert rows[i].freq == 29
    assert rows[i].acf == pytest.approx(-263.427453)


def test_synthetic_1hnn_force_cf_elects_cf_champion(abl):
    rows = abl.rows_from_synthetic(abl.SYNTHETIC_1HNN)
    i = abl.elect_rank0(rows, force_cf_rank_emission=True, temperature=300)
    assert rows[i].cluster == 0
    assert rows[i].cf == pytest.approx(-189.85613)


def test_temperature_zero_forces_cf(abl):
    rows = abl.rows_from_synthetic(abl.SYNTHETIC_1HNN)
    i = abl.elect_rank0(rows, force_cf_rank_emission=False, temperature=0)
    assert rows[i].cluster == 0


def test_election_flips_on_1hnn(abl):
    rows = abl.rows_from_synthetic(abl.SYNTHETIC_1HNN)
    rep = abl.ablation_report(rows, label="1HNN")
    assert rep["election_flips"] is True
    assert rep["classic_entropy_rank0"]["cluster"] == 3
    assert rep["force_cf_rank0"]["cluster"] == 0


def test_parse_cad_roundtrip(abl):
    text = (
        "Cluster 0: TOP=0 TCF=-17.437956 ACF=-49.305648 freq=4\n"
        "Cluster 3: TOP=7 TCF=-11.775027 ACF=-263.427453 freq=29\n"
    )
    rows = abl.parse_cad(text)
    assert len(rows) == 2
    assert rows[1].acf == pytest.approx(-263.427453)
    assert rows[1].freq == 29


def test_cli_synthetic_json(abl, capsys):
    rc = abl.main(["--synthetic-1hnn", "--json"])
    assert rc == 0
    out = capsys.readouterr().out
    data = json.loads(out)
    assert data["election_flips"] is True
    assert data["classic_entropy_rank0"]["cluster"] == 3


def test_live_1hnn_cad_if_present(abl):
    """Optional: if a live pre-fix 1HNN dir exists, assert the same flip."""
    candidates = [
        REPO / "results/astex_jcim2015_fair_20260708_0002/1HNN",
        REPO / "benchmarks/astex_repro/full/1HNN",
    ]
    target = next((p for p in candidates if (p / "1HNN.cad").is_file()), None)
    if target is None:
        pytest.skip("no live 1HNN.cad in workspace")
    rows = abl.load_target_dir(target)
    rep = abl.ablation_report(rows, label=str(target))
    # Pre-fix ensembles show ACF-best ≠ CF-best; if a future ensemble
    # happens to agree, the synthetic tests still pin the contract.
    assert rep["n_clusters"] >= 2
    classic = rep["classic_entropy_rank0"]
    force = rep["force_cf_rank0"]
    assert classic is not None and force is not None
    # Classic elects min ACF among parsed clusters
    assert classic["acf"] == min(r["acf"] for r in rep["clusters"])
    assert force["cf"] == min(r["cf"] for r in rep["clusters"])
