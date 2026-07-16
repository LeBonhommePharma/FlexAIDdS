#!/usr/bin/env python3
"""Unit tests for scripts/native_cf_oracle_gate.py."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"


def _load():
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    path = SCRIPTS / "native_cf_oracle_gate.py"
    spec = importlib.util.spec_from_file_location("native_cf_oracle_gate", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["native_cf_oracle_gate"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def oracle():
    return _load()


def test_native_competitive(oracle):
    # CF_native=-100, best_ga=-95 → native better → PASS
    res = oracle.evaluate_gate(
        cf_native=-100.0,
        best_ga_cf=-95.0,
        tolerance=0.0,
    )
    assert res.ok
    assert res.exit_code == oracle.EXIT_PASS
    assert res.ranking_forbidden is False


def test_native_pathology_1p62_style(oracle):
    # Crystal CF ≈ -112, GA best ≈ -657 → native not competitive
    res = oracle.evaluate_gate(
        cf_native=-112.0,
        best_ga_cf=-657.0,
        tolerance=0.0,
    )
    assert not res.ok
    assert res.exit_code == oracle.EXIT_FAIL_PATHOLOGY
    assert res.ranking_forbidden is True
    assert res.gap is not None and res.gap < 0


def test_tolerance_allows_near_miss(oracle):
    # native -100, best -102, tol 5 → PASS
    res = oracle.evaluate_gate(
        cf_native=-100.0,
        best_ga_cf=-102.0,
        tolerance=5.0,
    )
    assert res.ok


def test_sentinel_ini(oracle):
    res = oracle.evaluate_gate(
        cf_native=10000.0,
        best_ga_cf=-50.0,
        tolerance=0.0,
    )
    assert not res.ok
    assert res.exit_code == oracle.EXIT_MISSING_NATIVE
    assert res.ranking_forbidden is True


def test_missing_native(oracle):
    res = oracle.evaluate_gate(
        cf_native=None,
        best_ga_cf=-50.0,
        tolerance=0.0,
    )
    assert res.exit_code == oracle.EXIT_MISSING_NATIVE


def test_fixture_dir(oracle, tmp_path: Path):
    d = tmp_path / "1ABC"
    d.mkdir()
    (d / "1ABC_INI.pdb").write_text(
        "REMARK initial structure\nREMARK CF=-60.0\nEND\n"
    )
    (d / "1ABC_0.pdb").write_text(
        "REMARK optimized structure\nREMARK CF=-80.0\nEND\n"
    )
    (d / "1ABC_1.pdb").write_text(
        "REMARK optimized structure\nREMARK CF=-70.0\nEND\n"
    )
    (d / "result.csv").write_text(
        "pdb_id,best_score,cf_native\n1ABC,-80.0,-60.0\n"
    )
    res = oracle.run_gate(results=d, pdb_id="1ABC", tolerance=0.0)
    assert not res.ok
    assert res.exit_code == oracle.EXIT_FAIL_PATHOLOGY
    assert res.cf_native == pytest.approx(-60.0)
    assert res.best_ga_cf == pytest.approx(-80.0)


def test_cli_json(oracle, tmp_path: Path):
    d = tmp_path / "1ABC"
    d.mkdir()
    (d / "1ABC_INI.pdb").write_text("REMARK CF=-90.0\nEND\n")
    (d / "1ABC_0.pdb").write_text("REMARK CF=-85.0\nEND\n")
    out = tmp_path / "gate.json"
    rc = oracle.main(
        ["--results", str(d), "--pdb", "1ABC", "--json", str(out), "-q"]
    )
    assert rc == 0
    data = json.loads(out.read_text())
    assert data["ok"] is True
    assert data["ranking_forbidden"] is False


def test_is_sentinel(oracle):
    assert oracle.is_sentinel_cf(10000.0) is True
    assert oracle.is_sentinel_cf(-12.5) is False
    assert oracle.is_sentinel_cf(None) is True


def test_sentinel_ini_and_zero_csv(oracle, tmp_path: Path):
    """C0-style: INI CF=10000 and result.csv cf_native=0 must not count as native."""
    d = tmp_path / "1P62"
    d.mkdir()
    (d / "1P62_INI.pdb").write_text(
        "REMARK initial structure\nREMARK CF=10000.00000\nEND\n"
    )
    (d / "1P62_0.pdb").write_text("REMARK CF=-657.0\nEND\n")
    (d / "result.csv").write_text(
        "pdb_id,best_score,cf_native\n1P62,-657.0,0.0\n"
    )
    res = oracle.run_gate(results=d, pdb_id="1P62")
    assert not res.ok
    assert res.exit_code == oracle.EXIT_MISSING_NATIVE
    assert res.ranking_forbidden is True
