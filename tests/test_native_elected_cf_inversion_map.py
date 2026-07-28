"""Unit tests for scripts/native_elected_cf_inversion_map.py resolvers + classify."""
from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "native_elected_cf_inversion_map.py"


def _load_mod():
    spec = importlib.util.spec_from_file_location("inv_map", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["inv_map"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.skipif(not SCRIPT.is_file(), reason="inversion map script missing")
def test_classify_locked_miss_tied():
    m = _load_mod()
    # lower CF better; elected much better → LOCKED
    assert m.classify(cf_native=-10.0, cf_elected=-50.0, eps=0.5) == "SCORING-LOCKED"
    # native better → MISS
    assert m.classify(cf_native=-50.0, cf_elected=-10.0, eps=0.5) == "SEARCH-MISS"
    # within eps → TIED
    assert m.classify(cf_native=-10.0, cf_elected=-10.2, eps=0.5) == "TIED"


@pytest.mark.skipif(not SCRIPT.is_file(), reason="inversion map script missing")
def test_resolve_elected_prefers_result_csv_elected_path(tmp_path: Path):
    m = _load_mod()
    pdb = "1G9V"
    leaf = tmp_path / pdb
    leaf.mkdir()
    # wrong fallback that would mislead if preferred first
    (leaf / f"{pdb}_r0_0.pdb").write_text("R0\n", encoding="utf-8")
    elect = leaf / f"{pdb}_r9_0.pdb"
    elect.write_text("ELECT\n", encoding="utf-8")
    with (leaf / "result.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "elected_path", "rmsd_top1", "rmsd_bcr", "score_top1"],
        )
        w.writeheader()
        w.writerow(
            {
                "pdb_id": pdb,
                "elected_path": str(elect),
                "rmsd_top1": "11.0",
                "rmsd_bcr": "5.9",
                "score_top1": "-100.0",
            }
        )
    got = m.resolve_elected(tmp_path, pdb)
    assert got is not None
    assert got.resolve() == elect.resolve()


@pytest.mark.skipif(not SCRIPT.is_file(), reason="inversion map script missing")
def test_resolve_elected_appends_pdb_suffix(tmp_path: Path):
    m = _load_mod()
    pdb = "1P62"
    leaf = tmp_path / pdb
    leaf.mkdir()
    elect = leaf / f"{pdb}_r1_0.pdb"
    elect.write_text("E\n", encoding="utf-8")
    # path without suffix
    bare = str(elect.with_suffix(""))
    with (leaf / "result.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["elected_path"])
        w.writeheader()
        w.writerow({"elected_path": bare})
    got = m.resolve_elected(tmp_path, pdb)
    assert got is not None and got.is_file()
    assert got.name.endswith(".pdb")


@pytest.mark.skipif(not SCRIPT.is_file(), reason="inversion map script missing")
def test_load_result_metrics_pilot8_columns(tmp_path: Path):
    m = _load_mod()
    pdb = "1GPK"
    leaf = tmp_path / pdb
    leaf.mkdir()
    with (leaf / "result.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["rmsd_top1", "rmsd_bcr", "score_top1", "seed_echo", "elected_path"]
        )
        w.writeheader()
        w.writerow(
            {
                "rmsd_top1": "8.5",
                "rmsd_bcr": "4.4",
                "score_top1": "-900.1",
                "seed_echo": "0",
                "elected_path": "/x/y.pdb",
            }
        )
    meta = m.load_result_metrics(tmp_path, pdb)
    assert meta["rmsd_hungarian"] == 8.5
    assert meta["best_cluster_rmsd"] == 4.4
    assert meta["elected_cf"] == -900.1
    assert meta["seed_echo"] == 0.0
