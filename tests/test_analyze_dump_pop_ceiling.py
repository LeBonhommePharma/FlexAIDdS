"""Drive scripts/analyze_dump_pop_ceiling.py on real engine .pop.tsv schema."""
from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_dump_pop_ceiling.py"
FIXTURE = ROOT / "tests" / "fixtures" / "dump_pop" / "sample.pop.tsv"


def _load_mod():
    spec = importlib.util.spec_from_file_location("analyze_dump_pop_ceiling", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def test_fixture_has_engine_dump_columns():
    rows = list(csv.DictReader(FIXTURE.open(), delimiter="\t"))
    assert rows
    assert set(rows[0]) >= {
        "idx", "cluster", "rmsd_to_head", "rmsd_raw", "rmsd_sym",
        "cf_total", "cf_com", "cf_wal", "pose_id", "is_elected",
    }


def test_summarize_pop_uses_engine_columns_not_hardcoded_mins():
    mod = _load_mod()
    rows = mod.load_pop(FIXTURE)
    s = mod.summarize_pop(rows)
    # derived from fixture rows (if fixture changes, expectation follows data)
    assert s["n_dump"] == 3
    assert s["min_rmsd_raw"] == pytest.approx(3.1)
    assert s["min_rmsd_sym"] == pytest.approx(2.8)
    assert s["n_sym_le2"] == 0
    assert s["n_raw_le2"] == 0


def test_compare_search_wall_when_dump_near_bcr(tmp_path: Path):
    mod = _load_mod()
    # write a mini OUT
    code = "1XXX"
    d = tmp_path / code
    d.mkdir()
    (d / f"{code}.pop.tsv").write_text(FIXTURE.read_text())
    # BCR near min_sym 2.8
    (d / "result.csv").write_text(
        "pdb_id,best_cluster_rmsd,conditional_scanned_pool_ceiling,rmsd_to_crystal\n"
        "1XXX,2.85,2.85,6.0\n"
    )
    rec = mod.compare(mod.summarize_pop(mod.load_pop(d / f"{code}.pop.tsv")),
                      mod.load_emission_bcr(d / "result.csv"))
    assert rec["verdict"] == "SEARCH_WALL"


def test_compare_retention_when_sub2_in_dump(tmp_path: Path):
    mod = _load_mod()
    code = "1YYY"
    d = tmp_path / code
    d.mkdir()
    # one sub-2 sym row
    (d / f"{code}.pop.tsv").write_text(
        "idx\tcluster\trmsd_to_head\trmsd_raw\trmsd_sym\tcf_total\tcf_com\tcf_wal\tpose_id\tis_elected\n"
        "0\t0\t0\t1.5\t1.2\t-1\t-1\t0\t0\t0\n"
        "1\t0\t0\t4.0\t3.5\t-2\t-2\t0\t1\t1\n"
    )
    (d / "result.csv").write_text(
        "pdb_id,best_cluster_rmsd,conditional_scanned_pool_ceiling,rmsd_to_crystal\n"
        "1YYY,3.5,3.5,5.0\n"
    )
    rec = mod.compare(
        mod.summarize_pop(mod.load_pop(d / f"{code}.pop.tsv")),
        mod.load_emission_bcr(d / "result.csv"),
    )
    assert rec["verdict"] == "RETENTION"
    assert rec["n_sym_le2"] == 1


def test_cli_on_fixture_tree(tmp_path: Path):
    mod = _load_mod()
    code = "1ZZZ"
    d = tmp_path / code
    d.mkdir()
    (d / f"{code}.pop.tsv").write_text(FIXTURE.read_text())
    (d / "result.csv").write_text(
        "pdb_id,best_cluster_rmsd,conditional_scanned_pool_ceiling,rmsd_to_crystal\n"
        "1ZZZ,4.0,4.0,6.0\n"
    )
    rc = mod.main([str(tmp_path), "--codes", code])
    assert rc == 0
