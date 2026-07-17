"""Unit tests for scripts/bootstrap_3dsig_s_top10.py (fail-closed S_top10)."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bootstrap_3dsig_s_top10 as boot  # noqa: E402


def test_s_top10_any_of_ten():
    assert boot.s_top10([3.0, 1.5] + [5.0] * 8) is True
    assert boot.s_top10([2.0] + [5.0] * 9) is True  # inclusive ≤2.0
    assert boot.s_top10([2.01] + [5.0] * 9) is False
    assert boot.s_top10([None] * 10) is False


def test_extract_mode_rmsds_fail_closed_on_bcr_only(tmp_path: Path):
    row = {"rmsd_bcr": "1.1", "rmsd_top1": "1.2", "success_s1": "1"}
    with pytest.raises(boot.MissingModeRmsdError):
        boot.extract_mode_rmsds_from_row(row, require_mode_columns=True)


def test_extract_mode_rmsds_prefers_mode_columns():
    row = {f"mode_rmsd_{i}": str(3.0 + i * 0.1) for i in range(10)}
    row["mode_rmsd_3"] = "1.9"
    row["rmsd_bcr"] = "0.5"  # must not be used by s_top10 path alone
    modes = boot.extract_mode_rmsds_from_row(row)
    assert boot.s_top10(modes) is True
    assert modes[3] == pytest.approx(1.9)


def test_load_arm_dir_and_bootstrap(tmp_path: Path):
    # two successes, two fails → observed 0.5
    specs = {
        "AAAA": [1.5] + [5.0] * 9,
        "BBBB": [5.0] * 10,
        "CCCC": [5.0, 1.0] + [5.0] * 8,
        "DDDD": [5.0] * 10,
    }
    for pdb, modes in specs.items():
        d = tmp_path / pdb
        d.mkdir()
        path = d / "result.csv"
        fields = ["pdb_id"] + [f"mode_rmsd_{i}" for i in range(10)] + ["rmsd_bcr", "success_s_top10"]
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            row = {"pdb_id": pdb, "rmsd_bcr": "9.9", "success_s_top10": "0"}
            for i, m in enumerate(modes):
                row[f"mode_rmsd_{i}"] = str(m)
            row["success_s_top10"] = "1" if boot.s_top10(modes) else "0"
            w.writerow(row)

    cases = boot.load_arm_dir(tmp_path)
    assert cases["AAAA"] is True
    assert cases["BBBB"] is False
    assert cases["CCCC"] is True
    assert cases["DDDD"] is False

    summary = boot.summarize(cases, bootstraps=2000, seed=1, thresh=2.0, source="test")
    assert summary["n_cases"] == 4
    assert summary["n_hits"] == 2
    assert summary["observed_rate"] == pytest.approx(0.5)
    assert summary["bcr_not_used_as_s_top10"] is True
    assert 0.0 <= summary["bootstrap"]["median"] <= 1.0


def test_main_arm_dir_json(tmp_path: Path, capsys):
    d = tmp_path / "X111"
    d.mkdir()
    fields = [f"mode_rmsd_{i}" for i in range(10)]
    with (d / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow({f"mode_rmsd_{i}": "5.0" for i in range(10)})
    out = tmp_path / "out.json"
    rc = boot.main(["--arm-dir", str(tmp_path), "--bootstraps", "100", "--json-out", str(out)])
    assert rc == 0
    data = json.loads(out.read_text())
    assert data["n_hits"] == 0
    assert data["observed_rate"] == 0.0
