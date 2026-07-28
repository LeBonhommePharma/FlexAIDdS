"""Tests for scripts/build_pilot8_summary.py against real CSV schema rows."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "build_pilot8_summary.py"
MODE_COLS = [f"mode_rmsd_{i}" for i in range(10)]


def _write_result(path: Path, *, bcr: float, modes: list[float], s_top10: int = 0) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "arm": "A",
        "engine_sha": "deadbeef",
        "matrix_md5": "9dc93717dfed0698006d88dd6a9627bc",
        "pdb_id": path.parent.name,
        "rmsd_top1": f"{modes[0]:.4f}",
        "rmsd_bcr": f"{bcr:.4f}",
        "success_s_top10": str(s_top10),
        "seed_echo": "0",
        "native_pose_seeded": "0",
    }
    for i, m in enumerate(modes):
        row[f"mode_rmsd_{i}"] = f"{m:.4f}"
    fieldnames = list(row.keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerow(row)


@pytest.mark.skipif(not SCRIPT.is_file(), reason="build_pilot8_summary.py missing")
def test_build_pilot8_summary_schema_and_counts(tmp_path: Path) -> None:
    arm = tmp_path / "A"
    # 1GPK: schema ok, no success
    _write_result(
        arm / "1GPK" / "result.csv",
        bcr=5.3,
        modes=[6.8, 8.9, 7.3, 9.1, 8.5, 9.5, 9.6, 8.0, 9.4, 9.5],
        s_top10=0,
    )
    # 1P62: schema ok, BCR success
    _write_result(
        arm / "1P62" / "result.csv",
        bcr=1.5,
        modes=[1.2, 3.0, 4.0, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6],
        s_top10=1,
    )
    out = tmp_path / "summary.json"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--arm-dir",
            str(arm),
            "--panel",
            "1GPK",
            "1P62",
            "1MEH",
            "--json-out",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    summary = json.loads(out.read_text())
    assert summary["n_targets"] == 2
    assert summary["n_missing"] == 1
    assert summary["schema_ok"] is True
    assert summary["bcr_success"] == 1
    assert summary["s_top10_success"] == 1
    assert "1MEH" in summary["missing_pdbs"]
    # incomplete panel must science_hold so P3 cannot pass open P4 early
    assert summary["science_hold"] is True
    # mode columns present on parsed targets
    for t in summary["arms"][0]["targets"]:
        if t["status"] == "ok":
            assert t["schema_ok"] is True


@pytest.mark.skipif(not SCRIPT.is_file(), reason="build_pilot8_summary.py missing")
def test_build_pilot8_summary_schema_fail_missing_modes(tmp_path: Path) -> None:
    arm = tmp_path / "A"
    path = arm / "1GPK" / "result.csv"
    path.parent.mkdir(parents=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb_id", "rmsd_bcr", "rmsd_top1"])
        w.writeheader()
        w.writerow({"pdb_id": "1GPK", "rmsd_bcr": "1.0", "rmsd_top1": "1.0"})
    out = tmp_path / "summary.json"
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--arm-dir",
            str(arm),
            "--panel",
            "1GPK",
            "--json-out",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    summary = json.loads(out.read_text())
    assert summary["n_targets"] == 1
    assert summary["schema_ok"] is False
