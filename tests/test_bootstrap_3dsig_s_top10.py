#!/usr/bin/env python3
"""Unit tests for S_top10 bootstrap + FlexAID arm mode_rmsd emission.

Covers:
  - scripts/bootstrap_3dsig_s_top10.py (≤2.0, fail-closed, 10k-capable helpers)
  - scripts/parse_flexaid_arm_results.py (mode_rmsd_0..9 in emitted rank order)

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BOOT = ROOT / "scripts" / "bootstrap_3dsig_s_top10.py"
PARSE = ROOT / "scripts" / "parse_flexaid_arm_results.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_pose(
    path: Path,
    *,
    rmsd: float,
    cf: float,
    cf_app: float | None = None,
    sym: bool = True,
) -> None:
    lines = [
        "REMARK   FlexAID test pose",
        f"REMARK   CF={cf:.4f}",
    ]
    if cf_app is not None:
        lines.append(f"REMARK   CF.app={cf_app:.4f}")
    if sym:
        lines.append(
            f"REMARK   {rmsd:.4f} RMSD to ref. structure (symmetry corrected)"
        )
    else:
        lines.append(
            f"REMARK   {rmsd:.4f} RMSD to ref. structure (no symmetry correction)"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


# ── bootstrap helpers ────────────────────────────────────────────────────────


def test_s_top10_inclusive_threshold():
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    # Exactly 2.0 must count as success (≤ 2.0 claim contract)
    assert mod.s_top10([2.0], thresh=2.0) is True
    assert mod.s_top10([2.0001], thresh=2.0) is False
    assert mod.s_top10([1.9, 3.0], thresh=2.0) is True
    assert mod.s_top10([3.0, 4.0], thresh=2.0) is False
    assert mod.s_top10([None, None], thresh=2.0) is False
    assert mod.s_top10([], thresh=2.0) is False
    # Only first 10 matter
    vals = [3.0] * 10 + [0.5]
    assert mod.s_top10(vals, thresh=2.0) is False


def test_extract_mode_rmsds_prefer_mode_rmsd_keys():
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    row = {f"mode_rmsd_{i}": str(i * 0.5) for i in range(10)}
    row["rmsd_top1"] = "9.9"
    row["rmsd_bcr"] = "0.1"
    got = mod.extract_mode_rmsds_from_row(row)
    assert got[0] == 0.0
    assert got[4] == 2.0
    assert got[9] == 4.5


def test_extract_mode_rmsds_fail_closed_without_mode_cols():
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    row = {
        "pdb_id": "1ABC",
        "rmsd_top1": "1.2",
        "rmsd_bcr": "0.8",
        "success_s1": "1",
        "success_s3": "1",
    }
    with pytest.raises(mod.MissingModeRmsdError):
        mod.extract_mode_rmsds_from_row(row, require_mode_columns=True)


def test_load_arm_dir_fail_closed(tmp_path: Path):
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    d = tmp_path / "arm" / "1ABC"
    d.mkdir(parents=True)
    with (d / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["pdb_id", "rmsd_top1", "rmsd_bcr", "success_s1"]
        )
        w.writeheader()
        w.writerow(
            {
                "pdb_id": "1ABC",
                "rmsd_top1": "1.1",
                "rmsd_bcr": "0.5",
                "success_s1": "1",
            }
        )
    with pytest.raises(mod.MissingModeRmsdError):
        mod.load_arm_dir(tmp_path / "arm", strict=True)


def test_load_arm_dir_reads_mode_rmsd(tmp_path: Path):
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    d = tmp_path / "arm" / "1ABC"
    d.mkdir(parents=True)
    fields = ["pdb_id"] + [f"mode_rmsd_{i}" for i in range(10)]
    row = {"pdb_id": "1ABC"}
    # mode 0 bad, mode 3 good → S_top10 success without S1
    rmsds = [3.0, 4.0, 5.0, 1.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    for i, v in enumerate(rmsds):
        row[f"mode_rmsd_{i}"] = f"{v:.4f}"
    with (d / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(row)
    cases = mod.load_arm_dir(tmp_path / "arm", strict=True)
    assert "1ABC" in cases
    assert cases["1ABC"][3] == pytest.approx(1.5)
    assert mod.s_top10(cases["1ABC"]) is True
    # S1 would fail (mode 0 = 3.0)
    assert cases["1ABC"][0] == pytest.approx(3.0)


def _write_case(
    arm: Path, pdb: str, rmsds, *, n_poses: str = "500", restarts: str = "10"
) -> None:
    """Write one per-target result.csv. ``rmsds`` entries may be None → empty."""
    d = arm / pdb
    d.mkdir(parents=True)
    fields = (
        ["pdb_id"]
        + [f"mode_rmsd_{i}" for i in range(10)]
        + ["n_poses", "restarts_finished"]
    )
    row = {"pdb_id": pdb, "n_poses": n_poses, "restarts_finished": restarts}
    for i, v in enumerate(rmsds):
        row[f"mode_rmsd_{i}"] = "" if v is None else f"{v:.4f}"
    with (d / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(row)


def test_load_arm_dir_fail_closed_on_never_executed_case(tmp_path: Path):
    """A row that never ran must not be scored as an S_top10 miss."""
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case(arm, "1ABC", [3.0] * 10)
    _write_case(arm, "1XYZ", [None] * 10, n_poses="0", restarts="0")
    with pytest.raises(mod.UnrunCaseError):
        mod.load_arm_dir(arm, strict=True)


def test_load_arm_dir_allow_unrun_excludes_rather_than_fails(tmp_path: Path):
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case(arm, "1ABC", [1.5] + [9.0] * 9)
    _write_case(arm, "1XYZ", [None] * 10, n_poses="0", restarts="0")
    cases = mod.load_arm_dir(arm, strict=True, allow_unrun=True)
    # Excluded from N entirely — not counted as a failure.
    assert set(cases) == {"1ABC"}
    stats = mod.compute_s_top10_stats(cases, n_boot=64)
    assert stats["n_cases"] == 1
    assert stats["point"] == pytest.approx(1.0)


def test_empty_modes_without_execution_witness_still_counts_as_miss(tmp_path: Path):
    """Absent witness columns are not evidence the run failed to start."""
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm" / "1ABC"
    arm.mkdir(parents=True)
    fields = ["pdb_id"] + [f"mode_rmsd_{i}" for i in range(10)]
    row = {"pdb_id": "1ABC"}
    for i in range(10):
        row[f"mode_rmsd_{i}"] = ""
    with (arm / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(row)
    cases = mod.load_arm_dir(tmp_path / "arm", strict=True)
    assert cases["1ABC"] == [None] * 10
    assert mod.s_top10(cases["1ABC"]) is False


def test_cli_arm_dir_fail_closed_on_unrun(tmp_path: Path):
    """The CLI must refuse to print 0.0000 for an arm that never executed."""
    arm = tmp_path / "arm"
    _write_case(arm, "1ABC", [None] * 10, n_poses="0", restarts="0")
    proc = subprocess.run(
        [sys.executable, str(BOOT), "--arm-dir", str(arm), "--bootstraps", "10"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 2
    assert "never executed" in proc.stderr
    assert "0.0000" not in proc.stdout


def _write_case_with_budget(
    arm: Path, pdb: str, rmsds, *, evals_actual: str
) -> None:
    d = arm / pdb
    d.mkdir(parents=True)
    fields = (
        ["pdb_id"]
        + [f"mode_rmsd_{i}" for i in range(10)]
        + ["n_poses", "restarts_finished", "evals_actual", "protocol_claim_eligible"]
    )
    row = {
        "pdb_id": pdb,
        "n_poses": "500",
        "restarts_finished": "10",
        "evals_actual": evals_actual,
        # Deliberately asserts compliance; the guard must NOT trust it.
        "protocol_claim_eligible": "1",
    }
    for i, v in enumerate(rmsds):
        row[f"mode_rmsd_{i}"] = "" if v is None else f"{v:.4f}"
    with (d / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerow(row)


def test_budget_guard_off_by_default(tmp_path: Path):
    """No --intended-evals → unchanged behaviour, empty evals_actual scores."""
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case_with_budget(arm, "1ABC", [5.0] * 10, evals_actual="")
    cases = mod.load_arm_dir(arm, strict=True)
    assert set(cases) == {"1ABC"}


def test_budget_guard_refuses_unwitnessed_spend(tmp_path: Path):
    """evals_actual empty is not evidence of compliance, whatever the flag says."""
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case_with_budget(arm, "1ABC", [5.0] * 10, evals_actual="")
    with pytest.raises(mod.BudgetWitnessError):
        mod.load_arm_dir(arm, strict=True, intended_evals=2_000_000)


def test_budget_guard_refuses_truncated_run(tmp_path: Path):
    """The 8.25%-of-budget case measured on arm A must not be scored."""
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case_with_budget(arm, "1ABC", [5.0] * 10, evals_actual="165000")
    with pytest.raises(mod.BudgetWitnessError) as exc:
        mod.load_arm_dir(arm, strict=True, intended_evals=2_000_000)
    assert "8.2%" in str(exc.value)


def test_budget_guard_accepts_full_spend(tmp_path: Path):
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case_with_budget(arm, "1ABC", [1.5] + [9.0] * 9, evals_actual="2000000")
    cases = mod.load_arm_dir(arm, strict=True, intended_evals=2_000_000)
    assert mod.s_top10(cases["1ABC"]) is True


def test_budget_fraction_is_tunable(tmp_path: Path):
    """The 90% line is a knob, not a constant baked into the metric."""
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    arm = tmp_path / "arm"
    _write_case_with_budget(arm, "1ABC", [5.0] * 10, evals_actual="1000000")
    with pytest.raises(mod.BudgetWitnessError):
        mod.load_arm_dir(arm, strict=True, intended_evals=2_000_000)
    cases = mod.load_arm_dir(
        arm, strict=True, intended_evals=2_000_000, min_budget_fraction=0.5
    )
    assert set(cases) == {"1ABC"}


def test_cli_budget_guard_fail_closed(tmp_path: Path):
    arm = tmp_path / "arm"
    _write_case_with_budget(arm, "1ABC", [5.0] * 10, evals_actual="165000")
    proc = subprocess.run(
        [
            sys.executable,
            str(BOOT),
            "--arm-dir",
            str(arm),
            "--bootstraps",
            "10",
            "--intended-evals",
            "2000000",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 2
    assert "FIXED search effort" in proc.stderr
    assert "0.0000" not in proc.stdout


def test_bootstrap_median_deterministic():
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    # 2/4 success → point 0.5
    success = [True, True, False, False]
    a = mod.bootstrap_median(success, n_boot=1000, seed=42)
    b = mod.bootstrap_median(success, n_boot=1000, seed=42)
    assert a == b
    assert a["point"] == pytest.approx(0.5)
    assert a["n_cases"] == 4
    assert a["n_success"] == 2
    assert a["median"] is not None
    assert 0.0 <= a["p05"] <= a["median"] <= a["p95"] <= 1.0


def test_compute_s_top10_stats_json_roundtrip(tmp_path: Path):
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    cases = {
        "AAAA": [1.0] + [3.0] * 9,  # success
        "BBBB": [3.0] * 10,  # fail
        "CCCC": [None, None, 1.9] + [None] * 7,  # success via rank 2
    }
    stats = mod.compute_s_top10_stats(cases, n_boot=500, seed=1, thresh=2.0)
    assert stats["n_cases"] == 3
    assert stats["n_success"] == 2
    assert stats["point"] == pytest.approx(2 / 3)
    assert stats["thresh_op"] == "<="
    assert stats["per_case_success"]["AAAA"] is True
    assert stats["per_case_success"]["BBBB"] is False
    assert stats["per_case_success"]["CCCC"] is True
    out = tmp_path / "stats.json"
    out.write_text(json.dumps(stats))
    loaded = json.loads(out.read_text())
    assert loaded["metric"] == "S_top10"


def test_cli_arm_dir_fail_closed(tmp_path: Path):
    d = tmp_path / "arm" / "1ABC"
    d.mkdir(parents=True)
    with (d / "result.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "rmsd_top1", "rmsd_bcr"])
        w.writeheader()
        w.writerow({"pdb_id": "1ABC", "rmsd_top1": "0.5", "rmsd_bcr": "0.4"})
    proc = subprocess.run(
        [sys.executable, str(BOOT), "--arm-dir", str(tmp_path / "arm")],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 2
    assert "mode_rmsd" in proc.stderr.lower() or "fail-closed" in proc.stderr.lower()


def test_cli_json_cases_success(tmp_path: Path):
    cases_path = tmp_path / "cases.json"
    cases_path.write_text(
        json.dumps(
            {
                "cases": {
                    "1ABC": [0.5, 3.0, 4.0],
                    "2DEF": [3.0, 3.0, 3.0],
                }
            }
        )
    )
    out = tmp_path / "out.json"
    proc = subprocess.run(
        [
            sys.executable,
            str(BOOT),
            "--cases",
            str(cases_path),
            "--bootstraps",
            "200",
            "--json-out",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "S_top10" in proc.stdout
    data = json.loads(out.read_text())
    assert data["n_cases"] == 2
    assert data["n_success"] == 1
    assert data["thresh_A"] == 2.0


def test_rank_table_emitted_order_not_rmsd_sorted(tmp_path: Path):
    mod = _load("bootstrap_3dsig_s_top10", BOOT)
    table = tmp_path / "ranks.csv"
    # rank 0 has high RMSD, rank 1 has low — order must follow rank, not RMSD
    table.write_text(
        "pdb_id,rank,rmsd\n"
        "1ABC,1,0.5\n"
        "1ABC,0,5.0\n"
        "1ABC,2,4.0\n"
    )
    cases = mod.load_rank_table(table)
    assert cases["1ABC"][0] == pytest.approx(5.0)
    assert cases["1ABC"][1] == pytest.approx(0.5)
    assert mod.s_top10(cases["1ABC"]) is True  # via rank 1


# ── parse_flexaid_arm_results ────────────────────────────────────────────────


def test_parse_filename_cf_and_restart():
    mod = _load("parse_flexaid_arm_results", PARSE)
    p = Path("1HNN_0.pdb")
    assert mod.parse_pose_filename(p, "1HNN") == (None, 0)
    p = Path("1HNN_3.pdb")
    assert mod.parse_pose_filename(p, "1HNN") == (None, 3)
    p = Path("1HNN_r2_0.pdb")
    assert mod.parse_pose_filename(p, "1HNN") == (2, 0)
    p = Path("1HNN_r2_7.pdb")
    assert mod.parse_pose_filename(p, "1HNN") == (2, 7)
    p = Path("1HNN_10_0.pdb")  # FO dual-suffix minPts=10 rank=0
    assert mod.parse_pose_filename(p, "1HNN") == (None, 0)
    p = Path("1HNN_r1_10_3.pdb")  # restart FO dual
    assert mod.parse_pose_filename(p, "1HNN") == (1, 3)
    assert mod.parse_pose_filename(Path("1HNN_INI.pdb"), "1HNN") is None


def test_mode_rmsds_emitted_order_single_emission():
    mod = _load("parse_flexaid_arm_results", PARSE)
    # Simulate heads list: ranks 0,1,2 with rmsds; not sorted by RMSD
    heads = []
    for rank, rmsd, cf in [(0, 5.0, -10.0), (1, 1.0, -9.0), (2, 3.0, -8.0)]:
        meta = {
            "rmsd_sym": rmsd,
            "rmsd_nosym": None,
            "cf": cf,
            "cf_app": None,
        }
        heads.append((None, rank, Path(f"X_{rank}.pdb"), meta))
    got = mod.mode_rmsds_emitted_order(heads, elected_restart=None, n=10)
    assert got[0] == 5.0
    assert got[1] == 1.0
    assert got[2] == 3.0
    assert got[3] is None
    assert mod.success_s_top10(got) == 1
    # S1 fails (mode 0 = 5.0)
    assert not (got[0] is not None and got[0] <= 2.0)


def test_mode_rmsds_multi_restart_uses_elected():
    mod = _load("parse_flexaid_arm_results", PARSE)
    heads = []
    # restart 0: rank0 CF=-5, rmsd=4.0; rank1 rmsd=0.5
    heads.append(
        (
            0,
            0,
            Path("X_r0_0.pdb"),
            {"rmsd_sym": 4.0, "rmsd_nosym": None, "cf": -5.0, "cf_app": None},
        )
    )
    heads.append(
        (
            0,
            1,
            Path("X_r0_1.pdb"),
            {"rmsd_sym": 0.5, "rmsd_nosym": None, "cf": -4.0, "cf_app": None},
        )
    )
    # restart 1: rank0 CF=-20 (better), rmsd=3.0; rank1 rmsd=2.5
    heads.append(
        (
            1,
            0,
            Path("X_r1_0.pdb"),
            {"rmsd_sym": 3.0, "rmsd_nosym": None, "cf": -20.0, "cf_app": None},
        )
    )
    heads.append(
        (
            1,
            1,
            Path("X_r1_1.pdb"),
            {"rmsd_sym": 2.5, "rmsd_nosym": None, "cf": -19.0, "cf_app": None},
        )
    )
    best = mod.elect_best_rank0(heads)
    assert best is not None
    assert best[1] == 1  # restart 1
    got = mod.mode_rmsds_emitted_order(heads, elected_restart=1, n=10)
    assert got[0] == 3.0
    assert got[1] == 2.5
    # Must not pick restart 0's rank1 (0.5) — that would be cross-restart mix
    assert 0.5 not in [v for v in got if v is not None]


def test_parse_cli_writes_mode_rmsd_columns(tmp_path: Path):
    out = tmp_path / "1HNN"
    out.mkdir()
    # ranks 0..2 single emission; rank 1 is native-like
    _write_pose(out / "1HNN_0.pdb", rmsd=4.5, cf=-100.0, cf_app=-99.0)
    _write_pose(out / "1HNN_1.pdb", rmsd=1.2, cf=-90.0, cf_app=-89.0)
    _write_pose(out / "1HNN_2.pdb", rmsd=3.0, cf=-80.0, cf_app=-79.0)
    # noise file ignored
    _write_pose(out / "1HNN_INI.pdb", rmsd=0.0, cf=0.0)

    proc = subprocess.run(
        [
            sys.executable,
            str(PARSE),
            "--arm",
            "A",
            "--pdb",
            "1HNN",
            "--out-dir",
            str(out),
            "--matrix-md5",
            "72d7c7396702331d96ff12d18f831796",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    rows = list(csv.DictReader((out / "result.csv").open()))
    assert len(rows) == 1
    row = rows[0]
    assert row["pdb_id"] == "1HNN"
    assert "mode_rmsd_0" in row
    assert float(row["mode_rmsd_0"]) == pytest.approx(4.5)
    assert float(row["mode_rmsd_1"]) == pytest.approx(1.2)
    assert float(row["mode_rmsd_2"]) == pytest.approx(3.0)
    assert row["mode_rmsd_3"] == ""
    assert row["success_s1"] == "0"  # top1 = 4.5
    assert row["success_s_top10"] == "1"  # mode 1 = 1.2
    assert float(row["rmsd_top1"]) == pytest.approx(4.5)
    assert float(row["rmsd_bcr"]) == pytest.approx(1.2)

    # Bootstrap must accept this arm dir
    boot = _load("bootstrap_3dsig_s_top10", BOOT)
    cases = boot.load_arm_dir(tmp_path, strict=True)
    assert "1HNN" in cases
    assert boot.s_top10(cases["1HNN"]) is True


def test_parse_multi_restart_cli(tmp_path: Path):
    out = tmp_path / "2GBP"
    out.mkdir()
    # r0 rank0 worse score, good rmsd on rank1 — must not leak into elected modes
    _write_pose(out / "2GBP_r0_0.pdb", rmsd=5.0, cf=-50.0)
    _write_pose(out / "2GBP_r0_1.pdb", rmsd=0.8, cf=-40.0)
    # r1 elected by better CF
    _write_pose(out / "2GBP_r1_0.pdb", rmsd=1.5, cf=-200.0)
    _write_pose(out / "2GBP_r1_1.pdb", rmsd=2.5, cf=-180.0)

    proc = subprocess.run(
        [
            sys.executable,
            str(PARSE),
            "--arm",
            "B",
            "--pdb",
            "2GBP",
            "--out-dir",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    row = list(csv.DictReader((out / "result.csv").open()))[0]
    assert float(row["mode_rmsd_0"]) == pytest.approx(1.5)
    assert float(row["mode_rmsd_1"]) == pytest.approx(2.5)
    assert row["success_s1"] == "1"
    assert row["success_s_top10"] == "1"
    # BCR can still see 0.8 from other restart
    assert float(row["rmsd_bcr"]) == pytest.approx(0.8)
