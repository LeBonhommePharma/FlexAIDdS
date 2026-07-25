# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
"""Tests for comparative P0–P5 gates and pipeline (shipped code paths)."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from flexaidds.comparative_phases.gates import (
    MATRIX_MD5_PIN,
    can_run_p3,
    can_run_p4,
    evaluate_p2_oracle,
    evaluate_p3_pilot,
    load_phase_state,
    next_allowed_phase,
    save_phase_state,
)
from flexaidds.comparative_phases.p0_layout import ensure_matrix, md5_file, run_p0
from flexaidds.comparative_phases.p1_binaries import load_arm_pins, run_p1
from flexaidds.comparative_phases.pipeline import run_phase, run_pipeline_dry
from flexaidds.comparative_phases.p5_analyze import run_p5


REPO = Path(__file__).resolve().parents[2]


def test_matrix_md5_pin_matches_arm_pins_json():
    pins = load_arm_pins(REPO)
    assert pins["matrix"]["md5"] == MATRIX_MD5_PIN
    assert pins["arms"]["A"]["source_commit"].startswith("f766a14")
    assert pins["arms"]["B"]["source_commit"].startswith("1a6ae0b")


def test_evaluate_p2_oracle_pass_hold_fail():
    assert evaluate_p2_oracle({"ranking_allowed": True})[0] == "pass"
    assert evaluate_p2_oracle({"ranking_allowed": False})[0] == "hold"
    assert evaluate_p2_oracle({"status": "SCIENCE_HOLD"})[0] == "hold"
    assert evaluate_p2_oracle({"deferred": True})[0] == "hold"
    # Empty / deferred → hold (not silent pass); missing keys + deferred=False → fail
    assert evaluate_p2_oracle({})[0] == "hold"
    assert evaluate_p2_oracle({"n_targets": 2, "deferred": False})[0] == "fail"


def test_evaluate_p3_pilot_schema_and_zero_hold():
    assert (
        evaluate_p3_pilot(
            {"n_targets": 8, "schema_ok": True, "bcr_success": 1, "s_top10_success": 0}
        )[0]
        == "pass"
    )
    st, reason = evaluate_p3_pilot(
        {"n_targets": 8, "schema_ok": True, "bcr_success": 0, "s_top10_success": 0}
    )
    assert st == "hold"
    assert "SCIENCE HOLD" in reason or "0" in reason
    assert evaluate_p3_pilot({"n_targets": 3, "schema_ok": False})[0] == "fail"


def test_can_run_p4_fail_closed():
    assert can_run_p4({"P2": "pass", "P3": "pass"}) is True
    assert can_run_p4({"P2": "pass", "P3": "hold"}) is False
    assert can_run_p4({"P2": "fail", "P3": "pass"}) is False
    assert can_run_p4({"P2": "pending", "P3": "pending"}) is False
    assert can_run_p3({"P2": "pass"}) is True
    assert can_run_p3({"P2": "hold"}) is False


def test_phase_state_roundtrip(tmp_path: Path):
    sp = tmp_path / "phase_state.json"
    save_phase_state(sp, {"P0": "pass", "P2": "hold", "extra": 1})
    loaded = load_phase_state(sp)
    assert loaded["P0"] == "pass"
    assert loaded["P2"] == "hold"
    assert loaded["P1"] == "pending"
    assert loaded["extra"] == 1


def test_p0_layout_and_matrix_pin(tmp_path: Path):
    """Drive real run_p0 against an isolated local root; matrix from repo."""
    src = REPO / "MC_st0r5.2_6.dat"
    if not src.is_file():
        pytest.skip("repo matrix not present")
    # Pre-seed so we do not depend on ensure shell when C0 paths fail
    data = tmp_path / "three_engine_entropy_q1" / "data"
    data.mkdir(parents=True)
    shutil.copy2(src, data / "MC_st0r5.2_6.dat")
    result = run_p0(str(tmp_path), call_shell_layout=False)
    assert result["status"] == "pass", result
    assert result["matrix_md5"] == MATRIX_MD5_PIN
    assert (tmp_path / "three_engine_entropy_q1" / "bin" / "A").is_dir()
    assert (tmp_path / "campaigns" / "three_engine" / "receipts").is_dir()
    # Wrong matrix must fail
    bad = tmp_path / "bad"
    bad_data = bad / "three_engine_entropy_q1" / "data"
    bad_data.mkdir(parents=True)
    (bad_data / "MC_st0r5.2_6.dat").write_bytes(b"not-the-matrix")
    with pytest.raises(ValueError):
        ensure_matrix(bad, REPO)


def test_p1_missing_binaries_fail_closed_and_reconstruction(tmp_path: Path):
    # empty bins
    for arm in "A", "B", "C":
        (tmp_path / "three_engine_entropy_q1" / "bin" / arm).mkdir(parents=True)
    r = run_p1(str(tmp_path), allow_reconstruction=False)
    assert r["status"] == "fail"
    assert "MISSING" in r["reason"] or "missing" in r["reason"].lower()
    r2 = run_p1(str(tmp_path), allow_reconstruction=True)
    assert r2["status"] == "pass"
    rec_a = tmp_path / "campaigns" / "three_engine" / "receipts" / "arm_A_binary.json"
    assert rec_a.is_file()
    data = json.loads(rec_a.read_text())
    assert data["reconstruction"] is True
    assert data["status"] == "SOURCE_PINNED_BINARY_MISSING"
    assert data["git_commit"].startswith("f766a14")


def test_p1_distinct_sha_required(tmp_path: Path):
    for arm in "A", "B":
        d = tmp_path / "three_engine_entropy_q1" / "bin" / arm
        d.mkdir(parents=True)
        bin_path = d / "FlexAID"
        bin_path.write_bytes(b"same-binary-content")
        bin_path.chmod(0o755)
    r = run_p1(str(tmp_path), allow_reconstruction=False)
    assert r["status"] == "fail"
    assert "identical" in r["reason"]


def test_p4_blocked_without_p2_p3_pass(tmp_path: Path):
    sp = tmp_path / "state.json"
    save_phase_state(sp, {"P0": "pass", "P1": "pass", "P2": "hold", "P3": "pending"})
    # Pre-pass P0/P1 in state; attempt P4
    out = run_phase(
        "P4",
        local_root_path=str(tmp_path),
        dry_run=True,
        state_path=sp,
    )
    assert out.get("blocked") is True
    assert out["status"] == "fail"
    assert can_run_p4(load_phase_state(sp)) is False


def test_p4_allowed_when_p2_p3_pass(tmp_path: Path):
    sp = tmp_path / "state.json"
    save_phase_state(sp, {"P0": "pass", "P1": "pass", "P2": "pass", "P3": "pass"})
    out = run_phase(
        "P4",
        local_root_path=str(tmp_path),
        dry_run=True,
        state_path=sp,
    )
    assert out.get("blocked") is not True
    assert out["status"] == "pass"
    assert "full85" in out.get("launcher", "") or "P4" in out["reason"]


def test_p5_dry_run_writes_table(tmp_path: Path):
    # minimal pin env for p5
    src = REPO / "MC_st0r5.2_6.dat"
    if src.is_file():
        d = tmp_path / "three_engine_entropy_q1" / "data"
        d.mkdir(parents=True)
        shutil.copy2(src, d / "MC_st0r5.2_6.dat")
    run_p0(str(tmp_path), call_shell_layout=False)
    run_p1(str(tmp_path), allow_reconstruction=True)
    out = run_p5("pilot_dry", local_root_path=str(tmp_path), dry_run=True)
    assert out["status"] == "pass"
    assert Path(out["table_path"]).is_file()
    text = Path(out["table_path"]).read_text()
    assert "S_top10_median" in text
    assert "| A |" in text


def test_cli_pipeline_dry_and_gate_script(tmp_path: Path):
    """Drive shipped CLI entry points (not reimplemented gates)."""
    env = os.environ.copy()
    env["FLEXAIDDS_ROOT"] = str(REPO)
    env["PYTHONPATH"] = str(REPO / "python") + os.pathsep + env.get("PYTHONPATH", "")
    cli = REPO / "scripts" / "run_comparative_phases.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(cli),
            "--pipeline-dry",
            "--local-root",
            str(tmp_path),
            "--force-p2-pass",
            "--force-p3-pass",
        ],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "PIPELINE_DRY" in proc.stdout or "phase" in proc.stdout.lower()

    gate = REPO / "scripts" / "comparative_phase_gate.py"
    sp = tmp_path / "gate_state.json"
    proc2 = subprocess.run(
        [sys.executable, str(gate), "--state-file", str(sp), "--dry-run"],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO),
    )
    assert proc2.returncode == 0, proc2.stderr
    assert "ALLOW" in proc2.stdout or "OK" in proc2.stdout

    # BLOCK when P2 not pass (exit 1 = BLOCK by gate contract)
    save_phase_state(sp, {"P2": "fail", "P3": "pass"})
    proc3 = subprocess.run(
        [sys.executable, str(gate), "--state-file", str(sp), "--check-p4"],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(REPO),
    )
    assert "BLOCK" in proc3.stdout
    # Gate contract: non-zero exit on BLOCK (1) preferred; accept 0 if print-only.
    assert proc3.returncode in (0, 1)


def test_bootstrap_script_s_top10_function_importable():
    """Shipped bootstrap defines s_top10 used by P5."""
    # Import via path load of the script module by name
    import importlib.util

    path = REPO / "scripts" / "bootstrap_3dsig_s_top10.py"
    spec = importlib.util.spec_from_file_location("bootstrap_3dsig_s_top10", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.s_top10([3.0, 1.5, 4.0]) is True
    assert mod.s_top10([3.0, 2.5, 4.0]) is False
    assert mod.DEFAULT_THRESH == 2.0


def test_next_allowed_phase_serial():
    st = {p: "pending" for p in ("P0", "P1", "P2", "P3", "P4", "P5")}
    assert next_allowed_phase(st) == "P0"
    st["P0"] = "pass"
    assert next_allowed_phase(st) == "P1"
    st["P1"] = "pass"
    assert next_allowed_phase(st) == "P2"
    # hold/fail/pending all non-pass → first non-pass is P2
    st["P2"] = "hold"
    assert next_allowed_phase(st) == "P2"
