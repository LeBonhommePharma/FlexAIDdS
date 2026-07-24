#!/usr/bin/env python3
"""Tests for scripts/check_run_receipt.py scoring provenance checker."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_run_receipt.py"
PIN = "9dc93717dfed0698006d88dd6a9627bc"


def _load():
    spec = importlib.util.spec_from_file_location("check_run_receipt", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_missing_git_commit_fails():
    mod = _load()
    errs = mod.check_receipt(
        {"matrix_md5": PIN, "binary_sha256": "abc", "git_commit": ""},
        require_matrix_9dc9=True,
    )
    assert "missing_or_empty:git_commit" in errs


def test_ok_receipt():
    mod = _load()
    errs = mod.check_receipt(
        {
            "matrix_md5": PIN,
            "binary_sha256": "deadbeef",
            "git_commit": "6b81995b" + "0" * 32,
        },
        require_matrix_9dc9=True,
    )
    assert errs == []


def test_cli_fail_on_empty_commit(tmp_path: Path):
    rec = tmp_path / "RUN_RECEIPT.json"
    rec.write_text(
        json.dumps(
            {
                "matrix_md5": PIN,
                "binary_path": "/x/FlexAIDdS",
                "git_commit": "",
            }
        )
    )
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path), "--require-matrix-9dc9"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 1
    assert "git_commit" in proc.stdout


def test_cli_ok(tmp_path: Path):
    rec = tmp_path / "RUN_RECEIPT.json"
    rec.write_text(
        json.dumps(
            {
                "matrix_md5": PIN,
                "binary_sha256": "aa" * 32,
                "git_commit": "abcdef1234567890",
            }
        )
    )
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert "OK" in proc.stdout
