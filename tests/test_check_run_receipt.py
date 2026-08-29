#!/usr/bin/env python3
"""Tests for scripts/check_run_receipt.py scoring provenance checker."""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_run_receipt.py"
FIXTURES = Path(__file__).resolve().parent / "fixtures" / "receipts"
ENGINE_FIXTURE = FIXTURES / "engine_RUN_RECEIPT.json"
CAMPAIGN_FIXTURE = FIXTURES / "campaign_PREREGISTRATION.json"
PIN = "9dc93717dfed0698006d88dd6a9627bc"


def _load():
    spec = importlib.util.spec_from_file_location("check_run_receipt", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True
    )


# ── existing engine-dialect behaviour, unchanged ────────────────────────────


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
    proc = _run(str(tmp_path), "--require-matrix-9dc9")
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
    proc = _run(str(tmp_path))
    assert proc.returncode == 0
    assert "OK" in proc.stdout


# ── dialect detection ───────────────────────────────────────────────────────


def test_fixtures_exist():
    assert ENGINE_FIXTURE.is_file()
    assert CAMPAIGN_FIXTURE.is_file()


def test_detect_engine_dialect():
    mod = _load()
    data = json.loads(ENGINE_FIXTURE.read_text())
    assert mod.detect_dialect(data) == mod.DIALECT_ENGINE


def test_detect_campaign_dialect():
    mod = _load()
    data = json.loads(CAMPAIGN_FIXTURE.read_text())
    assert mod.detect_dialect(data) == mod.DIALECT_CAMPAIGN


def test_detect_legacy_provenance_is_engine_family():
    """provenance.json has no schema_version but is engine-written."""
    mod = _load()
    assert (
        mod.detect_dialect({"matrix_md5": PIN, "protocol_config": {}})
        == mod.DIALECT_ENGINE
    )


def test_detect_unknown_routes_to_engine_checks():
    """Bare dicts keep their historical behaviour — engine checks."""
    mod = _load()
    assert mod.detect_dialect({"matrix_md5": PIN}) == mod.DIALECT_UNKNOWN
    errs = mod.check_receipt({"matrix_md5": PIN, "binary_sha256": "x"})
    assert "missing_or_empty:git_commit" in errs


# ── the regression: campaign dialect must not fail the engine contract ──────


def test_campaign_dialect_does_not_demand_git_commit():
    mod = _load()
    data = json.loads(CAMPAIGN_FIXTURE.read_text())
    assert "git_commit" not in data
    errs = mod.check_receipt(data, require_matrix_9dc9=True)
    assert errs == [], f"campaign fixture should validate cleanly, got {errs}"


def test_campaign_dialect_cli_reports_dialect(tmp_path: Path):
    shutil.copy(CAMPAIGN_FIXTURE, tmp_path / "PREREGISTRATION.json")
    proc = _run(str(tmp_path), "--require-matrix-9dc9")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "campaign dialect" in proc.stdout
    assert "1jd0_ga_wal400_v1" in proc.stdout


def test_campaign_dialect_still_catches_its_own_missing_fields():
    mod = _load()
    data = json.loads(CAMPAIGN_FIXTURE.read_text())
    del data["frozen_utc"]
    errs = mod.check_receipt(data)
    assert "missing_or_empty:frozen_utc" in errs


def test_campaign_dialect_matrix_pin_enforced():
    mod = _load()
    data = json.loads(CAMPAIGN_FIXTURE.read_text())
    data["matrix_md5"] = "0" * 32
    errs = mod.check_receipt(data, require_matrix_9dc9=True)
    assert any(e.startswith("matrix_md5_not_9dc9") for e in errs)


# ── PREREGISTRATION.json resolution ─────────────────────────────────────────


def test_preregistration_resolves_from_directory(tmp_path: Path):
    """The rename must not turn a schema mismatch into FileNotFoundError."""
    shutil.copy(CAMPAIGN_FIXTURE, tmp_path / "PREREGISTRATION.json")
    proc = _run(str(tmp_path))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "FileNotFound" not in proc.stderr


def test_engine_receipt_wins_when_both_present(tmp_path: Path):
    shutil.copy(ENGINE_FIXTURE, tmp_path / "RUN_RECEIPT.json")
    shutil.copy(CAMPAIGN_FIXTURE, tmp_path / "PREREGISTRATION.json")
    mod = _load()
    assert mod.resolve_receipt_path(tmp_path).name == "RUN_RECEIPT.json"


def test_empty_dir_lists_every_accepted_name(tmp_path: Path):
    proc = _run(str(tmp_path))
    assert proc.returncode == 2
    for name in ("RUN_RECEIPT.json", "PREREGISTRATION.json", "provenance.json"):
        assert name in proc.stderr


# ── --require-engine-dialect ────────────────────────────────────────────────


def test_require_engine_dialect_rejects_campaign_clearly(tmp_path: Path):
    shutil.copy(CAMPAIGN_FIXTURE, tmp_path / "PREREGISTRATION.json")
    proc = _run(str(tmp_path), "--require-engine-dialect")
    assert proc.returncode == 1
    assert "campaign dialect, not engine dialect" in proc.stdout
    # exactly one message, not a pile of missing-key errors pointing the wrong way
    assert "git_commit" not in proc.stdout


def test_require_engine_dialect_accepts_engine(tmp_path: Path):
    shutil.copy(ENGINE_FIXTURE, tmp_path / "RUN_RECEIPT.json")
    proc = _run(str(tmp_path), "--require-engine-dialect", "--require-matrix-9dc9")
    assert proc.returncode == 0, proc.stdout + proc.stderr


# ── report ──────────────────────────────────────────────────────────────────


def test_json_out_records_dialect_and_source(tmp_path: Path):
    shutil.copy(CAMPAIGN_FIXTURE, tmp_path / "PREREGISTRATION.json")
    out = tmp_path / "report.json"
    proc = _run(str(tmp_path), "--json-out", str(out))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    rep = json.loads(out.read_text())
    assert rep["dialect"] == "campaign"
    assert rep["schema"] == "1jd0_ga_wal400_v1"
    assert rep["source"].endswith("PREREGISTRATION.json")
    assert rep["ok"] is True
