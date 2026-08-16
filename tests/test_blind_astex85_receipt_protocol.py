#!/usr/bin/env python3
"""Wave 4: receipt-gated blind Astex-85 protocol (no 85-target dock)."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "blind_astex85_receipt_protocol.py"
REPRO = ROOT / "scripts" / "reproduce_astex85.sh"
PIN = "72d7c7396702331d96ff12d18f831796"


def _load():
    spec = importlib.util.spec_from_file_location("blind_astex85_receipt_protocol", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _run(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )


def test_help_and_validate_defaults():
    help_p = _run(["--help"])
    assert help_p.returncode == 0
    assert "blind" in help_p.stdout.lower()
    val = _run(["validate-defaults"])
    assert val.returncode == 0
    assert "seed_elitism=0" in val.stdout
    assert PIN in val.stdout
    assert "%" not in val.stdout


def test_default_receipt_is_blind_not_oracle():
    mod = _load()
    rec = mod.default_blind_receipt(git_commit="abc", binary_path="/bin/FlexAIDdS")
    assert rec["n_targets"] == 85
    assert rec["native_pose_seeded"] == 0
    assert rec["seed_echo"] == 0
    assert rec["seed_elitism"] == 0
    assert rec["native_seed_frac"] == 0
    assert rec["matrix_md5"] == PIN
    assert mod.validate_blind_receipt(rec) == []
    try:
        mod.refuse_oracle_defaults("1", "0.90")
        raise AssertionError("oracle defaults must be refused")
    except mod.ProtocolError:
        pass


def test_claim_refuses_without_receipt(tmp_path: Path):
    proc = _run(["claim", "--dir", str(tmp_path)])
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "REFUSE" in combined
    assert "11.8%" not in combined
    assert "24.1%" not in combined
    assert "48.8%" not in combined


def test_claim_refuses_percent_without_csv(tmp_path: Path):
    proc = _run(
        [
            "emit",
            "--out",
            str(tmp_path),
            "--dry-run",
            "--git-commit",
            "deadbeef",
            "--binary-path",
            "/bin/FlexAIDdS",
        ]
    )
    assert proc.returncode == 0
    assert (tmp_path / "RUN_RECEIPT.json").is_file()
    assert "not launching an 85-target dock" in proc.stdout
    claim = _run(["claim", "--dir", str(tmp_path)])
    assert claim.returncode != 0
    combined = claim.stdout + claim.stderr
    assert "11.8%" not in combined
    assert "24.1%" not in combined
    assert "48.8%" not in combined
    assert "Not printing a success %" in combined


def test_claim_prints_s1_only_with_receipt_and_csv(tmp_path: Path):
    _run(
        [
            "emit",
            "--out",
            str(tmp_path),
            "--git-commit",
            "deadbeef",
            "--binary-path",
            "/bin/FlexAIDdS",
        ]
    )
    csv_path = tmp_path / "result.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["pdb_id", "rmsd_hungarian"])
        w.writeheader()
        for i in range(85):
            w.writerow({"pdb_id": f"{i:04d}", "rmsd_hungarian": "1.0" if i < 10 else "5.0"})
    claim = _run(["claim", "--dir", str(tmp_path)])
    assert claim.returncode == 0
    assert "10/85" in claim.stdout
    assert "11.8%" in claim.stdout


def test_oracle_emit_is_not_claimable(tmp_path: Path):
    emit = _run(
        [
            "emit",
            "--out",
            str(tmp_path),
            "--oracle-ceiling",
            "--git-commit",
            "deadbeef",
            "--binary-path",
            "/bin/FlexAIDdS",
        ]
    )
    assert emit.returncode == 0
    assert "ORACLE CEILING" in emit.stderr
    rec = json.loads((tmp_path / "RUN_RECEIPT.json").read_text(encoding="utf-8"))
    assert rec["seed_elitism"] == 1
    claim = _run(["claim", "--dir", str(tmp_path)])
    assert claim.returncode != 0
    combined = claim.stdout + claim.stderr
    assert "11.8%" not in combined
    assert "24.1%" not in combined
    assert "48.8%" not in combined


def test_reproduce_script_has_dry_run_and_stays_blind():
    text = REPRO.read_text(encoding="utf-8")
    assert "--dry-run" in text
    assert "SEED_ELITISM=0" in text
    assert "blind_astex85_receipt_protocol.py" in text
    assert "\nSEED_ELITISM=0\n" in text
    assert "\nSEED_ELITISM=1\n" not in text
