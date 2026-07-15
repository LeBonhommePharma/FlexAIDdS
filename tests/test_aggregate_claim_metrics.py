#!/usr/bin/env python3
"""Unit tests for scripts/aggregate_claim_metrics.py admission + S1/S2/S3 contract."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "aggregate_claim_metrics.py"
DEFAULT_PIN = "72d7c7396702331d96ff12d18f831796"

# Columns used by DatasetRunner-style result.csv claim rows
CSV_FIELDS = [
    "pdb_id",
    "rmsd_to_crystal",
    "rmsd_hungarian",
    "best_cluster_rmsd",
    "success",
    "success_rmsd",
    "pb_pass",
    "success_pb",
    "seed_echo",
    "native_pose_seeded",
    "protocol_claim_eligible",
    "matrix_md5",
]


def _load():
    spec = importlib.util.spec_from_file_location("aggregate_claim_metrics", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_campaign(tmp: Path, rows: list[dict], receipt_md5: str | None = DEFAULT_PIN) -> Path:
    """Write per-target result.csv tree + optional RUN_RECEIPT.json."""
    camp = tmp / "campaign"
    camp.mkdir(parents=True, exist_ok=True)
    for r in rows:
        pid = r["pdb_id"]
        d = camp / pid
        d.mkdir(exist_ok=True)
        path = d / "result.csv"
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
            w.writeheader()
            # fill missing keys
            full = {k: r.get(k, "") for k in CSV_FIELDS}
            w.writerow(full)
    if receipt_md5 is not None:
        (camp / "RUN_RECEIPT.json").write_text(
            json.dumps({"matrix_md5": receipt_md5, "run_id": "fixture"})
        )
    return camp


def _row(
    pdb_id: str,
    *,
    rmsd_h: float,
    bcr: float,
    pb: int = 0,
    seed_echo: int = 0,
    native_seeded: int = 0,
    claim: int = 1,
    matrix_md5: str = "",
    rmsd_c: float | None = None,
) -> dict:
    s1 = 1 if 0.0 <= rmsd_h <= 2.0 and seed_echo == 0 else 0
    return {
        "pdb_id": pdb_id,
        "rmsd_to_crystal": f"{(rmsd_c if rmsd_c is not None else rmsd_h + 0.5):.4f}",
        "rmsd_hungarian": f"{rmsd_h:.4f}",
        "best_cluster_rmsd": f"{bcr:.4f}",
        "success": str(s1),
        "success_rmsd": str(s1),
        "pb_pass": str(pb),
        "success_pb": str(1 if s1 and pb else 0),
        "seed_echo": str(seed_echo),
        "native_pose_seeded": str(native_seeded),
        "protocol_claim_eligible": str(claim),
        "matrix_md5": matrix_md5,
    }


# ── admission filter ─────────────────────────────────────────────────────────


def test_claim_filter_drops_seeded_rows(tmp_path: Path):
    mod = _load()
    rows = [
        _row("GOOD1", rmsd_h=1.2, bcr=0.8, pb=1),
        _row("SEED1", rmsd_h=0.5, bcr=0.4, pb=1, seed_echo=1, claim=0),
        _row("NAT1", rmsd_h=0.9, bcr=0.7, pb=1, native_seeded=1, claim=0),
        _row("GOOD2", rmsd_h=3.5, bcr=1.1, pb=0),  # S1 fail, S3 hit
    ]
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    assert pin == DEFAULT_PIN
    assert src == "RUN_RECEIPT.json"
    loaded = mod.load_campaign_rows(camp)
    report = mod.aggregate_rows(loaded, pin, src, str(camp))
    assert report["N_raw"] == 4
    assert report["N_claim"] == 2
    assert report["N_dropped"] == 2
    dropped_ids = {d["pdb_id"] for d in report["dropped_rows"]}
    assert dropped_ids == {"SEED1", "NAT1"}
    assert set(report["metrics"]["S1"]["ids"]) == {"GOOD1"}
    assert set(report["metrics"]["S3"]["ids"]) == {"GOOD1", "GOOD2"}


def test_matrix_md5_pin_filters_wrong_matrix(tmp_path: Path):
    mod = _load()
    rows = [
        _row("OK", rmsd_h=1.0, bcr=0.9, matrix_md5=DEFAULT_PIN),
        _row(
            "BADMAT",
            rmsd_h=1.0,
            bcr=0.9,
            matrix_md5="deadbeefdeadbeefdeadbeefdeadbeef",
        ),
    ]
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 1
    assert report["dropped_rows"][0]["pdb_id"] == "BADMAT"
    assert any("matrix_md5_mismatch" in r for r in report["dropped_rows"][0]["reasons"])


def test_s1_vs_s3_diverge_election_gap(tmp_path: Path):
    """S3 can exceed S1 when BCR finds a near-native the elector missed."""
    mod = _load()
    rows = [
        # S1+S2+S3
        _row("HIT", rmsd_h=1.5, bcr=0.9, pb=1),
        # election gap: BCR ≤2, elected >2
        _row("GAP1", rmsd_h=5.7, bcr=1.6, pb=0),
        _row("GAP2", rmsd_h=4.2, bcr=1.9, pb=0),
        # both fail
        _row("MISS", rmsd_h=6.0, bcr=3.5, pb=0),
        # S1 without PB → S2 fail
        _row("NOPB", rmsd_h=1.1, bcr=1.0, pb=0),
    ]
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 5
    assert report["metrics"]["S1"]["n"] == 2  # HIT, NOPB
    assert report["metrics"]["S2"]["n"] == 1  # HIT only
    assert report["metrics"]["S3"]["n"] == 4  # HIT, GAP1, GAP2, NOPB
    assert report["election_gap"]["n"] == 2
    assert set(report["election_gap"]["ids"]) == {"GAP1", "GAP2"}
    # Rates diverge
    assert report["metrics"]["S1"]["rate"] == pytest.approx(0.4)
    assert report["metrics"]["S3"]["rate"] == pytest.approx(0.8)
    assert report["metrics"]["S3"]["role"] == "diagnostic_only"
    assert report["headline"]["metric"] == "S1"


def test_headline_s3_rejected_without_diagnostic_flag():
    mod = _load()
    report = {
        "N_claim": 1,
        "metrics": {
            "S1": {"n": 0, "rate": 0.0, "definition": "s1", "role": "primary_headline"},
            "S2": {"n": 0, "rate": 0.0, "definition": "s2", "role": "secondary"},
            "S3": {
                "n": 1,
                "rate": 1.0,
                "definition": "s3",
                "role": "diagnostic_only",
            },
        },
    }
    _, err = mod.apply_headline(report, "s3", diagnostic_only=False)
    assert err == 2
    out, err2 = mod.apply_headline(report, "s3", diagnostic_only=True)
    assert err2 is None
    assert out["headline"]["metric"] == "S3"
    assert out["headline"].get("warning")


def test_cli_headline_s3_exits_nonzero(tmp_path: Path):
    rows = [_row("A", rmsd_h=5.0, bcr=1.0)]
    camp = _write_campaign(tmp_path, rows)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(camp), "--headline", "s3"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 2
    assert "CONTRACT VIOLATION" in proc.stderr


def test_cli_happy_path_json(tmp_path: Path):
    rows = [
        _row("P1", rmsd_h=1.0, bcr=0.5, pb=1),
        _row("P2", rmsd_h=4.0, bcr=1.5, pb=0),
        _row("SEED", rmsd_h=0.1, bcr=0.1, pb=1, seed_echo=1, claim=0),
    ]
    camp = _write_campaign(tmp_path, rows)
    out_json = tmp_path / "out.json"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), str(camp), "--json", str(out_json), "--quiet"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    report = json.loads(out_json.read_text())
    assert report["N_claim"] == 2
    assert report["metrics"]["S1"]["n"] == 1
    assert report["metrics"]["S3"]["n"] == 2
    assert report["metrics"]["S3"]["role"] == "diagnostic_only"
    assert report["matrix_md5_pin"] == DEFAULT_PIN


def test_flat_summary_csv(tmp_path: Path):
    mod = _load()
    rows = [
        _row("X", rmsd_h=1.8, bcr=1.2, pb=1),
        _row("Y", rmsd_h=2.5, bcr=2.1, pb=0),
    ]
    csv_path = tmp_path / "summary.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in CSV_FIELDS})
    # No receipt → default pin; rows omit matrix_md5 → campaign pin applies
    report = mod.aggregate_rows(
        mod.load_rows_from_csv(csv_path),
        DEFAULT_PIN,
        "default_pin",
        str(tmp_path),
    )
    assert report["N_claim"] == 2
    assert report["metrics"]["S1"]["n"] == 1
    assert report["metrics"]["S3"]["n"] == 1


def test_cli_csv_flag(tmp_path: Path):
    rows = [_row("Z", rmsd_h=0.5, bcr=0.4, pb=1)]
    csv_path = tmp_path / "astex_diverse_results.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerow({k: rows[0].get(k, "") for k in CSV_FIELDS})
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--csv", str(csv_path), "--quiet"],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    report = json.loads(proc.stdout)
    assert report["N_claim"] == 1
    assert report["metrics"]["S1"]["n"] == 1
    assert report["metrics"]["S2"]["n"] == 1


if __name__ == "__main__":
    # Allow bare `python3 tests/test_aggregate_claim_metrics.py`
    pytest.main([__file__, "-v"])


def test_missing_seed_columns_fail_closed(tmp_path: Path):
    """Fail-closed: omit seed columns → not claim-eligible."""
    mod = _load()
    r = _row("1AAA", rmsd_h=1.0, bcr=0.5)
    del r["seed_echo"]
    del r["native_pose_seeded"]
    camp = _write_campaign(tmp_path, [r])
    # rewrite without seed columns
    import csv
    path = camp / "1AAA" / "result.csv"
    fields = [c for c in CSV_FIELDS if c not in ("seed_echo", "native_pose_seeded")]
    row = {k: r.get(k, "") for k in fields}
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(row)
    rows = mod.load_campaign_rows(camp)
    pin, _ = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(rows, pin, "test", str(camp))
    assert rep["N_claim"] == 0
    assert rep["N_dropped"] == 1


def test_rmsd_top1_three_engine_schema(tmp_path: Path):
    mod = _load()
    camp = tmp_path / "camp"
    camp.mkdir()
    d = camp / "1BBB"
    d.mkdir()
    fields = ["pdb_id", "rmsd_top1", "rmsd_bcr", "seed_echo", "native_pose_seeded", "matrix_md5"]
    with (d / "result.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow({
            "pdb_id": "1BBB",
            "rmsd_top1": "1.2",
            "rmsd_bcr": "0.8",
            "seed_echo": "0",
            "native_pose_seeded": "0",
            "matrix_md5": "",
        })
    (camp / "RUN_RECEIPT.json").write_text('{"matrix_md5": "%s"}' % DEFAULT_PIN)
    rows = mod.load_campaign_rows(camp)
    pin, _ = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(rows, pin, "test", str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["S1"]["n"] == 1
    assert rep["metrics"]["S3"]["n"] == 1


def test_success_s1_flag_cannot_override_high_rmsd(tmp_path: Path):
    mod = _load()
    r = _row("1CCC", rmsd_h=5.0, bcr=5.0)
    r["success_s1"] = "1"
    r["success_rmsd"] = "1"
    camp = _write_campaign(tmp_path, [r])
    rows = mod.load_campaign_rows(camp)
    # inject success_s1 into loaded row via rewriting
    import csv
    path = camp / "1CCC" / "result.csv"
    fields = CSV_FIELDS + ["success_s1"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        full = {k: r.get(k, "") for k in CSV_FIELDS}
        full["success_s1"] = "1"
        w.writerow(full)
    rows = mod.load_campaign_rows(camp)
    pin, _ = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(rows, pin, "test", str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["S1"]["n"] == 0


def test_seed_echo_0_0_accepted():
    mod = _load()
    row = {"seed_echo": "0.0", "native_pose_seeded": "0.0", "matrix_md5": ""}
    assert mod._flag0(row, "seed_echo")
    assert mod._flag0(row, "native_pose_seeded")
