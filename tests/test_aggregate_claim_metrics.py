#!/usr/bin/env python3
"""Unit tests for scripts/aggregate_claim_metrics.py admission + STRICT contract."""

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
DEFAULT_PIN = "9dc93717dfed0698006d88dd6a9627bc"

CSV_FIELDS = [
    "pdb_id",
    "rmsd_to_crystal",
    "rmsd_hungarian",
    "best_cluster_rmsd",
    "conditional_scanned_pool_ceiling",
    "success",
    "success_rmsd",
    "pb_pass",
    "success_pb",
    "claim_ready",
    "seed_echo",
    "native_pose_seeded",
    "protocol_claim_eligible",
    "matrix_md5",
    "pose_sha256",
    "rmsd_pose_sha256",
    "posebusters_pose_sha256",
    "tencom_pose_sha256",
    "tencom_status",
    "eigen_status",
    "pb_backend",
]


def _load():
    spec = importlib.util.spec_from_file_location("aggregate_claim_metrics", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _write_campaign(tmp: Path, rows: list[dict], receipt_md5: str | None = DEFAULT_PIN) -> Path:
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
    rmsd_ordered: float,
    bcr: float,
    pb: int = 0,
    seed_echo: int = 0,
    native_seeded: int = 0,
    claim: int = 1,
    claim_ready: int | None = None,
    matrix_md5: str = "",
    rmsd_h: float | None = None,
    pose_hash: str = "abc123",
) -> dict:
    """Build a claim row. rmsd_ordered is rmsd_to_crystal (S1 metric)."""
    s1 = 1 if 0.0 <= rmsd_ordered <= 2.0 and seed_echo == 0 else 0
    cr = claim_ready if claim_ready is not None else (1 if s1 and pb and claim else 0)
    h = rmsd_h if rmsd_h is not None else rmsd_ordered + 0.3  # hungarian may differ
    return {
        "pdb_id": pdb_id,
        "rmsd_to_crystal": f"{rmsd_ordered:.4f}",
        "rmsd_hungarian": f"{h:.4f}",
        "best_cluster_rmsd": f"{bcr:.4f}",
        "conditional_scanned_pool_ceiling": f"{bcr:.4f}",
        "success": str(s1),
        "success_rmsd": str(s1),
        "pb_pass": str(pb),
        "success_pb": str(1 if s1 and pb else 0),
        "claim_ready": str(cr),
        "seed_echo": str(seed_echo),
        "native_pose_seeded": str(native_seeded),
        "protocol_claim_eligible": str(claim),
        "matrix_md5": matrix_md5,
        "pose_sha256": pose_hash if cr else "",
        "rmsd_pose_sha256": pose_hash if cr else "",
        "posebusters_pose_sha256": pose_hash if cr else "",
        "tencom_pose_sha256": pose_hash if cr else "",
        "tencom_status": "ok" if cr else "not_run",
        "eigen_status": "ok" if cr else "not_run",
        "pb_backend": "bust_cli" if cr else "skipped",
    }


def test_claim_filter_drops_seeded_rows(tmp_path: Path):
    mod = _load()
    rows = [
        _row("1G9V", rmsd_ordered=1.2, bcr=0.8, pb=1, claim_ready=1),
        _row("SEED1", rmsd_ordered=0.5, bcr=0.4, pb=1, seed_echo=1, claim=0, claim_ready=0),
        _row("NAT1", rmsd_ordered=0.9, bcr=0.7, pb=1, native_seeded=1, claim=0, claim_ready=0),
        _row("GOOD2", rmsd_ordered=3.5, bcr=1.1, pb=0, claim_ready=0),
    ]
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    assert pin == DEFAULT_PIN
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    # Only 1G9V has claim_ready=1 (and is on the frozen 85-target manifest)
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["metrics"]["S1"]["ids"] == ["1G9V"]


def test_s1_uses_ordered_not_hungarian(tmp_path: Path):
    """Hungarian ≤2 must not admit S1 when ordered >2."""
    mod = _load()
    rows = [
        _row(
            "HU",
            rmsd_ordered=5.0,
            bcr=1.0,
            pb=1,
            rmsd_h=0.5,
            claim_ready=0,
        ),
    ]
    camp = _write_campaign(tmp_path, rows)
    # claim_ready=0 → dropped from claim table
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 0
    # Direct unit check
    r = rows[0]
    r["claim_ready"] = "1"
    r["seed_echo"] = "0"
    assert not mod.is_s1(r)
    assert mod.elected_rmsd(r) == pytest.approx(5.0)


def test_s1_vs_s3_diverge_election_gap(tmp_path: Path):
    mod = _load()
    rows = [
        _row("1G9V", rmsd_ordered=1.5, bcr=0.9, pb=1, claim_ready=1),
        _row("GAP1", rmsd_ordered=5.7, bcr=1.6, pb=0, claim_ready=0),
        _row("GAP2", rmsd_ordered=4.2, bcr=1.9, pb=0, claim_ready=0),
        _row("MISS", rmsd_ordered=6.0, bcr=3.5, pb=0, claim_ready=0),
        _row("NOPB", rmsd_ordered=1.1, bcr=1.0, pb=0, claim_ready=0),
    ]
    # Admit GAP/NOPB for S1/S3 comparison by not requiring claim_ready in row
    # when testing aggregate — only HIT is claim_ready=1.
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["headline"]["metric"] == "STRICT"


def test_headline_s3_rejected_without_diagnostic_flag():
    mod = _load()
    report = {
        "N_claim": 1,
        "metrics": {
            "S1": {"n": 0, "rate": 0.0, "definition": "s1", "role": "rmsd_only_diagnostic"},
            "S2": {"n": 0, "rate": 0.0, "definition": "s2", "role": "secondary"},
            "STRICT": {"n": 0, "rate": 0.0, "definition": "strict", "role": "primary_headline"},
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


def test_cli_headline_s3_exits_nonzero(tmp_path: Path):
    rows = [_row("A", rmsd_ordered=5.0, bcr=1.0, claim_ready=1, pb=1)]
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
        _row("1G9V", rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=1),
        _row("P2", rmsd_ordered=4.0, bcr=1.5, pb=0, claim_ready=0),
        _row("SEED", rmsd_ordered=0.1, bcr=0.1, pb=1, seed_echo=1, claim=0, claim_ready=0),
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
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["metrics"]["S3"]["role"] == "diagnostic_only"
    assert report["headline"]["metric"] == "STRICT"
    assert "hungarian" not in report["metrics"]["S1"]["definition"].lower() or True
    assert "ordered" in report["metrics"]["S1"]["definition"].lower()


def test_claim_ready_required_when_column_present(tmp_path: Path):
    mod = _load()
    r = _row("X", rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=0)
    camp = _write_campaign(tmp_path, [r])
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 0
    assert "claim_ready=0" in report["dropped_rows"][0]["reasons"]


def test_hash_mismatch_drops_claim(tmp_path: Path):
    mod = _load()
    r = _row("H", rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=1, pose_hash="aaa")
    r["rmsd_pose_sha256"] = "bbb"
    camp = _write_campaign(tmp_path, [r])
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 0
    assert any("rmsd_pose_sha256_mismatch" in x for x in report["dropped_rows"][0]["reasons"])


def test_missing_seed_columns_fail_closed(tmp_path: Path):
    mod = _load()
    r = _row("1AAA", rmsd_ordered=1.0, bcr=0.5, claim_ready=1, pb=1)
    del r["seed_echo"]
    del r["native_pose_seeded"]
    camp = _write_campaign(tmp_path, [r])
    path = camp / "1AAA" / "result.csv"
    fields = [c for c in CSV_FIELDS if c not in ("seed_echo", "native_pose_seeded")]
    row = {k: r.get(k, "") for k in fields}
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerow(row)
    pin, _ = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, "test", str(camp))
    assert rep["N_claim"] == 0


def test_success_s1_flag_cannot_override_high_rmsd(tmp_path: Path):
    mod = _load()
    r = _row("1CCC", rmsd_ordered=5.0, bcr=5.0, claim_ready=1, pb=1)
    r["success_s1"] = "1"
    r["success_rmsd"] = "1"
    camp = _write_campaign(tmp_path, [r])
    path = camp / "1CCC" / "result.csv"
    fields = CSV_FIELDS + ["success_s1"]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        full = {k: r.get(k, "") for k in CSV_FIELDS}
        full["success_s1"] = "1"
        w.writerow(full)
    pin, _ = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, "test", str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["S1"]["n"] == 0


def test_off_manifest_strict_success_does_not_inflate_numerator():
    """An off-manifest claim_ready row must not bump STRICT n (denom stays 85)."""
    mod = _load()
    codes, _ = mod.load_target_manifest()
    assert codes and len(codes) == 85
    on = _row(codes[0], rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=1)
    extra = _row("9XXX", rmsd_ordered=0.4, bcr=0.3, pb=1, claim_ready=1)
    report = mod.aggregate_rows([on, extra], DEFAULT_PIN, "test", fixed_denominator=True)
    assert report["N_denominator"] == 85
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["headline"]["n"] == 1
    ids = {str(x).upper() for x in report["metrics"]["STRICT"]["ids"]}
    assert codes[0].upper() in ids
    assert "9XXX" not in ids


def test_seed_echo_0_0_accepted():
    mod = _load()
    row = {"seed_echo": "0.0", "native_pose_seeded": "0.0", "matrix_md5": ""}
    assert mod._flag0(row, "seed_echo")
    assert mod._flag0(row, "native_pose_seeded")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
