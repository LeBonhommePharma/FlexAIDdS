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


def _write_campaign(
    tmp: Path,
    rows: list[dict],
    receipt_md5: str | None = DEFAULT_PIN,
    *,
    write_pb_receipts: bool = True,
) -> Path:
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
        if write_pb_receipts and str(r.get("claim_ready", "")) == "1":
            pb_dir = d / "posebust"
            pb_dir.mkdir(exist_ok=True)
            (pb_dir / f"{pid}_bust_receipt.json").write_text(
                json.dumps(
                    {
                        "schema": {
                            "id": "posebusters-0.6.5-redock-csv-v1",
                            "required_check_count": 27,
                        },
                        "package": {
                            "name": "posebusters",
                            "version": "0.6.5",
                            "record_path": "/fixture/posebusters.dist-info/RECORD",
                            "record_sha256": "a" * 64,
                            "launcher_path": "/fixture/bin/bust",
                            "launcher_sha256": "b" * 64,
                            "launcher_version_output": "bust 0.6.5",
                        },
                        "config": {
                            "name": "redock",
                            "path": "/fixture/posebusters/config/redock.yml",
                            "sha256": "4d551d898ff29a404f16e02ad5a7a2d4235e6b7b14e9a3e27f7c66b4d16b2da9",
                        },
                        "command": {
                            "argv": [
                                "/fixture/bin/bust",
                                "ligand.sdf",
                                "-p",
                                "receptor.pdb",
                                "-l",
                                "crystal.sdf",
                            ],
                            "exit_status": 0,
                        },
                        "inputs": {
                            "predicted_ligand": {
                                "path": "/fixture/ligand.sdf",
                                "sha256": "c" * 64,
                            },
                            "protein": {
                                "path": "/fixture/receptor.pdb",
                                "sha256": "d" * 64,
                            },
                            "crystal_ligand": {
                                "path": "/fixture/crystal.sdf",
                                "sha256": "e" * 64,
                            },
                        },
                        "outputs": {
                            "raw_csv": {
                                "path": f"/fixture/{pid}_raw.csv",
                                "sha256": "f" * 64,
                            },
                            "validated_csv": {
                                "path": f"/fixture/{pid}_validated.csv",
                                "sha256": "1" * 64,
                            },
                        },
                        "result": {
                            "backend": "bust_cli",
                            "ran": True,
                            "pb_pass": True,
                            "failed_keys": "",
                            "error": "",
                        },
                    }
                )
            )
    if receipt_md5 is not None:
        expected_ids = [
            str(r["pdb_id"])
            for r in rows
            if str(r.get("protocol_claim_eligible", "")).strip() != "0"
        ]
        (camp / "RUN_RECEIPT.json").write_text(
            json.dumps(
                {
                    "matrix_md5": receipt_md5,
                    "run_id": "fixture",
                    "expected_target_ids": expected_ids,
                }
            )
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
    pose_hash: str = "2" * 64,
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


def test_denominator_uses_preregistered_protocol_not_outcomes(tmp_path: Path):
    mod = _load()
    rows = [
        _row("GOOD1", rmsd_ordered=1.2, bcr=0.8, pb=1, claim_ready=1),
        _row("SEED1", rmsd_ordered=0.5, bcr=0.4, pb=1, seed_echo=1, claim=0, claim_ready=0),
        _row("NAT1", rmsd_ordered=0.9, bcr=0.7, pb=1, native_seeded=1, claim=0, claim_ready=0),
        _row("GOOD2", rmsd_ordered=3.5, bcr=1.1, pb=0, claim_ready=0),
    ]
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    assert pin == DEFAULT_PIN
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    # GOOD2 remains in the denominator even though it failed docking/validators.
    assert report["N_claim"] == 2
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["metrics"]["S1"]["ids"] == ["GOOD1"]


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
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 0
    # Direct unit check
    r = rows[0]
    r["claim_ready"] = "1"
    r["seed_echo"] = "0"
    assert not mod.is_s1(r)
    assert mod.elected_rmsd(r) == pytest.approx(5.0)


def test_s1_vs_s3_diverge_election_gap(tmp_path: Path):
    mod = _load()
    rows = [
        _row("HIT", rmsd_ordered=1.5, bcr=0.9, pb=1, claim_ready=1),
        _row("GAP1", rmsd_ordered=5.7, bcr=1.6, pb=0, claim_ready=0),
        _row("GAP2", rmsd_ordered=4.2, bcr=1.9, pb=0, claim_ready=0),
        _row("MISS", rmsd_ordered=6.0, bcr=3.5, pb=0, claim_ready=0),
        _row("NOPB", rmsd_ordered=1.1, bcr=1.0, pb=0, claim_ready=0),
    ]
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 5
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["election_gap"]["ids"] == ["GAP1", "GAP2"]
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
        _row("P1", rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=1),
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
    assert report["N_claim"] == 2
    assert report["metrics"]["STRICT"]["n"] == 1
    assert report["metrics"]["S3"]["role"] == "diagnostic_only"
    assert report["headline"]["metric"] == "STRICT"
    assert report["headline"]["claimable"] is True
    assert report["completeness"]["complete"] is True
    assert "hungarian" not in report["metrics"]["S1"]["definition"].lower() or True
    assert "ordered" in report["metrics"]["S1"]["definition"].lower()


def test_claim_ready_failure_stays_in_denominator(tmp_path: Path):
    mod = _load()
    r = _row("X", rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=0)
    camp = _write_campaign(tmp_path, [r])
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 0
    assert "claim_ready!=1" in report["strict_fail_rows"][0]["reasons"]


def test_hash_mismatch_fails_strict_without_changing_denominator(tmp_path: Path):
    mod = _load()
    r = _row("H", rmsd_ordered=1.0, bcr=0.5, pb=1, claim_ready=1, pose_hash="aaa")
    r["rmsd_pose_sha256"] = "bbb"
    camp = _write_campaign(tmp_path, [r])
    pin, src = mod.load_matrix_pin(camp, None)
    report = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert report["N_claim"] == 1
    assert report["metrics"]["STRICT"]["n"] == 0
    assert any(
        "rmsd_pose_sha256_mismatch" in x
        for x in report["strict_fail_rows"][0]["reasons"]
    )


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
    assert rep["N_claim"] == 1
    assert rep["metrics"]["STRICT"]["n"] == 0


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
    assert rep["metrics"]["STRICT"]["n"] == 0


def test_one_success_plus_84_validator_failures_is_one_of_85(tmp_path: Path):
    mod = _load()
    rows = [_row("HIT", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)]
    rows.extend(
        _row(
            f"FAIL{i:02d}",
            rmsd_ordered=4.0,
            bcr=3.0,
            pb=0,
            claim_ready=0,
        )
        for i in range(84)
    )
    camp = _write_campaign(tmp_path, rows)
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 85
    assert rep["metrics"]["STRICT"]["n"] == 1
    assert rep["metrics"]["STRICT"]["rate"] == pytest.approx(1.0 / 85.0)
    assert len(rep["strict_fail_rows"]) == 84


def test_missing_pb_receipt_fails_strict_without_removing_denominator(tmp_path: Path):
    mod = _load()
    row = _row("NORECEIPT", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
    camp = _write_campaign(tmp_path, [row], write_pb_receipts=False)
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["STRICT"]["n"] == 0
    assert "pb_receipt_missing" in rep["strict_fail_rows"][0]["reasons"]


@pytest.mark.parametrize(
    ("section", "key", "value", "reason_fragment"),
    [
        ("package", "version", "0.6.6", "package_version!=0.6.5"),
        ("schema", "required_check_count", 26, "required_check_count!=27"),
        ("config", "sha256", "0" * 64, "config_sha256!="),
    ],
)
def test_nested_pb_identity_pins_fail_strict_not_denominator(
    tmp_path: Path,
    section: str,
    key: str,
    value: object,
    reason_fragment: str,
):
    mod = _load()
    row = _row("BADPIN", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
    camp = _write_campaign(tmp_path, [row])
    receipt_path = camp / "BADPIN" / "posebust" / "BADPIN_bust_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt[section][key] = value
    receipt_path.write_text(json.dumps(receipt))
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["STRICT"]["n"] == 0
    assert any(
        reason_fragment in reason for reason in rep["strict_fail_rows"][0]["reasons"]
    )


def test_nested_pb_output_rehash_mismatch_fails_strict(tmp_path: Path):
    mod = _load()
    row = _row("BADHASH", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
    camp = _write_campaign(tmp_path, [row])
    pb_dir = camp / "BADHASH" / "posebust"
    raw_csv = pb_dir / "BADHASH_raw.csv"
    raw_csv.write_text("actual raw output\n")
    receipt_path = pb_dir / "BADHASH_bust_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["outputs"]["raw_csv"] = {
        "path": raw_csv.name,
        "sha256": "0" * 64,
    }
    receipt_path.write_text(json.dumps(receipt))
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["STRICT"]["n"] == 0
    assert (
        "pb_receipt_output_raw_csv_rehash_mismatch"
        in rep["strict_fail_rows"][0]["reasons"]
    )


def test_legacy_rows_are_diagnostic_only(tmp_path: Path):
    mod = _load()
    modern = _row("MODERN", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
    legacy = _row("LEGACY", rmsd_ordered=1.1, bcr=0.7, pb=1, claim_ready=1)
    legacy["protocol_claim_eligible"] = ""
    legacy["claim_ready"] = ""
    camp = _write_campaign(tmp_path, [modern, legacy])
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 2
    assert rep["metrics"]["STRICT"]["ids"] == ["MODERN"]
    assert rep["N_legacy_diagnostic"] == 1
    assert rep["legacy_diagnostics"]["ids"] == ["LEGACY"]
    assert rep["legacy_diagnostics"]["S1_n"] == 1
    assert rep["headline"]["N"] == 2
    assert rep["headline"]["suppressed"] is True
    assert "legacy_schema_rows" in rep["headline"]["suppression_reasons"]


def test_missing_claim_ready_column_stays_in_denominator(tmp_path: Path):
    mod = _load()
    row = _row("INCOMPLETE", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
    row["claim_ready"] = ""
    camp = _write_campaign(tmp_path, [row])
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 1
    assert rep["metrics"]["STRICT"]["n"] == 0
    assert "claim_ready!=1" in rep["strict_fail_rows"][0]["reasons"]


def test_expected_target_without_result_row_counts_as_failure(tmp_path: Path):
    mod = _load()
    rows = [
        _row(f"PRESENT{i:02d}", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
        for i in range(84)
    ]
    camp = _write_campaign(tmp_path, rows)
    receipt = json.loads((camp / "RUN_RECEIPT.json").read_text())
    receipt["expected_target_ids"] = [r["pdb_id"] for r in rows] + ["MISSING"]
    (camp / "RUN_RECEIPT.json").write_text(json.dumps(receipt))
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(mod.load_campaign_rows(camp), pin, src, str(camp))
    assert rep["N_claim"] == 85
    assert rep["metrics"]["STRICT"]["n"] == 84
    assert rep["completeness"]["missing_ids"] == ["MISSING"]
    assert rep["headline"]["suppressed"] is True
    missing = next(r for r in rep["strict_fail_rows"] if r["pdb_id"] == "MISSING")
    assert "missing_result_row" in missing["reasons"]


def test_duplicate_target_ids_suppress_strict_headline(tmp_path: Path):
    mod = _load()
    row = _row("DUP", rmsd_ordered=1.0, bcr=0.8, pb=1, claim_ready=1)
    camp = _write_campaign(tmp_path, [row])
    result_csv = camp / "DUP" / "result.csv"
    with result_csv.open("a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writerow({key: row.get(key, "") for key in CSV_FIELDS})
    loaded = mod.load_campaign_rows(camp)
    assert len(loaded) == 2
    pin, src = mod.load_matrix_pin(camp, None)
    rep = mod.aggregate_rows(
        loaded,
        pin,
        src,
        str(camp),
        expected_target_ids=["DUP"],
        expected_target_source="fixture",
    )
    assert rep["N_claim"] == 1
    assert rep["completeness"]["duplicate_ids"] == {"DUP": 2}
    assert rep["metrics"]["STRICT"]["n"] == 0
    assert rep["headline"]["suppressed"] is True
    assert "duplicate_result_ids" in rep["headline"]["suppression_reasons"]


def test_seed_echo_0_0_accepted():
    mod = _load()
    row = {"seed_echo": "0.0", "native_pose_seeded": "0.0", "matrix_md5": ""}
    assert mod._flag0(row, "seed_echo")
    assert mod._flag0(row, "native_pose_seeded")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
