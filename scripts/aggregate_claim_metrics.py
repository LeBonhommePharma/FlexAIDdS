#!/usr/bin/env python3
"""Aggregate claim-table metrics under the admission + S1/S2/STRICT/S3 contract.

Normative contract:
  benchmarks/protocols/admission_metrics_contract.md

Claim denominator (pre-outcome):
  every expected/observed target unless explicitly preregistered ineligible;
  missing rows and legacy schema remain denominator failures

STRICT numerator (fail-closed):
  claim_ready == 1, no seed/oracle contamination, matrix pin, complete
  same-pose hashes, upstream PoseBusters pass + version/config receipt,
  tENCoM == ok, and Eigen == ok

Metrics (always separate):
  S1      ordered direct rmsd_to_crystal ≤ 2.0 Å  (RMSD-only diagnostic)
  S2      S1 ∧ pb_pass / success_pb
  STRICT  claim_ready == 1  ← primary headline
  S3      conditional scanned-pool ceiling ≤ 2.0 Å (diagnostic only; never any-pose)

S1 MUST use rmsd_to_crystal only — never rmsd_hungarian.
--headline s3 requires --diagnostic-only.

Usage:
  python3 scripts/aggregate_claim_metrics.py <campaign_dir> [--json out.json]
  python3 scripts/aggregate_claim_metrics.py --c0-full85
  python3 scripts/aggregate_claim_metrics.py <dir> --headline s1
  python3 scripts/aggregate_claim_metrics.py <dir> --headline s3 --diagnostic-only

Exit codes:
  0  complete, claimable report (N_claim > 0)
  1  empty or strict headline suppressed by incomplete target/schema evidence
  2  usage or contract violation
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable

# Production matrix pin used by three-engine / C0 full85 (see RUN_RECEIPT).
DEFAULT_MATRIX_MD5 = "9dc93717dfed0698006d88dd6a9627bc"
RMSD_SUCCESS_A = 2.0
C0_FULL85_REL = "campaigns/C0_full85_defined_cleft_nativeseed_forbidden"
PB_SCHEMA_ID = "posebusters-0.6.5-redock-csv-v1"
PB_REQUIRED_CHECK_COUNT = 27
PB_PACKAGE_NAME = "posebusters"
PB_PACKAGE_VERSION = "0.6.5"
PB_LAUNCHER_VERSION = "bust 0.6.5"
PB_CONFIG_NAME = "redock"
PB_CONFIG_SHA256 = "4d551d898ff29a404f16e02ad5a7a2d4235e6b7b14e9a3e27f7c66b4d16b2da9"

SUMMARY_CSV_NAMES = (
    "astex_diverse_results.csv",
    "astex_crossdock_85_results.csv",
    "results.csv",
    "summary.csv",
    "claim_summary.csv",
)


# ── helpers ──────────────────────────────────────────────────────────────────


def _f(row: dict[str, str], *keys: str) -> float:
    for k in keys:
        v = row.get(k)
        if v is None or v == "" or str(v).upper() == "NA":
            continue
        try:
            x = float(v)
            if math.isfinite(x):
                return x
        except (TypeError, ValueError):
            continue
    return float("nan")


def _truth(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip() in ("1", "True", "true", "YES", "yes")


def _flag0(row: dict[str, str], key: str) -> bool:
    """True when the admission flag is *explicitly* zero / false.

    Fail-closed: missing or blank keys fail admission (return False).
    Explicit non-zero also fails. Accept common zero spellings including "0.0".
    """
    if key not in row or row.get(key) is None:
        return False
    s = str(row.get(key, "")).strip()
    if s == "":
        return False
    return s in ("0", "0.0", "False", "false", "NO", "no")


def _is_sha256_hex(value: str) -> bool:
    digest = value.strip().lower()
    return len(digest) == 64 and all(c in "0123456789abcdef" for c in digest)


def _pdb_id(row: dict[str, str]) -> str:
    return str(row.get("pdb_id") or row.get("pdb") or row.get("target") or "?").strip()


def resolve_c0_full85_dir() -> Path:
    """Resolve C0 full85 campaign dir from FLEXAIDDS_RESULTS (or FLEXAIDDS_ICLOUD)."""
    results = os.environ.get("FLEXAIDDS_RESULTS", "").strip()
    if results:
        return Path(results) / "campaigns" / "C0_full85_defined_cleft_nativeseed_forbidden"
    icloud = os.environ.get("FLEXAIDDS_ICLOUD", "").strip()
    if icloud:
        # FLEXAIDDS_ICLOUD may be the benchmarks root or the nested astex_entropy path
        base = Path(icloud)
        candidates = [
            base / "results" / C0_FULL85_REL,
            base.parent / "results" / C0_FULL85_REL if base.name else None,
        ]
        for c in candidates:
            if c is not None and c.is_dir():
                return c
        return base / "results" / C0_FULL85_REL
    # Documented default relative to home iCloud Drive (not hardcoded user)
    home = Path.home()
    return (
        home
        / "Library"
        / "Mobile Documents"
        / "com~apple~CloudDocs"
        / "FlexAIDdS_benchmarks"
        / "results"
        / "campaigns"
        / "C0_full85_defined_cleft_nativeseed_forbidden"
    )


def _normalize_matrix_pin(md: str) -> str:
    s = str(md).strip().lower()
    if len(s) != 32 or any(c not in "0123456789abcdef" for c in s):
        raise ValueError(f"matrix_md5 pin must be 32 hex chars, got {md!r}")
    return s


def load_matrix_pin(campaign_dir: Path, cli_pin: str | None) -> tuple[str, str]:
    """Return (md5, source_label)."""
    if cli_pin:
        return _normalize_matrix_pin(cli_pin), "cli"
    for name in ("RUN_RECEIPT.json", "provenance.json"):
        p = campaign_dir / name
        if not p.is_file():
            continue
        try:
            data = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        md = data.get("matrix_md5")
        if md:
            return _normalize_matrix_pin(str(md)), name
    return _normalize_matrix_pin(DEFAULT_MATRIX_MD5), "default_pin"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _attach_pb_receipt(row: dict[str, str], target_dir: Path) -> dict[str, str]:
    """Attach the immutable upstream PoseBusters receipt to a result row.

    Receipt data is deliberately kept in private ``_pb_receipt_*`` keys so the
    input CSV schema remains unchanged. Missing, malformed, or incomplete
    receipts stay visible to the strict numerator instead of removing the row
    from the preregistered denominator.
    """
    enriched = dict(row)
    enriched["_pb_receipt_present"] = "0"
    receipt = target_dir / "posebust" / f"{_pdb_id(row)}_bust_receipt.json"
    if not receipt.is_file():
        return enriched
    try:
        data = json.loads(receipt.read_text())
    except (OSError, json.JSONDecodeError):
        enriched["_pb_receipt_malformed"] = "1"
        return enriched
    if not isinstance(data, dict):
        enriched["_pb_receipt_malformed"] = "1"
        return enriched
    enriched["_pb_receipt_present"] = "1"

    def obj(key: str) -> dict[str, Any]:
        value = data.get(key, {})
        return value if isinstance(value, dict) else {}

    schema = obj("schema")
    package = obj("package")
    config = obj("config")
    command = obj("command")
    inputs = obj("inputs")
    outputs = obj("outputs")
    result = obj("result")
    scalar_fields = {
        "schema_id": schema.get("id", ""),
        "schema_required_check_count": schema.get("required_check_count", ""),
        "package_name": package.get("name", ""),
        "package_version": package.get("version", ""),
        "package_record_path": package.get("record_path", ""),
        "package_record_sha256": package.get("record_sha256", ""),
        "package_launcher_path": package.get("launcher_path", ""),
        "package_launcher_sha256": package.get("launcher_sha256", ""),
        "package_launcher_version_output": package.get(
            "launcher_version_output", ""
        ),
        "config_name": config.get("name", ""),
        "config_path": config.get("path", ""),
        "config_sha256": config.get("sha256", ""),
        "command_exit_status": command.get("exit_status", ""),
        "result_backend": result.get("backend", ""),
        "result_ran": result.get("ran", ""),
        "result_pb_pass": result.get("pb_pass", ""),
    }
    argv = command.get("argv", [])
    enriched["_pb_receipt_command_argv"] = (
        json.dumps(argv) if isinstance(argv, list) else ""
    )
    for key, value in scalar_fields.items():
        enriched[f"_pb_receipt_{key}"] = "" if value is None else str(value)

    identity_objects: list[tuple[str, dict[str, Any]]] = []
    for name in ("predicted_ligand", "protein", "crystal_ligand"):
        value = inputs.get(name, {})
        identity_objects.append(
            (f"input_{name}", value if isinstance(value, dict) else {})
        )
    for name in ("raw_csv", "validated_csv"):
        value = outputs.get(name, {})
        identity_objects.append(
            (f"output_{name}", value if isinstance(value, dict) else {})
        )
    for name, identity in identity_objects:
        path_value = str(identity.get("path", "") or "").strip()
        hash_value = str(identity.get("sha256", "") or "").strip()
        enriched[f"_pb_receipt_{name}_path"] = path_value
        enriched[f"_pb_receipt_{name}_sha256"] = hash_value
        if name.startswith("output_") and path_value:
            output_path = Path(path_value)
            if not output_path.is_absolute():
                output_path = receipt.parent / output_path
            if output_path.is_file():
                try:
                    actual = _sha256_file(output_path)
                except OSError:
                    actual = ""
                enriched[f"_pb_receipt_{name}_actual_sha256"] = actual
                enriched[f"_pb_receipt_{name}_hash_matches"] = str(
                    bool(actual) and actual == hash_value
                )
    return enriched


def load_campaign_rows(out_dir: Path) -> list[dict[str, str]]:
    """Load per-target result.csv trees first, then flat summary CSVs."""
    rows: list[dict[str, str]] = []
    for rc in sorted(out_dir.glob("*/result.csv")):
        try:
            batch = list(csv.DictReader(rc.open(newline="")))
            for raw in batch:
                rows.append(_attach_pb_receipt(dict(raw), rc.parent))
        except OSError:
            continue
    if rows:
        return rows
    for name in SUMMARY_CSV_NAMES:
        p = out_dir / name
        if p.is_file():
            try:
                return [dict(r) for r in csv.DictReader(p.open(newline=""))]
            except OSError:
                continue
    # Also accept a single CSV path masquerading as "dir" (handled by caller)
    return rows


def load_rows_from_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def _target_ids_from_json(value: Any) -> list[str]:
    if isinstance(value, dict):
        for key in ("expected_target_ids", "target_ids", "targets"):
            if key in value:
                return _target_ids_from_json(value[key])
        return []
    if not isinstance(value, list):
        return []
    ids: list[str] = []
    for item in value:
        if isinstance(item, str):
            pid = item.strip()
        elif isinstance(item, dict):
            pid = str(
                item.get("pdb_id") or item.get("pdb") or item.get("target") or ""
            ).strip()
        else:
            pid = ""
        if pid:
            ids.append(pid)
    return ids


def load_expected_target_ids(
    campaign_dir: Path, manifest: Path | None
) -> tuple[list[str] | None, str | None]:
    """Load preregistered target IDs from an argument or campaign receipt."""
    if manifest is not None:
        if not manifest.is_file():
            raise ValueError(f"expected-target manifest not found: {manifest}")
        try:
            if manifest.suffix.lower() == ".json":
                ids = _target_ids_from_json(json.loads(manifest.read_text()))
            elif manifest.suffix.lower() == ".csv":
                with manifest.open(newline="") as fh:
                    batch = list(csv.DictReader(fh))
                ids = [_pdb_id(r) for r in batch if _pdb_id(r) != "?"]
            else:
                ids = [
                    token.strip()
                    for line in manifest.read_text().splitlines()
                    for token in line.split(",")
                    if token.strip() and not token.lstrip().startswith("#")
                ]
        except (OSError, json.JSONDecodeError, csv.Error) as exc:
            raise ValueError(f"cannot read expected-target manifest {manifest}: {exc}")
        if not ids:
            raise ValueError(f"expected-target manifest is empty: {manifest}")
        return ids, str(manifest)

    for name in ("RUN_RECEIPT.json", "provenance.json"):
        path = campaign_dir / name
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        ids = _target_ids_from_json(data)
        if ids:
            return ids, name
    return None, None


def elected_rmsd(row: dict[str, str]) -> float:
    """Ordered direct elected RMSD only (audit P0).

    Never use rmsd_hungarian for S1 / claim rates. Legacy three-engine
    `rmsd_top1` is accepted only when it is the ordered top-1 serial metric
    and `rmsd_to_crystal` is absent.
    """
    rc = _f(row, "rmsd_to_crystal")
    if math.isfinite(rc):
        return rc
    # Legacy engines without rmsd_to_crystal may emit ordered rmsd_top1 only.
    return _f(row, "rmsd_top1")


def is_s1(row: dict[str, str]) -> bool:
    """S1: ordered direct elected RMSD ≤ 2.0 Å (RMSD-only diagnostic).

    Finite ordered RMSD always wins over success_* flags so a stale flag
    cannot admit a high-RMSD pose. Hungarian is never consulted.
    """
    if _truth(row, "seed_echo"):
        return False
    rh = elected_rmsd(row)
    if math.isfinite(rh):
        return 0.0 <= rh <= RMSD_SUCCESS_A
    # No finite ordered RMSD: fall back to engine flags only
    if "success_s1" in row and str(row.get("success_s1", "")).strip() != "":
        return _truth(row, "success_s1")
    if "success_rmsd" in row and str(row.get("success_rmsd", "")).strip() != "":
        return _truth(row, "success_rmsd")
    if "success" in row and str(row.get("success", "")).strip() != "":
        return _truth(row, "success")
    return False


def is_s2(row: dict[str, str], s1: bool) -> bool:
    """S2: S1 ∧ PoseBusters pass on the same elected pose."""
    if not s1:
        return False
    if "success_pb" in row and str(row.get("success_pb", "")).strip() != "":
        return _truth(row, "success_pb") and s1
    if "pb_pass" in row and str(row.get("pb_pass", "")).strip() != "":
        return _truth(row, "pb_pass")
    return False


def is_s3(row: dict[str, str]) -> bool:
    """S3: conditional scanned-pool ceiling (diagnostic only; not any-pose)."""
    bc = _f(
        row,
        "conditional_scanned_pool_ceiling",
        "best_cluster_rmsd",
        "rmsd_bcr",
    )
    if math.isfinite(bc):
        return 0.0 <= bc <= RMSD_SUCCESS_A
    if "success_s3" in row and str(row.get("success_s3", "")).strip() != "":
        return _truth(row, "success_s3")
    return False


def _hash_receipts_ok(row: dict[str, str]) -> tuple[bool, list[str]]:
    """Require every validator to cite the same nonempty elected-pose hash."""
    reasons: list[str] = []
    pose = str(row.get("pose_sha256", "")).strip()
    if not pose:
        reasons.append("pose_sha256_missing")
    elif not _is_sha256_hex(pose):
        reasons.append("pose_sha256_invalid")
    for key in (
        "rmsd_pose_sha256",
        "posebusters_pose_sha256",
        "tencom_pose_sha256",
    ):
        value = str(row.get(key, "")).strip()
        if not value:
            reasons.append(f"{key}_missing")
            continue
        if not _is_sha256_hex(value):
            reasons.append(f"{key}_invalid")
        if pose and value != pose:
            reasons.append(f"{key}_mismatch")
    return (len(reasons) == 0, reasons)


def _pb_receipt_reasons(row: dict[str, str]) -> list[str]:
    """Validate the nested PoseBusters 0.6.5 redock identity receipt."""
    if not _truth(row, "_pb_receipt_present"):
        return ["pb_receipt_missing"]
    reasons: list[str] = []
    if _truth(row, "_pb_receipt_malformed"):
        reasons.append("pb_receipt_malformed")

    exact = {
        "schema_id": PB_SCHEMA_ID,
        "schema_required_check_count": str(PB_REQUIRED_CHECK_COUNT),
        "package_name": PB_PACKAGE_NAME,
        "package_version": PB_PACKAGE_VERSION,
        "package_launcher_version_output": PB_LAUNCHER_VERSION,
        "config_name": PB_CONFIG_NAME,
        "config_sha256": PB_CONFIG_SHA256,
        "result_backend": "bust_cli",
        "command_exit_status": "0",
    }
    for key, expected in exact.items():
        got = str(row.get(f"_pb_receipt_{key}", "")).strip()
        if got != expected:
            reasons.append(f"pb_receipt_{key}!={expected}")
    if not _truth(row, "_pb_receipt_result_ran"):
        reasons.append("pb_receipt_result_ran!=1")
    if not _truth(row, "_pb_receipt_result_pb_pass"):
        reasons.append("pb_receipt_result_pb_pass!=1")

    try:
        argv = json.loads(str(row.get("_pb_receipt_command_argv", "")))
    except json.JSONDecodeError:
        argv = []
    if not isinstance(argv, list) or not argv or any(not str(x).strip() for x in argv):
        reasons.append("pb_receipt_command_argv_missing")

    identities = (
        "package_record",
        "package_launcher",
        "config",
        "input_predicted_ligand",
        "input_protein",
        "input_crystal_ligand",
        "output_raw_csv",
        "output_validated_csv",
    )
    for name in identities:
        path = str(row.get(f"_pb_receipt_{name}_path", "")).strip()
        digest = str(row.get(f"_pb_receipt_{name}_sha256", "")).strip().lower()
        if not path:
            reasons.append(f"pb_receipt_{name}_path_missing")
        if not _is_sha256_hex(digest):
            reasons.append(f"pb_receipt_{name}_sha256_invalid")

    for name in ("output_raw_csv", "output_validated_csv"):
        match = str(row.get(f"_pb_receipt_{name}_hash_matches", "")).strip()
        if match and not _truth(row, f"_pb_receipt_{name}_hash_matches"):
            reasons.append(f"pb_receipt_{name}_rehash_mismatch")
    return reasons


def is_claim_ready(row: dict[str, str]) -> bool:
    """STRICT claim success: engine claim_ready when present."""
    if "claim_ready" in row and str(row.get("claim_ready", "")).strip() != "":
        return _truth(row, "claim_ready")
    return False


def strict_failure_reasons(row: dict[str, str], matrix_pin: str) -> list[str]:
    """Return fail-closed STRICT numerator reasons without changing eligibility."""
    reasons: list[str] = []
    if _truth(row, "_missing_result_row"):
        reasons.append("missing_result_row")
    duplicate_count = str(row.get("_duplicate_result_rows", "")).strip()
    if duplicate_count:
        reasons.append(f"duplicate_result_rows={duplicate_count}")
    if not _truth(row, "protocol_claim_eligible"):
        reasons.append("protocol_claim_eligible!=1")
    if not _flag0(row, "seed_echo"):
        reasons.append("seed_echo!=0")
    if not _flag0(row, "native_pose_seeded"):
        reasons.append("native_pose_seeded!=0")
    if not row_matrix_ok(row, matrix_pin):
        reasons.append(
            f"matrix_md5_mismatch(got={row.get('matrix_md5')!r}, pin={matrix_pin})"
        )
    if not is_claim_ready(row):
        reasons.append("claim_ready!=1")
    if not is_s1(row):
        reasons.append("ordered_rmsd_success!=1")
    if not _truth(row, "pb_pass"):
        reasons.append("pb_pass!=1")
    if str(row.get("pb_backend", "")).strip() != "bust_cli":
        reasons.append("pb_backend!=bust_cli")
    _, hash_reasons = _hash_receipts_ok(row)
    reasons.extend(hash_reasons)
    if str(row.get("tencom_status", "")).strip().lower() != "ok":
        reasons.append("tencom_status!=ok")
    if str(row.get("eigen_status", "")).strip().lower() != "ok":
        reasons.append("eigen_status!=ok")
    reasons.extend(_pb_receipt_reasons(row))
    return reasons


def is_strict_success(row: dict[str, str], matrix_pin: str) -> bool:
    """STRICT numerator; every required receipt is explicit and fail-closed."""
    return not strict_failure_reasons(row, matrix_pin)


def row_matrix_ok(row: dict[str, str], pin: str) -> bool:
    md = str(row.get("matrix_md5", "")).strip().lower()
    if not md:
        return True  # campaign-level pin applies
    return md == pin


def is_claim_eligible(row: dict[str, str], matrix_pin: str) -> tuple[bool, list[str]]:
    """Pre-outcome denominator gate for the preregistered claim table.

    ``matrix_pin`` remains in the signature for API compatibility; matrix and
    every docking/validator outcome are strict-numerator checks, never reasons
    to make a preregistered target vanish from the denominator.
    """
    del matrix_pin
    if _truth(row, "_expected_target"):
        return True, []
    value = str(row.get("protocol_claim_eligible", "")).strip()
    if value and not _truth(row, "protocol_claim_eligible"):
        return False, ["protocol_claim_eligible=0"]
    # A legacy/mixed-schema row is an observed preregistered target until proven
    # explicitly ineligible. It stays in N and fails STRICT for missing evidence.
    return True, []


def aggregate_rows(
    rows: Iterable[dict[str, str]],
    matrix_pin: str,
    matrix_pin_source: str,
    campaign_dir: str | None = None,
    expected_target_ids: Iterable[str] | None = None,
    expected_target_source: str | None = None,
) -> dict[str, Any]:
    all_rows = list(rows)
    claim: list[dict[str, str]] = []
    dropped: list[dict[str, Any]] = []
    if expected_target_ids is None and campaign_dir:
        expected_target_ids, discovered_source = load_expected_target_ids(
            Path(campaign_dir), None
        )
        if expected_target_source is None:
            expected_target_source = discovered_source

    rows_by_id: dict[str, list[dict[str, str]]] = {}
    for row in all_rows:
        rows_by_id.setdefault(_pdb_id(row), []).append(row)
    duplicate_ids = {
        pid: len(batch) for pid, batch in rows_by_id.items() if len(batch) > 1
    }

    expected_raw = (
        [str(pid).strip() for pid in expected_target_ids if str(pid).strip()]
        if expected_target_ids is not None
        else None
    )
    expected_duplicate_ids: dict[str, int] = {}
    expected_unique: list[str] | None = None
    missing_ids: list[str] = []
    extra_ids: list[str] = []
    if expected_raw is not None:
        expected_unique = list(dict.fromkeys(expected_raw))
        expected_duplicate_ids = {
            pid: expected_raw.count(pid)
            for pid in expected_unique
            if expected_raw.count(pid) > 1
        }
        expected_set = set(expected_unique)
        missing_ids = [pid for pid in expected_unique if pid not in rows_by_id]
        extra_candidates = [pid for pid in rows_by_id if pid not in expected_set]
        for pid in expected_unique:
            if pid in rows_by_id:
                row = dict(rows_by_id[pid][0])
            else:
                row = {
                    "pdb_id": pid,
                    "protocol_claim_eligible": "1",
                    "_missing_result_row": "1",
                }
            row["_expected_target"] = "1"
            if pid in duplicate_ids:
                row["_duplicate_result_rows"] = str(duplicate_ids[pid])
            claim.append(row)
        for pid in extra_candidates:
            ok, reasons = is_claim_eligible(rows_by_id[pid][0], matrix_pin)
            if ok:
                extra_ids.append(pid)
                reasons = ["not_in_expected_target_ids"]
            dropped.append({"pdb_id": pid, "reasons": reasons})
    else:
        for pid, batch in rows_by_id.items():
            row = dict(batch[0])
            if pid in duplicate_ids:
                row["_duplicate_result_rows"] = str(duplicate_ids[pid])
            ok, reasons = is_claim_eligible(row, matrix_pin)
            if ok:
                claim.append(row)
            else:
                dropped.append({"pdb_id": pid, "reasons": reasons})

    n = len(claim)
    s1_ids: list[str] = []
    s2_ids: list[str] = []
    strict_ids: list[str] = []
    s3_ids: list[str] = []
    election_gap_ids: list[str] = []
    s1_fail_ids: list[str] = []
    strict_fail_rows: list[dict[str, Any]] = []
    legacy_rows: list[dict[str, str]] = []

    for r in claim:
        if (
            not _truth(r, "_missing_result_row")
            and not str(r.get("protocol_claim_eligible", "")).strip()
        ):
            legacy_rows.append(r)

    for r in claim:
        pid = _pdb_id(r)
        s1 = is_s1(r)
        s2 = is_s2(r, s1)
        s3 = is_s3(r)
        strict_reasons = strict_failure_reasons(r, matrix_pin)
        strict = not strict_reasons
        if s1:
            s1_ids.append(pid)
        else:
            s1_fail_ids.append(pid)
        if s2:
            s2_ids.append(pid)
        if strict:
            strict_ids.append(pid)
        else:
            strict_fail_rows.append({"pdb_id": pid, "reasons": strict_reasons})
        if s3:
            s3_ids.append(pid)
        if s3 and not s1:
            election_gap_ids.append(pid)

    def rate(k: int) -> float:
        return (k / n) if n else 0.0

    headline_suppression_reasons: list[str] = []
    if expected_unique is None:
        headline_suppression_reasons.append("expected_target_ids_unavailable")
    if missing_ids:
        headline_suppression_reasons.append("missing_expected_targets")
    if duplicate_ids:
        headline_suppression_reasons.append("duplicate_result_ids")
    if expected_duplicate_ids:
        headline_suppression_reasons.append("duplicate_expected_target_ids")
    if extra_ids:
        headline_suppression_reasons.append("unexpected_result_ids")
    if legacy_rows:
        headline_suppression_reasons.append("legacy_schema_rows")

    report: dict[str, Any] = {
        "contract": "admission_metrics_contract",
        "contract_doc": "benchmarks/protocols/admission_metrics_contract.md",
        "campaign_dir": campaign_dir,
        "matrix_md5_pin": matrix_pin,
        "matrix_md5_pin_source": matrix_pin_source,
        "N_raw": len(all_rows),
        "N_claim": n,
        "N_dropped": len(dropped),
        "N_legacy_no_claim_ready": sum(
            1
            for r in legacy_rows
            if not str(r.get("claim_ready", "")).strip()
        ),
        "N_legacy_diagnostic": len(legacy_rows),
        "dropped_rows": dropped,
        "strict_fail_rows": strict_fail_rows,
        "completeness": {
            "verified": expected_unique is not None,
            "expected_target_source": expected_target_source,
            "N_expected": len(expected_unique) if expected_unique is not None else None,
            "N_observed_unique": len(rows_by_id),
            "missing_ids": missing_ids,
            "duplicate_ids": duplicate_ids,
            "expected_manifest_duplicate_ids": expected_duplicate_ids,
            "extra_ids": extra_ids,
            "complete": (
                expected_unique is not None
                and not missing_ids
                and not duplicate_ids
                and not expected_duplicate_ids
                and not extra_ids
            ),
        },
        "metrics": {
            "S1": {
                "definition": "ordered direct rmsd_to_crystal <= 2.0 A (never hungarian)",
                "role": "rmsd_only_diagnostic",
                "n": len(s1_ids),
                "rate": rate(len(s1_ids)),
                "ids": s1_ids,
            },
            "S2": {
                "definition": "S1 AND PoseBusters pass on same elected pose",
                "role": "rmsd_and_pb_diagnostic",
                "n": len(s2_ids),
                "rate": rate(len(s2_ids)),
                "ids": s2_ids,
            },
            "STRICT": {
                "definition": (
                    "claim_ready==1 with no seed/oracle contamination, complete "
                    "same-pose hashes, upstream PB pass/version/config receipt, "
                    "tENCoM==ok, and Eigen==ok"
                ),
                "role": "primary_headline",
                "n": len(strict_ids),
                "rate": rate(len(strict_ids)),
                "ids": strict_ids,
            },
            "S3": {
                "definition": (
                    "conditional_scanned_pool_ceiling / best_cluster_rmsd <= 2.0 A "
                    "(scanned heads/members only; NOT any-pose)"
                ),
                "role": "diagnostic_only",
                "n": len(s3_ids),
                "rate": rate(len(s3_ids)),
                "ids": s3_ids,
                "warning": (
                    "Do not report S3 as abstract / headline success. "
                    "Ceiling is conditional on scanned emission pool."
                ),
            },
        },
        "election_gap": {
            "definition": "S3=1 and S1=0 (scanned pool had near-native; elector missed)",
            "n": len(election_gap_ids),
            "rate": rate(len(election_gap_ids)),
            "ids": election_gap_ids,
        },
        "S1_fail_ids": s1_fail_ids,
        "legacy_diagnostics": {
            "definition": (
                "rows without explicit protocol_claim_eligible; retained in "
                "claim denominator but excluded from STRICT"
            ),
            "N": len(legacy_rows),
            "S1_n": sum(1 for r in legacy_rows if is_s1(r)),
            "S2_n": sum(1 for r in legacy_rows if is_s2(r, is_s1(r))),
            "S3_n": sum(1 for r in legacy_rows if is_s3(r)),
            "ids": [_pdb_id(r) for r in legacy_rows],
        },
        "headline": {
            "metric": "STRICT",
            "n": len(strict_ids),
            "N": n,
            "rate": rate(len(strict_ids)),
            "label": "strict successes / all preregistered protocol-eligible targets",
            "claimable": not headline_suppression_reasons,
            "suppressed": bool(headline_suppression_reasons),
            "suppression_reasons": headline_suppression_reasons,
            "schema_incomplete": len(legacy_rows) > 0,
        },
        "denominator": {
            "protocol_claim_eligible": 1,
            "outcome_fields_used": [],
        },
        "strict_numerator": {
            "seed_echo": 0,
            "native_pose_seeded": 0,
            "matrix_md5": matrix_pin,
            "claim_ready": 1,
            "pb_backend": "bust_cli",
            "pb_pass": 1,
            "pb_receipt_required": True,
            "tencom_status": "ok",
            "eigen_status": "ok",
        },
    }
    return report


def format_text_report(report: dict[str, Any]) -> str:
    m = report["metrics"]
    n = report["N_claim"]
    lines = [
        f"campaign_dir: {report.get('campaign_dir')}",
        f"matrix_md5_pin: {report['matrix_md5_pin']} (source={report['matrix_md5_pin_source']})",
        f"N_raw={report['N_raw']}  N_claim={n}  N_dropped={report['N_dropped']}",
        "completeness: "
        + ("verified" if report["completeness"]["complete"] else "INCOMPLETE/UNVERIFIED")
        + f" (source={report['completeness']['expected_target_source']})",
        "",
        f"STRICT (headline): {m['STRICT']['n']}/{n} = {100.0 * m['STRICT']['rate']:.2f}%  "
        f"[claim_ready]",
        f"S1 (RMSD-only):    {m['S1']['n']}/{n} = {100.0 * m['S1']['rate']:.2f}%  "
        f"[diagnostic; ordered rmsd_to_crystal]",
        f"S2 (S1∧PB):        {m['S2']['n']}/{n} = {100.0 * m['S2']['rate']:.2f}%",
        f"S3 (pool ceiling): {m['S3']['n']}/{n} = {100.0 * m['S3']['rate']:.2f}%  "
        f"[diagnostic; conditional scanned pool — NOT any-pose]",
        f"election_gap (S3∧¬S1): {report['election_gap']['n']}/{n} = "
        f"{100.0 * report['election_gap']['rate']:.2f}%",
    ]
    if report.get("headline", {}).get("suppressed"):
        lines.extend(
            [
                "",
                "STRICT HEADLINE SUPPRESSED: "
                + ", ".join(report["headline"].get("suppression_reasons", [])),
            ]
        )
    if report["N_dropped"]:
        lines.append("")
        lines.append("dropped (non-claim):")
        for d in report["dropped_rows"][:20]:
            lines.append(f"  {_pdb_id_from_drop(d)}: {', '.join(d['reasons'])}")
        if report["N_dropped"] > 20:
            lines.append(f"  ... +{report['N_dropped'] - 20} more")
    return "\n".join(lines) + "\n"


def _pdb_id_from_drop(d: dict[str, Any]) -> str:
    return str(d.get("pdb_id") or "?")


def apply_headline(
    report: dict[str, Any],
    headline: str,
    diagnostic_only: bool,
) -> tuple[dict[str, Any], int | None]:
    """Mutate report headline; return (report, error_exit_code_or_None)."""
    h = headline.strip().lower()
    if h in ("strict", "claim_ready", "claim-ready"):
        h = "strict"
    if h not in ("s1", "s2", "s3", "strict"):
        return report, 2
    if h == "s3" and not diagnostic_only:
        print(
            "CONTRACT VIOLATION: --headline s3 is not allowed as primary without "
            "--diagnostic-only. S3 is conditional scanned-pool ceiling only; "
            "use --headline strict (claim_ready) for claim success.",
            file=sys.stderr,
        )
        return report, 2
    if h in ("s1", "s2") and not diagnostic_only:
        # Allowed but must not be mistaken for claim success.
        pass

    key = "STRICT" if h == "strict" else h.upper()
    m = report["metrics"][key]
    strict_claimable = bool(report.get("headline", {}).get("claimable", False))
    suppression_reasons = list(
        report.get("headline", {}).get("suppression_reasons", [])
    )
    report["headline"] = {
        "metric": key,
        "n": m["n"],
        "N": report["N_claim"],
        "rate": m["rate"],
        "label": m["definition"],
        "role": m["role"],
        "diagnostic_only": h in ("s1", "s2", "s3") or diagnostic_only,
        "claimable": h == "strict" and strict_claimable and not diagnostic_only,
        "suppressed": h == "strict" and not strict_claimable,
        "suppression_reasons": suppression_reasons if h == "strict" else [],
    }
    if h == "s3":
        report["headline"]["warning"] = (
            "S3 is diagnostic only — conditional scanned-pool ceiling, not any-pose."
        )
    if h in ("s1", "s2"):
        report["headline"]["warning"] = (
            f"{key} is not the claim headline; prefer STRICT (claim_ready)."
        )
    return report, None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "campaign_dir",
        type=Path,
        nargs="?",
        default=None,
        help="Campaign directory with */result.csv and/or summary CSV",
    )
    ap.add_argument(
        "--c0-full85",
        action="store_true",
        help=f"Use $FLEXAIDDS_RESULTS/{C0_FULL85_REL}",
    )
    ap.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Aggregate a single summary CSV instead of a campaign tree",
    )
    ap.add_argument(
        "--matrix-md5",
        type=str,
        default=None,
        help=f"Matrix pin (default: receipt/provenance or {DEFAULT_MATRIX_MD5})",
    )
    ap.add_argument(
        "--expected-targets",
        type=Path,
        default=None,
        help=(
            "Preregistered target IDs as JSON, CSV, or newline/comma text. "
            "Otherwise expected_target_ids/target_ids/targets is read from RUN_RECEIPT.json."
        ),
    )
    ap.add_argument(
        "--headline",
        type=str,
        default="strict",
        choices=(
            "strict",
            "STRICT",
            "claim_ready",
            "s1",
            "s2",
            "s3",
            "S1",
            "S2",
            "S3",
        ),
        help="Primary headline (default strict/claim_ready). s3 requires --diagnostic-only.",
    )
    ap.add_argument(
        "--diagnostic-only",
        action="store_true",
        help="Allow --headline s3/s1/s2 as diagnostic headline (never abstract claim success).",
    )
    ap.add_argument("--json", type=Path, default=None, help="Write full JSON report")
    ap.add_argument(
        "--quiet",
        action="store_true",
        help="JSON only to stdout (no human text block)",
    )
    args = ap.parse_args(argv)

    n_sources = sum(
        [
            args.csv is not None,
            bool(args.c0_full85),
            args.campaign_dir is not None,
        ]
    )
    if n_sources > 1:
        print(
            "error: provide exactly one of campaign_dir, --csv, or --c0-full85",
            file=sys.stderr,
        )
        return 2

    campaign: Path | None = None
    rows: list[dict[str, str]] = []

    if args.csv is not None:
        if not args.csv.is_file():
            print(f"not a file: {args.csv}", file=sys.stderr)
            return 2
        rows = load_rows_from_csv(args.csv)
        campaign = args.csv.parent
    elif args.c0_full85:
        campaign = resolve_c0_full85_dir()
        if not campaign.is_dir():
            print(
                f"C0 full85 campaign not found: {campaign}\n"
                "Source ~/.flexaidds_env or set FLEXAIDDS_RESULTS.",
                file=sys.stderr,
            )
            return 2
        rows = load_campaign_rows(campaign)
    elif args.campaign_dir is not None:
        campaign = args.campaign_dir
        if campaign.is_file() and campaign.suffix.lower() == ".csv":
            rows = load_rows_from_csv(campaign)
            campaign = campaign.parent
        elif campaign.is_dir():
            rows = load_campaign_rows(campaign)
        else:
            print(f"not a directory or CSV: {campaign}", file=sys.stderr)
            return 2
    else:
        ap.print_help()
        print(
            "\nerror: provide campaign_dir, --csv, or --c0-full85",
            file=sys.stderr,
        )
        return 2

    assert campaign is not None
    try:
        pin, pin_src = load_matrix_pin(campaign, args.matrix_md5)
        expected_ids, expected_src = load_expected_target_ids(
            campaign, args.expected_targets
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    report = aggregate_rows(
        rows,
        matrix_pin=pin,
        matrix_pin_source=pin_src,
        campaign_dir=str(campaign.resolve()),
        expected_target_ids=expected_ids,
        expected_target_source=expected_src,
    )
    report, err = apply_headline(report, args.headline, args.diagnostic_only)
    if err is not None:
        return err

    if not args.quiet:
        print(format_text_report(report), end="")
        print(
            f"headline: {report['headline']['metric']} "
            f"{report['headline']['n']}/{report['headline']['N']} "
            f"= {100.0 * report['headline']['rate']:.2f}%"
        )

    text = json.dumps(report, indent=2)
    if args.json:
        args.json.write_text(text + "\n")
    if args.quiet:
        print(text)

    if report["headline"]["metric"] == "STRICT" and report["headline"]["suppressed"]:
        return 1
    return 0 if report["N_claim"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
