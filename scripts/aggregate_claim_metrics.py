#!/usr/bin/env python3
"""Aggregate target-level S1/S2/STRICT/S3 under METHODOLOGY.md section 0.3.

STRICT is recomputed from complete, internally consistent receipt fields. This
CSV tool does not authenticate validator execution or hash the underlying pose,
receptor, and raw validator artifacts. See evidence_level in every report.

One observation per target is the default. Repeated targets require a declared
--expected-seeds list; missing seeds fail and strictly more than half must pass.
Never mix arms or endpoints. Directory layouts must identify a single source;
use --csv to choose explicitly when both per-target and summary files exist.

The frozen roster is mandatory. --legacy-observed-denominator is available only
with --diagnostic-only and a diagnostic headline; it cannot produce STRICT rates.
Exit codes: 0 = strict receipt success (or diagnostic observations), 1 = no strict
receipt successes, 2 = input/contract error.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

DEFAULT_MATRIX_MD5 = "9dc93717dfed0698006d88dd6a9627bc"
RMSD_SUCCESS_A = 2.0
EXPECTED_PB_CHECKS = 27  # Versioned mandatory schema in LIB/PoseBust/BustCli.cpp.
SCORE_DELTA_TOLERANCE = 1e-4
C0_FULL85_REL = "campaigns/C0_full85_defined_cleft_nativeseed_forbidden"
SUMMARY_CSV_NAMES = ("astex_diverse_results.csv", "astex_crossdock_85_results.csv",
                     "results.csv", "summary.csv", "claim_summary.csv")
SYMMCORR_COL = "rmsd_symmcorr"
SYMMCORR_COLS = (SYMMCORR_COL, "rmsd_spyrmsd")
SHA256 = re.compile(r"[0-9a-fA-F]{64}\Z")


def _f(row: dict[str, str], *keys: str) -> float:
    for key in keys:
        value = row.get(key)
        if value is None or str(value).strip() in ("", "NA"):
            continue
        try:
            result = float(value)
        except (ValueError, TypeError):
            continue
        if math.isfinite(result):
            return result
    return float("nan")


def _truth(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip().lower() in ("1", "true", "yes")


def _flag0(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip().lower() in ("0", "0.0", "false", "no")


def _pdb_id(row: dict[str, str]) -> str:
    return str(row.get("pdb_id") or row.get("pdb") or row.get("target") or "?").strip().upper()


def _hash(row: dict[str, str], key: str) -> str:
    return str(row.get(key, "")).strip().lower()


def _valid_sha(value: str) -> bool:
    return bool(SHA256.fullmatch(value))


def resolve_c0_full85_dir() -> Path:
    # Live aggregation is local first. CloudDocs needs icloud_safe_io.py staging.
    results = os.environ.get("FLEXAIDDS_RESULTS", "").strip()
    if results:
        return Path(results) / C0_FULL85_REL
    local = Path(os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(Path.home() / "flexaidds_results")))
    return local / C0_FULL85_REL


def _require_local(path: Path) -> None:
    resolved = path.resolve()
    if "Mobile Documents" in resolved.parts or "com~apple~CloudDocs" in resolved.parts:
        raise ValueError("stage CloudDocs inputs locally with scripts/icloud_safe_io.py before aggregation")


def _normalize_matrix_pin(md: str) -> str:
    value = str(md).strip().lower()
    if not re.fullmatch(r"[0-9a-f]{32}", value):
        raise ValueError(f"matrix_md5 pin must be 32 hex characters, got {md!r}")
    return value


def load_matrix_pin(campaign_dir: Path, cli_pin: str | None) -> tuple[str, str]:
    _require_local(campaign_dir)
    pins: dict[str, str] = {}
    if cli_pin:
        pins["cli"] = _normalize_matrix_pin(cli_pin)
    for name in ("RUN_RECEIPT.json", "provenance.json"):
        path = campaign_dir / name
        if path.is_file():
            _require_local(path)
            data = json.loads(path.read_text())
            if not isinstance(data, dict):
                raise ValueError(f"{path}: receipt must be an object")
            if data.get("matrix_md5"):
                pins[name] = _normalize_matrix_pin(str(data["matrix_md5"]))
    if len(set(pins.values())) > 1:
        raise ValueError(f"conflicting matrix pins: {pins}")
    if pins:
        return next(iter(pins.values())), "+".join(pins)
    return _normalize_matrix_pin(DEFAULT_MATRIX_MD5), "default_expected_pin (row identity still required)"


def load_rows_from_csv(path: Path) -> list[dict[str, str]]:
    _require_local(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader, None)
        if not header or any(not key or key != key.strip() for key in header):
            raise ValueError(f"{path}: missing, blank, or whitespace-padded CSV header")
        if len(set(header)) != len(header):
            raise ValueError(f"{path}: duplicate CSV header")
        if not {"pdb_id", "pdb", "target"}.intersection(header):
            raise ValueError(f"{path}: missing target identity column")
        rows = []
        for values in reader:
            if len(values) != len(header):
                raise ValueError(f"{path}:{reader.line_num}: CSV row width {len(values)} != {len(header)}")
            row = dict(zip(header, values))
            identities = {str(row[k]).strip().upper() for k in ("pdb_id", "pdb", "target") if row.get(k)}
            if len(identities) != 1 or "?" in identities:
                raise ValueError(f"{path}:{reader.line_num}: missing or conflicting target identity")
            rows.append(row)
    return rows


def load_campaign_rows(out_dir: Path) -> list[dict[str, str]]:
    _require_local(out_dir)
    per_target = sorted(out_dir.glob("*/result.csv"))
    summaries = [out_dir / name for name in SUMMARY_CSV_NAMES if (out_dir / name).is_file()]
    if (per_target and summaries) or len(summaries) > 1:
        raise ValueError("ambiguous campaign sources; choose the intended source explicitly with --csv")
    sources = per_target or summaries
    rows = []
    for source in sources:
        batch = load_rows_from_csv(source)
        if source in per_target and len(batch) != 1:
            raise ValueError(f"{source}: per-target result.csv must contain exactly one observation")
        rows.extend(batch)
    return rows


def elected_rmsd_labelled(row: dict[str, str]) -> tuple[float, str]:
    value = _f(row, *SYMMCORR_COLS)
    if math.isfinite(value):
        return value, "symmcorr"
    value = _f(row, "rmsd_to_crystal")
    if math.isfinite(value):
        return value, "serial"
    return _f(row, "rmsd_top1"), "serial_legacy_top1"


def elected_rmsd(row: dict[str, str]) -> float:
    return elected_rmsd_labelled(row)[0]


def is_s1(row: dict[str, str]) -> bool:
    value = elected_rmsd(row)
    return _flag0(row, "seed_echo") and math.isfinite(value) and 0.0 <= value <= RMSD_SUCCESS_A


def load_symmcorr_sidecar(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    out = {}
    for row in load_rows_from_csv(path):
        if str(row.get("status", "")).strip() != "ok":
            continue
        value = _f(row, *SYMMCORR_COLS)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{path}: status=ok requires a finite, nonnegative symmetry-corrected RMSD")
        sha = _hash(row, "pose_sha256")
        if not _valid_sha(sha):
            raise ValueError(f"{path}: sidecar requires a valid pose_sha256")
        key = (_pdb_id(row), sha)
        if key in out:
            raise ValueError(f"{path}: duplicate sidecar target/pose identity {key}")
        row[SYMMCORR_COL] = str(value)
        out[key] = row
    return out


def join_symmcorr(rows: list[dict[str, str]], sidecar: dict) -> dict[str, Any]:
    # Accept the earlier direct-Python target-keyed API while enforcing hashes.
    index: dict[tuple[str, str], dict[str, str]] = {}
    targets = set()
    for rec in sidecar.values():
        pid = _pdb_id(rec)
        targets.add(pid)
        key = (pid, _hash(rec, "pose_sha256"))
        if key in index:
            raise ValueError(f"duplicate sidecar target/pose identity {key}")
        index[key] = rec
    joined = 0
    refused = []
    for row in rows:
        pid, sha = _pdb_id(row), _hash(row, "pose_sha256")
        if pid not in targets:
            continue
        rec = index.get((pid, sha))
        if not _valid_sha(sha) or rec is None or not _valid_sha(_hash(rec, "pose_sha256")):
            refused.append(pid)
            continue
        value = _f(rec, *SYMMCORR_COLS)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"invalid symmetry-corrected RMSD for {pid}")
        row[SYMMCORR_COL] = str(value)
        joined += 1
    return {"joined": joined, "refused_sha_mismatch": sorted(set(refused)), "sidecar_rows": len(sidecar)}


def is_s2(row: dict[str, str], s1: bool) -> bool:
    # Diagnostics are independent of STRICT's tENCoM/score gates, not immune
    # to contradictions in their own PB measurements.
    if "pb_ran" in row and not _truth(row, "pb_ran"):
        return False
    for key in ("pb_n_checks", "pb_n_pass"):
        if key in row and _f(row, key) != EXPECTED_PB_CHECKS:
            return False
    if "pb_n_fail" in row and _f(row, "pb_n_fail") != 0:
        return False
    if str(row.get("pb_failed_keys", "")).strip():
        return False
    if "pb_backend" in row and str(row["pb_backend"]).strip() != "bust_cli":
        return False
    pose = _hash(row, "pose_sha256")
    return (s1 and _truth(row, "pb_pass") and _valid_sha(pose)
            and _hash(row, "posebusters_pose_sha256") == pose)


def is_s3(row: dict[str, str]) -> bool:
    value = _f(row, "conditional_scanned_pool_ceiling", "best_cluster_rmsd", "rmsd_bcr")
    return math.isfinite(value) and 0.0 <= value <= RMSD_SUCCESS_A


def _hash_receipts_ok(row: dict[str, str]) -> tuple[bool, list[str]]:
    reasons = []
    pose = _hash(row, "pose_sha256")
    if not _valid_sha(pose):
        reasons.append("pose_sha256_missing_or_invalid")
    for key in ("rmsd_pose_sha256", "posebusters_pose_sha256", "tencom_pose_sha256"):
        value = _hash(row, key)
        if not _valid_sha(value):
            reasons.append(f"{key}_missing_or_invalid")
        elif value != pose:
            reasons.append(f"{key}_mismatch")
    if not _valid_sha(_hash(row, "posebusters_input_sha256")):
        reasons.append("posebusters_input_sha256_missing_or_invalid")
    return not reasons, reasons


def row_matrix_ok(row: dict[str, str], pin: str) -> bool:
    try:
        return _normalize_matrix_pin(str(row.get("matrix_md5", ""))) == pin
    except ValueError:
        return False


def is_claim_eligible(row: dict[str, str], matrix_pin: str) -> tuple[bool, list[str]]:
    """Protocol eligibility only: failures remain visible in diagnostic counts."""
    reasons = []
    for key in ("seed_echo", "native_pose_seeded"):
        if not _flag0(row, key):
            reasons.append(f"{key}!=0")
    if not _truth(row, "protocol_claim_eligible"):
        reasons.append("protocol_claim_eligible!=1")
    if not row_matrix_ok(row, matrix_pin):
        reasons.append("matrix_md5_missing_or_mismatch")
    if "native_pose_seed_fraction" in row:
        if _f(row, "native_pose_seed_fraction") != 0:
            reasons.append("native_pose_seed_fraction!=0")
    return not reasons, reasons


def strict_failure_reasons(row: dict[str, str], matrix_pin: str = DEFAULT_MATRIX_MD5) -> list[str]:
    _, reasons = is_claim_eligible(row, matrix_pin)
    # The producer flag is required as an attestation, never sufficient evidence.
    for key in ("claim_ready", "success_rmsd", "success_pb", "pb_pass", "score_pose_consistent"):
        if not _truth(row, key):
            reasons.append(f"{key}!=1")
    if "docking_completed" in row and not _truth(row, "docking_completed"):
        reasons.append("docking_completed!=1")
    if "docking_exit_code" in row and _f(row, "docking_exit_code") != 0:
        reasons.append("docking_exit_code!=0")
    # Current STRICT requires the full PB schema and actual Eigen output fields.
    # Missing older fields remain useful only for separately labelled diagnostics.
    for key in ("pb_n_checks", "pb_n_pass"):
        if _f(row, key) != EXPECTED_PB_CHECKS:
            reasons.append(f"{key}!={EXPECTED_PB_CHECKS}")
    if not _truth(row, "pb_ran"):
        reasons.append("pb_ran!=1")
    if _f(row, "pb_n_fail") != 0:
        reasons.append("pb_n_fail!=0")
    modes = _f(row, "eigen_n_modes")
    if not math.isfinite(modes) or modes <= 0 or not modes.is_integer():
        reasons.append("eigen_n_modes_not_positive_integer")
    if not math.isfinite(_f(row, "elected_H_vib")):
        reasons.append("elected_H_vib_not_finite")
    if "num_poses" in row:
        poses = _f(row, "num_poses")
        if not math.isfinite(poses) or poses <= 0 or not poses.is_integer():
            reasons.append("num_poses_not_positive_integer")
    if str(row.get("pb_failed_keys", "")).strip():
        reasons.append("pb_failed_keys_nonempty")
    if "rmsd_fail_reason" in row and str(row["rmsd_fail_reason"]).strip() != "none":
        reasons.append("rmsd_fail_reason!=none")
    serial = _f(row, "rmsd_to_crystal")
    if not math.isfinite(serial) or not 0 <= serial <= RMSD_SUCCESS_A:
        reasons.append("serial_rmsd_not_finite_or_outside_cutoff")
    delta = _f(row, "score_pose_delta")
    if not math.isfinite(delta) or abs(delta) > SCORE_DELTA_TOLERANCE:
        reasons.append("score_pose_delta_not_finite_or_outside_tolerance")
    for key, value in (("tencom_status", "ok"), ("eigen_status", "ok"), ("pb_backend", "bust_cli")):
        if str(row.get(key, "")).strip() != value:
            reasons.append(f"{key}!={value}")
    _, hash_reasons = _hash_receipts_ok(row)
    reasons.extend(hash_reasons)
    return reasons


def is_claim_ready(row: dict[str, str]) -> bool:
    return _truth(row, "claim_ready")


def is_strict_success(row: dict[str, str], matrix_pin: str = DEFAULT_MATRIX_MD5) -> bool:
    return not strict_failure_reasons(row, matrix_pin)


def load_target_manifest(path: Path | None = None) -> tuple[list[str], str]:
    here = Path(__file__).resolve().parent
    candidates = ([path] if path is not None else [
        here.parent / "benchmarks/protocols/astex85_target_manifest.json",
        here / "astex85_target_manifest.json",
    ])
    manifest_path = next((p for p in candidates if p.is_file()), None)
    if manifest_path is None:
        raise ValueError("frozen target manifest missing; supply --manifest (no implicit denominator fallback)")
    _require_local(manifest_path)
    data = json.loads(manifest_path.read_text())
    if (not isinstance(data, dict)
            or data.get("schema") != "flexaidds.astex.target_manifest/v1"
            or not isinstance(data.get("targets"), list)):
        raise ValueError("invalid target manifest schema")
    codes = data["targets"]
    if (not codes or any(not isinstance(c, str) or not re.fullmatch(r"[A-Z0-9]{4}", c) for c in codes)
            or len(set(codes)) != len(codes) or type(data.get("N")) is not int
            or data.get("N") != len(codes)):
        raise ValueError("manifest must contain unique uppercase four-character targets and matching N")
    codes = sorted(codes)
    digest = hashlib.sha256(",".join(codes).encode("ascii")).hexdigest()
    if data.get("sha256_of_sorted_codes") != digest:
        raise ValueError("target manifest digest mismatch (SHA256 of comma-joined sorted codes)")
    return codes, digest


def _seed_id(value: Any) -> str:
    """Canonical uint64 decimal identity, matching the engine's seed domain."""
    text = str(value).strip()
    if not text:
        return ""  # Legacy one-observation producers may omit seed metadata.
    if not re.fullmatch(r"[0-9]+", text) or int(text) > 2**64 - 1:
        raise ValueError(f"seed must be an unsigned 64-bit decimal integer, got {value!r}")
    return str(int(text))


def _observation_groups(rows: list[dict[str, str]], expected_seeds: list[str] | None,
                        arm: str | None) -> tuple[list[dict[str, str]], dict, dict]:
    selected = rows
    if arm is not None:
        selected = [r for r in rows if str(r.get("arm", "")).strip() == arm]
        if not selected:
            raise ValueError(f"no observations for selected arm {arm!r}")
    for field in ("arm", "endpoint"):
        labels = {str(r.get(field, "")).strip() for r in selected}
        if len(labels) > 1:
            raise ValueError(f"mixed {field} identities; select one explicit arm/endpoint source")
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen = set()
    for row in selected:
        pid = _pdb_id(row)
        if pid == "?":
            raise ValueError("missing target identity")
        key = (pid, _seed_id(row.get("seed", "")))
        if key in seen:
            raise ValueError(f"duplicate observation target/seed identity: {key}")
        seen.add(key)
        groups[pid].append(row)
    if expected_seeds is None:
        if any(len(group) != 1 for group in groups.values()):
            raise ValueError("repeated targets require an explicit --expected-seeds list")
        seed_labels = {_seed_id(r.get("seed", "")) for r in selected}
        if len(seed_labels) > 1:
            raise ValueError("mixed seed identities require an explicit --expected-seeds list")
        mode, required = "single_observation", 1
    else:
        expected_seeds = [_seed_id(s) for s in expected_seeds]
        if not expected_seeds or any(not s for s in expected_seeds) or len(set(expected_seeds)) != len(expected_seeds):
            raise ValueError("expected seeds must be a nonempty unique list")
        if any(_seed_id(r.get("seed", "")) not in expected_seeds for r in selected):
            raise ValueError("observation has missing or unexpected seed")
        mode, required = "majority_of_expected_seeds", len(expected_seeds) // 2 + 1
    meta = {"mode": mode, "expected_seeds": expected_seeds, "required_passes": required,
            "arm": next(iter({str(r.get('arm', '')).strip() for r in selected}), "") or "unspecified",
            "endpoint": next(iter({str(r.get('endpoint', '')).strip() for r in selected}), "") or "unspecified",
            "N_filtered_by_arm": len(rows) - len(selected), "missing_expected_observations": [],
            "missing_policy": "missing expected target/seed observations do not pass"}
    return selected, groups, meta


def aggregate_rows(rows: Iterable[dict[str, str]], matrix_pin: str, matrix_pin_source: str,
                   campaign_dir: str | None = None, fixed_denominator: bool = True,
                   *, manifest_path: Path | None = None, expected_seeds: list[str] | None = None,
                   arm: str | None = None) -> dict[str, Any]:
    all_rows = list(rows)
    matrix_pin = _normalize_matrix_pin(matrix_pin)
    selected, groups, aggregation = _observation_groups(all_rows, expected_seeds, arm)
    if fixed_denominator:
        codes, digest = load_target_manifest(manifest_path)
        denominator_source = f"frozen_manifest(N={len(codes)},sha={digest})"
    else:
        codes, digest = sorted(groups), None
        denominator_source = "observed_targets (explicit diagnostic legacy mode; STRICT unavailable)"
    roster = set(codes)
    denominator = len(codes)
    missing = sorted(roster - set(groups))
    dropped, strict_rows, eligible_rows = [], [], []
    verdicts: dict[str, list[dict]] = defaultdict(list)
    for row in selected:
        pid = _pdb_id(row)
        eligible, _ = is_claim_eligible(row, matrix_pin)
        reasons = strict_failure_reasons(row, matrix_pin)
        if pid not in roster:
            reasons.append("off_manifest")
        if reasons:
            dropped.append({"pdb_id": pid, "seed": row.get("seed", ""), "reasons": reasons})
        elif fixed_denominator:
            strict_rows.append(row)
        if eligible and pid in roster:
            eligible_rows.append(row)
        s1 = eligible and is_s1(row)
        verdicts[pid].append({"S1": s1, "S2": eligible and is_s2(row, s1),
                              "S3": eligible and is_s3(row), "STRICT": not reasons and fixed_denominator})
    required = aggregation["required_passes"]
    if expected_seeds is not None:
        for pid in codes:
            present = {_seed_id(r.get("seed", "")) for r in groups.get(pid, [])}
            for seed in aggregation["expected_seeds"]:
                if seed not in present:
                    aggregation["missing_expected_observations"].append({"pdb_id": pid, "seed": seed})
    definitions = {
        "S1": "finite elected in-place graph-symmetry RMSD <=2 A, or labelled serial fallback; protocol-eligible rows",
        "S2": "S1 AND pb_pass with matching elected-pose receipt hash (diagnostic)",
        "STRICT": "recomputed serial RMSD/PB/score/validator/protocol/matrix receipt checks; underlying artifacts not authenticated",
        "S3": "finite conditional scanned-pool ceiling <=2 A (diagnostic; not any-pose)",
    }
    metrics = {}
    for key in definitions:
        ids = sorted(pid for pid in roster if sum(v[key] for v in verdicts[pid]) >= required)
        rate = len(ids) / denominator if denominator else 0.0
        if key == "STRICT" and not fixed_denominator:
            rate = None
        metrics[key] = {"definition": definitions[key], "role": "primary_receipt_headline" if key == "STRICT" else "diagnostic_only",
                        "n": len(ids), "rate": rate, "ids": ids}
    metric_used = "+".join(sorted({elected_rmsd_labelled(r)[1] for r in eligible_rows
                                  if math.isfinite(elected_rmsd(r))})) or "none_finite"
    metrics["S1"]["metric_used"] = metric_used
    gap_ids = sorted(set(metrics["S3"]["ids"]) - set(metrics["S1"]["ids"]))
    symmcorr_gained = sorted({_pdb_id(r) for r in eligible_rows
                             if math.isfinite(_f(r, *SYMMCORR_COLS))
                             and 0 <= _f(r, *SYMMCORR_COLS) <= RMSD_SUCCESS_A
                             and not (math.isfinite(_f(r, "rmsd_to_crystal"))
                                      and 0 <= _f(r, "rmsd_to_crystal") <= RMSD_SUCCESS_A)})
    report = {"contract": "admission_metrics_contract/v2", "contract_doc": "METHODOLOGY.md#03-admission-identity-missingness-and-repair-evidence",
              "campaign_dir": campaign_dir, "matrix_md5_pin": matrix_pin, "matrix_md5_pin_source": matrix_pin_source,
              "evidence_level": "validated_receipt_fields", "artifacts_verified": False,
              "evidence_limit": "CSV consistency is not independent proof of validator execution, raw-artifact integrity, receptor condition, or scientific protocol compliance.",
              "N_raw": len(all_rows), "N_selected": len(selected), "N_claim": len(strict_rows),
              "N_protocol_eligible": len(eligible_rows), "N_denominator": denominator,
              "N_denominator_source": denominator_source, "manifest_sha256_of_sorted_codes": digest,
              "N_missing_from_manifest": len(missing), "missing_targets": missing,
              "off_manifest_targets": sorted(set(groups) - roster), "N_dropped": len(dropped), "dropped_rows": dropped,
              "N_legacy_no_claim_ready": sum(not str(r.get("claim_ready", "")).strip() for r in selected),
              "aggregation": aggregation, "metrics": metrics,
              "election_gap": {"definition": "target S3 majority passes and target S1 majority fails", "n": len(gap_ids),
                               "rate": len(gap_ids) / denominator if denominator else 0.0, "ids": gap_ids},
              "S1_fail_ids": sorted(roster - set(metrics["S1"]["ids"])),
              "symmcorr_delta": {"definition": "distinct targets with a measured serial-fail / symmcorr-pass observation (not a majority verdict)",
                                 "n": len(symmcorr_gained), "ids": symmcorr_gained},
              "strict_metric_inheritance": {"note": "STRICT recomputes the producer's serial gate; symmetry sidecars affect diagnostics only.",
                                            "s1_metric": metric_used, "strict_metric": "serial (recomputed)"},
              "legacy_diagnostic_mode": not fixed_denominator}
    report["headline"] = {"metric": "STRICT", "n": metrics["STRICT"]["n"], "N": denominator,
                          "rate": metrics["STRICT"]["rate"], "label": definitions["STRICT"]}
    return report


def apply_headline(report: dict[str, Any], headline: str, diagnostic_only: bool) -> tuple[dict[str, Any], int | None]:
    key = headline.upper().replace("-", "_")
    if key == "CLAIM_READY":
        key = "STRICT"
    if key not in report["metrics"] or (key == "S3" and not diagnostic_only):
        print("CONTRACT VIOLATION: S3 headline requires --diagnostic-only", file=sys.stderr)
        return report, 2
    if report.get("legacy_diagnostic_mode") and (key == "STRICT" or not diagnostic_only):
        print("CONTRACT VIOLATION: legacy denominator requires a diagnostic headline and --diagnostic-only", file=sys.stderr)
        return report, 2
    metric = report["metrics"][key]
    report["headline"] = {"metric": key, "n": metric["n"], "N": report["N_denominator"], "rate": metric["rate"],
                          "label": metric["definition"], "role": metric["role"],
                          "diagnostic_only": key != "STRICT" or diagnostic_only}
    return report, None


def format_text_report(report: dict[str, Any]) -> str:
    lines = [f"campaign_dir: {report['campaign_dir']}",
             f"evidence_level: {report['evidence_level']}; artifacts_verified=false",
             f"N_raw={report['N_raw']} N_claim={report['N_claim']} N_protocol_eligible={report['N_protocol_eligible']}",
             f"denominator: {report['N_denominator']} ({report['N_denominator_source']})",
             f"aggregation: {report['aggregation']['mode']}; arm={report['aggregation']['arm']}; endpoint={report['aggregation']['endpoint']}"]
    for key in ("STRICT", "S1", "S2", "S3"):
        metric = report["metrics"][key]
        rate = "unavailable" if metric["rate"] is None else f"{100 * metric['rate']:.2f}%"
        lines.append(f"{key}: {metric['n']}/{report['N_denominator']} = {rate}; {metric['definition']}")
    lines.append(report["evidence_limit"])
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("campaign_dir", type=Path, nargs="?")
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--c0-full85", action="store_true")
    parser.add_argument("--matrix-md5")
    parser.add_argument("--manifest", type=Path, help="Frozen target roster with N and validated sorted-code SHA256")
    parser.add_argument("--expected-seeds", help="Comma-separated prespecified seeds; target success requires a strict majority")
    parser.add_argument("--arm", help="Select this explicit arm value; report excluded row count")
    parser.add_argument("--headline", default="strict", choices=("strict", "STRICT", "claim_ready", "s1", "s2", "s3", "S1", "S2", "S3"))
    parser.add_argument("--diagnostic-only", action="store_true")
    parser.add_argument("--legacy-observed-denominator", action="store_true")
    parser.add_argument("--symmcorr", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    try:
        if sum([args.csv is not None, args.campaign_dir is not None, args.c0_full85]) != 1:
            raise ValueError("provide exactly one of campaign_dir, --csv, or --c0-full85")
        if args.legacy_observed_denominator and not args.diagnostic_only:
            raise ValueError("--legacy-observed-denominator requires --diagnostic-only")
        source = args.csv or (resolve_c0_full85_dir() if args.c0_full85 else args.campaign_dir)
        _require_local(source)
        if source.is_file():
            rows, campaign = load_rows_from_csv(source), source.parent
        elif source.is_dir() and args.csv is None:
            rows, campaign = load_campaign_rows(source), source
        else:
            raise ValueError(f"not a CSV or campaign directory: {source}")
        provenance = {"sidecar": None, "joined": 0, "refused_sha_mismatch": [], "sidecar_rows": 0}
        # Filter explicitly chosen arms before joining; an unrelated arm is not a pose mismatch.
        if args.arm is not None:
            side_rows = [r for r in rows if str(r.get("arm", "")).strip() == args.arm]
        else:
            side_rows = rows
        if args.symmcorr is not None:
            provenance = join_symmcorr(side_rows, load_symmcorr_sidecar(args.symmcorr))
            provenance["sidecar"] = str(args.symmcorr)
            if provenance["refused_sha_mismatch"]:
                raise ValueError("symmcorr pose_sha256 missing, invalid, or mismatched for: " + ", ".join(provenance["refused_sha_mismatch"]))
        pin, pin_source = load_matrix_pin(campaign, args.matrix_md5)
        report = aggregate_rows(rows, pin, pin_source, str(campaign.resolve()),
                                fixed_denominator=not args.legacy_observed_denominator,
                                manifest_path=args.manifest,
                                expected_seeds=args.expected_seeds.split(",") if args.expected_seeds is not None else None,
                                arm=args.arm)
        report["symmcorr"] = provenance
        report, error = apply_headline(report, args.headline, args.diagnostic_only)
        if error:
            return error
        output = json.dumps(report, indent=2, allow_nan=False)
        if args.json is not None:
            _require_local(args.json)
            args.json.write_text(output + "\n")
        if args.quiet:
            print(output)
        else:
            print(format_text_report(report), end="")
        return 0 if report["metrics"]["STRICT"]["n"] > 0 or (args.diagnostic_only and report["N_selected"] > 0) else 1
    except (ValueError, OSError, csv.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
