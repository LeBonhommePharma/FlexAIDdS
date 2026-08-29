#!/usr/bin/env python3
"""Aggregate claim-table metrics under the admission + S1/S2/STRICT/S3 contract.

Normative contract:
  benchmarks/protocols/admission_metrics_contract.md

Claim admission (fail-closed):
  seed_echo == 0, native_pose_seeded == 0, matrix pin, protocol_claim_eligible
  claim_ready == 1 (strict table) + PoseBusters/tENCoM/Eigen/hash receipts when present

Metrics (always separate):
  S1      symmetry-corrected elected-pose RMSD ≤ 2.0 Å  (RMSD-only diagnostic)
  S2      S1 ∧ pb_pass / success_pb
  STRICT  claim_ready == 1  ← primary headline
  S3      conditional scanned-pool ceiling ≤ 2.0 Å (diagnostic only; never any-pose)

S1 metric — read this before quoting a number (#365).
  "Success at 2.0 Å" used to have two producers with two meanings. This file
  MANDATED the serial column ("S1 MUST use rmsd_to_crystal only") while
  METHODOLOGY.md §claim, docs/swarm/2026-08-13/score_canonical.py and
  ops/gate_accuracy_rmsd.py all meant spyRMSD. Both were live. S1 now converges
  on the symmetry-corrected value, produced by exactly one implementation:

    scripts/rmsd_symmcorr.py  →  spyrmsd.rmsd.symmrmsd
    (METHODOLOGY.md §0 method 2 contract: crystal SDF bond block, heavy atoms,
     center=False, minimize=False — in-place, never superposed)

  Supply it with --symmcorr <sidecar.csv>. Rows are joined on pdb_id and, when
  both sides carry it, on pose_sha256 — so a sidecar from a different run
  cannot be silently attached to this claim table.

  `rmsd_to_crystal` KEEPS ITS MEANING: ordered/serial, identity atom mapping,
  as the engine documents (LIB/DatasetRunner.cpp). It is not redefined here.
  Without a sidecar S1 falls back to it, which is SAFE BUT CONSERVATIVE:
  symmetry correction minimises over graph automorphisms and the identity
  mapping is one of them, so symmcorr ≤ serial always. The serial gate can
  therefore under-count successes but never over-count them. The report says
  which metric produced the number; an unlabelled RMSD is not reportable.

  rmsd_hungarian remains BANNED for S1. Element-only Hungarian minimises over
  all same-element bijections — a superset of the chemically valid
  automorphisms — so it is over-permissive, not merely different: the repo
  measured it inflating the pool ceiling from 48.8% to 57.8%.
  Ordering, by construction: hungarian ≤ symmcorr ≤ serial.

--headline s3 requires --diagnostic-only.

Usage:
  python3 scripts/aggregate_claim_metrics.py <campaign_dir> [--json out.json]
  python3 scripts/aggregate_claim_metrics.py --c0-full85
  python3 scripts/aggregate_claim_metrics.py <dir> --headline s1
  python3 scripts/aggregate_claim_metrics.py <dir> --headline s3 --diagnostic-only

Exit codes:
  0  OK (N_claim > 0)
  1  no claim-eligible rows / empty campaign
  2  usage or contract violation
"""
from __future__ import annotations

import argparse
import csv
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


def load_campaign_rows(out_dir: Path) -> list[dict[str, str]]:
    """Load per-target result.csv trees first, then flat summary CSVs."""
    rows: list[dict[str, str]] = []
    for rc in sorted(out_dir.glob("*/result.csv")):
        try:
            batch = list(csv.DictReader(rc.open(newline="")))
            if batch:
                # One authoritative row per target dir (first row)
                rows.append(dict(batch[0]))
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


SYMMCORR_COL = "rmsd_symmcorr"
# The frozen nine-arm table (PER_POSE_CF_RMSD_NINE_ARMS.csv) and
# docs/swarm/2026-08-13/score_canonical.py already call this same quantity
# `rmsd_spyrmsd`. Accept both names rather than minting a rival one — the whole
# point of #365 is one quantity with one meaning.
SYMMCORR_COLS = (SYMMCORR_COL, "rmsd_spyrmsd")


def elected_rmsd_labelled(row: dict[str, str]) -> tuple[float, str]:
    """Elected-pose RMSD for S1, with the name of the metric that produced it.

    Preference order (#365):
      1. `rmsd_symmcorr`  — symmetry-corrected, spyrmsd graph automorphism,
         joined from scripts/rmsd_symmcorr.py. This is the claim metric.
      2. `rmsd_to_crystal` — ordered/serial identity mapping. Conservative
         fallback: symmcorr ≤ serial always, so this can only under-count.
      3. `rmsd_top1` — legacy three-engine ordered top-1, when neither is present.

    rmsd_hungarian is never consulted: it is over-permissive (see module
    docstring). Returns (value, metric_name); value is nan when none is finite.
    """
    rs = _f(row, *SYMMCORR_COLS)
    if math.isfinite(rs):
        return rs, "symmcorr"
    rc = _f(row, "rmsd_to_crystal")
    if math.isfinite(rc):
        return rc, "serial"
    # Legacy engines without rmsd_to_crystal may emit ordered rmsd_top1 only.
    return _f(row, "rmsd_top1"), "serial_legacy_top1"


def elected_rmsd(row: dict[str, str]) -> float:
    """Elected-pose RMSD used by S1. See elected_rmsd_labelled()."""
    return elected_rmsd_labelled(row)[0]


def is_s1(row: dict[str, str]) -> bool:
    """S1: symmetry-corrected elected RMSD ≤ 2.0 Å (RMSD-only diagnostic).

    A finite RMSD always wins over success_* flags so a stale flag cannot
    admit a high-RMSD pose. In particular the engine's `success_rmsd` column
    is a SERIAL gate; when a symmetry-corrected value is present it supersedes
    that flag rather than deferring to it — otherwise the engine's stricter
    definition would silently keep governing the claim (#365).

    Hungarian is never consulted.
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


def load_symmcorr_sidecar(path: Path) -> dict[str, dict[str, str]]:
    """Load scripts/rmsd_symmcorr.py output, keyed by pdb_id.

    Only rows with status == "ok" and a finite value are usable. Everything
    else is dropped here so a blank or errored row cannot be mistaken for a
    measurement.
    """
    out: dict[str, dict[str, str]] = {}
    with path.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            if str(rec.get("status", "")).strip() != "ok":
                continue
            val = ""
            for col in SYMMCORR_COLS:
                val = str(rec.get(col, "")).strip()
                if val:
                    break
            if val == "":
                continue
            rec = dict(rec)
            rec[SYMMCORR_COL] = val
            pid = str(rec.get("pdb_id", "")).strip().upper()
            if pid:
                out[pid] = dict(rec)
    return out


def join_symmcorr(
    rows: list[dict[str, str]], sidecar: dict[str, dict[str, str]]
) -> dict[str, Any]:
    """Attach `rmsd_symmcorr` to claim rows. Returns a provenance summary.

    The join requires pose identity: when both the claim row and the sidecar
    carry a pose_sha256 they must be equal, otherwise the sidecar value
    describes a different pose and is refused. This is the same same-pose
    discipline the engine applies via rmsd_pose_sha256 == pose_sha256.
    """
    joined = 0
    refused: list[str] = []
    for row in rows:
        pid = _pdb_id(row).upper()
        rec = sidecar.get(pid)
        if rec is None:
            continue
        row_sha = str(row.get("pose_sha256", "")).strip()
        side_sha = str(rec.get("pose_sha256", "")).strip()
        if row_sha and side_sha and row_sha != side_sha:
            refused.append(pid)
            continue
        row[SYMMCORR_COL] = rec["rmsd_symmcorr"]
        joined += 1
    return {
        "joined": joined,
        "refused_sha_mismatch": sorted(refused),
        "sidecar_rows": len(sidecar),
    }


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
    """When hash columns are present, require identity with pose_sha256."""
    reasons: list[str] = []
    pose = str(row.get("pose_sha256", "")).strip()
    if not pose:
        return True, reasons  # absent → not checked at row level
    for key in (
        "rmsd_pose_sha256",
        "posebusters_pose_sha256",
        "tencom_pose_sha256",
    ):
        if key not in row or str(row.get(key, "")).strip() == "":
            continue
        if str(row.get(key, "")).strip() != pose:
            reasons.append(f"{key}_mismatch")
    return (len(reasons) == 0, reasons)


def is_claim_ready(row: dict[str, str]) -> bool:
    """STRICT claim success: engine claim_ready when present."""
    if "claim_ready" in row and str(row.get("claim_ready", "")).strip() != "":
        return _truth(row, "claim_ready")
    return False


def is_strict_success(row: dict[str, str]) -> bool:
    """STRICT: claim_ready with receipts when available."""
    if not is_claim_ready(row):
        return False
    ok, _ = _hash_receipts_ok(row)
    if not ok:
        return False
    if "tencom_status" in row and str(row.get("tencom_status", "")).strip() != "":
        if str(row.get("tencom_status", "")).strip().lower() != "ok":
            return False
    if "eigen_status" in row and str(row.get("eigen_status", "")).strip() != "":
        if str(row.get("eigen_status", "")).strip().lower() != "ok":
            return False
    if "pb_backend" in row and str(row.get("pb_backend", "")).strip() != "":
        if str(row.get("pb_backend", "")).strip() != "bust_cli":
            return False
    return True


def row_matrix_ok(row: dict[str, str], pin: str) -> bool:
    md = str(row.get("matrix_md5", "")).strip().lower()
    if not md:
        return True  # campaign-level pin applies
    return md == pin


def is_claim_eligible(row: dict[str, str], matrix_pin: str) -> tuple[bool, list[str]]:
    """Admission gates for the claim table. Returns (ok, fail reasons).

    Strict admission requires claim_ready==1 when the column is present.
    """
    reasons: list[str] = []
    if not _flag0(row, "seed_echo"):
        reasons.append("seed_echo!=0")
    if not _flag0(row, "native_pose_seeded"):
        reasons.append("native_pose_seeded!=0")
    if not row_matrix_ok(row, matrix_pin):
        reasons.append(
            f"matrix_md5_mismatch(got={row.get('matrix_md5')!r}, pin={matrix_pin})"
        )
    if "protocol_claim_eligible" in row and str(
        row.get("protocol_claim_eligible", "")
    ).strip() != "":
        if not _truth(row, "protocol_claim_eligible"):
            reasons.append("protocol_claim_eligible=0")
    # Strict table: claim_ready required when column present
    if "claim_ready" in row and str(row.get("claim_ready", "")).strip() != "":
        if not _truth(row, "claim_ready"):
            reasons.append("claim_ready=0")
        else:
            ok_h, h_reasons = _hash_receipts_ok(row)
            if not ok_h:
                reasons.extend(h_reasons)
            if "tencom_status" in row and str(row.get("tencom_status", "")).strip() != "":
                if str(row.get("tencom_status", "")).strip().lower() != "ok":
                    reasons.append("tencom_status!=ok")
            if "eigen_status" in row and str(row.get("eigen_status", "")).strip() != "":
                if str(row.get("eigen_status", "")).strip().lower() != "ok":
                    reasons.append("eigen_status!=ok")
            if "pb_backend" in row and str(row.get("pb_backend", "")).strip() != "":
                if str(row.get("pb_backend", "")).strip() != "bust_cli":
                    reasons.append("pb_backend!=bust_cli")
    return (len(reasons) == 0, reasons)


def load_target_manifest() -> tuple[list[str], str] | tuple[None, None]:
    """Load the frozen pre-registered Astex-85 denominator manifest.

    Returns (sorted_upper_codes, sha256) or (None, None) if absent. When present,
    claim rates use the FIXED manifest count as denominator: a preregistered target
    that is absent from the campaign or dropped in admission counts as a FAILURE and
    is NEVER removed from the denominator. This is the P0 anti-inflation invariant.
    """
    here = Path(__file__).resolve().parent
    for cand in (
        here.parent / "benchmarks" / "protocols" / "astex85_target_manifest.json",
        here / "astex85_target_manifest.json",
    ):
        if cand.is_file():
            try:
                m = json.loads(cand.read_text())
                codes = [str(c).strip().upper() for c in m.get("targets", [])]
                if codes:
                    return sorted(set(codes)), str(m.get("sha256_of_sorted_codes", ""))
            except (ValueError, OSError):
                pass
    return None, None


def aggregate_rows(
    rows: Iterable[dict[str, str]],
    matrix_pin: str,
    matrix_pin_source: str,
    campaign_dir: str | None = None,
    fixed_denominator: bool = True,
) -> dict[str, Any]:
    all_rows = list(rows)
    claim: list[dict[str, str]] = []
    dropped: list[dict[str, Any]] = []

    for r in all_rows:
        ok, reasons = is_claim_eligible(r, matrix_pin)
        if ok:
            claim.append(r)
        else:
            dropped.append({"pdb_id": _pdb_id(r), "reasons": reasons})

    n = len(claim)

    # --- P0 fixed-denominator invariant ---------------------------------------
    # Claim rates must be reported over a FROZEN pre-registered target count, not
    # over the number of rows that happened to pass admission. Otherwise dropping
    # or losing hard targets mechanically inflates the rate. The denominator is
    # max(manifest N, observed distinct targets) so extra rows can never *shrink*
    # it below the preregistered set either.
    manifest_codes, manifest_sha = load_target_manifest()
    observed_ids = {_pdb_id(r).upper() for r in all_rows if _pdb_id(r)}
    if fixed_denominator and manifest_codes:
        denom = len(manifest_codes)
        missing_targets = sorted(set(manifest_codes) - observed_ids)
        denom_source = f"frozen_manifest(N={denom},sha={manifest_sha[:12]})"
        manifest_set: set[str] | None = {c.upper() for c in manifest_codes}
    else:
        denom = n
        missing_targets = []
        denom_source = "claim_eligible_rows(legacy)"
        manifest_set = None

    s1_ids: list[str] = []
    s2_ids: list[str] = []
    strict_ids: list[str] = []
    s3_ids: list[str] = []
    election_gap_ids: list[str] = []
    s1_fail_ids: list[str] = []
    n_legacy = 0

    for r in all_rows:
        if "claim_ready" not in r or str(r.get("claim_ready", "")).strip() == "":
            if _flag0(r, "seed_echo") and _flag0(r, "native_pose_seeded"):
                n_legacy += 1

    for r in claim:
        pid = _pdb_id(r)
        s1 = is_s1(r)
        s2 = is_s2(r, s1)
        s3 = is_s3(r)
        strict = is_strict_success(r)
        if s1:
            s1_ids.append(pid)
        else:
            s1_fail_ids.append(pid)
        if s2:
            s2_ids.append(pid)
        # STRICT numerator is ∩ the frozen 85-target manifest. Off-manifest
        # extras must not inflate n while the denominator stays 85.
        if strict and (
            manifest_set is None or pid.strip().upper() in manifest_set
        ):
            strict_ids.append(pid)
        if s3:
            s3_ids.append(pid)
        if s3 and not s1:
            election_gap_ids.append(pid)

    def rate(k: int) -> float:
        return (k / denom) if denom else 0.0

    # Targets whose RMSD verdict changes under the corrected metric. This is a
    # metric-only statement: it compares serial vs symmcorr on the SAME elected
    # pose and says nothing about admission. Direction is one-way by
    # construction (symmcorr <= serial), so this list can only ever contain
    # fail -> PASS moves; a pass -> fail entry would be a bug and is asserted
    # against in tests/test_rmsd_symmcorr.py.
    symmcorr_gained: list[str] = []
    for r in claim:
        rs = _f(r, *SYMMCORR_COLS)
        rc = _f(r, "rmsd_to_crystal")
        if not math.isfinite(rs):
            continue
        serial_pass = math.isfinite(rc) and 0.0 <= rc <= RMSD_SUCCESS_A
        if (0.0 <= rs <= RMSD_SUCCESS_A) and not serial_pass:
            symmcorr_gained.append(_pdb_id(r))

    # Which metric actually produced the S1 numbers on this table? An
    # unlabelled RMSD is not reportable (METHODOLOGY.md §0).
    _metrics_seen = sorted(
        {
            elected_rmsd_labelled(r)[1]
            for r in claim
            if math.isfinite(elected_rmsd_labelled(r)[0])
        }
    ) or ["none_finite"]
    s1_metric_used = "+".join(_metrics_seen)

    report: dict[str, Any] = {
        "contract": "admission_metrics_contract",
        "contract_doc": "benchmarks/protocols/admission_metrics_contract.md",
        "campaign_dir": campaign_dir,
        "matrix_md5_pin": matrix_pin,
        "matrix_md5_pin_source": matrix_pin_source,
        "N_raw": len(all_rows),
        "N_claim": n,
        "N_denominator": denom,
        "N_denominator_source": denom_source,
        "N_missing_from_manifest": len(missing_targets),
        "missing_targets": missing_targets,
        "N_dropped": len(dropped),
        "N_legacy_no_claim_ready": n_legacy,
        "dropped_rows": dropped,
        "metrics": {
            "S1": {
                "definition": (
                    "symmetry-corrected elected-pose RMSD <= 2.0 A "
                    "(spyrmsd graph automorphism, in-place; never hungarian). "
                    "Falls back to ordered rmsd_to_crystal when no symmcorr "
                    "sidecar is joined — conservative, since symmcorr <= serial."
                ),
                "metric_used": s1_metric_used,
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
                    "claim_ready==1 with PB + tENCoM/Eigen + hash receipts; "
                    "numerator ∩ frozen 85-target manifest"
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
        "symmcorr_delta": {
            "definition": (
                "targets whose elected-pose RMSD verdict moves fail->PASS when "
                "the serial metric is replaced by the symmetry-corrected one"
            ),
            "n": len(symmcorr_gained),
            "ids": symmcorr_gained,
            "direction_note": (
                "one-way by construction: symmetry correction minimises over "
                "graph automorphisms and the identity mapping is one of them, "
                "so symmcorr <= serial and no target can move PASS->fail"
            ),
        },
        "strict_metric_inheritance": {
            "note": (
                "STRICT is the engine's claim_ready column. In the engine "
                "claim_ready requires success_pb, and success_pb := "
                "success_rmsd AND pb_pass, where success_rmsd gates on the "
                "SERIAL rmsd_to_crystal (LIB/DatasetRunner.cpp). STRICT "
                "therefore still inherits the serial definition and is a "
                "CONSERVATIVE lower bound: it can under-count but never "
                "over-count. Repointing the engine gate is the remaining half "
                "of #365 and requires a rebuild."
            ),
            "s1_metric": s1_metric_used,
            "strict_metric": "serial (via engine success_pb)",
        },
        "headline": {
            "metric": "STRICT",
            "n": len(strict_ids),
            "N": denom,
            "rate": rate(len(strict_ids)),
            "label": "claim_ready strict success (rate over frozen 85-target denominator)",
        },
        "admission": {
            "seed_echo": 0,
            "native_pose_seeded": 0,
            "matrix_md5": matrix_pin,
            "claim_ready": 1,
        },
    }
    return report


def format_text_report(report: dict[str, Any]) -> str:
    m = report["metrics"]
    n = report.get("N_denominator", report["N_claim"])
    lines = [
        f"campaign_dir: {report.get('campaign_dir')}",
        f"matrix_md5_pin: {report['matrix_md5_pin']} (source={report['matrix_md5_pin_source']})",
        f"N_raw={report['N_raw']}  N_claim={report['N_claim']}  "
        f"N_denominator={n} ({report.get('N_denominator_source','')})  "
        f"N_dropped={report['N_dropped']}",
        (
            f"MISSING from frozen manifest (counted as failures): "
            f"{report['N_missing_from_manifest']} -> {', '.join(report['missing_targets'][:20])}"
            + (" ..." if report["N_missing_from_manifest"] > 20 else "")
        )
        if report.get("N_missing_from_manifest")
        else "all preregistered targets present",
        "",
        f"STRICT (headline): {m['STRICT']['n']}/{n} = {100.0 * m['STRICT']['rate']:.2f}%  "
        f"[claim_ready]",
        f"S1 (RMSD-only):    {m['S1']['n']}/{n} = {100.0 * m['S1']['rate']:.2f}%  "
        f"[diagnostic; metric={m['S1'].get('metric_used', 'unknown')}]",
        f"S2 (S1∧PB):        {m['S2']['n']}/{n} = {100.0 * m['S2']['rate']:.2f}%",
        f"S3 (pool ceiling): {m['S3']['n']}/{n} = {100.0 * m['S3']['rate']:.2f}%  "
        f"[diagnostic; conditional scanned pool — NOT any-pose]",
        f"election_gap (S3∧¬S1): {report['election_gap']['n']}/{n} = "
        f"{100.0 * report['election_gap']['rate']:.2f}%",
    ]
    _sd = report.get("symmcorr_delta")
    if _sd is not None:
        lines.append(
            f"symmcorr gain (fail->PASS): {_sd['n']}"
            + (f" -> {', '.join(_sd['ids'])}" if _sd["ids"] else "")
        )
    _sm = report.get("symmcorr")
    if _sm is not None and _sm.get("sidecar"):
        lines.append(
            f"symmcorr sidecar: {_sm['joined']}/{_sm['sidecar_rows']} rows joined"
        )
    elif report.get("metrics", {}).get("S1", {}).get("metric_used") == "serial":
        lines.append(
            "symmcorr sidecar: NONE — S1 used the serial metric "
            "(conservative; run scripts/rmsd_symmcorr.py and pass --symmcorr)"
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
    report["headline"] = {
        "metric": key,
        "n": m["n"],
        "N": report.get("N_denominator", report["N_claim"]),
        "rate": m["rate"],
        "label": m["definition"],
        "role": m["role"],
        "diagnostic_only": h in ("s1", "s2", "s3") or diagnostic_only,
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
    ap.add_argument(
        "--symmcorr",
        type=Path,
        default=None,
        help=(
            "Sidecar CSV from scripts/rmsd_symmcorr.py. Supplies the "
            "symmetry-corrected S1 metric (#365). Without it S1 falls back to "
            "the ordered serial column, which under-counts but never "
            "over-counts."
        ),
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
    symmcorr_provenance: dict[str, Any] = {
        "sidecar": None,
        "joined": 0,
        "refused_sha_mismatch": [],
        "sidecar_rows": 0,
    }
    if args.symmcorr is not None:
        if not args.symmcorr.is_file():
            print(f"not a file: {args.symmcorr}", file=sys.stderr)
            return 2
        sidecar = load_symmcorr_sidecar(args.symmcorr)
        symmcorr_provenance = join_symmcorr(rows, sidecar)
        symmcorr_provenance["sidecar"] = str(args.symmcorr)
        if symmcorr_provenance["refused_sha_mismatch"]:
            print(
                "error: symmcorr sidecar describes different poses for: "
                + ", ".join(symmcorr_provenance["refused_sha_mismatch"])
                + "\n(pose_sha256 mismatch — refusing to join a claim number "
                "to a pose it is not about)",
                file=sys.stderr,
            )
            return 2
    try:
        pin, pin_src = load_matrix_pin(campaign, args.matrix_md5)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    report = aggregate_rows(
        rows,
        matrix_pin=pin,
        matrix_pin_source=pin_src,
        campaign_dir=str(campaign.resolve()),
    )
    report["symmcorr"] = symmcorr_provenance
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

    return 0 if report["N_claim"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
