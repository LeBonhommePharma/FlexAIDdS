#!/usr/bin/env python3
"""
cf_ground_truth_audit.py — Phase-0 CF ground-truth audit for benchmark runs.

Verifies on disk whether failed targets are genuine CF scoring failures
(best/cluster CF beats native) or harness artifacts, before any scoring re-tune.

Usage:
    python3 scripts/cf_ground_truth_audit.py <result_dir> [--cf-delta 5.0]

Outputs:
    <result_dir>/cf_audit_report.json
    <result_dir>/cf_audit_summary.md

Copyright 2026 Le Bonhomme Pharma.  Apache-2.0.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import sys
from datetime import datetime, timezone

# Reuse failure taxonomy helpers from failure_classify.py
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from failure_classify import (  # noqa: E402
    DEFAULT_CF_DELTA,
    DEFAULT_LO_RMSD,
    SEED_ECHO_CF_TOL,
    _float,
    _fget,
    _harness_flags,
    _int,
    classify_one,
)

NATIVE_CF_RE = re.compile(
    r"\[NATIVE_CF\]\s+cf=([-\d.eE+]+)"
    r"(?:\s+breakdown=com:([-\d.eE+]+),wal:([-\d.eE+]+),"
    r"sas:([-\d.eE+]+),con:([-\d.eE+]+),hbond:([-\d.eE+]+))?"
)
REMARK_CF_RE = re.compile(r"REMARK CF=([-\d.eE+]+)")


def _rmsd_band(rmsd: float | None, lo_rmsd: float = DEFAULT_LO_RMSD) -> str:
    if rmsd is None or not math.isfinite(rmsd) or rmsd < 0:
        return "invalid"
    if rmsd < lo_rmsd:
        return "sub2"
    if rmsd < 2.0:
        return "sub2"  # defensive; lo_rmsd is usually 2.0
    if rmsd < 6.0:
        return "near_miss_2_6"
    if rmsd < 17.0:
        return "deep_6_17"
    return "catastrophic_17plus"


def _parse_native_cf(stderr_path: str) -> dict:
    out = {
        "cf_native_stderr": None,
        "native_breakdown": None,
        "native_cf_source": None,
    }
    if not os.path.isfile(stderr_path):
        return out
    try:
        with open(stderr_path) as fh:
            for line in fh:
                m = NATIVE_CF_RE.search(line)
                if not m:
                    continue
                out["cf_native_stderr"] = float(m.group(1))
                if m.group(2) is not None:
                    out["native_breakdown"] = {
                        "com": float(m.group(2)),
                        "wal": float(m.group(3)),
                        "sas": float(m.group(4)),
                        "con": float(m.group(5)),
                        "hbond": float(m.group(6)),
                    }
                out["native_cf_source"] = "stderr.log"
                break
    except OSError:
        pass
    return out


def _read_remark_cf(pdb_path: str) -> float | None:
    if not os.path.isfile(pdb_path):
        return None
    try:
        with open(pdb_path) as fh:
            for line in fh:
                m = REMARK_CF_RE.match(line.strip())
                if m:
                    v = float(m.group(1))
                    return v if math.isfinite(v) else None
    except OSError:
        pass
    return None


def _pose_pdb_candidates(target_dir: str, pdb_id: str, idx: int) -> list[str]:
    """Return candidate pose PDB paths for cluster/pose index `idx`."""
    pid = pdb_id.upper()
    pid_lo = pdb_id.lower()
    names = [
        f"{pid}_{idx}.pdb",
        f"{pid_lo}_{idx}.pdb",
        f"result_{idx}.pdb",
        os.path.join("r1", f"{pid}_{idx}.pdb"),
        os.path.join("r1", f"result_{idx}.pdb"),
    ]
    return [os.path.join(target_dir, n) for n in names]


def _discover_pose_pdbs(target_dir: str, pdb_id: str) -> list[tuple[int, str, float]]:
    """Return [(idx, path, cf), ...] for all emitted pose PDBs with REMARK CF."""
    pid = pdb_id.upper()
    pid_lo = pdb_id.lower()
    patterns = [
        os.path.join(target_dir, f"{pid}_*.pdb"),
        os.path.join(target_dir, f"{pid_lo}_*.pdb"),
        os.path.join(target_dir, "result_*.pdb"),
    ]
    poses: list[tuple[int, str, float]] = []
    seen_paths: set[str] = set()
    idx_re = re.compile(rf"(?:{re.escape(pid)}|{re.escape(pid_lo)}|result)_(\d+)\.pdb$", re.I)

    for pattern in patterns:
        for path in sorted(glob.glob(pattern)):
            if path in seen_paths:
                continue
            base = os.path.basename(path)
            if base.endswith("_INI.pdb"):
                continue
            m = idx_re.search(base)
            if not m:
                continue
            cf = _read_remark_cf(path)
            if cf is None:
                continue
            seen_paths.add(path)
            poses.append((int(m.group(1)), path, cf))
    return poses


def _audit_one(pdb_id: str, row: dict, target_dir: str, lo_rmsd: float, cf_delta: float) -> dict:
    rmsd = _float(row, ["rmsd_hungarian", "rmsd_to_crystal"], default=-1.0)
    oracle_rmsd = _float(row, ["best_cluster_rmsd"], default=rmsd)
    cf_native_csv = _float(row, ["cf_native"])
    best_score = _float(row, ["best_score", "best_cf"])
    best_cluster_idx = _int(row, ["best_cluster_idx"], default=0)

    stderr_path = os.path.join(target_dir, "stderr.log") if target_dir else ""
    native_info = _parse_native_cf(stderr_path)
    cf_native = cf_native_csv
    cf_native_mismatch = None
    if cf_native is not None and native_info["cf_native_stderr"] is not None:
        cf_native_mismatch = round(cf_native - native_info["cf_native_stderr"], 6)

    poses = _discover_pose_pdbs(target_dir, pdb_id) if target_dir else []
    pose_by_idx = {idx: (path, cf) for idx, path, cf in poses}
    min_pose_cf = min((cf for _, _, cf in poses), default=None)
    min_pose_idx = None
    if min_pose_cf is not None:
        for idx, _, cf in poses:
            if cf == min_pose_cf:
                min_pose_idx = idx
                break

    oracle_pose_cf = None
    oracle_pose_path = None
    if best_cluster_idx is not None:
        for cand in _pose_pdb_candidates(target_dir or "", pdb_id, best_cluster_idx):
            cf = _read_remark_cf(cand)
            if cf is not None:
                oracle_pose_cf = cf
                oracle_pose_path = cand
                break
        if oracle_pose_cf is None and best_cluster_idx in pose_by_idx:
            oracle_pose_path, oracle_pose_cf = pose_by_idx[best_cluster_idx]

    selected_pose_cf = best_score
    if selected_pose_cf is None and poses:
        # Fallback: rank-0 / lowest CF pose in ensemble
        selected_pose_cf = min_pose_cf

    cf_delta_selected = None
    cf_delta_oracle = None
    cf_delta_min_pose = None
    if cf_native is not None:
        if selected_pose_cf is not None:
            cf_delta_selected = round(selected_pose_cf - cf_native, 4)
        if oracle_pose_cf is not None:
            cf_delta_oracle = round(oracle_pose_cf - cf_native, 4)
        if min_pose_cf is not None:
            cf_delta_min_pose = round(min_pose_cf - cf_native, 4)

    scoring_proof = (
        cf_native is not None
        and selected_pose_cf is not None
        and selected_pose_cf < cf_native - cf_delta
    )
    oracle_scoring_proof = (
        cf_native is not None
        and oracle_pose_cf is not None
        and oracle_pose_cf < cf_native - cf_delta
    )

    harness = _harness_flags(row, target_dir, rmsd) if target_dir else []
    if cf_native_mismatch is not None and abs(cf_native_mismatch) > 0.05:
        harness.append("cf_native_csv_stderr_mismatch")
    # Oracle pose is near-native geometry but CF diverges from native_score path.
    if (
        oracle_rmsd is not None
        and oracle_rmsd < 0.5
        and cf_native is not None
        and oracle_pose_cf is not None
        and abs(oracle_pose_cf - cf_native) > 10.0
    ):
        harness.append("oracle_native_cf_divergence")

    failure = classify_one(pdb_id, row, target_dir, lo_rmsd, cf_delta)

    return {
        "pdb_id": pdb_id,
        "success": failure["success"],
        "rmsd_hungarian": round(rmsd, 4) if rmsd is not None else None,
        "best_cluster_rmsd": round(oracle_rmsd, 4) if oracle_rmsd is not None else None,
        "rmsd_band": _rmsd_band(rmsd, lo_rmsd),
        "failure_mode": failure.get("failure_mode"),
        "cf_native": cf_native,
        "cf_native_stderr": native_info["cf_native_stderr"],
        "cf_native_mismatch": cf_native_mismatch,
        "native_breakdown": native_info["native_breakdown"],
        "best_score": best_score,
        "selected_pose_cf": selected_pose_cf,
        "oracle_pose_cf": oracle_pose_cf,
        "oracle_pose_path": oracle_pose_path,
        "min_pose_cf": min_pose_cf,
        "min_pose_idx": min_pose_idx,
        "num_poses_with_cf": len(poses),
        "cf_delta_selected": cf_delta_selected,
        "cf_delta_oracle": cf_delta_oracle,
        "cf_delta_min_pose": cf_delta_min_pose,
        "cf_scoring_proof_selected": scoring_proof,
        "cf_scoring_proof_oracle": oracle_scoring_proof,
        "harness_flags": harness,
        "best_cluster_idx": best_cluster_idx,
        "seed_echo": _fget(row, ["seed_echo"]),
    }


def _load_rows(result_dir: str) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for csv_path in sorted(glob.glob(os.path.join(result_dir, "*/result.csv"))):
        pdb_id = os.path.basename(os.path.dirname(csv_path))
        try:
            with open(csv_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    rows[pdb_id] = row
                    break
        except OSError:
            pass
    if rows:
        return rows

    for pattern in (
        os.path.join(result_dir, "*results*.csv"),
        os.path.join(result_dir, "*summary*.csv"),
    ):
        for csv_path in sorted(glob.glob(pattern)):
            try:
                with open(csv_path, newline="") as fh:
                    reader = csv.DictReader(fh)
                    for row in reader:
                        pid = (row.get("pdb_id") or "").strip()
                        if pid:
                            rows[pid] = row
            except OSError:
                pass
        if rows:
            return rows
    return rows


def _load_provenance(result_dir: str) -> dict | None:
    path = os.path.join(result_dir, "provenance.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def _cohort_summary(targets: list[dict], band: str) -> dict:
    cohort = [t for t in targets if t["rmsd_band"] == band and not t["success"]]
    scoring_proven = [
        t for t in cohort
        if t["cf_scoring_proof_selected"] and not t["harness_flags"]
    ]
    harness = [t for t in cohort if t["harness_flags"]]
    return {
        "band": band,
        "count": len(cohort),
        "cf_scoring_proven": len(scoring_proven),
        "harness_suspect": len(harness),
        "targets": [t["pdb_id"] for t in cohort],
        "scoring_proven_ids": [t["pdb_id"] for t in scoring_proven],
        "harness_ids": [t["pdb_id"] for t in harness],
    }


def _build_markdown(report: dict) -> str:
    meta = report["meta"]
    targets = report["targets"]
    lines = [
        "# CF Ground-Truth Audit",
        "",
        f"**Result dir:** `{meta['result_dir']}`  ",
        f"**Generated:** {meta['generated_at']}  ",
        f"**Targets audited:** {meta['n_targets']} / 85 expected  ",
    ]
    if meta.get("git_commit"):
        lines.append(f"**Git commit:** `{meta['git_commit']}`  ")
    if meta.get("coverage_note"):
        lines.append(f"**Coverage:** {meta['coverage_note']}  ")
    lines.append("")

    n_success = sum(1 for t in targets if t["success"])
    n_fail = len(targets) - n_success
    lines += [
        "## Summary",
        "",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Success (RMSD < {meta['lo_rmsd']:.1f} Å) | {n_success}/{len(targets)} |",
        f"| Failed | {n_fail} |",
        f"| CF scoring proven (selected, Δ < -{meta['cf_delta']:.0f}) | "
        f"{sum(1 for t in targets if t['cf_scoring_proof_selected'] and not t['success'])} |",
        f"| Harness suspect | "
        f"{sum(1 for t in targets if t['harness_flags'] and not t['success'])} |",
        "",
    ]

    for band in ("deep_6_17", "near_miss_2_6", "catastrophic_17plus"):
        c = report["cohorts"].get(band, {})
        if c.get("count", 0) == 0:
            continue
        lines += [
            f"## Cohort: {band}",
            "",
            f"- Targets ({c['count']}): {', '.join(c['targets'])}",
            f"- CF scoring proven: {c['cf_scoring_proven']} — "
            f"{', '.join(c['scoring_proven_ids']) or '—'}",
            f"- Harness suspect: {c['harness_suspect']} — "
            f"{', '.join(c['harness_ids']) or '—'}",
            "",
        ]

    lines += [
        "## Per-target detail (failures)",
        "",
        "| PDB | RMSD | Band | Mode | cf_native | best_score | Δ_sel | "
        "oracle_CF | Δ_orc | Proof | Harness |",
        "|-----|------|------|------|-----------|------------|-------|"
        "-----------|-------|-------|---------|",
    ]
    for t in sorted(
        (x for x in targets if not x["success"]),
        key=lambda x: x.get("rmsd_hungarian") or 999,
    ):
        proof = "Y" if t["cf_scoring_proof_selected"] else "N"
        harness = ",".join(t["harness_flags"]) if t["harness_flags"] else "—"
        lines.append(
            f"| {t['pdb_id']} | {t['rmsd_hungarian']:.2f} | {t['rmsd_band']} | "
            f"{t.get('failure_mode') or '?'} | "
            f"{t['cf_native']:.2f} | {t['best_score']:.2f} | "
            f"{t['cf_delta_selected']:.2f} | "
            f"{t['oracle_pose_cf'] if t['oracle_pose_cf'] is not None else '—'} | "
            f"{t['cf_delta_oracle'] if t['cf_delta_oracle'] is not None else '—'} | "
            f"{proof} | {harness} |"
        )

    lines += [
        "",
        "## Phase-0 exit criterion",
        "",
    ]
    deep = report["cohorts"].get("deep_6_17", {})
    deep_unresolved = [
        t["pdb_id"] for t in targets
        if t["rmsd_band"] == "deep_6_17"
        and not t["success"]
        and not t["cf_scoring_proof_selected"]
        and not t["harness_flags"]
    ]
    if not deep.get("count"):
        lines.append(
            "No deep (6–17 Å) failures in this partial/full run — "
            "re-run full 85-target baseline to complete deep-cohort audit."
        )
    elif not deep_unresolved:
        lines.append(
            "PASS: Every deep failure has CF proof on disk or harness reclassification."
        )
    else:
        lines.append(
            "BLOCKED: Unresolved deep failures (no CF proof, not harness): "
            + ", ".join(deep_unresolved)
        )

    return "\n".join(lines) + "\n"


def audit_run(
    result_dir: str,
    lo_rmsd: float = DEFAULT_LO_RMSD,
    cf_delta: float = DEFAULT_CF_DELTA,
) -> dict:
    result_dir = os.path.expanduser(result_dir)
    rows = _load_rows(result_dir)
    if not rows:
        raise FileNotFoundError(f"No result data found in {result_dir!r}")

    provenance = _load_provenance(result_dir)
    targets = []
    for pdb_id, row in sorted(rows.items()):
        target_dir = os.path.join(result_dir, pdb_id)
        if not os.path.isdir(target_dir):
            target_dir = None
        targets.append(_audit_one(pdb_id, row, target_dir, lo_rmsd, cf_delta))

    cohorts = {
        band: _cohort_summary(targets, band)
        for band in ("near_miss_2_6", "deep_6_17", "catastrophic_17plus")
    }

    n_targets = len(targets)
    coverage_note = None
    if n_targets < 85:
        coverage_note = (
            f"PARTIAL RUN — {n_targets}/85 targets on disk; "
            "deep-cohort conclusions apply only to completed targets."
        )

    return {
        "meta": {
            "result_dir": result_dir,
            "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "lo_rmsd": lo_rmsd,
            "cf_delta": cf_delta,
            "seed_echo_cf_tol": SEED_ECHO_CF_TOL,
            "n_targets": n_targets,
            "git_commit": (provenance or {}).get("git_commit"),
            "binary_sha256": (provenance or {}).get("binary_sha256"),
            "coverage_note": coverage_note,
        },
        "cohorts": cohorts,
        "targets": targets,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="CF ground-truth audit for FlexAIDdS benchmark runs."
    )
    parser.add_argument("result_dir", help="Benchmark output root directory")
    parser.add_argument(
        "--lo-rmsd", type=float, default=DEFAULT_LO_RMSD,
        help=f"Success threshold in Å (default {DEFAULT_LO_RMSD})",
    )
    parser.add_argument(
        "--cf-delta", type=float, default=DEFAULT_CF_DELTA,
        help=f"CF gap for scoring-failure proof (default {DEFAULT_CF_DELTA})",
    )
    parser.add_argument(
        "--also-classify", action="store_true",
        help="Also run failure_classify.py on the same directory",
    )
    args = parser.parse_args(argv)

    result_dir = os.path.expanduser(args.result_dir)
    if not os.path.isdir(result_dir):
        print(f"ERROR: not a directory: {result_dir}", file=sys.stderr)
        return 1

    report = audit_run(result_dir, lo_rmsd=args.lo_rmsd, cf_delta=args.cf_delta)

    json_path = os.path.join(result_dir, "cf_audit_report.json")
    with open(json_path, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    print(f"  → {json_path}", file=sys.stderr)

    md_path = os.path.join(result_dir, "cf_audit_summary.md")
    with open(md_path, "w") as fh:
        fh.write(_build_markdown(report))
    print(f"  → {md_path}", file=sys.stderr)

    if args.also_classify:
        from failure_classify import classify_run, build_summary  # noqa: E402
        results, missing = classify_run(result_dir, lo_rmsd=args.lo_rmsd, cf_delta=args.cf_delta)
        summary = build_summary(results, args.lo_rmsd, 2.5, missing)
        print("\n".join(
            l.replace("(see JSON)", result_dir) if "(see JSON)" in l else l
            for l in summary
        ))

    # Exit 1 if deep failures lack proof and aren't harness
    deep_blocked = [
        t["pdb_id"] for t in report["targets"]
        if t["rmsd_band"] == "deep_6_17"
        and not t["success"]
        and not t["cf_scoring_proof_selected"]
        and not t["harness_flags"]
    ]
    if deep_blocked:
        print(
            f"WARNING: {len(deep_blocked)} deep failure(s) lack CF proof: "
            + ", ".join(deep_blocked),
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())