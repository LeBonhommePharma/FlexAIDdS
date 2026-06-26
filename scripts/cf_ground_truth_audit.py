#!/usr/bin/env python3
"""
cf_ground_truth_audit.py — verify CF(cluster) vs CF(native) on disk per target.

Cross-checks result.csv, stderr.log [NATIVE_CF], and REMARK CF= on emitted poses.
Assigns RMSD bands and harness red flags before any scoring re-tune.

Usage:
    python3 scripts/cf_ground_truth_audit.py <result_dir> [--lo-rmsd 2.0]

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

# Import shared helpers from failure_classify when available.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from failure_classify import (
        DEFAULT_CF_DELTA,
        DEFAULT_LO_RMSD,
        SEED_ECHO_CF_TOL,
        _float,
        _fget,
        _int,
        _load_per_target_csvs,
        _load_summary_csv,
        _seed_echo,
    )
except ImportError:
    DEFAULT_LO_RMSD = 2.0
    DEFAULT_CF_DELTA = 5.0
    SEED_ECHO_CF_TOL = 0.01

    def _fget(row, keys, default=None):
        for k in keys:
            v = row.get(k)
            if v is not None and str(v).strip() not in ("", "N/A", "nan"):
                return str(v).strip()
        return default

    def _float(row, keys, default=None):
        raw = _fget(row, keys)
        if raw is None:
            return default
        try:
            v = float(raw)
            return v if math.isfinite(v) else default
        except (ValueError, TypeError):
            return default

    def _int(row, keys, default=None):
        raw = _fget(row, keys)
        if raw is None:
            return default
        try:
            return int(float(raw))
        except (ValueError, TypeError):
            return default

    def _load_per_target_csvs(result_dir):
        rows = {}
        for csv_path in sorted(glob.glob(os.path.join(result_dir, "*/result.csv"))):
            pdb_id = os.path.basename(os.path.dirname(csv_path))
            try:
                with open(csv_path, newline="") as fh:
                    for row in csv.DictReader(fh):
                        rows[pdb_id] = row
                        break
            except OSError:
                pass
        return rows

    def _load_summary_csv(result_dir):
        candidates = sorted(
            glob.glob(os.path.join(result_dir, "*results*.csv"))
            + glob.glob(os.path.join(result_dir, "*summary*.csv"))
        )
        for csv_path in candidates:
            try:
                with open(csv_path, newline="") as fh:
                    rows = {}
                    for row in csv.DictReader(fh):
                        pid = (row.get("pdb_id") or "").strip()
                        if pid:
                            rows[pid] = row
                    if rows:
                        return rows
            except OSError:
                pass
        return {}

    def _seed_echo(row):
        return False

NATIVE_CF_RE = re.compile(
    r"\[NATIVE_CF\]\s+cf=([-\d.]+)(?:\s+breakdown=com:([-\d.]+),wal:([-\d.]+),"
    r"sas:([-\d.]+),con:([-\d.]+),hbond:([-\d.]+))?",
    re.IGNORECASE,
)
REMARK_CF_RE = re.compile(r"^REMARK\s+CF=([-\d.]+)", re.IGNORECASE)


def _rmsd_band(rmsd: float | None) -> str:
    if rmsd is None or rmsd < 0 or rmsd >= 998.0:
        return "invalid"
    if rmsd < 2.0:
        return "sub2"
    if rmsd < 6.0:
        return "near_miss_2_6"
    if rmsd < 17.0:
        return "deep_6_17"
    return "catastrophic_17plus"


def _parse_native_cf_stderr(target_dir: str) -> dict:
    out = {"cf_total": None, "breakdown": {}}
    log_path = os.path.join(target_dir, "stderr.log")
    if not os.path.isfile(log_path):
        return out
    try:
        with open(log_path) as fh:
            for line in fh:
                m = NATIVE_CF_RE.search(line)
                if m:
                    out["cf_total"] = float(m.group(1))
                    if m.group(2):
                        out["breakdown"] = {
                            "com": float(m.group(2)),
                            "wal": float(m.group(3)),
                            "sas": float(m.group(4)),
                            "con": float(m.group(5)),
                            "hbond": float(m.group(6)),
                        }
                    break
    except OSError:
        pass
    return out


def _parse_remark_cf(pdb_path: str) -> float | None:
    if not os.path.isfile(pdb_path):
        return None
    try:
        with open(pdb_path) as fh:
            for line in fh:
                m = REMARK_CF_RE.match(line.strip())
                if m:
                    return float(m.group(1))
    except OSError:
        pass
    return None


def _oracle_pose_cf(target_dir: str, cluster_idx: int | None) -> float | None:
    if cluster_idx is not None and cluster_idx >= 0:
        for name in (f"result_{cluster_idx}.pdb", f"result_{cluster_idx:02d}.pdb"):
            cf = _parse_remark_cf(os.path.join(target_dir, name))
            if cf is not None:
                return cf
    for pattern in ("result_*.pdb", "result_INI.pdb"):
        for path in sorted(glob.glob(os.path.join(target_dir, pattern))):
            cf = _parse_remark_cf(path)
            if cf is not None:
                return cf
    return None


def _harness_flags(row: dict, target_dir: str, rmsd: float | None) -> list[str]:
    flags = []
    count_delta = _int(row, ["count_delta"])
    if count_delta is not None and count_delta > 2:
        flags.append("v46_cache_count_delta")
    if rmsd is not None and rmsd >= 999.0:
        flags.append("rmsd_sentinel_999")
    cf_native = _float(row, ["cf_native"])
    best_score = _float(row, ["best_score", "best_cf"])
    if cf_native is not None and best_score is not None:
        if abs(best_score - cf_native) <= SEED_ECHO_CF_TOL and rmsd is not None and rmsd >= 6.0:
            flags.append("cf_match_high_rmsd_native_ic")
    if cf_native == 0.0 and best_score is not None and _int(row, ["num_poses"], 0) > 0:
        flags.append("v51_cf_native_zero")
    se_raw = _fget(row, ["seed_echo"])
    if se_raw in ("1", "true", "yes") and cf_native is not None and cf_native > 0:
        flags.append("seed_echo_clash_attractor")
    native_stderr = _parse_native_cf_stderr(target_dir)
    if cf_native is not None and native_stderr["cf_total"] is not None:
        if abs(cf_native - native_stderr["cf_total"]) > 1.0:
            flags.append("cf_native_stderr_mismatch")
    return flags


def _scoring_subtype(
    rmsd: float | None,
    oracle_rmsd: float | None,
    cf_best: float | None,
    cf_native: float | None,
    harness: list[str],
    lo_rmsd: float,
    cf_delta: float,
) -> str | None:
    if harness:
        return "harness_artifact"
    if rmsd is None or rmsd < 0 or rmsd < lo_rmsd:
        return None
    if oracle_rmsd is not None and oracle_rmsd < lo_rmsd and rmsd >= lo_rmsd:
        return "selection_miss"
    if cf_best is not None and cf_native is not None and cf_best < cf_native - cf_delta:
        if rmsd >= 6.0:
            return "CF_scoring_failure_deep"
        if 2.0 <= rmsd < 6.0:
            return "CF_scoring_failure_near"
        return "CF_false_minimum"
    if rmsd is not None and 2.0 <= rmsd < 6.0:
        if oracle_rmsd is not None and oracle_rmsd < lo_rmsd:
            return "CF_scoring_failure_near"
    return None


def audit_one(pdb_id: str, row: dict, result_dir: str, lo_rmsd: float, cf_delta: float) -> dict:
    target_dir = os.path.join(result_dir, pdb_id)
    rmsd = _float(row, ["rmsd_hungarian", "rmsd_to_crystal"], default=-1.0)
    oracle_rmsd = _float(row, ["best_cluster_rmsd"], default=rmsd)
    cf_native = _float(row, ["cf_native"])
    best_score = _float(row, ["best_score", "best_cf"])
    cluster_idx = _int(row, ["best_cluster_idx"])

    native_stderr = _parse_native_cf_stderr(target_dir) if os.path.isdir(target_dir) else {}
    oracle_cf = _oracle_pose_cf(target_dir, cluster_idx) if os.path.isdir(target_dir) else None

    cf_delta_sel = None
    cf_delta_oracle = None
    if cf_native is not None and best_score is not None:
        cf_delta_sel = round(best_score - cf_native, 4)
    if cf_native is not None and oracle_cf is not None:
        cf_delta_oracle = round(oracle_cf - cf_native, 4)

    harness = _harness_flags(row, target_dir, rmsd) if os.path.isdir(target_dir) else []
    band = _rmsd_band(rmsd)
    scoring_mode = _scoring_subtype(
        rmsd, oracle_rmsd, best_score, cf_native, harness, lo_rmsd, cf_delta
    )

    cf_proven_false_min = (
        scoring_mode in ("CF_scoring_failure_deep", "CF_false_minimum")
        and cf_delta_sel is not None
        and cf_delta_sel < -cf_delta
    )

    return {
        "pdb_id": pdb_id,
        "rmsd_hungarian": round(rmsd, 4) if rmsd >= 0 else rmsd,
        "best_cluster_rmsd": round(oracle_rmsd, 4) if oracle_rmsd is not None else None,
        "rmsd_band": band,
        "cf_native": cf_native,
        "best_score": best_score,
        "oracle_pose_cf": oracle_cf,
        "cf_delta_selected": cf_delta_sel,
        "cf_delta_oracle": cf_delta_oracle,
        "native_cf_stderr": native_stderr.get("cf_total"),
        "harness_flags": harness,
        "scoring_failure_mode": scoring_mode,
        "cf_false_minimum_proven": cf_proven_false_min,
        "success": rmsd >= 0 and rmsd < lo_rmsd,
    }


def build_markdown(report: dict, result_dir: str, lo_rmsd: float) -> str:
    targets = report["targets"]
    lines = [
        "# CF Ground-Truth Audit",
        "",
        f"**Result dir:** `{result_dir}`",
        f"**Targets:** {report['n_total']} | **Sub-2:** {report['n_sub2']} | "
        f"**Deep (6-17Å):** {report['n_deep']} | **Near-miss (2-6Å):** {report['n_near']}",
        "",
        "## Deep failures (6–17 Å)",
        "",
        "| PDB | RMSD | CF(native) | best_score | ΔCF | Oracle ΔCF | Mode | Harness |",
        "|-----|------|------------|------------|-----|------------|------|---------|",
    ]
    deep = [t for t in targets if t["rmsd_band"] == "deep_6_17"]
    def _fmt_cf(v):
        return f"{v:.2f}" if v is not None else "?"

    for t in sorted(deep, key=lambda x: x["rmsd_hungarian"]):
        lines.append(
            f"| {t['pdb_id']} | {t['rmsd_hungarian']:.2f} | "
            f"{_fmt_cf(t.get('cf_native'))} | {_fmt_cf(t.get('best_score'))} | "
            f"{t.get('cf_delta_selected', '?')} | {t.get('cf_delta_oracle', '?')} | "
            f"{t.get('scoring_failure_mode', '?')} | {','.join(t['harness_flags']) or '-'} |"
        )
    lines += [
        "",
        "## Near-misses (2–6 Å)",
        "",
        "| PDB | RMSD | Oracle | CF Δ | Mode |",
        "|-----|------|--------|------|------|",
    ]
    near = [t for t in targets if t["rmsd_band"] == "near_miss_2_6"]
    for t in sorted(near, key=lambda x: x["rmsd_hungarian"]):
        lines.append(
            f"| {t['pdb_id']} | {t['rmsd_hungarian']:.2f} | "
            f"{t.get('best_cluster_rmsd', '?'):.2f} | "
            f"{t.get('cf_delta_selected', '?')} | {t.get('scoring_failure_mode', '-')} |"
        )
    lines += [
        "",
        "## Harness suspects",
        "",
    ]
    harness = [t for t in targets if t["harness_flags"]]
    if harness:
        for t in harness:
            lines.append(f"- **{t['pdb_id']}**: {', '.join(t['harness_flags'])}")
    else:
        lines.append("_None detected._")
    lines.append("")
    return "\n".join(lines)


def audit_run(result_dir: str, lo_rmsd: float = DEFAULT_LO_RMSD, cf_delta: float = DEFAULT_CF_DELTA) -> dict:
    rows = _load_per_target_csvs(result_dir)
    if not rows:
        rows = _load_summary_csv(result_dir)
    if not rows:
        raise FileNotFoundError(f"No results found in {result_dir!r}")

    targets = [
        audit_one(pid, row, result_dir, lo_rmsd, cf_delta)
        for pid, row in sorted(rows.items())
    ]
    n_sub2 = sum(1 for t in targets if t["success"])
    n_deep = sum(1 for t in targets if t["rmsd_band"] == "deep_6_17")
    n_near = sum(1 for t in targets if t["rmsd_band"] == "near_miss_2_6")
    n_proven = sum(1 for t in targets if t["cf_false_minimum_proven"])
    n_harness = sum(1 for t in targets if t["harness_flags"])

    return {
        "result_dir": result_dir,
        "lo_rmsd": lo_rmsd,
        "cf_delta_threshold": cf_delta,
        "n_total": len(targets),
        "n_sub2": n_sub2,
        "n_deep": n_deep,
        "n_near": n_near,
        "n_cf_false_minimum_proven": n_proven,
        "n_harness_suspect": n_harness,
        "targets": targets,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description="CF ground-truth audit for benchmark runs.")
    parser.add_argument("result_dir", help="Benchmark output directory")
    parser.add_argument("--lo-rmsd", type=float, default=DEFAULT_LO_RMSD)
    parser.add_argument("--cf-delta", type=float, default=DEFAULT_CF_DELTA)
    args = parser.parse_args(argv)

    result_dir = os.path.expanduser(args.result_dir)
    if not os.path.isdir(result_dir):
        sys.exit(f"ERROR: not a directory: {result_dir}")

    report = audit_run(result_dir, lo_rmsd=args.lo_rmsd, cf_delta=args.cf_delta)

    json_path = os.path.join(result_dir, "cf_audit_report.json")
    with open(json_path, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    print(f"  → {json_path}", file=sys.stderr)

    md_path = os.path.join(result_dir, "cf_audit_summary.md")
    with open(md_path, "w") as fh:
        fh.write(build_markdown(report, result_dir, args.lo_rmsd))
    print(f"  → {md_path}", file=sys.stderr)

    print(
        f"Audit: {report['n_total']} targets | sub2={report['n_sub2']} | "
        f"deep={report['n_deep']} | near={report['n_near']} | "
        f"proven_false_min={report['n_cf_false_minimum_proven']} | "
        f"harness={report['n_harness_suspect']}"
    )


if __name__ == "__main__":
    main()