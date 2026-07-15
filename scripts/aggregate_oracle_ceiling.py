#!/usr/bin/env python3
"""Aggregate oracle-ceiling rates from a DatasetRunner campaign output directory.

Pinned ceiling metric (see goal baseline):
  ceiling_success ⇔ best_cluster_rmsd ≤ 2.0 Å
  ceiling_rate = count(ceiling_success) / N

Also reports top-1 success_rmsd / success_pb when columns exist.

Usage:
  python3 scripts/aggregate_oracle_ceiling.py <campaign_dir> [--json out.json]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def _f(row: dict, *keys: str) -> float:
    for k in keys:
        v = row.get(k)
        if v is None or v == "" or v == "NA":
            continue
        try:
            x = float(v)
            if math.isfinite(x):
                return x
        except ValueError:
            continue
    return float("nan")


def _truth(row: dict, key: str) -> bool:
    return str(row.get(key, "")).strip() in ("1", "True", "true", "YES", "yes")


def load_campaign_rows(out_dir: Path) -> list[dict]:
    rows: list[dict] = []
    # Prefer per-target result.csv (authoritative for resume runs)
    for rc in sorted(out_dir.glob("*/result.csv")):
        try:
            batch = list(csv.DictReader(rc.open()))
            if batch:
                rows.append(batch[0])
        except OSError:
            continue
    if rows:
        return rows
    # Fallback: flat results CSV
    for name in (
        "astex_diverse_results.csv",
        "astex_crossdock_85_results.csv",
        "results.csv",
    ):
        p = out_dir / name
        if p.is_file():
            return list(csv.DictReader(p.open()))
    return rows


def aggregate(out_dir: Path) -> dict:
    rows = load_campaign_rows(out_dir)
    N = len(rows)
    ceiling_ids: list[str] = []
    top1_ids: list[str] = []
    pb_ids: list[str] = []
    suc_pb_ids: list[str] = []
    fails_bcr: list[str] = []
    fails_top1: list[str] = []

    for r in rows:
        pid = r.get("pdb_id") or r.get("pdb") or "?"
        bc = _f(r, "best_cluster_rmsd")
        ok_bcr = math.isfinite(bc) and 0.0 <= bc <= 2.0

        if "success_rmsd" in r:
            ok_top1 = _truth(r, "success_rmsd")
        elif "success" in r:
            ok_top1 = _truth(r, "success")
        else:
            rh = _f(r, "rmsd_hungarian", "rmsd_to_crystal")
            ok_top1 = math.isfinite(rh) and 0.0 <= rh <= 2.0
            if _truth(r, "seed_echo"):
                ok_top1 = False

        ok_pb = _truth(r, "pb_pass") if "pb_pass" in r else False
        if "success_pb" in r:
            ok_spb = _truth(r, "success_pb")
        else:
            ok_spb = ok_top1 and ok_pb

        if ok_bcr:
            ceiling_ids.append(pid)
        else:
            fails_bcr.append(pid)
        if ok_top1:
            top1_ids.append(pid)
        else:
            fails_top1.append(pid)
        if ok_pb:
            pb_ids.append(pid)
        if ok_spb:
            suc_pb_ids.append(pid)

    rate = (len(ceiling_ids) / N) if N else 0.0
    return {
        "out_dir": str(out_dir.resolve()),
        "N": N,
        "ceiling_metric": "best_cluster_rmsd <= 2.0 Angstrom (any-pose BCR)",
        "ceiling_n": len(ceiling_ids),
        "ceiling_rate": rate,
        "top1_success_rmsd_n": len(top1_ids),
        "top1_success_rmsd_rate": (len(top1_ids) / N) if N else 0.0,
        "pb_pass_n": len(pb_ids),
        "success_pb_n": len(suc_pb_ids),
        "success_pb_rate": (len(suc_pb_ids) / N) if N else 0.0,
        "failed_BCR_ids": fails_bcr,
        "failed_top1_ids": fails_top1,
        "exceeds_90_ceiling": rate > 0.90,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("campaign_dir", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    if not args.campaign_dir.is_dir():
        print(f"not a directory: {args.campaign_dir}", file=sys.stderr)
        return 2
    report = aggregate(args.campaign_dir)
    text = json.dumps(report, indent=2)
    print(text)
    if args.json:
        args.json.write_text(text + "\n")
    return 0 if report["N"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
