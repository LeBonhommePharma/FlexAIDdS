#!/usr/bin/env python3
"""10k bootstrap median of S_top10 success rates (3Dsig 2017 metric).

Per case: success if min(RMSD among top-10 ranked modes) < 2.0 Å.
Headline: median of B bootstrap resamples of the case set (with replacement).

Input formats (auto-detected):
  - Directory of per-PDB result.csv with columns like top_k_rmsd / mode_rmsd_0..9
  - TSV/CSV with columns: pdb_id, rank, rmsd  (rank 0..9)
  - JSON list of {pdb_id, rmsds: [..]} 

Usage:
  python3 scripts/bootstrap_3dsig_s_top10.py --cases cases.json --bootstraps 10000
  python3 scripts/bootstrap_3dsig_s_top10.py --arm-dir path/to/3dsig_r10

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional


def s_top10(rmsds: List[float], thresh: float = 2.0) -> bool:
    vals = [r for r in rmsds[:10] if r is not None and math.isfinite(r) and r >= 0.0]
    if not vals:
        return False
    return min(vals) < thresh


def load_cases_json(path: Path) -> Dict[str, List[float]]:
    data = json.loads(path.read_text())
    out: Dict[str, List[float]] = {}
    if isinstance(data, dict) and "cases" in data:
        data = data["cases"]
    if isinstance(data, dict):
        for k, v in data.items():
            if isinstance(v, dict) and "rmsds" in v:
                out[str(k)] = [float(x) for x in v["rmsds"]]
            elif isinstance(v, list):
                out[str(k)] = [float(x) for x in v]
    elif isinstance(data, list):
        for row in data:
            pid = str(row.get("pdb_id") or row.get("pdb") or row["id"])
            out[pid] = [float(x) for x in row["rmsds"]]
    return out


def load_rank_table(path: Path) -> Dict[str, List[float]]:
    """pdb_id, rank, rmsd rows → top-10 rmsds sorted by rank."""
    by: Dict[str, List[tuple]] = {}
    with path.open(newline="") as f:
        # sniff delimiter
        sample = f.read(2048)
        f.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        r = csv.DictReader(f, delimiter=delim)
        for row in r:
            keys = {k.lower(): k for k in row}
            pid = row[keys.get("pdb_id") or keys.get("pdb") or keys.get("case") or list(row)[0]]
            rank = int(float(row[keys.get("rank") or list(row)[1]]))
            rmsd = float(row[keys.get("rmsd") or list(row)[2]])
            by.setdefault(str(pid), []).append((rank, rmsd))
    out: Dict[str, List[float]] = {}
    for pid, pairs in by.items():
        pairs.sort(key=lambda x: x[0])
        out[pid] = [r for _, r in pairs[:10]]
    return out


def load_arm_dir(arm_dir: Path) -> Dict[str, List[float]]:
    """Scan result.csv under arm_dir/<PDB>/ for rmsd columns if present."""
    out: Dict[str, List[float]] = {}
    for csv_path in sorted(arm_dir.rglob("result.csv")):
        try:
            with csv_path.open(newline="") as f:
                rows = list(csv.DictReader(f))
        except OSError:
            continue
        if not rows:
            continue
        row = rows[0]
        pdb = (
            row.get("pdb_id")
            or row.get("receptor_id")
            or csv_path.parent.name
        ).upper()[:4]
        rmsds: List[float] = []
        for i in range(10):
            for key in (f"mode_rmsd_{i}", f"rmsd_{i}", f"top{i}_rmsd", f"rank{i}_rmsd"):
                if key in row and row[key] not in ("", None, "NA"):
                    try:
                        v = float(row[key])
                        if math.isfinite(v) and v >= 0:
                            rmsds.append(v)
                            break
                    except ValueError:
                        pass
        # fallback: single elected + best_cluster only (incomplete for S_top10)
        if not rmsds:
            for key in ("rmsd_to_crystal", "rmsd_hungarian", "best_cluster_rmsd"):
                if key in row and row[key] not in ("", None, "NA"):
                    try:
                        v = float(row[key])
                        if math.isfinite(v) and v >= 0:
                            rmsds.append(v)
                    except ValueError:
                        pass
        if rmsds:
            out[pdb] = rmsds
    return out


def bootstrap_median(
    case_success: List[bool], n_boot: int, seed: int
) -> dict:
    rng = random.Random(seed)
    n = len(case_success)
    if n == 0:
        return {"n_cases": 0, "point": None, "median": None, "p05": None, "p95": None}
    point = sum(case_success) / n
    samples = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        samples.append(sum(case_success[i] for i in idx) / n)
    samples.sort()
    med = statistics.median(samples)
    p05 = samples[max(0, int(0.05 * n_boot) - 1)]
    p95 = samples[min(n_boot - 1, int(0.95 * n_boot))]
    return {
        "n_cases": n,
        "n_success": sum(case_success),
        "point": point,
        "median": med,
        "p05": p05,
        "p95": p95,
        "n_boot": n_boot,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=Path, help="JSON cases file")
    ap.add_argument("--rank-table", type=Path, help="CSV/TSV pdb,rank,rmsd")
    ap.add_argument("--arm-dir", type=Path, help="Campaign arm directory with result.csv")
    ap.add_argument("--bootstraps", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=20170715)
    ap.add_argument("--thresh", type=float, default=2.0)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args()

    cases: Dict[str, List[float]] = {}
    if args.cases:
        cases = load_cases_json(args.cases)
    elif args.rank_table:
        cases = load_rank_table(args.rank_table)
    elif args.arm_dir:
        cases = load_arm_dir(args.arm_dir)
    else:
        print("Provide --cases, --rank-table, or --arm-dir", file=sys.stderr)
        return 2

    success = [s_top10(v, args.thresh) for v in cases.values()]
    stats = bootstrap_median(success, args.bootstraps, args.seed)
    stats["metric"] = "S_top10"
    stats["thresh_A"] = args.thresh
    stats["n_pdbs_loaded"] = len(cases)
    stats["pdb_ids"] = sorted(cases.keys())

    print(
        f"S_top10 point={stats['point']:.4f}  "
        f"bootstrap_median={stats['median']:.4f}  "
        f"p05–p95=[{stats['p05']:.4f},{stats['p95']:.4f}]  "
        f"n={stats['n_cases']} success={stats['n_success']}"
    )
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(stats, indent=2) + "\n")
        print("wrote", args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
