#!/usr/bin/env python3
"""run_astex_diverse.py — Astex Diverse Set benchmark evaluation for FlexAIDdS.

Evaluates native-pose docking accuracy on the 85-complex Astex Diverse Set.
Reports top-1/top-3 success rates (RMSD < 2 Å), mean/median RMSD, and the
FlexAIDdS-specific entropy_rescue_rate metric.

Usage:
    export ASTEX_DATA=/path/to/astex_diverse_structures
    python run_astex_diverse.py --results /path/to/flexaids_results [--output report.json]
    python run_astex_diverse.py --results /path/to/results --csv benchmarks/astex_diverse/astex_diverse_set.csv

Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
SPDX-License-Identifier: Apache-2.0
"""

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


RMSD_THRESHOLD = 2.0


def load_astex_diverse_csv(csv_path: str) -> List[Dict]:
    """Load Astex Diverse Set entries from the canonical CSV.

    Returns list of dicts with keys: pdb_id, ligand_id, resolution_A,
    rmsd_threshold_A, reference.
    """
    entries = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            entries.append({
                "pdb_id": row["pdb_id"].strip(),
                "ligand_id": row["ligand_id"].strip(),
                "resolution_A": float(row["resolution_A"]),
                "rmsd_threshold_A": float(row["rmsd_threshold_A"]),
                "reference": row["reference"].strip(),
            })
    return entries


def load_docking_result(results_dir: str, pdb_id: str) -> Optional[Dict]:
    """Load FlexAIDdS docking result for one complex.

    Looks for binding_modes.json under results_dir/pdb_id/.
    Returns the parsed dict or None if not found.
    """
    result_path = Path(results_dir) / pdb_id / "binding_modes.json"
    if not result_path.is_file():
        return None
    with open(result_path) as fh:
        return json.load(fh)


def extract_rmsd_list(result: Dict) -> List[float]:
    """Extract RMSD values for top binding modes from a docking result dict."""
    rmsds = []
    for mode in result.get("binding_modes", []):
        rmsd = mode.get("best_pose_rmsd")
        if rmsd is not None:
            rmsds.append(float(rmsd))
    return rmsds


def extract_entropy_rescue(result: Dict) -> bool:
    """Return True if entropy re-ranking promoted a sub-2Å pose to rank 1."""
    modes = result.get("binding_modes", [])
    if not modes:
        return False
    top = modes[0]
    return bool(top.get("entropy_rescued", False))


def compute_docking_power(
    rmsd_results: Dict[str, List[float]],
    threshold: float = RMSD_THRESHOLD,
) -> Dict[str, float]:
    """Compute top-1 and top-3 docking success rates."""
    n = len(rmsd_results)
    if n == 0:
        return {"n": 0, "top1": 0.0, "top3": 0.0}

    top1 = sum(
        1 for rmsds in rmsd_results.values()
        if rmsds and rmsds[0] <= threshold
    )
    top3 = sum(
        1 for rmsds in rmsd_results.values()
        if any(r <= threshold for r in rmsds[:3])
    )
    return {
        "n": n,
        "top1": top1 / n,
        "top3": top3 / n,
    }


def compute_rmsd_stats(all_rmsds: List[float]) -> Dict[str, float]:
    """Compute mean and median RMSD across all complexes (best pose per complex)."""
    if not all_rmsds:
        return {"mean": float("nan"), "median": float("nan")}
    s = sorted(all_rmsds)
    n = len(s)
    mean = sum(s) / n
    mid = n // 2
    median = (s[mid] + s[mid - 1]) / 2.0 if n % 2 == 0 else s[mid]
    return {"mean": mean, "median": median}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Astex Diverse Set benchmark evaluation for FlexAIDdS"
    )
    parser.add_argument("--results", required=True, help="FlexAIDdS results directory")
    parser.add_argument(
        "--csv",
        default=None,
        help="Path to astex_diverse_set.csv (default: benchmarks/astex_diverse/astex_diverse_set.csv)",
    )
    parser.add_argument("--output", default="astex_diverse_report.json", help="Output JSON report")
    parser.add_argument("--threshold", type=float, default=RMSD_THRESHOLD, help="RMSD threshold in Å")
    args = parser.parse_args()

    # Resolve CSV path
    csv_path = args.csv
    if csv_path is None:
        repo_root = Path(__file__).parent.parent.parent.parent
        csv_path = str(repo_root / "benchmarks" / "astex_diverse" / "astex_diverse_set.csv")

    if not Path(csv_path).is_file():
        print(f"ERROR: CSV not found: {csv_path}")
        sys.exit(1)

    entries = load_astex_diverse_csv(csv_path)
    print(f"[Astex Diverse] Loaded {len(entries)} entries from CSV")
    print(f"[Astex Diverse] Results directory: {args.results}")
    print(f"[Astex Diverse] RMSD threshold: {args.threshold} Å")
    print()

    rmsd_by_pdb: Dict[str, List[float]] = {}
    best_rmsds: List[float] = []
    n_entropy_rescued = 0
    n_missing = 0
    per_complex = []

    for entry in entries:
        pdb = entry["pdb_id"]
        result = load_docking_result(args.results, pdb)

        if result is None:
            n_missing += 1
            per_complex.append({"pdb_id": pdb, "ligand_id": entry["ligand_id"], "status": "missing"})
            continue

        rmsds = extract_rmsd_list(result)
        rescued = extract_entropy_rescue(result)

        rmsd_by_pdb[pdb] = rmsds
        best = rmsds[0] if rmsds else float("nan")
        if not math.isnan(best):
            best_rmsds.append(best)
        if rescued:
            n_entropy_rescued += 1

        per_complex.append({
            "pdb_id": pdb,
            "ligand_id": entry["ligand_id"],
            "resolution_A": entry["resolution_A"],
            "best_rmsd": best,
            "all_rmsds": rmsds[:5],
            "success": best <= args.threshold,
            "entropy_rescued": rescued,
            "status": "ok",
        })

    n_evaluated = len(rmsd_by_pdb)
    dp = compute_docking_power(rmsd_by_pdb, args.threshold)
    rmsd_stats = compute_rmsd_stats(best_rmsds)
    entropy_rescue_rate = n_entropy_rescued / n_evaluated if n_evaluated else 0.0

    report = {
        "dataset": "Astex Diverse Set",
        "n_total": len(entries),
        "n_evaluated": n_evaluated,
        "n_missing": n_missing,
        "rmsd_threshold_A": args.threshold,
        "docking_power": dp,
        "rmsd_stats": rmsd_stats,
        "entropy_rescue_rate": entropy_rescue_rate,
        "per_complex": per_complex,
    }

    print(f"  N evaluated:          {n_evaluated} / {len(entries)}")
    print(f"  Top-1 success (<{args.threshold}Å): {dp['top1']:.1%}")
    print(f"  Top-3 success (<{args.threshold}Å): {dp['top3']:.1%}")
    print(f"  Mean RMSD:            {rmsd_stats['mean']:.2f} Å")
    print(f"  Median RMSD:          {rmsd_stats['median']:.2f} Å")
    print(f"  Entropy rescue rate:  {entropy_rescue_rate:.1%}")
    if n_missing:
        print(f"  Missing results:      {n_missing}")

    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"\n[Astex Diverse] Report written to {args.output}")


if __name__ == "__main__":
    main()
