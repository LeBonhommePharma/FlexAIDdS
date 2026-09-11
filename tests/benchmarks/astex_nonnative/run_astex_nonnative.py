#!/usr/bin/env python3
"""run_astex_nonnative.py — Astex Non-Native cross-docking benchmark for FlexAIDdS.

Evaluates cross-docking accuracy: ligands are docked into receptor conformations
crystallized with different compounds. Reports success rates at 2 Å and 3 Å,
mean/median RMSD, and the FlexAIDdS entropy_rescue_rate — the fraction of pairs
where entropy-weighted re-ranking corrects for receptor conformational mismatch.

Usage:
    python run_astex_nonnative.py --results /path/to/results [--output report.json]
    python run_astex_nonnative.py --results /path/to/results --csv benchmarks/astex_nonnative/astex_non_native_set.csv

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
from typing import Dict, List, Optional


RMSD_THRESHOLD_2A = 2.0
RMSD_THRESHOLD_3A = 3.0


def load_nonnative_csv(csv_path: str, include_invalid: bool = False) -> List[Dict]:
    """Load cross-docking pairs from the canonical CSV.

    Returns list of dicts with keys:
      target_pdb, ligand_pdb, ligand_id, target_name, ligand_pdb_name,
      target_uniprot, ligand_uniprot, pair_status

    Schema changed 2026-09-11 (dataset identity audit, C2/C4):

    - ``rmsd_threshold_A`` REMOVED. It was 2.0 on all 74 rows -- a protocol
      parameter, not a property of the pair -- and this module never used the
      column anyway: the threshold it reports comes from ``--threshold``.
      A constant column invites a reader to treat the default as a measurement.
    - ``target_name`` / ``ligand_id`` were REPLACED, not edited. Both were
      fabricated: ``target_name`` said CDK2 for the 1hwi/1hww pair (RCSB: HMG-CoA
      reductase and alpha-mannosidase II), and the stated ``ligand_id`` matched
      neither member of its pair on 0/74 rows. They now hold the RCSB polymer
      description of ``target_pdb`` and the cognate CCD code of ``ligand_pdb``.
    - ``pair_status`` ADDED. A cross-docking pair must be two structures of the
      SAME protein. Measured by UniProt accession overlap, 66 of the 74 pairs
      are two DIFFERENT proteins and cannot be cross-docked at all. Those rows
      are REFUSED here rather than silently entering a success-rate denominator.

    Pass ``include_invalid=True`` only to inspect the refused rows.
    """
    pairs: List[Dict] = []
    refused: Dict[str, int] = {}
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            status = (row.get("pair_status") or "").strip()
            rec = {
                "target_pdb": row["target_pdb"].strip(),
                "ligand_pdb": row["ligand_pdb"].strip(),
                "ligand_id": row["ligand_id"].strip(),
                "target_name": row["target_name"].strip(),
                "ligand_pdb_name": (row.get("ligand_pdb_name") or "").strip(),
                "target_uniprot": (row.get("target_uniprot") or "").strip(),
                "ligand_uniprot": (row.get("ligand_uniprot") or "").strip(),
                "pair_status": status,
            }
            if status != "SAME_PROTEIN" and not include_invalid:
                refused[status or "MISSING_PAIR_STATUS"] = refused.get(status or "MISSING_PAIR_STATUS", 0) + 1
                continue
            pairs.append(rec)
    if refused:
        total = sum(refused.values())
        print(
            f"[REFUSED] {total} of {total + len(pairs)} rows in {csv_path} are not "
            f"cross-docking pairs and were excluded: "
            + ", ".join(f"{k}={v}" for k, v in sorted(refused.items()))
        )
        print(
            "          A cross-docking pair must be two structures of the same "
            "protein. Running them would put a protein-A receptor against a "
            "protein-B ligand and report the result as a cross-docking failure."
        )
    if not pairs:
        raise SystemExit(
            f"ERROR: no valid cross-docking pairs in {csv_path} "
            f"(all {sum(refused.values())} rows refused). Refusing to report a "
            f"success rate over an empty denominator."
        )
    return pairs


def pair_id(target_pdb: str, ligand_pdb: str) -> str:
    return f"{target_pdb}_x_{ligand_pdb}"


def load_docking_result(results_dir: str, target_pdb: str, ligand_pdb: str) -> Optional[Dict]:
    """Load FlexAIDdS docking result for one cross-docking pair.

    Looks for binding_modes.json under results_dir/{target_pdb}_x_{ligand_pdb}/.
    """
    pid = pair_id(target_pdb, ligand_pdb)
    result_path = Path(results_dir) / pid / "binding_modes.json"
    if not result_path.is_file():
        return None
    with open(result_path) as fh:
        return json.load(fh)


def extract_best_rmsd(result: Dict) -> Optional[float]:
    """Return RMSD of the top-ranked pose from a docking result dict."""
    modes = result.get("binding_modes", [])
    if not modes:
        return None
    return modes[0].get("best_pose_rmsd")


def extract_entropy_rescue(result: Dict) -> bool:
    """Return True if entropy re-ranking promoted a sub-2Å pose to rank 1."""
    modes = result.get("binding_modes", [])
    if not modes:
        return False
    return bool(modes[0].get("entropy_rescued", False))


def compute_crossdock_metrics(rmsds: List[float]) -> Dict:
    """Compute cross-docking success rates and RMSD statistics."""
    if not rmsds:
        return {
            "n": 0,
            "success_rate_2A": 0.0,
            "success_rate_3A": 0.0,
            "mean_rmsd": float("nan"),
            "median_rmsd": float("nan"),
        }
    n = len(rmsds)
    s = sorted(rmsds)
    mid = n // 2
    median = (s[mid] + s[mid - 1]) / 2.0 if n % 2 == 0 else s[mid]
    return {
        "n": n,
        "success_rate_2A": sum(1 for r in rmsds if r <= RMSD_THRESHOLD_2A) / n,
        "success_rate_3A": sum(1 for r in rmsds if r <= RMSD_THRESHOLD_3A) / n,
        "mean_rmsd": sum(rmsds) / n,
        "median_rmsd": median,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Astex Non-Native cross-docking benchmark for FlexAIDdS"
    )
    parser.add_argument("--results", required=True, help="FlexAIDdS results directory")
    parser.add_argument(
        "--csv",
        default=None,
        help="Path to astex_non_native_set.csv (default: benchmarks/astex_nonnative/astex_non_native_set.csv)",
    )
    parser.add_argument("--output", default="astex_nonnative_report.json", help="Output JSON report")
    args = parser.parse_args()

    # Resolve CSV path
    csv_path = args.csv
    if csv_path is None:
        repo_root = Path(__file__).parent.parent.parent.parent
        csv_path = str(repo_root / "benchmarks" / "astex_nonnative" / "astex_non_native_set.csv")

    if not Path(csv_path).is_file():
        print(f"ERROR: CSV not found: {csv_path}")
        sys.exit(1)

    pairs = load_nonnative_csv(csv_path)
    print(f"[Astex Non-Native] Loaded {len(pairs)} cross-docking pairs from CSV")
    print(f"[Astex Non-Native] Results directory: {args.results}")
    print()

    all_rmsds: List[float] = []
    n_entropy_rescued = 0
    n_missing = 0
    per_pair = []

    # Group by target family for stratified reporting
    by_family: Dict[str, List[float]] = {}

    for pair in pairs:
        target = pair["target_pdb"]
        lig_src = pair["ligand_pdb"]
        family = pair["target_name"]

        result = load_docking_result(args.results, target, lig_src)
        if result is None:
            n_missing += 1
            per_pair.append({
                "pair_id": pair_id(target, lig_src),
                "target_pdb": target,
                "ligand_pdb": lig_src,
                "ligand_id": pair["ligand_id"],
                "target_name": family,
                "status": "missing",
            })
            continue

        best = extract_best_rmsd(result)
        rescued = extract_entropy_rescue(result)

        if best is not None and not math.isnan(best):
            all_rmsds.append(best)
            by_family.setdefault(family, []).append(best)
        if rescued:
            n_entropy_rescued += 1

        per_pair.append({
            "pair_id": pair_id(target, lig_src),
            "target_pdb": target,
            "ligand_pdb": lig_src,
            "ligand_id": pair["ligand_id"],
            "target_name": family,
            "best_rmsd": best,
            "success_2A": best is not None and best <= RMSD_THRESHOLD_2A,
            "success_3A": best is not None and best <= RMSD_THRESHOLD_3A,
            "entropy_rescued": rescued,
            "status": "ok",
        })

    n_evaluated = len(all_rmsds)
    metrics = compute_crossdock_metrics(all_rmsds)
    entropy_rescue_rate = n_entropy_rescued / n_evaluated if n_evaluated else 0.0

    # Per-family breakdown
    family_metrics = {
        fam: compute_crossdock_metrics(rmsds)
        for fam, rmsds in sorted(by_family.items())
    }

    report = {
        "dataset": "Astex Non-Native Cross-Docking Set",
        "n_total_pairs": len(pairs),
        "n_evaluated": n_evaluated,
        "n_missing": n_missing,
        "overall_metrics": metrics,
        "entropy_rescue_rate": entropy_rescue_rate,
        "by_target_family": family_metrics,
        "per_pair": per_pair,
    }

    print(f"  N evaluated:               {n_evaluated} / {len(pairs)}")
    print(f"  Success rate (<2.0 Å):     {metrics['success_rate_2A']:.1%}")
    print(f"  Success rate (<3.0 Å):     {metrics['success_rate_3A']:.1%}")
    print(f"  Mean RMSD:                 {metrics['mean_rmsd']:.2f} Å")
    print(f"  Median RMSD:               {metrics['median_rmsd']:.2f} Å")
    print(f"  Entropy rescue rate:       {entropy_rescue_rate:.1%}")
    if n_missing:
        print(f"  Missing results:           {n_missing}")

    if family_metrics:
        print("\n  Per-family breakdown:")
        for fam, fm in family_metrics.items():
            print(f"    {fam:<25} n={fm['n']:3d}  {fm['success_rate_2A']:.0%} <2Å  mean {fm['mean_rmsd']:.2f} Å")

    with open(args.output, "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(f"\n[Astex Non-Native] Report written to {args.output}")


if __name__ == "__main__":
    main()
