#!/usr/bin/env python3
"""Aggregate DP vs FO small-sim pilot metrics.

Metrics per target / arm:
  - S1: elected pose RMSD <= 2.0 Å (and not seed echo when columns exist)
  - BCR: any cluster-rep / pose RMSD <= 2.0 Å among top results
  - election_gap: BCR success but S1 fail
  - predictive_power (arm-level): among targets with BCR success, fraction with S1 success
    = P(elect native-like | native-like exists in cluster output)

Usage:
  python3 scripts/aggregate_dpfo_pilot.py --base-out <DPFO_pilot_base>
  python3 scripts/aggregate_dpfo_pilot.py --base-out <path> --json-out summary.json

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


def _f(row: dict, *keys: str) -> Optional[float]:
    for k in keys:
        if k in row and row[k] not in (None, "", "NA", "nan"):
            try:
                v = float(row[k])
                if math.isfinite(v):
                    return v
            except (TypeError, ValueError):
                pass
    return None


def _truth(row: dict, key: str) -> bool:
    v = str(row.get(key, "")).strip().lower()
    return v in ("1", "true", "yes", "t")


def load_result_csvs(arm_dir: Path) -> List[dict]:
    rows: List[dict] = []
    if not arm_dir.is_dir():
        return rows
    for csv_path in sorted(arm_dir.rglob("result.csv")):
        try:
            with csv_path.open(newline="") as f:
                r = csv.DictReader(f)
                for row in r:
                    row["_path"] = str(csv_path)
                    rows.append(row)
        except OSError:
            continue
    # also summary-level CSVs
    for name in ("summary.csv", "results.csv", "benchmark_summary.csv"):
        p = arm_dir / name
        if p.is_file():
            try:
                with p.open(newline="") as f:
                    for row in csv.DictReader(f):
                        row["_path"] = str(p)
                        rows.append(row)
            except OSError:
                pass
    return rows


def pdb_id(row: dict) -> str:
    for k in ("pdb_id", "receptor_id", "complex_id", "system", "name", "id"):
        if row.get(k):
            return str(row[k]).strip().upper()[:4]
    # from path
    p = Path(row.get("_path", ""))
    for part in p.parts[::-1]:
        if len(part) == 4 and part[0].isdigit():
            return part.upper()
    return "UNKN"


def score_row(row: dict) -> Dict[str, Any]:
    rmsd_s1 = _f(
        row,
        "rmsd",
        "elected_rmsd",
        "best_rmsd",
        "top_rmsd",
        "rmsd_top1",
        "S1_rmsd",
    )
    rmsd_bcr = _f(
        row,
        "best_cluster_rmsd",
        "bcr_rmsd",
        "min_rmsd",
        "rmsd_min",
        "lowest_rmsd",
        "best_pose_rmsd",
    )
    if rmsd_bcr is None:
        rmsd_bcr = rmsd_s1
    seed_echo = _truth(row, "seed_echo") if "seed_echo" in row else False
    s1 = (rmsd_s1 is not None and rmsd_s1 <= 2.0 and not seed_echo)
    bcr = (rmsd_bcr is not None and rmsd_bcr <= 2.0)
    return {
        "pdb_id": pdb_id(row),
        "rmsd_s1": rmsd_s1,
        "rmsd_bcr": rmsd_bcr,
        "s1": s1,
        "bcr": bcr,
        "election_gap": bool(bcr and not s1),
        "seed_echo": seed_echo,
        "path": row.get("_path"),
    }


def summarize_arm(arm: str, arm_dir: Path) -> Dict[str, Any]:
    rows = load_result_csvs(arm_dir)
    by_pdb: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        sc = score_row(row)
        pid = sc["pdb_id"]
        prev = by_pdb.get(pid)
        if prev is None:
            by_pdb[pid] = sc
        else:
            # keep best (lowest) S1 rmsd as elected if multiple
            a = prev.get("rmsd_s1")
            b = sc.get("rmsd_s1")
            if b is not None and (a is None or b < a):
                by_pdb[pid] = sc

    targets = sorted(by_pdb.values(), key=lambda x: x["pdb_id"])
    n = len(targets)
    n_s1 = sum(1 for t in targets if t["s1"])
    n_bcr = sum(1 for t in targets if t["bcr"])
    n_gap = sum(1 for t in targets if t["election_gap"])
    # predictive power: P(S1 | BCR)
    if n_bcr > 0:
        pred = sum(1 for t in targets if t["bcr"] and t["s1"]) / n_bcr
    else:
        pred = None

    return {
        "arm": arm,
        "arm_dir": str(arm_dir),
        "n_targets_with_results": n,
        "n_s1": n_s1,
        "n_bcr": n_bcr,
        "n_election_gap": n_gap,
        "s1_rate": (n_s1 / n) if n else None,
        "bcr_rate": (n_bcr / n) if n else None,
        "predictive_power_P_S1_given_BCR": pred,
        "targets": targets,
        "raw_row_count": len(rows),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-out", required=True, type=Path)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args()
    base: Path = args.base_out
    if not base.is_dir():
        print(f"WARN: base-out not found yet: {base}", file=sys.stderr)
        return 1

    arms = {}
    for arm in ("FO", "DP", "CF"):
        d = base / arm
        if d.is_dir():
            arms[arm] = summarize_arm(arm, d)

    # paired comparison if both FO and DP present
    paired = []
    if "FO" in arms and "DP" in arms:
        fo_map = {t["pdb_id"]: t for t in arms["FO"]["targets"]}
        dp_map = {t["pdb_id"]: t for t in arms["DP"]["targets"]}
        for pid in sorted(set(fo_map) | set(dp_map)):
            fo = fo_map.get(pid)
            dp = dp_map.get(pid)
            paired.append(
                {
                    "pdb_id": pid,
                    "FO_s1": None if not fo else fo["s1"],
                    "DP_s1": None if not dp else dp["s1"],
                    "FO_bcr": None if not fo else fo["bcr"],
                    "DP_bcr": None if not dp else dp["bcr"],
                    "FO_rmsd_s1": None if not fo else fo["rmsd_s1"],
                    "DP_rmsd_s1": None if not dp else dp["rmsd_s1"],
                    "FO_rmsd_bcr": None if not fo else fo["rmsd_bcr"],
                    "DP_rmsd_bcr": None if not dp else dp["rmsd_bcr"],
                    "winner_s1": (
                        "tie"
                        if fo and dp and fo["s1"] == dp["s1"]
                        else (
                            "DP"
                            if dp and dp["s1"] and (not fo or not fo["s1"])
                            else (
                                "FO"
                                if fo and fo["s1"] and (not dp or not dp["s1"])
                                else "neither"
                            )
                        )
                    ),
                }
            )

    report = {
        "base_out": str(base),
        "arms": arms,
        "paired_FO_vs_DP": paired,
        "headline": {
            arm: {
                "s1_rate": arms[arm]["s1_rate"],
                "bcr_rate": arms[arm]["bcr_rate"],
                "predictive_power_P_S1_given_BCR": arms[arm][
                    "predictive_power_P_S1_given_BCR"
                ],
                "n": arms[arm]["n_targets_with_results"],
            }
            for arm in arms
        },
    }

    out_path = args.json_out or (base / "DPFO_SUMMARY.json")
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["headline"], indent=2))
    print(f"wrote {out_path}")

    # markdown
    md = base / "DPFO_SUMMARY.md"
    lines = [
        "# DP vs FO clustering pilot (small sim)",
        "",
        f"Base: `{base}`",
        "",
        "| Arm | N | S1 rate | BCR rate | P(S1\\|BCR) predictive |",
        "|-----|---|---------|----------|----------------------|",
    ]
    for arm, h in report["headline"].items():
        s1 = h["s1_rate"]
        bcr = h["bcr_rate"]
        pred = h["predictive_power_P_S1_given_BCR"]
        lines.append(
            f"| {arm} | {h['n']} | "
            f"{'—' if s1 is None else f'{s1:.0%}'} | "
            f"{'—' if bcr is None else f'{bcr:.0%}'} | "
            f"{'—' if pred is None else f'{pred:.0%}'} |"
        )
    lines.append("")
    if paired:
        lines += [
            "## Paired targets",
            "",
            "| PDB | FO S1 | DP S1 | FO BCR | DP BCR | FO rmsd | DP rmsd |",
            "|-----|-------|-------|--------|--------|---------|---------|",
        ]
        for p in paired:
            lines.append(
                f"| {p['pdb_id']} | {p['FO_s1']} | {p['DP_s1']} | "
                f"{p['FO_bcr']} | {p['DP_bcr']} | "
                f"{p['FO_rmsd_s1']} | {p['DP_rmsd_s1']} |"
            )
    md.write_text("\n".join(lines) + "\n")
    print(f"wrote {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
