#!/usr/bin/env python3
"""Offline re-election over existing Astex campaign pose trees (Chunk 1).

Re-ranks emitted cluster poses under CF-primary + pathological-CF vetoes without
re-running the GA. Reports virtual success rates vs campaign top-1.

Usage:
  python3 scripts/offline_reelect_cf_primary.py benchmarks/astex_repro/full_v131
  python3 scripts/offline_reelect_cf_primary.py benchmarks/astex_repro/full_v130 --cf-floor -250
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from pathlib import Path


def parse_pose_cf(pdb: Path) -> float | None:
    try:
        text = pdb.read_text(errors="ignore")
    except OSError:
        return None
    m = re.search(r"REMARK\s+CF=\s*([-+eE0-9.]+)", text)
    if not m:
        m = re.search(r"REMARK CF=([-+eE0-9.]+)", text)
    if not m:
        return None
    try:
        v = float(m.group(1))
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def parse_freq(pdb: Path) -> int:
    try:
        text = pdb.read_text(errors="ignore")
    except OSError:
        return 1
    m = re.search(r"Frequency:\s*(\d+)", text)
    return int(m.group(1)) if m else 1


def elect_cf_primary(poses: list[tuple[Path, float, int]], cf_floor: float) -> Path | None:
    """Lowest non-pathological CF; veto CF>=0 if any CF<0 exists."""
    if not poses:
        return None
    sane = [(p, cf, f) for p, cf, f in poses if cf >= cf_floor]
    if not sane:
        sane = poses
    any_neg = any(cf < 0 for _, cf, _ in sane)
    cand = [(p, cf, f) for p, cf, f in sane if (not any_neg) or cf < 0]
    if not cand:
        cand = sane
    cand.sort(key=lambda t: (t[1], -t[2]))  # CF asc, freq desc
    return cand[0][0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("results_dir", type=Path, help="Campaign output dir (full_v130/full_v131/...)")
    ap.add_argument("--cf-floor", type=float, default=-250.0)
    ap.add_argument("--csv-out", type=Path, default=None)
    args = ap.parse_args()
    root: Path = args.results_dir
    if not root.is_dir():
        print(f"ERROR: not a directory: {root}", file=sys.stderr)
        return 2

    rows_out = []
    n = n_top1 = n_cf = n_bcr = 0
    for result_csv in sorted(root.glob("*/result.csv")):
        pdb = result_csv.parent.name
        camp = list(csv.DictReader(result_csv.open()))
        if not camp:
            continue
        r = camp[0]
        n += 1
        try:
            top1 = float(r.get("rmsd_hungarian") or r.get("rmsd_to_crystal") or 99)
            bcr = float(r.get("best_cluster_rmsd") or 99)
            top1_ok = top1 <= 2.0 and top1 >= 0
            bcr_ok = 0 <= bcr <= 2.0
        except ValueError:
            top1_ok = bcr_ok = False
            top1 = bcr = float("nan")
        if top1_ok:
            n_top1 += 1
        if bcr_ok:
            n_bcr += 1

        poses: list[tuple[Path, float, int]] = []
        for pdb_path in result_csv.parent.rglob("*.pdb"):
            name = pdb_path.name
            if "_INI" in name or "member" in name or "cleft" in name:
                continue
            # cluster heads: PDB_N.pdb or rX/PDB_N.pdb
            if not re.search(r"_\d+\.pdb$", name):
                continue
            cf = parse_pose_cf(pdb_path)
            if cf is None:
                continue
            poses.append((pdb_path, cf, parse_freq(pdb_path)))

        elected = elect_cf_primary(poses, args.cf_floor)
        # Without re-RMSD to crystal we report CF-primary CF only + whether BCR was ok.
        # Full RMSD recompute needs crystal SDF; optional enhancement.
        cf_e = parse_pose_cf(elected) if elected else None
        row = {
            "pdb_id": pdb,
            "campaign_success": int(top1_ok),
            "campaign_rmsd": f"{top1:.4f}" if top1 == top1 else "",
            "best_cluster_rmsd": f"{bcr:.4f}" if bcr == bcr else "",
            "bcr_le_2": int(bcr_ok),
            "n_poses_scanned": len(poses),
            "cf_primary_path": str(elected) if elected else "",
            "cf_primary_cf": f"{cf_e:.4f}" if cf_e is not None else "",
            "campaign_best_score": r.get("best_score", ""),
        }
        rows_out.append(row)
        # Virtual: if BCR ok, CF-primary *can* succeed if it picks that basin;
        # we conservatively count bcr as upper bound already.
        if bcr_ok:
            n_cf += 1  # placeholder: true CF-primary success needs RMSD of elected

    print(f"Campaign: {root}")
    print(f"Targets with result.csv: {n}")
    print(f"  campaign top-1 success (rmsd<=2): {n_top1}/{n}")
    print(f"  search ceiling BCR<=2:            {n_bcr}/{n}")
    print(f"  (BCR ceiling = upper bound if CF-primary recovers native cluster)")
    print(f"  poses scanned with REMARK CF across tree: "
          f"{sum(int(r['n_poses_scanned']) for r in rows_out)}")

    out = args.csv_out or (root / "offline_cf_primary_reelect.csv")
    if rows_out:
        with out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
            w.writeheader()
            w.writerows(rows_out)
        print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
