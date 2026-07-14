#!/usr/bin/env python3
"""Patch result.csv best_cluster_rmsd from on-disk poses via offline_bcr (CONECT extract).

Use when DatasetRunner writes RMSD=-1 (interrupted post-process) but pose PDBs exist.
Drives the real offline_bcr tool (LIB/PoseBust Loaders), not a reimplementation.

Usage:
  python3 scripts/patch_bcr_from_poses.py <campaign_dir> [--offline-bcr PATH] [--write]
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASTEX = ROOT / "benchmarks" / "astex_diverse" / "astex_diverse"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign_dir", type=Path)
    ap.add_argument(
        "--offline-bcr",
        type=Path,
        default=None,
        help="path to offline_bcr binary (default: search SCRATCH or build)",
    )
    ap.add_argument("--write", action="store_true", help="rewrite result.csv in place")
    args = ap.parse_args()
    tool = args.offline_bcr
    if tool is None:
        for cand in (
            Path(
                "/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/"
                "grok-goal-19923fdc9045/implementer/offline_bcr"
            ),
            ROOT / "build" / "offline_bcr",
        ):
            if cand.is_file():
                tool = cand
                break
    if tool is None or not tool.is_file():
        print("offline_bcr binary not found", file=sys.stderr)
        return 2

    camp = args.campaign_dir
    n_patch = 0
    for rcsv in sorted(camp.glob("*/result.csv")):
        pid = rcsv.parent.name
        crystal = ASTEX / pid / f"{pid}_ligand.sdf"
        if not crystal.is_file():
            print(f"SKIP {pid}: no crystal sdf")
            continue
        rows = list(csv.DictReader(rcsv.open()))
        if not rows:
            continue
        row = rows[0]
        try:
            bc = float(row.get("best_cluster_rmsd", -1))
        except ValueError:
            bc = -1.0
        if 0 <= bc <= 50:
            print(f"OK {pid}: already has BCR={bc}")
            continue
        try:
            out = subprocess.check_output(
                [str(tool), str(crystal), str(rcsv.parent), pid],
                text=True,
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError as e:
            print(f"FAIL {pid}: offline_bcr rc={e.returncode}")
            continue
        bcr_line = [ln for ln in out.splitlines() if ln.startswith("BCR")]
        if not bcr_line:
            print(f"FAIL {pid}: no BCR line")
            continue
        parts = bcr_line[0].split()
        bcr = float(parts[1])
        if bcr > 1e8:
            print(f"FAIL {pid}: invalid BCR {bcr}")
            continue
        print(f"PATCH {pid}: BCR {row.get('best_cluster_rmsd')} -> {bcr:.4f}")
        row["best_cluster_rmsd"] = f"{bcr:.4f}"
        n_patch += 1
        if args.write:
            with rcsv.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(row.keys()))
                w.writeheader()
                w.writerow(row)
    print(f"patched {n_patch} targets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
