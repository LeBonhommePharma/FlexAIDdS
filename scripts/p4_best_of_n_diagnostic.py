#!/usr/bin/env python3
"""
P4 — Best-of-N-clusters diagnostic (NO recompile).

For every Astex Diverse target in a FlexAIDdS results directory, compute the
Hungarian (optimal element-matched assignment) heavy-atom RMSD between each
emitted pose PDB (<CODE>_<N>.pdb) and the crystal ligand SDF, then take the
minimum across all emitted poses. This is the "oracle ceiling": the best RMSD a
perfect cluster-selector could have reached given the poses the GA actually
emitted.

Reported:
  - per-target best-of-N RMSD and which pose index achieved it
  - how many targets have min_RMSD < 2.0 A (the oracle ceiling)
  - distribution at <1.0, <1.5, <2.0, <3.0, <5.0 A
  - for comparison, the rank-0 pose RMSD (what the run actually reported)

Hungarian RMSD: build a cost matrix of squared inter-atomic distances between
the two heavy-atom sets, forbid cross-element matches (large cost), solve the
optimal assignment with scipy.optimize.linear_sum_assignment, RMSD = sqrt(mean).
This matches the semantics of the `rmsd_hungarian` column the runner emits.

Usage:
  python3 p4_best_of_n_diagnostic.py [RESULTS_DIR] [DATASET_DIR] [--csv OUT.csv]

Defaults:
  RESULTS_DIR = ~/flexaidds_results/v22_20260609_allfix
  DATASET_DIR = <repo>/benchmarks/astex_diverse/astex_diverse
"""
import os
import sys
import glob
import math
import argparse

import numpy as np
from scipy.optimize import linear_sum_assignment

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_RESULTS = os.path.expanduser("~/flexaidds_results/v22_20260609_allfix")
DEFAULT_DATASET = os.path.join(REPO, "benchmarks", "astex_diverse", "astex_diverse")

# Hydrogen-like element tokens FlexAID may emit (H tagged "Du" per pose conventions).
H_TOKENS = {"H", "D", "DU"}

def discover_pose_pdbs(tdir: str, code: str, *, recursive: bool = False, max_depth: int = 3):
    """Return persisted emission pose PDBs for *code* under *tdir*.

    Default (recursive=False): top-level CODE_[0-9]*.pdb only (legacy).
    recursive=True: rglob CODE_*.pdb up to max_depth under tdir, exclude INI,
    include elected_pose.pdb, dedupe by content SHA256.
    """
    import hashlib
    tdir_p = Path if False else __import__("pathlib").Path(tdir)
    if not recursive:
        poses = sorted(
            glob.glob(os.path.join(tdir, f"{code}_[0-9].pdb"))
            + glob.glob(os.path.join(tdir, f"{code}_[0-9][0-9].pdb"))
            + glob.glob(os.path.join(tdir, f"{code}_[0-9][0-9][0-9].pdb"))
        )
        return poses
    found = []
    for p in tdir_p.rglob(f"{code}_*.pdb"):
        try:
            rel = p.relative_to(tdir_p)
        except ValueError:
            continue
        if len(rel.parts) > max_depth:
            continue
        if "INI" in p.name.upper():
            continue
        found.append(p)
    ep = tdir_p / "elected_pose.pdb"
    if ep.is_file():
        found.append(ep)
    by_hash = {}
    for p in sorted(found, key=lambda x: str(x)):
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        by_hash.setdefault(h, str(p))
    return list(by_hash.values())



def parse_pdb_ligand(path):
    """Return (coords Nx3 float, elements list[str]) for the DOCKED ligand heavy atoms.

    The emitted pose PDB also carries crystallographic waters (HOH) and ions as
    HETATM records. The docked ligand is the HETATM block with serial >= 90000
    (the FlexAID 9000X convention, matching its CONECT records).
    """
    coords, elems = [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            try:
                serial = int(line[6:11])
            except ValueError:
                continue
            if serial < 90000:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            elem = line[76:78].strip().upper()
            if not elem:
                # Fall back to atom name (cols 13-16), strip digits.
                elem = "".join(c for c in line[12:16] if c.isalpha()).upper()[:2]
            if elem in H_TOKENS:
                continue
            coords.append((x, y, z))
            elems.append(elem)
    return np.asarray(coords, dtype=float), elems


def parse_sdf_ligand(path):
    """Return (coords Nx3 float, elements list[str]) for heavy atoms of first molecule."""
    coords, elems = [], []
    with open(path) as fh:
        lines = fh.readlines()
    # counts line is the 4th line (index 3) in V2000.
    if len(lines) < 4:
        return np.zeros((0, 3)), []
    try:
        natoms = int(lines[3][0:3])
    except ValueError:
        return np.zeros((0, 3)), []
    for i in range(4, 4 + natoms):
        line = lines[i]
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
        except (ValueError, IndexError):
            continue
        elem = line[31:34].strip().upper()
        if elem in H_TOKENS:
            continue
        coords.append((x, y, z))
        elems.append(elem)
    return np.asarray(coords, dtype=float), elems


def hungarian_rmsd(coordsA, elemsA, coordsB, elemsB):
    """Element-matched optimal-assignment heavy-atom RMSD. Returns None if shapes mismatch."""
    if coordsA.shape[0] == 0 or coordsB.shape[0] == 0:
        return None
    if coordsA.shape[0] != coordsB.shape[0]:
        return None
    n = coordsA.shape[0]
    # Squared-distance cost matrix.
    diff = coordsA[:, None, :] - coordsB[None, :, :]
    cost = np.einsum("ijk,ijk->ij", diff, diff)  # n x n squared distances
    # Forbid cross-element matches with a large penalty.
    BIG = 1.0e9
    ea = np.asarray(elemsA)
    eb = np.asarray(elemsB)
    mask = ea[:, None] != eb[None, :]
    cost = cost + mask * BIG
    row, col = linear_sum_assignment(cost)
    matched = cost[row, col]
    if np.any(matched >= BIG):
        # Element multiset mismatch — fall back to ordered RMSD (same extraction source).
        d = coordsA - coordsB
        return float(math.sqrt(np.mean(np.einsum("ij,ij->i", d, d))))
    return float(math.sqrt(matched.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir", nargs="?", default=DEFAULT_RESULTS)
    ap.add_argument("dataset_dir", nargs="?", default=DEFAULT_DATASET)
    ap.add_argument("--csv", default=None)
    ap.add_argument("--recursive", action="store_true",
                    help="include r*/ restart trees (depth<=3), excl INI, sha256-dedupe")
    ap.add_argument("--only-codes", default="",
                    help="comma-separated PDB codes to analyze (default: all)")
    args = ap.parse_args()

    results_dir = os.path.expanduser(args.results_dir)
    dataset_dir = os.path.expanduser(args.dataset_dir)

    codes = sorted(
        d for d in os.listdir(results_dir)
        if os.path.isdir(os.path.join(results_dir, d)) and d[0].isdigit()
    )
    if args.only_codes.strip():
        allow = {c.strip().upper() for c in args.only_codes.split(",") if c.strip()}
        codes = [c for c in codes if c.upper() in allow]

    rows = []
    for code in codes:
        tdir = os.path.join(results_dir, code)
        crystal = os.path.join(dataset_dir, code, f"{code}_ligand.sdf")
        if not os.path.exists(crystal):
            rows.append((code, None, None, None, "no_crystal_sdf"))
            continue
        cC, cE = parse_sdf_ligand(crystal)
        poses = discover_pose_pdbs(tdir, code, recursive=args.recursive)
        if not poses:
            rows.append((code, None, None, len(cC), "no_poses"))
            continue
        best_rmsd, best_idx, rank0_rmsd = None, None, None
        for p in poses:
            pC, pE = parse_pdb_ligand(p)
            r = hungarian_rmsd(pC, pE, cC, cE)
            if r is None:
                continue
            base = os.path.basename(p)
            idx = base[len(code) + 1:-4]
            if idx == "0":
                rank0_rmsd = r
            if best_rmsd is None or r < best_rmsd:
                best_rmsd, best_idx = r, idx
        rows.append((code, best_rmsd, best_idx, rank0_rmsd, f"{len(poses)}poses"))

    # Distribution.
    thresholds = [1.0, 1.5, 2.0, 3.0, 5.0]
    bins = {t: 0 for t in thresholds}
    valid = [r for r in rows if r[1] is not None]
    for r in valid:
        for t in thresholds:
            if r[1] < t:
                bins[t] += 1
    rank0_sub2 = sum(1 for r in rows if r[3] is not None and r[3] < 2.0)

    print(f"Results dir : {results_dir}")
    print(f"Dataset dir : {dataset_dir}")
    print(f"Targets     : {len(rows)}  (with valid RMSD: {len(valid)})")
    print()
    print(f"{'CODE':6} {'best_of_N':>10} {'pose':>5} {'rank0':>8}  note")
    for code, best, idx, rank0, note in rows:
        bs = f"{best:.4f}" if best is not None else "   --   "
        r0 = f"{rank0:.4f}" if rank0 is not None else "  --  "
        ix = idx if idx is not None else "-"
        print(f"{code:6} {bs:>10} {ix:>5} {r0:>8}  {note}")

    print()
    print("=== ORACLE CEILING (best-of-N over emitted poses) ===")
    n = len(valid)
    for t in thresholds:
        pct = 100.0 * bins[t] / n if n else 0.0
        print(f"  min_RMSD < {t:>3} A : {bins[t]:3d}/{n}  ({pct:5.1f}%)")
    print()
    print(f"  best-of-N sub-2A (oracle ceiling) : {bins[2.0]}/{n}")
    print(f"  rank-0 sub-2A (what run reported)  : {rank0_sub2}/{len(rows)}")
    gap = bins[2.0] - rank0_sub2
    print(f"  selection gap (oracle - rank0)     : {gap}  "
          f"(targets recoverable by better cluster selection alone)")

    if args.csv:
        import csv
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["code", "best_of_n_rmsd", "best_pose_idx", "rank0_rmsd", "note"])
            for row in rows:
                w.writerow(row)
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
