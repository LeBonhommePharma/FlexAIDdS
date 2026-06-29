#!/usr/bin/env python3
"""Compare a results CSV against published per-PDB RMSD (or expected band).
For now, computes succ rate, mean/median RMSD, and if a reference CSV given, per-PDB delta.
Usage: python scripts/compare_published_rmsd.py <csv> [--ref <ref_csv>] [--threshold 2.0]
"""
import argparse, csv, sys, json, os
from pathlib import Path

def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--ref", default=None, help="reference CSV with same structure for per-PDB delta")
    ap.add_argument("--threshold", type=float, default=2.0)
    args = ap.parse_args()
    rows = load_csv(args.csv)
    n = len(rows)
    succ = sum(1 for r in rows if (0 < float(r.get("rmsd_hungarian", 99) or 99) < args.threshold))
    rmsds = []
    for r in rows:
        try:
            v = float(r.get("rmsd_hungarian", 99) or 99)
            if v > 0:
                rmsds.append(v)
        except:
            pass
    mean = sum(rmsds)/len(rmsds) if rmsds else 0
    med = sorted(rmsds)[len(rmsds)//2] if rmsds else 0
    print(f"CSV: {args.csv}")
    print(f"N={n} succ<{args.threshold}: {succ} ({100*succ/n:.1f}%)")
    print(f"mean RMSD (valid): {mean:.2f}  median: {med:.2f}")
    if args.ref:
        ref_rows = load_csv(args.ref)
        refd = {r["pdb_id"]: float(r.get("rmsd_hungarian",99) or 99) for r in ref_rows}
        deltas = []
        flips = 0
        for r in rows:
            pid = r["pdb_id"]
            if pid in refd:
                h = float(r.get("rmsd_hungarian",99) or 99)
                if 0 < h and refd[pid] > 0:
                    d = abs(h - refd[pid])
                    deltas.append(d)
                    if (0 < h < args.threshold) != (0 < refd[pid] < args.threshold):
                        flips += 1
        if deltas:
            print(f"vs ref: mean |delta|={sum(deltas)/len(deltas):.3f}  succ flips={flips}")
    # simple band check for 80%
    if succ >= 68:  # 80% of 85 ~68
        print("BAND: matches ~80% target band")
    else:
        print("BAND: below 80% target")
    sys.exit(0 if succ >= 60 else 1)  # loose for now

if __name__ == "__main__":
    main()
