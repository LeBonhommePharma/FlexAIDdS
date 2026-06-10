#!/usr/bin/env python3
"""
A/B comparison report for FlexAIDdS targeted diagnostics (P10 matrix swap, P7 fine grid).

Reads the summary result.csv (or per-complex result.csv files) from two result
directories and prints a per-target side-by-side of rmsd_hungarian, best-of-N
cluster RMSD, and cf_native, plus how many cross the 2.0 A sub-2 line in each arm.

Usage:
  python3 ab_compare_report.py LABEL_A DIR_A LABEL_B DIR_B
"""
import sys, os, csv, glob


def load(d):
    rows = {}
    # Prefer per-complex result.csv (always present); fall back to a summary csv.
    files = glob.glob(os.path.join(d, "*", "result.csv"))
    if not files:
        files = glob.glob(os.path.join(d, "*.csv"))
    for f in files:
        try:
            with open(f) as fh:
                for r in csv.DictReader(fh):
                    if r.get("pdb_id"):
                        rows[r["pdb_id"]] = r
        except Exception:
            pass
    return rows


def fnum(r, k):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return None


def main():
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(1)
    la, da, lb, db = sys.argv[1:5]
    A, B = load(os.path.expanduser(da)), load(os.path.expanduser(db))
    codes = sorted(set(A) | set(B))

    print(f"{'CODE':6} | {la:>22} | {lb:>22} | delta")
    print(f"{'':6} | {'rmsdH  bestN  cf_nat':>22} | {'rmsdH  bestN  cf_nat':>22} | rmsdH")
    print("-" * 78)
    subA = subB = 0
    for c in codes:
        ra, rb = A.get(c, {}), B.get(c, {})
        ha, hb = fnum(ra, "rmsd_hungarian"), fnum(rb, "rmsd_hungarian")
        ba, bb = fnum(ra, "best_cluster_rmsd"), fnum(rb, "best_cluster_rmsd")
        ca, cb = fnum(ra, "cf_native"), fnum(rb, "cf_native")
        if ha is not None and ha < 2.0:
            subA += 1
        if hb is not None and hb < 2.0:
            subB += 1

        def cell(h, b, cf):
            hs = f"{h:5.2f}" if h is not None else "  -- "
            bs = f"{b:5.2f}" if b is not None else "  -- "
            cs = f"{cf:8.1f}" if cf is not None else "   --   "
            return f"{hs} {bs} {cs}"
        d = (f"{hb-ha:+.2f}" if (ha is not None and hb is not None) else " -- ")
        flag = ""
        if ha is not None and hb is not None:
            if ha >= 2.0 and hb < 2.0:
                flag = "  <-- CROSSED sub-2"
            elif ha < 2.0 and hb >= 2.0:
                flag = "  <-- LOST sub-2"
        print(f"{c:6} | {cell(ha,ba,ca):>22} | {cell(hb,bb,cb):>22} | {d:>6}{flag}")
    print("-" * 78)
    print(f"sub-2A (rmsd_hungarian top-1):  {la} = {subA}/{len(codes)}   "
          f"{lb} = {subB}/{len(codes)}")


if __name__ == "__main__":
    main()
