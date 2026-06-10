#!/usr/bin/env python3
# =============================================================================
# compare_all_astex.py — unified FlexAIDdS / rDock / Vina Astex-85 comparison
# =============================================================================
# Merges the three per-engine result CSVs into one side-by-side table and prints
# a summary. All engines are scored identically: TOP-1 emitted pose, element-
# matched Hungarian RMSD vs the crystal ligand, sub-2 A success.
#
# Inputs (any subset; missing engines just show blanks):
#   FlexAIDdS : astex_diverse_results.csv  (rmsd_to_crystal column)
#   rDock     : rdock_astex_results.csv    (from parse_rdock_results.py)
#   Vina      : vina_astex_results.csv     (from parse_vina_results.py)
#
# Output: astex_three_way.csv
# =============================================================================
import argparse
import csv
import os

HOME_RES = os.path.expanduser("~/flexaidds_benchmark_results")
THRESH = 2.0


def _load_rmsd(path, col="rmsd_to_crystal"):
    """pdb_id -> rmsd float (or None) from a results CSV."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            v = r.get(col, "")
            try:
                out[r["pdb_id"]] = float(v) if v not in ("", None) else None
            except ValueError:
                out[r["pdb_id"]] = None
    return out


def _succ(rmsd, thr):
    return int(rmsd is not None and rmsd < thr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flexaidds", default=os.path.join(HOME_RES, "astex_diverse_results.csv"))
    ap.add_argument("--rdock", default=os.path.join(HOME_RES, "rdock_astex", "rdock_astex_results.csv"))
    ap.add_argument("--vina", default=os.path.join(HOME_RES, "vina_astex", "vina_astex_results.csv"))
    ap.add_argument("--out", default=os.path.join(HOME_RES, "astex_three_way.csv"))
    ap.add_argument("--threshold", type=float, default=THRESH)
    args = ap.parse_args()

    fa = _load_rmsd(args.flexaidds)
    rd = _load_rmsd(args.rdock)
    vn = _load_rmsd(args.vina)

    codes = sorted(set(fa) | set(rd) | set(vn))
    if not codes:
        print("No engine CSVs found. Run the per-engine parsers first.")
        return

    thr = args.threshold
    fa_n = rd_n = vn_n = 0
    rows = []
    for c in codes:
        fr, rr, vr = fa.get(c), rd.get(c), vn.get(c)
        fs, rs, vs = _succ(fr, thr), _succ(rr, thr), _succ(vr, thr)
        fa_n += fs; rd_n += rs; vn_n += vs
        rows.append({
            "pdb_id": c,
            "flexaidds_rmsd": "" if fr is None else f"{fr:.4f}",
            "rdock_rmsd": "" if rr is None else f"{rr:.4f}",
            "vina_rmsd": "" if vr is None else f"{vr:.4f}",
            "flexaidds_success": fs,
            "rdock_success": rs,
            "vina_success": vs,
        })

    fields = ["pdb_id", "flexaidds_rmsd", "rdock_rmsd", "vina_rmsd",
              "flexaidds_success", "rdock_success", "vina_success"]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n = len(codes)
    def pct(x):
        return f"{x}/{n} ({100.0*x/n:.1f}%)" if n else "0/0"
    print(f"Astex-85 head-to-head (sub-{thr} A, top-1, Hungarian RMSD)")
    print(f"  FlexAIDdS : {pct(fa_n)}")
    print(f"  rDock     : {pct(rd_n)}")
    print(f"  Vina      : {pct(vn_n)}")
    print(f"  -> {args.out}")


if __name__ == "__main__":
    main()
