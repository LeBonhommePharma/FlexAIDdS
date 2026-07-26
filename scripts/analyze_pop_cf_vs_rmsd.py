#!/usr/bin/env python3
"""Label SCORING_PULL vs SEARCH_FAIL from engine DUMP_POP .pop.tsv columns.

Uses shipped columns: rmsd_sym, cf_total (and optional cf_com/cf_wal).
Does not reimplement CF or election.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def spearman(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 3:
        return float("nan")

    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    if dx < 1e-12 or dy < 1e-12:
        return float("nan")
    return num / (dx * dy)


def load_pop(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty pop dump: {path}")
    for col in ("rmsd_sym", "cf_total"):
        if col not in rows[0]:
            raise KeyError(f"missing engine column {col} in {path}")
    return rows


def diagnose(rows: list[dict], *, gross_rmsd_floor: float = 8.0) -> dict:
    rmsd = [float(r["rmsd_sym"]) for r in rows]
    cf = [float(r["cf_total"]) for r in rows]
    n = len(rows)
    i_rm = min(range(n), key=lambda i: rmsd[i])
    i_cf = min(range(n), key=lambda i: cf[i])  # most negative
    min_r = rmsd[i_rm]
    rmsd_at_best_cf = rmsd[i_cf]
    delta = rmsd_at_best_cf - min_r
    k = max(1, n // 10)
    order_cf = sorted(range(n), key=lambda i: cf[i])
    top_rm = sum(rmsd[i] for i in order_cf[:k]) / k
    bot_rm = sum(rmsd[i] for i in order_cf[-k:]) / k
    sp = spearman(cf, rmsd)

    if min_r > gross_rmsd_floor:
        if delta <= 1.0 and top_rm <= bot_rm + 1.0:
            label = "SEARCH_FAIL"
        elif delta >= 3.0 or top_rm > bot_rm + 2.0:
            label = "SCORING_PULL"
        else:
            label = "MIXED"
    else:
        if delta >= 2.0:
            label = "SCORING_PULL"
        elif delta <= 0.5:
            label = "SEARCH_NEAR_MISS"
        else:
            label = "MIXED"

    return {
        "n": n,
        "min_rmsd_sym": min_r,
        "cf_at_min_rmsd": cf[i_rm],
        "min_cf": cf[i_cf],
        "rmsd_at_min_cf": rmsd_at_best_cf,
        "delta_rmsd_bestCF_minus_minRMSD": delta,
        "spearman_cf_vs_rmsd": sp,
        "top10pct_cf_mean_rmsd": top_rm,
        "bottom10pct_cf_mean_rmsd": bot_rm,
        "label": label,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--codes", default="", help="comma codes; default: all */*.pop.tsv")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)
    out = args.out_dir
    if args.codes:
        codes = [c.strip() for c in args.codes.split(",") if c.strip()]
    else:
        codes = sorted({p.parent.name for p in out.glob("*/*.pop.tsv")})
    results = []
    for code in codes:
        pop = out / code / f"{code}.pop.tsv"
        if not pop.is_file():
            print(f"MISSING {pop}", file=sys.stderr)
            return 2
        rec = diagnose(load_pop(pop))
        rec["target"] = code
        results.append(rec)
        print(
            f"{code}\tn={rec['n']}\tmin_rmsd={rec['min_rmsd_sym']:.4f}\t"
            f"rmsd@bestCF={rec['rmsd_at_min_cf']:.4f}\t"
            f"spearman={rec['spearman_cf_vs_rmsd']:.4f}\tlabel={rec['label']}"
        )
    if args.json:
        print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
