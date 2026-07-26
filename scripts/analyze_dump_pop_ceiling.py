#!/usr/bin/env python3
"""Compare FLEXAIDDS_DUMP_POP .pop.tsv mins to emission BCR/S3 from result.csv.

Uses engine-written rmsd_raw/rmsd_sym columns (not a reimplementation of election).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def load_pop(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def summarize_pop(rows: list[dict]) -> dict:
    if not rows:
        raise ValueError("empty pop dump")
    for col in ("rmsd_raw", "rmsd_sym"):
        if col not in rows[0]:
            raise KeyError(f"missing column {col} in .pop.tsv (engine dump schema)")
    raw = [float(r["rmsd_raw"]) for r in rows]
    sym = [float(r["rmsd_sym"]) for r in rows]
    return {
        "n_dump": len(rows),
        "min_rmsd_raw": min(raw),
        "min_rmsd_sym": min(sym),
        "n_raw_le2": sum(1 for x in raw if x <= 2.0),
        "n_sym_le2": sum(1 for x in sym if x <= 2.0),
    }


def load_emission_bcr(result_csv: Path) -> dict:
    with result_csv.open(newline="") as f:
        row = next(csv.DictReader(f))
    bcr = float(row["best_cluster_rmsd"])
    s3 = float(row["conditional_scanned_pool_ceiling"])
    return {
        "best_cluster_rmsd": bcr,
        "conditional_scanned_pool_ceiling": s3,
        "elected_rmsd": float(row.get("rmsd_to_crystal") or "nan"),
    }


def compare(pop_summary: dict, emission: dict) -> dict:
    bcr = emission["best_cluster_rmsd"]
    delta = pop_summary["min_rmsd_sym"] - bcr
    if pop_summary["n_sym_le2"] > 0 and bcr > 2.0:
        verdict = "RETENTION"
    elif pop_summary["min_rmsd_sym"] + 0.25 < bcr:
        verdict = "MILD_RETENTION"
    elif abs(delta) <= 0.25:
        verdict = "SEARCH_WALL"
    else:
        verdict = "OTHER"
    return {
        **pop_summary,
        **emission,
        "delta_sym_vs_BCR": delta,
        "verdict": verdict,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir", type=Path, help="DatasetRunner OUT with CODE/result.csv + CODE.pop.tsv")
    ap.add_argument("--codes", default="", help="comma codes (default: discover *.pop.tsv)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)
    out = args.out_dir
    if args.codes:
        codes = [c.strip() for c in args.codes.split(",") if c.strip()]
    else:
        codes = sorted(p.parent.name for p in out.glob("*/*.pop.tsv"))
    results = []
    for code in codes:
        pop_path = out / code / f"{code}.pop.tsv"
        res_path = out / code / "result.csv"
        if not pop_path.is_file():
            print(f"MISSING_POP {code}", file=sys.stderr)
            return 2
        if not res_path.is_file():
            print(f"MISSING_RESULT {code}", file=sys.stderr)
            return 2
        rec = compare(summarize_pop(load_pop(pop_path)), load_emission_bcr(res_path))
        rec["target"] = code
        results.append(rec)
        print(
            f"{code}\tn_dump={rec['n_dump']}\tmin_sym={rec['min_rmsd_sym']:.4f}\t"
            f"BCR={rec['best_cluster_rmsd']:.4f}\tn_sym<=2={rec['n_sym_le2']}\t"
            f"verdict={rec['verdict']}"
        )
    if args.json:
        print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
