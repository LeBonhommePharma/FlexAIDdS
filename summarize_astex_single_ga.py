#!/usr/bin/env python3
"""Summarize Astex single-GA DatasetRunner results from local CSV evidence."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def find_aggregate(result_dir: Path) -> Path | None:
    for name in ("astex_crossdock_85_results.csv", "astex_diverse_results.csv"):
        path = result_dir / name
        if path.exists():
            return path
    matches = sorted(result_dir.glob("*_results.csv"))
    return matches[0] if matches else None


def collect_rows(result_dir: Path) -> tuple[list[dict[str, str]], str]:
    aggregate = find_aggregate(result_dir)
    if aggregate:
        return read_rows(aggregate), str(aggregate)
    per_target = sorted(result_dir.glob("*/result.csv"))
    rows: list[dict[str, str]] = []
    for path in per_target:
        rows.extend(read_rows(path))
    if rows:
        return rows, f"{result_dir}/*/result.csv ({len(per_target)} files)"
    raise FileNotFoundError(f"no Astex result CSVs under {result_dir}")


def as_float(row: dict[str, str], key: str, default: float = -1.0) -> float:
    try:
        value = row.get(key, "")
        if value == "" or value.upper() == "NA" or value.lower() == "nan":
            return default
        return float(value)
    except Exception:
        return default


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir")
    args = parser.parse_args()

    rows, source = collect_rows(Path(args.result_dir))
    valid = [r for r in rows if 0.0 <= as_float(r, "rmsd_hungarian") < 900.0]
    successes = [r for r in valid if as_float(r, "rmsd_hungarian") < 2.0]
    rmsds = [as_float(r, "rmsd_hungarian") for r in valid]

    print(f"source={source}")
    print(f"n_rows={len(rows)}")
    print(f"n_valid_rmsd={len(valid)}")
    print(f"successes={len(successes)}/{len(rows)}")
    print(f"success_rate={len(successes) / len(rows):.4f}" if rows else "success_rate=NA")
    print(f"mean_rmsd={statistics.mean(rmsds):.4f}" if rmsds else "mean_rmsd=NA")
    print(f"median_rmsd={statistics.median(rmsds):.4f}" if rmsds else "median_rmsd=NA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
