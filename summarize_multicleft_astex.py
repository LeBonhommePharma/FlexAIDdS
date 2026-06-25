#!/usr/bin/env python3
"""Summarize restored multi-cleft Astex runs.

Input is the DatasetRunner CSV where each row is one cleft variant
(`PDB__clfN`). The selected result is the lowest docking score per original
PDB. The RMSD-best row is reported separately as an oracle diagnostic only.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


DEFAULT_SINGLE_GA_SUMMARY = Path(
    "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v112_20260624_2130_oracle_full85/"
    "astex_diverse_summary.csv"
)
DEFAULT_SINGLE_GA_RESULTS = Path(
    "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v112_20260624_2130_oracle_full85/"
    "astex_diverse_results.csv"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        val = row.get(key, "")
        if val == "" or val.lower() == "nan":
            return default
        return float(val)
    except Exception:
        return default


def original_id(pdb_id: str) -> str:
    return pdb_id.split("__clf", 1)[0]


def summarize_single_ga_rows(rows: list[dict[str, str]]) -> dict[str, str]:
    valid = [r for r in rows if 0.0 <= f(r, "rmsd_hungarian", -1.0) < 900.0]
    successes = [r for r in valid if f(r, "rmsd_hungarian", -1.0) < 2.0]
    rmsds = [f(r, "rmsd_hungarian", -1.0) for r in valid]
    return {
        "total": str(len(rows)),
        "successful": str(len(successes)),
        "success_rate": f"{len(successes) / len(rows):.4f}" if rows else "0.0000",
        "mean_rmsd": f"{statistics.mean(rmsds):.4f}" if rmsds else "nan",
        "median_rmsd": f"{statistics.median(rmsds):.4f}" if rmsds else "nan",
    }


def find_results_csv(result_dir: Path) -> Path | None:
    preferred = [
        result_dir / "astex_crossdock_85_results.csv",
        result_dir / "astex_diverse_results.csv",
    ]
    for path in preferred:
        if path.exists():
            return path
    matches = sorted(result_dir.glob("*_results.csv"))
    return matches[0] if matches else None


def read_multicleft_rows(result_dir: Path, explicit_csv: str) -> tuple[list[dict[str, str]], str]:
    if explicit_csv:
        path = Path(explicit_csv)
        return read_rows(path), str(path)

    per_cleft = sorted(result_dir.glob("*/result.csv"))
    if per_cleft:
        rows: list[dict[str, str]] = []
        for path in per_cleft:
            rows.extend(read_rows(path))
        return rows, f"{result_dir}/*/result.csv ({len(per_cleft)} files)"

    aggregate = find_results_csv(result_dir)
    if aggregate:
        return read_rows(aggregate), str(aggregate)

    raise FileNotFoundError(f"no aggregate or per-cleft result.csv files under {result_dir}")


def summarize_multicleft(rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, float]]:
    groups: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(original_id(row["pdb_id"]), []).append(row)

    selected_rows: list[dict[str, str]] = []
    for code, items in sorted(groups.items()):
        selected = min(items, key=lambda r: f(r, "best_score", 1.0e9))
        valid_rmsd_items = [r for r in items if f(r, "rmsd_hungarian", -1.0) >= 0.0]
        oracle_best = min(valid_rmsd_items, key=lambda r: f(r, "rmsd_hungarian", 1.0e9)) if valid_rmsd_items else selected
        selected_rmsd = f(selected, "rmsd_hungarian", -1.0)
        oracle_rmsd = f(oracle_best, "rmsd_hungarian", -1.0) if valid_rmsd_items else -1.0
        selected_rows.append({
            "pdb_id": code,
            "selected_cleft": selected["pdb_id"],
            "selected_score": selected.get("best_score", ""),
            "selected_rmsd_hungarian": f"{selected_rmsd:.4f}",
            "selected_success": "1" if 0.0 <= selected_rmsd < 2.0 else "0",
            "oracle_best_cleft": oracle_best["pdb_id"],
            "oracle_best_rmsd_hungarian": f"{oracle_rmsd:.4f}",
            "n_clefts": str(len(items)),
        })

    rmsds = [float(r["selected_rmsd_hungarian"]) for r in selected_rows if float(r["selected_rmsd_hungarian"]) >= 0]
    successes = sum(int(r["selected_success"]) for r in selected_rows)
    metrics = {
        "total": float(len(selected_rows)),
        "successes": float(successes),
        "success_rate": successes / len(selected_rows) if selected_rows else 0.0,
        "mean_rmsd": statistics.mean(rmsds) if rmsds else -1.0,
        "median_rmsd": statistics.median(rmsds) if rmsds else -1.0,
        "oracle_successes": float(sum(
            1 for r in selected_rows
            if 0.0 <= float(r["oracle_best_rmsd_hungarian"]) < 2.0
        )),
    }
    return selected_rows, metrics


def read_current_single_ga_baseline(summary_path: Path, results_path: Path) -> dict[str, str]:
    if summary_path.exists():
        return read_rows(summary_path)[0]
    if results_path.exists():
        return summarize_single_ga_rows(read_rows(results_path))
    return {
        "total": "85",
        "successful": "38",
        "success_rate": "0.4471",
        "mean_rmsd": "3.5955",
        "median_rmsd": "2.1130",
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_comparison(
    path: Path,
    metrics: dict[str, float],
    source_csv: str,
    current_single_ga: dict[str, str],
    current_single_ga_results: Path,
) -> None:
    lines = [
        "# Astex Multi-Cleft Restoration Comparison",
        "",
        "| Campaign | Evidence | Cleft handling | Native/site handling | N | Successes <2 A | Success % | Mean RMSD | Median RMSD | Notes |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---|",
        "| Old 2017-2019 multi-cleft | evidence pending | independent GA per major cleft | evidence pending | TBD | TBD | TBD | TBD | TBD | Do not fill until old logs/table are recovered. |",
        (
            "| Current single-GA v112 | "
            f"{current_single_ga_results} | single GA | self-dock/oracle reference protocol | "
            f"{current_single_ga.get('total_systems', current_single_ga.get('total', '85'))} | "
            f"{current_single_ga.get('successful', '38')} | "
            f"{100.0 * float(current_single_ga.get('success_rate', '0')):.1f} | "
            f"{float(current_single_ga.get('mean_rmsd', 'nan')):.2f} | "
            f"{float(current_single_ga.get('median_rmsd', 'nan')):.2f} | Local current comparator. |"
        ),
        (
            "| Revived multi-cleft | "
            f"{source_csv} | one DatasetRunner entry per Get_Cleft sphere file; selected by best_score | autonomous; no oracle-site env | "
            f"{int(metrics['total'])} | "
            f"{int(metrics['successes'])} | "
            f"{100.0 * metrics['success_rate']:.1f} | "
            f"{metrics['mean_rmsd']:.2f} | "
            f"{metrics['median_rmsd']:.2f} | "
            f"Oracle-best diagnostic successes across clefts: {int(metrics['oracle_successes'])}. |"
        ),
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", required=True, help="DatasetRunner result directory")
    parser.add_argument("--results-csv", default="")
    parser.add_argument("--single-ga-summary", default=str(DEFAULT_SINGLE_GA_SUMMARY))
    parser.add_argument("--single-ga-results", default=str(DEFAULT_SINGLE_GA_RESULTS))
    parser.add_argument("--output-prefix", default="")
    args = parser.parse_args()

    result_dir = Path(args.result_dir)
    rows, source_csv = read_multicleft_rows(result_dir, args.results_csv)
    selected_rows, metrics = summarize_multicleft(rows)
    current_single_ga = read_current_single_ga_baseline(
        Path(args.single_ga_summary),
        Path(args.single_ga_results),
    )

    prefix = Path(args.output_prefix) if args.output_prefix else result_dir.parent / "multicleft_target"
    target_csv = prefix.with_suffix(".csv")
    comparison_md = prefix.parent / "multicleft_comparison.md"
    write_csv(target_csv, selected_rows)
    write_comparison(
        comparison_md,
        metrics,
        source_csv,
        current_single_ga,
        Path(args.single_ga_results),
    )

    print(f"source_csv={source_csv}")
    print(f"target_csv={target_csv}")
    print(f"comparison_md={comparison_md}")
    print(f"selected_successes={int(metrics['successes'])}/{int(metrics['total'])}")
    print(f"selected_success_rate={metrics['success_rate']:.4f}")
    print(f"oracle_best_successes={int(metrics['oracle_successes'])}/{int(metrics['total'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
