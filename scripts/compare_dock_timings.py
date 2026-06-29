#!/usr/bin/env python3
"""Compare harvested dock timings between pre- and post-optimization campaigns."""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "results/perf_swarm"


def normalize_job_key(raw: str) -> str:
    """Strip campaign wrapper so pre/post resume paths align."""
    markers = ("/results_resume_missing_198/", "/results/")
    for marker in markers:
        idx = raw.find(marker)
        if idx >= 0:
            return raw[idx + len(marker) :]
    return raw


def load_records(path: Path) -> dict[str, float]:
    data = json.loads(path.read_text())
    records = data.get("dock_timings_harvested", {}).get("records", [])
    out: dict[str, float] = {}
    for row in records:
        raw = row.get("job_key") or row.get("target") or ""
        key = normalize_job_key(raw)
        if key:
            out[key] = float(row["avg_ms_per_gen"])
    return out


def pct_delta(current: float, baseline: float) -> float:
    if baseline <= 0:
        return 0.0
    return 100.0 * (current - baseline) / baseline


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT / "post_p0_comparison.json")
    parser.add_argument("--label", type=str, default="pre_vs_post_p0")
    args = parser.parse_args()

    base = load_records(args.baseline)
    curr = load_records(args.current)
    common = sorted(set(base) & set(curr))

    if not common:
        base_vals = list(base.values())
        curr_vals = list(curr.values())
        if not base_vals or not curr_vals:
            print("No timings in baseline or current.", file=sys.stderr)
            return 1

        def cohort_stats(vals: list[float]) -> dict:
            vals_sorted = sorted(vals)
            return {
                "count": len(vals_sorted),
                "median_ms_per_gen": round(statistics.median(vals_sorted), 2),
                "mean_ms_per_gen": round(statistics.mean(vals_sorted), 2),
                "p25_ms_per_gen": round(vals_sorted[len(vals_sorted) // 4], 2),
                "p75_ms_per_gen": round(vals_sorted[3 * len(vals_sorted) // 4], 2),
                "min_ms_per_gen": round(min(vals_sorted), 2),
                "max_ms_per_gen": round(max(vals_sorted), 2),
            }

        b_stats = cohort_stats(base_vals)
        c_stats = cohort_stats(curr_vals)
        median_speedup = b_stats["median_ms_per_gen"] / c_stats["median_ms_per_gen"]
        report = {
            "label": args.label,
            "mode": "cohort_distribution",
            "baseline_file": str(args.baseline),
            "current_file": str(args.current),
            "note": "No paired job keys; comparing aggregate distributions (different target sets).",
            "baseline_cohort": b_stats,
            "current_cohort": c_stats,
            "median_speedup": round(median_speedup, 3),
            "median_delta_pct": round(
                100.0 * (c_stats["median_ms_per_gen"] - b_stats["median_ms_per_gen"])
                / b_stats["median_ms_per_gen"],
                2,
            ),
        }
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Mode: cohort distribution (no paired keys)")
        print(
            f"Baseline median: {b_stats['median_ms_per_gen']} ms/gen  "
            f"Current median: {c_stats['median_ms_per_gen']} ms/gen  "
            f"Speedup: {median_speedup:.2f}x"
        )
        print(f"Wrote {args.out}")
        return 0

    deltas: list[float] = []
    speedups: list[float] = []
    rows: list[dict] = []
    for key in common:
        b = base[key]
        c = curr[key]
        d = pct_delta(c, b)
        sp = b / c if c > 0 else 1.0
        deltas.append(d)
        speedups.append(sp)
        rows.append(
            {
                "job_key": key,
                "baseline_ms_per_gen": b,
                "current_ms_per_gen": c,
                "delta_pct": round(d, 2),
                "speedup": round(sp, 3),
            }
        )

    rows.sort(key=lambda r: r["delta_pct"])
    median_delta = statistics.median(deltas)
    mean_delta = statistics.mean(deltas)
    median_speedup = statistics.median(speedups)
    mean_speedup = statistics.mean(speedups)

    report = {
        "label": args.label,
        "baseline_file": str(args.baseline),
        "current_file": str(args.current),
        "matched_jobs": len(common),
        "summary": {
            "median_delta_pct": round(median_delta, 2),
            "mean_delta_pct": round(mean_delta, 2),
            "median_speedup": round(median_speedup, 3),
            "mean_speedup": round(mean_speedup, 3),
            "improved_count": sum(1 for d in deltas if d < 0),
            "regressed_count": sum(1 for d in deltas if d > 0),
        },
        "best_5": rows[:5],
        "worst_5": rows[-5:][::-1],
        "all_matched": rows,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")

    s = report["summary"]
    print(f"Matched jobs: {len(common)}")
    print(
        f"Median: {s['median_delta_pct']:+.1f}% ({s['median_speedup']:.2f}x)  "
        f"Mean: {s['mean_delta_pct']:+.1f}% ({s['mean_speedup']:.2f}x)"
    )
    print(
        f"Improved: {s['improved_count']}  Regressed: {s['regressed_count']}"
    )
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())