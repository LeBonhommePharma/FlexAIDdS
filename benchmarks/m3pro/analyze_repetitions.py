#!/usr/bin/env python3
"""analyze_repetitions.py — Bootstrap analysis for FlexAIDdS benchmark repetitions.

Loads N run-level CSV results per dataset, computes 95% CIs via bootstrap
(10,000 resamples), evaluates coefficient of variation, and generates a final
report recommending whether more runs are needed.

Usage:
    python benchmarks/m3pro/analyze_repetitions.py --results-dir RESULTS/tier2
    python benchmarks/m3pro/analyze_repetitions.py --results-dir RESULTS/tier2 --n-bootstrap 10000
    python benchmarks/m3ro/analyze_repetitions.py --results-dir RESULTS/tier2 --dataset astex

Apache-2.0 (c) 2026 NRGlab, Universite de Montreal
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


N_BOOTSTRAP_DEFAULT = 10_000
CI_LEVEL = 0.95
CV_SUFFICIENT = 0.05
CV_MODERATE = 0.10
MIN_RUNS = 30


def load_run_csv(csv_path: Path) -> dict[str, float]:
    """Load a per-target results CSV into run-level aggregate metrics.

    ``benchmark_datasets`` writes one row per target. The repetition analysis
    needs one scalar per run, so each numeric column is reduced to a run-level
    mean and the ``success`` column is reduced to both count and rate.
    """
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    metrics: dict[str, float] = {"total_systems": float(len(rows))}
    if not rows:
        return metrics

    numeric: dict[str, list[float]] = defaultdict(list)
    successes: list[float] = []

    for row in rows:
        for key, val in row.items():
            if key == "pdb_id":
                continue
            try:
                fval = float(val)
            except (ValueError, TypeError):
                continue
            if key == "success":
                successes.append(fval)
            else:
                numeric[key].append(fval)

    for key, values in numeric.items():
        if values:
            metrics[key] = float(np.mean(values))

    if successes:
        metrics["successful"] = float(np.sum(successes))
        metrics["success_rate"] = float(np.mean(successes))
    return metrics


def load_aggregate_metrics(csv_path: Path) -> dict[str, float]:
    """Load summary CSV with aggregate metrics (mean_rmsd, success_rate, etc.)."""
    metrics: dict[str, float] = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row.get("metric", row.get("Metric", ""))
            try:
                val = float(row.get("value", row.get("Value", 0)))
                metrics[key] = val
            except (ValueError, TypeError):
                pass
    return metrics


def collect_dataset_runs(dataset_dir: Path) -> list[dict[str, float]]:
    """Collect all run-level results from a dataset directory."""
    runs: list[dict[str, float]] = []

    if not dataset_dir.is_dir():
        return runs

    for run_dir in sorted(dataset_dir.iterdir()):
        if not run_dir.is_dir() or not run_dir.name.startswith("run"):
            continue

        csv_files = list(run_dir.glob("*_results.csv"))
        if csv_files:
            metrics = load_run_csv(csv_files[0])
            if metrics:
                runs.append(metrics)
                continue

        summary_files = list(run_dir.glob("*_summary.csv"))
        if summary_files:
            metrics = load_aggregate_metrics(summary_files[0])
            if metrics:
                runs.append(metrics)

    return runs


def bootstrap_ci(
    values: np.ndarray,
    n_bootstrap: int = N_BOOTSTRAP_DEFAULT,
    ci_level: float = CI_LEVEL,
) -> dict:
    """Bootstrap confidence interval for a set of scalar values."""
    n = len(values)
    if n < 2:
        return {
            "mean": float(np.mean(values)) if n == 1 else None,
            "std": 0.0,
            "ci_lo": None,
            "ci_hi": None,
            "cv": None,
            "n": n,
            "verdict": "insufficient_data",
        }

    rng = np.random.default_rng(42)
    boot_indices = rng.integers(0, n, size=(n_bootstrap, n))
    boot_means = np.array([np.mean(values[idx]) for idx in boot_indices])

    alpha = 1.0 - ci_level
    ci_lo = float(np.percentile(boot_means, 100 * alpha / 2))
    ci_hi = float(np.percentile(boot_means, 100 * (1 - alpha / 2)))
    mean_val = float(np.mean(values))
    std_val = float(np.std(values, ddof=1))
    cv = std_val / mean_val if abs(mean_val) > 1e-12 else float("inf")

    if cv < CV_SUFFICIENT:
        verdict = "sufficient"
    elif cv < CV_MODERATE:
        verdict = "consider_more"
    else:
        verdict = "need_more"

    return {
        "mean": mean_val,
        "std": std_val,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "ci_width": ci_hi - ci_lo,
        "cv": cv,
        "n": n,
        "n_bootstrap": n_bootstrap,
        "verdict": verdict,
    }


def compute_run_level_metrics(runs: list[dict[str, float]]) -> dict[str, np.ndarray]:
    """Transpose run-level dicts into metric→array-of-values."""
    all_keys: set[str] = set()
    for run in runs:
        all_keys.update(run.keys())

    metric_arrays: dict[str, np.ndarray] = {}
    for key in sorted(all_keys):
        vals = [run.get(key, np.nan) for run in runs]
        arr = np.array(vals, dtype=float)
        arr = arr[~np.isnan(arr)]
        if len(arr) > 0:
            metric_arrays[key] = arr

    return metric_arrays


def analyze_dataset(
    dataset_name: str,
    runs: list[dict[str, float]],
    n_bootstrap: int,
) -> dict:
    """Full bootstrap analysis for a single dataset."""
    metric_arrays = compute_run_level_metrics(runs)

    results: dict = {
        "dataset": dataset_name,
        "n_runs": len(runs),
        "metrics": {},
    }

    for metric_name, values in metric_arrays.items():
        ci_result = bootstrap_ci(values, n_bootstrap)
        results["metrics"][metric_name] = ci_result

    return results


def generate_markdown_report(
    all_results: list[dict],
    output_path: Path,
    n_bootstrap: int,
) -> None:
    """Generate a publication-ready Markdown report."""
    lines: list[str] = []
    lines.append("# FlexAIDdS Benchmark Repetition Analysis")
    lines.append("")
    lines.append("## Bootstrap Confidence Intervals (10,000 resamples)")
    lines.append("")

    for ds_result in all_results:
        ds_name = ds_result["dataset"]
        n_runs = ds_result["n_runs"]
        lines.append(f"### {ds_name} ({n_runs} runs)")
        lines.append("")
        lines.append("| Metric | Mean | 95% CI | CV | Verdict |")
        lines.append("|--------|------|--------|-----|---------|")

        for metric_name, ci in sorted(ds_result["metrics"].items()):
            if ci["mean"] is None:
                continue
            mean_str = f"{ci['mean']:.4f}"
            if ci["ci_lo"] is not None:
                ci_str = f"[{ci['ci_lo']:.4f}, {ci['ci_hi']:.4f}]"
            else:
                ci_str = "N/A"
            cv_str = f"{ci['cv']:.4f}" if ci["cv"] is not None else "N/A"
            verdict = ci["verdict"]
            lines.append(f"| {metric_name} | {mean_str} | {ci_str} | {cv_str} | {verdict} |")

        lines.append("")

    lines.append("## Run Sufficiency Summary")
    lines.append("")
    lines.append("| Dataset | Runs | Status | Recommendation |")
    lines.append("|---------|------|--------|----------------|")

    for ds_result in all_results:
        verdicts = [ci["verdict"] for ci in ds_result["metrics"].values()]
        if all(v == "sufficient" for v in verdicts):
            status = "PASS"
            rec = "No more runs needed"
        elif any(v == "need_more" for v in verdicts):
            status = "EXTEND"
            n_bad = sum(1 for v in verdicts if v == "need_more")
            rec = f"Add 20 more runs ({n_bad} metrics with CV>10%)"
        else:
            status = "CONSIDER"
            n_moderate = sum(1 for v in verdicts if v == "consider_more")
            rec = f"Consider 10 more runs ({n_moderate} metrics with CV 5-10%)"

        lines.append(f"| {ds_result['dataset']} | {ds_result['n_runs']} | {status} | {rec} |")

    lines.append("")
    lines.append("---")
    lines.append(f"*Generated with {n_bootstrap} bootstrap resamples per metric*")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    print(f"  Markdown report: {output_path}")


def generate_json_report(all_results: list[dict], output_path: Path) -> None:
    """Generate machine-readable JSON report."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"  JSON report: {output_path}")


def generate_campaign_status(
    all_results: list[dict],
    output_path: Path,
    min_success_rate: float = 0.70,
    max_cv_warn: float = 0.10,
    max_cv_fail: float = 0.15,
) -> dict:
    """
    Generate a hardened, machine-readable campaign status with quality gates.

    Thresholds are now configurable (refined in step A of the plan).
    """
    overall_status = "PASS"
    recommendations: list[str] = []
    worst_verdict = "sufficient"
    success_rates: list[float] = []
    max_cv_seen = 0.0

    for ds_result in all_results:
        for metric_name, ci in ds_result.get("metrics", {}).items():
            verdict = ci.get("verdict", "unknown")
            cv = ci.get("cv") or 0.0
            mean = ci.get("mean")

            if cv > max_cv_seen:
                max_cv_seen = cv

            # Use configurable CV thresholds + existing verdict
            if cv >= max_cv_fail or verdict == "need_more":
                worst_verdict = "need_more"
                overall_status = "FAIL"
                recommendations.append(
                    f"{ds_result['dataset']}/{metric_name}: FAIL (CV={cv:.3f} >= {max_cv_fail})"
                )
            elif cv >= max_cv_warn or verdict == "consider_more":
                if overall_status != "FAIL":
                    overall_status = "WARN"
                if worst_verdict not in ("need_more",):
                    worst_verdict = "consider_more"
                recommendations.append(
                    f"{ds_result['dataset']}/{metric_name}: WARN (CV={cv:.3f} >= {max_cv_warn})"
                )

            # Collect success_rate
            if metric_name == "success_rate" and mean is not None:
                success_rates.append(mean)

    if not recommendations:
        recommendations.append("All metrics look stable. Campaign quality is good.")

    # Aggregate success rates
    avg_success = None
    min_success = None
    if success_rates:
        avg_success = float(sum(success_rates) / len(success_rates))
        min_success = float(min(success_rates))

        # Apply success_rate threshold
        if min_success < min_success_rate and overall_status == "PASS":
            overall_status = "WARN"
            recommendations.append(
                f"Minimum success_rate {min_success:.3f} below threshold {min_success_rate}"
            )

    # Stability
    total_metrics = sum(len(ds.get("metrics", {})) for ds in all_results)
    sufficient_metrics = sum(
        1 for ds in all_results
        for m in ds.get("metrics", {}).values()
        if m.get("verdict") == "sufficient"
    )
    stability = round((sufficient_metrics / total_metrics) * 100, 1) if total_metrics > 0 else 0.0

    # Actionable
    if overall_status == "FAIL":
        action = "STOP or significantly extend the campaign before trusting results."
    elif overall_status == "WARN":
        action = "Proceed with caution. Consider adding more repetitions for key metrics."
    else:
        action = "Safe to proceed. Results are statistically stable."

    status = {
        "generated_at": __import__("datetime").datetime.utcnow().isoformat() + "Z",
        "overall_status": overall_status,
        "datasets_analyzed": len(all_results),
        "worst_verdict": worst_verdict,
        "stability_percent": stability,
        "thresholds": {
            "min_success_rate": min_success_rate,
            "max_cv_warn": max_cv_warn,
            "max_cv_fail": max_cv_fail,
        },
        "success_rate": {
            "average": avg_success,
            "minimum": min_success,
        },
        "max_cv_observed": round(max_cv_seen, 4),
        "recommendations": recommendations,
        "actionable": action,
        "requires_attention": overall_status != "PASS",
        "details": all_results,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(status, indent=2) + "\n")
    print(f"  Campaign status: {output_path}  (status={overall_status})")

    return status


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Bootstrap analysis of FlexAIDdS benchmark repetitions"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing per-dataset run subdirectories",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=N_BOOTSTRAP_DEFAULT,
        help=f"Number of bootstrap resamples (default: {N_BOOTSTRAP_DEFAULT}). "
             "Use lower values (1000-2000) for fast periodic quality-gate checks with --check-quality-every.",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        default=None,
        help="Analyze only this dataset (default: all)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for reports (default: results-dir/../analysis/)",
    )
    # P1 Quality Gate Thresholds (A: refined & configurable)
    parser.add_argument(
        "--min-success-rate",
        type=float,
        default=0.70,
        help="Minimum acceptable success_rate for PASS (default: 0.70)",
    )
    parser.add_argument(
        "--max-cv-warn",
        type=float,
        default=0.10,
        help="CV above this triggers WARN (default: 0.10)",
    )
    parser.add_argument(
        "--max-cv-fail",
        type=float,
        default=0.15,
        help="CV above this triggers FAIL (default: 0.15)",
    )

    args = parser.parse_args()

    results_dir: Path = args.results_dir
    if not results_dir.is_dir():
        print(f"ERROR: Results directory not found: {results_dir}", file=sys.stderr)
        return 1

    output_dir = args.output_dir or results_dir.parent / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset_dirs: list[Path] = []
    if args.dataset:
        target = results_dir / args.dataset
        if target.is_dir():
            dataset_dirs = [target]
        else:
            print(f"ERROR: Dataset directory not found: {target}", file=sys.stderr)
            return 1
    else:
        for d in sorted(results_dir.iterdir()):
            if d.is_dir():
                dataset_dirs.append(d)

    if not dataset_dirs:
        print(f"ERROR: No dataset directories found in {results_dir}", file=sys.stderr)
        return 1

    print(f"Bootstrap Analysis: {len(dataset_dirs)} datasets, {args.n_bootstrap} resamples")
    print("")

    all_results: list[dict] = []

    for ds_dir in dataset_dirs:
        ds_name = ds_dir.name
        runs = collect_dataset_runs(ds_dir)

        if not runs:
            print(f"  [{ds_name}] No run results found — skipping")
            continue

        print(f"  [{ds_name}] {len(runs)} runs loaded")

        ds_analysis = analyze_dataset(ds_name, runs, args.n_bootstrap)
        all_results.append(ds_analysis)

        for metric_name, ci in sorted(ds_analysis["metrics"].items()):
            if ci["mean"] is None:
                continue
            ci_str = f"[{ci['ci_lo']:.4f}, {ci['ci_hi']:.4f}]" if ci["ci_lo"] is not None else "N/A"
            print(
                f"    {metric_name}: {ci['mean']:.4f} {ci_str} CV={ci['cv']:.4f} ({ci['verdict']})"
            )

    if not all_results:
        print("No results to analyze.")
        return 1

    print("")
    generate_markdown_report(all_results, output_dir / "bootstrap_report.md", args.n_bootstrap)
    generate_json_report(all_results, output_dir / "bootstrap_report.json")

    # P1-1: Quality gate / automated decision making (refined thresholds)
    campaign_status_path = output_dir / "campaign_status.json"
    generate_campaign_status(
        all_results,
        campaign_status_path,
        min_success_rate=args.min_success_rate,
        max_cv_warn=args.max_cv_warn,
        max_cv_fail=args.max_cv_fail,
    )

    print("")
    print("Analysis complete.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
