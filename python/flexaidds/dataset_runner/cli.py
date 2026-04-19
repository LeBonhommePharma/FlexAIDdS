"""CLI entry point for the FlexAIDdS DatasetRunner.

Examples
--------
Run a tier-1 (PR sanity) benchmark on CASF-2016::

    python -m flexaidds.dataset_runner --dataset casf2016 --tier 1

Run all datasets at tier-2 with 4 MPI nodes::

    mpirun -n 4 python -m flexaidds.dataset_runner --all --tier 2 --distributed

Run a single metric on ITC-187::

    python -m flexaidds.dataset_runner --dataset itc187 --metric entropy_rescue_rate

Dry run to test the pipeline without actual docking::

    python -m flexaidds.dataset_runner --all --tier 1 --dry-run
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="flexaidds-benchmark",
        description="FlexAIDdS DatasetRunner — distributed molecular docking benchmarks",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # --- Dataset selection ---
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--dataset", "-d",
        metavar="SLUG",
        help="Run a single dataset by slug (e.g. casf2016, itc187).",
    )
    group.add_argument(
        "--all", "-a",
        action="store_true",
        help="Run all discovered datasets.",
    )

    # --- Tier ---
    p.add_argument(
        "--tier", "-t",
        type=int,
        choices=[1, 2],
        default=2,
        help=(
            "Benchmark tier: 1 = fast PR-sanity subset, "
            "2 = full comprehensive run (default: 2)."
        ),
    )

    # --- Metric filter ---
    p.add_argument(
        "--metric", "-m",
        metavar="NAME",
        default=None,
        help=(
            "Compute only this metric (e.g. entropy_rescue_rate, docking_power_top1). "
            "Default: all metrics defined in the dataset config."
        ),
    )

    # --- Distributed / parallel ---
    p.add_argument(
        "--distributed",
        action="store_true",
        help="Enable MPI-distributed execution (requires mpi4py; launch with mpirun).",
    )
    p.add_argument(
        "--nodes",
        type=int,
        default=1,
        metavar="N",
        help="Number of MPI nodes (informational; actual count is set by mpirun).",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Local worker processes for parallel target evaluation (default: 1).",
    )

    # --- I/O ---
    p.add_argument(
        "--datasets-dir",
        default=None,
        metavar="DIR",
        help=(
            "Directory containing dataset YAML configs. "
            "Default: shipped configs inside the flexaidds package."
        ),
    )
    p.add_argument(
        "--results-dir",
        default="results/benchmarks",
        metavar="DIR",
        help="Output directory for reports and per-dataset JSON (default: results/benchmarks).",
    )
    p.add_argument(
        "--data-dir",
        default=None,
        metavar="DIR",
        help=(
            "Root directory for cached dataset files. "
            "Overrides the FLEXAIDDS_BENCHMARK_DATA env variable."
        ),
    )

    # --- Engine ---
    p.add_argument(
        "--binary",
        default=None,
        metavar="PATH",
        help=(
            "Path to the FlexAID executable. "
            "Defaults to $FLEXAIDDS_BINARY env var, then PATH."
        ),
    )
    p.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        metavar="K",
        help="Simulation temperature in Kelvin (default: 300).",
    )

    # --- Bootstrap CIs ---
    p.add_argument(
        "--bootstrap",
        action="store_true",
        help="Compute 95%% bootstrap confidence intervals (slower).",
    )
    p.add_argument(
        "--n-bootstrap",
        type=int,
        default=5000,
        metavar="N",
        help="Number of bootstrap resamples when --bootstrap is set (default: 5000).",
    )

    # --- Misc ---
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Skip actual docking; use synthetic scores to test the pipeline.",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable DEBUG-level logging.",
    )
    p.add_argument(
        "--report-prefix",
        default=None,
        metavar="PATH",
        help=(
            "Save the final report to PATH.json and PATH.md. "
            "Default: results/benchmarks/report_<timestamp>."
        ),
    )

    return p


def main(argv: list[str] | None = None) -> int:
    """Entry point.  Returns 0 on success, non-zero on failure/regression."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    # Override data dir from CLI
    if args.data_dir:
        os.environ["FLEXAIDDS_BENCHMARK_DATA"] = args.data_dir

    # Import here so logging is already configured
    from .runner import DatasetRunner, BenchmarkReport, _git_sha, _runner_info

    runner_kwargs = dict(
        results_dir=args.results_dir,
        binary=args.binary,
        temperature=args.temperature,
        n_workers=args.workers,
        use_mpi=args.distributed,
        cache_dir=args.data_dir,
        bootstrap_ci=args.bootstrap,
        n_bootstrap=args.n_bootstrap,
        dry_run=args.dry_run,
    )
    if args.datasets_dir is not None:
        runner_kwargs["datasets_dir"] = args.datasets_dir

    runner = DatasetRunner(**runner_kwargs)

    # ----- Run -----
    if args.all:
        report = runner.run_all(
            tier=args.tier,
            distributed=args.distributed,
            n_nodes=args.nodes,
            metric_subset=[args.metric] if args.metric else None,
        )
    else:
        # Single dataset
        try:
            dr = runner.run_single(
                dataset_slug=args.dataset,
                tier=args.tier,
                metric=args.metric,
            )
        except FileNotFoundError as exc:
            logging.error("%s", exc)
            return 2

        # Build report for single-dataset case
        import datetime, socket
        report = BenchmarkReport(
            datasets=[dr],
            generated_at=datetime.datetime.utcnow().isoformat() + "Z",
            git_sha=_git_sha(),
            host=socket.gethostname(),
            runner_info=_runner_info(),
        )

    # Only root rank prints / saves
    if not runner._mpi_root:
        return 0

    # ----- Report output -----
    import datetime as _dt
    timestamp = _dt.datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    prefix = args.report_prefix or str(
        Path(args.results_dir) / f"report_{timestamp}"
    )
    json_path, md_path = report.save(prefix)
    print(f"\nReport saved:\n  JSON: {json_path}\n  Markdown: {md_path}\n")

    # Print markdown summary to stdout
    print(report.to_markdown())

    # Return non-zero if any regressions detected
    any_regression = any(
        any(dr.regression_flags.values())
        for dr in report.datasets
    )
    if any_regression:
        logging.warning("REGRESSION DETECTED — one or more metrics dropped below baseline.")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
