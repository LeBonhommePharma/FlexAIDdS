"""CLI entry point for the FlexAIDdS DatasetRunner.

THERMO/ITC (itc187, bindingdb_itc, scorpio) are automatically treated with:
- 298 K forced
- scoring+entropy_rescue metrics only
- full provenance (git_sha + binary + host + temp + command) always in outputs

Examples
--------
Run a tier-1 (PR sanity) benchmark on CASF-2016::

    python -m flexaidds.dataset_runner --dataset casf2016 --tier 1

ITC thermo dry-run (auto 298 K, correct metrics, provenance)::

    python -m flexaidds.dataset_runner --dataset itc187 --tier 1 --dry-run
    python -m flexaidds.dataset_runner --dataset itc187,bindingdb_itc --tier 1 --dry-run --output /tmp/out.json

Run all datasets at tier-2 with 4 MPI nodes::

    mpirun -n 4 python -m flexaidds.dataset_runner --all --tier 2 --distributed

Dry run to test the pipeline without actual docking::

    python -m flexaidds.dataset_runner --all --tier 1 --dry-run
"""

from __future__ import annotations

import argparse
import logging
import os
import shlex
import sys
from datetime import datetime
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
        help="Local worker threads for parallel target evaluation (default: 1).",
    )
    p.add_argument(
        "--omp-threads",
        type=int,
        default=None,
        metavar="N",
        help=(
            "OMP_NUM_THREADS to set per FlexAID subprocess. "
            "Default: 2 when --workers > 1, else 4. "
            "Override via FLEXAIDDS_OMP_THREADS env var."
        ),
    )
    p.add_argument(
        "--conc", "--concentration",
        type=float,
        default=1.0,
        metavar="M",
        dest="default_conc_M",
        help="Default concentration in M for grand canonical (P3; per-ligand from dataset yaml overrides).",
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
        default=None,
        metavar="DIR",
        help=(
            "Output directory for reports and per-dataset JSON. "
            "Default (omitted): prefers $FLEXAIDDS_ICLOUD/results/working/<ts> when set "
            "(or $FLEXAIDDS_RESULTS), else 'results/benchmarks'. "
            "Writes to working/ subdir for active runs (iCloud sync-safe pattern). "
            "Override always honored."
        ),
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
        help=(
            "Skip actual docking; use synthetic poses for pipeline smoke tests. "
            "docking_power_* is omitted (not real docking success rates); remaining "
            "metrics are also synthetic and must not be reported as production results."
        ),
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume a previous (or partially completed) run by skipping targets that already have per-entry result files. Enables fine-grained checkpointing and crash recovery.",
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

    # Resolve --results-dir preferring iCloud (FLEXAIDDS_ICLOUD > FLEXAIDDS_RESULTS).
    # Active runs default into working/ + timestamped dir (safer-than-safe for iCloud Drive).
    # Comments: iCloud sync can lag / produce placeholders / conflicted-copies during heavy concurrent writes.
    # Use working/ for live runs; promote final outputs via scripts/safe_archive_to_icoud.py --source ... --dest .../archived/
    # (the archiver does full SHA verification + atomic + never-evict-until-verified).
    # --results-dir override is always respected for full compatibility.
    if not getattr(args, "results_dir", None):
        icloud = os.environ.get("FLEXAIDDS_ICLOUD")
        results_env = os.environ.get("FLEXAIDDS_RESULTS")
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        if icloud:
            base = Path(icloud).expanduser() / "results" / "working"
            args.results_dir = str(base / f"flexaidds_bench_{ts}")
        elif results_env:
            base = Path(results_env).expanduser() / "working"
            args.results_dir = str(base / f"flexaidds_bench_{ts}")
        else:
            args.results_dir = "results/benchmarks"
        Path(args.results_dir).mkdir(parents=True, exist_ok=True)

    # Logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    logger = logging.getLogger(__name__)

    # Override data dir from CLI
    if args.data_dir:
        os.environ["FLEXAIDDS_BENCHMARK_DATA"] = args.data_dir

    # Import early for thermo constants and runner
    from .runner import DatasetRunner, BenchmarkReport, _git_sha, _runner_info, THERMO_DATASETS

    # -------------------------------------------------------------------------
    # Bulletproof ITC/thermo mode (handoff + validation prompt)
    # - Auto-force 298 K for itc187 / bindingdb_itc / scorpio (and sequences)
    # - Will be further restricted in runner for metrics + provenance
    # -------------------------------------------------------------------------
    if args.dataset:
        slugs = {s.strip().lower() for s in args.dataset.split(",")}
        if slugs & THERMO_DATASETS and args.temperature == 300.0:
            logger.info(
                "ITC/thermo dataset(s) %s detected — forcing --temperature 298.0 (reproducible as fuck)",
                args.dataset,
            )
            args.temperature = 298.0

    full_command = " ".join(shlex.quote(x) for x in sys.argv)
    runner_kwargs = dict(
        results_dir=args.results_dir,
        binary=args.binary,
        temperature=args.temperature,
        n_workers=args.workers,
        omp_threads=args.omp_threads,
        use_mpi=args.distributed,
        cache_dir=args.data_dir,
        bootstrap_ci=args.bootstrap,
        n_bootstrap=args.n_bootstrap,
        dry_run=args.dry_run,
        resume=args.resume,
        command_line=full_command,
        default_conc_M=getattr(args, 'default_conc_M', 1.0),
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
        import datetime as _datetime, socket
        report = BenchmarkReport(
            datasets=[dr],
            generated_at=_datetime.datetime.utcnow().isoformat() + "Z",
            git_sha=_git_sha(),
            host=socket.gethostname(),
            runner_info=_runner_info(),
            temperature=dr.temperature,
            binary=dr.binary,
            full_command=full_command,
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
