#!/usr/bin/env python3
"""
dataset_runner.py

Part of the flexaid-docking skill.

Convenient, high-quality wrapper around the FlexAIDδS DatasetRunner
with built-in safety, data checks, and helpful defaults.

This script follows the same professional patterns as the rest of the skill:
- Clear banners
- Integration with ensure_docking_data
- Support for --dry-run
- Helpful guidance and guardrails

Usage examples:
    # Basic tier-1 run on Astex Diverse (fast)
    python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --dataset astex_diverse --tier 1

    # Full distributed campaign
    mpirun -n 8 python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --all --tier 2 --distributed

    # Dry run to validate everything without docking
    python3 .grok/skills/flexaid-docking/scripts/dataset_runner.py \
        --all --tier 1 --dry-run
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from .ensure_docking_data import print_skill_banner  # reuse the nice banner style


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="flexaid-docking dataset-runner",
        description="Run FlexAIDδS DatasetRunner campaigns with skill-integrated safety and convenience.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fast sanity check on Astex Diverse
  python3 .../dataset_runner.py --dataset astex_diverse --tier 1

  # Full campaign with rich reports
  python3 .../dataset_runner.py --all --tier 2 --results-dir results/my_benchmark

  # Distributed run (launch via mpirun)
  mpirun -n 8 python3 .../dataset_runner.py --all --tier 2 --distributed --binary /path/to/FlexAIDδS

Always run ensure_docking_data.py first (or let this script remind you).
""",
    )

    p.add_argument("--dataset", "-d", help="Run a single dataset by slug (e.g. astex_diverse, casf2016, itc187)")
    p.add_argument("--all", "-a", action="store_true", help="Run all available datasets")
    p.add_argument("--tier", "-t", type=int, choices=[1, 2], default=2,
                   help="Benchmark tier (1=fast sanity, 2=full). Default: 2")
    p.add_argument("--metric", "-m", help="Run only a specific metric")
    p.add_argument("--distributed", action="store_true", help="Enable MPI distributed mode (launch with mpirun)")
    p.add_argument("--nodes", type=int, default=1, help="Informational: number of nodes")
    p.add_argument("--workers", type=int, default=1, help="Local parallel workers")
    p.add_argument("--binary", help="Path to FlexAIDδS binary (overrides FLEXAIDDS_BINARY)")
    p.add_argument("--results-dir", default="results/benchmarks",
                   help="Where to write reports (default: results/benchmarks)")
    p.add_argument("--dry-run", action="store_true",
                   help="Skip actual docking; useful for pipeline validation")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--ensure-data", action="store_true", default=True,
                   help="Automatically run ensure_docking_data.py first (default: True)")
    p.add_argument("--no-ensure-data", dest="ensure_data", action="store_false")

    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    print_skill_banner(verbose=args.verbose)

    if args.ensure_data:
        print("\n[Skill] Ensuring all critical runtime data is present before benchmarking...")
        ensure_script = Path(__file__).parent / "ensure_docking_data.py"
        subprocess.run([sys.executable, str(ensure_script)], check=False)

    # Build the command for the real module
    cmd = [sys.executable, "-m", "flexaidds.dataset_runner"]

    if args.dataset:
        cmd += ["--dataset", args.dataset]
    elif args.all:
        cmd += ["--all"]
    else:
        print("Error: You must specify either --dataset <slug> or --all")
        return 2

    cmd += ["--tier", str(args.tier)]

    if args.metric:
        cmd += ["--metric", args.metric]
    if args.distributed:
        cmd += ["--distributed"]
    if args.nodes and args.nodes > 1:
        cmd += ["--nodes", str(args.nodes)]
    if args.workers and args.workers > 1:
        cmd += ["--workers", str(args.workers)]
    if args.binary:
        cmd += ["--binary", args.binary]
    if args.results_dir:
        cmd += ["--results-dir", args.results_dir]
    if args.dry_run:
        cmd += ["--dry-run"]
    if args.verbose:
        cmd += ["--verbose"]

    print(f"\n[Skill] Launching DatasetRunner:")
    print("  " + " ".join(cmd))
    print()

    try:
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except KeyboardInterrupt:
        print("\n[Skill] Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
