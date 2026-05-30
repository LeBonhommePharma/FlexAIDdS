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
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Any

# Robust import for both "python -m" and direct execution
try:
    from .ensure_docking_data import print_skill_banner
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from ensure_docking_data import print_skill_banner


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
    p.add_argument("--package", "--create-package", action="store_true",
                   help="After the run, automatically create a professional, shareable validation package (zip + reproducibility manifest)")

    return p


def _get_file_sha256(path: Path) -> str:
    if not path.exists():
        return "missing"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def gather_reproducibility_metadata(args: argparse.Namespace, binary_path: str | None) -> Dict[str, Any]:
    """Collect comprehensive reproducibility information."""
    meta: Dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "skill_version": "2026-05 (flexaid-docking)",
    }

    # Git information (best effort)
    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).parents[4], stderr=subprocess.DEVNULL
        ).decode().strip()
        meta["git_sha"] = git_sha
    except Exception:
        meta["git_sha"] = "unavailable"

    # Binary information
    if binary_path:
        bin_path = Path(binary_path)
        meta["binary_path"] = str(bin_path)
        meta["binary_sha256"] = _get_file_sha256(bin_path)
        try:
            meta["binary_size_bytes"] = bin_path.stat().st_size
        except Exception:
            pass

    # Key data files (from skill data/)
    skill_data = Path(__file__).parent.parent / "data"
    key_data_files = ["MC_st0r5.2_6.dat", "AMINO.def", "Lovell_LIB.dat"]
    meta["data_files"] = {}
    for name in key_data_files:
        p = skill_data / name
        if p.exists():
            meta["data_files"][name] = _get_file_sha256(p)

    meta["command_line"] = " ".join(sys.argv)
    meta["environment"] = {k: v for k, v in os.environ.items() if k.startswith(("FLEXAID", "PYTHON", "OMP", "MPI"))}

    return meta


def create_validation_package(results_dir: Path, metadata: Dict[str, Any], package_name: str | None = None) -> Path:
    """Create a clean, professional, shareable validation package."""
    if package_name is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        package_name = f"flexaidds_validation_package_{timestamp}"

    package_dir = results_dir.parent / package_name
    package_dir.mkdir(parents=True, exist_ok=True)

    # Copy key artifacts
    if results_dir.exists():
        for item in results_dir.rglob("*"):
            if item.is_file():
                rel = item.relative_to(results_dir)
                dest = package_dir / "results" / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                import shutil
                shutil.copy2(item, dest)

    # Write reproducibility manifest
    with open(package_dir / "REPRODUCIBILITY_MANIFEST.json", "w") as f:
        json.dump(metadata, f, indent=2, default=str)

    # Human-readable summary
    with open(package_dir / "README.txt", "w") as f:
        f.write("FlexAIDδS Validation / Benchmark Package\n")
        f.write("========================================\n\n")
        f.write(f"Generated: {metadata.get('timestamp_utc')}\n")
        f.write(f"Git SHA:   {metadata.get('git_sha', 'N/A')}\n")
        f.write(f"Binary:    {metadata.get('binary_path', 'N/A')}\n")
        f.write(f"Binary SHA256: {metadata.get('binary_sha256', 'N/A')}\n\n")
        f.write("This package contains all necessary artifacts for audit and reproduction.\n")
        f.write("See REPRODUCIBILITY_MANIFEST.json for full details.\n")

    # Create zip
    import shutil
    zip_path = package_dir.with_suffix(".zip")
    shutil.make_archive(str(zip_path.with_suffix("")), 'zip', package_dir)

    return zip_path


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

    # Gather reproducibility metadata before the run
    binary_for_meta = args.binary
    metadata = gather_reproducibility_metadata(args, binary_for_meta)

    try:
        result = subprocess.run(cmd, check=False)

        if args.package and result.returncode == 0:
            results_dir = Path(args.results_dir)
            print("\n[Skill] Creating professional validation package...")
            package_path = create_validation_package(results_dir, metadata)
            print(f"[Skill] Validation package created: {package_path}")
            print("       This archive is suitable for internal audit, collaboration, or regulatory purposes.")

        return result.returncode
    except KeyboardInterrupt:
        print("\n[Skill] Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
