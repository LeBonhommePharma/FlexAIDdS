from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path
from typing import Any


DEFAULT_WORK_DIR = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy"


def count_csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(newline="") as fh:
        return max(0, sum(1 for _ in csv.reader(fh)) - 1)


def latest_summary(work_dir: Path) -> tuple[Path | None, dict[str, Any] | None]:
    summaries = sorted((work_dir / "orchestrator_runs").glob("*/orchestrator_summary.json"))
    if not summaries:
        return None, None
    path = summaries[-1]
    try:
        return path, json.loads(path.read_text())
    except json.JSONDecodeError:
        return path, {"error": "unparseable_summary"}


def active_processes() -> list[dict[str, str]]:
    patterns = ("benchmark_datasets", "benchmarks.astex_entropy", "vina", "rbdock", "rbcavity", "boltz")
    try:
        result = subprocess.run(
            ["ps", "-ax", "-o", "pid,ppid,stat,etime,command"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return []
    rows: list[dict[str, str]] = []
    for line in result.stdout.splitlines()[1:]:
        if not any(pattern in line for pattern in patterns):
            continue
        parts = line.split(None, 4)
        if len(parts) < 5 or "astex_entropy_status.py" in parts[4]:
            continue
        rows.append({"pid": parts[0], "ppid": parts[1], "stat": parts[2], "etime": parts[3], "command": parts[4]})
    return rows


def collect_status(work_dir: Path) -> dict[str, Any]:
    summary_path, summary = latest_summary(work_dir)
    pose_counts = {
        path.name: count_csv_rows(path)
        for path in sorted((work_dir / "poses").glob("*_poses.csv"))
    }
    rescored_counts = {
        str(path.relative_to(work_dir)): count_csv_rows(path)
        for path in sorted((work_dir / "rescored").glob("*/*/rescored_poses.csv"))
    }
    flex_command = work_dir / "flexaidds/native/flexaidds_command.txt"
    return {
        "work_dir": str(work_dir),
        "work_dir_exists": work_dir.exists(),
        "latest_summary_path": str(summary_path) if summary_path else "",
        "latest_summary": summary or {},
        "pose_counts": pose_counts,
        "rescored_counts": rescored_counts,
        "active_processes": active_processes(),
        "flexaidds_native_command": flex_command.read_text().strip() if flex_command.exists() else "",
    }


def print_text(status: dict[str, Any]) -> None:
    print(f"work_dir: {status['work_dir']}")
    print(f"work_dir_exists: {status['work_dir_exists']}")
    print(f"latest_summary: {status['latest_summary_path'] or 'none'}")
    if status["latest_summary"]:
        latest = status["latest_summary"]
        print(f"run_id: {latest.get('run_id', 'unknown')}")
        print(f"modes: {','.join(latest.get('modes', []))}")
        print(f"tools: {','.join(latest.get('tools', []))}")
        print(f"dry_run: {latest.get('dry_run')}")
    print("active_processes:")
    if status["active_processes"]:
        for proc in status["active_processes"]:
            print(f"  {proc['pid']} {proc['etime']} {proc['command']}")
    else:
        print("  none")
    print("pose_counts:")
    for name, count in status["pose_counts"].items():
        print(f"  {name}: {count}")
    print("rescored_counts:")
    for name, count in status["rescored_counts"].items():
        print(f"  {name}: {count}")
    if status["flexaidds_native_command"]:
        print("flexaidds_native_command:")
        print(f"  {status['flexaidds_native_command']}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Report Astex entropy benchmark status.")
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    status = collect_status(args.work_dir)
    if args.json:
        print(json.dumps(status, indent=2, sort_keys=True))
    else:
        print_text(status)


if __name__ == "__main__":
    main()
