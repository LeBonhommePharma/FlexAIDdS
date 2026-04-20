#!/usr/bin/env python3
"""Bonhomme Fleet Dashboard — FlexAIDdS Benchmark Campaign Monitor.

Production-grade terminal dashboard that monitors parallel benchmark campaigns
on an M3 Pro Mac. Auto-discovers fleet status JSON files, deep-scans result
directories for per-target progress, and displays a rich ANSI terminal UI.

Runners supported:
  - Claude Code (PID 96971): build/benchmark_datasets → fleet_status.json
  - opencode (GLM-5.1): build_opencode/benchmark_datasets → fleet_status_opencode.json

Usage:
  # Live dashboard, refresh every 30s (default)
  python3 benchmarks/m3pro/dashboard/fleet_monitor.py

  # Single snapshot
  python3 benchmarks/m3pro/dashboard/fleet_monitor.py --once

  # Machine-readable JSON for piping
  python3 benchmarks/m3pro/dashboard/fleet_monitor.py --json

  # Custom poll interval
  python3 benchmarks/m3pro/dashboard/fleet_monitor.py --interval 10

  # Override paths
  python3 benchmarks/m3pro/dashboard/fleet_monitor.py \\
      --results ~/.flexaidds_fast/results \\
      --icloud ~/Library/Mobile\\ Documents/com~apple~CloudDocs/FlexAIDdS

Environment variables:
  FLEXAIDDS_FAST_BASE  Override fast results path (default: ~/.flexaidds_fast)
  FLEXAIDDS_ICLOUD     Override iCloud FlexAIDdS path
  FLEXAIDDS_LOGS       Override logs path
  NO_COLOR             Disable ANSI colors (standard env var)

Apache-2.0 (c) 2026 NRGlab, Universite de Montreal
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from collections import OrderedDict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from glob import glob
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

__version__ = "1.0.0"

ICLOUD_DEFAULT = os.path.expanduser(
    "~/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS"
)
FAST_DEFAULT = os.path.expanduser("~/.flexaidds_fast")
RESULTS_SUBDIR = "results/tier2"
BOX_HORIZONTAL = "═"
BOX_VERTICAL = "║"
BOX_TL = "╔"
BOX_TR = "╗"
BOX_BL = "╚"
BOX_BR = "╝"
BOX_LJ = "╠"
BOX_RJ = "╣"
PROGRESS_FILLED = "█"
PROGRESS_EMPTY = "░"


class Colors:
    """ANSI color codes with auto-detection for piped output."""

    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[1;33m"
    BLUE = "\033[0;34m"
    MAGENTA = "\033[0;35m"
    CYAN = "\033[0;36m"
    WHITE = "\033[0;37m"
    BG_RED = "\033[41m"
    BG_GREEN = "\033[42m"

    def __init__(self, enabled: bool = True):
        self._enabled = enabled

    def _apply(self, code: str, text: str) -> str:
        if not self._enabled:
            return text
        return f"{code}{text}{self.RESET}"

    def bold(self, text: str) -> str:
        return self._apply(self.BOLD, text)

    def dim(self, text: str) -> str:
        return self._apply(self.DIM, text)

    def red(self, text: str) -> str:
        return self._apply(self.RED, text)

    def green(self, text: str) -> str:
        return self._apply(self.GREEN, text)

    def yellow(self, text: str) -> str:
        return self._apply(self.YELLOW, text)

    def blue(self, text: str) -> str:
        return self._apply(self.BLUE, text)

    def cyan(self, text: str) -> str:
        return self._apply(self.CYAN, text)

    def magenta(self, text: str) -> str:
        return self._apply(self.MAGENTA, text)


@dataclass
class TargetStatus:
    """Status of a single docking target (PDB code)."""
    name: str
    state: str = "queued"
    has_result: bool = False
    is_stuck: bool = False


@dataclass
class RunProgress:
    """Progress for a single run directory."""
    run_name: str
    done: int = 0
    in_progress: int = 0
    stuck: int = 0
    queued: int = 0
    total: int = 0
    targets: List[TargetStatus] = field(default_factory=list)


@dataclass
class DatasetProgress:
    """Progress for a dataset across all runs."""
    name: str
    total_runs: int = 0
    completed_runs: int = 0
    status: str = "pending"
    active_run: Optional[RunProgress] = None


@dataclass
class RunnerInfo:
    """Parsed information for one fleet runner."""
    name: str
    source_file: str
    active_dataset: str = "none"
    datasets: Dict[str, DatasetProgress] = field(default_factory=OrderedDict)
    active_workers: int = 0
    failed_chunks: int = 0
    eta_seconds: Optional[float] = None
    timestamp: Optional[str] = None
    raw_json: Optional[Dict[str, Any]] = None


@dataclass
class SystemStats:
    """System resource utilization."""
    cpu_percent: float = 0.0
    memory_pressure: str = "normal"
    thermal_state: str = "nominal"
    flexaidds_workers: int = 0


def should_use_color() -> bool:
    """Determine if ANSI colors should be used."""
    if os.environ.get("NO_COLOR"):
        return False
    if not sys.stdout.isatty():
        return False
    return True


def get_terminal_width() -> int:
    """Get terminal width with fallback to 80."""
    try:
        return os.get_terminal_size().columns
    except OSError:
        return 80


def safe_read_json(filepath: str) -> Optional[Dict[str, Any]]:
    """Read and parse a JSON file, handling partial writes and iCloud eviction.

    Returns None if the file cannot be read or parsed.
    """
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return None


def safe_read_file(filepath: str) -> Optional[str]:
    """Read a text file, handling iCloud eviction and missing files."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            return f.read()
    except OSError:
        return None


def safe_list_dir(path: str) -> List[str]:
    """List directory contents, returning empty list on failure."""
    try:
        return sorted(os.listdir(path))
    except OSError:
        return []


def safe_is_dir(path: str) -> bool:
    """Check if path is a directory, returning False on failure."""
    try:
        return os.path.isdir(path)
    except OSError:
        return False


def safe_is_file(path: str) -> bool:
    """Check if path is a file, returning False on failure."""
    try:
        return os.path.isfile(path)
    except OSError:
        return False


def discover_fleet_status_files(icloud_dir: str) -> List[str]:
    """Find all fleet_status*.json files in the iCloud directory.

    Args:
        icloud_dir: Path to the iCloud FlexAIDdS directory.

    Returns:
        Sorted list of absolute paths to fleet status JSON files.
    """
    pattern = os.path.join(icloud_dir, "fleet_status*.json")
    files = glob(pattern)
    return sorted(files)


def count_flexaidds_workers() -> int:
    """Count active FlexAIDdS worker processes.

    Uses `ps aux` to find processes matching FlexAIDdS/benchmark_datasets.
    Excludes the grep process itself.

    Returns:
        Number of active worker processes.
    """
    try:
        result = subprocess.run(
            ["ps", "aux"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode != 0:
            return 0
        count = 0
        for line in result.stdout.splitlines():
            lower = line.lower()
            if "flexaidds" in lower or "benchmark_datasets" in lower:
                if "grep" not in lower:
                    count += 1
        return count
    except (subprocess.TimeoutExpired, OSError):
        return 0


def get_system_stats() -> SystemStats:
    """Collect system resource statistics on macOS.

    Returns:
        SystemStats with CPU, memory, thermal, and worker info.
    """
    stats = SystemStats()

    try:
        result = subprocess.run(
            ["ps", "-A", "-o", "%cpu"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            total_cpu = 0.0
            for line in result.stdout.strip().splitlines()[1:]:
                try:
                    total_cpu += float(line.strip())
                except ValueError:
                    pass
            stats.cpu_percent = min(total_cpu, 100.0)
    except (subprocess.TimeoutExpired, OSError):
        pass

    try:
        result = subprocess.run(
            ["memory_pressure"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            output = result.stdout.lower()
            if "high" in output:
                stats.memory_pressure = "high"
            elif "medium" in output:
                stats.memory_pressure = "medium"
            else:
                stats.memory_pressure = "normal"
    except (subprocess.TimeoutExpired, OSError):
        pass

    try:
        thermal_path = "/Library/Application Support/HWMonitor/thermal_state"
        if safe_is_file(thermal_path):
            content = safe_read_file(thermal_path)
            if content:
                stats.thermal_state = content.strip().lower()
        else:
            stats.thermal_state = "nominal"
    except OSError:
        pass

    stats.flexaidds_workers = count_flexaidds_workers()
    return stats


def scan_target_dir(target_path: str) -> TargetStatus:
    """Determine the status of a single target directory.

    A target is:
      - "done" if it has a result.csv or non-INI .pdb output file
      - "stuck" if stdout.log contains [STUCK]
      - "in_progress" if it has dock_config.json but no result yet
      - "queued" if the directory doesn't exist yet

    Args:
        target_path: Absolute path to the target directory.

    Returns:
        TargetStatus with the determined state.
    """
    target_name = os.path.basename(target_path)
    status = TargetStatus(name=target_name)

    if not safe_is_dir(target_path):
        status.state = "queued"
        return status

    entries = safe_list_dir(target_path)

    has_result_csv = "result.csv" in entries
    has_dock_config = "dock_config.json" in entries
    has_stdout = "stdout.log" in entries

    has_output_pdb = False
    for entry in entries:
        if entry.endswith(".pdb") and not entry.endswith("_INI.pdb"):
            has_output_pdb = True
            break

    status.has_result = has_result_csv or has_output_pdb

    if has_stdout:
        stdout_path = os.path.join(target_path, "stdout.log")
        content = safe_read_file(stdout_path)
        if content and "[STUCK]" in content:
            status.is_stuck = True
            status.state = "stuck"
            return status

    if has_result_csv or has_output_pdb:
        status.state = "done"
    elif has_dock_config:
        status.state = "in_progress"
    else:
        status.state = "queued"

    return status


def scan_run_dir(run_path: str, known_total: Optional[int] = None) -> RunProgress:
    """Deep-scan a run directory for target-level progress.

    Args:
        run_path: Absolute path to the run directory (e.g., .../astex/run31/).
        known_total: Expected total number of targets (from fleet status).

    Returns:
        RunProgress with per-target breakdown.
    """
    run_name = os.path.basename(run_path)
    progress = RunProgress(run_name=run_name)

    entries = safe_list_dir(run_path)
    target_dirs = [
        e for e in entries
        if safe_is_dir(os.path.join(run_path, e))
        and not e.startswith(".")
        and e not in ("logs", "tmp")
    ]

    if not target_dirs and known_total and known_total > 0:
        progress.total = known_total
        progress.queued = known_total
        return progress

    for td in target_dirs:
        target_path = os.path.join(run_path, td)
        ts = scan_target_dir(target_path)
        progress.targets.append(ts)

        if ts.state == "done":
            progress.done += 1
        elif ts.state == "stuck":
            progress.stuck += 1
        elif ts.state == "in_progress":
            progress.in_progress += 1
        else:
            progress.queued += 1

    progress.total = known_total if known_total and known_total > len(target_dirs) else len(target_dirs)
    if progress.total < len(target_dirs):
        progress.total = len(target_dirs)

    return progress


def find_active_run(results_base: str, dataset: str) -> Optional[str]:
    """Find the most recently modified run directory for a dataset.

    Args:
        results_base: Base results path (e.g., ~/.flexaidds_fast/results/tier2).
        dataset: Dataset name (e.g., 'astex').

    Returns:
        Path to the most recent run directory, or None.
    """
    ds_dir = os.path.join(results_base, dataset)
    if not safe_is_dir(ds_dir):
        return None

    run_dirs = []
    for entry in safe_list_dir(ds_dir):
        full = os.path.join(ds_dir, entry)
        if safe_is_dir(full) and entry.startswith("run"):
            try:
                mtime = os.path.getmtime(full)
                run_dirs.append((mtime, full))
            except OSError:
                continue

    if not run_dirs:
        return None

    run_dirs.sort(key=lambda x: x[0], reverse=True)
    return run_dirs[0][1]


def estimate_dataset_total(results_dirs: List[str], dataset: str) -> int:
    """Estimate total targets for a dataset by scanning known runs.

    Looks across all result directories for the maximum number of
    targets seen in any single run for this dataset.

    Args:
        results_dirs: List of results base directories to scan.
        dataset: Dataset name.

    Returns:
        Estimated total number of targets (0 if unknown).
    """
    max_targets = 0
    for base in results_dirs:
        ds_dir = os.path.join(base, dataset)
        if not safe_is_dir(ds_dir):
            continue
        for run_entry in safe_list_dir(ds_dir):
            run_path = os.path.join(ds_dir, run_entry)
            if not safe_is_dir(run_path):
                continue
            count = len([
                e for e in safe_list_dir(run_path)
                if safe_is_dir(os.path.join(run_path, e))
                and not e.startswith(".")
            ])
            if count > max_targets:
                max_targets = count
    return max_targets


def parse_runner(
    filepath: str,
    results_dirs: List[str],
    runner_results_dir: Optional[str] = None,
) -> RunnerInfo:
    """Parse a fleet status JSON file and deep-scan result directories.

    Args:
        filepath: Path to the fleet_status*.json file.
        results_dirs: List of base result directories to scan.
        runner_results_dir: Specific results directory for this runner.

    Returns:
        RunnerInfo with parsed status and deep-scanned progress.
    """
    runner = RunnerInfo(
        name="unknown",
        source_file=filepath,
    )

    data = safe_read_json(filepath)
    if data is None:
        runner.name = _runner_name_from_file(filepath)
        return runner

    runner.raw_json = data
    runner.name = data.get("runner", _runner_name_from_file(filepath))
    runner.active_dataset = data.get("activeDataset", "none")
    runner.timestamp = data.get("timestamp")

    metrics = data.get("metrics", {})
    runner.failed_chunks = metrics.get("failedChunks", 0)
    eta_s = metrics.get("estimatedRemainingSeconds")
    if eta_s is not None:
        try:
            runner.eta_seconds = float(eta_s)
        except (ValueError, TypeError):
            pass

    campaign = data.get("campaign", {})
    for ds_name, ds_info in campaign.items():
        total = ds_info.get("total", 0)
        completed = ds_info.get("completed", 0)
        status = ds_info.get("status", "pending")

        ds_progress = DatasetProgress(
            name=ds_name,
            total_runs=total,
            completed_runs=completed,
            status=status,
        )

        scan_dir = runner_results_dir or (
            results_dirs[0] if results_dirs else None
        )

        if scan_dir and status in ("running", "complete"):
            tier2_base = os.path.join(scan_dir, RESULTS_SUBDIR) if not scan_dir.endswith(RESULTS_SUBDIR) else scan_dir
            active_run_path = find_active_run(tier2_base, ds_name)
            if active_run_path:
                known_total = estimate_dataset_total(
                    [tier2_base for _ in results_dirs] if not results_dirs else
                    [os.path.join(r, RESULTS_SUBDIR) if not r.endswith(RESULTS_SUBDIR) else r for r in results_dirs],
                    ds_name,
                )
                ds_progress.active_run = scan_run_dir(active_run_path, known_total)

        runner.datasets[ds_name] = ds_progress

    if runner_results_dir:
        tier2 = os.path.join(runner_results_dir, RESULTS_SUBDIR)
        if safe_is_dir(tier2):
            for ds_entry in safe_list_dir(tier2):
                ds_path = os.path.join(tier2, ds_entry)
                if not safe_is_dir(ds_path):
                    continue
                if ds_entry not in runner.datasets:
                    runs = [
                        e for e in safe_list_dir(ds_path)
                        if safe_is_dir(os.path.join(ds_path, e)) and e.startswith("run")
                    ]
                    if runs:
                        runner.datasets[ds_entry] = DatasetProgress(
                            name=ds_entry,
                            total_runs=len(runs),
                            completed_runs=0,
                            status="discovered",
                        )
                        active = find_active_run(tier2, ds_entry)
                        if active:
                            known_total = estimate_dataset_total(
                                [os.path.join(r, RESULTS_SUBDIR) if not r.endswith(RESULTS_SUBDIR) else r for r in results_dirs],
                                ds_entry,
                            )
                            runner.datasets[ds_entry].active_run = scan_run_dir(active, known_total)

    runner.active_workers = count_flexaidds_workers()
    return runner


def _runner_name_from_file(filepath: str) -> str:
    """Derive a runner name from the fleet status filename."""
    basename = os.path.basename(filepath)
    name = basename.replace("fleet_status", "").replace(".json", "").strip("_")
    if not name:
        return "unknown"
    return name


def format_eta(seconds: Optional[float]) -> str:
    """Format seconds into a human-readable ETA string.

    Args:
        seconds: Number of seconds, or None.

    Returns:
        Formatted string like '~2h 15m' or '---'.
    """
    if seconds is None or seconds <= 0:
        return "---"
    s = int(seconds)
    hours, remainder = divmod(s, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours > 0:
        return f"~{hours}h {minutes}m"
    elif minutes > 0:
        return f"~{minutes}m {secs}s"
    else:
        return f"~{secs}s"


def make_progress_bar(done: int, total: int, width: int = 20) -> str:
    """Create a text progress bar.

    Args:
        done: Completed count.
        total: Total count.
        width: Width of the bar in characters.

    Returns:
        String like '████████░░░░░░░░░░░░'
    """
    if total <= 0:
        return PROGRESS_EMPTY * width
    filled = int(width * done / total)
    filled = min(filled, width)
    return PROGRESS_FILLED * filled + PROGRESS_EMPTY * (width - filled)


def format_percent(done: int, total: int) -> str:
    """Format a percentage string."""
    if total <= 0:
        return "??%"
    return f"{100 * done // total}%"


def collect_runners(
    icloud_dir: str,
    results_dir: Optional[str] = None,
) -> List[RunnerInfo]:
    """Discover and parse all fleet runners.

    Determines runner identity from filename and maps to appropriate
    result directories:
      - opencode → ~/.flexaidds_fast (local SSD)
      - claude-code (or default fleet_status.json) → iCloud results

    Args:
        icloud_dir: Path to iCloud FlexAIDdS directory.
        results_dir: Override results directory.

    Returns:
        List of RunnerInfo objects.
    """
    fleet_files = discover_fleet_status_files(icloud_dir)
    fast_base = os.environ.get("FLEXAIDDS_FAST_BASE", FAST_DEFAULT)
    icloud_results = os.path.join(icloud_dir, "results")

    all_results_dirs = []
    if results_dir:
        all_results_dirs.append(results_dir)
    else:
        if os.path.isdir(os.path.join(fast_base, RESULTS_SUBDIR)):
            all_results_dirs.append(fast_base)
        if os.path.isdir(icloud_results):
            all_results_dirs.append(icloud_results)

    runners: List[RunnerInfo] = []

    for fpath in fleet_files:
        basename = os.path.basename(fpath)
        runner_results = None

        if "opencode" in basename.lower():
            runner_results = results_dir or fast_base
        else:
            runner_results = results_dir or icloud_results

        runner = parse_runner(fpath, all_results_dirs, runner_results)
        runners.append(runner)

    if not fleet_files and all_results_dirs:
        fallback_runner = RunnerInfo(
            name="scanner",
            source_file="(direct scan)",
        )
        for base in all_results_dirs:
            tier2 = os.path.join(base, RESULTS_SUBDIR)
            if not safe_is_dir(tier2):
                continue
            for ds_entry in safe_list_dir(tier2):
                ds_path = os.path.join(tier2, ds_entry)
                if not safe_is_dir(ds_path):
                    continue
                runs = [
                    e for e in safe_list_dir(ds_path)
                    if safe_is_dir(os.path.join(ds_path, e)) and e.startswith("run")
                ]
                if ds_entry not in fallback_runner.datasets:
                    active = find_active_run(tier2, ds_entry)
                    run_prog = None
                    if active:
                        known_total = estimate_dataset_total(
                            [os.path.join(r, RESULTS_SUBDIR) if not r.endswith(RESULTS_SUBDIR) else r for r in all_results_dirs],
                            ds_entry,
                        )
                        run_prog = scan_run_dir(active, known_total)
                    fallback_runner.datasets[ds_entry] = DatasetProgress(
                        name=ds_entry,
                        total_runs=len(runs),
                        status="discovered",
                        active_run=run_prog,
                    )
        fallback_runner.active_workers = count_flexaidds_workers()
        if fallback_runner.datasets:
            runners.append(fallback_runner)

    return runners


def render_dashboard(
    runners: List[RunnerInfo],
    system_stats: SystemStats,
    colors: Colors,
    width: Optional[int] = None,
) -> str:
    """Render the full terminal dashboard as a string.

    Args:
        runners: List of parsed runner information.
        system_stats: System resource statistics.
        colors: Colors instance for ANSI formatting.
        width: Terminal width (auto-detected if None).

    Returns:
        Multi-line string with the complete dashboard.
    """
    if width is None:
        width = get_terminal_width()
    width = max(width, 60)

    inner_width = width - 4
    lines: List[str] = []

    lines.append(f"{BOX_TL}{BOX_HORIZONTAL * inner_width}{BOX_TR}")

    title = " Bonhomme Fleet Dashboard \u2014 FlexAIDdS Benchmark Campaign "
    padded_title = title.center(inner_width)
    lines.append(f"{BOX_VERTICAL}{padded_title}{BOX_VERTICAL}")

    now = datetime.now().astimezone()
    tz_name = now.strftime("%Z") or "Local"
    ts_str = now.strftime(f"%Y-%m-%d %H:%M:%S {tz_name}")
    padded_ts = f" Updated: {ts_str} ".center(inner_width)
    lines.append(f"{BOX_VERTICAL}{padded_ts}{BOX_VERTICAL}")

    lines.append(f"{BOX_LJ}{BOX_HORIZONTAL * inner_width}{BOX_RJ}")

    for idx, runner in enumerate(runners):
        runner_name_raw = runner.name
        if "opencode" in runner_name_raw.lower() or "glm" in runner_name_raw.lower():
            runner_display = "opencode (GLM-5.1)"
            runner_line = f"  Runner: {runner_display}"
        elif "claude" in runner_name_raw.lower():
            runner_display = "claude-code (Opus 4.7)"
            runner_line = f"  Runner: {runner_display}"
        else:
            runner_display = runner_name_raw
            runner_line = f"  Runner: {runner_display}"

        runner_line = runner_line[:inner_width].ljust(inner_width)
        if colors._enabled:
            colored_line = runner_line.replace(f"Runner: {runner_display}", f"Runner: {colors.cyan(runner_display) if 'opencode' in runner_name_raw.lower() or 'glm' in runner_name_raw.lower() else colors.magenta(runner_display) if 'claude' in runner_name_raw.lower() else colors.bold(runner_display)}")
        else:
            colored_line = runner_line

        lines.append(f"{BOX_VERTICAL}{' ' * inner_width}{BOX_VERTICAL}")
        lines.append(f"{BOX_VERTICAL}{colored_line}{BOX_VERTICAL}")

        box_inner = inner_width - 6

        lines.append(f"{BOX_VERTICAL}  \u250c{'\u2500' * box_inner}\u2510  {BOX_VERTICAL}")

        if not runner.datasets:
            empty_msg = "  No datasets configured  ".center(box_inner)
            lines.append(f"{BOX_VERTICAL}  \u2502{empty_msg}\u2502  {BOX_VERTICAL}")
        else:
            for ds_name, ds_prog in runner.datasets.items():
                line_parts = _render_dataset_line(
                    ds_name, ds_prog, colors, box_inner
                )
                lines.append(f"{BOX_VERTICAL}  \u2502{line_parts}\u2502  {BOX_VERTICAL}")

        lines.append(f"{BOX_VERTICAL}  \u2514{'\u2500' * box_inner}\u2518  {BOX_VERTICAL}")

        workers = runner.active_workers
        failed = runner.failed_chunks
        eta = format_eta(runner.eta_seconds)

        stats_plain = f"  Workers: {workers} active  \u2502  Failed: {failed}  \u2502  ETA: {eta}"
        stats_plain = stats_plain[:inner_width].ljust(inner_width)

        if colors._enabled:
            w_str = colors.green(str(workers)) if workers > 0 else colors.dim(str(workers))
            f_str = colors.red(str(failed)) if failed > 0 else colors.green("0")
            stats_colored = f"  Workers: {w_str} active  \u2502  Failed: {f_str}  \u2502  ETA: {eta}"
            stats_colored = _pad_ansi(stats_colored, inner_width)
        else:
            stats_colored = stats_plain

        lines.append(f"{BOX_VERTICAL}{stats_colored}{BOX_VERTICAL}")

    sys_cpu = f"CPU: {system_stats.cpu_percent:.0f}%"
    mem_str = f"Mem: {system_stats.memory_pressure}"
    therm_str = f"Thermal: {system_stats.thermal_state}"
    workers_str = f"FlexAIDdS procs: {system_stats.flexaidds_workers}"

    sys_line = f"  System: {sys_cpu}  |  {mem_str}  |  {therm_str}  |  {workers_str}  "
    padded_sys = sys_line[:inner_width].ljust(inner_width)
    lines.append(f"{BOX_VERTICAL}{' ' * inner_width}{BOX_VERTICAL}")
    lines.append(f"{BOX_VERTICAL}{padded_sys}{BOX_VERTICAL}")

    lines.append(f"{BOX_VERTICAL}{'':^{inner_width}}{BOX_VERTICAL}")
    lines.append(f"{BOX_BL}{BOX_HORIZONTAL * inner_width}{BOX_BR}")

    return "\n".join(lines)


def _visible_len(text: str) -> int:
    """Return the visible (display) width of a string, ignoring ANSI escapes."""
    import re as _re
    stripped = _re.sub(r'\033\[[0-9;]*m', '', text)
    return len(stripped)


def _pad_ansi(text: str, width: int) -> str:
    """Right-pad an ANSI-colored string to a given visible width."""
    vis = _visible_len(text)
    if vis >= width:
        return text
    return text + ' ' * (width - vis)


def _render_dataset_line(
    ds_name: str,
    ds_prog: DatasetProgress,
    colors: Colors,
    box_inner: int,
) -> str:
    """Render a single dataset progress line.

    Args:
        ds_name: Dataset name (e.g., 'astex').
        ds_prog: DatasetProgress with run data.
        colors: Colors instance.
        box_inner: Inner width of the box.

    Returns:
        Formatted line string fitting within box_inner.
    """
    run_info = ds_prog.active_run
    if run_info and run_info.total > 0:
        done = run_info.done
        total = run_info.total
        active_run_name = run_info.run_name
    else:
        done = ds_prog.completed_runs
        total = ds_prog.total_runs
        active_run_name = "---"

    bar_width = 20
    bar = make_progress_bar(done, total, bar_width)
    pct = format_percent(done, total)

    if ds_prog.status == "running":
        bar = colors.green(bar)
        pct_str = colors.green(pct)
    elif ds_prog.status == "complete":
        bar = colors.bold(colors.green(bar))
        pct_str = colors.bold(colors.green(pct))
    elif ds_prog.status == "pending":
        bar = colors.dim(bar)
        pct_str = colors.dim(pct)
    else:
        pct_str = pct

    stuck_str = ""
    if run_info and run_info.stuck > 0:
        stuck_str = colors.red(f" stuck:{run_info.stuck}")

    plain_bar = make_progress_bar(done, total, bar_width)
    label_plain = f" {ds_name:<16s} {plain_bar} {active_run_name:>5s}  {done:>3}/{total:<3} {pct}{stuck_str.replace(colors.red(''), '').replace(colors._apply(colors.RED, ''), '')}"

    label_colored = f" {ds_name:<16s} {bar} {active_run_name:>5s}  {done:>3}/{total:<3} {pct_str}{stuck_str}"

    if _visible_len(label_plain) > box_inner:
        trimmed = label_plain[:box_inner]
        label_colored = trimmed
    else:
        label_colored = _pad_ansi(label_colored, box_inner)

    return label_colored


def render_json(runners: List[RunnerInfo], system_stats: SystemStats) -> str:
    """Render machine-readable JSON output.

    Args:
        runners: List of parsed runner information.
        system_stats: System resource statistics.

    Returns:
        JSON string with merged status.
    """
    output = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "system": {
            "cpu_percent": system_stats.cpu_percent,
            "memory_pressure": system_stats.memory_pressure,
            "thermal_state": system_stats.thermal_state,
            "flexaidds_workers": system_stats.flexaidds_workers,
        },
        "runners": [],
    }

    for runner in runners:
        r_data: Dict[str, Any] = {
            "name": runner.name,
            "source_file": runner.source_file,
            "active_dataset": runner.active_dataset,
            "active_workers": runner.active_workers,
            "failed_chunks": runner.failed_chunks,
            "eta_seconds": runner.eta_seconds,
            "timestamp": runner.timestamp,
            "datasets": {},
        }

        for ds_name, ds_prog in runner.datasets.items():
            ds_data: Dict[str, Any] = {
                "total_runs": ds_prog.total_runs,
                "completed_runs": ds_prog.completed_runs,
                "status": ds_prog.status,
            }
            if ds_prog.active_run:
                run = ds_prog.active_run
                ds_data["active_run"] = {
                    "name": run.run_name,
                    "done": run.done,
                    "in_progress": run.in_progress,
                    "stuck": run.stuck,
                    "queued": run.queued,
                    "total": run.total,
                }
            r_data["datasets"][ds_name] = ds_data

        output["runners"].append(r_data)

    return json.dumps(output, indent=2)


def build_argparser() -> argparse.ArgumentParser:
    """Build the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="fleet_monitor",
        description="Bonhomme Fleet Dashboard — FlexAIDdS Benchmark Campaign Monitor",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=30,
        help="Poll interval in seconds (default: 30)",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Single snapshot, then exit",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Machine-readable JSON output (for piping)",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable ANSI colors",
    )
    parser.add_argument(
        "--results",
        type=str,
        default=None,
        help="Override results directory",
    )
    parser.add_argument(
        "--icloud",
        type=str,
        default=None,
        help="Override iCloud status directory",
    )
    return parser


def run_single(
    icloud_dir: str,
    results_dir: Optional[str],
    json_output: bool,
    colors: Colors,
) -> str:
    """Collect a single snapshot and return the rendered output.

    Args:
        icloud_dir: iCloud FlexAIDdS directory.
        results_dir: Override results directory.
        json_output: If True, output JSON instead of dashboard.
        colors: Colors instance.

    Returns:
        Rendered string (dashboard or JSON).
    """
    runners = collect_runners(icloud_dir, results_dir)
    system_stats = get_system_stats()

    if json_output:
        return render_json(runners, system_stats)
    else:
        return render_dashboard(runners, system_stats, colors)


def main(argv: Optional[List[str]] = None) -> int:
    """Main entry point for the fleet monitor.

    Args:
        argv: Command-line arguments (defaults to sys.argv[1:]).

    Returns:
        Exit code (0 for success).
    """
    parser = build_argparser()
    args = parser.parse_args(argv)

    icloud_dir = args.icloud or os.environ.get(
        "FLEXAIDDS_ICLOUD", ICLOUD_DEFAULT
    )
    results_dir = args.results

    use_color = not args.no_color and should_use_color()
    colors = Colors(enabled=use_color)

    running = True

    def handle_sigint(signum: int, frame: Any) -> None:
        nonlocal running
        running = False

    signal.signal(signal.SIGINT, handle_sigint)

    first_iteration = True

    while running:
        try:
            output = run_single(icloud_dir, results_dir, args.json_output, colors)

            if not args.json_output and not first_iteration:
                num_lines = output.count("\n") + 1
                sys.stdout.write(f"\033[{num_lines}A")
            first_iteration = False

            sys.stdout.write(output + "\n")
            sys.stdout.flush()

        except Exception as exc:
            sys.stderr.write(f"Error collecting snapshot: {exc}\n")
            sys.stderr.flush()

        if args.once:
            break

        try:
            for _ in range(args.interval):
                if not running:
                    break
                time.sleep(1)
        except KeyboardInterrupt:
            running = False

    if not args.json_output and use_color:
        sys.stdout.write("\033[?25h")
        sys.stdout.flush()

    return 0


if __name__ == "__main__":
    sys.exit(main())
