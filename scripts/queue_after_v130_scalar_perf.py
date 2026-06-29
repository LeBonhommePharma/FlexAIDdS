#!/usr/bin/env python3
# queue_after_v130_scalar_perf.py — wait for v130 (+ v124), then scalar-only perf
#
# Watches big Astex campaigns until the computational queue is quiet, applies
# an optional settle delay, then launches launch_perf_scalar_quiet.py.
#
# Usage:
#   python3 scripts/queue_after_v130_scalar_perf.py [--poll 120] [--settle 90]
#
# Daemon (recommended):
#   python3 scripts/queue_after_v130_scalar_perf.py --daemon
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = Path(__file__).resolve().parents[1]
RESULTS = Path(
    os.environ.get(
        "FLEXAIDDS_RESULTS_ROOT",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results",
    )
)
PERF_ROOT = Path(
    os.environ.get(
        "FLEXAIDDS_PERF_VALIDATION_ROOT",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/perf_swarm_validation",
    )
)
QUEUE_LOG = RESULTS / "queue_after_v130_scalar_perf.log"
STATE_FILE = RESULTS / "queue_after_v130_scalar_perf.state.json"
V130_DIR = os.environ.get(
    "QUEUE_V130_DIR",
    "v130_20260629_0548_sulfo_expB_full85",
)
V124_OUTPUT = os.environ.get(
    "QUEUE_V124_OUTPUT",
    str(RESULTS / "multicleft_full_top3_knownsite_v124_15b536f8/results_resume_missing_198"),
)
LAUNCH_SCRIPT = REPO / "scripts" / "launch_perf_scalar_quiet.py"
TIMING_RE = re.compile(
    r"TIMING SUMMARY:\s+\d+\s+gens timed,\s+avg\s+([\d.]+)\s+ms/gen"
)


def log(msg: str) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(QUEUE_LOG, "a") as f:
        f.write(line + "\n")


def save_state(**kwargs) -> None:
    state = {}
    if STATE_FILE.exists():
        state = json.loads(STATE_FILE.read_text())
    state.update(kwargs)
    state["updated_at"] = datetime.now(timezone.utc).isoformat()
    STATE_FILE.write_text(json.dumps(state, indent=2) + "\n")


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def discover_pid_from_dir(run_dir: Path) -> int:
    if not run_dir.is_dir():
        return 0
    for name in ("benchmark.pid",):
        pf = run_dir / name
        if pf.is_file():
            try:
                pid = int(pf.read_text().strip().split()[0])
                if pid_alive(pid):
                    return pid
            except ValueError:
                pass
    prov = run_dir / "launch_provenance.json"
    if prov.is_file():
        try:
            pid = int(json.loads(prov.read_text()).get("pid") or 0)
            if pid_alive(pid):
                return pid
        except (json.JSONDecodeError, ValueError):
            pass
    return 0


def discover_v124_pid() -> int:
    env_pid = int(os.environ.get("QUEUE_V124_PID", "0") or "0")
    if env_pid > 0 and pid_alive(env_pid):
        return env_pid
    try:
        ps = subprocess.check_output(["pgrep", "-fl", "benchmark_datasets"], text=True)
    except subprocess.CalledProcessError:
        return 0
    for line in ps.splitlines():
        if "caffeinate" in line:
            continue
        if "multicleft_full_top3_knownsite_v124" in line or "results_resume_missing_198" in line:
            parts = line.split(None, 1)
            if parts:
                try:
                    return int(parts[0])
                except ValueError:
                    pass
    return 0


def active_foreign_benchmarks() -> list[str]:
    """benchmark_datasets PIDs excluding perf_swarm / this watcher."""
    try:
        ps = subprocess.check_output(["pgrep", "-fl", "benchmark_datasets"], text=True)
    except subprocess.CalledProcessError:
        return []
    lines = []
    for line in ps.strip().splitlines():
        if "caffeinate" in line:
            continue
        if "perf_swarm" in line or "tier1_scalar_quiet" in line or "tier1_paired" in line:
            continue
        lines.append(line.strip())
    return lines


def summarize_v130(run_dir: Path) -> dict:
    json_pairs = Path("/Users/lp.more/Projects/FlexAIDdS/benchmarks/datasets/benchmark_astex_native_85_v130.json")
    total = 85
    if json_pairs.is_file():
        total = len(json.loads(json_pairs.read_text())["pairs"])
    if not run_dir.is_dir():
        return {"run_dir": str(run_dir), "completed": 0, "total": total}
    done = [d.name for d in run_dir.iterdir() if d.is_dir() and (d / "result.csv").exists()]
    succ = 0
    for name in done:
        rf = run_dir / name / "result.csv"
        try:
            row = next(csv.DictReader(rf.open()))
            if row.get("success") == "1":
                succ += 1
        except (StopIteration, OSError):
            pass
    return {
        "run_dir": run_dir.name,
        "completed": len(done),
        "total": total,
        "success_result_csv": succ,
        "missing": max(0, total - len(done)),
    }


def harvest_scalar_timings(out_dir: Path) -> list[dict]:
    rows: list[dict] = []
    label = out_dir / "post_p0_scalar"
    if not label.is_dir():
        return rows
    for log in sorted(label.rglob("stderr.log")):
        try:
            text = log.read_text(errors="replace")
        except OSError:
            continue
        m = TIMING_RE.search(text)
        if m:
            rows.append(
                {
                    "target": log.parent.name,
                    "avg_ms_per_gen": float(m.group(1)),
                }
            )
    return rows


def launch_scalar_perf() -> int:
    log(f"Launching scalar quiet perf via {LAUNCH_SCRIPT}")
    proc = subprocess.run(
        [sys.executable, str(LAUNCH_SCRIPT), "--skip-build"],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if proc.stdout.strip():
        log(proc.stdout.strip())
    if proc.stderr.strip():
        log(f"scalar launch stderr: {proc.stderr.strip()}")
    if proc.returncode != 0:
        log(f"ERROR: scalar launch failed exit={proc.returncode}")
        return -1

    dirs = sorted(PERF_ROOT.glob("tier1_scalar_quiet_*"), key=lambda p: p.stat().st_mtime)
    if dirs:
        prov = dirs[-1] / "launch_provenance.json"
        if prov.is_file():
            pid = int(json.loads(prov.read_text()).get("pid") or 0)
            log(f"scalar perf queued pid={pid} output={dirs[-1].name}")
            return pid
    return 0


def wait_for_perf_completion(perf_pid: int, poll_s: int) -> Path | None:
    out_dirs = sorted(PERF_ROOT.glob("tier1_scalar_quiet_*"), key=lambda p: p.stat().st_mtime)
    out_dir = out_dirs[-1] if out_dirs else None
    while pid_alive(perf_pid):
        if out_dir:
            n = len(list((out_dir / "post_p0_scalar").glob("*/stderr.log"))) if (out_dir / "post_p0_scalar").is_dir() else 0
            log(f"scalar perf progress: {n}/5 targets")
        time.sleep(poll_s)
    return out_dir


def write_final_report(out_dir: Path | None) -> None:
    if not out_dir or not out_dir.is_dir():
        return
    timings = harvest_scalar_timings(out_dir)
    ms = [r["avg_ms_per_gen"] for r in timings]
    median = sorted(ms)[len(ms) // 2] if ms else None
    report = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "output_root": str(out_dir),
        "timings": timings,
        "median_ms_per_gen": median,
        "n_targets": len(timings),
        "complete": len(timings) == 5,
    }
    path = out_dir / "timing_report.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    log(f"Wrote {path} median_ms_per_gen={median}")


def run_watcher(poll_s: int, settle_s: int, also_v124: bool, no_launch: bool) -> int:
    v130_dir = RESULTS / V130_DIR
    v130_pid = discover_pid_from_dir(v130_dir)
    v124_pid = discover_v124_pid() if also_v124 else 0

    save_state(
        status="watching",
        poll_s=poll_s,
        settle_s=settle_s,
        v130_dir=V130_DIR,
        v130_pid=v130_pid,
        v124_pid=v124_pid,
    )
    log("Queue watcher v130→scalar_perf started")
    log(f"  v130_dir={V130_DIR} pid={v130_pid or 'done'}")
    if also_v124:
        log(f"  v124_pid={v124_pid or 'done'}")

    while True:
        v130_pid = discover_pid_from_dir(v130_dir) if pid_alive(v130_pid) else 0
        if also_v124:
            v124_pid = discover_v124_pid() if not v124_pid or pid_alive(v124_pid) else 0
            if v124_pid == 0:
                v124_pid = discover_v124_pid()

        foreign = active_foreign_benchmarks()
        busy = bool(v130_pid or v124_pid or foreign)

        if busy:
            parts = []
            if v130_pid:
                s = summarize_v130(v130_dir)
                parts.append(f"v130 {s['completed']}/{s['total']}")
            elif v130_dir.is_dir():
                parts.append("v130 done")
            if also_v124:
                parts.append(f"v124 pid={v124_pid or 'done'}")
            if foreign:
                parts.append(f"foreign_benchmarks={len(foreign)}")
            log("waiting: " + ", ".join(parts))
            time.sleep(poll_s)
            continue
        break

    log(f"Queue quiet — settling {settle_s}s before scalar perf")
    save_state(status="settling", settle_s=settle_s)
    time.sleep(settle_s)

    foreign = active_foreign_benchmarks()
    if foreign:
        log(f"WARN: foreign benchmarks appeared during settle ({len(foreign)}); waiting again")
        return run_watcher(poll_s, settle_s, also_v124, no_launch)

    log("v130 + queue clear — launching scalar-only tier-1 perf")
    save_state(status="queue_clear")

    perf_pid = -1
    if not no_launch:
        perf_pid = launch_scalar_perf()
        save_state(status="scalar_launched", perf_pid=perf_pid)
        if perf_pid > 0:
            out_dir = wait_for_perf_completion(perf_pid, poll_s)
            write_final_report(out_dir)
            save_state(status="scalar_done", perf_output=str(out_dir) if out_dir else None)
    else:
        log("Skipping scalar launch (--no-launch)")

    log("Queue watcher v130→scalar_perf exiting")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poll", type=int, default=120, help="seconds between checks")
    parser.add_argument(
        "--settle",
        type=int,
        default=90,
        help="seconds to wait after queue clears before launching perf",
    )
    parser.add_argument(
        "--no-v124",
        action="store_true",
        help="Only wait for v130, not v124 resume",
    )
    parser.add_argument("--no-launch", action="store_true")
    parser.add_argument(
        "--daemon",
        action="store_true",
        help="Run this watcher as a detached daemon",
    )
    args = parser.parse_args()

    if args.daemon:
        watcher_log = RESULTS / "queue_after_v130_scalar_perf_watcher"
        watcher_log.mkdir(parents=True, exist_ok=True)
        cmd = [
            "caffeinate",
            "-i",
            sys.executable,
            str(Path(__file__).resolve()),
            "--poll",
            str(args.poll),
            "--settle",
            str(args.settle),
        ]
        if args.no_v124:
            cmd.append("--no-v124")
        if args.no_launch:
            cmd.append("--no-launch")
        pid = launch_session_isolated(
            cmd,
            os.environ.copy(),
            str(watcher_log),
            cwd=str(REPO),
            stdout_log=str(watcher_log / "watcher_stdout.log"),
            stderr_log=str(watcher_log / "watcher_stderr.log"),
        )
        save_state(status="watcher_daemon_started", watcher_pid=pid)
        print(f"Watcher daemon pid={pid}")
        print(f"Log: {QUEUE_LOG}")
        print(f"State: {STATE_FILE}")
        return 0

    return run_watcher(args.poll, args.settle, not args.no_v124, args.no_launch)


if __name__ == "__main__":
    raise SystemExit(main())