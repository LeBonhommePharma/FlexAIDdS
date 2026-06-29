#!/usr/bin/env python3
# queue_v132_ablation_ladder.py — wait for v131, then run v132 ablation ladder serially
#
# Watches v131 full-85 + smoke-12 (and foreign benchmark_datasets) until quiet,
# applies settle delay, then launches each ablation step in order:
#   consensus_on → safe_binary → logsumexp_only → hbond_zero
#
# Usage:
#   python3 scripts/queue_v132_ablation_ladder.py --daemon
#   python3 scripts/queue_v132_ablation_ladder.py --poll 120 --settle 120
#   python3 scripts/queue_v132_ablation_ladder.py --steps consensus_on safe_binary
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from lib_launch import launch_session_isolated
from lib_v132_ablation import RESULTS, ladder_step_ids, step_by_id

REPO = SCRIPT_DIR.parent
LAUNCH_SCRIPT = SCRIPT_DIR / "launch_v132_ablation.py"
QUEUE_LOG = RESULTS / "queue_v132_ablation_ladder.log"
STATE_FILE = RESULTS / "queue_v132_ablation_ladder.state.json"

V131_FULL85_GLOB = "v131_*_r07_nofixb_full85"
V131_SMOKE_GLOB = "v131_*_smoke12_safe"
DEFAULT_STEPS = ladder_step_ids()


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


def latest_run_dir(glob_pat: str) -> Path | None:
    dirs = sorted(RESULTS.glob(glob_pat), key=lambda p: p.stat().st_mtime)
    return dirs[-1] if dirs else None


def active_foreign_benchmarks() -> list[str]:
    try:
        ps = subprocess.check_output(["pgrep", "-fl", "benchmark_datasets"], text=True)
    except subprocess.CalledProcessError:
        return []
    lines = []
    for line in ps.strip().splitlines():
        if "caffeinate" in line:
            continue
        if "v132_" in line or "queue_v132" in line:
            continue
        lines.append(line.strip())
    return lines


def summarize_full85(run_dir: Path) -> dict:
    total = 85
    if not run_dir.is_dir():
        return {"run_dir": str(run_dir), "completed": 0, "total": total, "success": 0}
    done = [d.name for d in run_dir.iterdir() if d.is_dir() and (d / "result.csv").exists()]
    succ = 0
    for name in done:
        try:
            row = next(csv.DictReader((run_dir / name / "result.csv").open()))
            if row.get("success") == "1":
                succ += 1
        except (StopIteration, OSError):
            pass
    return {
        "run_dir": run_dir.name,
        "completed": len(done),
        "total": total,
        "success": succ,
    }


def launch_step(step_id: str, skip_build: bool) -> tuple[int, Path | None]:
    log(f"Launching v132 step={step_id}")
    cmd = [sys.executable, str(LAUNCH_SCRIPT), step_id]
    if skip_build:
        cmd.append("--skip-build")
    proc = subprocess.run(cmd, cwd=str(REPO), capture_output=True, text=True)
    if proc.stdout.strip():
        log(proc.stdout.strip())
    if proc.stderr.strip():
        log(f"step={step_id} stderr: {proc.stderr.strip()}")
    if proc.returncode != 0:
        log(f"ERROR: step={step_id} launch failed exit={proc.returncode}")
        return -1, None

    dirs = sorted(RESULTS.glob(f"v132_*_{step_id}_full85"), key=lambda p: p.stat().st_mtime)
    out_dir = dirs[-1] if dirs else None
    pid = 0
    if out_dir and (out_dir / "launch_provenance.json").is_file():
        pid = int(json.loads((out_dir / "launch_provenance.json").read_text()).get("pid") or 0)
    log(f"step={step_id} queued pid={pid} output={out_dir.name if out_dir else '?'}")
    return pid, out_dir


def wait_for_run(pid: int, out_dir: Path | None, poll_s: int) -> bool:
    while pid_alive(pid):
        if out_dir and out_dir.is_dir():
            s = summarize_full85(out_dir)
            log(f"progress {out_dir.name}: {s['completed']}/{s['total']} ({s['success']} success)")
        time.sleep(poll_s)
    if out_dir and out_dir.is_dir():
        s = summarize_full85(out_dir)
        log(
            f"finished {out_dir.name}: {s['completed']}/{s['total']} "
            f"success={s['success']}"
        )
        return s["completed"] == s["total"]
    return True


def write_ladder_summary(results: list[dict]) -> Path:
    report = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "ladder": "v132_oracle_ablation",
        "steps": results,
        "audit_quote": (
            "Highest honest oracle next: restore consensus ON (step consensus_on), "
            "then ablate binary safe / logsumexp / hbond one knob at a time."
        ),
    }
    path = RESULTS / "v132_ablation_ladder_summary.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    log(f"Wrote ladder summary: {path}")
    return path


def run_watcher(
    poll_s: int,
    settle_s: int,
    steps: list[str],
    skip_build: bool,
    no_launch: bool,
    wait_v131_smoke: bool,
) -> int:
    v131_dir = latest_run_dir(V131_FULL85_GLOB)
    smoke_dir = latest_run_dir(V131_SMOKE_GLOB) if wait_v131_smoke else None
    v131_pid = discover_pid_from_dir(v131_dir) if v131_dir else 0
    smoke_pid = discover_pid_from_dir(smoke_dir) if smoke_dir else 0

    save_state(
        status="watching",
        poll_s=poll_s,
        settle_s=settle_s,
        steps=steps,
        v131_dir=v131_dir.name if v131_dir else None,
        v131_pid=v131_pid,
        smoke_dir=smoke_dir.name if smoke_dir else None,
        smoke_pid=smoke_pid,
    )
    log("v132 ablation ladder watcher started")
    log(f"  v131_full85={v131_dir.name if v131_dir else 'none'} pid={v131_pid or 'done'}")
    if wait_v131_smoke:
        log(f"  v131_smoke12={smoke_dir.name if smoke_dir else 'none'} pid={smoke_pid or 'done'}")
    log(f"  planned_steps={steps}")

    while True:
        v131_pid = discover_pid_from_dir(v131_dir) if v131_dir and pid_alive(v131_pid) else 0
        if not v131_pid and v131_dir:
            v131_pid = discover_pid_from_dir(v131_dir)
        smoke_pid = 0
        if wait_v131_smoke and smoke_dir:
            smoke_pid = discover_pid_from_dir(smoke_dir) if pid_alive(smoke_pid) else 0
            if not smoke_pid:
                smoke_pid = discover_pid_from_dir(smoke_dir)

        foreign = active_foreign_benchmarks()
        busy = bool(v131_pid or smoke_pid or foreign)
        if busy:
            parts = []
            if v131_dir:
                if v131_pid:
                    s = summarize_full85(v131_dir)
                    parts.append(f"v131 {s['completed']}/{s['total']}")
                else:
                    parts.append("v131 done")
            if wait_v131_smoke and smoke_dir:
                parts.append(f"smoke12 pid={smoke_pid or 'done'}")
            if foreign:
                parts.append(f"foreign_benchmarks={len(foreign)}")
            log("waiting: " + ", ".join(parts))
            time.sleep(poll_s)
            continue
        break

    log(f"Queue quiet — settling {settle_s}s before v132 ladder")
    save_state(status="settling", settle_s=settle_s)
    time.sleep(settle_s)

    if active_foreign_benchmarks():
        log("WARN: foreign benchmarks appeared during settle; restarting watch")
        return run_watcher(poll_s, settle_s, steps, skip_build, no_launch, wait_v131_smoke)

    if no_launch:
        log("Skipping ladder launch (--no-launch)")
        save_state(status="queue_clear_no_launch")
        return 0

    ladder_results: list[dict] = []
    for step_id in steps:
        step_by_id(step_id)  # validate
        save_state(status="launching", current_step=step_id)
        pid, out_dir = launch_step(step_id, skip_build)
        if pid <= 0:
            save_state(status="launch_failed", failed_step=step_id)
            return 1
        ok = wait_for_run(pid, out_dir, poll_s)
        summary = summarize_full85(out_dir) if out_dir else {}
        ladder_results.append(
            {
                "step_id": step_id,
                "output_dir": str(out_dir) if out_dir else None,
                "complete": ok,
                **summary,
            }
        )
        save_state(status="step_done", last_step=step_id, last_summary=summary)

    write_ladder_summary(ladder_results)
    save_state(status="ladder_done", ladder_results=ladder_results)
    log("v132 ablation ladder watcher exiting")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poll", type=int, default=120)
    parser.add_argument("--settle", type=int, default=120)
    parser.add_argument(
        "--steps",
        nargs="+",
        default=DEFAULT_STEPS,
        help=f"Ordered step ids (default: {' '.join(DEFAULT_STEPS)})",
    )
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--no-launch", action="store_true")
    parser.add_argument(
        "--no-wait-smoke12",
        action="store_true",
        help="Only wait for v131 full-85, not smoke-12",
    )
    parser.add_argument("--daemon", action="store_true")
    args = parser.parse_args()

    for step_id in args.steps:
        step_by_id(step_id)

    if args.daemon:
        watcher_log = RESULTS / "queue_v132_ablation_ladder_watcher"
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
            "--steps",
            *args.steps,
        ]
        if args.skip_build:
            cmd.append("--skip-build")
        if args.no_launch:
            cmd.append("--no-launch")
        if args.no_wait_smoke12:
            cmd.append("--no-wait-smoke12")
        pid = launch_session_isolated(
            cmd,
            os.environ.copy(),
            str(watcher_log),
            cwd=str(REPO),
            stdout_log=str(watcher_log / "watcher_stdout.log"),
            stderr_log=str(watcher_log / "watcher_stderr.log"),
        )
        save_state(status="watcher_daemon_started", watcher_pid=pid, steps=args.steps)
        print(f"v132 ladder watcher daemon pid={pid}")
        print(f"Log: {QUEUE_LOG}")
        print(f"State: {STATE_FILE}")
        return 0

    return run_watcher(
        args.poll,
        args.settle,
        args.steps,
        args.skip_build,
        args.no_launch,
        not args.no_wait_smoke12,
    )


if __name__ == "__main__":
    raise SystemExit(main())