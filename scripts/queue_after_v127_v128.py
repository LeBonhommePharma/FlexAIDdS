#!/usr/bin/env python3
# queue_after_v127_v128.py — wait for v127, summarize, launch v128 v50b repro
#
# Chains after queue_after_v124_v126.py (which launches v127 when v126 exits).
# Optionally pre-builds the efc4f5d worktree while v127 is still running.
#
# Usage:
#   python3 scripts/queue_after_v127_v128.py [--poll 120] [--no-prebuild]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path("/Users/lp.more/Projects/FlexAIDdS")
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
QUEUE_LOG = RESULTS / "queue_after_v127_v128.log"
STATE_FILE = RESULTS / "queue_after_v127_v128.state.json"
PREV_STATE = RESULTS / "queue_after_v124_v126.state.json"
BUILD_SCRIPT = REPO / "scripts" / "build_v128_repro.sh"
CANARY = ["1R55", "1G9V", "1OF6", "1T46", "1XOZ", "1Y6R"]
SUCCESS_GATE = 71


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
        state = json.load(open(STATE_FILE))
    state.update(kwargs)
    state["updated_at"] = datetime.now(timezone.utc).isoformat()
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)
        f.write("\n")


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def discover_v127_pid() -> int:
    env_pid = int(os.environ.get("QUEUE_V127_PID", "0") or "0")
    if env_pid > 0 and pid_alive(env_pid):
        return env_pid

    if PREV_STATE.exists():
        prev = json.load(open(PREV_STATE))
        v127_pid = int(prev.get("v127_pid") or 0)
        if v127_pid > 0 and pid_alive(v127_pid):
            return v127_pid

    v127_dirs = sorted(RESULTS.glob("v127_*_optB_full85"), key=lambda p: p.stat().st_mtime)
    for d in reversed(v127_dirs):
        prov = d / "launch_provenance.json"
        pid_file = d / "benchmark.pid"
        if prov.exists():
            pid = int(json.load(open(prov)).get("pid") or 0)
            if pid > 0 and pid_alive(pid):
                return pid
        if pid_file.exists():
            pid = int(pid_file.read_text().strip() or 0)
            if pid > 0 and pid_alive(pid):
                return pid
    return 0


def discover_v127_dir() -> Path | None:
    v127_dirs = sorted(RESULTS.glob("v127_*_optB_full85"), key=lambda p: p.stat().st_mtime)
    return v127_dirs[-1] if v127_dirs else None


def summarize_run(run_dir: Path) -> dict:
    json_pairs = REPO / "benchmarks/datasets/benchmark_astex_native_85.json"
    all_ids = [x["receptor_id"] for x in json.load(open(json_pairs))["pairs"]]
    if not run_dir.exists():
        return {"run_dir": str(run_dir), "completed": 0, "total": len(all_ids)}

    done = [d.name for d in run_dir.iterdir() if d.is_dir() and (d / "result.csv").exists()]
    succ = 0
    for d in done:
        rows = list(csv.DictReader(open(run_dir / d / "result.csv")))
        if rows and rows[0].get("success") == "1":
            succ += 1

    canary = {}
    for t in CANARY:
        rf = run_dir / t / "result.csv"
        if rf.exists():
            r = list(csv.DictReader(open(rf)))[0]
            canary[t] = {
                "rmsd": r.get("rmsd_to_crystal"),
                "success": r.get("success"),
                "pose_source": r.get("pose_source"),
            }

    return {
        "run_dir": run_dir.name,
        "completed": len(done),
        "total": len(all_ids),
        "success_result_csv": succ,
        "canary": canary,
        "missing": sorted(set(all_ids) - set(done)),
        "gate_pass": succ >= SUCCESS_GATE if len(done) == len(all_ids) else None,
    }


def prebuild_v128() -> None:
    if not BUILD_SCRIPT.exists():
        log(f"WARN: build script missing: {BUILD_SCRIPT}")
        return
    log(f"Pre-building v128 worktree via {BUILD_SCRIPT}")
    proc = subprocess.run(
        ["bash", str(BUILD_SCRIPT)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if proc.stdout.strip():
        log(proc.stdout.strip())
    if proc.returncode != 0:
        log(f"WARN: v128 prebuild failed exit={proc.returncode}")
        if proc.stderr.strip():
            log(f"  stderr: {proc.stderr.strip()}")
    else:
        log("v128 prebuild OK")


def launch_v128() -> int:
    script = REPO / "scripts" / "launch_v128_v50b_repro.py"
    log(f"Launching v128 via {script}")
    proc = subprocess.run(
        [sys.executable, str(script), "--skip-build"],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    log(proc.stdout.strip() or "(no stdout)")
    if proc.stderr.strip():
        log(f"v128 stderr: {proc.stderr.strip()}")
    if proc.returncode != 0:
        log("Retrying v128 launch with build step enabled")
        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(REPO),
            capture_output=True,
            text=True,
        )
        log(proc.stdout.strip() or "(no stdout)")
        if proc.stderr.strip():
            log(f"v128 stderr: {proc.stderr.strip()}")
        if proc.returncode != 0:
            log(f"ERROR: v128 launch failed exit={proc.returncode}")
            return -1

    v128_dirs = sorted(RESULTS.glob("v128_*_v50b_repro"), key=lambda x: x.stat().st_mtime)
    if v128_dirs:
        prov = v128_dirs[-1] / "launch_provenance.json"
        if prov.exists():
            pid = json.load(open(prov)).get("pid")
            log(f"v128 queued successfully pid={pid} output={v128_dirs[-1].name}")
            return int(pid or -1)
    return 0


def wait_for_v127_launch(poll_s: int) -> int:
    log("Waiting for v127 to be launched...")
    while True:
        pid = discover_v127_pid()
        if pid > 0:
            log(f"v127 detected pid={pid}")
            return pid
        if PREV_STATE.exists():
            prev = json.load(open(PREV_STATE))
            if prev.get("status") == "v127_launched":
                pid = int(prev.get("v127_pid") or 0)
                if pid > 0:
                    log(f"v127 launched (state file) pid={pid}")
                    return pid
        time.sleep(poll_s)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--poll", type=int, default=120, help="seconds between checks")
    parser.add_argument("--no-launch-v128", action="store_true")
    parser.add_argument("--no-prebuild", action="store_true")
    args = parser.parse_args()

    save_state(status="watching", poll_s=args.poll)
    log("Queue watcher v127→v128 started")

    v127_pid = discover_v127_pid()
    if v127_pid <= 0:
        v127_pid = wait_for_v127_launch(args.poll)

    v127_dir = discover_v127_dir()
    save_state(status="watching_v127", v127_pid=v127_pid, v127_dir=v127_dir.name if v127_dir else None)

    if not args.no_prebuild:
        prebuild_v128()

    while pid_alive(v127_pid):
        if v127_dir:
            s = summarize_run(v127_dir)
            log(f"v127 progress: {s['completed']}/{s['total']} success={s['success_result_csv']}")
        else:
            log(f"v127 pid={v127_pid} still running")
        time.sleep(args.poll)

    log("v127 finished — writing summary")
    if not v127_dir:
        v127_dir = discover_v127_dir()
    v127_sum = summarize_run(v127_dir) if v127_dir else {}
    if v127_sum:
        gate = v127_sum.get("gate_pass")
        log(
            f"v127 final: {v127_sum['completed']}/{v127_sum['total']} "
            f"success={v127_sum['success_result_csv']}"
        )
        if v127_sum.get("canary"):
            for t, c in sorted(v127_sum["canary"].items()):
                log(f"  v127 canary {t}: rmsd={c['rmsd']} ok={c['success']} src={c['pose_source']}")

    save_state(status="v127_done", v127_summary=v127_sum)

    v128_pid = -1
    if not args.no_launch_v128:
        v128_pid = launch_v128()
        save_state(status="v128_launched", v128_pid=v128_pid)
    else:
        log("Skipping v128 launch (--no-launch-v128)")

    log("Queue watcher v127→v128 exiting")
    return 0


if __name__ == "__main__":
    sys.exit(main())