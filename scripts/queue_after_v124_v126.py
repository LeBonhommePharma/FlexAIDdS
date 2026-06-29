#!/usr/bin/env python3
# queue_after_v124_v126.py — wait for v124 resume + v126, summarize, launch v127
#
# Usage:
#   python3 scripts/queue_after_v124_v126.py [--poll 120]
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
QUEUE_LOG = RESULTS / "queue_after_v124_v126.log"
STATE_FILE = RESULTS / "queue_after_v124_v126.state.json"

V124_DIR = "v124_full85_20260626_0413_consensus_guard"
V126_DIR = "v126_20260628_2347_optB_smoke"
CANARY = ["1R55", "1G9V", "1OF6", "1T46", "1XOZ", "1Y6R"]

# PIDs from the active campaigns (override via env if needed)
DEFAULT_PIDS = {
    "v126": int(os.environ.get("QUEUE_V126_PID", "13781")),
    "v124": int(os.environ.get("QUEUE_V124_PID", "48818")),
}


def log(msg: str) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(QUEUE_LOG, "a") as f:
        f.write(line + "\n")


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def summarize_run(run_dir: str) -> dict:
    p = RESULTS / run_dir
    json_pairs = REPO / "benchmarks/datasets/benchmark_astex_native_85.json"
    all_ids = [x["receptor_id"] for x in json.load(open(json_pairs))["pairs"]]
    done = [d.name for d in p.iterdir() if d.is_dir() and (d / "result.csv").exists()]
    succ = sent = None
    sc = p / "astex_diverse_results.csv"
    if sc.exists():
        rows = [r for r in csv.DictReader(open(sc)) if r.get("pdb_id", "").strip()]
        if rows:
            succ = sum(1 for r in rows if r.get("success") == "1")
            sent = sum(1 for r in rows if float(r.get("rmsd_to_crystal", -1) or -1) < 0)
    canary = {}
    for t in CANARY:
        rf = p / t / "result.csv"
        if rf.exists():
            r = list(csv.DictReader(open(rf)))[0]
            canary[t] = {
                "rmsd": r.get("rmsd_to_crystal"),
                "success": r.get("success"),
                "pose_source": r.get("pose_source"),
            }
    return {
        "run_dir": run_dir,
        "completed": len(done),
        "total": len(all_ids),
        "success_csv": succ,
        "sentinel_csv": sent,
        "canary": canary,
        "missing": sorted(set(all_ids) - set(done)),
    }


def save_state(**kwargs) -> None:
    state = {}
    if STATE_FILE.exists():
        state = json.load(open(STATE_FILE))
    state.update(kwargs)
    state["updated_at"] = datetime.now(timezone.utc).isoformat()
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)
        f.write("\n")


def launch_v127() -> int:
    script = REPO / "scripts" / "launch_v127_full85.py"
    log(f"Launching v127 via {script}")
    proc = subprocess.run(
        [sys.executable, str(script)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    log(proc.stdout.strip() or "(no stdout)")
    if proc.stderr.strip():
        log(f"v127 stderr: {proc.stderr.strip()}")
    if proc.returncode != 0:
        log(f"ERROR: v127 launch failed exit={proc.returncode}")
        return -1
    # Read pid from newest v127 provenance
    v127_dirs = sorted(RESULTS.glob("v127_*_optB_full85"), key=lambda x: x.stat().st_mtime)
    if v127_dirs:
        prov = v127_dirs[-1] / "launch_provenance.json"
        if prov.exists():
            pid = json.load(open(prov)).get("pid")
            log(f"v127 queued successfully pid={pid} output={v127_dirs[-1].name}")
            return int(pid or -1)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--poll", type=int, default=120, help="seconds between checks")
    parser.add_argument("--no-launch-v127", action="store_true")
    args = parser.parse_args()

    pids = dict(DEFAULT_PIDS)
    save_state(status="watching", pids=pids, poll_s=args.poll)

    log(f"Queue watcher started — waiting for PIDs v126={pids['v126']} v124={pids['v124']}")

    while pid_alive(pids["v126"]) or pid_alive(pids["v124"]):
        parts = []
        if pid_alive(pids["v126"]):
            s = summarize_run(V126_DIR)
            parts.append(f"v126 {s['completed']}/{s['total']}")
        else:
            parts.append("v126 done")
        if pid_alive(pids["v124"]):
            s = summarize_run(V124_DIR)
            parts.append(f"v124 {s['completed']}/{s['total']}")
        else:
            parts.append("v124 done")
        log("progress: " + ", ".join(parts))
        time.sleep(args.poll)

    log("Both campaigns finished — writing summaries")
    v124_sum = summarize_run(V124_DIR)
    v126_sum = summarize_run(V126_DIR)
    for label, s in [("v124", v124_sum), ("v126", v126_sum)]:
        rate = ""
        if s["success_csv"] is not None and s["completed"]:
            rate = f" success={s['success_csv']}/{s['completed']}"
        log(
            f"{label} final: {s['completed']}/{s['total']}{rate} "
            f"sentinel={s['sentinel_csv']} missing={len(s['missing'])}"
        )
        if s["canary"]:
            for t, c in sorted(s["canary"].items()):
                log(f"  {label} canary {t}: rmsd={c['rmsd']} ok={c['success']} src={c['pose_source']}")

    save_state(status="campaigns_done", v124_summary=v124_sum, v126_summary=v126_sum)

    v127_pid = -1
    if not args.no_launch_v127:
        v127_pid = launch_v127()
        save_state(status="v127_launched", v127_pid=v127_pid)
    else:
        log("Skipping v127 launch (--no-launch-v127)")

    log("Queue watcher exiting")
    return 0


if __name__ == "__main__":
    sys.exit(main())