#!/usr/bin/env python3
# monitor_campaigns.py — passive dashboard for v124 / v126 / queue / v127
#
# Usage:
#   python3 scripts/monitor_campaigns.py [--poll 60]
#
# Logs to: ~/Documents/PhD/Programs/FlexAIDdS/results/campaign_monitor.log
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import argparse
import csv
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path("/Users/lp.more/Projects/FlexAIDdS")
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
LOG = RESULTS / "campaign_monitor.log"
STATE = RESULTS / "campaign_monitor.state.json"

V124_DIR = "v124_full85_20260626_0413_consensus_guard"
V126_DIR = "v126_20260628_2347_optB_smoke"
CANARY = ["1R55", "1G9V", "1OF6", "1T46", "1XOZ", "1Y6R"]

WATCH_PIDS = {
    "v126": int(os.environ.get("MON_V126_PID", "13781")),
    "v124": int(os.environ.get("MON_V124_PID", "48818")),
    "queue": int(os.environ.get("MON_QUEUE_PID", "51275")),
}


def log(msg: str) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def count_run(run_dir: str) -> dict:
    p = RESULTS / run_dir
    if not p.is_dir():
        return {"exists": False}
    pairs = json.load(open(REPO / "benchmarks/datasets/benchmark_astex_native_85.json"))
    all_ids = [x["receptor_id"] for x in pairs["pairs"]]
    done = [d.name for d in p.iterdir() if d.is_dir() and (d / "result.csv").exists()]
    inprog = [
        d.name for d in p.iterdir()
        if d.is_dir() and d.name in all_ids
        and not (d / "result.csv").exists()
        and any(d.rglob("dock_config.json"))
    ]
    succ = sent = None
    sc = p / "astex_diverse_results.csv"
    if sc.exists():
        rows = [r for r in csv.DictReader(open(sc)) if r.get("pdb_id", "").strip()]
        if rows:
            succ = sum(1 for r in rows if r.get("success") == "1")
            sent = sum(1 for r in rows if float(r.get("rmsd_to_crystal", -1) or -1) < 0)
    canary_ok = 0
    for t in CANARY:
        rf = p / t / "result.csv"
        if rf.exists() and list(csv.DictReader(open(rf)))[0].get("success") == "1":
            canary_ok += 1
    return {
        "exists": True,
        "done": len(done),
        "total": len(all_ids),
        "inprog": inprog[:5],
        "success": succ,
        "sentinel": sent,
        "canary_ok": canary_ok,
    }


def latest_v127_dir() -> Path | None:
    dirs = sorted(RESULTS.glob("v127_*_optB_full85"), key=lambda x: x.stat().st_mtime)
    return dirs[-1] if dirs else None


def snapshot() -> dict:
    v127 = latest_v127_dir()
    qstate = {}
    if (RESULTS / "queue_after_v124_v126.state.json").exists():
        qstate = json.load(open(RESULTS / "queue_after_v124_v126.state.json"))
    return {
        "pids": {k: {"pid": v, "alive": alive(v)} for k, v in WATCH_PIDS.items()},
        "v124": count_run(V124_DIR),
        "v126": count_run(V126_DIR),
        "v127": count_run(v127.name) if v127 else {"exists": False},
        "v127_dir": v127.name if v127 else None,
        "queue_status": qstate.get("status"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--poll", type=int, default=60)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()

    log("Campaign monitor started")
    while True:
        s = snapshot()
        v124 = s["v124"]
        v126 = s["v126"]
        parts = []
        if v124.get("exists"):
            parts.append(
                f"v124 {v124['done']}/{v124['total']} ok={v124['success']} "
                f"canary={v124['canary_ok']}/6 active={','.join(v124['inprog']) or '-'}"
            )
        if v126.get("exists"):
            parts.append(
                f"v126 {v126['done']}/{v126['total']} ok={v126['success']} "
                f"canary={v126['canary_ok']}/6 active={','.join(v126['inprog']) or '-'}"
            )
        pid_bits = []
        for name, info in s["pids"].items():
            pid_bits.append(f"{name}={info['pid']}{'↑' if info['alive'] else '✗'}")
        parts.append("pids " + " ".join(pid_bits))
        if s.get("v127_dir"):
            v7 = s["v127"]
            parts.append(f"v127 {v7['done']}/{v7['total']} ({s['v127_dir']})")
        elif s.get("queue_status"):
            parts.append(f"queue={s['queue_status']}")
        line = " | ".join(parts)
        log(line)
        s["logged_at"] = datetime.now(timezone.utc).isoformat()
        with open(STATE, "w") as f:
            json.dump(s, f, indent=2)
            f.write("\n")
        if args.once:
            break
        time.sleep(args.poll)
    return 0


if __name__ == "__main__":
    sys.exit(main())