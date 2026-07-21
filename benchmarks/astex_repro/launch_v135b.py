#!/usr/bin/env python3
"""v135b — oracle-ceiling run with GA param fix via FlexAIDdS wrapper.

v134 failed (5/85 = 5.9%) because benchmark_datasets oracle-ceiling mode
hardcodes broken GA params causing population collapse:
  sharing_alpha=2.28571 (should be 4)
  num_chromosomes=1750   (should be 1000)
  num_generations=2000   (using 875)
  vct_dist_weight_r0=4   (should be 7) ← partially fixed via pre-gen; wrapper also covers it

Fix: FLEXAIDDS_BINARY points to flexaidds_v135_wrapper.sh, which patches
the rN/dock_config.json AFTER the runner generates it but BEFORE FlexAIDdS
reads it, then execs the real build_v135/FlexAIDdS binary.
"""

import json
import os
import subprocess
import sys
from datetime import datetime

ROOT    = "/Users/lp.more/Projects/FlexAIDdS"
REPRO   = f"{ROOT}/benchmarks/astex_repro"
RESULTS = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"

TIMESTAMP  = datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT_DIR = f"{RESULTS}/v135_full85_clean_{TIMESTAMP}"
PID_FILE   = f"{RESULTS}/v135_full85_clean_{TIMESTAMP}_launcher.pid"
LOG_FILE   = f"{RESULTS}/v135_full85_clean_{TIMESTAMP}_launcher.log"

REAL_BINARY = f"{ROOT}/build_v135/FlexAIDdS"
WRAPPER     = f"{REPRO}/flexaidds_v135_wrapper.sh"
RUNNER      = f"{REPRO}/engine/benchmark_datasets"
ORACLE_DIR  = f"{ROOT}/benchmarks/astex_diverse/astex_diverse"


def main():
    # ── Pre-flight checks ────────────────────────────────────────────────────
    for path, label in [(REAL_BINARY, "binary"), (WRAPPER, "wrapper"), (RUNNER, "runner")]:
        if not os.path.isfile(path):
            sys.exit(f"FATAL: {label} not found: {path}")

    # ── Record git HEAD and binary SHA ───────────────────────────────────────
    git_head   = subprocess.check_output(["git", "-C", ROOT, "rev-parse", "HEAD"], text=True).strip()
    binary_sha = subprocess.check_output(["shasum", "-a", "256", REAL_BINARY], text=True).split()[0]

    # ── Create output dir and write provenance ───────────────────────────────
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    provenance = {
        "dataset":           "astex_diverse",
        "version":           "v135b",
        "binary_path":       REAL_BINARY,
        "binary_sha256":     binary_sha,
        "git_commit":        git_head,
        "oracle_site_dir":   ORACLE_DIR,
        "oracle_site_dir_set": True,
        "wrapper":           WRAPPER,
        "fix_notes": (
            "flexaidds_v135_wrapper.sh patches rN/dock_config.json before exec: "
            "sharing_alpha=4 (was 2.28571), num_chromosomes=1000 (was 1750), "
            "num_generations=875, vct_dist_weight_r0=7 (was 4)"
        ),
    }
    with open(os.path.join(OUTPUT_DIR, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=2)

    # ── Report to caller before fork ─────────────────────────────────────────
    print(f"[v135b] git_commit : {git_head}", flush=True)
    print(f"[v135b] binary_sha : {binary_sha}", flush=True)
    print(f"[v135b] output_dir : {OUTPUT_DIR}", flush=True)
    print(f"[v135b] pid_file   : {PID_FILE}", flush=True)
    print(f"[v135b] wrapper    : {WRAPPER}", flush=True)
    sys.stdout.flush()

    # ── Build command (identical to v134 invocation) ─────────────────────────
    cmd = [
        RUNNER,
        "--benchmark",   "astex",
        "--mode",        "oracle-ceiling",
        "--threads",     "5",
        "--omp-threads", "2",
        "--output",      OUTPUT_DIR,
    ]

    env = {
        **os.environ,
        "PATH":                    "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
        "FLEXAIDDS_BINARY":        WRAPPER,        # ← wrapper, not real binary
        "FLEXAIDDS_DATA_DIR":      f"{REPRO}/engine",
        "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,
        # FLEXAIDDS_IGNORE_CACHE absent → runner uses pre-built rN configs when present
    }

    # ── Double-fork / os.setsid() daemon launch ───────────────────────────────
    if os.fork() > 0:
        os._exit(0)

    os.setsid()

    if os.fork() > 0:
        os._exit(0)

    os.chdir(REPRO)
    log_fh  = open(LOG_FILE, "ab", buffering=0)
    os.dup2(log_fh.fileno(), 1)
    os.dup2(log_fh.fileno(), 2)
    devnull = open(os.devnull, "rb")
    os.dup2(devnull.fileno(), 0)

    with open(PID_FILE, "w") as fh:
        fh.write(str(os.getpid()) + "\n")

    os.execve(cmd[0], cmd, env)


if __name__ == "__main__":
    main()
