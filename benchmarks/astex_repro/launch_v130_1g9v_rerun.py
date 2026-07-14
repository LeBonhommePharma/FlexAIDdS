#!/usr/bin/env python3
"""Requeue 1G9V with HEM-stripped receptor (commit 225060b55).

v130 processed 1G9V first (before the HEM fix) → RMSD=12.0 due to
iron-porphyrin HEM cofactor in apo PDB causing CF=+343 at crystal coords.
Fix: 1G9V_apo.pdb now contains chains A+B+C only (3261 ATOM, 0 HETATM).

Runs immediately (parallel to v130 which is past 1G9V in sequential order).
IGNORE_CACHE=1 forces overwrite of the poisoned result.csv.
Output lands in full_v130/1G9V/ — same dir used by the main v130 run.
"""

import os

ROOT  = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"

LOG_FILE   = f"{REPRO}/v130_1g9v_rerun.log"
PID_FILE   = f"{REPRO}/v130_1g9v_rerun.pid"
CODE_FILE  = f"{REPRO}/1g9v_only.txt"

ENV = {
    **os.environ,
    "PATH":                         "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY":             f"{REPRO}/engine/FlexAIDdS",
    "FLEXAIDDS_DATA_DIR":           f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR":   f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR":    f"{ROOT}/benchmarks/astex_diverse/astex_diverse",
    "FLEXAIDDS_RESTARTS":           "5",
    "FLEXAIDDS_PARALLEL_RESTARTS":  "0",
    "FLEXAIDDS_IGNORE_CACHE":       "1",  # overwrite poisoned 1G9V result
}


def main():
    with open(CODE_FILE, "w") as f:
        f.write("1G9V\n")

    CMD = [
        f"{REPRO}/engine/benchmark_datasets",
        "--benchmark", "astex",
        "--cache",     f"{ROOT}/benchmarks/astex_diverse",
        "--output",    f"{REPRO}/full_v130",    # same output dir as main v130 run
        "--mode",      "autonomous",
        "--threads",   "1",
        "--omp-threads", "1",
        "--ga-population", "1000",
        "--ga-generations", "2000",
        "--job-timeout-seconds", "10800",
        "--only-codes", CODE_FILE,
    ]

    if os.fork() > 0:
        os._exit(0)
    os.setsid()
    if os.fork() > 0:
        os._exit(0)

    os.chdir(REPRO)
    log = open(LOG_FILE, "ab", buffering=0)
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    devnull = open(os.devnull, "rb")
    os.dup2(devnull.fileno(), 0)

    with open(PID_FILE, "w") as f:
        f.write(str(os.getpid()) + "\n")

    os.execve(CMD[0], CMD, ENV)


if __name__ == "__main__":
    main()
