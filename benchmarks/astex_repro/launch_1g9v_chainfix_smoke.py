#!/usr/bin/env python3
"""1G9V chain-fix smoke — runs beside live v130 (does not touch full_v130).

Validates A+B+C receptor prep + AUTONOMOUS multimer prune (f847044d8).
Uses fixed benchmark_datasets binary; 1 restart only to limit memory contention.
"""

import os

ROOT = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"
PID_FILE = f"{REPRO}/smoke_1g9v_abc.pid"
LOG_FILE = f"{REPRO}/smoke_1g9v_abc.log"
CODES_FILE = f"{REPRO}/smoke_1g9v_only.txt"
OUT_DIR = f"{REPRO}/smoke_1g9v_abc"
# Fixed runner (prune + ligand-centroid); live v130 keeps old binary inode.
RUNNER = f"{REPRO}/engine/benchmark_datasets.v130_1g9v_chainfix"

ENV = {
    **os.environ,
    "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY": f"{REPRO}/engine/FlexAIDdS",
    "FLEXAIDDS_DATA_DIR": f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR": f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR": f"{ROOT}/benchmarks/astex_diverse/astex_diverse",
    "FLEXAIDDS_RESTARTS": "1",  # smoke: one restart (v130 uses 5)
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    "FLEXAIDDS_IGNORE_CACHE": "1",
}


def main():
    with open(CODES_FILE, "w") as f:
        f.write("1G9V\n")

    cmd = [
        RUNNER,
        "--benchmark", "astex",
        "--cache", f"{ROOT}/benchmarks/astex_diverse",
        "--output", OUT_DIR,
        "--mode", "autonomous",
        "--threads", "1",
        "--omp-threads", "1",
        "--ga-population", "1000",
        "--ga-generations", "2000",
        "--job-timeout-seconds", "7200",
        "--only-codes", CODES_FILE,
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

    os.execve(cmd[0], cmd, ENV)


if __name__ == "__main__":
    main()
