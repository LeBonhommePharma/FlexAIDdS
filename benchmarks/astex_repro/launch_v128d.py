#!/usr/bin/env python3
"""v128d launcher — same config as v128c but FLEXAIDDS_RESTARTS=5 (matches v127
restart depth for direct comparison). Double-fork/setsid isolation survives
terminal teardown. Resumes cached targets in full_v128/ (17 cached + clean 1HQ2)."""
import os

ROOT = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"
PID_FILE = f"{REPRO}/v128d.pid"
LOG_FILE = f"{REPRO}/v128d.log"

ENV = {
    **os.environ,
    "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY": f"{REPRO}/engine/FlexAIDdS",
    "FLEXAIDDS_DATA_DIR": f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR": f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR": f"{ROOT}/benchmarks/astex_diverse/astex_diverse",
    "FLEXAIDDS_RESTARTS": "5",
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    "FLEXAIDDS_NO_SEC": "1",
}

CMD = [
    f"{REPRO}/engine/benchmark_datasets",
    "--benchmark", "astex",
    "--cache", f"{ROOT}/benchmarks/astex_diverse",
    "--output", f"{REPRO}/full_v128",
    "--mode", "autonomous",
    "--threads", "1",
    "--omp-threads", "1",
    "--ga-population", "1000",
    "--ga-generations", "2000",
    "--job-timeout-seconds", "10800",
]


def main():
    # First fork: parent returns to shell.
    if os.fork() > 0:
        os._exit(0)
    # New session leader (detach from controlling terminal).
    os.setsid()
    # Second fork: grandchild can never re-acquire a controlling terminal.
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
