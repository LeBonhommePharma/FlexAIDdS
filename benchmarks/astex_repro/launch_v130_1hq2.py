#!/usr/bin/env python3
"""v130-1HQ2 — standalone 1HQ2 job queued after v130 (PID 98837) finishes.

1HQ2 was excluded from v130's 84-code run with the comment "8 flex bonds;
~3600s per restart". Opus audit showed the old timeout (wall_time=10810s,
0 poses) was caused by three now-eliminated factors:
  1. coarse_init enabled   (+224s/restart at OMP=1)
  2. n_gen=4000            (double v130's 2000)
  3. PGO-instrumented binary (8-10x overhead per basic block)

With the fixed v130 config (non-PGO d30652373fa6, coarse_init=false, gen=2000):
  Regime A (SEC ~gen300): 5 × 139s = 695s  ≈ 12 min
  Regime B (full 2000 gen): 5 × 581s = 2905s ≈ 48 min
Both are comfortably inside the 10800s (3h) timeout.

Pop budget for 1HQ2 (n_flex=8):
  pop_effective = 1000 × (8/4.0) = 2000
  n_gen         = 2000  (pop-scaling only, NOT gen-scaling)
  total_evals   = 4,000,000
  sharing_alpha = 4.0 × 1000/2000 = 2.0
  n_genes       = 12   (8 dihedral + 4 rigid-body; below budget_scale ≥14)

Caveats:
  - May still fail on docking quality (r0=7 false-minimum landscape).
  - Result added to full_v130/ alongside the 84-target run for unified scoring.
"""

import os
import time
import sys

ROOT  = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"

V130_PID  = 98837          # PID of the running 84-target v130 job
POLL_S    = 60             # check every 60 s
PID_FILE  = f"{REPRO}/v130_1hq2.pid"
LOG_FILE  = f"{REPRO}/v130_1hq2.log"

ENV = {
    **os.environ,
    "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY":           f"{REPRO}/engine/FlexAIDdS",
    "FLEXAIDDS_DATA_DIR":         f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR": f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR":  f"{ROOT}/benchmarks/astex_diverse/astex_diverse",
    "FLEXAIDDS_RESTARTS":         "5",
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    "FLEXAIDDS_IGNORE_CACHE":     "1",
    # NO_SEC absent: SEC joint gate (CF stagnant + allele_H < 0.300) is correct
}

CMD = [
    f"{REPRO}/engine/benchmark_datasets",
    "--benchmark",         "astex",
    "--cache",             f"{ROOT}/benchmarks/astex_diverse",
    "--output",            f"{REPRO}/full_v130",   # same dir as 84-target run
    "--mode",              "autonomous",
    "--threads",           "1",
    "--omp-threads",       "1",
    "--ga-population",     "1000",
    "--ga-generations",    "2000",
    "--job-timeout-seconds", "10800",
    "--only-codes",        "1HQ2",
]


def pid_alive(pid: int) -> bool:
    """Return True if process with given PID exists."""
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def wait_for_v130() -> None:
    """Block until v130 (PID 98837) exits, logging progress every poll."""
    print(f"[queue] Waiting for v130 PID {V130_PID} to finish …", flush=True)
    while pid_alive(V130_PID):
        time.sleep(POLL_S)
    print(f"[queue] v130 PID {V130_PID} is gone — launching 1HQ2.", flush=True)


def main() -> None:
    # ── Phase 1: wait (runs in the intermediate child before 2nd fork) ──────
    # Double-fork so the queue daemon is fully detached from this shell.
    if os.fork() > 0:
        os._exit(0)          # parent exits immediately
    os.setsid()
    if os.fork() > 0:
        os._exit(0)          # session leader exits

    # ── Phase 2: redirect I/O ───────────────────────────────────────────────
    os.chdir(REPRO)
    log = open(LOG_FILE, "ab", buffering=0)
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    devnull = open(os.devnull, "rb")
    os.dup2(devnull.fileno(), 0)

    with open(PID_FILE, "w") as f:
        f.write(str(os.getpid()) + "\n")

    # ── Phase 3: wait for v130 then exec ────────────────────────────────────
    wait_for_v130()
    os.execve(CMD[0], CMD, ENV)


if __name__ == "__main__":
    main()
