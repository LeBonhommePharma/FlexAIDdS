#!/usr/bin/env python3
"""v130 — coarse_init disabled + non-PGO binary + binary contamination scrub.

Key fixes vs v129 (which produced 0/84 due to OOM + PGO overhead):

1. coarse_init DISABLED: commit 5498fee7c added coarse_init AFTER v127's 91.8%
   record. At OMP=1, 25 seeds × 16 orientations × ~70ms/pose = ~28s, but that
   is with 8 OMP threads. Single-threaded it runs ~224s per restart before GA
   even starts. 5 restarts × 85 targets × 224s = 27 hours of pure overhead.
   DatasetRunner.cpp now hardcodes enabled=false (v130 comment).

2. Non-PGO binary: engine/FlexAIDdS was built with FLEXAIDS_PGO_GENERATE=ON
   (Stage-1 instrumentation, 8-10× overhead per basic block). Rebuilt with
   PGO_GENERATE=OFF and PGO_USE= (clean optimized binary). SHA d30652373fa6.

3. Binary contamination scrub: benchmark_datasets (SHA c2d3679387ab) includes
   commit 60bd43c32 (fix/per-restart-binary-contamination): pose .pdb files
   from prior benchmark runs are scrubbed before IGNORE_CACHE re-runs, so
   stale poses from a different binary cannot contaminate RMSD results.

4. IGNORE_CACHE=1 retained: bypasses any stale result.csv from v129 (0/84
   sentinel run). All 84 targets re-run from scratch.

Timing (estimated):
  - Setup per restart (OMP=1): ~61s (grid computation, no coarse_init)
  - GA per-eval (OMP=1): ~130 µs → 1000 pop × 2000 gen = 260s
  - Regime A (SEC early exit at gen ~300): 61 + 39 = 100s per restart
  - 5 restarts Regime A: ~500s ≈ 8.3 min per target
  - 85 targets mostly Regime A: ~12h overnight run

Baseline: v127 = 78/85 = 91.8% (oracle-ceiling). v129 = 0/84 (OOM + PGO).
Expected: ≥ v127 (v130 should be identical in search quality to v127 modulo
the binary contamination scrub and SEC improvements from 62a216549).

Binaries: FlexAIDdS d30652373fa6, benchmark_datasets c2d3679387ab.
"""

import os

ROOT = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"
PID_FILE = f"{REPRO}/v130.pid"
LOG_FILE = f"{REPRO}/v130.log"
CODES_FILE = f"{REPRO}/astex84_no1hq2.txt"

# 84 codes = all 85 Astex Diverse minus 1HQ2 (8 flex bonds; 5 restarts at
# full 2000 gens exceed 3h timeout even with SEC active for rigid-ish cases).
CODES_84 = [
    "1G9V","1GM8","1GPK","1HNN","1HP0","1IA1","1IGJ","1J3J","1JD0","1JJE",
    "1K3U","1KE5","1KZK","1L2S","1L7F","1LPZ","1M2Z","1MEH","1MQ6","1N1M",
    "1N2J","1N2V","1N46","1NAV","1OF1","1OF6","1OPK","1OQ5","1OWE","1P2Y",
    "1P62","1PMN","1Q1G","1Q41","1Q4G","1R1H","1R55","1R58","1R9O","1S19",
    "1S3V","1SG0","1SJ0","1SQ5","1T40","1T46","1T9B","1TT1","1TW6","1TZ8",
    "1U1C","1U4D","1UML","1UNL","1UOU","1V0P","1V48","1V4S","1VCJ","1W1P",
    "1W2G","1X8X","1XM6","1XOZ","1Y6B","1Y6R","1YGC","1YQY","1YV3","1YVF",
    "1YWR","1Z95","2BM2","2BR1","2BSM","2BYS","2C3I","2CET","2CGR","2D3U",
    "2GBP","2HB1","2HR7","2J62",
]

ENV = {
    **os.environ,
    "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY":           f"{REPRO}/engine/FlexAIDdS",
    "FLEXAIDDS_DATA_DIR":         f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR": f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR":  f"{ROOT}/benchmarks/astex_diverse/astex_diverse",
    "FLEXAIDDS_RESTARTS":         "5",
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    # NO_SEC intentionally absent: SEC joint gate (CF stagnant + allele_H<0.300)
    # is correct after 62a216549 CF-stagnation fix.
    "FLEXAIDDS_IGNORE_CACHE":     "1",  # force re-run; ignore stale v129 results
}


def main():
    # Write codes to file (avoids ENAMETOOLONG in exec argv)
    with open(CODES_FILE, "w") as f:
        f.write("\n".join(CODES_84) + "\n")

    CMD = [
        f"{REPRO}/engine/benchmark_datasets",
        "--benchmark", "astex",
        "--cache", f"{ROOT}/benchmarks/astex_diverse",
        "--output", f"{REPRO}/full_v130",   # fresh directory, no stale restarts
        "--mode", "autonomous",
        "--threads", "1",
        "--omp-threads", "1",
        "--ga-population", "1000",
        "--ga-generations", "2000",
        "--job-timeout-seconds", "10800",
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

    os.execve(CMD[0], CMD, ENV)


if __name__ == "__main__":
    main()
