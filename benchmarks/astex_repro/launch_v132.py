#!/usr/bin/env python3
"""v132 — clean restart with parallel chunk fixes (CF-primary election).

Killed polluted/partial v131. Stack:

  Chunk 1: offline re-election script (scripts/offline_reelect_cf_primary.py)
  Chunk 2: CF-primary pool election + pathological CF floor + consensus guard
           (default ON; FLEXAIDDS_CF_PRIMARY=0 restores Z+H)
  Chunk 3: AUTONOMOUS multimer prune, 1G9V A+B+C, site preflight, chain catalog
  Chunk 5: FLEXAIDDS_POP_CAP=2000 (default in code)

  + CF-stagnation, pop-scale, sharing_alpha inverse, IGNORE_CACHE, coarse_init=false
  + Non-PGO binaries; fresh output full_v132
"""

import os

ROOT = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"
PID_FILE = f"{REPRO}/v132.pid"
LOG_FILE = f"{REPRO}/v132.log"
CODES_FILE = f"{REPRO}/astex84_no1hq2.txt"

CODES_84 = [
    "1G9V", "1GM8", "1GPK", "1HNN", "1HP0", "1IA1", "1IGJ", "1J3J", "1JD0", "1JJE",
    "1K3U", "1KE5", "1KZK", "1L2S", "1L7F", "1LPZ", "1M2Z", "1MEH", "1MQ6", "1N1M",
    "1N2J", "1N2V", "1N46", "1NAV", "1OF1", "1OF6", "1OPK", "1OQ5", "1OWE", "1P2Y",
    "1P62", "1PMN", "1Q1G", "1Q41", "1Q4G", "1R1H", "1R55", "1R58", "1R9O", "1S19",
    "1S3V", "1SG0", "1SJ0", "1SQ5", "1T40", "1T46", "1T9B", "1TT1", "1TW6", "1TZ8",
    "1U1C", "1U4D", "1UML", "1UNL", "1UOU", "1V0P", "1V48", "1V4S", "1VCJ", "1W1P",
    "1W2G", "1X8X", "1XM6", "1XOZ", "1Y6B", "1Y6R", "1YGC", "1YQY", "1YV3", "1YVF",
    "1YWR", "1Z95", "2BM2", "2BR1", "2BSM", "2BYS", "2C3I", "2CET", "2CGR", "2D3U",
    "2GBP", "2HB1", "2HR7", "2J62",
]

ENV = {
    **os.environ,
    "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY": f"{REPRO}/engine/FlexAIDdS",
    "FLEXAIDDS_DATA_DIR": f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR": f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR": f"{ROOT}/benchmarks/astex_diverse/astex_diverse",
    "FLEXAIDDS_RESTARTS": "5",
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    "FLEXAIDDS_IGNORE_CACHE": "1",
    "FLEXAIDDS_CF_PRIMARY": "1",
    "FLEXAIDDS_CF_FLOOR": "-250",
    "FLEXAIDDS_POP_CAP": "2000",
}


def main():
    with open(CODES_FILE, "w") as f:
        f.write("\n".join(CODES_84) + "\n")

    cmd = [
        f"{REPRO}/engine/benchmark_datasets",
        "--benchmark", "astex",
        "--cache", f"{ROOT}/benchmarks/astex_diverse",
        "--output", f"{REPRO}/full_v132",
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

    os.execve(cmd[0], cmd, ENV)


if __name__ == "__main__":
    main()
