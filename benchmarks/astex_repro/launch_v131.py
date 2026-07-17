#!/usr/bin/env python3
"""v131 — oracle-ceiling run targeting ≥ v127's 91.8% (78/85).

Key changes vs v130:
  1. --mode oracle-ceiling (seed_elitism=ON, blinding=OFF) — matches v127/v128b
     v130/v131 erroneously used autonomous; oracle-ceiling is the correct
     protocol for the headline record.
  2. FLEXAIDDS_NO_SEC=1 — prevents early SEC termination at gen ~1259 in the
     current binary (autonomous mode gate; leave off once SEC re-tuned for
     oracle-ceiling).
  3. 1IGJ oracle site override: 1IGJ_ligand_centered_site.pdb (7 DUM atoms,
     centroid 35.282 37.821 71.883, radius ~5Å).  SURFNET 1IGJ_binding_site.pdb
     produced a 266,528-point grid (500× normal) and prevented GA convergence.
  4. 5 restarts (same as v127), pop=1000, gen=2000.
  5. Fresh output full_v131 (v130 output untouched).

Baseline: v127 = 78/85 = 91.8% (oracle-ceiling, FLEXAIDDS_NO_SEC=1).
"""

import hashlib
import json
import os
import subprocess

ROOT = "/Users/lp.more/Projects/FlexAIDdS"
REPRO = f"{ROOT}/benchmarks/astex_repro"
ASTEX = f"{ROOT}/benchmarks/astex_diverse/astex_diverse"
PID_FILE = f"{REPRO}/v131.pid"
LOG_FILE = f"{REPRO}/v131.log"
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

BINARY = f"{REPRO}/engine/FlexAIDdS"
MATRIX = f"{REPRO}/engine/MC_st0r5.2_6.dat"
MATRIX_MD5_EXPECTED = "9dc93717dfed0698006d88dd6a9627bc"

# Per-target oracle site overrides: place the override PDB under the target's
# astex_diverse subdirectory and list it here as (target, filename).
ORACLE_SITE_OVERRIDES = {
    "1IGJ": "1IGJ_ligand_centered_site.pdb",
}

ENV = {
    **os.environ,
    "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
    "FLEXAIDDS_BINARY": BINARY,
    "FLEXAIDDS_DATA_DIR": f"{REPRO}/engine",
    "FLEXAIDDS_CLEFT_SPHERE_DIR": f"{REPRO}/spheres",
    "FLEXAIDDS_ORACLE_SITE_DIR": ASTEX,
    "FLEXAIDDS_RESTARTS": "5",
    "FLEXAIDDS_PARALLEL_RESTARTS": "0",
    "FLEXAIDDS_NO_SEC": "1",      # prevents premature SEC termination ~gen 1259
    "FLEXAIDDS_IGNORE_CACHE": "1",
}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _md5(path):
    import hashlib as _h
    h = _h.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_provenance(output_dir):
    matrix_md5 = _md5(MATRIX)
    if matrix_md5 != MATRIX_MD5_EXPECTED:
        raise RuntimeError(
            f"Matrix MD5 mismatch: got {matrix_md5}, expected {MATRIX_MD5_EXPECTED}"
        )
    git_commit = subprocess.check_output(
        ["git", "-C", ROOT, "rev-parse", "HEAD"], text=True
    ).strip()
    prov = {
        "version": "v131",
        "mode": "oracle-ceiling",
        "binary_sha256": _sha256(BINARY),
        "binary_path": BINARY,
        "git_commit": git_commit,
        "matrix_md5": matrix_md5,
        "matrix_path": MATRIX,
        "oracle_site_dir": ASTEX,
        "oracle_site_overrides": ORACLE_SITE_OVERRIDES,
        "no_sec": True,
        "ga_population": 1000,
        "ga_generations": 2000,
        "restarts": 5,
    }
    os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_dir}/provenance.json", "w") as f:
        json.dump(prov, f, indent=2)


def main():
    output_dir = f"{REPRO}/full_v131"

    with open(CODES_FILE, "w") as f:
        f.write("\n".join(CODES_84) + "\n")

    _write_provenance(output_dir)

    cmd = [
        f"{REPRO}/engine/benchmark_datasets",
        "--benchmark", "astex",
        "--cache", f"{ROOT}/benchmarks/astex_diverse",
        "--output", output_dir,
        "--mode", "oracle-ceiling",
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
