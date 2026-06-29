#!/usr/bin/env python3
# launch_v131_full85.py — v109-chase run: r0=7, consensus OFF, full fixes bundled
#
# Prescription (Claude/Grok consensus):
#   r0=7 (v109) + HEAD bug fixes (H-bond/VCT, sulfo remap, wal cap) +
#   consensus OFF (avoid Fix-B selector stacking on flat landscapes) +
#   pipeline holo (1G9V apo A+B+C, 1TW6_holo A+B) + 1HNN expB site.
#
# Success gate: beat v109 record 80/85 (94.1%). Forecast: 82-83/85 if clean.
#
# Note: logsumexp boltzmann_composite remains in HEAD binary (no runtime off-switch).
# CONSENSUS_SCORER=0 restores v124-default consensus path; primary v109 delta is r0=7.
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = "/Users/lp.more/Projects/FlexAIDdS"
BUILD = f"{REPO}/build_lto"
BINARY_SRC = f"{BUILD}/FlexAIDdS"
BINARY = "/tmp/FlexAIDdS_v131"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85_v131.json"
HNN_LIGAND_SITE = f"{ORACLE_DIR}/1HNN/1HNN_ligand_centered_site.pdb"
TW6_HOLO = f"{ORACLE_DIR}/1TW6/1TW6_holo.pdb"
RESULTS_DIR = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"

TAG = datetime.datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT = f"{RESULTS_DIR}/v131_{TAG}_r07_nofixb_full85"
CACHE = f"{RESULTS_DIR}/cache_v131_r07_nofixb"
PROV = f"{OUTPUT}/launch_provenance.json"

REF_V109_DIR = "v109_20260626_tier1_consensus5r"
REF_V130_DIR = "v130_20260629_0548_sulfo_expB_full85"
HOLO_SMOKE_DIR = "holo_smoke_20260629_0658_1G9V_1TW6"
MIN_COMMIT = "04ff1735"
SUCCESS_GATE = 80

BENCH_THREADS = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")

ENV_SNAPSHOT_KEYS = (
    "FLEXAIDDS_BINARY",
    "FLEXAIDDS_ORACLE_SITE_DIR",
    "FLEXAIDDS_RESTARTS",
    "FLEXAIDDS_PARALLEL_RESTARTS",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL",
    "FLEXAIDDS_CONSENSUS_SCORER",
    "FLEXAIDDS_SEED_ELITISM",
    "FLEXAIDDS_N_ELITE",
    "FLEXAIDDS_BUDGET_SCALE",
    "FLEXAIDDS_SOFTCORE_WAL",
    "FLEXAIDDS_SOFTCORE_FLOOR",
    "FLEXAIDDS_T_HOT",
    "FLEXAIDDS_NATIVE_SEED_FRAC",
    "FLEXAIDDS_VCT_R0",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
    "FLEXAIDDS_DATA_DIR",
    "FLEXAIDDS_BENCH_CACHE",
    "FLEXAIDDS_ALLOW_CONCURRENT",
)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    for p in (
        BINARY_SRC,
        RUNNER,
        ORACLE_DIR,
        JSON_PAIRS,
        HNN_LIGAND_SITE,
        TW6_HOLO,
        f"{DATA_DIR}/MC_st0r5.2_6.dat",
    ):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    git_commit = subprocess.check_output(
        ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
    ).strip().split()[0]

    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", MIN_COMMIT, "HEAD"], cwd=REPO
    ).returncode != 0:
        sys.exit(f"ERROR: HEAD missing sulfonamide fix {MIN_COMMIT}")

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

    native = json.load(open(JSON_PAIRS))
    assert len(native["pairs"]) == 85
    hnn = next(p for p in native["pairs"] if p["receptor_id"] == "1HNN")
    tw6 = next(p for p in native["pairs"] if p["receptor_id"] == "1TW6")
    assert hnn["oracle_site_pdb"] == HNN_LIGAND_SITE
    assert tw6["receptor_pdb"] == TW6_HOLO

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY":                BINARY,
        "FLEXAIDDS_BUILD":                 BUILD,
        "FLEXAIDDS_REPO":                  REPO,
        "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
        "FLEXAIDDS_RESTARTS":              "5",
        "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
        "FLEXAIDDS_CONSENSUS_SCORER":      "0",
        "FLEXAIDDS_SEED_ELITISM":          "1",
        "FLEXAIDDS_N_ELITE":               "1",
        "FLEXAIDDS_BUDGET_SCALE":          "1",
        "FLEXAIDDS_SOFTCORE_WAL":          "1",
        "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
        "FLEXAIDDS_T_HOT":                 "500",
        "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
        "FLEXAIDDS_VCT_R0":                "7",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
        "FLEXAIDDS_DATA_DIR":              DATA_DIR,
        "FLEXAIDDS_ALLOW_CONCURRENT":      "1",
        "FLEXAIDDS_BENCH_CACHE":           CACHE,
        "OMP_WAIT_POLICY":                 "passive",
        "OMP_PLACES":                      "cores",
        "OMP_PROC_BIND":                   "spread",
    })
    for k in (
        "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
        "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
        "FLEXAIDDS_VCT_NORM",
        "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_RING_FLEX",
        "FLEXAIDDS_THERMO", "FLEXAIDDS_HVIB",
        "FLEXAIDDS_PRIORITY_TARGETS", "FLEXAIDDS_FREQSEL",
    ):
        env.pop(k, None)

    cmd = [
        "caffeinate", "-i",
        RUNNER,
        "--benchmark",           f"crossdock_json:{JSON_PAIRS}",
        "--output",              OUTPUT,
        "--threads",             BENCH_THREADS,
        "--temperature",         "298",
        "--job-timeout-seconds", "7200",
        "--cache",               CACHE,
        "--mode",                "oracle-ceiling",
    ]

    os.makedirs(OUTPUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)

    print("\nLaunching v131 full-85 — r0=7, consensus OFF, holo pipeline + sulfo + expB")
    print(f"  commit    : {git_commit}")
    print(f"  output    : {OUTPUT}")
    print(f"  gate      : >={SUCCESS_GATE}/85 (beat v109)")
    print(f"  threads   : {BENCH_THREADS}")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":       "v131_r07_nofixb_full85",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    git_commit,
        "description": (
            "v109-chase full-85: r0=7 + consensus OFF + HEAD bug fixes + sulfo remap + "
            "1HNN expB site + 1TW6_holo receptor. Target >=80/85."
        ),
        "success_gate":  f"{SUCCESS_GATE}/85",
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs":     JSON_PAIRS,
        "output_dir":     OUTPUT,
        "cache_dir":      CACHE,
        "pid":            child_pid,
        "reference_runs": {
            "v109_record": REF_V109_DIR,
            "v130_prior": REF_V130_DIR,
            "holo_smoke": HOLO_SMOKE_DIR,
        },
        "protocol_delta_vs_v130": {
            "FLEXAIDDS_VCT_R0": "7 (was 4)",
            "FLEXAIDDS_CONSENSUS_SCORER": "0 (was 1)",
            "1TW6_receptor": "1TW6_holo.pdb (was 1TW6.pdb)",
        },
        "forecast": {
            "conservative": "79/85",
            "target": "82-83/85",
            "record": "80/85 v109",
        },
        "audit_notes": {
            "logsumexp_still_active": (
                "HEAD binary includes Option B logsumexp boltzmann_composite; "
                "no runtime disable. CONSENSUS_SCORER=0 is the actionable no-Fix-B lever."
            ),
            "holo_smoke_validated": "1G9V 0.9Å + 1TW6 1.4Å in holo_smoke_0658",
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(PROV, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\nv131 launched pid={child_pid}")
    print(f"  prov: {PROV}")
    return child_pid


if __name__ == "__main__":
    main()