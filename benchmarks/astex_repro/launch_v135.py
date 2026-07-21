#!/usr/bin/env python3
"""v135 — clean oracle-ceiling run fixing v134's GA collapse.

Root cause of v134 failure (5/85 = 5.9%):
  engine/benchmark_datasets (Jul 14 binary, SHA c35413a2) oracle-ceiling mode
  hardcodes bad GA params that cause population diversity collapse:
    sharing_alpha=2.28571 (should be 4)   → niche penalty too flat, all pop converges
    vct_dist_weight_r0=4   (should be 7)   → wrong CF landscape
    num_chromosomes=1750   (should be 1000) → oversized pop with wrong alpha

v135 fix:
  Pre-generate all 85 oracle-ceiling dock_configs with CORRECT params BEFORE
  launching benchmark_datasets. Without FLEXAIDDS_IGNORE_CACHE, the runner
  uses pre-built configs rather than regenerating from its broken defaults.

GA params restored to v133-era values:
  sharing_alpha=4, vct_dist_weight_r0=7, num_chromosomes=1000, num_generations=875

C++ binary: build_v135/FlexAIDdS (git HEAD 2a60f65132, SHA 36cbede1...)
  (Dead I_ES code was LTO-eliminated; binary identical to v134's C++ output
   but the GA param fix is what matters for correctness.)
"""

import json
import os
import subprocess
import sys
from datetime import datetime

ROOT       = "/Users/lp.more/Projects/FlexAIDdS"
REPRO      = f"{ROOT}/benchmarks/astex_repro"
RESULTS    = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"
TIMESTAMP  = datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT_DIR = f"{RESULTS}/v135_full85_clean_{TIMESTAMP}"
PID_FILE   = f"{RESULTS}/v135_full85_clean_{TIMESTAMP}_launcher.pid"
LOG_FILE   = f"{RESULTS}/v135_full85_clean_{TIMESTAMP}_launcher.log"

BINARY     = f"{ROOT}/build_v135/FlexAIDdS"
RUNNER     = f"{REPRO}/engine/benchmark_datasets"
LIGAND_CACHE = "/Users/lp.more/.flexaidds/benchmarks/astex_diverse"
ORACLE_DIR   = f"{ROOT}/benchmarks/astex_diverse/astex_diverse"

ASTEX_85 = [
    "1G9V", "1GM8", "1GPK", "1HNN", "1HP0", "1HQ2", "1IA1", "1IGJ", "1J3J", "1JD0",
    "1JJE", "1K3U", "1KE5", "1KZK", "1L2S", "1L7F", "1LPZ", "1M2Z", "1MEH", "1MQ6",
    "1N1M", "1N2J", "1N2V", "1N46", "1NAV", "1OF1", "1OF6", "1OPK", "1OQ5", "1OWE",
    "1P2Y", "1P62", "1PMN", "1Q1G", "1Q41", "1Q4G", "1R1H", "1R55", "1R58", "1R9O",
    "1S19", "1S3V", "1SG0", "1SJ0", "1SQ5", "1T40", "1T46", "1T9B", "1TT1", "1TW6",
    "1TZ8", "1U1C", "1U4D", "1UML", "1UNL", "1UOU", "1V0P", "1V48", "1V4S", "1VCJ",
    "1W1P", "1W2G", "1X8X", "1XM6", "1XOZ", "1Y6B", "1Y6R", "1YGC", "1YQY", "1YV3",
    "1YVF", "1YWR", "1Z95", "2BM2", "2BR1", "2BSM", "2BYS", "2C3I", "2CET", "2CGR",
    "2D3U", "2GBP", "2HB1", "2HR7", "2J62",
]

assert len(ASTEX_85) == 85, f"Expected 85 targets, got {len(ASTEX_85)}"


def make_dock_config(target: str) -> dict:
    """Oracle-ceiling dock_config with corrected GA params."""
    return {
        "flexibility": {
            "force_rigid": False,
            "intramolecular": True,
            "permeability": 0.9,
            "soft_wall_cutoff": 0.4,
            "receptor_rotamer_prep": False,
        },
        "optimization": {
            "grid_spacing": 0.375,
        },
        "scoring": {
            "normalize_area": True,
            "vct_dist_weight_r0": 7,          # FIXED: was 4 in v134
            "vct_normalize_contacts": False,
            "hbond_enabled": True,
            "hbond_search_enabled": True,
            "hbond_rank_enabled": False,
            "metal_coord_enabled": True,
            "sas_weight": 1.0,
            "tencom_weight": 0.0,
            "vct_entropy_weight": 0,
        },
        "seeding": {
            "mif_enabled": True,
        },
        "reference_ligand": {
            "file": f"{LIGAND_CACHE}/{target}/{target}_ligand.sdf",
            "seed_fraction": 0.9,             # oracle-ceiling: 90% pop seeded
            "pose_seed_enabled": True,         # oracle-ceiling: use crystal pose seed
            "k_nearest": 10,
        },
        "coarse_init": {
            "enabled": False,                  # oracle-ceiling: skip, we have oracle seed
            "grid_step": 3.0,
            "n_seeds": 25,
            "n_orientations": 16,
        },
        "thermodynamics": {
            "temperature": 300,
            "clustering_algorithm": "CF",
            "cluster_rmsd": 2.0,
            "classic_entropy_ranking": True,
            "force_cf_rank_emission": False,
        },
        "ga": {
            "num_chromosomes": 1000,           # FIXED: was 1750 in v134
            "num_generations": 875,
            "crossover_rate": 0.8,
            "mutation_rate": 0.03,
            "diversity_monitoring": True,
            "adaptive": True,
            "adaptive_k": [0.95, 0.1, 1.0, 0.05],
            "sharing_alpha": 4,                # FIXED: was 2.28571 in v134
            "boom_inject_interval": 100,
            "boom_inject_fraction": 1,         # oracle-ceiling: refill from oracle seed
            "n_elite": 1,
            "fitness_model": "SMFREE",
        },
    }


def main():
    # ── Pre-flight checks ────────────────────────────────────────────────────
    if not os.path.isfile(BINARY):
        sys.exit(f"FATAL: binary not found: {BINARY}")
    if not os.path.isfile(RUNNER):
        sys.exit(f"FATAL: runner not found: {RUNNER}")

    # ── Record git HEAD and binary SHA ───────────────────────────────────────
    git_head = subprocess.check_output(
        ["git", "-C", ROOT, "rev-parse", "HEAD"], text=True
    ).strip()
    binary_sha = subprocess.check_output(
        ["shasum", "-a", "256", BINARY], text=True
    ).split()[0]

    # ── Create output dir and pre-generate 85 dock_configs ───────────────────
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for target in ASTEX_85:
        tdir = os.path.join(OUTPUT_DIR, target)
        os.makedirs(tdir, exist_ok=True)
        with open(os.path.join(tdir, "dock_config.json"), "w") as fh:
            json.dump(make_dock_config(target), fh, indent=2)

    # ── Write provenance ─────────────────────────────────────────────────────
    provenance = {
        "dataset": "astex_diverse",
        "version": "v135",
        "binary_path": BINARY,
        "binary_sha256": binary_sha,
        "git_commit": git_head,
        "oracle_site_dir": ORACLE_DIR,
        "oracle_site_dir_set": True,
        "fix_notes": (
            "Pre-generated oracle-ceiling dock_configs fix v134 GA collapse: "
            "sharing_alpha=4 (was 2.28571), vct_dist_weight_r0=7 (was 4), "
            "num_chromosomes=1000 (was 1750), num_generations=875"
        ),
    }
    with open(os.path.join(OUTPUT_DIR, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=2)

    # ── Report to caller before fork ─────────────────────────────────────────
    print(f"[v135] git_commit  : {git_head}", flush=True)
    print(f"[v135] binary_sha  : {binary_sha}", flush=True)
    print(f"[v135] output_dir  : {OUTPUT_DIR}", flush=True)
    print(f"[v135] pid_file    : {PID_FILE}", flush=True)
    print(f"[v135] configs_written: {len(ASTEX_85)}", flush=True)
    sys.stdout.flush()

    # ── Build command ─────────────────────────────────────────────────────────
    cmd = [
        RUNNER,
        "--benchmark", "astex",
        "--mode",      "oracle-ceiling",
        "--threads",   "5",
        "--omp-threads", "2",
        "--output",    OUTPUT_DIR,
    ]

    env = {
        **os.environ,
        "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
        "FLEXAIDDS_BINARY":       BINARY,
        "FLEXAIDDS_DATA_DIR":     f"{REPRO}/engine",
        "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,
        # FLEXAIDDS_IGNORE_CACHE intentionally ABSENT → runner uses pre-built configs
    }

    # ── Double-fork / os.setsid() daemon launch ───────────────────────────────
    if os.fork() > 0:
        os._exit(0)        # parent exits; grandchild is fully detached

    os.setsid()

    if os.fork() > 0:
        os._exit(0)        # intermediate exits

    # Grandchild (daemon): redirect I/O and exec
    os.chdir(REPRO)
    log_fh = open(LOG_FILE, "ab", buffering=0)
    os.dup2(log_fh.fileno(), 1)
    os.dup2(log_fh.fileno(), 2)
    devnull = open(os.devnull, "rb")
    os.dup2(devnull.fileno(), 0)

    with open(PID_FILE, "w") as fh:
        fh.write(str(os.getpid()) + "\n")

    os.execve(cmd[0], cmd, env)


if __name__ == "__main__":
    main()
