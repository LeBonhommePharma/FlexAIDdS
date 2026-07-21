#!/usr/bin/env python3
"""v136 — proper blind docking benchmark (autonomous mode, seeding OFF).

v135 was an oracle-ceiling run (seed_elitism=ON, blinding=OFF) where
pose_seed_enabled=True / seed_fraction=0.9 fed the crystal ligand directly
into the GA population, making it a ceiling measurement rather than a
blind prediction benchmark.

v136 fixes this by:
  1. Running the benchmark_datasets runner in '--mode autonomous'
     (seed_elitism=OFF, blinding=ON at the runner level).
  2. Pre-generating all 85 dock_configs with:
       pose_seed_enabled = false   ← no crystal pose seed
       seed_fraction      = 0.0    ← no oracle fraction
       coarse_init.enabled = true  ← GA must initialize blindly via grid
       boom_inject_fraction = 0    ← no oracle refill during GA cycles
  3. Keeping FLEXAIDDS_ORACLE_SITE_DIR set so the runner can compute
     RMSD against crystal poses for evaluation (read-only, not for seeding).

GA params: identical to v135 (sharing_alpha=4, vct_dist_weight_r0=7,
num_chromosomes=1000, num_generations=875) — isolates seeding as the
only variable change between v135 and v136.

Binding site files audit (pre-launch):
  All 85 targets have non-degenerate _binding_site.pdb files (32–381 atoms).
  No rcut/sphere_radius parameter exists in the JSON config format;
  binding site extent is derived from _binding_site.pdb atom coordinates
  by the runner. No binding site issues found.

C++ binary: build/FlexAIDdS (SHA 36cbede1…) — same as v135, no rebuild.
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
OUTPUT_DIR = f"{RESULTS}/v136_full85_blind_{TIMESTAMP}"
PID_FILE   = f"{RESULTS}/v136_full85_blind_{TIMESTAMP}_launcher.pid"
LOG_FILE   = f"{RESULTS}/v136_full85_blind_{TIMESTAMP}_launcher.log"

BINARY       = f"{ROOT}/build/FlexAIDdS"
RUNNER       = f"{REPRO}/engine/benchmark_datasets"
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
    """Blind-docking dock_config: seeding OFF, coarse_init ON.

    Key differences from v135 oracle-ceiling config:
      - pose_seed_enabled : false  (was true)
      - seed_fraction     : 0.0   (was 0.9)
      - coarse_init.enabled: true  (was false) — blind GA init via grid
      - boom_inject_fraction: 0   (was 1)      — no oracle refill
    All GA quality params (sharing_alpha, vct_dist_weight_r0, etc.) unchanged.
    """
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
            "vct_dist_weight_r0": 7,
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
            "seed_fraction": 0.0,          # BLIND: no oracle fraction (was 0.9)
            "pose_seed_enabled": False,    # BLIND: no crystal pose seed (was True)
            "k_nearest": 10,
        },
        "coarse_init": {
            "enabled": True,               # BLIND: grid init required (was False)
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
            "num_chromosomes": 1000,
            "num_generations": 875,
            "crossover_rate": 0.8,
            "mutation_rate": 0.03,
            "diversity_monitoring": True,
            "adaptive": True,
            "adaptive_k": [0.95, 0.1, 1.0, 0.05],
            "sharing_alpha": 4,
            "boom_inject_interval": 100,
            "boom_inject_fraction": 0,     # BLIND: no oracle refill (was 1)
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

    # ── Verify all 85 binding site PDB files are present and non-degenerate ─
    missing_sites, degenerate_sites = [], []
    for target in ASTEX_85:
        bsite = os.path.join(ORACLE_DIR, target, f"{target}_binding_site.pdb")
        if not os.path.isfile(bsite):
            missing_sites.append(target)
        else:
            with open(bsite) as fh:
                natoms = sum(1 for l in fh if l.startswith(("ATOM", "HETATM")))
            if natoms <= 3:
                degenerate_sites.append(f"{target}(n={natoms})")
    if missing_sites:
        sys.exit(f"FATAL: missing _binding_site.pdb for: {missing_sites}")
    if degenerate_sites:
        sys.exit(f"FATAL: degenerate binding sites (≤3 atoms): {degenerate_sites}")

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
        "version": "v136",
        "mode": "autonomous",       # REQUIRED for blind: oracle-ceiling anchors IC frame to crystal pose regardless of pose_seed_enabled
        "binary_path": BINARY,
        "binary_sha256": binary_sha,
        "git_commit": git_head,
        "oracle_site_dir": ORACLE_DIR,
        "oracle_site_dir_set": True,
        "blind_notes": (
            "Blind docking benchmark: pose_seed_enabled=false, seed_fraction=0.0, "
            "coarse_init.enabled=true, boom_inject_fraction=0. "
            "Runner mode=autonomous (required: oracle-ceiling injects [ORACLE-SEED] and anchors "
            "IC frame to crystal pose at runner level, bypassing pose_seed_enabled=false). "
            "autonomous mode uses per-target adaptive GA params (runner-calibrated for blind docking). "
            "Oracle site dir retained for RMSD evaluation only. "
            "All GA quality params identical to v135."
        ),
    }
    with open(os.path.join(OUTPUT_DIR, "provenance.json"), "w") as fh:
        json.dump(provenance, fh, indent=2)

    # ── Report to caller before fork ─────────────────────────────────────────
    print(f"[v136] git_commit      : {git_head}", flush=True)
    print(f"[v136] binary_sha      : {binary_sha}", flush=True)
    print(f"[v136] output_dir      : {OUTPUT_DIR}", flush=True)
    print(f"[v136] pid_file        : {PID_FILE}", flush=True)
    print(f"[v136] configs_written : {len(ASTEX_85)}", flush=True)
    print(f"[v136] mode            : autonomous (IC-frame blind; oracle-ceiling anchors IC to crystal regardless of config)", flush=True)
    print(f"[v136] pose_seed       : OFF (seed_fraction=0.0, pose_seed_enabled=false)", flush=True)
    sys.stdout.flush()

    # ── Build command (autonomous = blind; IC frame NOT anchored to crystal) ─
    cmd = [
        RUNNER,
        "--benchmark", "astex",
        "--mode",      "autonomous",      # REQUIRED: oracle-ceiling logs [ORACLE-SEED] and anchors IC frame to crystal regardless of dock_config
        "--threads",   "5",
        "--omp-threads", "2",
        "--output",    OUTPUT_DIR,
    ]

    env = {
        **os.environ,
        "PATH": "/opt/homebrew/bin:" + os.environ.get("PATH", ""),
        "FLEXAIDDS_BINARY":          BINARY,
        "FLEXAIDDS_DATA_DIR":        f"{REPRO}/engine",
        "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,  # for RMSD eval only
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
