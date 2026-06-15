#!/usr/bin/env python3
# launch_vH_guard.py — 9-target guard benchmark for virtual-H binary.
#
# vH = virtual-H architecture: H positions reconstructed from standard
# bond geometry (DonorGeom recipe stored per atom, read from live coords
# at scoring time). N.am uses VHG_AMIDE (planar external bisector) —
# correct in-plane direction, no false minima from scalar-0.3 fallback.
# Angular discrimination replaces blanket suppression from v57e.
#
# Changes from v57e (flexaid.h + hbond_potential.h + top.cpp):
#   - atom_struct gains vH_kind, vH_n, vH_nbr[2] (recipe fields)
#   - DonorGeom enum + build_virtual_H() + assign_virtual_h_geometry()
#   - donor_angle_term(): explicit H → virtual H → 0.0 (then 0.3 fallback)
#   - N.am restored to conservative_implicit_h_count (donor bit set again)
#
# Guard targets: 9-target subset used as regression gate before full-85.
# Must beat v50b (5/9) to qualify. Projected: 7/9.
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
import os, sys, subprocess, hashlib, json, datetime, shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO       = "/Users/lp.more/Projects/FlexAIDdS"
BUILD      = f"{REPO}/build_lto"
BINARY_SRC = f"{BUILD}/FlexAIDdS"
BINARY     = "/tmp/FlexAIDdS_vH"
RUNNER     = f"{BUILD}/benchmark_datasets"
DATA_DIR   = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85.json"
OUTPUT     = os.path.expanduser(
    "~/flexaidds_results/vH_guard_" +
    datetime.datetime.now().strftime("%Y%m%d_%H%M")
)
PROV_FILE  = f"{OUTPUT}/provenance.json"

# 9-target guard set — validated regression gate
GUARD_CODES = "1JD0,1MEH,1R55,1S3V,1SJ0,1X8X,1XM6,1XOZ,2D3U"

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

# ── Pre-flight ────────────────────────────────────────────────────────────────
for p in (BINARY_SRC, RUNNER, ORACLE_DIR, JSON_PAIRS, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing required path: {p}")

git_commit = subprocess.check_output(
    ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
).strip().split()[0]

shutil.copy2(BINARY_SRC, BINARY)
os.chmod(BINARY, 0o755)

engine_sha = sha256(BINARY)
matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

print(f"vH guard binary   : {engine_sha[:16]}...")
print(f"git commit        : {git_commit}")
print(f"guard targets     : {GUARD_CODES}")
print(f"output            : {OUTPUT}")

# ── Environment — same stack as v57e ─────────────────────────────────────────
env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY":                BINARY,
    "FLEXAIDDS_BUILD":                 BUILD,
    "FLEXAIDDS_REPO":                  REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
    "FLEXAIDDS_RESTARTS":              "3",
    "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
    "FLEXAIDDS_CONSENSUS_SCORER":      "1",
    "FLEXAIDDS_SEED_ELITISM":          "1",
    "FLEXAIDDS_N_ELITE":               "1",
    "FLEXAIDDS_BUDGET_SCALE":          "1",
    "FLEXAIDDS_SOFTCORE_WAL":          "1",
    "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
    "FLEXAIDDS_T_HOT":                 "500",
    "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
    "FLEXAIDDS_DATA_DIR":              DATA_DIR,
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
    "OMP_WAIT_POLICY":                 "passive",
    "OMP_PLACES":                      "cores",
    "OMP_PROC_BIND":                   "spread",
})
for k in (
    "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
    "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
    "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
    "FLEXAIDDS_RING_FLEX",
):
    env.pop(k, None)

# ── Command ───────────────────────────────────────────────────────────────────
cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark",           f"crossdock_json:{JSON_PAIRS}",
    "--output",              OUTPUT,
    "--threads",             "9",
    "--temperature",         "298",
    "--job-timeout-seconds", "1800",
    "--only-codes",          GUARD_CODES,
]

# ── Launch ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)

    print(f"\nLaunching vH guard benchmark (9 targets, 3 restarts) ...")
    print(f"  N.am donor:   VHG_AMIDE (planar external bisector, tracks GA moves)")
    print(f"  N.3 donor:    VHG_SP3_2NBR / VHG_SP3_1NBR (tetrahedral)")
    print(f"  O.3/S.3:      VHG_HYDROXYL (canonical 104.5° bend)")
    print(f"  H-bond cap:   -2.0 per pair")
    print(f"  sigma_angle:  30° (tuning knob; may need widening to 40-50°)")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":      "vH_guard",
        "launched_at":  datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":   git_commit,
        "description": (
            "Virtual-H architecture: DonorGeom recipe per atom, reconstructed from "
            "live heavy-neighbor coords at scoring time (tracks GA moves). N.am uses "
            "VHG_AMIDE (planar external bisector), N.3 uses SP3 tetrahedral, O.3/S.3 "
            "uses HYDROXYL canonical bend. angular discrimination replaces v57e blanket "
            "suppression. 3-restart consensus. Binary: FlexAIDdS_vH."
        ),
        "binary":        BINARY,
        "binary_sha256": engine_sha,
        "matrix_md5":    matrix_md5,
        "guard_codes":   GUARD_CODES,
        "oracle_site_dir": ORACLE_DIR,
        "output_dir":    OUTPUT,
        "pid":           child_pid,
        "sigma_angle_note": "sigma=30 from FA_Global default; if 1JD0/1S3V still fail consider 40-50",
    }
    with open(PROV_FILE, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\n✓ vH guard launched")
    print(f"  pid:     {child_pid}")
    print(f"  output:  {OUTPUT}")
    print(f"  monitor: tail -f {OUTPUT}/1JD0/stdout.log")
    print(f"  results: for t in 1JD0 1MEH 1R55 1S3V 1SJ0 1X8X 1XM6 1XOZ 2D3U; do")
    print(f"             awk -F, 'NR==2{{print FILENAME, $3}}' {OUTPUT}/$t/result.csv")
    print(f"           done")
