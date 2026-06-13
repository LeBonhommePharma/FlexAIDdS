#!/usr/bin/env python3
# launch_v40_clean.py — clean-oracle launcher for the v40 Astex Diverse 85 run.
#
# Purpose: a from-scratch, provenance-recording launch that CANNOT repeat the
# v35–v39 failure mode (0 oracle binding sites → pose-blinding → <5% success).
# The runner now hard-aborts on a 0-oracle Astex run (DatasetRunner Fix 1), and
# this launcher additionally sets FLEXAIDDS_ORACLE_SITE_DIR explicitly and
# records the exact binary / matrix / commit it launched with.
#
# Template: scripts/launch_v34_ctrl.py (same env contract, same threading).
# Differences from v34_ctrl:
#   * NO FLEXAIDDS_DATA_DIR override — uses the in-tree (default) MC_st0r5.2_6.dat
#   * FLEXAIDDS_SOFTCORE_FLOOR=0.5  (softcore_floor gate, per HEAD 8c0c840)
#   * Records binary SHA256 + matrix MD5 + git commit at launch time to
#     <OUTPUT>/launch_provenance.json (the runner also drops its own
#     provenance.json once it resolves the binary — see DatasetRunner Fix 4)
#   * nohup caffeinate -i launch pattern
#
# DOES NOT LAUNCH ANYTHING ON IMPORT. Run explicitly:  python3 scripts/launch_v40_clean.py
#
import os, sys, signal, subprocess, hashlib, json, datetime

REPO       = "/Users/lp.more/Projects/FlexAIDdS"
BUILD      = f"{REPO}/build_lto"
BINARY     = f"{BUILD}/FlexAIDdS"
RUNNER     = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"   # same as launch_v34_ctrl.py
MATRIX     = f"{REPO}/MC_st0r5.2_6.dat"                          # operative scoring matrix (cwd=REPO)
CODES_CSV  = f"{REPO}/benchmarks/astex_diverse/astex_diverse_set.csv"

DATESTAMP  = datetime.date.today().strftime("%Y%m%d")
OUTPUT     = os.path.expanduser(f"~/flexaidds_results/v40_{DATESTAMP}_cleanoracle")
CODES_FILE = f"{OUTPUT}/v40_codes_85.txt"
PIDFILE    = f"{OUTPUT}/v40_clean.pid"
STDOUT_LOG = f"{OUTPUT}/v40_clean_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v40_clean_stderr.log"
PROV_FILE  = f"{OUTPUT}/launch_provenance.json"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def git_commit():
    try:
        return subprocess.run(
            ["git", "-C", REPO, "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return ""


def make_codes_file(csv_path, out_path):
    codes = []
    with open(csv_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            pdb_id = line.split(",")[0].strip()
            if pdb_id:
                codes.append(pdb_id)
    with open(out_path, "w") as f:
        f.write("\n".join(codes) + "\n")
    return len(codes)


def main():
    # Preflight: every required input must exist before we fork.
    for p in (BINARY, RUNNER, ORACLE_DIR, CODES_CSV, MATRIX):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing {p}")

    # Guard: never start a second benchmark over a live one (shared-cache corruption).
    ps = subprocess.run(["pgrep", "-x", "benchmark_datasets"],
                        capture_output=True, text=True)
    if ps.stdout.strip():
        sys.exit(f"ERROR: benchmark_datasets already running (pids {ps.stdout.split()}) — abort")

    os.makedirs(OUTPUT, exist_ok=True)

    # ── Record launch-time provenance ────────────────────────────────────
    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(MATRIX)
    commit     = git_commit()
    provenance = {
        "run":              f"v40_{DATESTAMP}_cleanoracle",
        "launched_at":      datetime.datetime.now().isoformat(timespec="seconds"),
        "git_commit":       commit,
        "binary_path":      BINARY,
        "binary_sha256":    engine_sha,
        "runner_path":      RUNNER,
        "runner_sha256":    runner_sha,
        "matrix_path":      MATRIX,
        "matrix_md5":       matrix_md5,
        "oracle_site_dir":  ORACLE_DIR,
        "softcore_floor":   "0.5",
    }
    with open(PROV_FILE, "w") as f:
        json.dump(provenance, f, indent=2)
        f.write("\n")
    print("launch provenance:")
    print(json.dumps(provenance, indent=2))
    print(f"  → {PROV_FILE}")

    n_codes = make_codes_file(CODES_CSV, CODES_FILE)
    print(f"codes: {n_codes} targets → {CODES_FILE}")

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY":          BINARY,
        "FLEXAIDDS_BUILD":           BUILD,
        "FLEXAIDDS_REPO":            REPO,
        "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,       # ← prevents the v35–v39 0-oracle failure
        "FLEXAIDDS_SOFTCORE_FLOOR":  "0.5",
        "FLEXAIDDS_RESTARTS":        "3",
        "FLEXAIDDS_SEED_ELITISM":    "1",
        "FLEXAIDDS_N_ELITE":         "1",
        "FLEXAIDDS_BUDGET_SCALE":    "1",
        "FLEXAIDDS_SOFTCORE_WAL":    "1",
        "FLEXAIDDS_T_HOT":           "500",
        "OMP_WAIT_POLICY":           "passive",
        "OMP_PLACES":                "cores",
        "OMP_PROC_BIND":             "spread",
    })
    # Clear any stale diagnostic/ablation knobs and the matrix-swap override so
    # this run uses the in-tree default matrix.
    for k in ("FLEXAIDDS_DATA_DIR", "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
              "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
              "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
              "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
              "FLEXAIDDS_RING_FLEX"):
        env.pop(k, None)

    cmd = [
        "nohup", "caffeinate", "-i",
        RUNNER,
        "--benchmark",           "astex",
        "--only-codes",          CODES_FILE,
        "--output",              OUTPUT,
        "--threads",             "5",
        "--omp-threads",         "2",
        "--temperature",         "298",
        "--job-timeout-seconds", "5400",
    ]

    # ── nohup caffeinate -i launch, detached, PID recorded ────────────────
    os.chdir(REPO)
    with open(STDOUT_LOG, "w") as out, open(STDERR_LOG, "w") as err:
        p = subprocess.Popen(cmd, stdout=out, stderr=err, env=env,
                             start_new_session=True)
    with open(PIDFILE, "w") as pf:
        pf.write(str(p.pid) + "\n")

    print(f"launched v40 clean-oracle run → {OUTPUT}")
    print(f"  engine SHA: {engine_sha}")
    print(f"  matrix MD5: {matrix_md5}")
    print(f"  monitor:    tail -f {STDOUT_LOG}")
    print(f"  pid:        {PIDFILE}  ({p.pid})")
    print(f"  est. wall time: ~12.5 h")


if __name__ == "__main__":
    main()
