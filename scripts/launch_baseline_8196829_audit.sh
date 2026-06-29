#!/usr/bin/env bash
# launch_baseline_8196829_audit.sh — build + run published-commit baseline for CF audit.
set -euo pipefail

WT="/Users/lp.more/.grok/worktrees/projects-flexaidds/opus-48-baseline-8196829"
OUT="${HOME}/flexaidds_results/baseline_8196829_audit"
BUILD="${WT}/build_baseline_audit"
JSON="${WT}/benchmarks/datasets/benchmark_astex_native_85.json"
ORACLE="${WT}/benchmarks/astex_diverse/astex_diverse"
LOG="${OUT}/baseline_audit_runner.log"

mkdir -p "${OUT}"
exec > >(tee -a "${LOG}") 2>&1

echo "[baseline_8196829] started $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "[baseline_8196829] worktree=${WT}"
echo "[baseline_8196829] output=${OUT}"

cd "${WT}"
git log --oneline -1

cmake -B "${BUILD}" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build "${BUILD}" -j"$(sysctl -n hw.ncpu 2>/dev/null || nproc)" --target benchmark_datasets FlexAIDdS

export FLEXAIDDS_BINARY="${BUILD}/FlexAIDdS"
export FLEXAIDDS_BUILD="${BUILD}"
export FLEXAIDDS_REPO="${WT}"
export FLEXAIDDS_DATA_DIR="${BUILD}"
export FLEXAIDDS_ORACLE_SITE_DIR="${ORACLE}"
export FLEXAIDDS_THERMO=1
export FLEXAIDDS_T_EFF=0.596
export FLEXAIDDS_TENCOM_SCALE=1.0
export FLEXAIDDS_RESTARTS=7
export FLEXAIDDS_PARALLEL_RESTARTS=1
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
export FLEXAIDDS_CONSENSUS_SCORER=1
export FLEXAIDDS_SEED_ELITISM=1
export FLEXAIDDS_N_ELITE=1
export FLEXAIDDS_BUDGET_SCALE=1
export FLEXAIDDS_SOFTCORE_WAL=1
export FLEXAIDDS_SOFTCORE_FLOOR=0.5
export FLEXAIDDS_T_HOT=500
export FLEXAIDDS_NATIVE_SEED_FRAC=0.90
export FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1

"${BUILD}/benchmark_datasets" \
  --benchmark "crossdock_json:${JSON}" \
  --mode oracle-ceiling \
  --output "${OUT}" \
  --threads 4 \
  --omp-threads 2 \
  --temperature 298 \
  --job-timeout-seconds 7200

SCIENCE="/Users/lp.more/.grok/worktrees/projects-flexaidds/opus-48-science-fixes"
python3 "${SCIENCE}/scripts/cf_ground_truth_audit.py" "${OUT}" || true
python3 "${SCIENCE}/scripts/failure_classify.py" "${OUT}"
# Aggregate CSV for downstream compare (if not already emitted)
if [[ ! -f "${OUT}/astex_crossdock_85_results.csv" ]]; then
  python3 - "${OUT}" <<'PY'
import csv, glob, os, sys
out = sys.argv[1]
rows = []
for p in sorted(glob.glob(os.path.join(out, "*/result.csv"))):
    with open(p, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
            break
if rows:
    path = os.path.join(out, "astex_crossdock_85_results.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"[baseline_8196829] wrote {path} ({len(rows)} targets)")
PY
fi

echo "[baseline_8196829] finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"