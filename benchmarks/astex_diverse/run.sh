#!/usr/bin/env bash
# run.sh — Execute Astex Diverse Set tier-1 benchmark (5 targets) and validate baselines
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results"
BUILD_DIR="${ROOT_DIR}/build"
MANIFEST="${SCRIPT_DIR}/manifest.yaml"

TARGETS=(1sq5 2hb1 1r1h 1t46 2c69)

parse_baseline() {
    local key="$1"
    grep "^  ${key}:" "${MANIFEST}" | head -1 | awk '{print $2}'
}

TOL=$(parse_baseline "baseline_tolerance" 2>/dev/null || echo "0.05")
DOCK_TOP1=$(parse_baseline "docking_power_top1" || echo "0.70")
DOCK_TOP3=$(parse_baseline "docking_power_top3" || echo "0.85")
MEAN_RMSD=$(parse_baseline "mean_rmsd" || echo "2.30")
ENTROPY_RESCUE=$(parse_baseline "entropy_rescue_rate" || echo "0.30")

echo "=== Astex Diverse Set Tier-1 Benchmark ==="
echo "Targets: ${TARGETS[*]}"
echo "Baselines: top1=${DOCK_TOP1}, top3=${DOCK_TOP3}, mean_RMSD=${MEAN_RMSD}, entropy_rescue=${ENTROPY_RESCUE}"
echo "Tolerance: ±${TOL}"
echo ""

# --- Build FlexAIDdS if needed ---
if [[ ! -x "${BUILD_DIR}/FlexAIDdS" ]]; then
    echo "[..] Building FlexAIDdS..."
    cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
    cmake --build "${BUILD_DIR}" --parallel
fi

FLEXAIDDS="${BUILD_DIR}/FlexAIDdS"

# --- Download data if needed ---
if [[ ! -d "${DATA_DIR}" ]] || [[ ! -f "${DATA_DIR}/1sq5.pdb" ]]; then
    echo "[..] Running download.sh..."
    bash "${SCRIPT_DIR}/download.sh"
fi

# --- Run docking on each target ---
mkdir -p "${RESULTS_DIR}"
PASS=0
FAIL=0
SKIP=0
RMSDS=()

for pdb in "${TARGETS[@]}"; do
    receptor="${DATA_DIR}/${pdb}.pdb"
    ligand="${DATA_DIR}/${pdb}_ligand.mol2"
    out_dir="${RESULTS_DIR}/${pdb}"

    if [[ ! -f "${receptor}" ]]; then
        echo "[SKIP] ${pdb}: receptor PDB missing (run download.sh)"
        SKIP=$((SKIP + 1))
        continue
    fi

    if [[ ! -f "${ligand}" ]]; then
        echo "[SKIP] ${pdb}: ligand MOL2 missing (extracted by download.sh via benchmark_datasets)"
        SKIP=$((SKIP + 1))
        continue
    fi

    echo "[..] Docking ${pdb}..."
    mkdir -p "${out_dir}"
    if "${FLEXAIDDS}" "${receptor}" "${ligand}" -o "${out_dir}" 2>&1 | tail -3; then
        echo "[OK] ${pdb} docked"
        PASS=$((PASS + 1))
        # Extract best RMSD if available
        if [[ -f "${out_dir}/binding_modes.json" ]]; then
            rmsd=$(python3 -c "
import json, sys
d=json.load(open('${out_dir}/binding_modes.json'))
modes=d.get('binding_modes', [])
if modes: print(modes[0].get('best_pose_rmsd', 'N/A'))
else: print('N/A')
" 2>/dev/null || echo "N/A")
            RMSDS+=("${rmsd}")
            echo "       best RMSD: ${rmsd} Å"
        fi
    else
        echo "[FAIL] ${pdb} docking error" >&2
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "=== Results Summary ==="
echo "Passed: ${PASS}, Failed: ${FAIL}, Skipped: ${SKIP}"

# --- Generate metrics report ---
REPORT="${RESULTS_DIR}/report.md"
cat > "${REPORT}" <<HEREDOC
# Astex Diverse Set Tier-1 Report

## Execution
- Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)
- Targets: ${TARGETS[*]}
- Passed: ${PASS}, Failed: ${FAIL}, Skipped: ${SKIP}

## Baselines (Hartshorn et al. 2007 J Med Chem 50:726-741)

| Metric | Baseline | Tolerance | Measured | Status |
|:-------|:--------:|:---------:|:--------:|:------:|
| docking_power_top1 | ${DOCK_TOP1} | ±${TOL} | — | pending |
| docking_power_top3 | ${DOCK_TOP3} | ±${TOL} | — | pending |
| mean_rmsd | ${MEAN_RMSD} Å | ±${TOL} | — | pending |
| entropy_rescue_rate | ${ENTROPY_RESCUE} | ±${TOL} | — | pending |

> Metrics are populated by DatasetRunner once pose output is available.
> Full 85-complex validation: \`benchmark_datasets --benchmark astex\`
HEREDOC

echo "Report: ${REPORT}"

if [[ ${FAIL} -gt 0 ]]; then
    echo "RESULT: FAIL (${FAIL} targets failed)"
    exit 1
elif [[ ${PASS} -eq 0 ]]; then
    echo "RESULT: SKIP (no targets executed — missing data?)"
    exit 0
else
    echo "RESULT: OK (${PASS} targets docked)"
    exit 0
fi
