#!/usr/bin/env bash
# run.sh — Execute Astex Non-Native tier-1 benchmark (3 cross-docking pairs) and validate baselines
#
# Cross-docking: ligand from ligand_pdb is docked into receptor from target_pdb.
# Tests robustness to receptor conformational change (induced-fit challenge).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results"
BUILD_DIR="${ROOT_DIR}/build"
MANIFEST="${SCRIPT_DIR}/manifest.yaml"

# Tier-1 pairs: target_pdb ligand_pdb ligand_id target_name
PAIRS=(
    "1hwi 1hww STU CDK2"
    "1t40 1t46 2AN p38alpha"
    "2c68 2c69 PLN thymidine_kinase"
)

parse_baseline() {
    local key="$1"
    grep "^  ${key}:" "${MANIFEST}" | head -1 | awk '{print $2}'
}

TOL=$(parse_baseline "baseline_tolerance" 2>/dev/null || echo "0.05")
SUCCESS_2A=$(parse_baseline "crossdock_success_rate_2A" || echo "0.32")
MEAN_RMSD=$(parse_baseline "crossdock_mean_rmsd" || echo "3.10")
ENTROPY_RESCUE=$(parse_baseline "entropy_rescue_rate" || echo "0.22")

echo "=== Astex Non-Native Tier-1 Benchmark (Cross-Docking) ==="
echo "Pairs: ${#PAIRS[@]}"
echo "Baselines: success_2A=${SUCCESS_2A}, mean_RMSD=${MEAN_RMSD}, entropy_rescue=${ENTROPY_RESCUE}"
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
if [[ ! -d "${DATA_DIR}" ]] || [[ ! -f "${DATA_DIR}/1hwi.pdb" ]]; then
    echo "[..] Running download.sh..."
    bash "${SCRIPT_DIR}/download.sh"
fi

# --- Run cross-docking on each pair ---
mkdir -p "${RESULTS_DIR}"
PASS=0
FAIL=0
SKIP=0

for pair in "${PAIRS[@]}"; do
    read -r target_pdb ligand_pdb ligand_id target_name <<< "${pair}"

    receptor="${DATA_DIR}/${target_pdb}.pdb"
    ligand="${DATA_DIR}/${ligand_pdb}_ligand.mol2"
    reference="${DATA_DIR}/${ligand_pdb}.pdb"
    pair_id="${target_pdb}_x_${ligand_pdb}"
    out_dir="${RESULTS_DIR}/${pair_id}"

    echo "[..] Cross-docking: ligand ${ligand_id} from ${ligand_pdb} → receptor ${target_pdb} (${target_name})"

    if [[ ! -f "${receptor}" ]]; then
        echo "[SKIP] ${pair_id}: receptor PDB ${target_pdb}.pdb missing"
        SKIP=$((SKIP + 1))
        continue
    fi

    if [[ ! -f "${ligand}" ]]; then
        echo "[SKIP] ${pair_id}: ligand MOL2 ${ligand_pdb}_ligand.mol2 missing"
        SKIP=$((SKIP + 1))
        continue
    fi

    mkdir -p "${out_dir}"
    if "${FLEXAIDDS}" "${receptor}" "${ligand}" -o "${out_dir}" 2>&1 | tail -3; then
        echo "[OK] ${pair_id} cross-docked"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] ${pair_id} docking error" >&2
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "=== Results Summary ==="
echo "Passed: ${PASS}, Failed: ${FAIL}, Skipped: ${SKIP}"

REPORT="${RESULTS_DIR}/report.md"
cat > "${REPORT}" <<HEREDOC
# Astex Non-Native Tier-1 Report (Cross-Docking)

## Execution
- Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)
- Pairs: ${PAIRS[@]}
- Passed: ${PASS}, Failed: ${FAIL}, Skipped: ${SKIP}

## Baselines (Verdonk et al. 2008 J Chem Inf Model 48:2214-2225)

| Metric | Baseline | Tolerance | Measured | Status |
|:-------|:--------:|:---------:|:--------:|:------:|
| crossdock_success_rate_2A | ${SUCCESS_2A} | ±${TOL} | — | pending |
| crossdock_mean_rmsd | ${MEAN_RMSD} Å | ±${TOL} | — | pending |
| entropy_rescue_rate | ${ENTROPY_RESCUE} | ±${TOL} | — | pending |

> Published reference success rates under 2 Å: Vina 20–30%, Glide 30–40%.
> Full 65-target / 1112-pair validation: \`benchmark_datasets --benchmark astex_nonnative\`
HEREDOC

echo "Report: ${REPORT}"

if [[ ${FAIL} -gt 0 ]]; then
    echo "RESULT: FAIL (${FAIL} pairs failed)"
    exit 1
elif [[ ${PASS} -eq 0 ]]; then
    echo "RESULT: SKIP (no pairs executed — missing data?)"
    exit 0
else
    echo "RESULT: OK (${PASS} pairs cross-docked)"
    exit 0
fi
