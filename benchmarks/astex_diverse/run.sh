#!/usr/bin/env bash
# run.sh — Execute Astex Diverse Set tier-1 benchmark (5 targets) and validate baselines
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results"
BUILD_DIR="${ROOT_DIR}/build"
MANIFEST="${SCRIPT_DIR}/manifest.yaml"
SEED="${FLEXAID_SEED:-42}"

TARGETS=(1sq5 2hb1 1r1h 1t46 2c69)

compute_sha256() {
    local file="$1"
    if [[ ! -f "${file}" ]]; then
        return 1
    fi
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "${file}" | awk '{print $1}'
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "${file}" | awk '{print $1}'
    else
        python3 - "${file}" <<'PY'
import hashlib
import sys
from pathlib import Path
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
    fi
}

resolve_binary() {
    if [[ -n "${FLEXAIDDS_BINARY:-}" && -x "${FLEXAIDDS_BINARY}" ]]; then
        printf '%s\n' "${FLEXAIDDS_BINARY}"
        return 0
    fi
    if [[ -n "${FLEXAID_BINARY:-}" && -x "${FLEXAID_BINARY}" ]]; then
        printf '%s\n' "${FLEXAID_BINARY}"
        return 0
    fi
    if [[ -x "${BUILD_DIR}/FlexAIDdS" ]]; then
        printf '%s\n' "${BUILD_DIR}/FlexAIDdS"
        return 0
    fi
    return 1
}

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
echo "Seed: ${SEED}"
echo ""

# --- Resolve / build FlexAIDdS if needed ---
FLEXAIDDS="$(resolve_binary || true)"
if [[ -z "${FLEXAIDDS}" ]]; then
    echo "[..] Building FlexAIDdS..."
    cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
    cmake --build "${BUILD_DIR}" --parallel
    FLEXAIDDS="${BUILD_DIR}/FlexAIDdS"
fi

if [[ ! -x "${FLEXAIDDS}" ]]; then
    echo "[FAIL] FlexAIDdS binary not found after build/override resolution" >&2
    exit 1
fi

export FLEXAIDDS_BINARY="${FLEXAIDDS}"
export FLEXAID_SEED="${SEED}"

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
    pdb_up="$(printf '%s' "${pdb}" | tr '[:lower:]' '[:upper:]')"
    prepared_dir="${DATA_DIR}/astex_diverse/${pdb_up}"

    receptor="${prepared_dir}/${pdb_up}_apo.pdb"
    if [[ ! -f "${receptor}" ]]; then
        receptor="${DATA_DIR}/${pdb}_apo.pdb"
    fi
    if [[ ! -f "${receptor}" ]]; then
        receptor="${DATA_DIR}/${pdb}.pdb"
    fi

    ligand="${prepared_dir}/${pdb_up}_ligand.sdf"
    if [[ ! -f "${ligand}" ]]; then
        ligand="${DATA_DIR}/${pdb}_ligand.sdf"
    fi
    if [[ ! -f "${ligand}" ]]; then
        ligand="${DATA_DIR}/${pdb}_ligand.mol2"
    fi

    out_dir="${RESULTS_DIR}/${pdb}"

    if [[ ! -f "${receptor}" ]]; then
        echo "[SKIP] ${pdb}: receptor PDB missing (run download.sh)"
        SKIP=$((SKIP + 1))
        continue
    fi

    if [[ ! -f "${ligand}" ]]; then
        echo "[SKIP] ${pdb}: ligand SDF/MOL2 missing (extracted by download.sh via benchmark_datasets)"
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

PROVENANCE_TARGETS="$(IFS=,; echo "${TARGETS[*]}")"
export PROVENANCE_TARGETS
export PROVENANCE_BINARY="${FLEXAIDDS}"
export PROVENANCE_MANIFEST="${MANIFEST}"
export PROVENANCE_DATA_DIR="${DATA_DIR}"
export PASS_COUNT="${PASS}"
export FAIL_COUNT="${FAIL}"
export SKIP_COUNT="${SKIP}"
python3 - <<'PY' > "${RESULTS_DIR}/provenance.json"
import hashlib
import json
import os
import time
from pathlib import Path

def sha256(path):
    if not path:
        return None
    p = Path(path)
    if not p.is_file():
        return None
    return hashlib.sha256(p.read_bytes()).hexdigest()

record = {
    "benchmark": "astex_diverse_tier1",
    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "seed": int(os.environ.get("FLEXAID_SEED", "42") or 42),
    "targets": [t for t in os.environ.get("PROVENANCE_TARGETS", "").split(",") if t],
    "data_dir": os.environ.get("PROVENANCE_DATA_DIR"),
    "manifest": os.environ.get("PROVENANCE_MANIFEST"),
    "binary": {
        "path": os.environ.get("PROVENANCE_BINARY"),
        "sha256": sha256(os.environ.get("PROVENANCE_BINARY")),
    },
    "counts": {
        "passed": int(os.environ.get("PASS_COUNT", "0") or 0),
        "failed": int(os.environ.get("FAIL_COUNT", "0") or 0),
        "skipped": int(os.environ.get("SKIP_COUNT", "0") or 0),
    },
}
print(json.dumps(record, indent=2, sort_keys=False))
PY

# --- Generate metrics report ---
REPORT="${RESULTS_DIR}/report.md"
cat > "${REPORT}" <<HEREDOC
# Astex Diverse Set Tier-1 Report

## Execution
- Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)
- Targets: ${TARGETS[*]}
- Seed: ${SEED}
- Binary: ${FLEXAIDDS}
- Binary SHA256: $(compute_sha256 "${FLEXAIDDS}")
- Manifest: ${MANIFEST}
- Passed: ${PASS}, Failed: ${FAIL}, Skipped: ${SKIP}

## Baselines (manifest + published JCIM comparator)

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
