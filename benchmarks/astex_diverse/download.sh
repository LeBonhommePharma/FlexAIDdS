#!/usr/bin/env bash
# download.sh — Acquire Astex Diverse Set tier-1 structures from RCSB
#
# Downloads receptor PDB files. Ligand extraction uses benchmark_datasets
# (--prepare-only flag) which calls extract_ligand() from HETATM records and
# writes apo receptor + ligand SDF into the cache tree.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
BUILD_DIR="${ROOT_DIR}/build"

TARGETS=(1sq5 2hb1 1r1h 1t46 2c69)
RCSB_BASE="https://files.rcsb.org/download"

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

mkdir -p "${DATA_DIR}"

echo "=== Astex Diverse Set Tier-1 Data Download ==="
echo "Targets: ${TARGETS[*]}"
echo ""

for pdb in "${TARGETS[@]}"; do
    pdb_file="${DATA_DIR}/${pdb}.pdb"
    if [[ -f "${pdb_file}" ]]; then
        echo "[OK] ${pdb}.pdb already exists"
    else
        echo "[..] Downloading ${pdb}.pdb from RCSB..."
        if curl -fsSL "${RCSB_BASE}/${pdb}.pdb" -o "${pdb_file}"; then
            echo "[OK] ${pdb}.pdb ($(wc -l < "${pdb_file}") lines)"
        else
            echo "[FAIL] Could not download ${pdb}.pdb" >&2
            exit 1
        fi
    fi
done

# Prepare ligands and apo receptors using benchmark_datasets if available.
if [[ -x "${BUILD_DIR}/benchmark_datasets" ]]; then
    echo ""
    echo "[..] Preparing apo receptors and ligands via benchmark_datasets..."
    FLEXAID_RESOLVED_BINARY="$(resolve_binary || true)"
    if [[ -n "${FLEXAID_RESOLVED_BINARY}" ]]; then
        export FLEXAIDDS_BINARY="${FLEXAID_RESOLVED_BINARY}"
    fi
    "${BUILD_DIR}/benchmark_datasets" \
        --benchmark astex \
        --prepare-only \
        --cache "${DATA_DIR}" \
        --output "${DATA_DIR}"

    for pdb in "${TARGETS[@]}"; do
        pdb_up="$(printf '%s' "${pdb}" | tr '[:lower:]' '[:upper:]')"
        src_dir="${DATA_DIR}/astex_diverse/${pdb_up}"
        if [[ -f "${src_dir}/${pdb_up}_apo.pdb" ]]; then
            cp "${src_dir}/${pdb_up}_apo.pdb" "${DATA_DIR}/${pdb}_apo.pdb"
        fi
        if [[ -f "${src_dir}/${pdb_up}_ligand.sdf" ]]; then
            cp "${src_dir}/${pdb_up}_ligand.sdf" "${DATA_DIR}/${pdb}_ligand.sdf"
        fi
        if [[ -f "${src_dir}/${pdb_up}_ligand.mol2" ]]; then
            cp "${src_dir}/${pdb_up}_ligand.mol2" "${DATA_DIR}/${pdb}_ligand.mol2"
        fi
    done

    PROVENANCE_TARGETS="$(IFS=,; echo "${TARGETS[*]}")"
    export PROVENANCE_TARGETS
    export PROVENANCE_BINARY="${FLEXAIDDS_BINARY:-}"
    export PROVENANCE_BINARY_SHA256="$(compute_sha256 "${FLEXAIDDS_BINARY:-}" 2>/dev/null || true)"
    export DATA_DIR BUILD_DIR
    python3 - <<'PY' > "${DATA_DIR}/prepare_provenance.json"
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
    "data_dir": os.environ.get("DATA_DIR"),
    "build_dir": os.environ.get("BUILD_DIR"),
    "targets": [t for t in os.environ.get("PROVENANCE_TARGETS", "").split(",") if t],
    "binary": {
        "path": os.environ.get("PROVENANCE_BINARY"),
        "sha256": sha256(os.environ.get("PROVENANCE_BINARY")),
    },
    "benchmark_datasets": str(Path(os.environ.get("BUILD_DIR", "")) / "benchmark_datasets"),
}
print(json.dumps(record, indent=2, sort_keys=False))
PY
else
    echo ""
    echo "[INFO] benchmark_datasets not built — build with ENABLE_BENCHMARK_DATASETS=ON"
    echo "       then re-run this script to extract ligands automatically."
    echo "       Alternatively, place {pdb}_apo.pdb and {pdb}_ligand.sdf files in ${DATA_DIR}/"
fi

echo ""
echo "Data directory: ${DATA_DIR}"
echo "Done."
