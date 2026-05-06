#!/usr/bin/env bash
# download.sh — Acquire Astex Diverse Set tier-1 structures from RCSB
#
# Downloads receptor PDB files. Ligand MOL2 extraction uses benchmark_datasets
# (--prepare-only flag) which calls extract_ligand() from HETATM records.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
BUILD_DIR="${ROOT_DIR}/build"

TARGETS=(1sq5 2hb1 1r1h 1t46 2c69)
RCSB_BASE="https://files.rcsb.org/download"

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

# Extract ligands using benchmark_datasets if available
if [[ -x "${BUILD_DIR}/benchmark_datasets" ]]; then
    echo ""
    echo "[..] Extracting ligands via benchmark_datasets..."
    for pdb in "${TARGETS[@]}"; do
        ligand_out="${DATA_DIR}/${pdb}_ligand.mol2"
        if [[ ! -f "${ligand_out}" ]]; then
            "${BUILD_DIR}/benchmark_datasets" \
                --benchmark astex \
                --prepare-only \
                --cache "${DATA_DIR}" \
                --output "${DATA_DIR}" 2>/dev/null || true
            # Rename extracted file to expected location if needed
            if [[ -f "${DATA_DIR}/${pdb}/${pdb}_ligand.mol2" ]]; then
                cp "${DATA_DIR}/${pdb}/${pdb}_ligand.mol2" "${ligand_out}"
            fi
        fi
    done
else
    echo ""
    echo "[INFO] benchmark_datasets not built — build with ENABLE_BENCHMARK_DATASETS=ON"
    echo "       then re-run this script to extract ligands automatically."
    echo "       Alternatively, place {pdb}_ligand.mol2 files in ${DATA_DIR}/"
fi

echo ""
echo "Data directory: ${DATA_DIR}"
echo "Done."
