#!/usr/bin/env bash
# download.sh — Acquire Astex Non-Native tier-1 structures from RCSB
#
# Downloads all PDB files needed for cross-docking pairs:
#   - Receptor PDBs (target_pdb)
#   - Ligand source PDBs (ligand_pdb, for ligand extraction)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
BUILD_DIR="${ROOT_DIR}/build"

# All unique PDB codes for tier-1 pairs (target + ligand source)
PDBS=(1hwi 1hww 1t40 1t46 2c68 2c69)
RCSB_BASE="https://files.rcsb.org/download"

mkdir -p "${DATA_DIR}"

echo "=== Astex Non-Native Tier-1 Data Download ==="
echo "PDB codes: ${PDBS[*]}"
echo ""

for pdb in "${PDBS[@]}"; do
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

# Extract ligands from ligand-source PDBs
LIGAND_PDBS=(1hww 1t46 2c69)
if [[ -x "${BUILD_DIR}/benchmark_datasets" ]]; then
    echo ""
    echo "[..] Extracting ligands via benchmark_datasets --prepare-only..."
    "${BUILD_DIR}/benchmark_datasets" \
        --benchmark astex_nonnative \
        --prepare-only \
        --cache "${DATA_DIR}" \
        --output "${DATA_DIR}" 2>/dev/null || true
    # Copy extracted ligands to flat layout expected by run.sh
    for lpdb in "${LIGAND_PDBS[@]}"; do
        src="${DATA_DIR}/${lpdb}/${lpdb}_ligand.mol2"
        dst="${DATA_DIR}/${lpdb}_ligand.mol2"
        if [[ -f "${src}" ]] && [[ ! -f "${dst}" ]]; then
            cp "${src}" "${dst}"
        fi
    done
else
    echo ""
    echo "[INFO] benchmark_datasets not built — ligand extraction skipped."
    echo "       Build with ENABLE_BENCHMARK_DATASETS=ON then re-run."
    echo "       Or place {pdb}_ligand.mol2 files manually in ${DATA_DIR}/"
fi

echo ""
echo "Data directory: ${DATA_DIR}"
echo "Done."
