#!/usr/bin/env bash
# Direct FlexAIDdS runner — bypasses benchmark_datasets, lower RAM footprint
# Usage: run_direct.sh <pdb_code> <cache_dir> <out_dir> [seed]
# Cache layout: <cache_dir>/<PDB>/<PDB>.cif + <PDB>_ligand.sdf
set -euo pipefail

PDB=${1:?need PDB code}
CACHE=${2:?need cache dir}
OUT=${3:?need out dir}
SEED=${4:-42}

REPO="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${REPO}/build/FlexAIDdS"
WRK="${REPO}/WRK"

ENTRY="${CACHE}/${PDB}"
REC="${ENTRY}/${PDB}.cif"
LIG="${ENTRY}/${PDB}_ligand.sdf"

# Skip if no structure
[[ ! -f "$REC" ]] && { echo "SKIP $PDB: no receptor"; exit 0; }
[[ ! -f "$LIG" ]] && { echo "SKIP $PDB: no ligand";   exit 0; }

# Skip if already done (has any output pose)
OUTD="${OUT}/${PDB}"
[[ -f "${OUTD}/result.csv" ]] && { echo "SKIP $PDB: done"; exit 0; }

mkdir -p "$OUTD"

"$BIN" "$REC" "$LIG" \
    --data-dir "$WRK" \
    --folded \
    -c "${OUTD}/dock_config.json" \
    -o "${OUTD}/${PDB}" \
    > "${OUTD}/stdout.log" 2>"${OUTD}/stderr.log"

EXIT=$?

# Write minimal result.csv
H_FINAL=$(grep "H_final" "${OUTD}/stdout.log" 2>/dev/null | tail -1 | grep -oE '[0-9.]+' | tail -1 || echo "0")
F_ENERGY=$(grep "Free energy F" "${OUTD}/stdout.log" 2>/dev/null | tail -1 | grep -oE '[-0-9.]+' | tail -1 || echo "0")
SUCCESS=$([[ $EXIT -eq 0 ]] && echo 1 || echo 0)

cat > "${OUTD}/result.csv" << CSV
pdb_id,predicted_dG,shannon_entropy,success,exit_code
${PDB},${F_ENERGY},${H_FINAL},${SUCCESS},${EXIT}
CSV

echo "DONE $PDB exit=$EXIT H=$H_FINAL dG=$F_ENERGY"
