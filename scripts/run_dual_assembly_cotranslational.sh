#!/usr/bin/env bash
# run_dual_assembly_cotranslational.sh — DualAssembly cotranslational docking driver
#
# Wraps the `dual_assembly` C++ binary (built with -DENABLE_DUAL_ASSEMBLY_TOOL=ON)
# and parses its CSV output. The MVP backend is synthetic; --real-ga fails closed
# until the FlexAID GA hookup lands.
#
# Usage:
#   scripts/run_dual_assembly_cotranslational.sh \
#       --target-pdb <protofibril.pdb> \
#       --sequence <FASTA-1-letter or @file.fasta> \
#       [--monomer-pdb <monomer.pdb>] \
#       [--checkpoint-interval 10] \
#       [--sim-c-interval 5] \
#       [--monomer-conc-M 1e-6] \
#       [--temperature 310.15] \
#       [--threads 6] \
#       [--no-reciprocal-controls] \
#       [--no-sim-c] \
#       [--output-csv cotranslational_trajectory.csv] \
#       [--nascent-pdb-dir .]
#
# Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
set -euo pipefail

# ─── Defaults ────────────────────────────────────────────────────────────────
TARGET_PDB=""
SEQUENCE=""
MONOMER_PDB=""
CHECKPOINT_INTERVAL=10
SIM_C_INTERVAL=5
MONOMER_CONC_M="1e-6"
TEMPERATURE="310.15"
THREADS=6
RECIPROCAL_CONTROLS=1
SIM_C_ENABLED=1
OUTPUT_CSV="cotranslational_trajectory.csv"
NASCENT_PDB_DIR="."

# ─── arg parse ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-pdb)          TARGET_PDB="$2"; shift 2 ;;
        --sequence)            SEQUENCE="$2"; shift 2 ;;
        --monomer-pdb)         MONOMER_PDB="$2"; shift 2 ;;
        --checkpoint-interval) CHECKPOINT_INTERVAL="$2"; shift 2 ;;
        --sim-c-interval)      SIM_C_INTERVAL="$2"; shift 2 ;;
        --monomer-conc-M)      MONOMER_CONC_M="$2"; shift 2 ;;
        --temperature)         TEMPERATURE="$2"; shift 2 ;;
        --threads)             THREADS="$2"; shift 2 ;;
        --no-reciprocal-controls) RECIPROCAL_CONTROLS=0; shift ;;
        --no-sim-c)            SIM_C_ENABLED=0; shift ;;
        --output-csv)          OUTPUT_CSV="$2"; shift 2 ;;
        --nascent-pdb-dir)     NASCENT_PDB_DIR="$2"; shift 2 ;;
        -h|--help)
            sed -n '1,30p' "$0"
            exit 0
            ;;
        *) echo "unknown flag: $1" >&2; exit 2 ;;
    esac
done

# ─── Validate ────────────────────────────────────────────────────────────────
if [[ -z "$TARGET_PDB" || -z "$SEQUENCE" ]]; then
    echo "ERROR: --target-pdb and --sequence are required" >&2
    sed -n '1,30p' "$0" >&2
    exit 2
fi
if [[ ! -f "$TARGET_PDB" ]]; then
    echo "ERROR: target-pdb not readable: $TARGET_PDB" >&2
    exit 2
fi

# ─── Locate the binary ────────────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN=""
for candidate in \
    "$REPO_ROOT/build/dual_assembly" \
    "$REPO_ROOT/build/Release/dual_assembly" \
    "$REPO_ROOT/build/Debug/dual_assembly" \
    "$(command -v dual_assembly || true)"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then
        BIN="$candidate"
        break
    fi
done

if [[ -z "$BIN" ]]; then
    cat >&2 <<'EOF'
ERROR: dual_assembly binary not found.
Build it first:
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_DUAL_ASSEMBLY_TOOL=ON
  cmake --build build --target dual_assembly -j 6
EOF
    exit 1
fi

mkdir -p "$NASCENT_PDB_DIR"

# ─── Construct CLI args ──────────────────────────────────────────────────────
ARGS=(
    --target-pdb "$TARGET_PDB"
    --sequence   "$SEQUENCE"
    --checkpoint-interval "$CHECKPOINT_INTERVAL"
    --sim-c-interval      "$SIM_C_INTERVAL"
    --monomer-conc-M      "$MONOMER_CONC_M"
    --temperature         "$TEMPERATURE"
    --threads             "$THREADS"
    --output-csv          "$OUTPUT_CSV"
    --nascent-pdb-dir     "$NASCENT_PDB_DIR"
)

if [[ -n "$MONOMER_PDB" ]]; then
    if [[ ! -f "$MONOMER_PDB" ]]; then
        echo "ERROR: monomer-pdb not readable: $MONOMER_PDB" >&2
        exit 2
    fi
    ARGS+=( --monomer-pdb "$MONOMER_PDB" )
fi
[[ $RECIPROCAL_CONTROLS -eq 0 ]] && ARGS+=( --no-reciprocal-controls )
[[ $SIM_C_ENABLED       -eq 0 ]] && ARGS+=( --no-sim-c )

# OMP threads
export OMP_NUM_THREADS="$THREADS"

echo "[dual_assembly] binary:     $BIN"
echo "[dual_assembly] target-pdb: $TARGET_PDB"
echo "[dual_assembly] sequence:   ${SEQUENCE:0:60}$([[ ${#SEQUENCE} -gt 60 ]] && echo '…')"
echo "[dual_assembly] output:     $OUTPUT_CSV"
echo "[dual_assembly] threads:    $THREADS"

"$BIN" "${ARGS[@]}"

# ─── Summary ─────────────────────────────────────────────────────────────────
if [[ -f "$OUTPUT_CSV" ]]; then
    N_ROWS=$(($(wc -l < "$OUTPUT_CSV") - 1))
    echo "[dual_assembly] wrote $N_ROWS checkpoints to $OUTPUT_CSV"
fi
