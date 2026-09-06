#!/usr/bin/env bash
# Isolated two-engine parity per METHODOLOGY.md section 1.
# All scientific inputs and both source/build identities are explicit.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: gate_parity.sh ENGINE_A ENGINE_B \
  --baseline-source DIR --baseline-build DIR \
  --candidate-source DIR --candidate-build DIR \
  --receptor FILE --ligand FILE --config FILE --data-dir DIR --out NEW_DIR

No input or matrix fallback is inferred from a mutable checkout. NEW_DIR must
not exist. Engines with the same basename receive distinct P0/P1 output trees.
EOF
}
if [[ ${1:-} == --help ]]; then usage; exit 0; fi
if [[ $# -lt 2 ]]; then usage >&2; exit 2; fi
engine_a=$1; engine_b=$2; shift 2
baseline_source= baseline_build= candidate_source= candidate_build=
receptor= ligand= config= data_dir= output=
while [[ $# -gt 0 ]]; do
  if [[ $# -lt 2 ]]; then usage >&2; exit 2; fi
  case "$1" in
    --baseline-source) baseline_source=$2 ;;
    --baseline-build) baseline_build=$2 ;;
    --candidate-source) candidate_source=$2 ;;
    --candidate-build) candidate_build=$2 ;;
    --receptor) receptor=$2 ;;
    --ligand) ligand=$2 ;;
    --config) config=$2 ;;
    --data-dir) data_dir=$2 ;;
    --out) output=$2 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
  shift 2
done
for value in "$baseline_source" "$baseline_build" "$candidate_source" "$candidate_build" \
             "$receptor" "$ligand" "$config" "$data_dir" "$output"; do
  if [[ -z "$value" ]]; then usage >&2; exit 2; fi
done
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-python3}
if [[ -e "$output" ]]; then printf 'Output already exists: %s\n' "$output" >&2; exit 2; fi
mkdir -p -- "$(dirname -- "$output")"
mkdir -- "$output"
"$python_bin" "$script_dir/engine_repro_gate.py" run-one --label P0 \
  --engine "$engine_a" --source-dir "$baseline_source" --build-dir "$baseline_build" \
  --receptor "$receptor" --ligand "$ligand" --config "$config" --data-dir "$data_dir" \
  --out "$output/P0" --omp-threads 1 --parallel-reproduce off --require-gen0
"$python_bin" "$script_dir/engine_repro_gate.py" run-one --label P1 \
  --engine "$engine_b" --source-dir "$candidate_source" --build-dir "$candidate_build" \
  --receptor "$receptor" --ligand "$ligand" --config "$config" --data-dir "$data_dir" \
  --out "$output/P1" --omp-threads 1 --parallel-reproduce off --require-gen0
"$python_bin" "$script_dir/engine_repro_gate.py" compare --kind parity \
  --runs "$output/P0" "$output/P1" --json "$output/parity.json"
