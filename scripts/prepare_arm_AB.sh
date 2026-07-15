#!/usr/bin/env bash
# Prepare FlexAID work trees for arms A / B0 / B (ProcessLigand + CONFIG + ga).
#
# Usage:
#   scripts/prepare_arm_AB.sh                 # pilot8, all three arms
#   scripts/prepare_arm_AB.sh --pdb 1GPK
#   scripts/prepare_arm_AB.sh --full85
#   scripts/prepare_arm_AB.sh --dry-run
set -euo pipefail

# shellcheck disable=SC1090
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"

ICLOUD_DEFAULT="${HOME}/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$ICLOUD_DEFAULT}"
Q="${QUEUE_ROOT:-${FLEXAIDDS_QUEUE_ROOT:-$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1}}"
export FLEXAIDDS_QUEUE_ROOT="$Q"

REPO="${FLEXAIDDS_ROOT:-}"
if [[ -z "$REPO" ]]; then
  REPO="$(cd "$(dirname "$0")/.." && pwd)"
fi

if [[ -z "${FLEXAIDDS_PROCESSLIGAND:-}" ]]; then
  for m in "$REPO"/.venv-processligand/lib/python*/site-packages/processligandpy/bin/ProcessLigand; do
    if [[ -x "$m" ]]; then
      export FLEXAIDDS_PROCESSLIGAND="$m"
      break
    fi
  done
fi

ARGS=(--queue-root "$Q" --work-root "$Q/work" --arms A,B0,B --pop 1000 --gen 6000 --restarts 5)
MODE=pilot8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) ARGS+=(--dry-run); shift ;;
    --force) ARGS+=(--force); shift ;;
    --pdb) ARGS+=(--pdb "$2"); MODE=single; shift 2 ;;
    --full85) MODE=full; shift ;;
    --arms) ARGS+=(--arms "$2"); shift 2 ;;
    *) echo "Unknown: $1" >&2; exit 2 ;;
  esac
done

if [[ "$MODE" == "pilot8" ]]; then
  ARGS+=(--pilot8)
fi

echo "ProcessLigand=${FLEXAIDDS_PROCESSLIGAND:-unset}"
exec python3 "$REPO/scripts/generate_flexaid_inp.py" "${ARGS[@]}"
