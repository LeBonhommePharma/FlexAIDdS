#!/usr/bin/env bash
# Serial 3Dsig red-pair: A → B0 → B  (LOCAL-FIRST; sync iCloud later)
#
# Deck knobs: pop=1000 gen=2000 (2e6) restarts=10 · matrix 72d7
# Live I/O: ~/flexaidds_results (APFS). No CloudDocs during prepare/dock.
# Later: bash scripts/sync_three_engine_local_to_icloud.sh
#
# Usage:
#   bash scripts/run_3dsig_red_pair_serial.sh
#   bash scripts/run_3dsig_red_pair_serial.sh --from B0
#   bash scripts/run_3dsig_red_pair_serial.sh --only A
#   bash scripts/run_3dsig_red_pair_serial.sh --dry-run
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_local_first_benchmark_storage.sh"

FROM="A"
ONLY=""
DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="${2:-A}"; shift 2 ;;
    --only) ONLY="${2:-}"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n '1,18p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

export FLEXAID_POP="${FLEXAID_POP:-1000}"
export FLEXAID_GEN="${FLEXAID_GEN:-2000}"
export FLEXAID_RESTARTS="${FLEXAID_RESTARTS:-10}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export FLEXAIDDS_PARALLEL_RESTARTS="${FLEXAIDDS_PARALLEL_RESTARTS:-0}"

MATRIX_PIN="72d7c7396702331d96ff12d18f831796"
MAT="$FLEXAIDDS_QUEUE_ROOT/data/MC_st0r5.2_6.dat"
[[ -f "$MAT" ]] || MAT="$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat"
GOT=$(md5 -q "$MAT" 2>/dev/null || true)
[[ "$GOT" == "$MATRIX_PIN" ]] || { echo "FAIL: matrix MD5 '$GOT' != $MATRIX_PIN" >&2; exit 90; }
echo "OK: matrix md5=$GOT path=$MAT"

# Refuse live docks
if pgrep -f '/bin/[ABC]/FlexAID|/bin/C/FlexAIDdS|/bin/C/benchmark_datasets' >/dev/null 2>&1; then
  echo "REFUSE: live dock process" >&2
  pgrep -lf 'FlexAID|FlexAIDdS|benchmark_datasets' || true
  exit 92
fi

LOGDIR="$FLEXAIDDS_LOCAL_LOGDIR"
mkdir -p "$LOGDIR" "$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine"
TS=$(date -u +%Y-%m-%dT%H:%M:%SZ)
cat >"$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/3DSIG_RED_PAIR_PRIORITY.md" <<EOF
# 3Dsig red-pair (local-first)

**Updated:** $TS  
**Live OUT:** \$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/3dsig_r10/  
**Work:** \$FLEXAID_WORK_ROOT (local APFS)  
**iCloud:** deferred — \`bash scripts/sync_three_engine_local_to_icloud.sh\`  
**C0:** SUSPENDED  
**Knobs:** pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS  
EOF

run_arm() {
  local arm="$1"
  export FLEXAID_ARM_OUT="$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/$arm/3dsig_r10"
  mkdir -p "$FLEXAID_ARM_OUT"
  echo "=========================================="
  echo "3DSIG red-pair arm=$arm (LOCAL)"
  echo "OUT=$FLEXAID_ARM_OUT"
  echo "WORK=$FLEXAID_WORK_ROOT"
  echo "Q=$FLEXAIDDS_QUEUE_ROOT"
  echo "pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS"
  echo "=========================================="
  if (( DRY )); then
    bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" "$arm" --dry-run
    return 0
  fi
  # Prep flags: FLEXAID_FORCE_PREP=1 → --force; FLEXAID_SKIP_PREP=1 → --no-prepare
  extra=()
  [[ "${FLEXAID_FORCE_PREP:-0}" == "1" ]] && extra+=(--force)
  [[ "${FLEXAID_SKIP_PREP:-0}" == "1" ]] && extra+=(--no-prepare)
  NOHUP=0 bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" "$arm" "${extra[@]}"
}

ARMS=(A B0 B)
if [[ -n "$ONLY" ]]; then ARMS=("$ONLY")
else
  case "$FROM" in
    A) ARMS=(A B0 B) ;;
    B0) ARMS=(B0 B) ;;
    B) ARMS=(B) ;;
    *) echo "bad --from $FROM" >&2; exit 2 ;;
  esac
fi

CHAIN_LOG="$LOGDIR/run_3dsig_red_pair_serial.log"
CHAIN_PID="$LOGDIR/run_3dsig_red_pair_serial.pid"
if [[ -f "$CHAIN_PID" ]] && kill -0 "$(cat "$CHAIN_PID")" 2>/dev/null; then
  echo "REFUSE: serial chain already live pid=$(cat "$CHAIN_PID")" >&2
  exit 92
fi

if (( DRY )); then
  for a in "${ARMS[@]}"; do run_arm "$a"; done
  echo "DRY-RUN complete"
  exit 0
fi

{
  echo "$$" >"$CHAIN_PID"
  echo "[$TS] START local-first arms=${ARMS[*]} pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS"
  for a in "${ARMS[@]}"; do
    run_arm "$a" || echo "WARN: arm $a non-zero"
  done
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] COMPLETE red-pair serial (local)"
  # Best-effort thin iCloud sync (never fails the chain)
  bash "$ROOT/scripts/sync_three_engine_local_to_icloud.sh" --campaign 3dsig_r10 >>"$CHAIN_LOG" 2>&1 || true
  rm -f "$CHAIN_PID"
} >>"$CHAIN_LOG" 2>&1 &

echo "STARTED 3Dsig red-pair serial (LOCAL-FIRST) pid=$! log=$CHAIN_LOG"
echo "  arms: ${ARMS[*]}"
echo "  OUT:  $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/3dsig_r10/"
echo "  sync later: bash scripts/sync_three_engine_local_to_icloud.sh"
