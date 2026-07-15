#!/usr/bin/env bash
# Serial 3Dsig red-pair campaign: A (FlexAID CF) → B0 (master CF) → B (entropy).
# Deck knobs: pop=1000 gen=2000 (2e6 evals) restarts=10 · matrix 72d7 · no dual-launch.
#
# Usage:
#   bash scripts/run_3dsig_red_pair_serial.sh              # A then B0 then B
#   bash scripts/run_3dsig_red_pair_serial.sh --from B0    # resume from B0
#   bash scripts/run_3dsig_red_pair_serial.sh --only A     # single arm
#   bash scripts/run_3dsig_red_pair_serial.sh --dry-run
#
# OUT namespace (separate from old pilot8 R=5):
#   $FLEXAIDDS_RESULTS/campaigns/three_engine/{A,B0,B}/3dsig_r10/
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_icloud_benchmark_storage.sh"

Q="${FLEXAIDDS_QUEUE_ROOT}"
LOGDIR="$Q/logs"
mkdir -p "$LOGDIR"

FROM="A"
ONLY=""
DRY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="${2:-A}"; shift 2 ;;
    --only) ONLY="${2:-}"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help)
      sed -n '1,20p' "$0"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

# Deck knobs
export FLEXAID_POP="${FLEXAID_POP:-1000}"
export FLEXAID_GEN="${FLEXAID_GEN:-2000}"
export FLEXAID_RESTARTS="${FLEXAID_RESTARTS:-10}"
export NOHUP="${NOHUP:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export FLEXAIDDS_PARALLEL_RESTARTS="${FLEXAIDDS_PARALLEL_RESTARTS:-0}"

MATRIX_PIN="72d7c7396702331d96ff12d18f831796"
MAT="$Q/data/MC_st0r5.2_6.dat"
GOT=$(md5 -q "$MAT" 2>/dev/null || true)
if [[ "$GOT" != "$MATRIX_PIN" ]]; then
  echo "FAIL: matrix MD5 '$GOT' != $MATRIX_PIN ($MAT)" >&2
  exit 90
fi
echo "OK: matrix md5=$GOT"

# Refuse dual-launch with live C0 / FlexAIDdS / FlexAID docks
refuse_if_live() {
  local label="$1" pidfile="$2"
  if [[ -f "$pidfile" ]]; then
    local pid
    pid=$(tr -d ' \n' <"$pidfile" || true)
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
      echo "REFUSE: $label still live pid=$pid ($pidfile)" >&2
      exit 92
    fi
  fi
}
refuse_if_live "C0_claim" "$LOGDIR/C0_claim_clean.pid"
refuse_if_live "C0_full85" "$LOGDIR/C0_full85.pid"

if pgrep -f '/bin/C/FlexAIDdS|/bin/C/benchmark_datasets' >/dev/null 2>&1; then
  echo "REFUSE: live FlexAIDdS/DatasetRunner dock process" >&2
  pgrep -lf 'FlexAIDdS|benchmark_datasets' || true
  exit 92
fi
if pgrep -f '/bin/[AB]/]/FlexAID' >/dev/null 2>&1; then
  echo "REFUSE: live FlexAID dock process" >&2
  pgrep -lf 'FlexAID' || true
  exit 92
fi

# Clear stale C0 locks (dead pids only)
for L in "$LOGDIR/C0_claim_clean.lock" "$LOGDIR/C0_full85.lock"; do
  if [[ -f "$L" ]]; then
    pid=$(tr -d ' \n' <"$L" || true)
    if [[ -z "$pid" ]] || ! kill -0 "$pid" 2>/dev/null; then
      rm -f "$L"
      echo "removed stale lock $(basename "$L")"
    fi
  fi
done

# Priority markers
TS=$(date -u +%Y-%m-%dT%H:%M:%SZ)
mkdir -p "$FLEXAIDDS_RESULTS/campaigns/three_engine"
cat >"$FLEXAIDDS_RESULTS/campaigns/three_engine/3DSIG_RED_PAIR_PRIORITY.md" <<EOF
# 3Dsig red-pair priority

**Updated:** $TS  
**Status:** C0 STOPPED/SUSPENDED · A→B0→B PRIORITY  
**Protocol:** docs/implementation/3dsig_red_pair_protocol.md  
**Knobs:** pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS (2e6 evals/sim)  
**Metric:** S_top10 RMSD<2Å · 10k bootstrap median (not S1/BCR/PB for 2017 bars)  
**Matrix:** $MATRIX_PIN  
EOF
# queue status hint
if [[ -d "$Q" ]]; then
  cat >"$Q/3DSIG_RED_PAIR_ACTIVE.md" <<EOF
# 3Dsig red-pair ACTIVE

C0 claim/legacy: **SUSPENDED** (do not relaunch without operator).  
Serial arms: **A → B0 → B** under deck knobs.  
Started: $TS
EOF
fi

run_arm() {
  local arm="$1"
  export FLEXAID_ARM_OUT="$FLEXAIDDS_RESULTS/campaigns/three_engine/$arm/3dsig_r10"
  mkdir -p "$FLEXAID_ARM_OUT"
  echo "=========================================="
  echo "3DSIG red-pair arm=$arm OUT=$FLEXAID_ARM_OUT"
  echo "pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS"
  echo "=========================================="
  if (( DRY )); then
    bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" "$arm" --dry-run
    return 0
  fi
  # Foreground serial chain: wait for each arm before next
  # NOHUP=0 for serial wait inside this script; outer nohup wraps whole chain
  NOHUP=0 bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" "$arm"
}

ARMS=(A B0 B)
if [[ -n "$ONLY" ]]; then
  ARMS=("$ONLY")
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
  echo "[$TS] START arms=${ARMS[*]} pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS"
  for a in "${ARMS[@]}"; do
    run_arm "$a" || echo "WARN: arm $a returned non-zero"
  done
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] COMPLETE red-pair serial"
  rm -f "$CHAIN_PID"
} >>"$CHAIN_LOG" 2>&1 &

echo "STARTED 3Dsig red-pair serial pid=$! log=$CHAIN_LOG"
echo "  arms: ${ARMS[*]}"
echo "  OUT:  \$FLEXAIDDS_RESULTS/campaigns/three_engine/{A,B0,B}/3dsig_r10/"
echo "  protocol: docs/implementation/3dsig_red_pair_protocol.md"
