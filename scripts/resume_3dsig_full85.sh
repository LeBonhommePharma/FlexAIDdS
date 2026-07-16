#!/usr/bin/env bash
# Resume / continue / multi-restart extension for 3Dsig full85 red-pair.
#
# First pass (R=1) is launched by: bash scripts/run_3dsig_red_pair_full85.sh
# This helper:
#   - refuses dual-launch if full85 chain or FlexAID --legacy is live
#   - resumes arms that skip existing result.csv (unless --force)
#   - optionally raises FLEXAID_RESTARTS for a second pass (does not wipe R0)
#
# Usage:
#   bash scripts/resume_3dsig_full85.sh --status
#   bash scripts/resume_3dsig_full85.sh --from B0
#   bash scripts/resume_3dsig_full85.sh --only B --restarts 10
#   FLEXAID_FORCE_PREP=0 bash scripts/resume_3dsig_full85.sh --from A
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_local_first_benchmark_storage.sh"

FROM="A"
ONLY=""
STATUS=0
RESTARTS="${FLEXAID_RESTARTS:-1}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="${2:-A}"; shift 2 ;;
    --only) ONLY="${2:-}"; shift 2 ;;
    --restarts) RESTARTS="${2:-1}"; shift 2 ;;
    --status) STATUS=1; shift ;;
    -h|--help) sed -n '1,20p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

export FLEXAID_POP="${FLEXAID_POP:-1000}"
export FLEXAID_GEN="${FLEXAID_GEN:-2000}"
export FLEXAID_RESTARTS="$RESTARTS"
export FLEXAID_CAMPAIGN="${FLEXAID_CAMPAIGN:-3dsig_full85_r1}"
export FLEXAID_FORCE_PREP="${FLEXAID_FORCE_PREP:-0}"  # resume: do not force wipe prep
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SOFTBETA_ELECTION=0

CAMPAIGN="$FLEXAID_CAMPAIGN"
BASE="$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine"

count_done() {
  local arm="$1"
  local n=0
  local d="$BASE/$arm/$CAMPAIGN"
  [[ -d "$d" ]] || { echo 0; return; }
  for f in "$d"/*/result.csv; do
    [[ -f "$f" ]] && n=$((n + 1))
  done
  echo "$n"
}

echo "=== full85 resume status (campaign=$CAMPAIGN) ==="
for a in A B0 B; do
  echo "  $a: result.csv count=$(count_done "$a") / 85  OUT=$BASE/$a/$CAMPAIGN"
done
if [[ -f "$FLEXAIDDS_LOCAL_LOGDIR/run_3dsig_red_pair_full85.pid" ]]; then
  pid=$(cat "$FLEXAIDDS_LOCAL_LOGDIR/run_3dsig_red_pair_full85.pid" 2>/dev/null || true)
  if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
    echo "  chain LIVE pid=$pid"
  else
    echo "  chain pid file stale/absent"
  fi
fi
pgrep -lf 'FlexAID --legacy' 2>/dev/null | head -3 || echo "  no FlexAID --legacy"

if (( STATUS )); then
  exit 0
fi

# Resume = re-invoke full85 launcher (skip completed via result.csv)
export FLEXAID_RESTARTS
export FLEXAID_CAMPAIGN
export FLEXAID_FORCE_PREP
args=()
[[ -n "$ONLY" ]] && args+=(--only "$ONLY") || args+=(--from "$FROM")
echo "Resuming full85: ${args[*]} R=$FLEXAID_RESTARTS FORCE_PREP=$FLEXAID_FORCE_PREP"
exec bash "$ROOT/scripts/run_3dsig_red_pair_full85.sh" "${args[@]}"
