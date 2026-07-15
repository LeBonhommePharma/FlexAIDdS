#!/usr/bin/env bash
# DP vs FO clustering pilot — 8 random Astex Diverse cases, SMALL GA budget.
# Isolated OUT namespaces (never writes into C0 claim OUT).
#
# Protocol (small-sim):
#   pop=200  gen=50  T=298  EVAL_SCALE=-1 (fixed budget)  seed OFF
#   arms: FO and DP (matched seeds via engine ga.seed in dock_config if fixed later)
#   mode: defined-cleft-redock
#
# Usage:
#   source ~/.flexaidds_env
#   bash scripts/run_DPFO_pilot8_small.sh [--dry-run] [--nowait] [--arm FO|DP|both]
#
# --nowait: do not wait for live FlexAIDdS claim processes (still uses nice 19)
# default: wait until no FlexAIDdS binary is running before each target (yield to C0)
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_icloud_benchmark_storage.sh" 2>/dev/null || true
# shellcheck disable=SC1091
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"

Q="${FLEXAIDDS_QUEUE_ROOT:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/queues/three_engine_entropy_q1}"
R="${FLEXAIDDS_RESULTS:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/results}"
MANIFEST="${DPFO_MANIFEST:-$Q/inputs/astex_dpfo_pilot8_random.json}"
RUNNER="${DPFO_RUNNER:-$Q/bin/C/benchmark_datasets}"
BINARY="${DPFO_BINARY:-$Q/bin/C/FlexAIDdS}"
LOGDIR="$Q/logs"
STAMP="20260715"
BASE_OUT="${DPFO_OUT_BASE:-$R/campaigns/DPFO_pilot8_small_g50_p200_${STAMP}}"

POP="${DPFO_POP:-200}"
GEN="${DPFO_GEN:-50}"
TEMP="${DPFO_TEMP:-298}"
TIMEOUT="${DPFO_TIMEOUT:-1800}"

DRY=0
NOWAIT=0
ARM="both"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1 ;;
    --nowait) NOWAIT=1 ;;
    --arm) ARM="${2:?}"; shift ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
  shift
done

mkdir -p "$LOGDIR" "$BASE_OUT"
LOCK="$LOGDIR/DPFO_pilot8_small.lock"
PIDF="$LOGDIR/DPFO_pilot8_small.pid"
LOG="$LOGDIR/DPFO_pilot8_small.log"

export FLEXAIDDS_BINARY="$BINARY"
export FLEXAIDDS_DATA_DIR="${FLEXAIDDS_DATA_DIR:-$Q/data}"
export FLEXAIDDS_RESTARTS=2
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_VCT_R0=4
export FLEXAIDDS_SHARING_ALPHA=4
export SHARING_ALPHA=4
# Fixed small budget — do not pop×DoF scale for this pilot
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=-1
export EVAL_SCALE_DIHEDRAL=-1
export FLEXAIDDS_BUDGET_SCALE=1
export FLEXAIDDS_NATIVE_SEED_FRAC=0
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_POSEBUSTERS_BIN="${FLEXAIDDS_POSEBUSTERS_BIN:-$ROOT/.venv-posebusters/bin/bust}"
unset FLEXAIDDS_FORCE_SEED 2>/dev/null || true
unset FLEXAIDDS_USE_DP 2>/dev/null || true

echo "=== DP/FO SMALL-SIM PILOT8 ===" | tee -a "$LOG"
echo "MANIFEST=$MANIFEST" | tee -a "$LOG"
echo "BASE_OUT=$BASE_OUT" | tee -a "$LOG"
echo "pop=$POP gen=$GEN T=$TEMP EVAL_SCALE=-1 ARM=$ARM NOWAIT=$NOWAIT" | tee -a "$LOG"
echo "RUNNER=$RUNNER BINARY=$BINARY" | tee -a "$LOG"

[[ -x "$RUNNER" ]] || { echo "FAIL: runner not executable: $RUNNER" >&2; exit 1; }
[[ -x "$BINARY" ]] || { echo "FAIL: binary not executable: $BINARY" >&2; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "FAIL: manifest missing: $MANIFEST" >&2; exit 1; }

if [[ -f "$LOCK" ]] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
  echo "REFUSE: pilot already running pid $(cat "$LOCK")" >&2
  exit 92
fi

# Refuse if OUT collides with claim
if [[ "$BASE_OUT" == *"C0_full85_claim"* ]]; then
  echo "REFUSE: pilot OUT must not be claim OUT" >&2
  exit 91
fi

if (( DRY )); then
  echo "DRY-RUN OK"
  exit 0
fi

echo $$ >"$LOCK"
echo $$ >"$PIDF"

wait_for_slot() {
  if (( NOWAIT )); then
    return 0
  fi
  local n=0
  while pgrep -f '/bin/C/FlexAIDdS|bin/C/FlexAID' >/dev/null 2>&1; do
    n=$((n + 1))
    if (( n % 12 == 1 )); then
      echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] yield: claim FlexAIDdS still running; waiting 30s…" | tee -a "$LOG"
    fi
    sleep 30
  done
  # brief settle so claim can respawn next target
  sleep 5
  if pgrep -f '/bin/C/FlexAIDdS|bin/C/FlexAID' >/dev/null 2>&1; then
    wait_for_slot
  fi
}

run_arm() {
  local cluster="$1"
  local out="$BASE_OUT/${cluster}"
  mkdir -p "$out"
  local arm_log="$LOGDIR/DPFO_pilot8_${cluster}.log"
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] START arm=$cluster OUT=$out" | tee -a "$LOG" "$arm_log"

  # FO/DP via CLI --clustering (not FLEXAIDDS_USE_DP alone; FO needs explicit FO)
  unset FLEXAIDDS_USE_DP 2>/dev/null || true
  if [[ "$cluster" == "DP" ]]; then
    export FLEXAIDDS_USE_DP=1
  fi

  wait_for_slot

  # shellcheck disable=SC2086
  nice -n 19 env \
    FLEXAIDDS_BINARY="$BINARY" \
    FLEXAIDDS_DATA_DIR="$FLEXAIDDS_DATA_DIR" \
    FLEXAIDDS_EVAL_SCALE_DIHEDRAL=-1 \
    EVAL_SCALE_DIHEDRAL=-1 \
    FLEXAIDDS_BUDGET_SCALE=1 \
    FLEXAIDDS_NATIVE_SEED_FRAC=0 \
    FLEXAIDDS_SEED_ELITISM=0 \
    FLEXAIDDS_RESTARTS=2 \
    FLEXAIDDS_USE_DP="${FLEXAIDDS_USE_DP:-}" \
    "$RUNNER" \
      --benchmark "crossdock_json:$MANIFEST" \
      --mode defined-cleft-redock \
      --output "$out/" \
      --threads 1 \
      --omp-threads 1 \
      --ga-population "$POP" \
      --ga-generations "$GEN" \
      --temperature "$TEMP" \
      --clustering "$cluster" \
      --job-timeout-seconds "$TIMEOUT" \
      --force \
    >>"$arm_log" 2>&1

  local rc=$?
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] END arm=$cluster rc=$rc" | tee -a "$LOG" "$arm_log"
  return $rc
}

# Provenance
python3 - <<PY
import hashlib, json, time, os
from pathlib import Path
out = Path(r"""$BASE_OUT""")
out.mkdir(parents=True, exist_ok=True)
manifest = Path(r"""$MANIFEST""")
runner = Path(r"""$RUNNER""")
binary = Path(r"""$BINARY""")

def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

ids = [p["receptor_id"] for p in json.loads(manifest.read_text())["pairs"]]
receipt = {
    "campaign": "DPFO_pilot8_small_g50_p200_${STAMP}",
    "purpose": "DP vs FO clustering success rate + predictive power (small sim)",
    "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "manifest": str(manifest),
    "manifest_sha256": sha256(manifest),
    "pdb_ids": ids,
    "seed_sample": 20260715,
    "protocol": {
        "pop": int("$POP"),
        "gen": int("$GEN"),
        "temperature": int("$TEMP"),
        "EVAL_SCALE_DIHEDRAL": -1,
        "mode": "defined-cleft-redock",
        "arms": ["FO", "DP"],
        "metrics": ["S1_elected_rmsd<=2", "BCR_best_cluster_rmsd<=2", "election_gap", "predictive_power"],
    },
    "runner_sha256": sha256(runner),
    "binary_sha256": sha256(binary),
    "note": "Does not touch C0 claim OUT. Yields to live FlexAIDdS unless --nowait.",
}
(out / "PROVENANCE.json").write_text(json.dumps(receipt, indent=2) + "\n")
print("wrote", out / "PROVENANCE.json")
print("targets", ids)
PY

rc_all=0
case "$ARM" in
  FO|fo) run_arm FO || rc_all=$? ;;
  DP|dp) run_arm DP || rc_all=$? ;;
  both)
    # FO first (default CF→FO path), then DP
    run_arm FO || rc_all=$?
    run_arm DP || rc_all=$?
    ;;
  *) echo "bad --arm $ARM" >&2; exit 2 ;;
esac

# Aggregate when both arms have any result.csv
python3 "$ROOT/scripts/aggregate_dpfo_pilot.py" --base-out "$BASE_OUT" 2>>"$LOG" || true

rm -f "$LOCK"
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] PILOT DONE rc=$rc_all" | tee -a "$LOG"
exit $rc_all
