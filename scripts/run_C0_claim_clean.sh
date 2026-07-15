#!/usr/bin/env bash
# Clean claim C0 full85 — gen=2000, pop base 1000 + EVAL_SCALE=1 (pop×DoF), seed OFF.
# New iCloud OUT only (never dual-launch dirty residual OUT).
#
# Usage:
#   source ~/.flexaidds_env
#   source scripts/use_icloud_benchmark_storage.sh
#   bash scripts/run_C0_claim_clean.sh [--dry-run]
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_icloud_benchmark_storage.sh"
# shellcheck disable=SC1091
source "$ROOT/scripts/require_icloud_out.sh"
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"

Q="${FLEXAIDDS_QUEUE_ROOT}"
R="${FLEXAIDDS_RESULTS}"
CAMPAIGN_ID="C0_full85_claim_g2000_popmod_20260715"
OUT="${C0_CLAIM_OUT:-$R/campaigns/$CAMPAIGN_ID}"
MANIFEST="${C0_CLAIM_MANIFEST:-$Q/inputs/astex_native_85.json}"
LOGDIR="$Q/logs"
LOCK="$LOGDIR/C0_claim_clean.lock"
PIDF="$LOGDIR/C0_claim_clean.pid"
LOG="$LOGDIR/C0_claim_clean.log"
RUNNER="${C0_CLAIM_RUNNER:-$Q/bin/C/benchmark_datasets}"
BINARY="${C0_CLAIM_BINARY:-$Q/bin/C/FlexAIDdS}"

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

require_icloud_out "$OUT" || exit 91
mkdir -p "$OUT" "$LOGDIR"

# Claim knobs (publish freeze)
export FLEXAIDDS_BINARY="$BINARY"
export FLEXAIDDS_DATA_DIR="${FLEXAIDDS_DATA_DIR:-$Q/data}"
export FLEXAIDDS_RESTARTS=5
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_VCT_R0=4
export FLEXAIDDS_SHARING_ALPHA=4
export SHARING_ALPHA=4
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
export EVAL_SCALE_DIHEDRAL=1
export FLEXAIDDS_BUDGET_SCALE=1
export FLEXAIDDS_NATIVE_SEED_FRAC=0
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_POSEBUSTERS_BIN="${FLEXAIDDS_POSEBUSTERS_BIN:-$ROOT/.venv-posebusters/bin/bust}"
unset FLEXAIDDS_FORCE_SEED 2>/dev/null || true

POP=1000
GEN=2000
TEMP=298

echo "=== C0 CLEAN CLAIM LAUNCH ==="
echo "OUT=$OUT"
echo "MANIFEST=$MANIFEST"
echo "RUNNER=$RUNNER"
echo "BINARY=$BINARY"
echo "pop=$POP gen=$GEN T=$TEMP R=5 EVAL_SCALE=1 BUDGET_SCALE=1 seed=OFF"
echo "max_results (engine)=50 (kBenchmarkPoseLimit)"

[[ -x "$RUNNER" ]] || { echo "FAIL: runner not executable: $RUNNER" >&2; exit 1; }
[[ -x "$BINARY" ]] || { echo "FAIL: binary not executable: $BINARY" >&2; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "FAIL: manifest missing: $MANIFEST" >&2; exit 1; }

# Dual-launch guard: refuse if this OUT or claim lock is live
if [[ -f "$LOCK" ]] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
  echo "REFUSE: claim C0 already running pid $(cat "$LOCK")" >&2
  exit 92
fi
if [[ -f "$Q/logs/C0_full85.pid" ]] && kill -0 "$(cat "$Q/logs/C0_full85.pid" 2>/dev/null)" 2>/dev/null; then
  echo "REFUSE: legacy C0 residual still live — stop it first" >&2
  exit 92
fi

if (( DRY )); then
  echo "DRY-RUN OK"
  exit 0
fi

# Provenance receipt
python3 - <<PY
import hashlib, json, time, os
from pathlib import Path

def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

out = Path(r"""$OUT""")
q = Path(r"""$Q""")
runner = Path(r"""$RUNNER""")
binary = Path(r"""$BINARY""")
mat = q / "data" / "MC_st0r5.2_6.dat"
receipt = {
    "run_id": r"""$CAMPAIGN_ID""",
    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "output": str(out),
    "storage": "icloud_drive",
    "mode": "defined-cleft-redock",
    "temperature_K": $TEMP,
    "pop": $POP,
    "gen": $GEN,
    "restarts": 5,
    "native_seed_frac": 0,
    "seed_elitism": 0,
    "eval_scale_dihedral": 1,
    "budget_scale": 1,
    "max_results": 50,
    "matrix_md5": hashlib.md5(mat.read_bytes()).hexdigest() if mat.is_file() else None,
    "binary_sha256": sha256(binary) if binary.is_file() else None,
    "runner_sha256": sha256(runner) if runner.is_file() else None,
    "manifest": r"""$MANIFEST""",
    "claim_path": True,
    "supersedes": "C0_full85_defined_cleft_nativeseed_forbidden (dirty EVAL_SCALE=-1 residual)",
    "env": {
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL": "1",
        "FLEXAIDDS_BUDGET_SCALE": "1",
        "FLEXAIDDS_NATIVE_SEED_FRAC": "0",
        "FLEXAIDDS_SEED_ELITISM": "0",
        "FLEXAIDDS_PARALLEL_RESTARTS": "0",
        "OMP_NUM_THREADS": "1",
    },
}
(out / "RUN_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
(q / "provenance_run_C0_claim_clean.json").write_text(json.dumps(receipt, indent=2) + "\n")
print("wrote", out / "RUN_RECEIPT.json")
PY

export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0

nohup caffeinate -i -s "$RUNNER" \
  --benchmark "crossdock_json:$MANIFEST" \
  --mode defined-cleft-redock \
  --output "$OUT/" \
  --threads 1 \
  --omp-threads 1 \
  --ga-population "$POP" \
  --ga-generations "$GEN" \
  --temperature "$TEMP" \
  --job-timeout-seconds 10800 \
  >>"$LOG" 2>&1 &
echo $! | tee "$LOCK" > "$PIDF"
# also point legacy names at claim run so monitors find it
echo "$(cat "$PIDF")" > "$Q/logs/C0_full85.pid"
echo "$(cat "$PIDF")" > "$Q/logs/C0_full85.lock"

echo "STARTED claim C0 pid=$(cat "$PIDF") log=$LOG out=$OUT"
sleep 2
kill -0 "$(cat "$PIDF")" && echo LIVE || echo DEAD
tail -5 "$LOG" || true
