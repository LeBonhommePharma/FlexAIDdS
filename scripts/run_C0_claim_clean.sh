#!/usr/bin/env bash
# Clean claim C0 full85 — gen=2000, pop base 1000 + EVAL_SCALE=1 (pop×DoF), seed OFF.
#
# DEFAULT (anti-hang): live GA OUT + logs + binaries on local APFS outside CloudDocs.
# iCloud holds the durable mirror (result.csv + thin metadata) via claim_icloud_sync_loop.sh.
#
# Usage:
#   bash scripts/run_C0_claim_clean.sh [--dry-run] [--icloud-out]
#
#   --icloud-out   legacy: write OUT directly under iCloud (NOT recommended; hangs)
#
# Env overrides: C0_CLAIM_OUT, C0_CLAIM_RUNNER, C0_CLAIM_BINARY, FLEXAIDDS_CLAIM_LOCAL=0
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/claim_local_staging_paths.sh"
# shellcheck disable=SC1091
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"

DRY=0
FORCE_ICLOUD_OUT=0
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --icloud-out) FORCE_ICLOUD_OUT=1 ;;
  esac
done

# Local staging ON by default (fix CloudDocs hang during GA + re-election).
USE_LOCAL=1
if [[ "${FLEXAIDDS_CLAIM_LOCAL:-1}" == "0" ]] || (( FORCE_ICLOUD_OUT )); then
  USE_LOCAL=0
fi

CAMPAIGN_ID="${C0_CAMPAIGN_ID}"
Q="${FLEXAIDDS_QUEUE_ROOT}"
MANIFEST="${C0_CLAIM_MANIFEST}"

if (( USE_LOCAL )); then
  # Ensure local-first layout (campaigns/, logs/, pins/, three_engine_entropy_q1/…)
  bash "$ROOT/scripts/ensure_local_first_layout.sh" || {
    echo "FAIL: ensure_local_first_layout.sh" >&2
    exit 1
  }
  OUT="${C0_CLAIM_OUT:-$C0_CLAIM_LOCAL_OUT}"
  LOGDIR="${C0_CLAIM_LOCAL_LOGDIR}"
  RUNNER="${C0_CLAIM_RUNNER}"
  BINARY="${C0_CLAIM_BINARY}"
  DATA_DIR="${C0_CLAIM_DATA_DIR}"
  ICLOUD_MIRROR="${C0_CLAIM_ICLOUD_OUT}"
  STORAGE_LABEL="local_apfs_mirror_icloud"
else
  # shellcheck disable=SC1091
  source "$ROOT/scripts/use_icloud_benchmark_storage.sh"
  # shellcheck disable=SC1091
  source "$ROOT/scripts/require_icloud_out.sh"
  OUT="${C0_CLAIM_OUT:-$FLEXAIDDS_RESULTS/campaigns/$CAMPAIGN_ID}"
  LOGDIR="${Q}/logs"
  RUNNER="${C0_CLAIM_RUNNER:-$Q/bin/C/benchmark_datasets}"
  BINARY="${C0_CLAIM_BINARY:-$Q/bin/C/FlexAIDdS}"
  DATA_DIR="${FLEXAIDDS_DATA_DIR:-$Q/data}"
  ICLOUD_MIRROR=""
  STORAGE_LABEL="icloud_drive_direct"
  require_icloud_out "$OUT" || exit 91
fi

LOCK="$LOGDIR/C0_claim_clean.lock"
PIDF="$LOGDIR/C0_claim_clean.pid"
LOG="$LOGDIR/C0_claim_clean.log"

mkdir -p "$OUT" "$LOGDIR"

# Prefer local data if present when local mode
if (( USE_LOCAL )) && [[ ! -d "$DATA_DIR" ]]; then
  DATA_DIR="${Q}/data"
fi

# Claim knobs (publish freeze)
export FLEXAIDDS_BINARY="$BINARY"
export FLEXAIDDS_DATA_DIR="$DATA_DIR"
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
# 3Dsig Shannon ranking path (G̃=H̃−T·S̃) — engine default is OFF until Astex pilot;
# this claim launch explicitly opts in. LEGACY_ZH=1 rolls back to legacy ZH / min-CF.
export FLEXAIDDS_ELECTION_SHANNON_F="${FLEXAIDDS_ELECTION_SHANNON_F:-1}"
export FLEXAIDDS_ELECTION_SOFT_T="${FLEXAIDDS_ELECTION_SOFT_T:-298}"
unset FLEXAIDDS_FORCE_SEED 2>/dev/null || true
unset FLEXAIDDS_ELECTION_LEGACY_ZH 2>/dev/null || true

POP=1000
GEN=2000
TEMP=298

echo "=== C0 CLEAN CLAIM LAUNCH ==="
echo "storage=$STORAGE_LABEL USE_LOCAL=$USE_LOCAL"
echo "OUT=$OUT"
[[ -n "$ICLOUD_MIRROR" ]] && echo "ICLOUD_MIRROR=$ICLOUD_MIRROR"
echo "LOGDIR=$LOGDIR"
echo "MANIFEST=$MANIFEST"
echo "RUNNER=$RUNNER"
echo "BINARY=$BINARY"
echo "DATA_DIR=$DATA_DIR"
echo "pop=$POP gen=$GEN T=$TEMP R=5 EVAL_SCALE=1 BUDGET_SCALE=1 seed=OFF"
echo "election: SHANNON_F=${FLEXAIDDS_ELECTION_SHANNON_F} SOFT_T=${FLEXAIDDS_ELECTION_SOFT_T}"
echo "max_results (engine)=50 (kBenchmarkPoseLimit)"

[[ -x "$RUNNER" ]] || { echo "FAIL: runner not executable: $RUNNER" >&2; exit 1; }
[[ -x "$BINARY" ]] || { echo "FAIL: binary not executable: $BINARY" >&2; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "FAIL: manifest missing: $MANIFEST" >&2; exit 1; }
[[ -d "$DATA_DIR" ]] || { echo "FAIL: data dir missing: $DATA_DIR" >&2; exit 1; }

# Dual-launch guard: refuse if this OUT lock is live (local or iCloud pid files)
is_live_pidfile() {
  local f="$1" p
  [[ -f "$f" ]] || return 1
  p=$(tr -d ' \n' <"$f" 2>/dev/null || true)
  [[ -n "$p" ]] && kill -0 "$p" 2>/dev/null
}

if is_live_pidfile "$LOCK"; then
  echo "REFUSE: claim C0 already running pid $(cat "$LOCK")" >&2
  exit 92
fi
# Also refuse if another claim worker still holds the old iCloud lock
if is_live_pidfile "${Q}/logs/C0_claim_clean.pid" || is_live_pidfile "${Q}/logs/C0_full85.pid"; then
  echo "REFUSE: another claim residual still live (iCloud pid file)" >&2
  exit 92
fi
if is_live_pidfile "${C0_CLAIM_LOCAL_LOGDIR}/C0_claim_clean.pid"; then
  if [[ "$(cat "${C0_CLAIM_LOCAL_LOGDIR}/C0_claim_clean.pid" 2>/dev/null)" != "$(cat "$LOCK" 2>/dev/null || true)" ]]; then
    # only if different lock path still live
    if [[ "$LOCK" != "${C0_CLAIM_LOCAL_LOGDIR}/C0_claim_clean.lock" ]] && \
       is_live_pidfile "${C0_CLAIM_LOCAL_LOGDIR}/C0_claim_clean.pid"; then
      echo "REFUSE: local claim already running" >&2
      exit 92
    fi
  fi
fi

if (( DRY )); then
  echo "DRY-RUN OK"
  exit 0
fi

# Clear stale (dead) locks so we can relaunch
rm -f "$LOCK" "$PIDF" 2>/dev/null || true
rm -f "${Q}/logs/C0_claim_clean.lock" "${Q}/logs/C0_claim_clean.pid" \
      "${Q}/logs/C0_full85.lock" "${Q}/logs/C0_full85.pid" 2>/dev/null || true

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
runner = Path(r"""$RUNNER""")
binary = Path(r"""$BINARY""")
data = Path(r"""$DATA_DIR""")
mat = data / "MC_st0r5.2_6.dat"
if not mat.is_file():
    mat = Path(r"""$Q""") / "data" / "MC_st0r5.2_6.dat"
receipt = {
    "run_id": r"""$CAMPAIGN_ID""",
    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "output": str(out),
    "storage": r"""$STORAGE_LABEL""",
    "icloud_mirror": r"""$ICLOUD_MIRROR""" or None,
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
    "anti_hang": bool(int("""$USE_LOCAL""")),
    "supersedes": "C0_full85_defined_cleft_nativeseed_forbidden (dirty EVAL_SCALE=-1 residual)",
    "env": {
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL": "1",
        "FLEXAIDDS_BUDGET_SCALE": "1",
        "FLEXAIDDS_NATIVE_SEED_FRAC": "0",
        "FLEXAIDDS_SEED_ELITISM": "0",
        "FLEXAIDDS_PARALLEL_RESTARTS": "0",
        "FLEXAIDDS_ELECTION_SHANNON_F": os.environ.get("FLEXAIDDS_ELECTION_SHANNON_F", "1"),
        "FLEXAIDDS_ELECTION_SOFT_T": os.environ.get("FLEXAIDDS_ELECTION_SOFT_T", "298"),
        "OMP_NUM_THREADS": "1",
    },
}
(out / "RUN_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
# thin local provenance only (never block on iCloud write)
logdir = Path(r"""$LOGDIR""")
logdir.mkdir(parents=True, exist_ok=True)
(logdir / "provenance_run_C0_claim_clean.json").write_text(json.dumps(receipt, indent=2) + "\n")
print("wrote", out / "RUN_RECEIPT.json")
PY

export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0

# Truncate log for this launch (keep prior as .prev if large)
if [[ -f "$LOG" ]] && [[ $(stat -f%z "$LOG" 2>/dev/null || echo 0) -gt 100000 ]]; then
  mv "$LOG" "${LOG}.prev.$(date +%Y%m%d_%H%M%S)" 2>/dev/null || true
fi

nohup caffeinate -i -s env \
  FLEXAIDDS_BINARY="$BINARY" \
  FLEXAIDDS_DATA_DIR="$DATA_DIR" \
  FLEXAIDDS_RESTARTS=5 \
  FLEXAIDDS_PARALLEL_RESTARTS=0 \
  FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1 \
  EVAL_SCALE_DIHEDRAL=1 \
  FLEXAIDDS_BUDGET_SCALE=1 \
  FLEXAIDDS_NATIVE_SEED_FRAC=0 \
  FLEXAIDDS_SEED_ELITISM=0 \
  FLEXAIDDS_ELECTION_SHANNON_F="${FLEXAIDDS_ELECTION_SHANNON_F}" \
  FLEXAIDDS_ELECTION_SOFT_T="${FLEXAIDDS_ELECTION_SOFT_T}" \
  OMP_NUM_THREADS=1 \
  OMP_WAIT_POLICY=passive \
  "$RUNNER" \
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

# Mirror pid to iCloud queue logs only if quick (non-blocking best-effort)
{
  echo "$(cat "$PIDF")" > "${Q}/logs/C0_claim_clean.pid" 2>/dev/null || true
  echo "$(cat "$PIDF")" > "${Q}/logs/C0_full85.pid" 2>/dev/null || true
  cp "$PIDF" "${Q}/logs/C0_claim_clean.lock" 2>/dev/null || true
} &

echo "STARTED claim C0 pid=$(cat "$PIDF") log=$LOG out=$OUT"
sleep 3
if kill -0 "$(cat "$PIDF")" 2>/dev/null; then
  echo LIVE
else
  echo DEAD
fi
tail -8 "$LOG" 2>/dev/null || true
