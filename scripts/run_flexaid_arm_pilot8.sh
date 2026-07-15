#!/usr/bin/env bash
# Run FlexAID arm A / B0 / B (pilot8 or full85 Astex Diverse).
#
# Default storage (anti-hang): local APFS via use_local_first_benchmark_storage.sh
#   OUT:  $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/...
#   work: $FLEXAIDDS_LOCAL_QUEUE/work/{A,B0,B}/...
# Sync later: bash scripts/sync_three_engine_local_to_icloud.sh
#
# Force live iCloud OUT (legacy): FLEXAIDDS_FORCE_ICLOUD_OUT=1
#
# Usage:
#   scripts/run_flexaid_arm_pilot8.sh A            # real pilot8
#   scripts/run_flexaid_arm_pilot8.sh B0 --dry-run
#   scripts/run_flexaid_arm_pilot8.sh B --pdb 1GPK # single target
#   scripts/run_flexaid_arm_pilot8.sh A --smoke    # tiny pop/gen smoke (1 target)
#   scripts/run_flexaid_arm_pilot8.sh A --full85   # Astex Diverse N=85 (3Dsig full)
#
# Does NOT touch C0 full85 outputs. Separate lock per arm/panel (pilot8 vs full85).
set -euo pipefail

# shellcheck disable=SC1090
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"

_ROOT_FOR_STORAGE="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
if [[ "${FLEXAIDDS_FORCE_ICLOUD_OUT:-0}" == "1" ]]; then
  source "$_ROOT_FOR_STORAGE/scripts/use_icloud_benchmark_storage.sh"
else
  source "$_ROOT_FOR_STORAGE/scripts/use_local_first_benchmark_storage.sh"
fi
# shellcheck disable=SC1091
source "$_ROOT_FOR_STORAGE/scripts/require_icloud_out.sh"

Q="${QUEUE_ROOT:-${FLEXAIDDS_QUEUE_ROOT}}"
export QUEUE_ROOT="$Q"
export FLEXAIDDS_QUEUE_ROOT="$Q"

ARM="${1:-}"
if [[ -z "$ARM" || "$ARM" == -* ]]; then
  echo "Usage: $0 <A|B0|B> [--dry-run|--smoke|--full85|--pdb ID|--force|--no-prepare]" >&2
  exit 2
fi
shift || true

DRY=0
SMOKE=0
FULL85=0
FORCE=0
NOPREP=0
ONLY_PDB=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --smoke) SMOKE=1; shift ;;
    --full85) FULL85=1; shift ;;
    --force) FORCE=1; shift ;;
    --no-prepare) NOPREP=1; shift ;;
    --pdb) ONLY_PDB="${2:-}"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if (( SMOKE && FULL85 )); then
  echo "ERROR: --smoke and --full85 are mutually exclusive" >&2
  exit 2
fi

case "$ARM" in
  A|B0|B) ;;
  *) echo "ERROR: arm must be A, B0, or B (got $ARM)" >&2; exit 2 ;;
esac

BIN_ARM="B"
[[ "$ARM" == "A" ]] && BIN_ARM="A"
BINARY="$Q/bin/$BIN_ARM/FlexAID"
MATRIX="$Q/data/MC_st0r5.2_6.dat"
MATRIX_PIN="72d7c7396702331d96ff12d18f831796"

REPO="${FLEXAIDDS_ROOT:-}"
if [[ -z "$REPO" ]]; then
  if git -C "${HOME}/Projects/FlexAIDdS" rev-parse --show-toplevel >/dev/null 2>&1; then
    REPO="$(git -C "${HOME}/Projects/FlexAIDdS" rev-parse --show-toplevel)"
  else
    REPO="$(cd "$(dirname "$0")/.." && pwd)"
  fi
fi
export FLEXAIDDS_ROOT="$REPO"
GEN="$REPO/scripts/generate_flexaid_inp.py"
PARSE="$REPO/scripts/parse_flexaid_arm_results.py"

PANEL="pilot8"
(( FULL85 )) && PANEL="full85"
(( SMOKE )) && PANEL="smoke"

OUT="${FLEXAID_ARM_OUT:-}"
if [[ -z "$OUT" ]]; then
  if (( SMOKE )); then
    OUT="$FLEXAIDDS_RESULTS/campaigns/three_engine/$ARM/smoke"
  elif (( FULL85 )); then
    OUT="$FLEXAIDDS_RESULTS/campaigns/three_engine/$ARM/3dsig_full85_r10"
  else
    OUT="$FLEXAIDDS_RESULTS/campaigns/three_engine/$ARM/pilot8"
  fi
fi
# Smoke must NEVER share pilot8/full85 OUT (skip would keep tiny-GA result.csv).
if (( SMOKE )) && [[ -z "${FLEXAID_ARM_OUT:-}" ]]; then
  OUT="$FLEXAIDDS_RESULTS/campaigns/three_engine/$ARM/smoke"
fi

WORK_ROOT="${FLEXAID_WORK_ROOT:-$Q/work}"
# Logs on local when local-first; else queue logs (may be iCloud)
LOGDIR="${FLEXAIDDS_LOCAL_LOGDIR:-$Q/logs}"
mkdir -p "$OUT" "$LOGDIR" "$WORK_ROOT" "$FLEXAIDDS_RESULTS/campaigns/three_engine/$ARM"

require_icloud_out "$OUT" || exit 91

# Prefer local binary if queue is on CloudDocs but local staging has the arm binary
_LOCAL_BIN="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/three_engine_entropy_q1/bin/$BIN_ARM/FlexAID"
if [[ ! -x "$BINARY" && -x "$_LOCAL_BIN" ]]; then
  BINARY="$_LOCAL_BIN"
elif [[ "$BINARY" == *"/Mobile Documents/"* && -x "$_LOCAL_BIN" ]]; then
  BINARY="$_LOCAL_BIN"
  echo "OK: using local FlexAID binary (avoid CloudDocs): $BINARY"
fi
_LOCAL_MAT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/three_engine_entropy_q1/data/MC_st0r5.2_6.dat"
if [[ -f "$_LOCAL_MAT" ]]; then
  MATRIX="$_LOCAL_MAT"
fi

echo "=== FlexAID arm $ARM $PANEL ==="
echo "Q=$Q"
echo "OUT=$OUT"
echo "BINARY=$BINARY"
echo "WORK_ROOT=$WORK_ROOT"
echo "LOGDIR=$LOGDIR"
echo "ALLOW_LOCAL_OUT=${FLEXAIDDS_ALLOW_LOCAL_OUT:-0}"

[[ -x "$BINARY" ]] || { echo "FAIL: not executable: $BINARY" >&2; exit 1; }
[[ -f "$MATRIX" ]] || { echo "FAIL: matrix missing: $MATRIX" >&2; exit 1; }
GOT=$(md5 -q "$MATRIX")
[[ "$GOT" == "$MATRIX_PIN" ]] || { echo "FAIL: matrix MD5 $GOT != $MATRIX_PIN" >&2; exit 90; }
echo "OK: matrix md5=$GOT"
[[ -f "$GEN" ]] || { echo "FAIL: missing $GEN" >&2; exit 1; }

if [[ -z "${FLEXAIDDS_PROCESSLIGAND:-}" ]]; then
  for m in "$REPO"/.venv-processligand/lib/python*/site-packages/processligandpy/bin/ProcessLigand "$Q/bin/ProcessLigand"; do
    if [[ -x "$m" ]]; then
      export FLEXAIDDS_PROCESSLIGAND="$m"
      break
    fi
  done
fi
if [[ -z "${FLEXAIDDS_PROCESSLIGAND:-}" ]]; then
  echo "FAIL: ProcessLigand not found. Install processligand-py or set FLEXAIDDS_PROCESSLIGAND" >&2
  exit 1
fi
echo "OK: ProcessLigand=$FLEXAIDDS_PROCESSLIGAND"

LOCK="$LOGDIR/run_${ARM}_${PANEL}.lock"
if [[ -f "$LOCK" ]] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
  echo "REFUSE: arm $ARM $PANEL already running pid $(cat "$LOCK")" >&2
  exit 92
fi
if [[ -f "$LOGDIR/C0_full85.lock" ]] && kill -0 "$(cat "$LOGDIR/C0_full85.lock" 2>/dev/null)" 2>/dev/null; then
  echo "NOTE: C0 full85 is live (pid $(cat "$LOGDIR/C0_full85.lock")) — proceeding for arm $ARM only (separate OUT)."
fi

PILOT8=(1G9V 1GPK 1MEH 1P62 1Q4G 1R9O 1T40 2BYS)
if [[ -n "$ONLY_PDB" ]]; then
  TARGETS=("${ONLY_PDB^^}")
elif (( SMOKE )); then
  TARGETS=(1GPK)
elif (( FULL85 )); then
  # All prepared Astex Diverse inputs under queue (must be 85 for deck claim)
  TARGETS=()
  while IFS= read -r _tid; do
    [[ -n "$_tid" ]] && TARGETS+=("$_tid")
  done < <(
    # Local queue only — never walk CloudDocs. basename dirs under inputs/astex_diverse.
    for _d in "$Q/inputs/astex_diverse"/*/; do
      [[ -d "$_d" ]] || continue
      basename "${_d%/}"
    done | LC_ALL=C sort
  )
  if (( ${#TARGETS[@]} == 0 )); then
    echo "FAIL: no inputs under $Q/inputs/astex_diverse" >&2
    exit 1
  fi
  if (( ${#TARGETS[@]} != 85 )); then
    echo "WARN: expected 85 Astex Diverse targets, found ${#TARGETS[@]} under $Q/inputs/astex_diverse" >&2
  fi
else
  TARGETS=("${PILOT8[@]}")
fi

POP="${FLEXAID_POP:-1000}"
GENS="${FLEXAID_GEN:-2000}"
RESTARTS="${FLEXAID_RESTARTS:-5}"
if (( SMOKE )); then
  POP="${FLEXAID_SMOKE_POP:-20}"
  GENS="${FLEXAID_SMOKE_GEN:-5}"
  RESTARTS="${FLEXAID_SMOKE_RESTARTS:-1}"
  echo "SMOKE mode: pop=$POP gen=$GENS R=$RESTARTS targets=${TARGETS[*]}"
fi

echo "targets_n=${#TARGETS[@]} pop=$POP gen=$GENS R=$RESTARTS"
if (( ${#TARGETS[@]} <= 12 )); then
  echo "targets=${TARGETS[*]}"
else
  _last=$((${#TARGETS[@]} - 1))
  echo "targets=${TARGETS[0]} ... ${TARGETS[$_last]} (n=${#TARGETS[@]})"
fi

if (( DRY )); then
  echo "DRY-RUN: prepare plan only (no FlexAID exec)"
  PREP_ARGS=(--queue-root "$Q" --work-root "$WORK_ROOT" --arms "$ARM" --pop "$POP" --gen "$GENS" --restarts "$RESTARTS" --dry-run)
  if [[ -n "$ONLY_PDB" ]]; then
    PREP_ARGS+=(--pdb "$ONLY_PDB")
  elif (( SMOKE )); then
    PREP_ARGS+=(--pdb 1GPK)
  elif (( FULL85 )); then
    : # all targets under inputs/astex_diverse (no --pilot8)
  else
    PREP_ARGS+=(--pilot8)
  fi
  python3 "$GEN" "${PREP_ARGS[@]}"
  echo "DRY-RUN complete (no launch)"
  exit 0
fi

if (( NOPREP == 0 )); then
  PREP_ARGS=(--queue-root "$Q" --work-root "$WORK_ROOT" --arms "$ARM" --pop "$POP" --gen "$GENS" --restarts "$RESTARTS")
  (( FORCE )) && PREP_ARGS+=(--force)
  if [[ ${#TARGETS[@]} -eq 1 ]]; then
    PREP_ARGS+=(--pdb "${TARGETS[0]}")
  elif (( FULL85 )); then
    : # generate for full queue inputs/astex_diverse
  else
    PREP_ARGS+=(--pilot8)
  fi
  python3 "$GEN" "${PREP_ARGS[@]}"
fi

python3 - <<PY
import hashlib, json, time
from pathlib import Path
out = Path(r'''$OUT''')
q = Path(r'''$Q''')
binp = Path(r'''$BINARY''').resolve()
mat = Path(r'''$MATRIX''')
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()
receipt = {
    "run_id": f"flexaid_{'''$ARM'''}_{'''$PANEL'''}",
    "arm": '''$ARM''',
    "panel": '''$PANEL''',
    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "output": str(out),
    "storage": "local_first_or_icloud",
    "pop": int('''$POP'''),
    "gen": int('''$GENS'''),
    "restarts": int('''$RESTARTS'''),
    "targets": '''${TARGETS[*]}'''.split(),
    "n_targets": int('''${#TARGETS[@]}'''),
    "matrix_md5": hashlib.md5(mat.read_bytes()).hexdigest(),
    "binary_sha256": sha(binp),
    "binary": str(binp),
    "smoke": bool(int('''$SMOKE''')),
    "full85": bool(int('''$FULL85''')),
    "seed_note": "STRTSEED written; staged binaries may use time(0)",
    "cli_mode": "--legacy classic CONFIG+ga",
}
(out / "RUN_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
(q / f"provenance_run_{'''$ARM'''}_{'''$PANEL'''}.json").write_text(json.dumps(receipt, indent=2) + "\n")
print("wrote", out / "RUN_RECEIPT.json")
PY

LOG="$LOGDIR/run_${ARM}_${PANEL}.log"
run_one() {
  local pdb="$1"
  local wdir="$WORK_ROOT/$ARM/$pdb"
  local odir="$OUT/$pdb"
  mkdir -p "$odir"
  if [[ -f "$odir/result.csv" && "$FORCE" -eq 0 ]]; then
    # Refuse skip when existing receipt is smoke-budget (prevents contamination)
    if [[ -f "$odir/RUN_RECEIPT.json" ]] && grep -q '"smoke"[[:space:]]*:[[:space:]]*true' "$odir/RUN_RECEIPT.json" 2>/dev/null; then
      echo "REFUSE SKIP $ARM/$pdb: existing result is smoke-budget — use --force or remove $odir"
      return 1
    fi
    if (( ! SMOKE )) && [[ -f "$odir/meta.json" ]] && grep -q '"smoke"[[:space:]]*:[[:space:]]*true' "$odir/meta.json" 2>/dev/null; then
      echo "REFUSE SKIP $ARM/$pdb: meta marks smoke — use --force or remove $odir"
      return 1
    fi
    echo "SKIP $ARM/$pdb (result.csv exists)"
    return 0
  fi
  if [[ ! -f "$wdir/meta.json" ]]; then
    echo "FAIL: missing prepared work $wdir" >&2
    return 1
  fi

  local r=0
  local n_ok=0
  while (( r < RESTARTS )); do
    local rdir="$wdir/restart_$r"
    local prefix="$odir/${pdb}_r${r}"
    echo "=== $ARM $pdb restart $r ==="
    # Staged A/B binaries are FlexAIDdS-unified CLI; classic 3-file mode needs --legacy
    echo "CMD: $BINARY --legacy $rdir/CONFIG.inp $rdir/ga.inp $prefix"
    if ! (
      cd "$rdir"
      caffeinate -i -s "$BINARY" --legacy "$rdir/CONFIG.inp" "$rdir/ga.inp" "$prefix"
    ); then
      echo "WARN: FlexAID exit non-zero for $ARM $pdb r$r" >&2
    else
      n_ok=$((n_ok + 1))
    fi
    r=$((r + 1))
  done

  if [[ -f "$PARSE" ]]; then
    python3 "$PARSE" --arm "$ARM" --pdb "$pdb" --out-dir "$odir" --work-dir "$wdir" \
      --matrix-md5 "$MATRIX_PIN" --binary "$BINARY" || true
  fi
  echo "DONE $ARM/$pdb restarts_ok=$n_ok/$RESTARTS"
}

if [[ "${NOHUP:-0}" == "1" ]]; then
  {
    echo "$$" >"$LOCK"
    for pdb in "${TARGETS[@]}"; do
      run_one "$pdb" || true
    done
    rm -f "$LOCK"
  } >>"$LOG" 2>&1 &
  echo "STARTED background pid=$! log=$LOG out=$OUT"
else
  echo "$$" >"$LOCK"
  trap 'rm -f "$LOCK"' EXIT
  for pdb in "${TARGETS[@]}"; do
    run_one "$pdb" || true
  done
  echo "COMPLETE arm=$ARM panel=$PANEL out=$OUT log=$LOG"
fi
