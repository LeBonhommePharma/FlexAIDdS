#!/usr/bin/env bash
# Serial 3Dsig red-pair FULL Astex Diverse (N=85): A then B0 then B (LOCAL-FIRST)
#
# PREPARE-ONLY while pilot8 is live: this launcher REFUSES to start if the
# pilot8 chain or any FlexAID --legacy dock is running.
#
# Deck knobs (same as pilot8 3dsig_r10):
#   pop=1000 gen=2000 (2e6 evals) restarts=10
#   matrix MD5 72d7c7396702331d96ff12d18f831796
#   CLI: FlexAID --legacy
#   Arm order: A (CF) -> B0 (master CF) -> B (TEMPER 21 + single FO MinPts)
# Live OUT namespace (separate from pilot8):
#   $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/3dsig_full85_r10/
#
# Usage:
#   bash scripts/run_3dsig_red_pair_full85.sh --dry-run   # plan only (safe while pilot8 live)
#   bash scripts/run_3dsig_red_pair_full85.sh --check     # coverage + refuse status; no launch
#   bash scripts/run_3dsig_red_pair_full85.sh             # start full85 AFTER pilot8 is dead
#   bash scripts/run_3dsig_red_pair_full85.sh --from B0
#   bash scripts/run_3dsig_red_pair_full85.sh --only A
#
# Protocol: docs/implementation/3dsig_red_pair_protocol.md
# Parent:   docs/implementation/3dsig_shannon_ranking.md
#
# macOS note: written for /bin/bash 3.2 (no mapfile / no associative arrays).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_local_first_benchmark_storage.sh"

FROM="A"
ONLY=""
DRY=0
CHECK_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --from) FROM="${2:-A}"; shift 2 ;;
    --only) ONLY="${2:-}"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --check) CHECK_ONLY=1; shift ;;
    -h|--help) sed -n '1,30p' "$0"; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

export FLEXAID_POP="${FLEXAID_POP:-1000}"
export FLEXAID_GEN="${FLEXAID_GEN:-2000}"
export FLEXAID_RESTARTS="${FLEXAID_RESTARTS:-1}"  # first pass R=1; raise later to resume multi-restart
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export FLEXAIDDS_PARALLEL_RESTARTS="${FLEXAIDDS_PARALLEL_RESTARTS:-0}"

CAMPAIGN="${FLEXAID_CAMPAIGN:-3dsig_full85_r1}"
MATRIX_PIN="72d7c7396702331d96ff12d18f831796"
MAT="$FLEXAIDDS_QUEUE_ROOT/data/MC_st0r5.2_6.dat"
[[ -f "$MAT" ]] || MAT="$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat"
GOT=$(md5 -q "$MAT" 2>/dev/null || true)
[[ "$GOT" == "$MATRIX_PIN" ]] || { echo "FAIL: matrix MD5 '$GOT' != $MATRIX_PIN" >&2; exit 90; }
echo "OK: matrix md5=$GOT path=$MAT"

# --- Input coverage (Astex Diverse N=85) via python (bash 3.2 safe) ----------
INP_DIR="$FLEXAIDDS_QUEUE_ROOT/inputs/astex_diverse"
EXPECTED_YAML="$ROOT/benchmarks/datasets/astex_diverse.yaml"
COV_JSON=$(
  python3 - "$INP_DIR" "$EXPECTED_YAML" <<'PY'
import json, re, sys
from pathlib import Path

inp = Path(sys.argv[1])
yaml_path = Path(sys.argv[2])
present = sorted(p.name for p in inp.iterdir() if p.is_dir()) if inp.is_dir() else []
expected = []
if yaml_path.is_file():
    ids = re.findall(r"^\s*-\s*[\"']?([0-9][A-Za-z0-9]{3})[\"']?\s*$", yaml_path.read_text(), re.M)
    expected = sorted({i.upper() for i in ids})
exp_set = set(expected) if expected else set(present)
pres_set = set(present)
missing = sorted(exp_set - pres_set)
extra = sorted(pres_set - exp_set) if expected else []
incomplete = []
for pdb in present:
    d = inp / pdb
    apo = (d / f"{pdb}_apo.pdb").is_file()
    lig = (d / f"{pdb}_ligand.sdf").is_file() or (d / f"{pdb}_ligand.mol2").is_file()
    site = (
        (d / f"{pdb}_site.pdb").is_file()
        or (d / f"{pdb}_binding_site.pdb").is_file()
    )
    if not (apo and lig and site):
        incomplete.append(
            f"{pdb}(apo={int(apo)} lig={int(lig)} site={int(site)})"
        )
print(
    json.dumps(
        {
            "present": len(present),
            "expected": len(expected) if expected else len(present),
            "missing": missing,
            "extra": extra,
            "incomplete": incomplete,
            "present_ids": present,
        }
    )
)
PY
)
N_PRESENT=$(python3 -c 'import json,sys; print(json.loads(sys.argv[1])["present"])' "$COV_JSON")
N_EXPECTED=$(python3 -c 'import json,sys; print(json.loads(sys.argv[1])["expected"])' "$COV_JSON")
N_MISSING=$(python3 -c 'import json,sys; print(len(json.loads(sys.argv[1])["missing"]))' "$COV_JSON")
N_INCOMPLETE=$(python3 -c 'import json,sys; print(len(json.loads(sys.argv[1])["incomplete"]))' "$COV_JSON")
MISSING_STR=$(python3 -c 'import json,sys; print(" ".join(json.loads(sys.argv[1])["missing"]))' "$COV_JSON")
INCOMPLETE_STR=$(python3 -c 'import json,sys; print(" ".join(json.loads(sys.argv[1])["incomplete"]))' "$COV_JSON")

echo "OK: queue inputs present=$N_PRESENT expected=$N_EXPECTED missing=$N_MISSING incomplete=$N_INCOMPLETE"
if [[ "$N_MISSING" != "0" ]]; then
  echo "MISSING PDB ids ($N_MISSING): $MISSING_STR" >&2
fi
if [[ "$N_INCOMPLETE" != "0" ]]; then
  echo "INCOMPLETE inputs ($N_INCOMPLETE): $INCOMPLETE_STR" >&2
fi

# --- Refuse guards (hard) ----------------------------------------------------
# Exit 92 = REFUSE (live conflict). --dry-run / --check still report status.
refuse_reasons=""

add_refuse() {
  if [[ -z "$refuse_reasons" ]]; then
    refuse_reasons="$1"
  else
    refuse_reasons="${refuse_reasons}
$1"
  fi
}

# 1) pilot8 serial chain (pid file)
PILOT_PIDF="$FLEXAIDDS_LOCAL_LOGDIR/run_3dsig_red_pair_serial.pid"
if [[ -f "$PILOT_PIDF" ]]; then
  _pp=$(cat "$PILOT_PIDF" 2>/dev/null || true)
  if [[ -n "$_pp" ]] && kill -0 "$_pp" 2>/dev/null; then
    add_refuse "pilot8 serial chain live pid=$_pp ($PILOT_PIDF)"
  fi
fi
# process-name detection (pid file can be stale while children live)
if pgrep -f 'scripts/run_3dsig_red_pair_serial\.sh' >/dev/null 2>&1; then
  add_refuse "pilot8 run_3dsig_red_pair_serial.sh process live"
fi
# pilot8 arm runners (without --full85)
if pgrep -lf 'scripts/run_flexaid_arm_pilot8\.sh' 2>/dev/null | grep -v -- '--full85' | grep -q .; then
  _line=$(pgrep -lf 'scripts/run_flexaid_arm_pilot8\.sh' 2>/dev/null | grep -v -- '--full85' | head -1)
  add_refuse "pilot8 arm runner live: $_line"
fi
for _lk in \
  "$FLEXAIDDS_LOCAL_LOGDIR/run_A_pilot8.lock" \
  "$FLEXAIDDS_LOCAL_LOGDIR/run_B0_pilot8.lock" \
  "$FLEXAIDDS_LOCAL_LOGDIR/run_B_pilot8.lock"
do
  if [[ -f "$_lk" ]]; then
    _lp=$(cat "$_lk" 2>/dev/null || true)
    if [[ -n "$_lp" ]] && kill -0 "$_lp" 2>/dev/null; then
      add_refuse "pilot8 lock live $_lk pid=$_lp"
    fi
  fi
done

# 2) Any FlexAID --legacy dock
if pgrep -f 'FlexAID --legacy' >/dev/null 2>&1; then
  add_refuse "FlexAID --legacy process(es) live"
  pgrep -lf 'FlexAID --legacy' 2>/dev/null | head -5 || true
fi

# 3) Self full85 chain already live
FULL_PIDF="$FLEXAIDDS_LOCAL_LOGDIR/run_3dsig_red_pair_full85.pid"
if [[ -f "$FULL_PIDF" ]]; then
  _fp=$(cat "$FULL_PIDF" 2>/dev/null || true)
  if [[ -n "$_fp" ]] && kill -0 "$_fp" 2>/dev/null; then
    add_refuse "full85 serial chain already live pid=$_fp"
  fi
fi

# 4) Disk headroom soft note
AVAIL_GB=$(df -g "$FLEXAIDDS_LOCAL_ROOT" 2>/dev/null | awk 'NR==2{print $4}')
if [[ -n "${AVAIL_GB:-}" ]] && [[ "$AVAIL_GB" -lt 15 ]]; then
  echo "WARN: low free space on local root: ${AVAIL_GB}G (want >=15G headroom for full85 poses)" >&2
fi

N_REFUSE=0
if [[ -n "$refuse_reasons" ]]; then
  N_REFUSE=$(printf '%s\n' "$refuse_reasons" | grep -c . || true)
  echo "REFUSE STATUS: $N_REFUSE conflict(s)" >&2
  printf '%s\n' "$refuse_reasons" | while IFS= read -r r; do
    [[ -n "$r" ]] && echo "  - $r" >&2
  done
else
  echo "OK: no pilot8 / FlexAID --legacy conflicts detected"
fi

if (( CHECK_ONLY )); then
  echo "CHECK-ONLY complete (no launch)"
  echo "  OUT namespace: \$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/$CAMPAIGN/"
  echo "  knobs: pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS OMP=$OMP_NUM_THREADS"
  if [[ "$N_MISSING" != "0" || "$N_INCOMPLETE" != "0" ]]; then
    exit 91
  fi
  if [[ -n "$refuse_reasons" ]]; then
    exit 92
  fi
  exit 0
fi

# Hard refuse for real launch (dry-run still allowed to exercise prep plan)
if (( DRY == 0 )) && [[ -n "$refuse_reasons" ]]; then
  echo "REFUSE: cannot start full85 while pilot8/legacy docks are live (exit 92)" >&2
  echo "  Wait for pilot8 A -> B0 -> B to finish, then re-run." >&2
  echo "  Safe now: bash scripts/run_3dsig_red_pair_full85.sh --dry-run" >&2
  exit 92
fi

if (( DRY == 0 )) && [[ "$N_MISSING" != "0" || "$N_INCOMPLETE" != "0" ]]; then
  echo "REFUSE: input coverage incomplete (missing=$N_MISSING incomplete=$N_INCOMPLETE)" >&2
  exit 91
fi

if (( DRY )) && [[ "$N_MISSING" != "0" || "$N_INCOMPLETE" != "0" ]]; then
  echo "WARN: dry-run continues despite coverage gaps" >&2
fi

LOGDIR="$FLEXAIDDS_LOCAL_LOGDIR"
mkdir -p "$LOGDIR" "$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine"
TS=$(date -u +%Y-%m-%dT%H:%M:%SZ)
STATUS_MD="$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/3DSIG_RED_PAIR_FULL85.md"
{
  echo "# 3Dsig red-pair FULL85 (local-first)"
  echo
  echo "**Updated:** $TS"
  echo "**Status:** prepared — launch only when pilot8 is dead"
  echo "**Live OUT:** \$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/$CAMPAIGN/"
  echo "**Work:** \$FLEXAID_WORK_ROOT (local APFS)"
  echo "**iCloud:** deferred — bash scripts/sync_three_engine_local_to_icloud.sh --campaign $CAMPAIGN"
  echo "**C0:** SUSPENDED (do not dual-launch)"
  echo "**Knobs:** pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS OMP=$OMP_NUM_THREADS"
  echo "**Matrix:** $MATRIX_PIN"
  echo "**CLI:** FlexAID --legacy"
  echo "**Arm B:** TEMPER 21 + CLUSTA FO + single literature MinPts (no ladder)"
  echo "**Serial only on ~18 GiB Mac**"
  echo "**Inputs:** present=$N_PRESENT expected=$N_EXPECTED missing=$N_MISSING"
} >"$STATUS_MD"
echo "OK: wrote $STATUS_MD"

run_arm() {
  local arm="$1"
  export FLEXAID_ARM_OUT="$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/$arm/$CAMPAIGN"
  mkdir -p "$FLEXAID_ARM_OUT"
  echo "=========================================="
  echo "3DSIG red-pair FULL85 arm=$arm (LOCAL)"
  echo "OUT=$FLEXAID_ARM_OUT"
  echo "WORK=$FLEXAID_WORK_ROOT"
  echo "Q=$FLEXAIDDS_QUEUE_ROOT"
  echo "pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS"
  echo "=========================================="
  if (( DRY )); then
    bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" "$arm" --full85 --dry-run
    return 0
  fi
  extra=(--full85)
  # Default force re-prep so SHARESCL 10 + clean apo land (FLEXAID_FORCE_PREP=0 to skip)
  if [[ "${FLEXAID_FORCE_PREP:-1}" == "1" ]]; then extra+=(--force); fi
  [[ "${FLEXAID_SKIP_PREP:-0}" == "1" ]] && extra+=(--no-prepare)
  NOHUP=0 bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" "$arm" "${extra[@]}"
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

CHAIN_LOG="$LOGDIR/run_3dsig_red_pair_full85.log"
CHAIN_PID="$LOGDIR/run_3dsig_red_pair_full85.pid"

if (( DRY )); then
  for a in "${ARMS[@]}"; do run_arm "$a"; done
  echo "DRY-RUN complete (FULL85 not launched)"
  echo "  would OUT: $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/$CAMPAIGN/"
  echo "  arms: ${ARMS[*]}"
  exit 0
fi

if [[ -f "$CHAIN_PID" ]] && kill -0 "$(cat "$CHAIN_PID")" 2>/dev/null; then
  echo "REFUSE: full85 serial chain already live pid=$(cat "$CHAIN_PID")" >&2
  exit 92
fi

{
  echo "$$" >"$CHAIN_PID"
  echo "[$TS] START full85 local-first arms=${ARMS[*]} pop=$FLEXAID_POP gen=$FLEXAID_GEN R=$FLEXAID_RESTARTS campaign=$CAMPAIGN"
  for a in "${ARMS[@]}"; do
    run_arm "$a" || echo "WARN: arm $a non-zero"
  done
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] COMPLETE red-pair full85 serial (local)"
  bash "$ROOT/scripts/sync_three_engine_local_to_icloud.sh" --campaign "$CAMPAIGN" >>"$CHAIN_LOG" 2>&1 || true
  rm -f "$CHAIN_PID"
} >>"$CHAIN_LOG" 2>&1 &

echo "STARTED 3Dsig red-pair FULL85 serial (LOCAL-FIRST) pid=$! log=$CHAIN_LOG"
echo "  arms: ${ARMS[*]}"
echo "  OUT:  $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B}/$CAMPAIGN/"
echo "  sync later: bash scripts/sync_three_engine_local_to_icloud.sh --campaign $CAMPAIGN"
