#!/usr/bin/env bash
# Serial 3Dsig red-pair FULL Astex Diverse (N=85): A then B0 then B (LOCAL-FIRST)
#
# PREPARE-ONLY while pilot8 is live: this launcher REFUSES to start if the
# pilot8 chain or any FlexAID --legacy dock is running.
#
# Deck knobs (same as pilot8 3dsig_r10):
#   pop=1000 gen=2000 (2e6 evals) restarts=10
#   matrix MD5 9dc93717dfed0698006d88dd6a9627bc (repo/baseline-validated)
#   CLI: FlexAID --legacy + SOFTWA 0.40 + sphere prune + metal-near-ligand
#   Arm order (default): A → B  (B0 skip unless FLEXAID_INCLUDE_B0=1)
#   Arm C FO@298K only after native CF oracle PASS
# Live OUT namespace:
#   $FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/{A,B0,B,C}/$CAMPAIGN/
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
# Claim-comparable default: R=10. Diagnostic R=1 only if FLEXAID_DIAGNOSTIC=1.
# FLEXAID_CLAIM_MODE=1 adds hard gates (oracle PASS + A≠B binary split).
CLAIM_MODE="${FLEXAID_CLAIM_MODE:-0}"
DIAGNOSTIC="${FLEXAID_DIAGNOSTIC:-0}"
if [[ "$DIAGNOSTIC" == "1" ]]; then
  export FLEXAID_RESTARTS="${FLEXAID_RESTARTS:-1}"
else
  # Recovery / claim path: multi-restart
  export FLEXAID_RESTARTS="${FLEXAID_RESTARTS:-10}"
fi
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export FLEXAIDDS_PARALLEL_RESTARTS="${FLEXAIDDS_PARALLEL_RESTARTS:-0}"
# Softβ always OFF unless operator explicitly overrides after oracle PASS
export FLEXAIDDS_SOFTBETA_ELECTION="${FLEXAIDDS_SOFTBETA_ELECTION:-0}"
# Soft-core wall on classic path (SOFTWA in CONFIG.inp)
export FLEXAIDDS_SOFT_WALL="${FLEXAIDDS_SOFT_WALL:-0.40}"
# Sphere prune default for new preps (oversized LOCCLF); 0 disables
export FLEXAIDDS_SPHERE_MAX="${FLEXAIDDS_SPHERE_MAX:-2500}"
export FLEXAIDDS_SPHERE_MAX_DIST="${FLEXAIDDS_SPHERE_MAX_DIST:-8.0}"
# Keep catalytic metals within 4 Å of crystal ligand
export FLEXAIDDS_KEEP_METALS_NEAR_LIGAND="${FLEXAIDDS_KEEP_METALS_NEAR_LIGAND:-4.0}"

CAMPAIGN="${FLEXAID_CAMPAIGN:-3dsig_full85_r10_cf_fix}"
MATRIX_PIN="9dc93717dfed0698006d88dd6a9627bc"
MAT="$FLEXAIDDS_QUEUE_ROOT/data/MC_st0r5.2_6.dat"
[[ -f "$MAT" ]] || MAT="$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat"
GOT=$(md5 -q "$MAT" 2>/dev/null || true)
[[ "$GOT" == "$MATRIX_PIN" ]] || { echo "FAIL: matrix MD5 '$GOT' != $MATRIX_PIN" >&2; exit 90; }
echo "OK: matrix md5=$GOT path=$MAT"

# --- Claim / science gates (serial only; no dual-launch) ---------------------
# B0 is NOT a scientific control when bin/A SHA == bin/B SHA (deterministic twin).
BIN_ROOT="${FLEXAIDDS_LOCAL_ROOT}/three_engine_entropy_q1/bin"
if [[ "$CLAIM_MODE" == "1" ]]; then
  if (( FLEXAID_RESTARTS < 10 )); then
    echo "FAIL: FLEXAID_CLAIM_MODE=1 requires FLEXAID_RESTARTS>=10 (got $FLEXAID_RESTARTS)" >&2
    exit 93
  fi
  if ! python3 "$ROOT/scripts/stage_three_engine_bins.py" --check-claim-split --dest "$BIN_ROOT"; then
    echo "FAIL: claim mode needs distinct historical A vs master B binaries" >&2
    exit 94
  fi
  ORACLE_STATUS="${FLEXAID_ORACLE_STATUS:-$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/${CAMPAIGN}_oracle_status.json}"
  if [[ ! -f "$ORACLE_STATUS" ]]; then
    echo "FAIL: claim mode needs native CF oracle status at $ORACLE_STATUS" >&2
    echo "  Run: python3 scripts/run_panel_native_cf_oracle.py --out-dir ... --status-out $ORACLE_STATUS" >&2
    exit 95
  fi
  if ! python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); sys.exit(0 if s.get("ranking_allowed") else 1)' "$ORACLE_STATUS"; then
    echo "FAIL: native CF oracle ranking_allowed=false — Softβ/arm-C/claim ranking FORBIDDEN" >&2
    exit 96
  fi
  echo "OK: claim mode gates (R>=10, binary split, oracle PASS)"
else
  if [[ "$DIAGNOSTIC" == "1" ]]; then
    echo "OK: diagnostic mode (R=$FLEXAID_RESTARTS) — not 3Dsig-claimable"
  else
    echo "OK: recovery path R=$FLEXAID_RESTARTS Softβ=OFF — set FLEXAID_CLAIM_MODE=1 after oracle PASS for claim gates"
  fi
  # Warn if A==B (B0 twin)
  if [[ -x "$BIN_ROOT/A/FlexAID" && -x "$BIN_ROOT/B/FlexAID" ]]; then
    if ! python3 "$ROOT/scripts/stage_three_engine_bins.py" --check-claim-split --dest "$BIN_ROOT" 2>/dev/null; then
      echo "WARN: bin A SHA == bin B SHA → B0 is a deterministic twin of A (not science control)" >&2
    fi
  fi
fi
# Refuse Softβ unless oracle allows
if [[ "${FLEXAIDDS_SOFTBETA_ELECTION}" != "0" && "${FLEXAIDDS_SOFTBETA_ELECTION}" != "false" ]]; then
  ORACLE_STATUS="${FLEXAID_ORACLE_STATUS:-$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/${CAMPAIGN}_oracle_status.json}"
  if [[ ! -f "$ORACLE_STATUS" ]] || ! python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); sys.exit(0 if s.get("softbeta_allowed") else 1)' "$ORACLE_STATUS" 2>/dev/null; then
    echo "FAIL: Softβ election requested but oracle softbeta_allowed is false/missing" >&2
    exit 97
  fi
fi

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

# Science chain (serial, no dual-launch):
#   default: A (historical pin) → B (master TEMPER21 FO)
#   B0: optional only — NOT a claim control when A≡B binary; set FLEXAID_INCLUDE_B0=1
#   C: FO@298K only after native CF oracle PASS; set FLEXAID_INCLUDE_ARM_C=1
#   FLEXAID_ARMS="A,B0,B,C" overrides (still enforces C oracle gate)
INCLUDE_B0="${FLEXAID_INCLUDE_B0:-0}"
INCLUDE_C="${FLEXAID_INCLUDE_ARM_C:-0}"
ORACLE_STATUS="${FLEXAID_ORACLE_STATUS:-$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/${CAMPAIGN}_oracle_status.json}"

oracle_allows_c() {
  [[ -f "$ORACLE_STATUS" ]] && python3 -c 'import json,sys; s=json.load(open(sys.argv[1])); sys.exit(0 if s.get("arm_c_fo298_allowed") else 1)' "$ORACLE_STATUS" 2>/dev/null
}

if [[ -n "${FLEXAID_ARMS:-}" ]]; then
  IFS=',' read -r -a ARMS <<< "$FLEXAID_ARMS"
elif [[ -n "$ONLY" ]]; then
  ARMS=("$ONLY")
else
  case "$FROM" in
    A)
      ARMS=(A)
      if [[ "$INCLUDE_B0" == "1" ]]; then ARMS+=(B0); fi
      ARMS+=(B)
      if [[ "$INCLUDE_C" == "1" ]]; then ARMS+=(C); fi
      ;;
    B0) ARMS=(B0 B); [[ "$INCLUDE_C" == "1" ]] && ARMS+=(C) ;;
    B) ARMS=(B); [[ "$INCLUDE_C" == "1" ]] && ARMS+=(C) ;;
    C) ARMS=(C) ;;
    *) echo "bad --from $FROM" >&2; exit 2 ;;
  esac
fi

# Enforce arm-C oracle gate whenever C is in the chain
for _a in "${ARMS[@]}"; do
  if [[ "$_a" == "C" ]]; then
    if [[ "${FLEXAIDDS_ALLOW_ARM_C:-0}" != "1" ]] && ! oracle_allows_c; then
      echo "FAIL: arm C requires oracle arm_c_fo298_allowed=true in $ORACLE_STATUS" >&2
      echo "  (or FLEXAIDDS_ALLOW_ARM_C=1 diagnostic-only). Softβ/FO@298K blocked." >&2
      exit 98
    fi
  fi
done
echo "OK: arm order ${ARMS[*]} (INCLUDE_B0=$INCLUDE_B0 INCLUDE_C=$INCLUDE_C CLAIM_MODE=$CLAIM_MODE R=$FLEXAID_RESTARTS SOFT_WALL=$FLEXAIDDS_SOFT_WALL)"

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
