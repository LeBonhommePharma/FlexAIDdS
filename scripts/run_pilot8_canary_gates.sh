#!/usr/bin/env bash
# P0 canary gates for pilot8 failure recovery (1P62 / 1T40 by default).
#
# Runs (per target):
#   1) ligand integrity on work LIG_ref (+ INI if present)
#   2) optional clean-apo dry-run report on work apo
#   3) native CF oracle gate when results/INI available
#
# Does NOT launch docking. Read-only except optional --json-dir reports.
#
# Usage:
#   bash scripts/run_pilot8_canary_gates.sh
#   bash scripts/run_pilot8_canary_gates.sh --arm B0 --pdb 1P62,1T40
#   WORK_ROOT=$HOME/flexaidds_results/three_engine_entropy_q1/work \\
#     bash scripts/run_pilot8_canary_gates.sh --results-root $HOME/flexaidds_results/.../C0_pilot
#
# Exit: 0 only if all integrity gates pass AND (if results present) oracle
# decisions are reported. Oracle pathology returns non-zero unless --report-only.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1090
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env" || true

ARM="${ARM:-B0}"
PDBS="${PDBS:-1P62,1T40}"
WORK_ROOT="${WORK_ROOT:-${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/three_engine_entropy_q1/work}"
RESULTS_ROOT="${RESULTS_ROOT:-}"
JSON_DIR=""
REPORT_ONLY=0
REQUIRE_INI=0
MAX_BOND="${MAX_BOND:-3.0}"
TOLERANCE="${TOLERANCE:-0.0}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --arm) ARM="${2:-}"; shift 2 ;;
    --pdb) PDBS="${2:-}"; shift 2 ;;
    --work-root) WORK_ROOT="${2:-}"; shift 2 ;;
    --results-root) RESULTS_ROOT="${2:-}"; shift 2 ;;
    --json-dir) JSON_DIR="${2:-}"; shift 2 ;;
    --report-only) REPORT_ONLY=1; shift ;;
    --require-ini) REQUIRE_INI=1; shift ;;
    --max-bond) MAX_BOND="${2:-}"; shift 2 ;;
    --tolerance) TOLERANCE="${2:-}"; shift 2 ;;
    -h|--help)
      sed -n '2,25p' "$0"
      exit 0
      ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

VAL_LIG="$ROOT/scripts/validate_ligand_integrity.py"
ORACLE="$ROOT/scripts/native_cf_oracle_gate.py"
CLEAN="$ROOT/scripts/clean_target_apo.py"

for f in "$VAL_LIG" "$ORACLE" "$CLEAN"; do
  [[ -f "$f" ]] || { echo "FAIL: missing $f" >&2; exit 2; }
done

IFS=',' read -r -a TARGETS <<< "$PDBS"
FAIL=0

echo "=== pilot8 canary gates ==="
echo "arm=$ARM work_root=$WORK_ROOT"
echo "targets=${TARGETS[*]}"
[[ -n "$RESULTS_ROOT" ]] && echo "results_root=$RESULTS_ROOT"
[[ -n "$JSON_DIR" ]] && mkdir -p "$JSON_DIR"

for pdb in "${TARGETS[@]}"; do
  pdb="${pdb^^}"
  pdb="${pdb// /}"
  [[ -z "$pdb" ]] && continue
  w="$WORK_ROOT/$ARM/$pdb"
  echo ""
  echo "── $pdb  work=$w"

  if [[ ! -d "$w" ]]; then
    echo "FAIL: work dir missing: $w"
    FAIL=1
    continue
  fi

  # 1) Ligand integrity (LIG_ref; INI if --require-ini)
  jarg=()
  if [[ -n "$JSON_DIR" ]]; then
    jarg=(--json "$JSON_DIR/${pdb}_ligand_integrity.json")
  fi
  ini_arg=()
  (( REQUIRE_INI )) && ini_arg=(--require-ini)
  set +e
  python3 "$VAL_LIG" --work "$w" --max-bond "$MAX_BOND" "${ini_arg[@]}" "${jarg[@]}"
  rc_lig=$?
  set -e
  if (( rc_lig != 0 )); then
    echo "ligand_integrity exit=$rc_lig"
    FAIL=1
  fi

  # 2) Clean-apo dry-run on work apo (informational)
  apo=""
  for cand in "$w/${pdb}_apo.pdb" "$w/TARGET.inp.pdb"; do
    if [[ -f "$cand" ]]; then
      apo="$cand"
      break
    fi
  done
  if [[ -n "$apo" ]]; then
    echo "clean_target_apo dry-run on $apo"
    python3 "$CLEAN" "$apo" --dry-run || true
  fi

  # 3) Native CF oracle when results available
  res=""
  if [[ -n "$RESULTS_ROOT" && -d "$RESULTS_ROOT/$pdb" ]]; then
    res="$RESULTS_ROOT/$pdb"
  elif [[ -d "$w" ]] && { [[ -f "$w/result.csv" ]] || compgen -G "$w/*_0.pdb" >/dev/null; }; then
    res="$w"
  fi

  if [[ -n "$res" ]]; then
    jarg2=()
    if [[ -n "$JSON_DIR" ]]; then
      jarg2=(--json "$JSON_DIR/${pdb}_native_cf_oracle.json")
    fi
    set +e
    python3 "$ORACLE" --work "$w" --results "$res" --pdb "$pdb" \
      --tolerance "$TOLERANCE" "${jarg2[@]}"
    rc_or=$?
    set -e
    echo "native_cf_oracle exit=$rc_or"
    if (( rc_or != 0 && REPORT_ONLY == 0 )); then
      FAIL=1
    fi
  else
    echo "NOTE: no results dir for $pdb — oracle deferred (prep integrity only)."
    echo "      After FlexAID: $ORACLE --work $w --results <OUT>/$pdb"
  fi
done

echo ""
if (( FAIL )); then
  echo "CANARY GATES: FAIL"
  exit 1
fi
echo "CANARY GATES: PASS (or oracle deferred without results)"
exit 0
