#!/usr/bin/env bash
# wave3_preflight.sh — Wave 3 sampling / BCR-raiser preflight (no live dock).
#
# Checks:
#   1) Production matrix MC_st0r5.2_6.dat MD5 == 9dc93717dfed0698006d88dd6a9627bc
#   2) Seed-off / Softβ-off env echoes (claim hygiene)
#   3) Prints documented next commands (does NOT launch GA)
#
# Usage:
#   bash scripts/wave3_preflight.sh
#   FLEXAIDDS_LOCAL_ROOT=/path/to/local bash scripts/wave3_preflight.sh
#
# Exit codes:
#   0 — matrix pin OK (env hygiene printed even if seed-off vars unset)
#   1 — matrix missing or MD5 mismatch
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"

EXPECTED_MATRIX_MD5="9dc93717dfed0698006d88dd6a9627bc"
MATRIX_NAME="MC_st0r5.2_6.dat"

fail() {
  echo "WAVE3_PREFLIGHT status=fail" >&2
  echo "ERROR: $*" >&2
  exit 1
}

echo "=== wave3_preflight (no dock) ==="
echo "FLEXAIDDS_ROOT=$FLEXAIDDS_ROOT"

# Optional layout ensure (idempotent; local APFS only)
if [[ -f "$ROOT/scripts/ensure_local_first_layout.sh" ]]; then
  # shellcheck disable=SC1091
  bash "$ROOT/scripts/ensure_local_first_layout.sh" || true
fi

LOCAL="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
export FLEXAIDDS_LOCAL_ROOT="$LOCAL"
MATRIX_DST="$LOCAL/three_engine_entropy_q1/data/$MATRIX_NAME"

# Resolve matrix: live local path preferred; fall back to repo root copy
MATRIX_PATH=""
if [[ -f "$MATRIX_DST" ]]; then
  MATRIX_PATH="$MATRIX_DST"
elif [[ -f "$ROOT/$MATRIX_NAME" ]]; then
  MATRIX_PATH="$ROOT/$MATRIX_NAME"
  echo "WARN: using repo matrix at $MATRIX_PATH (live path missing: $MATRIX_DST)"
elif [[ -f "$FLEXAIDDS_ROOT/$MATRIX_NAME" ]]; then
  MATRIX_PATH="$FLEXAIDDS_ROOT/$MATRIX_NAME"
  echo "WARN: using FLEXAIDDS_ROOT matrix at $MATRIX_PATH"
else
  fail "matrix not found at $MATRIX_DST or $ROOT/$MATRIX_NAME"
fi

if command -v md5 >/dev/null 2>&1; then
  GOT_MD5="$(md5 -q "$MATRIX_PATH")"
elif command -v md5sum >/dev/null 2>&1; then
  GOT_MD5="$(md5sum "$MATRIX_PATH" | awk '{print $1}')"
else
  fail "neither md5 nor md5sum available"
fi

echo "matrix_path=$MATRIX_PATH"
echo "matrix_md5=$GOT_MD5"
echo "matrix_md5_expected=$EXPECTED_MATRIX_MD5"

if [[ "$GOT_MD5" != "$EXPECTED_MATRIX_MD5" ]]; then
  fail "matrix MD5 mismatch (got $GOT_MD5, expected $EXPECTED_MATRIX_MD5) — refuse 72d7 packing fork"
fi
echo "OK matrix pin 9dc9"

# ── Seed-off / claim hygiene env echoes (informational; does not setenv) ─────
echo ""
echo "=== claim hygiene env (echo only) ==="
echo "FLEXAIDDS_SEED_ELITISM=${FLEXAIDDS_SEED_ELITISM-<unset>}"
echo "FLEXAIDDS_NATIVE_SEED_FRAC=${FLEXAIDDS_NATIVE_SEED_FRAC-<unset>}"
echo "FLEXAIDDS_SOFTBETA_ELECTION=${FLEXAIDDS_SOFTBETA_ELECTION-<unset>}"
echo "FLEXAIDDS_ELECTION_SHANNON_F=${FLEXAIDDS_ELECTION_SHANNON_F-<unset>}"
echo "FLEXAIDDS_COM_BURIAL_CAP=${FLEXAIDDS_COM_BURIAL_CAP-<unset>}"
echo "FLEXAIDDS_BOOM_FRAC=${FLEXAIDDS_BOOM_FRAC-<unset>}"
echo "FLEXAIDDS_SHARING_ALPHA=${FLEXAIDDS_SHARING_ALPHA-<unset>}"
echo "FLEXAIDDS_COARSE_GRID_STEP=${FLEXAIDDS_COARSE_GRID_STEP-<unset>}"
echo "FLEXAIDDS_SEED_BASE=${FLEXAIDDS_SEED_BASE-<unset>}"

SEED_OK=1
if [[ -n "${FLEXAIDDS_SEED_ELITISM:-}" && "${FLEXAIDDS_SEED_ELITISM}" != "0" && "${FLEXAIDDS_SEED_ELITISM}" != "false" ]]; then
  echo "WARN: FLEXAIDDS_SEED_ELITISM is ON — claim path requires seed-off (set 0)"
  SEED_OK=0
fi
if [[ -n "${FLEXAIDDS_NATIVE_SEED_FRAC:-}" ]]; then
  # bash arithmetic: non-numeric → treat as non-zero warning path
  if [[ "${FLEXAIDDS_NATIVE_SEED_FRAC}" != "0" && "${FLEXAIDDS_NATIVE_SEED_FRAC}" != "0.0" ]]; then
    echo "WARN: FLEXAIDDS_NATIVE_SEED_FRAC=${FLEXAIDDS_NATIVE_SEED_FRAC} — claim path requires 0"
    SEED_OK=0
  fi
fi
if [[ -n "${FLEXAIDDS_SOFTBETA_ELECTION:-}" && "${FLEXAIDDS_SOFTBETA_ELECTION}" != "0" && "${FLEXAIDDS_SOFTBETA_ELECTION}" != "false" ]]; then
  echo "WARN: Softβ S1 ON — Wave 3 primary lever is sampling, not Softβ default ON"
fi
if [[ -n "${FLEXAIDDS_COM_BURIAL_CAP:-}" ]]; then
  echo "WARN: FLEXAIDDS_COM_BURIAL_CAP is set — REJECT as product default (see WAVE3 plan §3)"
fi

if [[ "$SEED_OK" -eq 1 ]]; then
  echo "OK seed-off hygiene (vars unset or explicitly off)"
else
  echo "WARN seed-off hygiene incomplete — fix before claim-style pilot"
fi

echo ""
echo "=== baseline (documented; not re-measured here) ==="
echo "campaign_ref=v_autonomous_20260724_160919"
echo "genuine_baseline=20/79=25.3%  BCR=22/79=27.8%  election_gap~2  seed_echo=0"
echo "goal_anchor_JCIM_top1=45.2% (goal metric only)"
echo "plan=docs/implementation/WAVE3_SAMPLING_BCR_PLAN.md"

echo ""
echo "=== next commands (documented; not executed) ==="
cat <<'EOF'
# 1) Optional canary (no dock):
#    bash scripts/run_pilot8_canary_gates.sh --pdb 1P62,1T40 --report-only
#
# 2) Comparative pipeline dry (no dock):
#    PYTHONPATH=$PWD/python python3 scripts/run_comparative_phases.py --pipeline-dry
#
# 3) Claim-style pilot only when operator authorizes (serial, seed-off, 9dc9):
#    export FLEXAIDDS_SEED_ELITISM=0 FLEXAIDDS_NATIVE_SEED_FRAC=0
#    export FLEXAIDDS_SOFTBETA_ELECTION=0
#    # micro-set e.g. 1P62,1T40,1G9V — see WAVE3_SAMPLING_BCR_PLAN.md §4
#    # Prefer local-first launchers (run_C0_claim_clean / DatasetRunner subset).
#    # P2 native_cf_oracle_gate must PASS before claim full85.
#    # No dual full85. No Softβ as primary BCR lever.
EOF

echo ""
echo "WAVE3_PREFLIGHT status=pass matrix_md5=$GOT_MD5"
exit 0
