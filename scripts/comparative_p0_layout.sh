#!/usr/bin/env bash
# comparative_p0_layout.sh — Phase 0 (Pins & storage) for the comparative pipeline.
#
# Implements docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md Phase 0:
#   - Local-first APFS dirs under $FLEXAIDDS_LOCAL_ROOT
#   - Matrix pin MC_st0r5.2_6.dat MD5 9dc93717dfed0698006d88dd6a9627bc
#
# Idempotent. Does not run docks. Does not touch CloudDocs for live I/O.
#
# Usage:
#   bash scripts/comparative_p0_layout.sh
#   FLEXAIDDS_LOCAL_ROOT=/tmp/foo bash scripts/comparative_p0_layout.sh
#
# Exit codes:
#   0 — layout OK and matrix MD5 matches (prints PHASE=P0 status=pass)
#   1 — layout or matrix failure (prints PHASE=P0 status=fail)
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"

EXPECTED_MATRIX_MD5="9dc93717dfed0698006d88dd6a9627bc"
MATRIX_NAME="MC_st0r5.2_6.dat"

fail() {
  echo "PHASE=P0 status=fail" >&2
  echo "ERROR: $*" >&2
  exit 1
}

pass() {
  echo "PHASE=P0 status=pass"
  exit 0
}

# Base local-first layout (campaigns/, logs/, pins/, three_engine_entropy_q1/{bin/C,data,inputs})
if [[ ! -f "$ROOT/scripts/ensure_local_first_layout.sh" ]]; then
  fail "missing scripts/ensure_local_first_layout.sh under FLEXAIDDS_ROOT=$FLEXAIDDS_ROOT"
fi
# shellcheck disable=SC1091
# ensure_local_first_layout sources claim_local_staging_paths and mkdir -p base dirs.
bash "$ROOT/scripts/ensure_local_first_layout.sh" || fail "ensure_local_first_layout.sh failed"

LOCAL="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
TEQ1="$LOCAL/three_engine_entropy_q1"
MATRIX_DST="$TEQ1/data/$MATRIX_NAME"

# Comparative dirs beyond the base ensure_local_first_layout set
mkdir -p \
  "$TEQ1/bin/A" \
  "$TEQ1/bin/B" \
  "$TEQ1/bin/C" \
  "$TEQ1/data" \
  "$TEQ1/inputs" \
  "$LOCAL/campaigns/three_engine/A" \
  "$LOCAL/campaigns/three_engine/B0" \
  "$LOCAL/campaigns/three_engine/B" \
  "$LOCAL/campaigns/three_engine/C" \
  "$LOCAL/campaigns/three_engine/analysis" \
  "$LOCAL/logs"

REQUIRED_DIRS=(
  "$TEQ1/bin/A"
  "$TEQ1/bin/B"
  "$TEQ1/bin/C"
  "$TEQ1/data"
  "$TEQ1/inputs"
  "$LOCAL/campaigns/three_engine/A"
  "$LOCAL/campaigns/three_engine/B0"
  "$LOCAL/campaigns/three_engine/B"
  "$LOCAL/campaigns/three_engine/C"
  "$LOCAL/campaigns/three_engine/analysis"
  "$LOCAL/logs"
)

echo "=== comparative_p0_layout ==="
echo "FLEXAIDDS_ROOT=$FLEXAIDDS_ROOT"
echo "FLEXAIDDS_LOCAL_ROOT=$LOCAL"
echo "dirs:"
for d in "${REQUIRED_DIRS[@]}"; do
  if [[ -d "$d" ]]; then
    echo "  OK  $d"
  else
    fail "missing directory: $d"
  fi
done

# Matrix: copy from repo if missing on the live local path
if [[ ! -f "$MATRIX_DST" ]]; then
  SRC=""
  if [[ -f "$ROOT/$MATRIX_NAME" ]]; then
    SRC="$ROOT/$MATRIX_NAME"
  elif [[ -f "$FLEXAIDDS_ROOT/$MATRIX_NAME" ]]; then
    SRC="$FLEXAIDDS_ROOT/$MATRIX_NAME"
  else
    fail "matrix missing at $MATRIX_DST and no source at $ROOT/$MATRIX_NAME or \$FLEXAIDDS_ROOT/$MATRIX_NAME"
  fi
  echo "copying matrix: $SRC -> $MATRIX_DST"
  cp "$SRC" "$MATRIX_DST"
fi

if [[ ! -f "$MATRIX_DST" ]]; then
  fail "matrix still missing after copy: $MATRIX_DST"
fi

# MD5 pin (macOS: md5 -q; Linux: md5sum)
if command -v md5 >/dev/null 2>&1; then
  GOT_MD5="$(md5 -q "$MATRIX_DST")"
elif command -v md5sum >/dev/null 2>&1; then
  GOT_MD5="$(md5sum "$MATRIX_DST" | awk '{print $1}')"
else
  fail "neither md5 nor md5sum available to verify matrix"
fi

echo "matrix_path=$MATRIX_DST"
echo "matrix_md5=$GOT_MD5"
echo "matrix_md5_expected=$EXPECTED_MATRIX_MD5"

if [[ "$GOT_MD5" != "$EXPECTED_MATRIX_MD5" ]]; then
  fail "matrix MD5 mismatch (got $GOT_MD5, expected $EXPECTED_MATRIX_MD5)"
fi

echo "OK comparative Phase 0 layout + matrix pin under $LOCAL"
pass
