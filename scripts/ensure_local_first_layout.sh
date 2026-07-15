#!/usr/bin/env bash
# ensure_local_first_layout.sh — create local-first APFS layout for claims/ops.
#
# Idempotent. Does not touch CloudDocs. Safe to run before every claim launch.
#
# Creates under $FLEXAIDDS_LOCAL_ROOT (default ~/flexaidds_results):
#   campaigns/
#   logs/ops
#   logs/ops_monitor
#   logs/C0_claim
#   pins/materialize
#   three_engine_entropy_q1/{bin/C,data,inputs}
#
# Usage:
#   bash scripts/ensure_local_first_layout.sh
#   FLEXAIDDS_LOCAL_ROOT=/tmp/foo bash scripts/ensure_local_first_layout.sh
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"

# shellcheck disable=SC1091
source "$ROOT/scripts/claim_local_staging_paths.sh"

LOCAL="${FLEXAIDDS_LOCAL_ROOT}"
mkdir -p \
  "$LOCAL/campaigns" \
  "$LOCAL/logs/ops" \
  "$LOCAL/logs/ops_monitor" \
  "$LOCAL/logs/C0_claim" \
  "$LOCAL/pins/materialize" \
  "$LOCAL/three_engine_entropy_q1/bin/C" \
  "$LOCAL/three_engine_entropy_q1/data" \
  "$LOCAL/three_engine_entropy_q1/inputs"

# Campaign OUT + claim logdir (from claim_local_staging_paths)
mkdir -p "${C0_CLAIM_LOCAL_OUT}" "${C0_CLAIM_LOCAL_LOGDIR}"

echo "=== ensure_local_first_layout ==="
echo "FLEXAIDDS_ROOT=$FLEXAIDDS_ROOT"
echo "FLEXAIDDS_LOCAL_ROOT=$LOCAL"
echo "C0_CAMPAIGN_ID=${C0_CAMPAIGN_ID}"
echo "C0_CLAIM_LOCAL_OUT=${C0_CLAIM_LOCAL_OUT}"
echo "C0_CLAIM_LOCAL_LOGDIR=${C0_CLAIM_LOCAL_LOGDIR}"
echo "C0_CLAIM_ICLOUD_OUT=${C0_CLAIM_ICLOUD_OUT}"
echo "FLEXAIDDS_LOCAL_QUEUE=${FLEXAIDDS_LOCAL_QUEUE}"
echo "C0_CLAIM_RUNNER=${C0_CLAIM_RUNNER}"
echo "C0_CLAIM_BINARY=${C0_CLAIM_BINARY}"
echo "C0_CLAIM_DATA_DIR=${C0_CLAIM_DATA_DIR}"
echo "dirs:"
for d in \
  campaigns \
  logs/ops \
  logs/ops_monitor \
  logs/C0_claim \
  pins/materialize \
  three_engine_entropy_q1/bin/C \
  three_engine_entropy_q1/data \
  three_engine_entropy_q1/inputs
do
  p="$LOCAL/$d"
  if [[ -d "$p" ]]; then
    echo "  OK  $p"
  else
    echo "  MISSING  $p" >&2
    exit 1
  fi
done
echo "OK layout ready under $LOCAL"
exit 0
