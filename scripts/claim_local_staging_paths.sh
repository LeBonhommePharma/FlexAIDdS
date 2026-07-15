#!/usr/bin/env bash
# claim_local_staging_paths.sh — resolve local work + iCloud mirror paths.
#
# Intermediate anti-hang policy:
#   - Live GA I/O runs on local APFS outside CloudDocs (FLEXAIDDS_LOCAL_ROOT).
#   - iCloud remains the durable mirror for completed result.csv + thin metadata.
#   - A background loop (claim_icloud_sync_loop.sh) rsyncs local → iCloud every N min.
#
# Source this file; do not execute.
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
# shellcheck disable=SC2034

: "${FLEXAIDDS_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"
export FLEXAIDDS_QUEUE_ROOT="${FLEXAIDDS_QUEUE_ROOT:-$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1}"

# Prefer local queue staging for binaries when present (avoids CloudDocs fileprovider).
_LOCAL_Q="${FLEXAIDDS_LOCAL_ROOT}/three_engine_entropy_q1"
if [[ -x "$_LOCAL_Q/bin/C/FlexAIDdS" && -x "$_LOCAL_Q/bin/C/benchmark_datasets" ]]; then
  export FLEXAIDDS_LOCAL_QUEUE="$_LOCAL_Q"
else
  export FLEXAIDDS_LOCAL_QUEUE="${FLEXAIDDS_QUEUE_ROOT}"
fi

export C0_CAMPAIGN_ID="${C0_CAMPAIGN_ID:-C0_full85_claim_g2000_popmod_20260715}"
export C0_CLAIM_LOCAL_OUT="${C0_CLAIM_LOCAL_OUT:-$FLEXAIDDS_LOCAL_ROOT/campaigns/$C0_CAMPAIGN_ID}"
export C0_CLAIM_ICLOUD_OUT="${C0_CLAIM_ICLOUD_OUT:-$FLEXAIDDS_ICLOUD/results/campaigns/$C0_CAMPAIGN_ID}"
export C0_CLAIM_LOCAL_LOGDIR="${C0_CLAIM_LOCAL_LOGDIR:-$FLEXAIDDS_LOCAL_ROOT/logs/C0_claim}"
export C0_CLAIM_RUNNER="${C0_CLAIM_RUNNER:-$FLEXAIDDS_LOCAL_QUEUE/bin/C/benchmark_datasets}"
export C0_CLAIM_BINARY="${C0_CLAIM_BINARY:-$FLEXAIDDS_LOCAL_QUEUE/bin/C/FlexAIDdS}"
export C0_CLAIM_DATA_DIR="${C0_CLAIM_DATA_DIR:-$FLEXAIDDS_LOCAL_QUEUE/data}"
export C0_CLAIM_MANIFEST="${C0_CLAIM_MANIFEST:-$FLEXAIDDS_QUEUE_ROOT/inputs/astex_native_85.json}"

# Sync cadence (minutes) for the background loop
export C0_CLAIM_SYNC_EVERY_MIN="${C0_CLAIM_SYNC_EVERY_MIN:-5}"
# Stall: log mtime older than this (minutes) + CPU ~0 → ALERT
export C0_CLAIM_STALL_MIN="${C0_CLAIM_STALL_MIN:-20}"
