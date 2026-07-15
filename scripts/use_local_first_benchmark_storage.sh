#!/usr/bin/env bash
# use_local_first_benchmark_storage.sh — live GA I/O on local APFS; iCloud = delayed mirror.
#
# Policy (anti-hang):
#   - OUT, work, logs, binaries, matrix, inputs: under $FLEXAIDDS_LOCAL_ROOT (default ~/flexaidds_results)
#   - iCloud ($FLEXAIDDS_ICLOUD) is durable archive only — rsync later via
#       scripts/sync_three_engine_local_to_icloud.sh
#   - Never require CloudDocs for a live dock or prepare step
#
# Source before any pilot / 3Dsig red-pair launch:
#   source scripts/use_local_first_benchmark_storage.sh
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
# shellcheck disable=SC2034

_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$_ROOT}"

export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"

# Local queue (binaries + data + inputs already staged)
export FLEXAIDDS_LOCAL_QUEUE="${FLEXAIDDS_LOCAL_QUEUE:-$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1}"
export FLEXAIDDS_QUEUE_ROOT="${FLEXAIDDS_QUEUE_ROOT:-$FLEXAIDDS_LOCAL_QUEUE}"

# Live results on local APFS (NOT CloudDocs)
export FLEXAIDDS_RESULTS="${FLEXAIDDS_RESULTS:-$FLEXAIDDS_LOCAL_ROOT}"
export FLEXAID_WORK_ROOT="${FLEXAID_WORK_ROOT:-$FLEXAIDDS_LOCAL_QUEUE/work}"
export FLEXAIDDS_LOCAL_LOGDIR="${FLEXAIDDS_LOCAL_LOGDIR:-$FLEXAIDDS_LOCAL_ROOT/logs/three_engine}"

# Allow run_flexaid_arm_pilot8 to skip require_icloud_out
export FLEXAIDDS_ALLOW_LOCAL_OUT=1

# iCloud mirror targets (sync later only)
export FLEXAIDDS_ICLOUD_RESULTS="${FLEXAIDDS_ICLOUD_RESULTS:-$FLEXAIDDS_ICLOUD/results}"
export FLEXAIDDS_ICLOUD_QUEUE="${FLEXAIDDS_ICLOUD_QUEUE:-$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1}"

mkdir -p \
  "$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine" \
  "$FLEXAID_WORK_ROOT" \
  "$FLEXAIDDS_LOCAL_LOGDIR" \
  "$FLEXAIDDS_LOCAL_QUEUE/bin/A" \
  "$FLEXAIDDS_LOCAL_QUEUE/bin/B" \
  "$FLEXAIDDS_LOCAL_QUEUE/data" \
  "$FLEXAIDDS_LOCAL_QUEUE/inputs" \
  2>/dev/null || true

# Refuse if someone still pointed RESULTS at CloudDocs while local-first is active
case "${FLEXAIDDS_RESULTS}" in
  *"/Mobile Documents/com~apple~CloudDocs/"*)
    echo "WARN: FLEXAIDDS_RESULTS was on iCloud ($FLEXAIDDS_RESULTS) — forcing local" >&2
    export FLEXAIDDS_RESULTS="$FLEXAIDDS_LOCAL_ROOT"
    ;;
esac

echo "OK local-first storage:"
echo "  FLEXAIDDS_LOCAL_ROOT=$FLEXAIDDS_LOCAL_ROOT"
echo "  FLEXAIDDS_QUEUE_ROOT=$FLEXAIDDS_QUEUE_ROOT"
echo "  FLEXAIDDS_RESULTS=$FLEXAIDDS_RESULTS"
echo "  FLEXAID_WORK_ROOT=$FLEXAID_WORK_ROOT"
echo "  FLEXAIDDS_LOCAL_LOGDIR=$FLEXAIDDS_LOCAL_LOGDIR"
echo "  iCloud mirror (later): $FLEXAIDDS_ICLOUD_RESULTS"
echo "  FLEXAIDDS_ALLOW_LOCAL_OUT=$FLEXAIDDS_ALLOW_LOCAL_OUT"
