#!/usr/bin/env bash
# use_icloud_benchmark_storage.sh — force production benchmark I/O onto iCloud Drive.
#
# Source this before any campaign launch:
#   source scripts/use_icloud_benchmark_storage.sh
#   # or: source ~/.flexaidds_env && source "$FLEXAIDDS_ROOT/scripts/use_icloud_benchmark_storage.sh"
#
# Policy (non-negotiable for production / claim runs):
#   - All result.csv, elected poses, RUN_RECEIPT, logs that define a campaign
#     MUST live under iCloud Drive:
#       $HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/
#   - Local ~/flexaidds_results is archive / binary staging ONLY (not new claim OUT).
#   - Mach-O binaries may stay on local disk (queue bin/ symlinks) to avoid sync corruption.
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
# shellcheck disable=SC2034

_ICLOUD_HOME="${HOME}/Library/Mobile Documents/com~apple~CloudDocs"
_ICLOUD_BENCH="${FLEXAIDDS_ICLOUD:-$_ICLOUD_HOME/FlexAIDdS_benchmarks}"

if [[ ! -d "$_ICLOUD_HOME" ]]; then
  echo "ERROR: iCloud Drive path missing: $_ICLOUD_HOME" >&2
  echo "Enable iCloud Drive and ensure CloudDocs is available." >&2
  return 91 2>/dev/null || exit 91
fi

export FLEXAIDDS_ICLOUD="$_ICLOUD_BENCH"
export FLEXAIDDS_RESULTS="${FLEXAIDDS_RESULTS:-$FLEXAIDDS_ICLOUD/results}"
export FLEXAIDDS_QUEUE_ROOT="${FLEXAIDDS_QUEUE_ROOT:-$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1}"
export FLEXAIDDS_WORKING="${FLEXAIDDS_WORKING:-$FLEXAIDDS_RESULTS/working}"

# Force RESULTS under iCloud even if a stale env pointed elsewhere
case "$FLEXAIDDS_RESULTS" in
  *"/Mobile Documents/com~apple~CloudDocs/"*) ;;
  *)
    echo "WARN: FLEXAIDDS_RESULTS was not on iCloud ($FLEXAIDDS_RESULTS) — resetting" >&2
    export FLEXAIDDS_RESULTS="$FLEXAIDDS_ICLOUD/results"
    ;;
esac

mkdir -p \
  "$FLEXAIDDS_ICLOUD" \
  "$FLEXAIDDS_RESULTS/campaigns" \
  "$FLEXAIDDS_RESULTS/working" \
  "$FLEXAIDDS_RESULTS/archive" \
  "$FLEXAIDDS_QUEUE_ROOT/logs" \
  "$FLEXAIDDS_QUEUE_ROOT/work" 2>/dev/null || true

echo "OK iCloud storage:"
echo "  FLEXAIDDS_ICLOUD=$FLEXAIDDS_ICLOUD"
echo "  FLEXAIDDS_RESULTS=$FLEXAIDDS_RESULTS"
echo "  FLEXAIDDS_QUEUE_ROOT=$FLEXAIDDS_QUEUE_ROOT"
