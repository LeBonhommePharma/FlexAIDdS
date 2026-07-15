#!/usr/bin/env bash
# sync_three_engine_local_to_icloud.sh — hang-safe local → iCloud mirror for A/B0/B OUT.
#
# Mirrors completed result.csv + thin metadata only (default).
# Never recursive find on CloudDocs. rsync with timeout when available.
#
# Usage:
#   bash scripts/sync_three_engine_local_to_icloud.sh
#   bash scripts/sync_three_engine_local_to_icloud.sh --with-poses
#   bash scripts/sync_three_engine_local_to_icloud.sh --arm A --campaign 3dsig_r10
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_local_first_benchmark_storage.sh"

WITH_POSES=0
DRY=0
ARM_FILTER=""
CAMPAIGN="${FLEXAID_SYNC_CAMPAIGN:-3dsig_r10}"
for a in "$@"; do
  case "$a" in
    --with-poses) WITH_POSES=1 ;;
    --dry-run) DRY=1 ;;
    --arm) shift; ARM_FILTER="${1:-}"; shift || true ;;
    --campaign) shift; CAMPAIGN="${1:-3dsig_r10}"; shift || true ;;
  esac
done

LOCAL_BASE="$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine"
REMOTE_BASE="$FLEXAIDDS_ICLOUD_RESULTS/campaigns/three_engine"
LOGDIR="$FLEXAIDDS_LOCAL_LOGDIR"
mkdir -p "$LOGDIR" 2>/dev/null || true
LOG="$LOGDIR/sync_three_engine_local_to_icloud.log"

log() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" | tee -a "$LOG"; }

RSYNC=(rsync -a --protect-args)
if command -v gtimeout >/dev/null 2>&1; then TO=(gtimeout 120)
elif command -v timeout >/dev/null 2>&1; then TO=(timeout 120)
else TO=(); fi

run_rsync() {
  if (( DRY )); then
    log "DRY: rsync $*"
    return 0
  fi
  if ((${#TO[@]})); then
    "${TO[@]}" "${RSYNC[@]}" "$@" 2>>"$LOG" || log "WARN: rsync exit=$? for $*"
  else
    "${RSYNC[@]}" "$@" 2>>"$LOG" || log "WARN: rsync failed for $*"
  fi
}

ARMS=(A B0 B)
if [[ -n "$ARM_FILTER" ]]; then
  ARMS=("$ARM_FILTER")
fi

log "START campaign=$CAMPAIGN arms=${ARMS[*]} with_poses=$WITH_POSES"
total=0
for arm in "${ARMS[@]}"; do
  LOCAL="$LOCAL_BASE/$arm/$CAMPAIGN"
  REMOTE="$REMOTE_BASE/$arm/$CAMPAIGN"
  if [[ ! -d "$LOCAL" ]]; then
    log "SKIP no local dir $LOCAL"
    continue
  fi
  mkdir -p "$REMOTE" 2>/dev/null || true
  for name in RUN_RECEIPT.json provenance.json; do
    [[ -f "$LOCAL/$name" ]] && run_rsync "$LOCAL/$name" "$REMOTE/$name"
  done
  for rc in "$LOCAL"/*/result.csv; do
    [[ -f "$rc" ]] || continue
    tdir=$(dirname "$rc")
    tid=$(basename "$tdir")
    case "$tid" in *_incomplete_*|backup_*|.*) continue ;; esac
    mkdir -p "$REMOTE/$tid" 2>/dev/null || true
    run_rsync "$rc" "$REMOTE/$tid/result.csv"
    for extra in validator_provenance.json RUN_RECEIPT.json meta.json elected_pose.pdb; do
      [[ -f "$tdir/$extra" ]] && run_rsync "$tdir/$extra" "$REMOTE/$tid/$extra"
    done
    if (( WITH_POSES )); then
      # only small ranked heads if present (no recursive chrom dumps)
      for p in "$tdir"/*.pdb; do
        [[ -f "$p" ]] || continue
        sz=$(stat -f%z "$p" 2>/dev/null || echo 0)
        # skip huge files
        if (( sz > 0 && sz < 5000000 )); then
          run_rsync "$p" "$REMOTE/$tid/$(basename "$p")"
        fi
      done
    fi
    total=$((total + 1))
    log "synced $arm/$tid"
  done
done
log "DONE synced_targets=$total"
echo "synced_targets=$total log=$LOG"
