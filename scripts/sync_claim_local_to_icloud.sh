#!/usr/bin/env bash
# sync_claim_local_to_icloud.sh — thin, hang-safe local → iCloud mirror for claim OUT.
#
# Only copies:
#   - */result.csv (completed targets)
#   - RUN_RECEIPT.json, provenance.json
#   - elected* / summary* if present (small)
# Never recursive pose dumps by default (use --with-poses).
# Uses rsync with timeouts; never `find` over CloudDocs trees.
#
# Usage:
#   bash scripts/sync_claim_local_to_icloud.sh [--with-poses] [--dry-run]
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/claim_local_staging_paths.sh"

WITH_POSES=0
DRY=0
for a in "$@"; do
  case "$a" in
    --with-poses) WITH_POSES=1 ;;
    --dry-run) DRY=1 ;;
  esac
done

LOCAL="${C0_CLAIM_LOCAL_OUT}"
REMOTE="${C0_CLAIM_ICLOUD_OUT}"
STAMP_DIR="${C0_CLAIM_LOCAL_LOGDIR}"
mkdir -p "$STAMP_DIR" "$REMOTE" 2>/dev/null || true

if [[ ! -d "$LOCAL" ]]; then
  echo "NO_LOCAL: $LOCAL"
  exit 0
fi

ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
echo "=== sync_claim_local_to_icloud $ts ==="
echo "LOCAL=$LOCAL"
echo "REMOTE=$REMOTE"

# Count completed via glob only (no find)
n=0
for f in "$LOCAL"/*/result.csv; do
  [[ -f "$f" ]] || continue
  n=$((n + 1))
done
echo "local_result_csv=$n"

RSYNC=(rsync -a --protect-args)
# Prefer not to hang forever if iCloud stalls
if command -v gtimeout >/dev/null 2>&1; then
  TO=(gtimeout 120)
elif command -v timeout >/dev/null 2>&1; then
  TO=(timeout 120)
else
  TO=()
fi

run_rsync() {
  if (( DRY )); then
    echo "DRY: ${RSYNC[*]} $*"
    return 0
  fi
  if ((${#TO[@]})); then
    "${TO[@]}" "${RSYNC[@]}" "$@" || {
      ec=$?
      echo "WARN: rsync exit=$ec (timeout or error) for: $*" >&2
      return 0
    }
  else
    "${RSYNC[@]}" "$@" || {
      echo "WARN: rsync failed for: $*" >&2
      return 0
    }
  fi
}

# Root thin files
for name in RUN_RECEIPT.json provenance.json; do
  if [[ -f "$LOCAL/$name" ]]; then
    run_rsync "$LOCAL/$name" "$REMOTE/$name"
  fi
done

# Per-target result.csv (+ optional thin extras)
synced=0
for rc in "$LOCAL"/*/result.csv; do
  [[ -f "$rc" ]] || continue
  tdir=$(dirname "$rc")
  tid=$(basename "$tdir")
  # Skip incomplete archives
  case "$tid" in
    *_incomplete_*|backup_*|.*) continue ;;
  esac
  mkdir -p "$REMOTE/$tid" 2>/dev/null || true
  run_rsync "$rc" "$REMOTE/$tid/result.csv"
  # thin optional files
  for extra in summary.csv elected.pdb elected_pose.pdb S1.pdb result.json; do
    if [[ -f "$tdir/$extra" ]]; then
      run_rsync "$tdir/$extra" "$REMOTE/$tid/$extra"
    fi
  done
  if (( WITH_POSES )); then
    # poses only if explicitly requested; still timeout-bound
    run_rsync --include='*/' --include='*.pdb' --include='*.mcf' --exclude='*' \
      "$tdir/" "$REMOTE/$tid/"
  fi
  synced=$((synced + 1))
done

# Status stamp (local + thin iCloud)
status_local="$STAMP_DIR/last_sync_status.txt"
{
  echo "utc=$ts"
  echo "local_result_csv=$n"
  echo "synced_targets=$synced"
  echo "local=$LOCAL"
  echo "remote=$REMOTE"
  echo "with_poses=$WITH_POSES"
} >"$status_local"

if (( ! DRY )); then
  # Tiny write to iCloud (single file, not a tree walk)
  if ((${#TO[@]})); then
    "${TO[@]}" cp "$status_local" "$REMOTE/LAST_LOCAL_SYNC.txt" 2>/dev/null || \
      echo "WARN: could not write LAST_LOCAL_SYNC.txt to iCloud" >&2
  else
    cp "$status_local" "$REMOTE/LAST_LOCAL_SYNC.txt" 2>/dev/null || true
  fi
fi

echo "OK synced_targets=$synced"
exit 0
