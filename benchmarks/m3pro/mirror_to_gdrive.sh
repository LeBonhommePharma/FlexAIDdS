#!/usr/bin/env bash
# mirror_to_gdrive.sh — DEPRECATED / REMOVED
#
# Google Drive mirror support has been removed per user request.
# All FlexAIDdS m3pro benchmarking now uses iCloud 2TB exclusively.
#
# This script is kept as a placeholder. It is now a safe no-op.
#
# Apache-2.0 (c) 2026 NRGlab, Universite de Montreal

echo "[INFO] Google Drive mirror support has been removed."
echo "       m3pro benchmarks use iCloud 2TB as the only storage."
echo "       (File kept for future reference; no action taken.)"
exit 0

# ─── Logging ─────────────────────────────────────────────────────────────────

TIMESTAMP="$(date +%Y%m%dT%H%M%S)"
LOG_DIR="$FLEXAIDDS_ICLOUD/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/mirror_${TIMESTAMP}.log"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG_FILE"; }

log "Mirror started: iCloud -> Google Drive"
log "Source:      $FLEXAIDDS_ICLOUD"
log "Destination: $FLEXAIDDS_GDRIVE"

# ─── Sync directories ───────────────────────────────────────────────────────

DIRS_TO_SYNC=("benchmark_data" "results" "logs")
FAILURES=0
TOTAL_BYTES=0

for dir in "${DIRS_TO_SYNC[@]}"; do
    SRC="$FLEXAIDDS_ICLOUD/$dir/"
    DST="$FLEXAIDDS_GDRIVE/$dir/"

    if [[ ! -d "$SRC" ]]; then
        log "SKIP: $dir/ (source does not exist)"
        continue
    fi

    mkdir -p "$DST"

    log "SYNC: $dir/ ..."
    SYNC_START=$(date +%s)

    # rsync with:
    #   -a  archive mode (preserves permissions, timestamps)
    #   -v  verbose (for logging)
    #   -z  compress during transfer (helps over cloud FS)
    #   --delete  remove files from dst that don't exist in src
    #   --stats   show transfer statistics
    if rsync_output=$(rsync -avz --delete --stats "$SRC" "$DST" 2>&1); then
        SYNC_END=$(date +%s)
        DURATION=$((SYNC_END - SYNC_START))

        # Extract bytes transferred from rsync stats
        bytes=$(echo "$rsync_output" | grep -o 'Total transferred file size: [0-9,]*' | grep -o '[0-9,]*' | tr -d ',' || echo "0")
        TOTAL_BYTES=$((TOTAL_BYTES + ${bytes:-0}))

        log "  OK: $dir/ synced in ${DURATION}s (${bytes:-0} bytes)"
    else
        FAILURES=$((FAILURES + 1))
        log "  FAIL: $dir/ rsync returned $?"
        echo "$rsync_output" >> "$LOG_FILE"
    fi
done

# ─── Summary ─────────────────────────────────────────────────────────────────

log ""
log "Mirror complete: ${#DIRS_TO_SYNC[@]} dirs attempted, $FAILURES failures"
log "Total bytes synced: $TOTAL_BYTES"
log "Log: $LOG_FILE"

if [[ $FAILURES -eq ${#DIRS_TO_SYNC[@]} ]]; then
    log "STATUS: TOTAL FAILURE"
    exit 2
elif [[ $FAILURES -gt 0 ]]; then
    log "STATUS: PARTIAL FAILURE ($FAILURES/${#DIRS_TO_SYNC[@]} dirs failed)"
    exit 1
else
    log "STATUS: SUCCESS"
    exit 0
fi
