#!/bin/bash
# Close the launcher's check-then-act window so two shard workers cannot enter
# the same target directory.
#
# WHY THIS IS NEEDED
#   launch_wall_campaign.sh line 94 gates on `[ -f "$d/result.csv" ]`, but the
#   engine writes its result one level deeper at $d/$T/result.csv. The gate has
#   therefore NEVER fired. With one worker per arm that was invisible; with
#   overlapping shards it lets a second worker re-enter a finished target.
#
# WHAT IT DOES, AND WHAT IT REFUSES TO DO
#   - a target dir with no $d/result.csv and no engine result  -> write a CLAIM
#     sentinel, so a second worker's gate fires and it skips.
#   - a target dir with an engine result but no gate file       -> hardlink the
#     engine's own file up to the gate path (same bytes, not a re-derivation).
#   - a FINISHED cell still holding a CLAIM sentinel            -> replace the
#     sentinel with the engine's real result.
#   It never overwrites a real result, never deletes anything, and never
#   fabricates a receipt for a cell that has not finished.
#
# THE HAZARD IT INTRODUCES, AND THE DETECTOR FOR IT
#   A leftover CLAIM makes an unfinished target look complete. Detect with:
#   CELL_PROVENANCE.json absent while $d/result.csv starts with "CLAIM,".
#   shard_collision_ledger.py reports exactly that condition.
set -u
O="${1:?usage: claim_daemon.sh <arm_dir> <seconds>}"
DUR="${2:-34200}"
LOG="$(dirname "$O")/claim_daemon.log"
echo "$(date -u +%FT%TZ) DAEMON-START dir=$O duration=${DUR}s" >> "$LOG"
END=$(( $(date +%s) + DUR ))
while [ "$(date +%s)" -lt "$END" ]; do
  for d in "$O"/*/; do
    [ -d "$d" ] || continue
    t=$(basename "$d"); [ ${#t} -eq 4 ] || continue
    if [ ! -f "$d/result.csv" ]; then
      real=$(find "$d" -mindepth 2 -name result.csv 2>/dev/null | head -1)
      if [ -n "${real:-}" ]; then
        ln "$real" "$d/result.csv" 2>/dev/null || cp "$real" "$d/result.csv"
        echo "$(date -u +%FT%TZ) PROMOTE $t" >> "$LOG"
      else
        printf 'CLAIM,%s,%s,claim_daemon\n' "$t" "$(date -u +%FT%TZ)" > "$d/result.csv"
        echo "$(date -u +%FT%TZ) CLAIM $t" >> "$LOG"
      fi
    elif head -1 "$d/result.csv" 2>/dev/null | grep -q '^CLAIM,'; then
      if [ -f "$d/CELL_PROVENANCE.json" ]; then
        real=$(find "$d" -mindepth 2 -name result.csv 2>/dev/null | head -1)
        if [ -n "${real:-}" ]; then
          cp "$real" "$d/result.csv"
          echo "$(date -u +%FT%TZ) REPLACE-SENTINEL $t" >> "$LOG"
        fi
      fi
    fi
  done
  sleep 3
done
echo "$(date -u +%FT%TZ) DAEMON-EXIT" >> "$LOG"
