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
#   - a nested result whose own pdb_id is NOT this target                -> REFUSE
#     to promote it, and log REFUSE-PROMOTE. `find | head -1` picked the first
#     nested result.csv with no id check, so a run written into another
#     target's directory would be promoted into the wrong authoritative slot.
#     Observed: campaigns/wall_paired_85/{c1_on,c1_off}/2HR7/ each hold a row
#     whose pdb_id is 2J62 (CELL_PROVENANCE.json records target=2J62 while
#     RUN_RECEIPT.json records output=.../2HR7). 2HR7's own run timed out at
#     7200 s with zero poses; 2J62's scored 1.2996 A. Promoting by path there
#     converts a timeout into a sub-2 A success.
#
# THE HAZARD IT INTRODUCES, AND THE DETECTOR FOR IT
#   A leftover CLAIM makes an unfinished target look complete. Detect with:
#   CELL_PROVENANCE.json absent while $d/result.csv starts with "CLAIM,".
#   shard_collision_ledger.py reports exactly that condition.
set -u
O="${1:?usage: claim_daemon.sh <arm_dir> <seconds>}"
DUR="${2:-34200}"
LOG="$(dirname "$O")/claim_daemon.log"

# First nested result.csv whose own pdb_id equals $1 (the target dir name).
# Empty if none matches -- the caller must then refuse to promote.
matching_nested_result() {
  local want p id
  want=$(echo "$1" | tr '[:lower:]' '[:upper:]')
  while IFS= read -r p; do
    id=$(awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i=="pdb_id") c=i; next}
                  NR==2 && c{gsub(/^[ \t"]+|[ \t"]+$/,"",$c); print toupper($c); exit}' "$p" 2>/dev/null)
    if [ -n "$id" ] && [ "$id" = "$want" ]; then printf '%s\n' "$p"; return 0; fi
    [ -n "$id" ] && echo "$(date -u +%FT%TZ) REFUSE-PROMOTE $1: $p carries pdb_id=$id" >> "$LOG"
  done < <(find "$2" -mindepth 2 -name result.csv 2>/dev/null)
  return 1
}

echo "$(date -u +%FT%TZ) DAEMON-START dir=$O duration=${DUR}s" >> "$LOG"
END=$(( $(date +%s) + DUR ))
while [ "$(date +%s)" -lt "$END" ]; do
  for d in "$O"/*/; do
    [ -d "$d" ] || continue
    t=$(basename "$d"); [ ${#t} -eq 4 ] || continue
    if [ ! -f "$d/result.csv" ]; then
      real=$(matching_nested_result "$t" "$d" || true)
      if [ -n "${real:-}" ]; then
        ln "$real" "$d/result.csv" 2>/dev/null || cp "$real" "$d/result.csv"
        echo "$(date -u +%FT%TZ) PROMOTE $t" >> "$LOG"
      else
        printf 'CLAIM,%s,%s,claim_daemon\n' "$t" "$(date -u +%FT%TZ)" > "$d/result.csv"
        echo "$(date -u +%FT%TZ) CLAIM $t" >> "$LOG"
      fi
    elif head -1 "$d/result.csv" 2>/dev/null | grep -q '^CLAIM,'; then
      if [ -f "$d/CELL_PROVENANCE.json" ]; then
        real=$(matching_nested_result "$t" "$d" || true)
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
