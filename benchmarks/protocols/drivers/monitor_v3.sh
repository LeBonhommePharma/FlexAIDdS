#!/usr/bin/env bash
# Campaign monitor v3. Every wake condition EXITS (a background cell only
# notifies on exit -- v1's P2 branch logged and looped, so it never woke).
# v2 died with exit -1 at 20:17Z with no output, one minute after a foreground
# stop_child on an unrelated cell; v3 lives in a file so the loop is not tied to
# an in-kernel heredoc.
set -u
R=/Users/lp.more/flexaidds_results; Q=$(cat "$R/state/requeue_dir")
LOG=$R/state/monitor_final.log
deadline=$(( $(date +%s) + 12*3600 )); last=-1
echo "MONITOR v3 ARMED $(date -u +%FT%TZ) batch=$(basename "$Q")" >> "$LOG"
while [ "$(date +%s)" -lt "$deadline" ]; do
  n=$(( $(find "$Q" -name DONE 2>/dev/null | wc -l | tr -d ' ') - 2 ))
  p2=$(find "$Q/P2_nosec_flex" -name DONE 2>/dev/null | wc -l | tr -d ' ')
  p3=$(find "$Q/P3_inducedfit" -name DONE 2>/dev/null | wc -l | tr -d ' ')
  g_rc=$(find "$Q" -name DONE -exec grep -l 'rc=[1-9]' {} + 2>/dev/null | wc -l | tr -d ' ')
  g_op=$(find "$Q" -name DONE -exec grep -l 'optres_ok=0 ' {} + 2>/dev/null | wc -l | tr -d ' ')
  g_pz=$(find "$Q" -name DONE -exec grep -l 'poses=0 ' {} + 2>/dev/null | wc -l | tr -d ' ')
  g_sv=$(find "$Q" -type f -newermt '-15 minutes' -print0 2>/dev/null | xargs -0 grep -laiE 'exit 13[0-9]|segmentation|bus error|SIGBUS|SIGSEGV|SIGABRT' 2>/dev/null | wc -l | tr -d ' ')
  fresh=$(find "$Q" -type f -newermt '-600 seconds' 2>/dev/null | wc -l | tr -d ' ')
  load=$(uptime | sed 's/.*load average[s]*: //' | awk '{print $1}')
  disk=$(df -m /System/Volumes/Data | awk 'NR==2{printf "%.1f", $4/1024}')
  line="$(date -u +%H:%MZ) cells=$n/345 P2=$p2/84 P3=$p3/9 guards rc=$g_rc optres=$g_op poses0=$g_pz sig15m=$g_sv fresh600s=$fresh load=$load disk=${disk}G"
  echo "$line" >> "$LOG"
  if [ "$g_rc" -gt 0 ] || [ "$g_op" -gt 0 ] || [ "$g_pz" -gt 0 ] || [ "$g_sv" -gt 0 ]; then echo "WAKE=GUARD_TRIPPED $line"; exit 0; fi
  if [ "$(df -m /System/Volumes/Data | awk 'NR==2{print int($4/1024)}')" -lt 10 ]; then echo "WAKE=DISK_LOW $line"; exit 0; fi
  if [ "$fresh" -eq 0 ] && [ "$n" -eq "$last" ] && [ ! -f "$Q/DONE" ]; then echo "WAKE=STALL $line"; exit 0; fi
  if [ "$p2" -ge 84 ] && [ -z "${P2_ACK:-}" ]; then echo "WAKE=P2_COMPLETE $line"; exit 0; fi
  if [ -f "$Q/DONE" ]; then
    echo "WAKE=BATCH_DONE $(cat "$Q/DONE")   $line"
    echo "  waiting 600s for the chains to fire..."; sleep 600
    for l in chain_water chain_queue3 chain_recovery; do p=$R/state/$l.log
      printf "  %-16s age %5ss  last: %s\n" "$l" "$(( $(date +%s) - $(stat -f %m "$p" 2>/dev/null || echo 0) ))" "$(tail -1 "$p" 2>/dev/null | cut -c1-80)"; done
    echo "  new batch dirs (15 min): $(find "$R" -maxdepth 1 -type d -newermt '-15 minutes' 2>/dev/null | xargs -n1 basename 2>/dev/null | grep -vE '^(state|cache.*|flexaidds_results)$' | tr '\n' ' ')"
    echo "  box_queue: $(ls "$R/state/box_queue/" 2>/dev/null | tr '\n' ' ')"
    exit 0; fi
  last=$n; sleep 600
done
echo "WAKE=DEADLINE $(tail -1 "$LOG")"
