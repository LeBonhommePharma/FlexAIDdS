#!/bin/bash
# FlexAIDdS OPS monitor — run this on the Mac alongside external coding agents.
# Emits a machine-readable log the OPS agent (Claude Science) reads to audit.
# Usage:
#   ops/monitor.sh [interval_s] [logfile] [pattern]
#     interval_s : sample period (default 30)
#     logfile    : output (default /tmp/flexaidds_ops_monitor.log)
#     pattern    : egrep of process comm/args to track (default docking+agents)
set -u
INT="${1:-30}"
LOG="${2:-/tmp/flexaidds_ops_monitor.log}"
PAT="${3:-FlexAIDdS|benchmark_datasets|claude|codex|grok|node .*claude-code}"
REPO="/Users/lp.more/Projects/FlexAIDdS"
echo "# FlexAIDdS OPS monitor started $(date -u +%FT%TZ) interval=${INT}s pattern='${PAT}'" >> "$LOG"
while true; do
  TS=$(date -u +%FT%TZ)
  # --- machine resources ---
  LOAD=$(sysctl -n vm.loadavg 2>/dev/null | tr -d '{}'); [ -z "$LOAD" ] && LOAD=$(uptime 2>/dev/null | sed -E 's/.*load averages?: //')
  # free + inactive pages -> GiB (page size 16384 on Apple Silicon)
  PS=$(vm_stat 2>/dev/null)
  PGSIZE=$(echo "$PS" | awk '/page size of/{gsub(/[^0-9]/,"",$0);print}')
  PGSIZE=${PGSIZE:-16384}
  FREE=$(echo "$PS" | awk '/Pages free/{gsub(/\./,"",$3);print $3}')
  INACT=$(echo "$PS" | awk '/Pages inactive/{gsub(/\./,"",$3);print $3}')
  AVAIL_GB=$(awk -v f="${FREE:-0}" -v i="${INACT:-0}" -v p="$PGSIZE" 'BEGIN{printf "%.2f",(f+i)*p/1073741824}')
  echo "RES ${TS} load=${LOAD} avail_gb=${AVAIL_GB}" >> "$LOG"
  # --- tracked processes (PID RSS_MB %CPU ELAPSED CMD) ---
  ps -axo pid,rss,pcpu,etime,comm,args 2>/dev/null \
    | grep -E "$PAT" | grep -v 'grep -E' \
    | while read -r pid rss pcpu etime comm rest; do
        rssmb=$(awk -v r="$rss" 'BEGIN{printf "%.0f", r/1024}')
        echo "PROC ${TS} pid=${pid} rss_mb=${rssmb} cpu=${pcpu} elapsed=${etime} comm=${comm}" >> "$LOG"
      done
  # --- repo HEAD + dirty state (so OPS sees what each agent committed) ---
  HEAD=$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null)
  BR=$(git -C "$REPO" rev-parse --abbrev-ref HEAD 2>/dev/null)
  DIRTY=$(git -C "$REPO" status --porcelain 2>/dev/null | grep -vE '^\?\?' | wc -l | tr -d ' ')
  echo "GIT ${TS} branch=${BR} head=${HEAD} tracked_dirty=${DIRTY}" >> "$LOG"
  echo "---" >> "$LOG"
  sleep "$INT"
done
