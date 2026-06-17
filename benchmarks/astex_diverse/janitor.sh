#!/usr/bin/env bash
# Disk janitor for the Astex Diverse benchmark run.
# - Caps per-complex stderr.log growth (memory: a v7 run once produced 50 MB/complex
#   and nearly filled the disk). Truncates any stderr.log > 5 MB to its last 200 lines.
# - Monitors free disk; if it drops below the hard floor, kills the benchmark to
#   protect the system, then exits.
set -u

CACHE_DIR="$HOME/.flexaidds/benchmarks/astex_diverse"
JLOG="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/janitor.log"
HARD_FLOOR_MB=1200          # kill benchmark if free space drops below this
STDERR_CAP_BYTES=$((5*1024*1024))
PIDFILE="/tmp/astex_bench.pid"

echo "[janitor] started $(date)" >> "$JLOG"

while true; do
  # 1. Cap oversized stderr.log files
  if [[ -d "$CACHE_DIR" ]]; then
    while IFS= read -r f; do
      [[ -z "$f" ]] && continue
      tail -n 200 "$f" > "$f.tmp" 2>/dev/null && mv "$f.tmp" "$f" 2>/dev/null
      echo "[janitor] $(date '+%H:%M:%S') capped $f" >> "$JLOG"
    done < <(find "$CACHE_DIR" -name 'stderr.log' -size +${STDERR_CAP_BYTES}c 2>/dev/null)
  fi

  # 2. Disk guard
  FREE_MB=$(df -m "$HOME" | tail -1 | awk '{print $4}')
  if [[ -n "$FREE_MB" && "$FREE_MB" -lt "$HARD_FLOOR_MB" ]]; then
    echo "[janitor] $(date '+%H:%M:%S') CRITICAL free=${FREE_MB}MB < ${HARD_FLOOR_MB}MB — killing benchmark" >> "$JLOG"
    if [[ -f "$PIDFILE" ]]; then
      kill "$(cat "$PIDFILE")" 2>/dev/null
      pkill -P "$(cat "$PIDFILE")" 2>/dev/null
    fi
    pkill -f 'benchmark_datasets --benchmark astex' 2>/dev/null
    echo "[janitor] exiting after emergency stop" >> "$JLOG"
    exit 1
  fi

  sleep 30
done
