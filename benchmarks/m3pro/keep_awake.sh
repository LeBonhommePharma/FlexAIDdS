#!/usr/bin/env bash
# keep_awake.sh — Persistent macOS blocker for long benchmarks
# Prevents display sleep, system sleep, idle sleep, disk sleep.
# Run once at start of benchmark campaign.
# Safe: just kills the caffeinate when done or on logout.
#
# Usage:
#   chmod +x benchmarks/m3pro/keep_awake.sh
#   nohup ./benchmarks/m3pro/keep_awake.sh &
#   # or
#   ./benchmarks/m3pro/keep_awake.sh start
#   ./benchmarks/m3pro/keep_awake.sh stop
#   ./benchmarks/m3pro/keep_awake.sh status

set -euo pipefail
PIDFILE=/tmp/caffeinate_bench_master.pid
LOG=/tmp/caffeinate_master.log

start() {
  if [ -f "$PIDFILE" ] && ps -p "$(cat "$PIDFILE")" > /dev/null 2>&1; then
    echo "Already running: $(cat $PIDFILE)"
    ps -p "$(cat $PIDFILE)" -o pid,etime,command
    return 0
  fi
  echo "Starting persistent caffeinate -dimsu ..."
  nohup caffeinate -dimsu >> "$LOG" 2>&1 &
  echo $! > "$PIDFILE"
  sleep 1
  if ps -p "$(cat "$PIDFILE")" > /dev/null 2>&1; then
    echo "Started. pid=$(cat $PIDFILE)"
    echo "Log: $LOG"
    pmset -g assertions | grep -E 'PreventUserIdle|PreventSystem' | head -3 || true
  else
    echo "Failed to start caffeinate" >&2
    exit 1
  fi
}

stop() {
  if [ -f "$PIDFILE" ]; then
    PID=$(cat "$PIDFILE")
    if ps -p "$PID" > /dev/null 2>&1; then
      kill "$PID" && echo "Stopped caffeinate $PID"
    else
      echo "No process for $PID"
    fi
    rm -f "$PIDFILE"
  else
    echo "No pidfile"
  fi
  # Also clean any lingering plain caffeinate if you want (dangerous if others use)
  # pkill -f 'caffeinate -dimsu' 2>/dev/null || true
}

status() {
  if [ -f "$PIDFILE" ] && ps -p "$(cat "$PIDFILE")" > /dev/null 2>&1; then
    echo "Master running:"
    ps -p "$(cat "$PIDFILE")" -o pid,etime,command
  else
    echo "Master NOT running"
  fi
  echo ""
  echo "Current sleep state:"
  pmset -g | grep -E 'sleep |displaysleep '
  echo ""
  echo "Active caffeinate assertions:"
  pmset -g assertions | grep -c caffeinate || echo 0
}

case "${1:-start}" in
  start) start ;;
  stop)  stop ;;
  status) status ;;
  *) echo "Usage: $0 [start|stop|status]"; exit 1 ;;
esac
