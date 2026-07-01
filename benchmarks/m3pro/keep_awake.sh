#!/usr/bin/env bash
# keep_awake.sh — Persistent macOS blocker for long benchmarks
# Prevents display sleep, system sleep, idle sleep, disk sleep + screen lock / logout.
# Uses:
#   - caffeinate -dimsu : core assertions (Prevent*Idle*Sleep, PreventSystemSleep)
#   - user pinger loop : repeated short -u assertions to keep "UserIsActive" high
#     This defeats lock screen timers and idle-logout even when askForPassword=0.
# Run once at start of benchmark campaign.
# Safe: just kills the caffeinates when done (or on logout / reboot they die).
#
# Usage:
#   chmod +x benchmarks/m3pro/keep_awake.sh
#   ./benchmarks/m3pro/keep_awake.sh start     # ensure master + pinger
#   ./benchmarks/m3pro/keep_awake.sh stop
#   ./benchmarks/m3pro/keep_awake.sh status
#
# Also safe to run from nohup at boot of a long campaign.

set -euo pipefail
PIDFILE=/tmp/caffeinate_bench_master.pid
LOG=/tmp/caffeinate_master.log
PINGER_PIDFILE=/tmp/caffeinate_user_pinger.pid
PINGER_LOG=/tmp/caffeinate_pinger.log

start_pinger() {
  if [ -f "$PINGER_PIDFILE" ] && ps -p "$(cat "$PINGER_PIDFILE")" > /dev/null 2>&1; then
    echo "User pinger already: $(cat "$PINGER_PIDFILE")"
    return 0
  fi
  echo "Starting user-activity pinger (repeated caffeinate -u to hold UserIsActive) ..."
  nohup bash -c '
    while true; do
      caffeinate -u -t 30 >/dev/null 2>&1 || true
      sleep 25
    done
  ' >> "$PINGER_LOG" 2>&1 &
  echo $! > "$PINGER_PIDFILE"
  sleep 0.5
  if ps -p "$(cat "$PINGER_PIDFILE")" > /dev/null 2>&1; then
    echo "Pinger pid=$(cat "$PINGER_PIDFILE")  Log: $PINGER_LOG"
  else
    echo "WARNING: pinger failed to start"
  fi
}

start() {
  ensured=0

  # Master -dimsu (core sleep prevention)
  if [ -f "$PIDFILE" ] && ps -p "$(cat "$PIDFILE")" > /dev/null 2>&1; then
    echo "Master already running: $(cat "$PIDFILE")  ($(ps -p "$(cat "$PIDFILE")" -o etime= 2>/dev/null))"
  else
    echo "Starting persistent caffeinate -dimsu (display/system/disk) ..."
    nohup caffeinate -dimsu >> "$LOG" 2>&1 &
    echo $! > "$PIDFILE"
    ensured=1
  fi

  # Pinger: keeps UserIsActive refreshed → blocks screen lock / auto-logout
  start_pinger

  # Extra: lock out screensaver + password prompt (prevents "locking off")
  echo "Enforcing screensaver/lock prevention (idleTime=0, no password prompt)..."
  defaults -currentHost write com.apple.screensaver idleTime -int 0 2>/dev/null || true
  defaults -currentHost write com.apple.screensaver askForPassword -int 0 2>/dev/null || true
  defaults -currentHost write com.apple.screensaver askForPasswordDelay -int 0 2>/dev/null || true
  killall ScreenSaverEngine 2>/dev/null || true

  # Note: for even stronger base (powernap/standby/sleep timers), run once in a shell with sudo:
  #   sudo pmset sleep 0 displaysleep 0 disksleep 0 powernap 0 standby 0

  if [ "$ensured" = "1" ]; then
    sleep 1
    if ps -p "$(cat "$PIDFILE")" > /dev/null 2>&1; then
      echo "Master started. pid=$(cat "$PIDFILE")  Log: $LOG"
      pmset -g assertions | grep -E 'PreventUserIdle|PreventSystem' | head -3 || true
    else
      echo "Failed to start master caffeinate" >&2
    fi
  fi
}

stop_pinger() {
  if [ -f "$PINGER_PIDFILE" ]; then
    PID=$(cat "$PINGER_PIDFILE")
    if ps -p "$PID" > /dev/null 2>&1; then
      kill "$PID" 2>/dev/null && echo "Stopped pinger $PID"
    fi
    rm -f "$PINGER_PIDFILE"
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
  stop_pinger
  # Also clean any lingering plain caffeinate if you want (dangerous if others use)
  # pkill -f 'caffeinate -dimsu' 2>/dev/null || true
}

status() {
  echo "=== Master (dimsu) ==="
  if [ -f "$PIDFILE" ] && ps -p "$(cat "$PIDFILE")" > /dev/null 2>&1; then
    ps -p "$(cat "$PIDFILE")" -o pid,etime,command
  else
    echo "Master NOT running"
  fi

  echo ""
  echo "=== User pinger ==="
  if [ -f "$PINGER_PIDFILE" ] && ps -p "$(cat "$PINGER_PIDFILE")" > /dev/null 2>&1; then
    ps -p "$(cat "$PINGER_PIDFILE")" -o pid,etime,command
  else
    echo "User pinger NOT running"
  fi

  echo ""
  echo "=== Current sleep/display state ==="
  pmset -g | grep -E 'sleep |displaysleep |disksleep '
  echo ""
  echo "=== Active caffeinate assertions ==="
  pmset -g assertions | grep -c caffeinate || echo 0
  echo ""
  echo "=== Key assertions (from pmset) ==="
  pmset -g assertions | grep -E 'PreventUserIdle|PreventSystem|UserIsActive' | head -6 || true
}

case "${1:-start}" in
  start) start ;;
  stop)  stop ;;
  status) status ;;
  *) echo "Usage: $0 [start|stop|status]"; exit 1 ;;
esac
