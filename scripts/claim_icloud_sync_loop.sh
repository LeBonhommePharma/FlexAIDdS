#!/usr/bin/env bash
# claim_icloud_sync_loop.sh — hang-safe ops + local→iCloud sync every N minutes.
#
# Intermediate iCloud hang mitigation:
#   1) Live claim I/O is on local disk (see run_C0_claim_local.sh).
#   2) This loop periodically mirrors completed result.csv to iCloud.
#   3) Writes hang-safe status without recursive CloudDocs walks / find.
#   4) Detects stall (0% CPU + stale log) and ALERTs; never dual-launches.
#
# Usage:
#   bash scripts/claim_icloud_sync_loop.sh              # foreground loop
#   bash scripts/claim_icloud_sync_loop.sh --once       # single tick
#   nohup bash scripts/claim_icloud_sync_loop.sh &      # background
#
# Env:
#   C0_CLAIM_SYNC_EVERY_MIN  (default 5)
#   C0_CLAIM_STALL_MIN       (default 20)
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/claim_local_staging_paths.sh"

ONCE=0
[[ "${1:-}" == "--once" ]] && ONCE=1

EVERY_MIN="${C0_CLAIM_SYNC_EVERY_MIN:-5}"
STALL_MIN="${C0_CLAIM_STALL_MIN:-20}"
LOGDIR="${C0_CLAIM_LOCAL_LOGDIR}"
mkdir -p "$LOGDIR"
LOOP_LOG="$LOGDIR/sync_ops_loop.log"
PIDF="$LOGDIR/sync_ops_loop.pid"
STATUS_MD="$LOGDIR/OPS_STATUS_local.md"
# Thin mirror on iCloud (single file — never a tree)
ICLOUD_STATUS="${C0_CLAIM_ICLOUD_OUT}/OPS_STATUS_agent.md"
# Also keep a copy under iCloud queue logs if writable
QLOG_STATUS="${FLEXAIDDS_QUEUE_ROOT}/logs/benchmark_ops_latest.md"

echo $$ >"$PIDF"

log() {
  local line="[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
  echo "$line" | tee -a "$LOOP_LOG"
}

# Resolve claim parent PID from local or iCloud pid files (no pgrep -f self-match)
claim_pid() {
  local f p
  for f in \
    "$LOGDIR/C0_claim_clean.pid" \
    "${FLEXAIDDS_QUEUE_ROOT}/logs/C0_claim_clean.pid" \
    "${FLEXAIDDS_QUEUE_ROOT}/logs/C0_full85.pid"
  do
    if [[ -f "$f" ]]; then
      p=$(tr -d ' \n' <"$f" 2>/dev/null || true)
      if [[ -n "$p" ]] && kill -0 "$p" 2>/dev/null; then
        echo "$p"
        return 0
      fi
    fi
  done
  # Fallback: scan process table for local OUT path only
  local local_out="${C0_CLAIM_LOCAL_OUT}"
  ps -ax -o pid=,command= 2>/dev/null | while read -r pid cmd; do
    case "$cmd" in
      *benchmark_datasets*"$local_out"*)
        echo "$pid"
        return 0
        ;;
    esac
  done
  return 1
}

flexaidds_child() {
  ps -ax -o pid=,%cpu=,rss=,etime=,command= 2>/dev/null | while read -r line; do
    case "$line" in
      *bin/C/FlexAIDdS*|*bin/C/FlexAID\ *|*FlexAIDdS\ *)
        # skip self
        case "$line" in *claim_icloud_sync_loop*) continue ;; esac
        echo "$line"
        ;;
    esac
  done | head -5
}

count_result_csv() {
  local root="$1" n=0 f
  [[ -d "$root" ]] || { echo 0; return; }
  for f in "$root"/*/result.csv; do
    [[ -f "$f" ]] || continue
    n=$((n + 1))
  done
  echo "$n"
}

log_mtime_age_min() {
  local logf="$1"
  if [[ ! -f "$logf" ]]; then
    echo 99999
    return
  fi
  local mt now
  mt=$(stat -f %m "$logf" 2>/dev/null || echo 0)
  now=$(date +%s)
  echo $(( (now - mt) / 60 ))
}

cpu_of() {
  local pid="$1"
  ps -p "$pid" -o %cpu= 2>/dev/null | tr -d ' ' || echo 0
}

tick() {
  local ts pid cpu age n_local n_cloud logf child_lines alert
  ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  logf="${LOGDIR}/C0_claim_clean.log"
  [[ -f "$logf" ]] || logf="${FLEXAIDDS_QUEUE_ROOT}/logs/C0_claim_clean.log"

  n_local=$(count_result_csv "${C0_CLAIM_LOCAL_OUT}")
  # Cloud count can hang if iCloud is wedged — bound with soft failure
  n_cloud="?"
  if command -v timeout >/dev/null 2>&1; then
    n_cloud=$(timeout 30 bash -c "source '$ROOT/scripts/claim_local_staging_paths.sh'; c=0; for f in \"\$C0_CLAIM_ICLOUD_OUT\"/*/result.csv; do [[ -f \"\$f\" ]] && c=\$((c+1)); done; echo \$c" 2>/dev/null || echo "?")
  else
    n_cloud=$(count_result_csv "${C0_CLAIM_ICLOUD_OUT}" 2>/dev/null || echo "?")
  fi

  pid=$(claim_pid || true)
  cpu="n/a"
  age=$(log_mtime_age_min "$logf")
  alert=""
  child_lines=$(flexaidds_child || true)

  if [[ -n "${pid:-}" ]]; then
    cpu=$(cpu_of "$pid")
    # Stall heuristic: live, ~0 CPU, log stale
    if awk "BEGIN{exit !($cpu < 0.5)}"; then
      if (( age >= STALL_MIN )); then
        alert="STALL: pid=$pid cpu=$cpu log_age_min=$age (>=${STALL_MIN})"
      fi
    fi
  else
    alert="NO_LIVE_CLAIM"
  fi

  # Sync local → iCloud (completed only)
  bash "$ROOT/scripts/sync_claim_local_to_icloud.sh" >>"$LOOP_LOG" 2>&1 || true

  # Hang-safe status markdown (local first)
  {
    echo "# Claim ops (local staging) — $ts"
    echo
    echo "- **local_OUT**: \`${C0_CLAIM_LOCAL_OUT}\`"
    echo "- **icloud_mirror**: \`${C0_CLAIM_ICLOUD_OUT}\`"
    echo "- **N_local**: ${n_local}/85"
    echo "- **N_icloud**: ${n_cloud}/85"
    echo "- **claim_pid**: ${pid:-none}  cpu=${cpu}  log_age_min=${age}"
    echo "- **sync_every_min**: ${EVERY_MIN}"
    echo "- **stall_threshold_min**: ${STALL_MIN}"
    if [[ -n "$alert" ]]; then
      echo
      echo "**ALERT:** $alert"
    fi
    echo
    echo "## FlexAIDdS / runner (ps snapshot)"
    echo '```'
    if [[ -n "${pid:-}" ]]; then
      ps -p "$pid" -o pid=,%cpu=,rss=,etime=,command= 2>/dev/null || true
    fi
    echo "${child_lines:-"(no FlexAIDdS child)"}"
    echo '```'
    echo
    echo "## Log tail (local preferred)"
    echo '```'
    if [[ -f "$logf" ]]; then
      tail -20 "$logf" 2>/dev/null || true
    else
      echo "(no log)"
    fi
    echo '```'
    echo
    echo "CF = scoring proxy a.u. Shannon election G̃=H̃−T·S̃ on DatasetRunner."
    echo "Live I/O is local; iCloud is periodic mirror only (intermediate hang fix)."
  } >"$STATUS_MD"

  # Thin push of status to iCloud (timeout-bound, single file)
  if command -v timeout >/dev/null 2>&1; then
    timeout 20 cp "$STATUS_MD" "$ICLOUD_STATUS" 2>/dev/null || true
    timeout 20 cp "$STATUS_MD" "$QLOG_STATUS" 2>/dev/null || true
  else
    cp "$STATUS_MD" "$ICLOUD_STATUS" 2>/dev/null || true
    cp "$STATUS_MD" "$QLOG_STATUS" 2>/dev/null || true
  fi

  if [[ -n "$alert" ]]; then
    log "ALERT $alert N_local=$n_local"
  else
    log "OK pid=${pid:-none} cpu=$cpu N_local=$n_local N_icloud=$n_cloud log_age=$age"
  fi
}

log "START loop every=${EVERY_MIN}m stall=${STALL_MIN}m once=$ONCE"
tick
if (( ONCE )); then
  log "ONCE done"
  exit 0
fi

while true; do
  sleep $(( EVERY_MIN * 60 ))
  tick || log "tick error (continuing)"
done
