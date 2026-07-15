#!/usr/bin/env bash
# reap_hung_icloud_walkers.sh — kill ONLY stuck CloudDocs walkers; never dockers.
#
# Targets hung agent/ops processes that block on iCloud FileProvider:
#   find / mdfind / rg / python …rglob… under Mobile Documents
#
# NEVER kills:
#   FlexAIDdS, FlexAID, benchmark_datasets, caffeinate claim workers, dock PIDs
#
# Usage:
#   bash scripts/reap_hung_icloud_walkers.sh           # report + reap
#   bash scripts/reap_hung_icloud_walkers.sh --dry-run # report only
#   MIN_AGE_SEC=120 bash scripts/reap_hung_icloud_walkers.sh
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1
MIN_AGE_SEC="${MIN_AGE_SEC:-90}"

# etime is [[dd-]hh:]mm:ss — rough age in seconds
etime_to_sec() {
  local e="$1" d=0 h=0 m=0 s=0
  if [[ "$e" == *-* ]]; then
    d="${e%%-*}"
    e="${e#*-}"
  fi
  IFS=: read -r a b c <<<"$e"
  if [[ -n "${c:-}" ]]; then
    h=$a; m=$b; s=$c
  else
    m=$a; s=$b
  fi
  echo $((10#$d * 86400 + 10#$h * 3600 + 10#$m * 60 + 10#$s))
}

is_protected() {
  local cmd="$1"
  case "$cmd" in
    *bin/C/FlexAIDdS*|*bin/C/FlexAID*|*bin/A/FlexAID*|*bin/B/FlexAID*) return 0 ;;
    *benchmark_datasets*) return 0 ;;
    *caffeinate*FlexAIDdS*|*caffeinate*benchmark_datasets*) return 0 ;;
    *run_C0_claim*|*claim_icloud_sync_loop*) return 0 ;;
    *cmake*|*ninja*|*xcodebuild*) return 0 ;;
  esac
  return 1
}

is_icloud_walker() {
  local cmd="$1"
  case "$cmd" in
    *Mobile\ Documents*|*com~apple~CloudDocs*)
      case "$cmd" in
        *find\ *|*mdfind*|*rglob*|*Path.rglob*|*os.walk*|*bfs\ *|*/usr/bin/find*|*rg\ -*CloudDocs*|*python*CloudDocs*)
          return 0
          ;;
      esac
      # long-running md5/shasum/openssl on CloudDocs
      case "$cmd" in
        *md5*|*openssl\ md5*|*shasum*|*cmp\ *)
          return 0
          ;;
      esac
      ;;
  esac
  return 1
}

reaped=0
skipped=0
echo "=== reap_hung_icloud_walkers dry=$DRY min_age=${MIN_AGE_SEC}s ==="

# pid, etime, command
while IFS= read -r line; do
  pid=$(echo "$line" | awk '{print $1}')
  etime=$(echo "$line" | awk '{print $2}')
  cmd=$(echo "$line" | sed -E 's/^[[:space:]]*[0-9]+[[:space:]]+[^[:space:]]+[[:space:]]+//')
  [[ -n "$pid" ]] || continue
  if is_protected "$cmd"; then
    skipped=$((skipped + 1))
    continue
  fi
  if ! is_icloud_walker "$cmd"; then
    continue
  fi
  age=$(etime_to_sec "$etime" 2>/dev/null || echo 0)
  if (( age < MIN_AGE_SEC )); then
    echo "SKIP young age=${age}s pid=$pid ${cmd:0:100}"
    continue
  fi
  echo "REAP age=${age}s pid=$pid ${cmd:0:120}"
  if (( ! DRY )); then
    kill -TERM "$pid" 2>/dev/null || true
    sleep 0.5
    kill -KILL "$pid" 2>/dev/null || true
  fi
  reaped=$((reaped + 1))
done < <(ps -ax -o pid=,etime=,command= 2>/dev/null)

echo "reaped=$reaped protected_seen=$skipped"
exit 0
