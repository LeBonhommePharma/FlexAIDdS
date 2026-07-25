#!/usr/bin/env bash
# safe_free_verified_icloud_duplicates.sh — free local bulk ONLY after iCloud verify.
#
# Safety contract:
#   1) Never free whole live agent roots (~/.claude, ~/.codex, ~/.grok).
#   2) Never replace homes with CloudDocs symlinks.
#   3) Never free without a proof line: remote mirror exists + non-empty (shallow).
#   4) No recursive find over Mobile Documents.
#   5) Default is --dry-run; pass --execute to actually delete freeable paths.
#
# Modes of free (both require agent_homes mirror verify for that agent):
#   --free-regenerable   free regenerable cache subdirs (thin backup excludes them;
#                        live config remains; parent agent mirror must be non-empty
#                        OR archive_batch already holds complete agent tree)
#   --free-flexaidds-dup free local flexaidds_results campaign/result trees only when
#                        the matching archive/iCloud path has result.csv (one-level)
#
# Usage:
#   bash scripts/safe_free_verified_icloud_duplicates.sh --dry-run --free-regenerable
#   bash scripts/safe_free_verified_icloud_duplicates.sh --execute --free-regenerable
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"
export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
PATHS_PY="$ROOT/scripts/agent_icloud_paths.py"

DRY=1
FREE_REGEN=0
FREE_FLEX=0
AGENTS="claude,codex,grok,claude_app"
ARCHIVE_BATCH="${FLEXAIDDS_ICLOUD}/archived_from_ssd/archive_batch_20260725T095624Z"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --execute) DRY=0; shift ;;
    --free-regenerable) FREE_REGEN=1; shift ;;
    --free-flexaidds-dup) FREE_FLEX=1; shift ;;
    --agents) AGENTS="${2:-$AGENTS}"; shift 2 ;;
    --archive-batch) ARCHIVE_BATCH="${2:-}"; shift 2 ;;
    -h|--help)
      sed -n '2,28p' "$0" | sed 's/^# \?//'
      exit 0
      ;;
    *) echo "Unknown: $1" >&2; exit 2 ;;
  esac
done

if (( ! FREE_REGEN && ! FREE_FLEX )); then
  FREE_REGEN=1
fi

ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
echo "=== safe_free_verified_icloud_duplicates $ts ==="
echo "dry_run=$DRY free_regenerable=$FREE_REGEN free_flexaidds_dup=$FREE_FLEX"
echo "agent_homes=${FLEXAIDDS_ICLOUD}/agent_homes"
echo "archive_batch=$ARCHIVE_BATCH"

AGENT_ROOT="${FLEXAIDDS_ICLOUD}/agent_homes"

# Shallow verify: remote dir exists and has ≥1 top-level entry
verify_remote_nonempty() {
  local remote="$1"
  [[ -d "$remote" ]] || return 1
  local n
  n=$(ls -1A "$remote" 2>/dev/null | wc -l | tr -d ' ')
  [[ "${n:-0}" -ge 1 ]]
}

# Map agent_id -> remote_name via python
remote_name_for() {
  local aid="$1"
  python3 "$PATHS_PY" --print-map --agents "$aid" --icloud "$FLEXAIDDS_ICLOUD" 2>/dev/null \
    | python3 -c 'import json,sys; d=json.load(sys.stdin); print(d["pairs"][0]["dest"].split("/")[-1] if d.get("pairs") else "")' \
    2>/dev/null || true
}

agent_mirror_ok() {
  local aid="$1"
  local rname remote archive_remote
  rname=$(remote_name_for "$aid")
  [[ -n "$rname" ]] || return 1
  remote="$AGENT_ROOT/$rname"
  archive_remote="$ARCHIVE_BATCH/$rname"
  if verify_remote_nonempty "$remote"; then
    echo "VERIFY_OK agent=$aid remote=$remote"
    return 0
  fi
  # Archive batch complete copy counts as durable proof (claude case)
  if verify_remote_nonempty "$archive_remote"; then
    echo "VERIFY_OK agent=$aid archive=$archive_remote"
    return 0
  fi
  echo "VERIFY_FAIL agent=$aid remote=$remote archive=$archive_remote"
  return 1
}

freed=0
skipped=0
bytes_note=""

if (( FREE_REGEN )); then
  echo "--- regenerable free candidates ---"
  free_json=$(python3 "$PATHS_PY" --print-freeable --agents "$AGENTS" --icloud "$FLEXAIDDS_ICLOUD")
  while IFS=$'\t' read -r agent_id relative local agent_root remote_name; do
    [[ -n "$agent_id" ]] || continue
    # Never free the agent root itself
    if [[ "$local" == "$agent_root" || "$local" == "$agent_root/" ]]; then
      echo "REFUSE free whole agent root: $local"
      skipped=$((skipped + 1))
      continue
    fi
    if ! agent_mirror_ok "$agent_id"; then
      echo "SKIP no durable proof for $agent_id — not freeing $local"
      skipped=$((skipped + 1))
      continue
    fi
    if [[ ! -e "$local" ]]; then
      echo "SKIP absent: $local"
      skipped=$((skipped + 1))
      continue
    fi
    # Refuse if path is CloudDocs symlink
    if [[ -L "$local" ]]; then
      t=$(readlink "$local" 2>/dev/null || true)
      case "$t" in
        *"Mobile Documents"*|*"com~apple~CloudDocs"*)
          echo "REFUSE CloudDocs symlink: $local"
          skipped=$((skipped + 1))
          continue
          ;;
      esac
    fi
    proof="agent_homes_or_archive non-empty for $agent_id (remote_name=$remote_name)"
    if (( DRY )); then
      echo "WOULD_FREE local=$local proof=$proof"
      freed=$((freed + 1))
    else
      echo "FREE local=$local proof=$proof"
      rm -rf "$local"
      echo "FREED local=$local"
      freed=$((freed + 1))
    fi
  done < <(
    echo "$free_json" | python3 -c '
import json, sys
data = json.load(sys.stdin)
for row in data.get("freeable", []):
    print("\t".join([
        row["agent_id"], row["relative"], row["local"],
        row["agent_local_root"], row["remote_name"],
    ]))
'
  )
fi

if (( FREE_FLEX )); then
  echo "--- flexaidds_results local dups (one-level campaigns) ---"
  LOCAL_CAMP="${FLEXAIDDS_LOCAL_ROOT}/campaigns"
  # Prefer archive flexaidds_results and iCloud results/campaigns
  IC_CAMP="${FLEXAIDDS_ICLOUD}/results/campaigns"
  ARCH_RES="${ARCHIVE_BATCH}/flexaidds_results"
  if [[ -d "$LOCAL_CAMP" ]]; then
    for d in "$LOCAL_CAMP"/*; do
      [[ -d "$d" ]] || continue
      name=$(basename "$d")
      # Protect skeleton / active names lightly — only free if remote has result.csv one-level
      remote=""
      if [[ -f "$IC_CAMP/$name/result.csv" ]] || [[ -d "$IC_CAMP/$name" ]]; then
        # one-level: any target result.csv
        has=0
        for rc in "$IC_CAMP/$name"/*/result.csv; do
          [[ -f "$rc" ]] && has=1 && break
        done
        [[ -f "$IC_CAMP/$name/result.csv" ]] && has=1
        if (( has )); then
          remote="$IC_CAMP/$name"
        fi
      fi
      if [[ -z "$remote" && -d "$ARCH_RES" ]]; then
        # archive layout may nest campaigns differently — only free if exact name dir exists with entries
        if [[ -d "$ARCH_RES/$name" ]]; then
          n=$(ls -1A "$ARCH_RES/$name" 2>/dev/null | wc -l | tr -d ' ')
          if [[ "${n:-0}" -ge 1 ]]; then
            remote="$ARCH_RES/$name"
          fi
        fi
      fi
      if [[ -z "$remote" ]]; then
        echo "SKIP no iCloud proof for local campaign $name"
        skipped=$((skipped + 1))
        continue
      fi
      if (( DRY )); then
        echo "WOULD_FREE local=$d proof=$remote"
        freed=$((freed + 1))
      else
        echo "FREE local=$d proof=$remote"
        # Prefer freeing large pose trees if present, not entire campaign if active
        # Only free completed target subdirs that have result.csv on both sides
        for td in "$d"/*; do
          [[ -d "$td" ]] || continue
          tid=$(basename "$td")
          if [[ -f "$td/result.csv" && -f "$remote/$tid/result.csv" ]]; then
            # Free bulky pose dumps only, keep result.csv local for resume
            for bulk in poses dock_config elected_pool chroms; do
              if [[ -d "$td/$bulk" ]]; then
                echo "FREE bulk $td/$bulk proof=$remote/$tid/result.csv"
                rm -rf "$td/$bulk"
              fi
            done
          fi
        done
        freed=$((freed + 1))
      fi
    done
  else
    echo "SKIP no local campaigns dir $LOCAL_CAMP"
  fi
fi

# Post-check: live agent roots still local dirs
echo "--- live agent roots after free ---"
for p in "$HOME/.claude" "$HOME/.codex" "$HOME/.grok"; do
  if [[ -L "$p" ]]; then
    echo "FAIL symlink $p -> $(readlink "$p")"
  elif [[ -d "$p" ]]; then
    echo "OK local_dir $p"
  else
    echo "NOTE missing $p"
  fi
done

echo "OK free_actions=$freed skipped=$skipped dry_run=$DRY"
exit 0
