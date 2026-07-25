#!/usr/bin/env bash
# sync_home_dots_to_icloud.sh — hang-safe local→iCloud mirror for allowlisted ~/.*
#
# Mirrors meaningful home-dot configs (not agent trees) under:
#   $FLEXAIDDS_ICLOUD/home_dots/
# Named missing paths (.venv, .env) are SKIP-logged, never invented.
# Ephemeral noise (.zcompdump*, histories) is not uploaded.
#
# Agent homes (.claude/.codex/.grok) use sync_agent_homes_to_icloud.sh.
#
# Usage:
#   bash scripts/sync_home_dots_to_icloud.sh --dry-run
#   bash scripts/sync_home_dots_to_icloud.sh
#   bash scripts/sync_home_dots_to_icloud.sh --print-inventory
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"
PATHS_PY="$ROOT/scripts/agent_icloud_paths.py"
DRY=0
PRINT_INV=0
TIMEOUT_SEC="${HOME_DOT_SYNC_TIMEOUT:-120}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --print-inventory) PRINT_INV=1; shift ;;
    -h|--help)
      sed -n '2,20p' "$0" | sed 's/^# \?//'
      exit 0
      ;;
    *) echo "Unknown: $1" >&2; exit 2 ;;
  esac
done

if [[ ! -f "$PATHS_PY" ]]; then
  echo "ERROR: missing $PATHS_PY" >&2
  exit 1
fi

DEST_ROOT="${FLEXAIDDS_ICLOUD}/home_dots"
ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
echo "=== sync_home_dots_to_icloud $ts ==="
echo "DEST_ROOT=$DEST_ROOT dry_run=$DRY"

if (( PRINT_INV )); then
  python3 "$PATHS_PY" --print-inventory --icloud "$FLEXAIDDS_ICLOUD"
  exit 0
fi

# Log named MISSING (.venv/.env) explicitly
python3 "$PATHS_PY" --print-inventory --icloud "$FLEXAIDDS_ICLOUD" | python3 -c '
import json, sys
data = json.load(sys.stdin)
for it in data.get("items", []):
    if it["name"] in (".venv", ".env", ".claude", ".codex", ".grok"):
        print("NAMED {status}: {name} — {reason}".format(**it))
'

map_json=$(python3 "$PATHS_PY" --print-home-dots-map --icloud "$FLEXAIDDS_ICLOUD" --dest-root "$DEST_ROOT")
echo "--- home_dots map ---"
echo "$map_json" | python3 -c '
import json, sys
data = json.load(sys.stdin)
print("dest_root=", data.get("dest_root"))
for row in data.get("pairs", []):
    print("  [{id}] {src} -> {dst} exists={ex}".format(
        id=row["agent_id"], src=row["source"], dst=row["dest"],
        ex=row.get("source_exists")))
'

if command -v gtimeout >/dev/null 2>&1; then TO=(gtimeout "$TIMEOUT_SEC")
elif command -v timeout >/dev/null 2>&1; then TO=(timeout "$TIMEOUT_SEC")
else TO=(); fi

RSYNC=(rsync -a --protect-args)
synced=0
skipped=0

while IFS=$'\t' read -r agent_id source dest kind; do
  [[ -n "$agent_id" ]] || continue
  if [[ ! -e "$source" ]]; then
    echo "SKIP missing: $source"
    skipped=$((skipped + 1))
    continue
  fi
  if (( DRY )); then
    if [[ -d "$source" ]]; then
      echo "DRY: rsync -a ${source}/ ${dest}/"
    else
      echo "DRY: rsync -a ${source} ${dest}"
    fi
    synced=$((synced + 1))
    continue
  fi
  mkdir -p "$(dirname "$dest")" 2>/dev/null || {
    echo "WARN: mkdir failed for $(dirname "$dest")" >&2
    continue
  }
  if [[ -d "$source" ]]; then
    mkdir -p "$dest" 2>/dev/null || true
    cmd=("${RSYNC[@]}" "${source}/" "${dest}/")
  else
    # file: ensure parent exists; dest is file path under home_dots
    cmd=("${RSYNC[@]}" "$source" "$dest")
  fi
  if ((${#TO[@]})); then
    "${TO[@]}" "${cmd[@]}" || {
      echo "WARN: rsync exit=$? for $source -> $dest" >&2
      continue
    }
  else
    "${cmd[@]}" || {
      echo "WARN: rsync failed for $source -> $dest" >&2
      continue
    }
  fi
  echo "OK: $source -> $dest"
  synced=$((synced + 1))
done < <(
  echo "$map_json" | python3 -c '
import json, sys
from pathlib import Path
data = json.load(sys.stdin)
for row in data.get("pairs", []):
    src = Path(row["source"])
    kind = "dir" if src.is_dir() else "file"
    print("\t".join([row["agent_id"], row["source"], row["dest"], kind]))
'
)

stamp_dir="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/logs/home_dots_sync"
mkdir -p "$stamp_dir" 2>/dev/null || true
{
  echo "utc=$ts"
  echo "synced=$synced"
  echo "skipped=$skipped"
  echo "dest_root=$DEST_ROOT"
} >"$stamp_dir/last_sync_status.txt" 2>/dev/null || true

echo "OK synced_home_dots=$synced skipped=$skipped"
exit 0
