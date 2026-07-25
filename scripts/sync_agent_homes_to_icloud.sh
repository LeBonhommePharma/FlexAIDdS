#!/usr/bin/env bash
# sync_agent_homes_to_icloud.sh — hang-safe local↔iCloud backup/restore for agent homes.
#
# PRIMARY (save to iCloud): backup mode — local agent homes → CloudDocs mirror.
# SECONDARY (restore): restore mode — iCloud archive/mirror → local homes
#   (matches the seed rsync direction for Claude from archive_batch layouts).
#
# Agents covered:
#   claude          ~/.claude                         → …/agent_homes/dot_claude
#   claude_app      ~/Library/Application Support/Claude → …/Application_Support_Claude
#   claude_science  ~/.claude-science (selective)     → …/dot_claude_science
#   codex           ~/.codex                          → …/dot_codex
#   grok            ~/.grok                           → …/dot_grok
#
# Safety contract (AGENTS.md local-first / thin-iCloud):
#   - Never replace live agent homes with CloudDocs symlinks.
#   - Default backup is thin (excludes caches, conda, vm_bundles, large runtimes).
#   - rsync is timeout-wrapped; never `find` / rglob over Mobile Documents.
#   - Does not delete local sources. Restore never uses --delete by default.
#
# Usage:
#   bash scripts/sync_agent_homes_to_icloud.sh --dry-run
#   bash scripts/sync_agent_homes_to_icloud.sh --backup
#   bash scripts/sync_agent_homes_to_icloud.sh --backup --agents claude,codex,grok
#   bash scripts/sync_agent_homes_to_icloud.sh --restore --archive-batch \
#     "$FLEXAIDDS_ICLOUD/archived_from_ssd/archive_batch_20260725T095624Z"
#   bash scripts/sync_agent_homes_to_icloud.sh --print-restore-cmds --archive-batch …
#   bash scripts/sync_agent_homes_to_icloud.sh --print-seed-restore --archive-batch …
#
# Env:
#   FLEXAIDDS_ICLOUD   durable CloudDocs root (default …/CloudDocs/FlexAIDdS_benchmarks)
#   FLEXAIDDS_ROOT     repo root (auto-detected)
#   AGENT_SYNC_TIMEOUT seconds per rsync (default 180)
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"

PATHS_PY="$ROOT/scripts/agent_icloud_paths.py"
MODE="backup"
DRY=0
FULL=0
AGENTS="all"
ARCHIVE_BATCH=""
DEST_ROOT=""
PRINT_RESTORE=0
PRINT_SEED=0
PRINT_MAP=0
TIMEOUT_SEC="${AGENT_SYNC_TIMEOUT:-180}"

usage() {
  sed -n '2,40p' "$0" | sed 's/^# \?//'
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --backup) MODE="backup"; shift ;;
    --restore) MODE="restore"; shift ;;
    --dry-run) DRY=1; shift ;;
    --full) FULL=1; shift ;;
    --agents) AGENTS="${2:-all}"; shift 2 ;;
    --archive-batch) ARCHIVE_BATCH="${2:-}"; shift 2 ;;
    --dest-root) DEST_ROOT="${2:-}"; shift 2 ;;
    --print-restore-cmds) PRINT_RESTORE=1; shift ;;
    --print-seed-restore) PRINT_SEED=1; shift ;;
    --print-map) PRINT_MAP=1; shift ;;
    -h|--help) usage ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

if [[ ! -f "$PATHS_PY" ]]; then
  echo "ERROR: missing $PATHS_PY" >&2
  exit 1
fi

# Resolve dest root
if [[ -n "$ARCHIVE_BATCH" ]]; then
  DEST_ROOT="$ARCHIVE_BATCH"
elif [[ -z "$DEST_ROOT" ]]; then
  DEST_ROOT="${FLEXAIDDS_ICLOUD}/agent_homes"
fi

py_args=(python3 "$PATHS_PY" --icloud "$FLEXAIDDS_ICLOUD" --agents "$AGENTS")
# Restore matches seed rsync -a (full archive → local). Thin excludes apply
# only to backup so large reinstallable trees are not pushed to iCloud.
if (( FULL )) || [[ "$MODE" == "restore" ]]; then
  py_args+=(--full)
fi
if [[ -n "$ARCHIVE_BATCH" ]]; then
  py_args+=(--archive-batch "$ARCHIVE_BATCH")
else
  py_args+=(--dest-root "$DEST_ROOT")
fi

if (( PRINT_SEED )); then
  if [[ -z "$ARCHIVE_BATCH" ]]; then
    # Default to the known seed batch if present, else print with placeholder name
    _seed_default="${FLEXAIDDS_ICLOUD}/archived_from_ssd/archive_batch_20260725T095624Z"
    ARCHIVE_BATCH="$_seed_default"
    py_args=(python3 "$PATHS_PY" --icloud "$FLEXAIDDS_ICLOUD" --archive-batch "$ARCHIVE_BATCH")
  fi
  echo "=== seed Claude restore (iCloud archive → local homes) ==="
  echo "LOCAL agent homes stay real directories (never CloudDocs symlinks)."
  "${py_args[@]}" --print-seed-restore
  exit 0
fi

if (( PRINT_RESTORE )); then
  if [[ -z "$ARCHIVE_BATCH" && -z "$DEST_ROOT" ]]; then
    echo "ERROR: --print-restore-cmds needs --archive-batch or --dest-root" >&2
    exit 2
  fi
  echo "=== restore rsync commands (iCloud → local) ==="
  echo "mode=restore dest_root=$DEST_ROOT"
  "${py_args[@]}" --print-restore-cmds
  exit 0
fi

if (( PRINT_MAP )); then
  "${py_args[@]}" --print-map --mode "$MODE"
  exit 0
fi

ts=$(date -u +%Y-%m-%dT%H:%M:%SZ)
echo "=== sync_agent_homes_to_icloud $ts ==="
echo "mode=$MODE dry_run=$DRY full=$FULL agents=$AGENTS"
echo "FLEXAIDDS_ICLOUD=$FLEXAIDDS_ICLOUD"
echo "DEST_ROOT=$DEST_ROOT"
echo "timeout_sec=$TIMEOUT_SEC"

# Timeout wrapper (same pattern as sync_claim_local_to_icloud.sh)
if command -v gtimeout >/dev/null 2>&1; then
  TO=(gtimeout "$TIMEOUT_SEC")
elif command -v timeout >/dev/null 2>&1; then
  TO=(timeout "$TIMEOUT_SEC")
else
  TO=()
fi

RSYNC=(rsync -a --protect-args)
# Never delete local or remote by default (safe mirror).
# Operator can re-run with manual rsync if they need --delete.

map_json=$("${py_args[@]}" --print-map --mode "$MODE")
echo "--- path map ---"
echo "$map_json" | python3 -c '
import json, sys
data = json.load(sys.stdin)
print("dest_root=", data.get("dest_root"))
print("icloud=", data.get("icloud"))
for row in data.get("pairs", []):
    print("  [{id}] {src}  ->  {dst}".format(
        id=row["agent_id"], src=row["source"], dst=row["dest"]
    ))
    for ex in row.get("excludes") or []:
        print("      exclude: {0}".format(ex))
'

# Guard: refuse to treat live homes as CloudDocs symlinks (backup sources)
check_not_cloud_symlink() {
  local p="$1"
  [[ -L "$p" ]] || return 0
  local t
  t=$(readlink "$p" 2>/dev/null || true)
  case "$t" in
    *"Mobile Documents"*|*"com~apple~CloudDocs"*)
      echo "ERROR: refusing CloudDocs symlink live home: $p -> $t" >&2
      echo "  Live agent homes must stay on local APFS (AGENTS.md)." >&2
      return 1
      ;;
  esac
  return 0
}

run_one_rsync() {
  local src="$1" dst="$2"
  shift 2 || true
  local rsync_cmd=("${RSYNC[@]}")
  # Iterate remaining args as exclude patterns; empty $@ is fine under set -u.
  local ex
  while [[ $# -gt 0 ]]; do
    ex="$1"
    shift
    [[ -n "$ex" ]] && rsync_cmd+=(--exclude="$ex")
  done

  if (( DRY )); then
    echo "DRY: ${rsync_cmd[*]} ${src}/ ${dst}/"
    return 0
  fi

  mkdir -p "$dst" 2>/dev/null || {
    echo "WARN: could not mkdir $dst (iCloud may be stalled)" >&2
    return 0
  }

  if ((${#TO[@]})); then
    "${TO[@]}" "${rsync_cmd[@]}" "${src}/" "${dst}/" || {
      ec=$?
      echo "WARN: rsync exit=$ec (timeout or error) for $src -> $dst" >&2
      return 0
    }
  else
    "${rsync_cmd[@]}" "${src}/" "${dst}/" || {
      echo "WARN: rsync failed for $src -> $dst" >&2
      return 0
    }
  fi
  echo "OK: $src -> $dst"
}

# Parse pairs and run rsync without deep CloudDocs walks
synced=0
skipped=0
while IFS=$'\t' read -r agent_id source dest excludes_csv; do
  [[ -n "$agent_id" ]] || continue

  if [[ "$MODE" == "backup" ]]; then
    check_not_cloud_symlink "$source" || exit 3
    if [[ ! -d "$source" ]]; then
      echo "SKIP missing source: $source ($agent_id)"
      skipped=$((skipped + 1))
      continue
    fi
  else
    # restore: source is under iCloud; dest is local home
    check_not_cloud_symlink "$dest" || exit 3
    if [[ ! -d "$source" ]]; then
      echo "SKIP missing archive dir: $source ($agent_id)"
      skipped=$((skipped + 1))
      continue
    fi
  fi

  # Split excludes on | (empty for full/restore — set -u safe)
  local_excludes=()
  if [[ -n "$excludes_csv" ]]; then
    IFS='|' read -r -a local_excludes <<<"$excludes_csv"
  fi

  if ((${#local_excludes[@]})); then
    run_one_rsync "$source" "$dest" "${local_excludes[@]}"
  else
    run_one_rsync "$source" "$dest"
  fi
  synced=$((synced + 1))
done < <(
  echo "$map_json" | python3 -c '
import json, sys
data = json.load(sys.stdin)
for row in data.get("pairs", []):
    excl = "|".join(row.get("excludes") or [])
    print("\t".join([row["agent_id"], row["source"], row["dest"], excl]))
'
)

# Local stamp only (never deep write trees on CloudDocs for status)
LOCAL_STAMP_DIR="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/logs/agent_homes_sync"
mkdir -p "$LOCAL_STAMP_DIR" 2>/dev/null || true
stamp="$LOCAL_STAMP_DIR/last_sync_status.txt"
{
  echo "utc=$ts"
  echo "mode=$MODE"
  echo "dry_run=$DRY"
  echo "full=$FULL"
  echo "agents=$AGENTS"
  echo "dest_root=$DEST_ROOT"
  echo "synced_pairs=$synced"
  echo "skipped_pairs=$skipped"
  echo "icloud=$FLEXAIDDS_ICLOUD"
} >"$stamp" 2>/dev/null || true

echo "OK synced_pairs=$synced skipped_pairs=$skipped stamp=$stamp"
echo "NOTE: live agent homes remain local directories; iCloud is durable mirror only."
exit 0
