#!/usr/bin/env bash
# local_control_plane_probe.sh — prove agent control plane is not blocked by iCloud.
#
# Runs ONLY on local APFS paths. Never walks Mobile Documents / CloudDocs.
# Use after FileProvider storms or when cells hang at the 600s spawn cap.
#
# Usage:
#   bash scripts/local_control_plane_probe.sh
#   bash scripts/local_control_plane_probe.sh --out /path/to/log
#
# Exit 0 if local git + $HOME/.claude listing finish under WALL_SEC (default 45).
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

WALL_SEC="${WALL_SEC:-45}"
OUT="${1:-}"
if [[ "${1:-}" == "--out" ]]; then
  OUT="${2:?}"
fi

REPO="$(cd "$(dirname "$0")/.." && pwd)"
run() {
  python3 - "$WALL_SEC" "$@" <<'PY'
import subprocess, sys, time
wall = float(sys.argv[1])
cmd = sys.argv[2:]
t0 = time.time()
try:
    r = subprocess.run(cmd, timeout=wall, capture_output=True, text=True)
    dt = time.time() - t0
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)
    print(f"[probe] exit={r.returncode} wall_s={dt:.3f} cmd={' '.join(cmd)}", file=sys.stderr)
    sys.exit(0 if r.returncode == 0 and dt < wall else 1)
except subprocess.TimeoutExpired:
    print(f"[probe] TIMEOUT {wall}s: {' '.join(cmd)}", file=sys.stderr)
    sys.exit(124)
PY
}

log() { echo "$*"; }

{
  log "=== local_control_plane_probe $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  log "repo=$REPO wall=${WALL_SEC}s"
  log "--- df ---"
  df -h "$HOME" | tail -1
  log "--- .claude local? ---"
  if [[ -L "$HOME/.claude" ]]; then
    log "WARN: $HOME/.claude is a symlink -> $(readlink "$HOME/.claude")"
  else
    log "OK: $HOME/.claude is not a symlink"
  fi
  run /bin/ls "$HOME/.claude" >/dev/null
  log "--- git log -1 ---"
  run /usr/bin/git -C "$REPO" log -1 --oneline
  log "--- fileproviderd (informational) ---"
  ps aux 2>/dev/null | grep -E '[f]ileproviderd|[b]ird ' | head -5 || true
  log "=== PASS: local control plane responsive under ${WALL_SEC}s ==="
} 2>&1 | if [[ -n "$OUT" ]]; then tee "$OUT"; else cat; fi
