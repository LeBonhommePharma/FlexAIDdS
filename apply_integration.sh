#!/usr/bin/env bash
# apply_integration.sh — FlexAIDdS CMA-ES integration on-ramp (CHUNK 3)
#
# Idempotent installer for the opt-in CMA-ES search backend wiring.
# Run from a clean tree (or re-run after partial apply); second run is a no-op.
#
# What it does:
#   1. Copies LIB/cmaes_search.{cpp,h} from SWARM_ARTIFACTS / relative chunk paths
#      (when present — does not invent the adapter; chunk1 owns that TU)
#   2. Inserts cmaes_search.cpp into LIB/CMakeLists.txt FLEXAID_CORE_SOURCES
#      immediately after gaboom.cpp (if missing)
#   3. Applies the FLEXAIDDS_SEARCH=cmaes branch in LIB/top.cpp between
#      // FLEXAIDDS_CMAES_BEGIN / END markers (if missing)
#      Branch matches chunk1 API: cmaes_run_dock / cmaes_fill_chromosomes /
#      cmaes_write_trace_csv (population, max_evals, write_trace, n_evals).
#   4. Ensures analysis/ exists (collapse_fingerprint.py lands via chunk4)
#   5. Installs CMAES_INTEGRATION.md to repo root when present next to this script
#   6. Smoke-compiles LIB/cmaes_search.cpp with proper -I paths when the TU exists
#
# Usage:
#   bash .swarm/cmaes/chunk3_onramp/artifacts/apply_integration.sh
#   SWARM_ARTIFACTS=/path/to/chunk1/artifacts bash apply_integration.sh
#   FLEXAIDDS_ROOT=/path/to/repo bash apply_integration.sh
#   APPLY_INTEGRATION_SKIP_SMOKE=1 bash apply_integration.sh
#
# Copyright 2026 Le Bonhomme Pharma / Louis-Philippe Morency / NRGlab
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

resolve_repo_root() {
  if [[ -n "${FLEXAIDDS_ROOT:-}" && -d "${FLEXAIDDS_ROOT}/LIB" ]]; then
    (cd "$FLEXAIDDS_ROOT" && pwd)
    return 0
  fi
  local d="$SCRIPT_DIR"
  while [[ "$d" != "/" ]]; do
    if [[ -f "$d/AGENTS.md" && -d "$d/LIB" && -f "$d/LIB/CMakeLists.txt" ]]; then
      (cd "$d" && pwd)
      return 0
    fi
    d="$(dirname "$d")"
  done
  if command -v git >/dev/null 2>&1; then
    local gt
    gt="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel 2>/dev/null || true)"
    if [[ -n "$gt" && -d "$gt/LIB" ]]; then
      (cd "$gt" && pwd)
      return 0
    fi
  fi
  return 1
}

REPO_ROOT="$(resolve_repo_root)" || {
  echo "FAIL: cannot resolve FlexAIDdS repo root (set FLEXAIDDS_ROOT)" >&2
  exit 1
}
export FLEXAIDDS_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0
SKIP_COUNT=0

log()  { printf '%s\n' "$*"; }
ok()   { PASS_COUNT=$((PASS_COUNT + 1)); log "  [OK]   $*"; }
warn() { WARN_COUNT=$((WARN_COUNT + 1)); log "  [WARN] $*"; }
skip() { SKIP_COUNT=$((SKIP_COUNT + 1)); log "  [SKIP] $*"; }
fail() { FAIL_COUNT=$((FAIL_COUNT + 1)); log "  [FAIL] $*" >&2; }

MARKER_BEGIN="// FLEXAIDDS_CMAES_BEGIN"
MARKER_END="// FLEXAIDDS_CMAES_END"
INCLUDE_BEGIN="// FLEXAIDDS_CMAES_INCLUDE_BEGIN"
INCLUDE_END="// FLEXAIDDS_CMAES_INCLUDE_END"
CMAKE_BEGIN="# FLEXAIDDS_CMAES_CMAKE_BEGIN"
CMAKE_END="# FLEXAIDDS_CMAES_CMAKE_END"

log "=== FlexAIDdS CMA-ES apply_integration.sh ==="
log "REPO_ROOT=$REPO_ROOT"
log "SCRIPT_DIR=$SCRIPT_DIR"

# ─── 1. Locate & copy adapter sources ────────────────────────────────────────
log ""
log "[1/5] Adapter sources (LIB/cmaes_search.{cpp,h})"

find_artifact() {
  local base="$1"
  local candidates=()

  if [[ -n "${SWARM_ARTIFACTS:-}" ]]; then
    candidates+=("$SWARM_ARTIFACTS/$base")
    candidates+=("$SWARM_ARTIFACTS/LIB/$base")
  fi

  candidates+=(
    "$REPO_ROOT/.swarm/cmaes/chunk1_adapter/artifacts/LIB/$base"
    "$REPO_ROOT/.swarm/cmaes/chunk1_adapter/artifacts/$base"
    "$REPO_ROOT/.swarm/cmaes/chunk2_wiring/artifacts/LIB/$base"
    "$REPO_ROOT/.swarm/cmaes/chunk2_wiring/artifacts/$base"
    "$SCRIPT_DIR/LIB/$base"
    "$SCRIPT_DIR/$base"
    "$REPO_ROOT/LIB/$base"
  )

  local c
  for c in "${candidates[@]}"; do
    if [[ -f "$c" ]]; then
      printf '%s\n' "$c"
      return 0
    fi
  done
  return 1
}

copy_if_needed() {
  local base="$1"
  local dest="$REPO_ROOT/LIB/$base"

  if [[ -f "$dest" ]]; then
    ok "already present: LIB/$base"
    return 0
  fi

  local src=""
  if src="$(find_artifact "$base")"; then
    cp -f "$src" "$dest"
    ok "copied LIB/$base  ←  $src"
    return 0
  fi

  warn "LIB/$base not found (set SWARM_ARTIFACTS or land chunk1 adapter first)"
  return 0
}

copy_if_needed "cmaes_search.h"
copy_if_needed "cmaes_search.cpp"

# ─── 2. CMake FLEXAID_CORE_SOURCES insertion after gaboom.cpp ────────────────
log ""
log "[2/5] LIB/CMakeLists.txt — FLEXAID_CORE_SOURCES"

CMAKE_FILE="$REPO_ROOT/LIB/CMakeLists.txt"
if [[ ! -f "$CMAKE_FILE" ]]; then
  fail "missing $CMAKE_FILE"
else
  if grep -qE '^\s*cmaes_search\.cpp\s*$' "$CMAKE_FILE" \
     || grep -qF "$CMAKE_BEGIN" "$CMAKE_FILE"; then
    ok "cmaes_search.cpp already listed in FLEXAID_CORE_SOURCES"
  else
    if ! grep -qE '^\s*gaboom\.cpp\s*$' "$CMAKE_FILE"; then
      fail "gaboom.cpp not found in $CMAKE_FILE — cannot insert after it"
    else
      python3 - "$CMAKE_FILE" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
out = []
inserted = False
block = (
    "# FLEXAIDDS_CMAES_CMAKE_BEGIN\n"
    "    cmaes_search.cpp\n"
    "# FLEXAIDDS_CMAES_CMAKE_END\n"
)
for line in lines:
    out.append(line)
    if (not inserted) and line.strip() == "gaboom.cpp":
        out.append(block)
        inserted = True
path.write_text("".join(out), encoding="utf-8")
sys.exit(0 if inserted else 2)
PY
      if grep -qE 'cmaes_search\.cpp' "$CMAKE_FILE"; then
        ok "inserted cmaes_search.cpp after gaboom.cpp in FLEXAID_CORE_SOURCES"
      else
        fail "failed to insert cmaes_search.cpp into LIB/CMakeLists.txt"
      fi
    fi
  fi
fi

# ─── 3. top.cpp FLEXAIDDS_SEARCH branch ──────────────────────────────────────
log ""
log "[3/5] LIB/top.cpp — FLEXAIDDS_SEARCH=cmaes branch"

TOP_FILE="$REPO_ROOT/LIB/top.cpp"
if [[ ! -f "$TOP_FILE" ]]; then
  fail "missing $TOP_FILE"
else
  if grep -qF "$INCLUDE_BEGIN" "$TOP_FILE"; then
    ok "cmaes_search.h include already present (markers)"
  else
    python3 - "$TOP_FILE" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
include_block = (
    "// FLEXAIDDS_CMAES_INCLUDE_BEGIN\n"
    "#include \"cmaes_search.h\"\n"
    "// FLEXAIDDS_CMAES_INCLUDE_END\n"
)
needle = '#include "GAContext.h"\n'
if needle in text:
    text = text.replace(needle, needle + include_block, 1)
else:
    needle2 = '#include "gaboom.h"\n'
    if needle2 in text:
        text = text.replace(needle2, needle2 + include_block, 1)
    else:
        text = include_block + text
path.write_text(text, encoding="utf-8")
PY
    if grep -qF "$INCLUDE_BEGIN" "$TOP_FILE"; then
      ok "added #include \"cmaes_search.h\" with markers"
    else
      fail "failed to insert cmaes_search.h include into top.cpp"
    fi
  fi

  if grep -qF "$MARKER_BEGIN" "$TOP_FILE"; then
    ok "FLEXAIDDS_SEARCH branch already present (markers)"
  else
    # Apply branch via external helper (avoids bash/heredoc quoting issues).
    # Helper lives next to this script for re-runs from clean trees.
    HELPER="$SCRIPT_DIR/_patch_top_cmaes.py"
    if [[ ! -f "$HELPER" ]]; then
      fail "missing $HELPER (chunk3 patch helper)"
    else
      if python3 "$HELPER" "$TOP_FILE"; then
        if grep -qF "$MARKER_BEGIN" "$TOP_FILE" && grep -qF "$MARKER_END" "$TOP_FILE"; then
          ok "applied FLEXAIDDS_SEARCH branch with FLEXAIDDS_CMAES_BEGIN/END markers"
        else
          fail "patch helper ran but markers missing in top.cpp"
        fi
      else
        fail "failed to apply FLEXAIDDS_SEARCH branch to top.cpp"
      fi
    fi
  fi
fi

# ─── 4. analysis/ directory ──────────────────────────────────────────────────
log ""
log "[4/5] analysis/ directory"

ANALYSIS_DIR="$REPO_ROOT/analysis"
if [[ -d "$ANALYSIS_DIR" ]]; then
  ok "analysis/ already exists"
else
  mkdir -p "$ANALYSIS_DIR"
  if [[ ! -f "$ANALYSIS_DIR/.gitkeep" ]]; then
    cat > "$ANALYSIS_DIR/.gitkeep" <<'KEEP'
# Placeholder so analysis/ is present before collapse_fingerprint.py lands (chunk4).
KEEP
  fi
  ok "created analysis/"
fi

DOC_SRC="$SCRIPT_DIR/CMAES_INTEGRATION.md"
DOC_DST="$REPO_ROOT/CMAES_INTEGRATION.md"
if [[ -f "$DOC_SRC" ]]; then
  if [[ -f "$DOC_DST" ]] && cmp -s "$DOC_SRC" "$DOC_DST"; then
    ok "CMAES_INTEGRATION.md already installed at repo root (identical)"
  else
    cp -f "$DOC_SRC" "$DOC_DST"
    ok "installed CMAES_INTEGRATION.md → repo root"
  fi
else
  skip "CMAES_INTEGRATION.md not next to script (docs stay in swarm artifacts)"
fi

# ─── 5. Compile smoke ───────────────────────────────────────────────────────
log ""
log "[5/5] Compile smoke (cmaes_search.cpp)"

if [[ "${APPLY_INTEGRATION_SKIP_SMOKE:-0}" == "1" ]]; then
  skip "APPLY_INTEGRATION_SKIP_SMOKE=1 — not compiling"
elif [[ ! -f "$REPO_ROOT/LIB/cmaes_search.cpp" ]]; then
  skip "LIB/cmaes_search.cpp absent — wiring only; land chunk1 adapter then re-run smoke"
else
  CXX="${CXX:-}"
  if [[ -z "$CXX" ]]; then
    if command -v g++ >/dev/null 2>&1; then
      CXX="g++"
    elif command -v c++ >/dev/null 2>&1; then
      CXX="c++"
    else
      CXX=""
    fi
  fi
  if [[ -z "$CXX" ]]; then
    skip "no C++ compiler found — wiring applied; compile later via cmake"
  else
    SMOKE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/cmaes_smoke.XXXXXX")"
    set +e
    "$CXX" -std=c++23 -Wall -Wextra -c \
      -I"$REPO_ROOT/LIB" \
      -I"$REPO_ROOT/LIB/tENCoM" \
      -I"$REPO_ROOT/LIB/ShannonThermoStack" \
      -I"$REPO_ROOT/LIB/LigandRingFlex" \
      -I"$REPO_ROOT/LIB/NATURaL" \
      -I"$REPO_ROOT/LIB/ChiralCenter" \
      -I"$REPO_ROOT/LIB/CavityDetect" \
      -I"$REPO_ROOT/LIB/PTMAttachment" \
      "$REPO_ROOT/LIB/cmaes_search.cpp" \
      -o "$SMOKE_DIR/cmaes_search.o" \
      >"$SMOKE_DIR/smoke.log" 2>&1
    rc=$?
    set -e
    if [[ $rc -eq 0 ]]; then
      ok "smoke compile PASS ($CXX -std=c++23 -c LIB/cmaes_search.cpp with -I LIB…)"
    else
      warn "smoke compile exited $rc — log: $SMOKE_DIR/smoke.log (wiring still applied; full link = VALIDATION B)"
      mkdir -p "$ANALYSIS_DIR"
      cp -f "$SMOKE_DIR/smoke.log" "$ANALYSIS_DIR/cmaes_smoke_compile.log" 2>/dev/null || true
    fi
    rm -rf "$SMOKE_DIR"
  fi
fi

# ─── Summary ────────────────────────────────────────────────────────────────
log ""
log "=== Summary ==="
log "  PASS=$PASS_COUNT  WARN=$WARN_COUNT  SKIP=$SKIP_COUNT  FAIL=$FAIL_COUNT"
log "  ic2cf.cpp / gaboom.cpp: never modified by this script (additive surface only)"
log "  Re-run is a no-op when markers + sources are already integrated."

if [[ $FAIL_COUNT -gt 0 ]]; then
  log "FAIL"
  exit 1
fi

log "PASS"
exit 0
