#!/usr/bin/env bash
# format_check.sh — non-destructive style check for FlexAIDdS.
#
# Does NOT reformat the tree. Exits non-zero if checked files would change
# under the repo .clang-format.
#
# Usage:
#   scripts/format_check.sh                 # default allowlist (modern modules)
#   scripts/format_check.sh --all-new       # untracked + staged C++ only
#   scripts/format_check.sh path/to/file.cpp
#   scripts/format_check.sh --list          # print default allowlist
#   scripts/format_check.sh --tidy          # optional clang-tidy on allowlist
#
# Mass reformat of LIB/ is intentionally out of scope for this PR series.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Default allowlist: newer modules written in 4-space modern style.
# Expand carefully; do not dump gaboom.cpp / Vcontacts.cpp here until
# those files are deliberately modernized.
DEFAULT_ALLOWLIST=(
  "LIB/statmech.cpp"
  "LIB/statmech.h"
  "LIB/encom.cpp"
  "LIB/encom.h"
  "LIB/UnifiedHardwareDispatch.h"
  "LIB/VoronoiCFBatch.h"
  "python/bindings/core_bindings.cpp"
)

CLANG_FORMAT_BIN="${CLANG_FORMAT:-clang-format}"
CLANG_TIDY_BIN="${CLANG_TIDY:-clang-tidy}"
DO_TIDY=0
MODE="allowlist"
EXTRA_PATHS=()

usage() {
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --list)
      printf '%s\n' "${DEFAULT_ALLOWLIST[@]}"
      exit 0
      ;;
    --all-new) MODE="new"; shift ;;
    --tidy) DO_TIDY=1; shift ;;
    --) shift; EXTRA_PATHS+=("$@"); break ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *) EXTRA_PATHS+=("$1"); shift ;;
  esac
done

if ! command -v "$CLANG_FORMAT_BIN" >/dev/null 2>&1; then
  echo "ERROR: clang-format not found (set CLANG_FORMAT=... or install llvm)." >&2
  exit 127
fi

collect_new_files() {
  # Untracked + staged/modified C/C++ under LIB/, src/, tests/, python/bindings/
  {
    git ls-files --others --exclude-standard -- 'LIB/**' 'src/**' 'tests/**' 'python/bindings/**' 2>/dev/null || true
    git diff --name-only --diff-filter=ACMR HEAD -- 'LIB/**' 'src/**' 'tests/**' 'python/bindings/**' 2>/dev/null || true
    git diff --cached --name-only --diff-filter=ACMR -- 'LIB/**' 'src/**' 'tests/**' 'python/bindings/**' 2>/dev/null || true
  } | sort -u | while read -r f; do
    case "$f" in
      *.c|*.cc|*.cpp|*.cxx|*.h|*.hh|*.hpp|*.cu|*.cuh|*.mm|*.m) echo "$f" ;;
    esac
  done
}

FILES=()
if [[ ${#EXTRA_PATHS[@]} -gt 0 ]]; then
  FILES=("${EXTRA_PATHS[@]}")
elif [[ "$MODE" == "new" ]]; then
  while IFS= read -r f; do
    [[ -n "$f" ]] && FILES+=("$f")
  done < <(collect_new_files)
else
  FILES=("${DEFAULT_ALLOWLIST[@]}")
fi

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No files to check."
  exit 0
fi

echo "clang-format: $($CLANG_FORMAT_BIN --version | head -1)"
echo "Checking ${#FILES[@]} file(s) against .clang-format (dry-run)..."

failed=0
checked=0
missing=0
for f in "${FILES[@]}"; do
  if [[ ! -f "$f" ]]; then
    echo "  SKIP (missing): $f"
    missing=$((missing + 1))
    continue
  fi
  checked=$((checked + 1))
  # --dry-run --Werror available in clang-format >= 10
  if "$CLANG_FORMAT_BIN" --dry-run --Werror "$f" >/dev/null 2>&1; then
    echo "  OK  $f"
  else
    echo "  DIFF $f  (would reformat; run: $CLANG_FORMAT_BIN -i $f)"
    failed=$((failed + 1))
  fi
done

echo
echo "Summary: checked=$checked missing=$missing would_reformat=$failed"

if [[ "$DO_TIDY" -eq 1 ]]; then
  if ! command -v "$CLANG_TIDY_BIN" >/dev/null 2>&1; then
    echo "WARN: clang-tidy not found; skipping --tidy" >&2
  elif [[ ! -f build/compile_commands.json && ! -f compile_commands.json ]]; then
    echo "WARN: no compile_commands.json (configure with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON); skipping tidy" >&2
  else
    COMPDB="$ROOT/build"
    [[ -f compile_commands.json ]] && COMPDB="$ROOT"
    echo "Running clang-tidy on allowlist (advisory)..."
    for f in "${FILES[@]}"; do
      [[ -f "$f" ]] || continue
      case "$f" in
        *.h|*.hpp) continue ;;  # tidy needs a TU; headers via main files
      esac
      "$CLANG_TIDY_BIN" -p "$COMPDB" "$f" || true
    done
  fi
fi

if [[ "$failed" -gt 0 ]]; then
  echo
  echo "format_check: FAIL ($failed file(s) differ from .clang-format)"
  echo "This script never rewrites files. Format intentionally with clang-format -i."
  exit 1
fi

echo "format_check: PASS"
exit 0
