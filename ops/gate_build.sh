#!/bin/bash
# METHODOLOGY.md §0 — build an engine variant. Usage: gate_build.sh <target> [out_copy]
set -e
export PATH="/opt/homebrew/bin:$PATH"
REPO="/Users/lp.more/Projects/FlexAIDdS"; cd "$REPO/build"
TGT="${1:-FlexAIDdS}"
/opt/homebrew/bin/cmake --build . --target "$TGT" -j4
echo "built $TGT md5=$(md5 -q $TGT 2>/dev/null || md5sum $TGT | cut -d' ' -f1)"
[ -n "${2:-}" ] && cp "$TGT" "$2" && echo "copied -> $2"
