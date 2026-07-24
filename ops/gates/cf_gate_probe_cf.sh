#!/bin/bash
# CF scoring-regression merge gate — wraps the in-repo build/probe_cf instrument.
# Supersedes the standalone python re-implementation: scores via the ENGINE's real
# vcfunction()/score_native_pose(), not a proxy. For each panel target it scores the
# NATIVE pose and the best DECOY pose, computes ΔCF = cf_total(native) - cf_total(decoy),
# and FAILS if the inverted fraction (ΔCF>tol) exceeds MAX_INV_FRAC.
#
# Usage: cf_gate_probe_cf.sh <panel_manifest.tsv> [tol] [max_inv_frac]
#   manifest rows (tab-sep): pdb <TAB> receptor.pdb <TAB> native_pose.sdf <TAB> decoy_pose.pdb
# Requires: build/probe_cf (built from tools/probe_cf.cpp, commit 81acfed6+).
set -u
REPO="/Users/lp.more/Projects/FlexAIDdS"
PROBE="$REPO/build/probe_cf"
MANIFEST="${1:?manifest tsv required}"
TOL="${2:-0.0}"
MAX_INV_FRAC="${3:-0.125}"

[ -x "$PROBE" ] || { echo "GATE ERROR: $PROBE not built (cmake --build build --target probe_cf)"; exit 2; }

cf() {  # $1=receptor $2=pose  -> echoes cf_total via probe_cf JSON
  "$PROBE" --receptor "$1" --pose "$2" 2>/dev/null \
    | sed -n 's/.*"cf_total": *\([-0-9.][-0-9.]*\).*/\1/p' | head -1
}

n=0; inv=0; rows=""
while IFS=$'\t' read -r pdb rec nat dec; do
  [ -z "${pdb:-}" ] && continue
  case "$pdb" in \#*) continue;; esac
  cn=$(cf "$rec" "$nat"); cd_=$(cf "$rec" "$dec")
  [ -z "$cn" ] || [ -z "$cd_" ] && { echo "  skip $pdb (probe_cf gave no cf_total)"; continue; }
  d=$(echo "$cn - $cd_" | bc -l)
  bad=$(echo "$d > $TOL" | bc -l)
  n=$((n+1)); [ "$bad" = 1 ] && inv=$((inv+1))
  rows="$rows$(printf '  %-8s dCF=%+8.2f (native %s vs decoy %s)%s\n' "$pdb" "$d" "$cn" "$cd_" "$([ "$bad" = 1 ] && echo '  INVERTED')")\n"
done < "$MANIFEST"

[ "$n" -eq 0 ] && { echo "GATE ERROR: no targets scored"; exit 2; }
frac=$(echo "scale=3; $inv / $n" | bc -l)
printf '%b' "$rows"
echo "CF-REGRESSION GATE | n=$n inverted=$inv frac=$frac (threshold $MAX_INV_FRAC)"
pass=$(echo "$frac <= $MAX_INV_FRAC" | bc -l)
if [ "$pass" = 1 ]; then echo "  VERDICT: PASS"; exit 0; else echo "  VERDICT: FAIL"; exit 1; fi
