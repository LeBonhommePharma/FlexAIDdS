#!/usr/bin/env bash
# ============================================================================
# CONFOUND PROBE — separates "newer build fixed it" from "gen=300 too short"
#
# WHY. determinism_probe_20260906_045811 found the campaign's fixed-seed
# divergence DOES NOT REPRODUCE: on both 1P2Y and 2D3U, same-seed repeats are
# byte-identical at threads=1 AND threads=3, flag on AND off. But that probe
# differs from the campaign on TWO axes at once, so the null is uninterpretable:
#
#   campaign (diverged)   engine d7916d136f89aa70   POP=1000  GEN=1000
#   determinism probe     engine 74b407a450a361a0   POP=1000  GEN=300
#
# This probe holds the ENGINE at the campaign's own binary (which survived at
# prodgate_040617/bin/FlexAIDdS) and varies ONLY the generation count, so the
# two explanations separate:
#
#   diverges at GEN=1000 but not 300  -> budget-dependent; my probe was short
#   diverges at BOTH                  -> build 74b407a4 fixed it; find the commit
#   identical at BOTH                 -> cause is OUTSIDE the engine (harness),
#                                        and the P1-vs-P2 divergence needs a
#                                        different explanation entirely
#
# Target 2D3U: its exact-pose audit was 0/24000 inconsistent, so its CF was
# already reproducible, yet its ELECTED pose moved between P1 and P2. That
# makes it the informative one -- the Voronoi jitter cannot be its cause.
# 1P2Y added at GEN=1000 only, as the high-audit-inconsistency counterpart.
#
# threads=3 and flag-off throughout: those are the campaign's own settings.
# Not varied here -- the determinism probe already showed both are inert.
#
# STOP LEVER: touch $B/STOP.
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
BIN="$R/prodgate_040617/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
STAMP=$(date -u +%Y%m%d_%H%M%S)
B="$R/confound_probe_$STAMP"
mkdir -p "$B/run" "$B/tmp"
echo "$B" > "$R/state/confound_probe_dir"
export TMPDIR="$B/tmp"

exec > >(tee -a "$B/driver.log") 2>&1
echo "=== CONFOUND PROBE  $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
RSHA=$(shasum -a256 "$BIN/benchmark_datasets" 2>/dev/null | awk '{print $1}')
echo "  engine $ESHA"
echo "  runner ${RSHA:-MISSING}"

# GATE: this must be the CAMPAIGN's engine, or the probe answers nothing.
if [ "${ESHA:0:16}" != "d7916d136f89aa70" ]; then
  echo "FATAL: engine is ${ESHA:0:16}, expected d7916d136f89aa70 (the campaign's)."
  echo "  Varying generations on the WRONG binary cannot separate the confound."
  exit 1
fi
if [ -z "${RSHA:-}" ]; then
  echo "FATAL: no benchmark_datasets beside that engine at $BIN"
  exit 1
fi
echo "  GATE PASS — running the campaign's own binary."

posehash(){
  for f in $(find "$1" -name '*.pdb' 2>/dev/null | grep -v _INI | sort); do
    grep '^\(ATOM\|HETATM\)' "$f" | cut -c31-54
  done | shasum -a256 | cut -c1-16
}
nposes(){ find "$1" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0; }
electedsha(){
  local c; c=$(find "$1" -name 'astex_diverse_results.csv' 2>/dev/null | head -1)
  [ -n "${c:-}" ] && awk -F, 'NR==2{print substr($20,1,16)}' "$c" || echo NA
}

# cell <outdir> <target> <generations> <rep>
cell(){
  local O="$1" T="$2" GEN="$3" REP="$4"
  if [ -f "$O/DONE" ]; then return 0; fi
  mkdir -p "$O"
  local t0=$(date +%s) rc=0
  ( cd "$O" && env OMP_NUM_THREADS=3 \
      FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
      FLEXAIDDS_RESTARTS=3 FLEXAIDDS_NO_SEC=1 FLEXAIDDS_SOFTBETA_ELECTION=1 \
      "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
        --only-codes "$T" --output "$O" --cache "$CACHE" \
        --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
        --threads 1 --omp-threads 3 --ga-population 1000 --ga-generations "$GEN" \
        --job-timeout-seconds 28800 > "$O/stdout.log" 2>&1 ) || rc=$?
  local t1=$(date +%s)
  local aud; aud=$(grep -haoE 'inconsistent=[0-9]+/[0-9]+ max_delta=[0-9.]+' "$O/stdout.log" 2>/dev/null | tr '\n' ';')
  printf 'rc=%s target=%s gen=%s rep=%s poses=%s allhash=%s elected=%s wall_s=%s audit=%s at=%s\n' \
    "$rc" "$T" "$GEN" "$REP" "$(nposes "$O")" "$(posehash "$O")" "$(electedsha "$O")" \
    "$((t1-t0))" "${aud:-none}" "$(date -u +%FT%TZ)" > "$O/DONE"
}

i=0
for spec in "2D3U 300" "2D3U 1000" "1P2Y 1000"; do
  set -- $spec; T=$1; GEN=$2
  for REP in A B; do
    if [ -f "$B/STOP" ]; then echo "  STOP sentinel — halting"; break 3; fi
    i=$((i+1))
    O="$B/run/${T}_g${GEN}_${REP}"
    cell "$O" "$T" "$GEN" "$REP"
    printf "  [%s/6] %-6s gen=%-5s rep=%s  %s\n" "$i" "$T" "$GEN" "$REP" \
      "$(tr '\n' ' ' < "$O/DONE" | cut -c1-110)"
  done
done

echo ""
echo "=== VERDICT: same-seed repeat identity on the CAMPAIGN binary ==="
printf "  %-6s %-6s %-18s %-18s %s\n" target gen all-poses elected repeats
for spec in "2D3U 300" "2D3U 1000" "1P2Y 1000"; do
  set -- $spec; T=$1; GEN=$2
  A="$B/run/${T}_g${GEN}_A/DONE"; Bf="$B/run/${T}_g${GEN}_B/DONE"
  [ -f "$A" ] && [ -f "$Bf" ] || { printf "  %-6s %-6s (incomplete)\n" "$T" "$GEN"; continue; }
  ha=$(grep -o 'allhash=[0-9a-f]*' "$A" | cut -d= -f2); hb=$(grep -o 'allhash=[0-9a-f]*' "$Bf" | cut -d= -f2)
  ea=$(grep -o 'elected=[0-9a-zA-Z]*' "$A" | cut -d= -f2); eb=$(grep -o 'elected=[0-9a-zA-Z]*' "$Bf" | cut -d= -f2)
  printf "  %-6s %-6s %-18s %-18s all=%s elected=%s\n" "$T" "$GEN" "$ha" "$ea" \
    "$([ "$ha" = "$hb" ] && echo same || echo DIFFER)" \
    "$([ "$ea" = "$eb" ] && echo same || echo DIFFER)"
done
echo ""
echo "  2D3U DIFFER at gen=1000 and same at gen=300 -> budget-dependent."
echo "  2D3U DIFFER at BOTH                         -> build 74b407a4 fixed it."
echo "  2D3U SAME at both                           -> cause is outside the engine."
echo "COMPLETE $(date -u +%FT%TZ)" > "$B/DONE"
echo "=== CONFOUND PROBE COMPLETE  $(date -u +%FT%TZ) ==="
