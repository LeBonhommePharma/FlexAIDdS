#!/usr/bin/env bash
# ============================================================================
# CONTENTION PROBE — is the campaign's fixed-seed divergence caused by
# PROCESS-LEVEL CONTENTION rather than by anything inside one engine run?
#
# WHAT IS ALREADY MEASURED (2026-09-06):
#   determinism_probe_20260906_045811  engine 74b407a4, gen=300, SEQUENTIAL
#     12/12 conditions byte-identical; omp-threads 1 vs 3 INERT on all 6 pairs
#   confound_probe_20260906_060321     engine d7916d13, gen=300 and 1000, SEQUENTIAL
#     3/3 conditions byte-identical; 1P2Y gen=1000 -> b8eea006d19ac629 twice
#   campaign requeue_20260904_050430   engine d7916d13, gen=1000, WINDOW=3
#     1P2Y and 2D3U elected DIFFERENT poses between P1 and P2 at the same seed
#
# WHAT RULES OUT EVERYTHING EXCEPT CONTENTION:
#   The two divergent campaign cells had a BYTE-IDENTICAL emitted config
#   (dock_config.json sha 3bb0ead6f6b40734, 0 of 59 keys differing), the same
#   seed, the same binary, and even the SAME co-scheduled neighbours
#   (1OQ5, 1OWE). Their walls differed by 12 s on identical work.
#
# THE AXIS NEITHER PROBE TESTED:
#   Both probes ran ONE engine process at a time on an idle box (0 backgrounded
#   dispatches in either driver). The campaign ran WINDOW=3, i.e. 3 concurrent
#   engine processes x 3 OMP threads = 9 threads on 11 cores. Varying
#   --omp-threads inside one process on a quiet machine does NOT reproduce the
#   thread scheduling, work-stealing order or timing of 3 processes competing.
#   So "the engine is reproducible" is established only IN ISOLATION.
#
# DESIGN: 3 targets x {isolated, concurrent} x 2 repeats.
#   ISO: each target run alone, sequentially            6 cells
#   CON: all 3 targets launched simultaneously, twice   6 cells (2 waves of 3)
#   The concurrent arm tests 3 targets at once, so it gives THREE independent
#   tests of the hypothesis per wave rather than one.
#
# PREDICTION (falsifiable both ways):
#   ISO repeats identical AND CON repeats DIFFER  -> contention is the mechanism
#   both identical                                -> contention is NOT it either,
#                                                    and the campaign divergence
#                                                    has a cause still unfound
#   ISO also differs                              -> my confound probe's null was
#                                                    luck, not reproducibility
#
# THE LOAD-BEARING GATE: a contention experiment whose cells did not actually
# overlap in time is VACUOUS. This driver records per-cell start/end epochs and
# computes the real pairwise overlap in the CON waves, then refuses to report a
# verdict if the overlap is zero. Assuming '&' worked is exactly the class of
# error that made an earlier driver silently serial.
#
# STOP LEVER: touch $B/STOP.
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
BIN="$R/prodgate_040617/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
GEN=1000
TARGETS="1P2Y 1OQ5 1OWE"
STAMP=$(date -u +%Y%m%d_%H%M%S)
PREV=$(ls -d "$R"/contention_probe_* 2>/dev/null | sort | tail -1)
if [ -n "${PREV:-}" ] && [ -d "$PREV/run" ]; then B="$PREV"; echo "  REUSING batch $(basename "$B")"; else B="$R/contention_probe_$STAMP"; fi
mkdir -p "$B/run" "$B/tmp"
echo "$B" > "$R/state/contention_probe_dir"
export TMPDIR="$B/tmp"

exec > >(tee -a "$B/driver.log") 2>&1
echo "=== CONTENTION PROBE  $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
RSHA=$(shasum -a256 "$BIN/benchmark_datasets" | awk '{print $1}')
echo "  engine $ESHA"
echo "  runner $RSHA"
if [ "${ESHA:0:16}" != "d7916d136f89aa70" ]; then
  echo "FATAL: engine ${ESHA:0:16} != d7916d136f89aa70 (the campaign's). Wrong binary."
  exit 1
fi
echo "  GATE PASS — campaign binary.  targets: $TARGETS  gen=$GEN"

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

# cell <outdir> <target>   — records START and END epochs for the overlap gate
cell(){
  local O="$1" T="$2"
  if [ -f "$O/DONE" ]; then return 0; fi
  mkdir -p "$O"
  local t0=$(date +%s) rc=0
  echo "$t0" > "$O/START"
  ( cd "$O" && env OMP_NUM_THREADS=3 \
      FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
      FLEXAIDDS_RESTARTS=3 FLEXAIDDS_NO_SEC=1 FLEXAIDDS_SOFTBETA_ELECTION=1 \
      "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
        --only-codes "$T" --output "$O" --cache "$CACHE" \
        --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
        --threads 1 --omp-threads 3 --ga-population 1000 --ga-generations "$GEN" \
        --job-timeout-seconds 28800 > "$O/stdout.log" 2>&1 ) || rc=$?
  local t1=$(date +%s)
  printf 'rc=%s target=%s poses=%s allhash=%s elected=%s start=%s end=%s wall_s=%s load1=%s at=%s\n' \
    "$rc" "$T" "$(nposes "$O")" "$(posehash "$O")" "$(electedsha "$O")" \
    "$t0" "$t1" "$((t1-t0))" \
    "$(uptime | sed 's/.*load average[s]*: //' | awk '{print $1}' | tr -d ,)" \
    "$(date -u +%FT%TZ)" > "$O/DONE"
}

# ARM ORDER REVERSED vs run_contention_probe.sh: CON runs FIRST.
# The ISO arm is largely redundant -- confound_probe_20260906_060321 already
# measured 1P2Y isolated at gen=1000 TWICE, byte-identical (b8eea006d19ac629).
# CON is the decisive arm and costs ~25 min against ISO's ~38, so CON first
# delivers the verdict sooner. cell() skips any cell with a DONE, so cells
# already completed by the previous launch are reused, not redone.
# This is a NEW FILE, not an edit of a script that may still be draining.

# ---- ARM CON: all three simultaneously, twice -------------------------------
echo ""
echo "=== ARM CON (3 engine processes launched simultaneously) ==="
for REP in A B; do
  [ -f "$B/STOP" ] && { echo "  STOP sentinel"; break; }
  echo "  --- wave $REP: launching $TARGETS together at $(date -u +%H:%M:%SZ) ---"
  for T in $TARGETS; do
    cell "$B/run/CON_${T}_${REP}" "$T" &
  done
  wait
  for T in $TARGETS; do
    printf "  CON %-6s %s  %s\n" "$T" "$REP" "$(tr '\n' ' ' < "$B/run/CON_${T}_${REP}/DONE" | cut -c1-104)"
  done
done

# ---- ARM ISO: one process at a time -----------------------------------------
echo ""
echo "=== ARM ISO (sequential, one engine process) ==="
for REP in A B; do
  for T in $TARGETS; do
    [ -f "$B/STOP" ] && { echo "  STOP sentinel"; break 3; }
    O="$B/run/ISO_${T}_${REP}"
    cell "$O" "$T"
    printf "  ISO %-6s %s  %s\n" "$T" "$REP" "$(tr '\n' ' ' < "$O/DONE" | cut -c1-104)"
  done
done

# ---- THE OVERLAP GATE: did CON actually run concurrently? -------------------
echo ""
echo "=== OVERLAP GATE (a contention probe that did not overlap is vacuous) ==="
python3 - "$B" "$TARGETS" <<'PY'
import os, sys, itertools
B = sys.argv[1]; targets = sys.argv[2].split()
def rec(p):
    if not os.path.exists(p): return None
    return dict(t.split("=",1) for t in open(p).read().split() if "=" in t)
ok = True
for REP in ("A","B"):
    spans = {}
    for T in targets:
        d = rec(f"{B}/run/CON_{T}_{REP}/DONE")
        if d: spans[T] = (int(d["start"]), int(d["end"]), int(d["wall_s"]))
    if len(spans) < 2:
        print(f"  wave {REP}: fewer than 2 cells complete — no overlap measurable"); ok = False; continue
    print(f"  wave {REP}:")
    for T,(s,e,w) in spans.items():
        print(f"    {T:<6} {w:>5}s  span {s}..{e}")
    mins = []
    for a,b in itertools.combinations(spans,2):
        s1,e1,_ = spans[a]; s2,e2,_ = spans[b]
        ov = max(0, min(e1,e2) - max(s1,s2))
        frac = ov / max(1, min(e1-s1, e2-s2))
        mins.append(ov)
        print(f"    overlap {a}/{b}: {ov}s ({frac:.0%} of the shorter cell)")
    if min(mins) <= 0:
        print(f"    *** wave {REP} DID NOT OVERLAP — verdict for this wave is VOID"); ok = False
print(f"\n  CONCURRENCY ACHIEVED: {ok}")
if not ok:
    print("  -> do NOT read the CON verdict below as a contention result.")
PY

# ---- VERDICT ---------------------------------------------------------------
echo ""
echo "=== VERDICT: same-seed repeat identity, isolated vs concurrent ==="
printf "  %-4s %-6s %-18s %-18s %-9s %s\n" arm target allhash_A allhash_B trajectory elected
for ARM in ISO CON; do
  for T in $TARGETS; do
    A="$B/run/${ARM}_${T}_A/DONE"; Bf="$B/run/${ARM}_${T}_B/DONE"
    if [ -f "$A" ] && [ -f "$Bf" ]; then
      ha=$(grep -o 'allhash=[0-9a-f]*' "$A" | cut -d= -f2); hb=$(grep -o 'allhash=[0-9a-f]*' "$Bf" | cut -d= -f2)
      ea=$(grep -o 'elected=[0-9a-zA-Z]*' "$A" | cut -d= -f2); eb=$(grep -o 'elected=[0-9a-zA-Z]*' "$Bf" | cut -d= -f2)
      printf "  %-4s %-6s %-18s %-18s %-9s %s\n" "$ARM" "$T" "$ha" "$hb" \
        "$([ "$ha" = "$hb" ] && echo SAME || echo DIFFER)" \
        "$([ "$ea" = "$eb" ] && echo same || echo DIFFER)"
    else
      printf "  %-4s %-6s (incomplete)\n" "$ARM" "$T"
    fi
  done
done
echo ""
echo "  ISO all SAME + CON any DIFFER  ->  CONTENTION is the mechanism."
echo "  both all SAME                  ->  contention is not it; cause still unfound."
echo "  ISO any DIFFER                 ->  the earlier sequential nulls were luck."
echo "COMPLETE $(date -u +%FT%TZ)" > "$B/DONE"
echo "=== CONTENTION PROBE COMPLETE  $(date -u +%FT%TZ) ==="
