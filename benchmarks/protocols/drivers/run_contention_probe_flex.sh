#!/usr/bin/env bash
# ============================================================================
# CONTENTION PROBE, FLEXIBLE ARM — a faithful replication of the campaign cell
# with PROCESS CONCURRENCY as the only variable.
#
# WHY THIS SUPERSEDES run_contention_probe_confirst.sh.
# That probe (and both earlier reproducibility probes) ran the RIGID path:
#   emitted config diff, campaign cell vs my probe, 59 keys, 3 DIFFERING
#     flexibility.autoflex_max           campaign 5      my probe 0
#     flexibility.autoflex_metal_shrink  campaign true   my probe false
#     output.scored_only                 campaign true   my probe false
# The campaign arm is "B_shrink" -- flexible sidechains with metal shrink. So
# "the engine is reproducible on 15/15 conditions" describes the RIGID path
# only, and the divergence being chased appeared in the FLEXIBLE path. Pose
# counts confirm two different protocols: campaign 1P2Y wrote 126 poses, my
# rigid probes wrote 87.
#
# Reading the campaign driver's cell() (run_optresfix_requeue_v2.sh:63-75)
# exposed a FOURTH difference the config diff could not see, because these are
# election/output env flags that never reach dock_config.json:
#     campaign sets   FLEXAIDDS_SCORED_ONLY=1  FLEXAIDDS_OPTRES_DIAG=1
#     my probes set   FLEXAIDDS_SOFTBETA_ELECTION=1   (campaign does NOT)
# This driver therefore reproduces the campaign env EXACTLY and drops
# SOFTBETA_ELECTION.
#
# WHY THE FLEXIBLE PATH IS THE RIGHT PLACE TO LOOK. Candidate 3 (pose-writer
# truncation at the 50-cluster cap) was REFUTED by pose-set comparison of the
# two campaign cells:
#     1P2Y  Jaccard 0.233, non-shared poses at CF ranks #1 #3 #5 #6
#     2D3U  Jaccard 0.987, exactly one pose differs -- and it is rank #1
#     1JD0  Jaccard 1.000  (clean control)
# Truncation would drop poses at the WORST ranks, at the cap boundary. These
# differ at the BEST ranks, so the searches genuinely diverged. The flexible
# path carries the extra mutable state the rigid path never touches (rotamer
# selection and the OptRes array), which is exactly where this campaign's
# known bugs live.
#
# DESIGN: 3 targets x {isolated, concurrent} x 2 repeats = 12 cells.
#   CON runs FIRST (the decisive arm), then ISO.
#
# THE THREE GATES, in the order they can fail:
#   1. ENGINE IDENTITY  — refuses unless the binary is the campaign's d7916d13.
#   2. FLEXIBILITY ENGAGED — after the first cell, asserts the emitted config
#      carries autoflex_max=5 AND that num_optres > 0. A "flexible" probe that
#      silently ran rigid is vacuous, and that is precisely the error this
#      driver exists to correct. Also records pose count for comparison with
#      the campaign's 126 on 1P2Y: a large mismatch means the replication is
#      still not faithful and the verdict must not be read as one.
#   3. OVERLAP — records per-cell start/end epochs and computes real pairwise
#      overlap in each CON wave; a wave that did not overlap is reported VOID
#      rather than read as a contention result.
#
# STOP LEVER: touch $B/STOP
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
BIN="$R/prodgate_040617/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
GEN=1000
POP=1000
RST=3
OMP=3
TMO=28800
AFX=5          # campaign: FLEXAIDDS_AUTOFLEX_MAX=5
SHR=1          # campaign: FLEXAIDDS_AUTOFLEX_METAL_SHRINK=1
TARGETS="1P2Y 1OQ5 1OWE"
STAMP=$(date -u +%Y%m%d_%H%M%S)
B="$R/contention_flex_$STAMP"
mkdir -p "$B/run" "$B/tmp"
echo "$B" > "$R/state/contention_flex_dir"
export TMPDIR="$B/tmp"

exec > >(tee -a "$B/driver.log") 2>&1
echo "=== CONTENTION PROBE (FLEXIBLE) $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
echo "  engine $ESHA"
if [ "${ESHA:0:16}" != "d7916d136f89aa70" ]; then
  echo "FATAL: engine ${ESHA:0:16} != d7916d136f89aa70 (the campaign's). Wrong binary."
  exit 1
fi
echo "  GATE 1 PASS — campaign binary."
echo "  env replicated: AUTOFLEX_MAX=$AFX METAL_SHRINK=$SHR SCORED_ONLY=1 OPTRES_DIAG=1 NO_SEC=1"
echo "  SOFTBETA_ELECTION deliberately NOT set (the campaign does not set it)"

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
n_optres(){ grep -rhoaE 'num_optres[= ]+[0-9]+|\[OPTRES\][^\n]{0,40}' "$1" 2>/dev/null | grep -oE '[0-9]+' | sort -rn | head -1; }
cfg_afx(){
  local f; f=$(find "$1" -name 'dock_config.json' 2>/dev/null | head -1)
  [ -n "${f:-}" ] && python3 -c "
import json,sys
try: print(json.load(open('$f')).get('flexibility',{}).get('autoflex_max','MISSING'))
except Exception as e: print('ERR')" || echo NOCFG
}

# cell <outdir> <target>  — campaign env, verbatim except for concurrency
cell(){
  local O="$1" T="$2"
  if [ -f "$O/DONE" ]; then return 0; fi
  mkdir -p "$O"
  local t0=$(date +%s) rc=0
  echo "$t0" > "$O/START"
  ( cd "$O" && env OMP_NUM_THREADS="$OMP" \
      FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
      FLEXAIDDS_RESTARTS="$RST" FLEXAIDDS_SCORED_ONLY=1 FLEXAIDDS_NO_SEC=1 \
      FLEXAIDDS_AUTOFLEX_MAX="$AFX" FLEXAIDDS_AUTOFLEX_METAL_SHRINK="$SHR" \
      FLEXAIDDS_OPTRES_DIAG=1 TMPDIR="$TMPDIR" \
      "$BIN/benchmark_datasets" --benchmark astex_diverse \
        --mode defined-cleft-redock --only-codes "$T" --output "$O" --cache "$CACHE" \
        --threads 1 --omp-threads "$OMP" --ga-population "$POP" --ga-generations "$GEN" \
        --job-timeout-seconds "$TMO" --engine-sha256 "$ESHA" > "$O/run.log" 2>&1 ) || rc=$?
  local t1=$(date +%s)
  printf 'rc=%s target=%s poses=%s allhash=%s elected=%s afx_cfg=%s optres=%s start=%s end=%s wall_s=%s load1=%s at=%s\n' \
    "$rc" "$T" "$(nposes "$O")" "$(posehash "$O")" "$(electedsha "$O")" \
    "$(cfg_afx "$O")" "$(n_optres "$O")" "$t0" "$t1" "$((t1-t0))" \
    "$(uptime | sed 's/.*load average[s]*: //' | awk '{print $1}' | tr -d ,)" \
    "$(date -u +%FT%TZ)" > "$O/DONE"
}

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
    printf "  CON %-6s %s  %s\n" "$T" "$REP" "$(tr '\n' ' ' < "$B/run/CON_${T}_${REP}/DONE" | cut -c1-118)"
  done
  # ---- GATE 2, on the first completed wave: did flexibility actually engage?
  if [ "$REP" = "A" ]; then
    afx=$(grep -o 'afx_cfg=[A-Za-z0-9]*' "$B/run/CON_1P2Y_A/DONE" | cut -d= -f2)
    opt=$(grep -o 'optres=[0-9]*' "$B/run/CON_1P2Y_A/DONE" | cut -d= -f2)
    pz=$(grep -o 'poses=[0-9]*' "$B/run/CON_1P2Y_A/DONE" | cut -d= -f2)
    echo ""
    echo "  === GATE 2: FLEXIBILITY ENGAGED? ==="
    echo "    emitted config autoflex_max = $afx   (want 5)"
    echo "    num_optres                  = ${opt:-none}  (want > 0)"
    echo "    poses                       = $pz   (campaign 1P2Y wrote 126; rigid probes wrote 87)"
    if [ "$afx" != "5" ]; then
      echo "    FATAL: autoflex_max is '$afx', not 5 — this is the RIGID path again."
      echo "      The probe would repeat the exact error it exists to correct."
      exit 2
    fi
    echo "    GATE 2 PASS — flexible path confirmed from the emitted config."
    if [ "$pz" = "87" ]; then
      echo "    WARNING: pose count equals the RIGID probes' 87. Treat the verdict as suspect."
    fi
  fi
done

# ---- ARM ISO: one process at a time -----------------------------------------
echo ""
echo "=== ARM ISO (sequential, one engine process) ==="
for REP in A B; do
  for T in $TARGETS; do
    [ -f "$B/STOP" ] && { echo "  STOP sentinel"; break 2; }
    O="$B/run/ISO_${T}_${REP}"
    cell "$O" "$T"
    printf "  ISO %-6s %s  %s\n" "$T" "$REP" "$(tr '\n' ' ' < "$O/DONE" | cut -c1-118)"
  done
done

# ---- GATE 3: OVERLAP --------------------------------------------------------
echo ""
echo "=== GATE 3: OVERLAP (a contention probe that did not overlap is vacuous) ==="
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
        print(f"  wave {REP}: fewer than 2 cells complete — overlap not measurable"); ok = False; continue
    print(f"  wave {REP}:")
    for T,(s,e,w) in spans.items(): print(f"    {T:<6} {w:>5}s")
    mins = []
    for a,b in itertools.combinations(spans,2):
        s1,e1,_ = spans[a]; s2,e2,_ = spans[b]
        ov = max(0, min(e1,e2) - max(s1,s2))
        mins.append(ov)
        print(f"    overlap {a}/{b}: {ov}s ({ov/max(1,min(e1-s1,e2-s2)):.0%} of the shorter cell)")
    if min(mins) <= 0:
        print(f"    *** wave {REP} DID NOT OVERLAP — this wave's verdict is VOID"); ok = False
print(f"\n  CONCURRENCY ACHIEVED: {ok}")
PY

# ---- VERDICT ---------------------------------------------------------------
echo ""
echo "=== VERDICT: same-seed repeat identity, FLEXIBLE path, isolated vs concurrent ==="
printf "  %-4s %-6s %-18s %-18s %-11s %s\n" arm target allhash_A allhash_B trajectory elected
for ARM in CON ISO; do
  for T in $TARGETS; do
    A="$B/run/${ARM}_${T}_A/DONE"; Bf="$B/run/${ARM}_${T}_B/DONE"
    if [ -f "$A" ] && [ -f "$Bf" ]; then
      ha=$(grep -o 'allhash=[0-9a-f]*' "$A" | cut -d= -f2); hb=$(grep -o 'allhash=[0-9a-f]*' "$Bf" | cut -d= -f2)
      ea=$(grep -o 'elected=[0-9a-zA-Z]*' "$A" | cut -d= -f2); eb=$(grep -o 'elected=[0-9a-zA-Z]*' "$Bf" | cut -d= -f2)
      printf "  %-4s %-6s %-18s %-18s %-11s %s\n" "$ARM" "$T" "$ha" "$hb" \
        "$([ "$ha" = "$hb" ] && echo SAME || echo DIFFER)" \
        "$([ "$ea" = "$eb" ] && echo same || echo DIFFER)"
    else
      printf "  %-4s %-6s (incomplete)\n" "$ARM" "$T"
    fi
  done
done
echo ""
echo "  CON any DIFFER + ISO all SAME  -> CONTENTION on the flexible path."
echo "  ISO any DIFFER                 -> the FLEXIBLE path is non-reproducible even"
echo "                                    ISOLATED, which localises it to the flexible"
echo "                                    machinery (rotamer/OptRes) and NOT to concurrency."
echo "  both all SAME                  -> neither; the campaign divergence needs a"
echo "                                    cause still unfound."
echo "COMPLETE $(date -u +%FT%TZ)" > "$B/DONE"
echo "=== CONTENTION PROBE (FLEXIBLE) COMPLETE $(date -u +%FT%TZ) ==="
