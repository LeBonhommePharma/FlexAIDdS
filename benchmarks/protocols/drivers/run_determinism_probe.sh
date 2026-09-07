#!/usr/bin/env bash
# ============================================================================
# DETERMINISM PROBE — what makes the engine non-reproducible at fixed seed?
#
# MEASURED FACT THAT MOTIVATES THIS (requeue_20260904_050430, OMP=3):
#   2 of 84 elected poses differ between two same-seed runs (1P2Y, 2D3U).
#   The engine's own exact-pose audit (top.cpp:3425) reports the GA's stored
#   score disagreeing with a master re-score in 435/756 restarts, median
#   max_delta 40.9 CF units.
#
# HYPOTHESIS (code evidence, Vcontacts.cpp:851-855):
#   On a Voronoi hull degeneracy the engine perturbs atom coordinates by a
#   +-0.005 A jitter drawn from lazy_thread_rng (RngSeed.h:181). That RNG is
#   thread_local and its draw sequence advances by however many degeneracies
#   each thread happens to process -> depends on OpenMP scheduling -> CF is
#   not reproducible when omp-threads > 1.
#
# THE FIX THAT ALREADY EXISTS, DEFAULT OFF -- AND THE RIGHT GATE FOR IT:
#   Vcontacts.cpp:841  if (flexaids_rng::voronoi_keyed_jitter_enabled()) {
#   Vcontacts.cpp:842    // Independent of FLEXAIDDS_RNG_STREAM_FIX (D6).
#   Vcontacts.cpp:847-849  keyed_jitter(pose_atom_identity(...), k)
#   Vcontacts.cpp:850-855  else -> legacy vc_dist(lazy_thread_rng(0x0C0A11))
#   RngSeed.h:160-162  voronoi_keyed_jitter_enabled() =
#                        env_bool("FLEXAIDDS_VORONOI_KEYED_JITTER", false)
#
#   CORRECTED (review finding, and verified from the two sources above): the
#   schedule-independent keyed jitter is selected by
#   FLEXAIDDS_VORONOI_KEYED_JITTER, NOT by FLEXAIDDS_RNG_STREAM_FIX. The first
#   version of this script set RNG_STREAM_FIX=1 in its ON arm, which flips only
#   the per-stream generator map inside lazy_thread_rng and leaves the call site
#   on the legacy else-branch -- so all 24 cells would have run the SAME jitter
#   and the probe could not have tested its own hypothesis. Binary check
#   confirms both names are present (engine 1 hit each), so the gate is live.
#
#   THE VARIED VARIABLE IS FLEXAIDDS_VORONOI_KEYED_JITTER, AND ONLY IT.
#   FLEXAIDDS_RNG_STREAM_FIX is deliberately left unset in BOTH arms: it governs
#   a different consumer (the generator map), so varying it here would confound
#   the jitter measurement. In the OFF arm the variable is not set AT ALL rather
#   than set to 0, because tonight's NO_SEC finding showed "=0" is not a safe
#   way to express off (that one is presence-tested; this one is env_bool, but
#   the safe convention is worth keeping).
#
# DESIGN: 3 targets x omp-threads{1,3} x flag{off,on} x 2 repeats = 24 cells.
#   1P2Y  audit inconsistent 2256/2269 (99.4%) -- hits the degeneracy branch
#         constantly, and its elected pose DID diverge. Strongest signal.
#   2D3U  audit inconsistent 0/24000 -- never hits it, CF is reproducible,
#         YET its elected pose still diverged. If this one still diverges at
#         threads=1 with the flag on, there is a SECOND source in the search
#         path and the jitter is not the whole story.
#   1JD0  control, audit 7/820.
#
# PREDICTIONS (any of these failing falsifies the hypothesis):
#   threads=1, keyed off -> repeats IDENTICAL (single thread = fixed draw order)
#   threads=3, keyed off -> repeats DIFFER      (schedule-dependent draw order)
#   threads=3, keyed ON  -> repeats IDENTICAL   (jitter from atom identity)
#   keyed ON vs OFF at threads=1 -> if the POSES differ, the gate changes
#         physics and not only reproducibility, which decides whether it can
#         ever be default-on. Measured, not assumed.
#
# Budget gen=300 not 1000: divergence appears in the first generations, and
# this keeps the probe ~3 h instead of ~10 h. Determinism is an A-vs-B
# identity at IDENTICAL settings, so a shorter budget does not invalidate it.
#
# STOP LEVER (the thing run_water_ablation.sh lacks): touch $B/STOP and the
# loop exits at the next cell boundary without writing a fake receipt.
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
WA=$(ls -d "$R"/water_ablation_* 2>/dev/null | sort | tail -1)   # reuse its built binaries
BIN="$WA/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
STAMP=$(date -u +%Y%m%d_%H%M%S)
B="$R/determinism_probe_$STAMP"
mkdir -p "$B/run" "$B/tmp"
echo "$B" > "$R/state/determinism_probe_dir"
export TMPDIR="$B/tmp"

exec > >(tee -a "$B/driver.log") 2>&1
echo "=== DETERMINISM PROBE  $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
RSHA=$(shasum -a256 "$BIN/benchmark_datasets" | awk '{print $1}')
echo "  engine $ESHA"
echo "  runner $RSHA"
echo "  cache  $CACHE"
echo "  sites  $SITES"

# GATE: the flag must exist in the binary, else the probe cannot vary its own
# variable and would produce a null that means nothing.
nf=$(strings "$BIN/FlexAIDdS" 2>/dev/null | grep -c FLEXAIDDS_VORONOI_KEYED_JITTER || true)
echo "  FLEXAIDDS_VORONOI_KEYED_JITTER in engine: $nf hit(s)   <- the gate under test"
echo "  FLEXAIDDS_RNG_STREAM_FIX in engine:       $(strings "$BIN/FlexAIDdS" 2>/dev/null | grep -c FLEXAIDDS_RNG_STREAM_FIX || true) hit(s)   (unset in both arms, different consumer)"
if [ "${nf:-0}" -lt 1 ]; then
  echo "FATAL: FLEXAIDDS_VORONOI_KEYED_JITTER is absent from this binary, so the"
  echo "  keyed-jitter branch at Vcontacts.cpp:841 cannot be reached and the ON"
  echo "  arm would be identical to the OFF arm. The probe would be vacuous."
  exit 1
fi

# all-poses coordinate hash (search trajectory), excluding the _INI seed pose
posehash(){
  for f in $(find "$1" -name '*.pdb' 2>/dev/null | grep -v _INI | sort); do
    grep '^\(ATOM\|HETATM\)' "$f" | cut -c31-54
  done | shasum -a256 | cut -c1-16
}
nposes(){ find "$1" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0; }
# elected pose hash, as the engine itself reports it (CSV col 20)
electedsha(){
  local c; c=$(find "$1" -name 'astex_diverse_results.csv' 2>/dev/null | head -1)
  [ -n "${c:-}" ] && awk -F, 'NR==2{print substr($20,1,16)}' "$c" || echo NA
}

# cell <outdir> <target> <threads> <flagmode:off|on> <rep>
cell(){
  local O="$1" T="$2" TH="$3" FM="$4" REP="$5"
  if [ -f "$O/DONE" ]; then return 0; fi
  mkdir -p "$O"
  local t0=$(date +%s) rc=0
  if [ "$FM" = "keyed" ]; then
    ( cd "$O" && env OMP_NUM_THREADS="$TH" \
        FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
        FLEXAIDDS_RESTARTS=3 FLEXAIDDS_NO_SEC=1 FLEXAIDDS_SOFTBETA_ELECTION=1 \
        FLEXAIDDS_VORONOI_KEYED_JITTER=1 \
        "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
          --only-codes "$T" --output "$O" --cache "$CACHE" \
          --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
          --threads 1 --omp-threads "$TH" --ga-population 1000 --ga-generations 300 \
          --job-timeout-seconds 28800 > "$O/stdout.log" 2>&1 ) || rc=$?
  else
    # FLEXAIDDS_VORONOI_KEYED_JITTER DELIBERATELY NOT SET — the legacy
    # thread-local jitter branch (Vcontacts.cpp:850-855) is what runs here.
    ( cd "$O" && env OMP_NUM_THREADS="$TH" \
        FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
        FLEXAIDDS_RESTARTS=3 FLEXAIDDS_NO_SEC=1 FLEXAIDDS_SOFTBETA_ELECTION=1 \
        "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
          --only-codes "$T" --output "$O" --cache "$CACHE" \
          --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
          --threads 1 --omp-threads "$TH" --ga-population 1000 --ga-generations 300 \
          --job-timeout-seconds 28800 > "$O/stdout.log" 2>&1 ) || rc=$?
  fi
  local t1=$(date +%s)
  # did the variable actually reach the child? recorded, not assumed
  local seen; seen=$(grep -c 'VORONOI_KEYED_JITTER\|keyed_jitter' "$O/stdout.log" 2>/dev/null || true)
  local aud; aud=$(grep -haoE 'inconsistent=[0-9]+/[0-9]+ max_delta=[0-9.]+' "$O/stdout.log" 2>/dev/null | tr '\n' ';')
  printf 'rc=%s target=%s threads=%s flag=%s rep=%s poses=%s allhash=%s elected=%s flagseen=%s wall_s=%s audit=%s at=%s\n' \
    "$rc" "$T" "$TH" "$FM" "$REP" "$(nposes "$O")" "$(posehash "$O")" "$(electedsha "$O")" \
    "${seen:-0}" "$((t1-t0))" "${aud:-none}" "$(date -u +%FT%TZ)" > "$O/DONE"
}

i=0
for T in 1P2Y 2D3U 1JD0; do
  for TH in 1 3; do
    for FM in off keyed; do
      for REP in A B; do
        if [ -f "$B/STOP" ]; then echo "  STOP sentinel present — halting at cell boundary"; break 4; fi
        i=$((i+1))
        O="$B/run/${T}_t${TH}_${FM}_${REP}"
        cell "$O" "$T" "$TH" "$FM" "$REP"
        printf "  [%2s/24] %-6s t=%s flag=%-3s rep=%s  %s\n" "$i" "$T" "$TH" "$FM" "$REP" \
          "$(tr '\n' ' ' < "$O/DONE" | cut -c1-110)"
      done
    done
  done
done

echo ""
echo "=== VERDICT: repeat-identity per (target, threads, flag) ==="
printf "  %-6s %-8s %-5s %-10s %-10s %s\n" target threads flag all-poses elected identical
for T in 1P2Y 2D3U 1JD0; do
  for TH in 1 3; do
    for FM in off keyed; do
      A="$B/run/${T}_t${TH}_${FM}_A/DONE"; Bf="$B/run/${T}_t${TH}_${FM}_B/DONE"
      [ -f "$A" ] && [ -f "$Bf" ] || { printf "  %-6s %-8s %-5s %s\n" "$T" "$TH" "$FM" "(incomplete)"; continue; }
      ha=$(grep -o 'allhash=[0-9a-f]*' "$A" | cut -d= -f2); hb=$(grep -o 'allhash=[0-9a-f]*' "$Bf" | cut -d= -f2)
      ea=$(grep -o 'elected=[0-9a-fA-Z]*' "$A" | cut -d= -f2); eb=$(grep -o 'elected=[0-9a-fA-Z]*' "$Bf" | cut -d= -f2)
      sa="DIFFER"; [ "$ha" = "$hb" ] && sa="same"
      se="DIFFER"; [ "$ea" = "$eb" ] && se="same"
      printf "  %-6s %-8s %-5s %-10s %-10s all=%s elected=%s\n" "$T" "$TH" "$FM" "${ha:0:8}/${hb:0:8}" "${ea:0:8}/${eb:0:8}" "$sa" "$se"
    done
  done
done
echo ""
echo "  Read it as: threads=1 flag=off SAME + threads=3 flag=off DIFFER + threads=3"
echo "  flag=on SAME  =>  hypothesis confirmed, the flag is the fix."
echo "  2D3U still DIFFER at threads=1 flag=on  =>  a SECOND source exists."
echo "COMPLETE $(date -u +%FT%TZ)" > "$B/DONE"
echo "=== DETERMINISM PROBE COMPLETE  $(date -u +%FT%TZ) ==="
