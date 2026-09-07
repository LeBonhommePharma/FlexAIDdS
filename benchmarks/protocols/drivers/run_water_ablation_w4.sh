#!/usr/bin/env bash
# ============================================================================
# WATER ABLATION, 3-WIDE — a NEW FILE, not an edit of run_water_ablation.sh
# (that one may still be executing; never edit a script a process is reading).
#
# WHY THIS EXISTS. The original driver runs cells strictly sequentially:
#   for ARM; do while read T; do cell "$O" "$T"; done < ROSTER; done
# with --omp-threads 3, so ~3 of 11 cores are busy and 168 cells take ~40 h.
# MEASURED: per-cell wall is IDENTICAL to the requeue's (median ratio 1.00 on
# the 6 overlapping targets), and the requeue ran 3 cells CONCURRENTLY. So
# one-at-a-time buys nothing and costs ~27 h. 3-wide -> ~13 h for the same work.
#
# RESUMES BY DESIGN. cell() returns early when $O/DONE exists, so pointing this
# at the existing batch keeps every completed cell. Pass the batch dir as $1:
#   bash run_water_ablation_w3.sh /Users/lp.more/flexaidds_results/water_ablation_20260906_022153
#
# NO acquire_box_slot CALL. This is launched by hand, so it must not take or
# release a box slot -- in particular it must NEVER write box_queue/*.done,
# which is what would release queue3's 100-cell run.
#
# STOP LEVER: touch $B/STOP -> the loop stops dispatching at the next window
# boundary. In-flight cells finish and write real receipts.
#
# GATE ADDED: FLEXAIDDS_KEEP_STRUCTURAL_WATERS appears 0 times in the FlexAIDdS
# engine binary and once in benchmark_datasets, i.e. the runner is expected to
# translate it into the emitted dock_config. That is the same shape as the
# NO_SEC defect, so this driver PROVES the arms differ before spending 13 h:
# it diffs W0 vs W2 emitted config and receptor atom counts on one target and
# refuses to run the arms if they are identical.
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
B="${1:-$(ls -d "$R"/water_ablation_* 2>/dev/null | sort | tail -1)}"
[ -d "$B" ] || { echo "FATAL: no batch dir at '$B'"; exit 1; }
BIN="$B/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
ROSTER="$R/state/astex85_codes_84.txt"
WINDOW=3
mkdir -p "$B/run" "$B/tmp" "$B/gates"
export TMPDIR="$B/tmp"

exec > >(tee -a "$B/driver_w3.log") 2>&1
echo "=== WATER ABLATION 3-WIDE  $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
RSHA=$(shasum -a256 "$BIN/benchmark_datasets" | awk '{print $1}')
echo "  engine $ESHA"
echo "  runner $RSHA"
echo "  resuming: $(find "$B/run" -name DONE 2>/dev/null | wc -l | tr -d ' ') cells already complete"

posehash(){
  for f in $(find "$1" -name '*.pdb' 2>/dev/null | grep -v _INI | sort); do
    grep '^\(ATOM\|HETATM\)' "$f" | cut -c31-54
  done | shasum -a256 | cut -c1-16
}
nposes(){ find "$1" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0; }

# cell <outdir> <target> [extra env assignments...]
# ============================================================================
# RESUME GUARD (added after a MEASURED failure, 2026-09-06).
#
# WHAT WENT WRONG. run_water_ablation.sh was killed mid-cell, leaving
# run/W0/1igj with RUN_RECEIPT.json + provenance.json already written but NO
# DONE. cell() skips only when DONE exists, so w3 re-ran INTO that dirty
# directory. DatasetRunner saw an existing receipt in --output and resumed the
# WHOLE DATASET instead of honouring --only-codes: that one cell accumulated 77
# target directories and 9,000+ poses over 13 h, wasting the box and blocking
# the W2 arm behind it.
#
# TWO GATES, because one of them would have caught it in the first minute:
#   PRE  : any cell dir without DONE is RENAMED ASIDE with a timestamp before
#          the run (renamed, never deleted -- these dirs carry provenance).
#   POST : a finished cell whose output holds more than ONE target directory is
#          marked contaminated (rc=90) instead of being trusted. --only-codes
#          names exactly one target, so >1 dir means the filter was not honoured.
# ============================================================================
prescrub(){
  local TS; TS=$(date -u +%Y%m%d_%H%M%S)
  local n=0
  for ARM in W0 W2; do
    [ -d "$B/run/$ARM" ] || continue
    for d in "$B/run/$ARM"/*; do
      [ -d "$d" ] || continue
      case "$(basename "$d")" in *.dirty_*) continue ;; esac
      if [ ! -f "$d/DONE" ]; then
        np=$(find "$d" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0)
        nd=$(find "$d" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
        mv "$d" "${d}.dirty_$TS"
        echo "  PRESCRUB renamed aside (no DONE, $np poses, $nd target dirs): $ARM/$(basename "$d")"
        n=$((n+1))
      fi
    done
  done
  echo "  PRESCRUB: $n incomplete cell dir(s) renamed aside; resume is now clean"
}

cell(){
  local O="$1" T="$2"; shift 2
  if [ -f "$O/DONE" ]; then return 0; fi
  mkdir -p "$O"
  local t0=$(date +%s) rc=0
  ( cd "$O" && env OMP_NUM_THREADS=3 \
      FLEXAIDDS_ORACLE_SITE_DIR="$SITES" \
      FLEXAIDDS_SEED_BASE=12345 \
      FLEXAIDDS_RESTARTS=3 \
      FLEXAIDDS_NO_SEC=1 \
      FLEXAIDDS_SOFTBETA_ELECTION=1 \
      TMPDIR="$TMPDIR" \
      "$@" \
      "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
        --only-codes "$T" --output "$O" --cache "$CACHE" \
        --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
        --threads 1 --omp-threads 3 --ga-population 1000 --ga-generations 1000 \
        --job-timeout-seconds 28800 > "$O/stdout.log" 2>&1 ) || rc=$?
  local t1=$(date +%s)
  # POST GATE: --only-codes names ONE target, so >1 target dir means the filter
  # was not honoured and this cell's output is not what it claims to be.
  local ndir; ndir=$(find "$O" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | grep -vcE 'posebust|pose_ledger' || echo 0)
  if [ "${ndir:-0}" -gt 1 ]; then
    echo "  *** CONTAMINATED: $O holds $ndir target dirs (expected 1) — marking rc=90"
    rc=90
  fi
  local sn; sn=$(grep -c 'SENTINEL' "$O/stdout.log" 2>/dev/null || true)
  local vd; vd=$(grep -c 'VCF-DIAG' "$O/stdout.log" 2>/dev/null || true)
  local to; to=$(grep -cE 'timed out|SIGKILL|killed after' "$O/stdout.log" 2>/dev/null || true)
  printf 'rc=%s target=%s seed=12345 engine=%s runner=%s poses=%s hash=%s sentinel=%s vcf_diag=%s restart_timeout=%s wall_s=%s at=%s\n' \
    "$rc" "$T" "${ESHA:0:16}" "${RSHA:0:16}" "$(nposes "$O")" "$(posehash "$O")" \
    "${sn:-0}" "${vd:-0}" "${to:-0}" "$((t1-t0))" "$(date -u +%FT%TZ)" > "$O/DONE"
}

# ---- GATE: do the two arms actually DIFFER at the engine? -------------------
# Refuses to spend 13 h on two arms that turn out to be one arm (the NO_SEC
# lesson). Compares emitted config + receptor HOH counts on a single target.
G="$B/gates/w3_armdiff"
mkdir -p "$G"
echo ""
echo "=== ARM-DIFFERENCE GATE (1JD0, config + water counts) ==="
for ARM in W0 W2; do
  O="$G/$ARM"; mkdir -p "$O"
  if [ "$ARM" = "W0" ]; then
    ( cd "$O" && env OMP_NUM_THREADS=2 FLEXAIDDS_ORACLE_SITE_DIR="$SITES" \
        FLEXAIDDS_SEED_BASE=12345 FLEXAIDDS_RESTARTS=1 FLEXAIDDS_NO_SEC=1 \
        FLEXAIDDS_KEEP_STRUCTURAL_WATERS=0 FLEXAIDDS_WRITE_FLEXED_RECEPTOR=1 \
        "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
          --only-codes 1JD0 --output "$O" --cache "$CACHE" \
          --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
          --threads 1 --omp-threads 2 --ga-population 200 --ga-generations 50 \
          --job-timeout-seconds 3600 > "$O/stdout.log" 2>&1 ) || true
  else
    ( cd "$O" && env OMP_NUM_THREADS=2 FLEXAIDDS_ORACLE_SITE_DIR="$SITES" \
        FLEXAIDDS_SEED_BASE=12345 FLEXAIDDS_RESTARTS=1 FLEXAIDDS_NO_SEC=1 \
        FLEXAIDDS_WRITE_FLEXED_RECEPTOR=1 \
        "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
          --only-codes 1JD0 --output "$O" --cache "$CACHE" \
          --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
          --threads 1 --omp-threads 2 --ga-population 200 --ga-generations 50 \
          --job-timeout-seconds 3600 > "$O/stdout.log" 2>&1 ) || true
  fi
  cfg=$(find "$O" -name 'dock_config.json' | head -1)
  hoh=$(find "$O" -name '*flexed*receptor*.pdb' -o -name '*receptor*.pdb' 2>/dev/null | head -1)
  nh=$([ -n "${hoh:-}" ] && grep -c 'HOH' "$hoh" 2>/dev/null || echo NA)
  printf "  %-3s cfg=%s  receptor_HOH=%s  protein_block=%s\n" "$ARM" \
    "$([ -n "${cfg:-}" ] && shasum -a256 "$cfg" | cut -c1-12 || echo none)" "$nh" \
    "$([ -n "${cfg:-}" ] && grep -c '"protein"' "$cfg" 2>/dev/null || echo NA)"
  echo "$nh" > "$G/$ARM.hoh"
done
h0=$(cat "$G/W0.hoh" 2>/dev/null || echo NA); h2=$(cat "$G/W2.hoh" 2>/dev/null || echo NA)
echo "  W0 HOH=$h0   W2 HOH=$h2"
if [ "$h0" = "$h2" ]; then
  echo "FATAL: the two arms load the SAME number of waters -> they are ONE arm."
  echo "  Refusing to spend ~13 h on a duplicate. Fix the water plumbing first."
  echo "  (This is the NO_SEC failure mode: a flag that never reaches the engine.)"
  exit 1
fi
echo "  GATE PASS — the arms differ at the engine."

# ---- RESUME GUARD, before any cell runs -----------------------------------
prescrub

# ---- THE ARMS, 3-WIDE ------------------------------------------------------
for ARM in W0 W2; do
  echo ""
  echo "=== ARM $ARM  $(date -u +%FT%TZ) ==="
  # SLIDING WINDOW, bash-3.2-SAFE.
  #
  # An earlier version of this loop used `wait -n`, which is a bash>=4.3
  # builtin option. /bin/bash on this box is 3.2.57 and is the ONLY bash
  # present, so `wait -n` returns "invalid option" (rc=2), the `|| wait`
  # fallback then blocked on ALL children while `running` was decremented by
  # one, pinning `running` at WINDOW-1 and making every dispatch after the
  # first triple join immediately -- STRICTLY SERIAL, i.e. exactly the 40 h
  # pathology this file exists to remove. Found by review, then reproduced:
  #   9 jobs x sleep 2, WINDOW=3   ->  wait -n form 14 s, this form 6 s
  # `bash -n` cannot catch an unsupported builtin option, so this idiom is
  # verified BEHAVIOURALLY (timing) and not by syntax check alone.
  #
  # `jobs -pr` (running background jobs of the current shell) works in 3.2.
  # This is a SLIDING window, not a barrier: a slot refills as soon as any one
  # cell finishes. The requeue used a barrier (groups of WINDOW), which is why
  # a slow target like 1UML held two idle slots for ~30 min. `done < "$ROSTER"`
  # is a redirect rather than a pipe, so the loop body stays in this shell and
  # `jobs`/`$!` refer to these children.
  i=0
  while read -r T; do
    [ -z "$T" ] && continue
    if [ -f "$B/STOP" ]; then echo "  STOP sentinel — no further dispatch"; break 2; fi
    i=$((i+1))
    O="$B/run/$ARM/$T"
    if [ -f "$O/DONE" ]; then printf "  [%s %3s/84] %-6s already complete\n" "$ARM" "$i" "$T"; continue; fi
    while [ "$(jobs -pr 2>/dev/null | wc -l | tr -d ' ')" -ge "$WINDOW" ]; do sleep 5; done
    if [ "$ARM" = "W0" ]; then
      cell "$O" "$T" FLEXAIDDS_KEEP_STRUCTURAL_WATERS=0 &
    else
      cell "$O" "$T" &
    fi
    printf "  [%s %3s/84] %-6s dispatched pid=%s  live=%s/%s\n" "$ARM" "$i" "$T" "$!" \
      "$(jobs -pr 2>/dev/null | wc -l | tr -d ' ')" "$WINDOW"
  done < "$ROSTER"
  wait
  echo "  ARM $ARM complete: $(find "$B/run/$ARM" -name DONE 2>/dev/null | wc -l | tr -d ' ')/84"
done

echo ""
echo "=== GUARDS ACROSS ALL RECEIPTS ==="
for k in 'rc=[1-9]' 'poses=0 ' 'sentinel=[1-9]' 'restart_timeout=[1-9]'; do
  printf "  %-24s %s\n" "$k" "$(find "$B/run" -name DONE -exec grep -l "$k" {} + 2>/dev/null | wc -l | tr -d ' ')"
done
echo "  cells: $(find "$B/run" -name DONE 2>/dev/null | wc -l | tr -d ' ')/168"
echo "COMPLETE_W3 $(date -u +%FT%TZ)" > "$B/DONE_W3"
echo "=== WATER ABLATION 3-WIDE COMPLETE  $(date -u +%FT%TZ) ==="
