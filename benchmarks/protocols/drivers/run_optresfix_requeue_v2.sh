#!/usr/bin/env bash
# =============================================================================
# run_optresfix_requeue_v2.sh -- re-run every FLEXIBLE-arm result invalidated by the
# OptRes capacity-index defect (fixed in bf85c89d).
#
# WHY. build_rotamers wrote each flexible residue's OptRes at FA->optres[MIN_OPTRES-1]
# (CAPACITY) while update_optres reads 0..num_optres-1 (COUNT). reserve_optres()
# decoupled them, so every flexible side chain went unmapped. MEASURED on 1MQ6 seed
# 12345: ligand coord hash 4a12f6e6495babb4 -> 74339c1b01575132, min CF -97.84604 ->
# -240.58793. The ENDPOINT MOVED, so flexible-arm results are void.
#
# RIGID ARMS ARE NOT RE-RUN. num_optres=1 for a rigid cell so there was nothing to
# mismap; the pre-launch gate PROVED the fixed engine reproduces stored rigid cells
# bit-for-bit. That proof is what makes reusing stored rigid arms legitimate.
#
# NO DELETE VERBS. Stale cell dirs are RENAMED ASIDE with a timestamp.
# =============================================================================
set -u

R="${FLEXAIDDS_RESULTS:-$HOME/flexaidds_results}"
case "$R" in /*) ;; *) echo "FATAL: results root not absolute: $R" >&2; exit 2;; esac
TMPDIR="$R/tmp"; mkdir -p "$TMPDIR"; export TMPDIR
case "$TMPDIR" in /*) ;; *) echo "FATAL: TMPDIR not absolute" >&2; exit 2;; esac
: > "$TMPDIR/.writable" || { echo "FATAL: TMPDIR not writable: $TMPDIR" >&2; exit 2; }

CACHE="$R/cache_v2"                 # PARENT of astex_diverse. Passing the dataset dir
                                    # itself is the defect behind the additive-docking
                                    # incident; do not "simplify" this.
SITES="$R/astex85_sites_clean"
BIN_SRC="$(cat "$R/state/prodgate_dir")/bin"
ESHA="$(shasum -a256 "$BIN_SRC/FlexAIDdS" | awk '{print $1}')"
COMMIT="$(cd /Users/lp.more/Projects/FlexAIDdS && GIT_CONFIG_GLOBAL=/dev/null git rev-parse HEAD)"
for p in "$R" "$CACHE" "$SITES" "$BIN_SRC"; do
  [ -d "$p" ] || { echo "FATAL: missing $p" >&2; exit 2; }
done

QUEUE="$R/requeue_$(date -u +%Y%m%d_%H%M%S)"
mkdir -p "$QUEUE"; echo "$QUEUE" > "$R/state/requeue_dir"
STATUS="$R/REQUEUE_STATUS.txt"
POP=1000; GEN=1000; RST=3; WINDOW=3; OMP=3; TMO=28800; NOSEC=0
JOBS="P0_gates P1_flex_3seed P2_nosec_flex P3_inducedfit"

log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$QUEUE/queue.log"; }

status(){
  {
    echo "REQUEUE $(date -u +%FT%TZ)"
    echo "queue   $QUEUE"
    echo "engine  $(echo "$ESHA" | cut -c1-16)   commit $(echo "$COMMIT" | cut -c1-12)"
    for j in $JOBS; do printf '  %-26s %s\n' "$j" "$(cat "$QUEUE/$j.state" 2>/dev/null || echo pending)"; done
    b=$(cat "$QUEUE/current_batch" 2>/dev/null || echo none)
    if [ -d "${b:-/nonexistent}" ]; then
      echo "current $(basename "$b")"
      echo "  cells    $(find "$b/run" -mindepth 3 -name DONE 2>/dev/null | wc -l | tr -d ' ')"
      echo "  unmapped $(find "$b/run" -mindepth 3 -name DONE -exec grep -l 'optres_ok=0' {} + 2>/dev/null | wc -l | tr -d ' ')"
      echo "  rc!=0    $(find "$b/run" -mindepth 3 -name DONE -exec grep -l 'rc=[1-9]' {} + 2>/dev/null | wc -l | tr -d ' ')"
    fi
    echo "load    $(uptime | sed 's/.*load averages*: //')"
    echo "disk    $(df -m /System/Volumes/Data | awk 'NR==2{printf "%.1f GB", $4/1024}')"
  } > "$STATUS.tmp" && mv "$STATUS.tmp" "$STATUS"
}

cell(){  # arm target seed afx shrink outdir
  arm=$1; t=$2; sb=$3; afx=$4; shr=$5; O=$6
  mkdir -p "$O"
  t0=$(date +%s)
  ( cd "$O" && env OMP_NUM_THREADS="$OMP" \
      FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE="$sb" \
      FLEXAIDDS_RESTARTS="$RST" FLEXAIDDS_SCORED_ONLY=1 FLEXAIDDS_NO_SEC="$NOSEC" \
      FLEXAIDDS_AUTOFLEX_MAX="$afx" FLEXAIDDS_AUTOFLEX_METAL_SHRINK="$shr" \
      FLEXAIDDS_OPTRES_DIAG=1 TMPDIR="$TMPDIR" \
      "$BIN_SRC/benchmark_datasets" --benchmark astex_diverse \
        --mode defined-cleft-redock --only-codes "$t" --output "$O" --cache "$CACHE" \
        --threads 1 --omp-threads "$OMP" --ga-population "$POP" --ga-generations "$GEN" \
        --job-timeout-seconds "$TMO" --engine-sha256 "$ESHA" > "$O/run.log" 2>&1 )
  rc=$?; t1=$(date +%s)
  pz=$(find "$O" -name '*.pdb' 2>/dev/null | grep -vce '_INI\.pdb$' || true); pz=${pz:-0}
  rr=$(find "$O" -name 'astex_diverse_results.csv' -exec sh -c 'tail -n +2 "$1" | grep -c .' _ {} \; 2>/dev/null | head -1); rr=${rr:-0}
  sent=$(grep -rIlaF 'SKIPPED' "$O" 2>/dev/null | wc -l | tr -d ' ')
  tmo_hit=$(grep -rIlaE 'timed out|SIGKILL|killed after' "$O" 2>/dev/null | wc -l | tr -d ' ')
  nopt=$(grep -rIhoa 'num_optres=[0-9]*' "$O" 2>/dev/null | head -1 | cut -d= -f2); nopt=${nopt:-0}
  withopt=$(grep -rIhoa 'atoms_with_optres=[0-9]*' "$O" 2>/dev/null | head -1 | cut -d= -f2); withopt=${withopt:-0}
  prot=$(grep -rIhoa 'protein=[0-9]*' "$O" 2>/dev/null | head -1 | cut -d= -f2); prot=${prot:-0}
  ligo=$(grep -rIhoa 'ligand=[0-9]*' "$O" 2>/dev/null | head -1 | cut -d= -f2); ligo=${ligo:-0}
  unmapped=$(grep -rIlaF 'OPTRES] WARNING' "$O" 2>/dev/null | wc -l | tr -d ' ')
  if [ "$afx" -gt 0 ]; then
    if [ "$prot" -gt 0 ] && [ "$unmapped" -eq 0 ]; then optok=1; else optok=0; fi
  else
    if [ "$nopt" -eq 1 ] && [ "$unmapped" -eq 0 ]; then optok=1; else optok=0; fi
  fi
  printf 'rc=%s arm=%s target=%s seed=%s afx=%s shrink=%s engine=%s commit=%s poses=%s result_rows=%s sentinel=%s timeout_hit=%s num_optres=%s atoms_with_optres=%s protein=%s ligand=%s optres_ok=%s nosec=%s pop=%s gen=%s restarts=%s wall_s=%s at=%s\n' \
    "$rc" "$arm" "$t" "$sb" "$afx" "$shr" "$(echo "$ESHA" | cut -c1-16)" "$(echo "$COMMIT" | cut -c1-12)" \
    "$pz" "$rr" "$sent" "$tmo_hit" "$nopt" "$withopt" "$prot" "$ligo" "$optok" \
    "$NOSEC" "$POP" "$GEN" "$RST" "$((t1-t0))" "$(date -u +%FT%TZ)" > "$O/DONE"
}

runjob(){  # jobname arms seeds targets_file afx shrink nosec pop gen
  job=$1; arms=$2; seeds=$3; tf=$4; afx=$5; shr=$6; NOSEC=$7; POP=$8; GEN=$9
  BATCH="$QUEUE/$job"
  mkdir -p "$BATCH/run"; echo "$BATCH" > "$QUEUE/current_batch"
  cp "$tf" "$BATCH/targets.tsv"
  echo "running $(date -u +%FT%TZ)" > "$QUEUE/$job.state"
  log "JOB $job START arms=$arms seeds=$seeds afx=$afx targets=$(grep -c . "$BATCH/targets.tsv")"
  n=0; PIDS=""
  for arm in $arms; do
    for sb in $seeds; do
      while read -r t rest; do
        [ -z "${t:-}" ] && continue
        case "$t" in \#*) continue;; esac
        O="$BATCH/run/${arm}_s${sb}/$t"
        [ -f "$O/DONE" ] && continue
        [ -d "$O" ] && mv "$O" "${O}.stale.$(date -u +%H%M%S)"
        cell "$arm" "$t" "$sb" "$afx" "$shr" "$O" &
        PIDS="$PIDS $!"
        n=$((n+1))
        if [ $(( n % WINDOW )) -eq 0 ]; then
          for p in $PIDS; do wait "$p" || log "WARN child $p nonzero"; done
          PIDS=""; status
        fi
      done < "$BATCH/targets.tsv"
    done
  done
  for p in $PIDS; do wait "$p" || log "WARN child $p nonzero"; done
  tot=$(find "$BATCH/run" -mindepth 3 -name DONE 2>/dev/null | wc -l | tr -d ' ')
  bad=$(find "$BATCH/run" -mindepth 3 -name DONE -exec grep -l 'rc=[1-9]' {} + 2>/dev/null | wc -l | tr -d ' ')
  unm=$(find "$BATCH/run" -mindepth 3 -name DONE -exec grep -l 'optres_ok=0' {} + 2>/dev/null | wc -l | tr -d ' ')
  if [ "$unm" -eq 0 ] && [ "$bad" -eq 0 ] && [ "$tot" -gt 0 ]; then
    echo "VERDICT=PASS cells=$tot" > "$BATCH/VERDICT"
    echo "done $(date -u +%FT%TZ) cells=$tot" > "$QUEUE/$job.state"
    log "JOB $job PASS cells=$tot"; status; return 0
  fi
  echo "VERDICT=FAIL cells=$tot rc_nonzero=$bad optres_unmapped=$unm" > "$BATCH/VERDICT"
  echo "FAILED $(date -u +%FT%TZ) cells=$tot rc=$bad unmapped=$unm" > "$QUEUE/$job.state"
  log "JOB $job FAIL cells=$tot rc_nonzero=$bad unmapped=$unm -- CHAIN HALTED"
  status; return 1
}

# GUARD: P1 and P2 must differ in NO_SEC or P2 is not a control. Checked from the
# invocation lines of this very file so it cannot drift from the comments.
# The seed list is a quoted argument, so it splits into a DIFFERENT number of awk
# fields for P1 (three seeds) than for P2 (one) -- a positional field from the
# raw line reads the wrong column. Strip the trailing '|| exit 1' first; the
# remaining tail is always '<afx> <shrink> <nosec> <pop> <gen>', so nosec is the
# THIRD field from the end and that is stable across both invocations.
nsarg(){ grep -E "^runjob $1 " "$0" | sed 's/ *|| *exit.*//' | awk '{print $(NF-2)}'; }
p1ns=$(nsarg P1_flex_3seed)
p2ns=$(nsarg P2_nosec_flex)
if [ "$p1ns" = "$p2ns" ]; then
  echo "FATAL: P1 and P2 both pass NO_SEC=$p1ns -- P2 would duplicate P1 seed 12345." >&2
  echo "       P1 must be 0 (reproduces astex85_full) and P2 must be 1 (equal-budget)." >&2
  exit 2
fi
echo "GUARD ok: P1 NO_SEC=$p1ns  P2 NO_SEC=$p2ns"
for j in $JOBS; do echo pending > "$QUEUE/$j.state"; done
status
log "QUEUE START engine=$(echo "$ESHA" | cut -c1-16) commit=$(echo "$COMMIT" | cut -c1-12)"

echo "running $(date -u +%FT%TZ)" > "$QUEUE/P0_gates.state"
G="$QUEUE/P0_gates"; mkdir -p "$G"; echo "$G" > "$QUEUE/current_batch"
pass=1
LIGNAME=$(grep -m1 -E '^HETATM|^ATOM' "$CACHE/astex_diverse/1TW6/1TW6.pdb" 2>/dev/null | cut -c18-20)
echo "G1 reference_ligand_probe 1TW6 first_resname=${LIGNAME:-UNKNOWN}" >> "$G/GATES"
echo "G2 engine_sha=$ESHA commit=$COMMIT" >> "$G/GATES"
for spec in "flexON 5 1" "rigidOFF 0 0"; do
  set -- $spec
  nm=$1; a=$2; s=$3
  POP=300; GEN=200; NOSEC=0
  cell "canary" 1MQ6 12345 "$a" "$s" "$G/canary_$nm"
  ok=$(grep -o 'optres_ok=[01]' "$G/canary_$nm/DONE" | cut -d= -f2)
  pr=$(grep -o 'protein=[0-9]*' "$G/canary_$nm/DONE" | cut -d= -f2)
  echo "G3 canary_$nm optres_ok=${ok:-0} protein_atoms_mapped=${pr:-0}" >> "$G/GATES"
  if [ "$nm" = flexON ] && { [ "${ok:-0}" -ne 1 ] || [ "${pr:-0}" -le 0 ]; }; then pass=0; fi
  if [ "$nm" = rigidOFF ] && [ "${ok:-0}" -ne 1 ]; then pass=0; fi
done
if [ "$pass" -eq 1 ]; then
  echo "VERDICT=PASS" > "$G/VERDICT"; echo "done $(date -u +%FT%TZ)" > "$QUEUE/P0_gates.state"
  log "P0 gates PASS"
else
  echo "VERDICT=FAIL" > "$G/VERDICT"; echo "FAILED $(date -u +%FT%TZ)" > "$QUEUE/P0_gates.state"
  log "P0 gates FAIL -- nothing else runs"; status; exit 1
fi

# P1 is the headline claim: the flexibility contrast at three seeds, paired against
#    the STORED rigid arms (proven bit-identical by the pre-launch gate).
# NOSEC=0 -- P1 reproduces astex85_full_20260830_212437, MEASURED to have run with
# early exits ENABLED (0 "early exits DISABLED" banners in its B_shrink_s12345 cells,
# vs 252 in the nosec batch). P1 and P2 differing ONLY in this arg is what makes P2 a
# control rather than a duplicate; they were both 1 in v1, which would have run the
# headline arm in the control configuration and made P2 redundant.
runjob P1_flex_3seed  "B_shrink" "12345 777777 999999" "$R/state/astex85_codes_84.txt"   5 1 0 1000 1000 || exit 1
# P2 is the equal-budget control (NO_SEC=1), single seed, same roster.
runjob P2_nosec_flex  "B_shrink" "12345"               "$R/state/astex85_codes_84.txt"   5 1 1 1000 1000 || exit 1
# P3 is the induced-fit readings. The optres fix makes FLEXAIDDS_SCORED_ONLY write the
#    flexed side chains (measured: 67 atoms over GLN/GLU/TRP/TYR + ligand on 1MQ6), so
#    the displacement metric no longer needs a whole-receptor dump.
runjob P3_inducedfit  "B_shrink" "12345"               "$R/state/inducedfit_targets.txt"  5 1 1 1000 1000 || exit 1

# DELIBERATELY NOT IN THIS QUEUE, and why -- a job that cannot express its own
# experiment is worse than an absent one:
#   cap/ceiling flexible cells  need FLEXAIDDS_MAX_RESULTS=100; cell() has no budget
#                               parameter, so it would silently run at the default 50
#                               and NOT reproduce the cap experiment.
#   CF-vs-FO flexible cells     need the clustering algorithm per arm; cell() has no
#                               algo parameter, so as first drafted P4 was BYTE-FOR-BYTE
#                               the same invocation as P3 under a different job name --
#                               exactly the mislabelling this project already paid for.
# Both are secondary analyses. They get a driver that takes those parameters, not a
# driver that quietly drops them.

echo "COMPLETE $(date -u +%FT%TZ)" > "$QUEUE/DONE"
log "QUEUE COMPLETE"
status
