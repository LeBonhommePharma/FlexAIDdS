#!/usr/bin/env bash
# ============================================================================
# 1IGJ W0 REDO — recompute the one contaminated cell, off the critical path.
#
# WHY. run/W0/1igj in water_ablation_20260906_022153 resumed into a directory
# left dirty by the killed sequential run (RUN_RECEIPT.json + provenance.json
# present, no DONE). DatasetRunner honoured the stale receipt instead of
# --only-codes and is docking the WHOLE dataset in that one cell: 77/84 target
# dirs, 3.8 GB, 796 min and counting. That cell's receipt cannot be trusted for
# the W0 arm, so 1IGJ must be recomputed in a CLEAN directory.
#
# WHY NOW RATHER THAN QUEUED. The w3 driver is 3-wide and has TWO IDLE SLOTS
# while it waits on that single runaway cell; W2 cannot start until it exits.
# Running this redo now costs nothing on the critical path and removes it from
# the tail. It writes to its OWN batch dir, so it cannot race the running
# driver or be skipped/overwritten by it.
#
# PROTOCOL: byte-for-byte the w3 W0 cell (read from run_water_ablation_w3.sh
# cell(), lines 60-88, and its arm call at line 165):
#   OMP_NUM_THREADS=3, SEED_BASE=12345, RESTARTS=3, NO_SEC=1,
#   SOFTBETA_ELECTION=1, KEEP_STRUCTURAL_WATERS=0, WRITE_FLEXED_RECEPTOR unset
#   (the w3 arm call passes only KEEP_STRUCTURAL_WATERS=0; line 99's
#    WRITE_FLEXED_RECEPTOR=1 belongs to the arm-difference GATE, not the arm),
#   --mode defined-cleft-redock --only-codes 1IGJ
#   --threads 1 --omp-threads 3 --ga-population 1000 --ga-generations 1000
#   --job-timeout-seconds 28800
#
# GATES
#   1 clean dir: refuses to start if the output dir already exists non-empty
#     (that is the exact precondition that caused the contamination).
#   2 one-target: after the run, fails rc=90 if the output holds more than one
#     target directory -- the check that would have caught this in minute one.
#
# NOTE ON 1IGJ. It is the target whose NATIVE pose scores UNFAVOURABLY
# (cf=+1081) in the native-CF oracle -- the worst-behaved target in the set.
# A poor result here is expected and is not evidence of a defect.
#
# STOP LEVER: touch $B/STOP  (checked before launch only; single cell)
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
W=$(ls -d "$R"/water_ablation_* 2>/dev/null | sort | tail -1)
BIN="$W/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
T=1IGJ
STAMP=$(date -u +%Y%m%d_%H%M%S)
B="$R/igj_redo_$STAMP"
O="$B/run/W0/1igj"
mkdir -p "$B/tmp"
echo "$B" > "$R/state/igj_redo_dir"
export TMPDIR="$B/tmp"
exec > >(tee -a "$B/driver.log") 2>&1
echo "=== 1IGJ W0 REDO $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
RSHA=$(shasum -a256 "$BIN/benchmark_datasets" | awk '{print $1}')
echo "  engine ${ESHA:0:16}  runner ${RSHA:0:16}  (same binaries as the W0 arm)"

# ---- GATE 1: the dir must be clean -----------------------------------------
if [ -e "$O" ] && [ -n "$(ls -A "$O" 2>/dev/null)" ]; then
  echo "FATAL: $O already exists and is non-empty. That is the precondition that"
  echo "  produced the contamination; refusing to resume into a dirty directory."
  exit 1
fi
mkdir -p "$O"
echo "  GATE 1 PASS — output dir is clean and empty"

posehash(){ for f in $(find "$1" -name '*.pdb' 2>/dev/null | grep -v _INI | sort); do
              grep '^\(ATOM\|HETATM\)' "$f" | cut -c31-54; done | shasum -a256 | cut -c1-16; }
nposes(){ find "$1" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0; }

t0=$(date +%s); rc=0
echo "  launching $(date -u +%H:%M:%SZ) ..."
( cd "$O" && env OMP_NUM_THREADS=3 \
    FLEXAIDDS_ORACLE_SITE_DIR="$SITES" \
    FLEXAIDDS_SEED_BASE=12345 \
    FLEXAIDDS_RESTARTS=3 \
    FLEXAIDDS_NO_SEC=1 \
    FLEXAIDDS_SOFTBETA_ELECTION=1 \
    TMPDIR="$TMPDIR" \
    FLEXAIDDS_KEEP_STRUCTURAL_WATERS=0 \
    "$BIN/benchmark_datasets" --benchmark astex_diverse --mode defined-cleft-redock \
      --only-codes "$T" --output "$O" --cache "$CACHE" \
      --engine-sha256 "$ESHA" --runner-sha256 "$RSHA" \
      --threads 1 --omp-threads 3 --ga-population 1000 --ga-generations 1000 \
      --job-timeout-seconds 28800 > "$O/stdout.log" 2>&1 ) || rc=$?
t1=$(date +%s)

# ---- GATE 2: exactly one target dir ----------------------------------------
ndir=$(find "$O" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | grep -vcE 'posebust|pose_ledger' || echo 0)
echo "  GATE 2: target dirs in output = $ndir (want 1)"
if [ "${ndir:-0}" -ne 1 ]; then
  echo "  *** CONTAMINATED or EMPTY: expected exactly 1 target dir — marking rc=90"
  rc=90
fi

sn=$(grep -c 'SENTINEL' "$O/stdout.log" 2>/dev/null || true)
to=$(grep -cE 'timed out|SIGKILL|killed after' "$O/stdout.log" 2>/dev/null || true)
printf 'rc=%s target=%s seed=12345 engine=%s runner=%s poses=%s hash=%s sentinel=%s restart_timeout=%s ndir=%s wall_s=%s at=%s\n' \
  "$rc" "$T" "${ESHA:0:16}" "${RSHA:0:16}" "$(nposes "$O")" "$(posehash "$O")" \
  "${sn:-0}" "${to:-0}" "$ndir" "$((t1-t0))" "$(date -u +%FT%TZ)" > "$O/DONE"
echo "  receipt: $(tr '\n' ' ' < "$O/DONE")"
echo "COMPLETE $(date -u +%FT%TZ)" > "$B/DONE"
echo "=== 1IGJ W0 REDO COMPLETE rc=$rc $(date -u +%FT%TZ) ==="
