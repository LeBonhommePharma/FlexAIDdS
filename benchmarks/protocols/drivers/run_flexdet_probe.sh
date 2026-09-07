#!/usr/bin/env bash
# ============================================================================
# FLEX DETERMINISM PROBE — the decisive test of the amplification mechanism
#
# WHAT IS MEASURED SO FAR (contention_flex_20260906_091252, 1P2Y, omp=3, 4 runs)
#   INI poses            IDENTICAL in all 4 runs  -> setup/seeding deterministic
#   restart 0            jaccard 0.976, 40/41 shared, ONE cluster differs,
#                        pose counts 41/40/41/41    -> a TINY perturbation
#   restart 1            jaccard 0.000, ZERO shared -> COMPLETELY different
#   restart 2            jaccard 0.138
#   best CF              -254.37 / -230.86 / -231.03 / -231.41
#   restarts run SEQUENTIALLY in one process:
#     "[PARALLEL-RESTART] launch window = 1 concurrent restart(s)"
#   concurrency ruled out: ISO and CON diverge identically, and the RIGID path
#   is 15/15 + 12/12 reproducible across omp 1 and 3.
#
# THE HYPOTHESIS THIS TESTS
#   A near-degenerate tie in restart 0 resolves differently run-to-run, which
#   (a) changes one cluster and (b) consumes a DIFFERENT NUMBER OF RNG DRAWS.
#   Because restarts share one process and one RNG stream, restart 1 then
#   starts from a shifted stream and diverges completely. Divergence is
#   therefore AMPLIFIED, not independently random per restart.
#   The candidate symmetry-breaker is the Voronoi degeneracy jitter
#   (Vcontacts.cpp:850-855), whose draws come from a THREAD_LOCAL RNG, so which
#   thread handles which degeneracy -- and hence the draw sequence -- depends on
#   OMP scheduling. 1P2Y is the target whose exact-pose audit reported
#   2256/2269 restarts inconsistent, i.e. it hits that branch constantly.
#
# WHY THIS IS NOT A REPEAT OF THE EARLIER THREAD TEST. That one ran the RIGID
# path (autoflex_max=0) and found thread count inert -- correctly, for rigid.
# This runs the campaign's FLEXIBLE configuration.
#
# ARMS (1P2Y only; 1OQ5/1OWE are already reproducible and serve as controls)
#   A  omp=1  jitter unset   -> if REPRODUCIBLE, the mechanism is thread
#                               scheduling in the jitter path, and single-thread
#                               is a workaround
#   B  omp=3  jitter=keyed   -> if REPRODUCIBLE, the schedule-independent keyed
#                               jitter is the FIX (note: it is also a physics
#                               change; it moved the elected pose on 2/3 targets)
#   The existing omp=3 / jitter-unset runs (4 of them, all divergent) are the
#   negative control; this probe does not re-run them.
#
# PREDICTION TABLE, written before the run:
#   A same, B same  -> jitter path confirmed; keyed jitter is the fix
#   A same, B diff  -> thread scheduling matters but keyed jitter is not enough
#   A diff, B same  -> keyed jitter fixes it for a reason other than threading
#   A diff, B diff  -> neither; the symmetry-breaker is elsewhere (next suspect:
#                      the clustering tie-break itself, which needs source work)
#
# GATES
#   1 engine identity: refuses unless the binary is the campaign's d7916d13
#   2 flexibility engaged: emitted config must carry autoflex_max=5
#   3 per-restart reporting: the verdict is computed PER RESTART (r0/r1/r2), not
#     only on the pooled pose set, because r0's single-cluster difference is the
#     signal and pooling hides it
#
# STOP LEVER: touch $B/STOP
# ============================================================================
set -u
R=/Users/lp.more/flexaidds_results
BIN="$R/prodgate_040617/bin"
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
T=1P2Y
GEN=1000; POP=1000; RST=3; TMO=28800
AFX=5; SHR=1
STAMP=$(date -u +%Y%m%d_%H%M%S)
B="$R/flexdet_$STAMP"
mkdir -p "$B/run" "$B/tmp"
echo "$B" > "$R/state/flexdet_dir"
export TMPDIR="$B/tmp"
exec > >(tee -a "$B/driver.log") 2>&1
echo "=== FLEX DETERMINISM PROBE $(date -u +%FT%TZ)  batch=$(basename "$B") ==="

ESHA=$(shasum -a256 "$BIN/FlexAIDdS" | awk '{print $1}')
if [ "${ESHA:0:16}" != "d7916d136f89aa70" ]; then
  echo "FATAL: engine ${ESHA:0:16} != d7916d136f89aa70 (campaign). Wrong binary."; exit 1
fi
echo "  GATE 1 PASS — campaign binary ${ESHA:0:16}"

posehash(){ for f in $(find "$1" -name '*.pdb' 2>/dev/null | grep -v _INI | sort); do
              grep '^\(ATOM\|HETATM\)' "$f" | cut -c31-54; done | shasum -a256 | cut -c1-16; }
nposes(){ find "$1" -name '*.pdb' 2>/dev/null | grep -vc _INI || echo 0; }
electedsha(){ local c; c=$(find "$1" -name 'astex_diverse_results.csv' 2>/dev/null | head -1)
              [ -n "${c:-}" ] && awk -F, 'NR==2{print substr($20,1,16)}' "$c" || echo NA; }
cfg_afx(){ local f; f=$(find "$1" -name 'dock_config.json' 2>/dev/null | head -1)
           [ -n "${f:-}" ] && python3 -c "
import json
try: print(json.load(open('$f')).get('flexibility',{}).get('autoflex_max','MISSING'))
except Exception: print('ERR')" || echo NOCFG; }

# cell <outdir> <omp> <jittermode>
cell(){
  local O="$1" OMP="$2" JIT="$3"
  [ -f "$O/DONE" ] && return 0
  mkdir -p "$O"
  local t0=$(date +%s) rc=0
  if [ "$JIT" = "keyed" ]; then
    ( cd "$O" && env OMP_NUM_THREADS="$OMP" \
        FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
        FLEXAIDDS_RESTARTS="$RST" FLEXAIDDS_SCORED_ONLY=1 FLEXAIDDS_NO_SEC=1 \
        FLEXAIDDS_AUTOFLEX_MAX="$AFX" FLEXAIDDS_AUTOFLEX_METAL_SHRINK="$SHR" \
        FLEXAIDDS_OPTRES_DIAG=1 FLEXAIDDS_VORONOI_KEYED_JITTER=1 TMPDIR="$TMPDIR" \
        "$BIN/benchmark_datasets" --benchmark astex_diverse \
          --mode defined-cleft-redock --only-codes "$T" --output "$O" --cache "$CACHE" \
          --threads 1 --omp-threads "$OMP" --ga-population "$POP" --ga-generations "$GEN" \
          --job-timeout-seconds "$TMO" --engine-sha256 "$ESHA" > "$O/run.log" 2>&1 ) || rc=$?
  else
    # FLEXAIDDS_VORONOI_KEYED_JITTER DELIBERATELY NOT SET (legacy thread-local branch)
    ( cd "$O" && env OMP_NUM_THREADS="$OMP" \
        FLEXAIDDS_ORACLE_SITE_DIR="$SITES" FLEXAIDDS_SEED_BASE=12345 \
        FLEXAIDDS_RESTARTS="$RST" FLEXAIDDS_SCORED_ONLY=1 FLEXAIDDS_NO_SEC=1 \
        FLEXAIDDS_AUTOFLEX_MAX="$AFX" FLEXAIDDS_AUTOFLEX_METAL_SHRINK="$SHR" \
        FLEXAIDDS_OPTRES_DIAG=1 TMPDIR="$TMPDIR" \
        "$BIN/benchmark_datasets" --benchmark astex_diverse \
          --mode defined-cleft-redock --only-codes "$T" --output "$O" --cache "$CACHE" \
          --threads 1 --omp-threads "$OMP" --ga-population "$POP" --ga-generations "$GEN" \
          --job-timeout-seconds "$TMO" --engine-sha256 "$ESHA" > "$O/run.log" 2>&1 ) || rc=$?
  fi
  local t1=$(date +%s)
  printf 'rc=%s omp=%s jitter=%s poses=%s allhash=%s elected=%s afx_cfg=%s wall_s=%s load1=%s at=%s\n' \
    "$rc" "$OMP" "$JIT" "$(nposes "$O")" "$(posehash "$O")" "$(electedsha "$O")" \
    "$(cfg_afx "$O")" "$((t1-t0))" \
    "$(uptime | sed 's/.*load average[s]*: //' | awk '{print $1}' | tr -d ,)" \
    "$(date -u +%FT%TZ)" > "$O/DONE"
}

first=1
for spec in "omp1_off 1 off" "omp3_keyed 3 keyed"; do
  set -- $spec
  LBL=$1; OMP=$2; JIT=$3
  for REP in A B; do
    [ -f "$B/STOP" ] && { echo "  STOP sentinel"; break 2; }
    O="$B/run/${LBL}_${REP}"
    echo "  --- ${LBL}_${REP}: omp=$OMP jitter=$JIT  $(date -u +%H:%M:%SZ) ---"
    cell "$O" "$OMP" "$JIT"
    echo "    $(tr '\n' ' ' < "$O/DONE" | cut -c1-118)"
    if [ "$first" -eq 1 ]; then
      a=$(grep -o 'afx_cfg=[A-Za-z0-9]*' "$O/DONE" | cut -d= -f2)
      echo "    GATE 2: autoflex_max=$a (want 5)"
      [ "$a" != "5" ] && { echo "    FATAL: rigid path again — probe would be vacuous."; exit 2; }
      echo "    GATE 2 PASS"
      first=0
    fi
  done
done

echo ""
echo "=== VERDICT: PER-RESTART, because r0's single-cluster difference is the signal ==="
python3 - "$B" <<'PY'
import glob, os, re, hashlib, sys, itertools
B = sys.argv[1]
def sha(f):
    h = hashlib.sha256()
    for l in open(f, errors="ignore"):
        if l.startswith(("ATOM","HETATM")): h.update(l[30:54].encode())
    return h.hexdigest()[:12]
def cf(f):
    for l in open(f, errors="ignore"):
        if l.startswith("REMARK"):
            m = re.search(r"CF\s*=\s*(-?[\d.]+)", l)
            if m: return float(m.group(1))
    return None
def per_restart(cell):
    out = {}
    for f in glob.glob(f"{B}/run/{cell}/1P2Y/**/*.pdb", recursive=True):
        if "_INI" in os.path.basename(f): continue
        m = re.search(r"/r(\d+)/", f)
        out.setdefault(int(m.group(1)) if m else 0, {})[sha(f)] = cf(f)
    return out
def inihash(cell):
    ini = sorted(glob.glob(f"{B}/run/{cell}/**/*_INI.pdb", recursive=True))
    return hashlib.sha256("".join(sha(x) for x in ini).encode()).hexdigest()[:12] if ini else "-"
print(f"  {'arm':<12} {'INI':<14} {'restart':<8} {'nA':>4} {'nB':>4} {'shared':>7} {'jaccard':>8} verdict")
for lbl in ("omp1_off","omp3_keyed"):
    ca, cb = f"{lbl}_A", f"{lbl}_B"
    if not (os.path.exists(f"{B}/run/{ca}/DONE") and os.path.exists(f"{B}/run/{cb}/DONE")):
        print(f"  {lbl:<12} incomplete"); continue
    A, Bd = per_restart(ca), per_restart(cb)
    ini_same = inihash(ca) == inihash(cb)
    allsame = True
    for r in sorted(set(A) | set(Bd)):
        sa, sb = set(A.get(r, {})), set(Bd.get(r, {}))
        j = len(sa & sb)/max(1, len(sa | sb))
        if j != 1.0: allsame = False
        print(f"  {lbl:<12} {inihash(ca)+('=' if ini_same else '!'):<14} r{r:<7} "
              f"{len(sa):>4} {len(sb):>4} {len(sa&sb):>7} {j:>8.3f} "
              f"{'IDENTICAL' if j == 1.0 else 'DIVERGED'}")
    bestA = min([v for d in A.values() for v in d.values() if v is not None] or [float('nan')])
    bestB = min([v for d in Bd.values() for v in d.values() if v is not None] or [float('nan')])
    print(f"  {lbl:<12} POOLED: {'REPRODUCIBLE' if allsame else 'NON-REPRODUCIBLE'}"
          f"   best CF {bestA:.4f} vs {bestB:.4f}")
print()
print("  reference (omp=3, jitter unset, 4 runs): r0 j=0.976, r1 j=0.000, r2 j=0.138 -> NON-REPRODUCIBLE")
PY
echo "COMPLETE $(date -u +%FT%TZ)" > "$B/DONE"
echo "=== FLEX DETERMINISM PROBE COMPLETE $(date -u +%FT%TZ) ==="
