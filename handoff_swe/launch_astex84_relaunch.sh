#!/bin/bash
# Astex-84 relaunch — post-#403 (deterministic cleft grid), ΔG̃ mode ranking.
#
# Differences from the killed campaign astex84_cleftsort_20260809_124322:
#   - FLEXAIDDS_CLEFT_SORT=1        REMOVED  (imposes a THIRD canonical probe order,
#                                             neither serial nor thread-arrival; changes
#                                             1-thread poses. #403 made it unnecessary.)
#   - FLEXAIDDS_SOFTBETA_ELECTION=1 ADDED    (elect by G̃ = H̃ − T·S̃, S̃ = −Σ p ln p,
#                                             per Morency ISMB/ECCB 2017 3Dsig slide 12.
#                                             Default is min-CF = purely enthalpic.)
#   - binaries RESTAGED from main            (the harness runs its own copied binaries;
#                                             rebuilding the repo does not reach them)
#
# See handoff_swe/HANDOFF_CLAUDE_SCIENCE_benchmark_relaunch_20260809.md
set -euo pipefail

REPO=${FLEXAIDDS_REPO:-/Users/lp.more/Projects/FlexAIDdS}
CACHE=${CACHE:-/Users/lp.more/flexaidds_results/cache_v2}
SITE=${SITE:-/Users/lp.more/flexaidds_results/astex85_sites_clean}
STAMP=$(date +%Y%m%d_%H%M%S)
OUT=${OUT:-/Users/lp.more/flexaidds_results/astex84_dG_${STAMP}}
CODES=${CODES:-}

mkdir -p "$OUT/bin" "$OUT/run"

# ── Build from main ──────────────────────────────────────────────────────────
cd "$REPO"
BRANCH=$(git rev-parse --abbrev-ref HEAD)
SHA=$(git rev-parse --short HEAD)
echo "[$(date +%T)] building $BRANCH @ $SHA"
cmake --build build_tests --target FlexAIDdS benchmark_datasets -j "$(sysctl -n hw.ncpu)"

# ── FAIL-CLOSED PRECONDITIONS ────────────────────────────────────────────────
# The campaign is only interpretable if these hold. Check before burning 5 hours.
fail () { echo "PRECONDITION FAILED: $*" >&2; exit 1; }

grep -q "buckets\[ii\]" "$REPO/LIB/CleftDetector.cpp" \
  || fail "#403 bucket merge absent from LIB/CleftDetector.cpp — cleft grid is nondeterministic"
grep -q "omp critical" "$REPO/LIB/CleftDetector.cpp" \
  && fail "the racy omp-critical probe merge is back in CleftDetector.cpp"
[ "$(grep -c 'election_mode,consensus_count,rank0_demoted' "$REPO/LIB/DatasetRunner.cpp")" -eq 2 ] \
  || fail "#404 CSV header fix missing — results CSV columns will shift by 3"
[ -z "${FLEXAIDDS_CLEFT_SORT:-}" ] \
  || fail "FLEXAIDDS_CLEFT_SORT is set in the environment — unset it"

cp "$REPO/build_tests/FlexAIDdS" "$REPO/build_tests/benchmark_datasets" "$OUT/bin/"
# stage the runtime data files the engine resolves relative to its binary
cp "$REPO"/build_tests/*.dat "$REPO"/build_tests/*.def "$OUT/bin/" 2>/dev/null || true
# fail closed: the engine resolves AMINO.def / NUCLEOTIDES.def next to its binary and
# aborts with exit code 8 before docking if they are absent — which yields a 13-second
# "rc=0" campaign with 84 result rows, 0 poses and a 0.0% rate that looks like a result.
for f in AMINO.def NUCLEOTIDES.def MC_st0r5.2_6.dat; do
  [ -s "$OUT/bin/$f" ] || fail "runtime data file $f missing from $OUT/bin — engine would abort (exit 8)"
done

# ── Campaign environment ─────────────────────────────────────────────────────
export FLEXAID_SEED=12345
export OMP_NUM_THREADS=3
export FLEXAIDDS_ORACLE_SITE_DIR="$SITE"
export FLEXAIDDS_RESTARTS=10
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=fixed
export FLEXAIDDS_BUDGET_SCALE=0
export FLEXAIDDS_SOFTBETA_ELECTION=1     # ΔG̃ = ΔH̃ − TΔS̃ mode ranking
unset  FLEXAIDDS_CLEFT_SORT              # explicit: never set for this protocol

{
  echo "campaign  = $(basename "$OUT")"
  echo "repo_sha  = $SHA  ($BRANCH)"
  echo "engine_md5= $(md5 -q "$OUT/bin/FlexAIDdS")"
  echo "election  = SOFTBETA (G=H-T*S)   cleft_sort = OFF   omp = 3   restarts = 10"
} | tee "$OUT/provenance.txt"

# Input integrity: hash before and after, compare at the end.
shasum -a 256 "$CACHE"/astex_diverse/*/*_apo.pdb "$CACHE"/astex_diverse/*/*_ligand.sdf \
  > "$OUT/inputs_at_launch.txt"

ONLY=()
[ -n "$CODES" ] && ONLY=(--only-codes "$CODES")

echo "[$(date +%T)] launching astex84 | seed 12345 | omp=3 | R=10 | 1000x1000 | ΔG̃ election"
set +e
"$OUT/bin/benchmark_datasets" --benchmark astex --mode defined-cleft-redock \
  "${ONLY[@]}" --threads 2 --omp-threads 3 \
  --ga-generations 1000 --ga-population 1000 \
  --cache "$CACHE" --output "$OUT/run" --job-timeout-seconds 7200 \
  > "$OUT/claim.log" 2>&1
rc=$?
set -e

shasum -a 256 "$CACHE"/astex_diverse/*/*_apo.pdb "$CACHE"/astex_diverse/*/*_ligand.sdf \
  > "$OUT/inputs_at_end.txt"
if diff -q "$OUT/inputs_at_launch.txt" "$OUT/inputs_at_end.txt" >/dev/null; then
  echo "[$(date +%T)] INPUT INTEGRITY: OK"
else
  echo "[$(date +%T)] INPUT INTEGRITY: VIOLATED"
  diff "$OUT/inputs_at_launch.txt" "$OUT/inputs_at_end.txt" || true
fi
echo "[$(date +%T)] benchmark_datasets rc=$rc"
echo "results: $OUT/run"
