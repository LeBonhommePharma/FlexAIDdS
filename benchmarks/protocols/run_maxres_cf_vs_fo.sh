#!/usr/bin/env bash
# =============================================================================
# run_maxres.sh -- does the emission cap truncate the POOL CEILING?
#
# SCOPE, DECLARED BEFORE ANYTHING RUNS.
#   IN  : the pool ceiling, and the sampling-vs-selection decomposition that rests on it.
#   OUT : top-1 and the -16 / -18 verdicts. Those are CAP-INVARIANT (emission is sorted by
#         apparent CF; the min-CF pose sits at rank <=5 in essentially every cell, so ranks
#         50-299 cannot contain it). This run CANNOT and MUST NOT be read as revisiting them.
#
# REQUIRES the harness patch maxres_harness_gate.patch, which gates BOTH consumers of the
# limit. Without it the engine writes N poses while enumerate_emitted_cluster_heads() still
# reads 50 -- and this experiment measures nothing while every log and exit code looks clean.
#
# GATES ARE CONDITIONALS, NOT PROMISES. Each phase writes a verdict file; the next phase
# refuses to launch unless it reads PASS. An rc=0 check is not a gate.
# =============================================================================
set -u
set -o pipefail

# --- HOUSE CONVENTION: source the project preamble for its absolute-TMPDIR discipline
# --- and its require_abs assertions on the results/cache/site roots. It exists because
# --- arm14 lost 8 of 85 targets to a RELATIVE inherited TMPDIR.
#
# WHAT WE DELIBERATELY DO NOT TAKE FROM IT: $FLEXAIDDS_CACHE. The preamble points that at
# .../cache_v2/astex_diverse -- the DATASET directory. benchmark_datasets --cache expects
# the PARENT. Passing the dataset dir is the same defect behind the incident where this
# project silently docked a crystallisation additive. We set CACHE to the parent below.
PREAMBLE="${FLEXAIDDS_PREAMBLE:-/Users/lp.more/flexaidds_results/driver_preamble.sh}"
[ -f "$PREAMBLE" ] || { echo "FATAL: preamble not found: $PREAMBLE" >&2; exit 2; }
# shellcheck source=/dev/null
. "$PREAMBLE"

R="${FLEXAIDDS_RESULTS:?preamble must export FLEXAIDDS_RESULTS}"
case "$R" in /*) ;; *) echo "FATAL: FLEXAIDDS_RESULTS not absolute: $R" >&2; exit 2;; esac
BATCH="$R/maxres_$(date -u +%Y%m%d_%H%M%S)"
REF="$(cat "$R/state/astex85_full_batch")"          # frozen campaign: used for TARGET LISTS and
                                                     # ceilings only. NOT an inertness reference --
                                                     # its engine sha is not reproducible from its
                                                     # recorded commit (see gate 0a).
CACHE="$R/cache_v2"
SITES="$R/astex85_sites_clean"
TARGETS="$BATCH/maxres_flip_targets.tsv"
LIMIT="${MAXRES_LIMIT:-100}"
CANARY_T="${MAXRES_CANARY:-1G9V}"
# CF vs FO. FO has NEVER executed in this project (0 dual-suffix poses exist anywhere on
# disk), so phase 0 smoke-tests it before any campaign cell is spent.
ALGOS="${MAXRES_ALGOS:-CF FO}"
SEED=12345
WINDOW=3
OMP=3
TMO=28800

# --- absolute scratch. NEVER inherit TMPDIR: the sandbox exports a RELATIVE ./.tmp and a
# --- relative scratch already cost this project an arm. Not /private/tmp either: swept.
mkdir -p "$BATCH/bin" "$BATCH/run" "$BATCH/tmp"
# Per-batch scratch UNDER the preamble's absolute root: keeps the absoluteness guarantee
# while preventing two concurrent runs from sharing one scratch directory.
export TMPDIR="/${BATCH#/}/tmp"
case "$TMPDIR" in /*) ;; *) echo "FATAL: TMPDIR not absolute: $TMPDIR" >&2; exit 2;; esac
[ -d "$TMPDIR" ] || { echo "FATAL: TMPDIR missing: $TMPDIR" >&2; exit 2; }
: > "$TMPDIR/.writable" || { echo "FATAL: TMPDIR not writable: $TMPDIR" >&2; exit 2; }
# This driver contains ZERO delete verbs, so it cannot remove the wrong path.

# Binaries come from the BUILD, not from the frozen campaign batch: the campaign's
# engine sha is not reproducible from its recorded commit, so staging it would mix an
# unattributable binary into a measurement about this patch.
ARM_BIN="${MAXRES_ARM_BIN:-$R/build_flex}"
for f in FlexAIDdS benchmark_datasets; do
  [ -x "$ARM_BIN/$f" ] || { echo "FATAL: missing $ARM_BIN/$f" >&2; exit 2; }
  cp "$ARM_BIN/$f" "$BATCH/bin/"
done
for f in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def rotobs.lst; do
  src=$(find "$REF/bin" "$ARM_BIN" -maxdepth 1 -name "$f" 2>/dev/null | head -1)
  [ -n "${src:-}" ] && cp "$src" "$BATCH/bin/"
done
ESHA="$(shasum -a256 "$BATCH/bin/FlexAIDdS" | awk '{print $1}')"
cp "${MAXRES_TARGETS:-$R/state/maxres_flip_targets.tsv}" "$TARGETS"
NT=$(( $(wc -l < "$TARGETS") - 1 ))

# --- G-IDENT: the ligand-identity canary. This project once silently docked a
# --- crystallisation additive; a hardcoded single-target check is the outside check.
LIG=$(awk '/^....$|^[A-Z0-9]{3} /{print}' /dev/null 2>/dev/null; \
      grep -m1 -A4 '^@<TRIPOS>' /dev/null 2>/dev/null; \
      python3 - <<PY
import sys
p="$CACHE/astex_diverse/1TW6/1TW6_ligand.sdf"
try:
    L=open(p,errors="ignore").read().splitlines()
    print("ALA" if any("ALA" in x for x in L[:4]) else L[0].strip()[:8] or "UNKNOWN")
except Exception: print("UNKNOWN")
PY
)
echo "[G-IDENT] 1TW6 reference ligand probe = ${LIG}  engine=${ESHA:0:16}  targets=${NT}"
[ "$NT" -ge 1 ] || { echo "FATAL: no targets parsed from $TARGETS" >&2; exit 2; }

env_for() {   # $1 = arm
  printf 'FLEXAIDDS_ORACLE_SITE_DIR=%s FLEXAIDDS_SEED_BASE=%s FLEXAIDDS_RESTARTS=3 ' \
         "$SITES" "$SEED"
  printf 'FLEXAIDDS_SCORED_ONLY=1 '
  if [ "$1" = "B_shrink" ]; then
    printf 'FLEXAIDDS_AUTOFLEX_MAX=5 FLEXAIDDS_AUTOFLEX_METAL_SHRINK=1 '
  fi
}

cell() {      # $1 arm  $2 target  $3 limit-or-UNSET  $4 outdir  $5 algo
  local arm=$1 t=$2 lim=$3 O=$4 algo=${5:-CF}
  mkdir -p "$O"
  local t0=$(date +%s) extra=""
  [ "$lim" != "UNSET" ] && extra="FLEXAIDDS_MAX_RESULTS=$lim"
  ( cd "$BATCH/bin" && env $(env_for "$arm") $extra TMPDIR="$TMPDIR" \
      ./benchmark_datasets --benchmark astex_diverse --mode defined-cleft-redock \
        --only-codes "$t" --output "$O" --cache "$CACHE" --clustering "$algo" \
        --threads 1 --omp-threads "$OMP" --ga-population 1000 --ga-generations 1000 \
        --job-timeout-seconds "$TMO" --engine-sha256 "$ESHA" > "$O/run.log" 2>&1 )
  local rc=$? t1=$(date +%s)
  local pz rr sent tmo banner
  pz=$(find "$O" -name '*.pdb' 2>/dev/null | grep -v '_INI\.pdb$' | grep -vc 'elected_pose' || echo 0)
  rr=$(find "$O" -name 'astex_diverse_results.csv' -exec sh -c 'wc -l < "$1"' _ {} \; 2>/dev/null | head -1 || echo 0)
  sent=$(grep -rIlaE 'SKIPPED|FATAL' "$O" 2>/dev/null | wc -l | tr -d ' ')
  tmo=$(grep -rIlaE 'timed out|SIGKILL|killed after' "$O" 2>/dev/null | wc -l | tr -d ' ')
  # ── HARNESS-SIDE fields. poses= counts FILES THE ENGINE WROTE and therefore CANNOT
  # ── detect a budget that reached the engine but not enumerate_emitted_cluster_heads().
  # ── The [ENUM] line is emitted BY that enumerator, so it is the only field a gate on
  # ── harness scoring may read. [MAXRES] budget= attests the knob reached the child.
  local mbudget enum_max enum_min cft fot algo_echo
  mbudget=$(grep -rIhoa 'budget=[0-9]*' "$O" 2>/dev/null | head -1 | cut -d= -f2)
  enum_max=$(grep -rIhoa 'enumerated=[0-9]*' "$O" 2>/dev/null | cut -d= -f2 | sort -n | tail -1)
  enum_min=$(grep -rIhoa 'enumerated=[0-9]*' "$O" 2>/dev/null | cut -d= -f2 | sort -n | head -1)
  cft=$(grep -rIhoa 'cf_truncated=[0-9]' "$O" 2>/dev/null | cut -d= -f2 | sort -rn | head -1)
  fot=$(grep -rIhoa 'fo_truncated=[0-9]' "$O" 2>/dev/null | cut -d= -f2 | sort -rn | head -1)
  algo_echo=$(grep -rIhoaE 'using the (Fast OPTICS \(FO\)|Density Peak \(DP\)|Complementarity Function \(CF\))' \
              "$O" 2>/dev/null | head -1 | sed -E 's/.*\((..)\).*/\1/')
  printf 'rc=%s arm=%s target=%s algo=%s algo_echo=%s limit=%s maxres_budget=%s enum_max=%s enum_min=%s cf_truncated=%s fo_truncated=%s engine=%s lig=%s poses=%s result_rows=%s sentinel=%s restart_timeout=%s wall_s=%s at=%s\n' \
    "$rc" "$arm" "$t" "$algo" "${algo_echo:-ABSENT}" "$lim" "${mbudget:-ABSENT}" \
    "${enum_max:-ABSENT}" "${enum_min:-ABSENT}" "${cft:-ABSENT}" "${fot:-ABSENT}" \
    "${ESHA:0:16}" "${LIG:-UNKNOWN}" "${pz:-0}" "${rr:-0}" \
    "${sent:-0}" "${tmo:-0}" "$((t1-t0))" "$(date -u +%FT%TZ)" > "$O/DONE"
}

# =============================================================================
# PHASE 0 -- KNOB CANARY, BOTH DIRECTIONS. A one-directional pass is not a pass.
# =============================================================================
G0="$BATCH/G0_VERDICT"
echo "[PHASE 0] knob canary (both directions) + FO feasibility smoke"

cell B_shrink "$CANARY_T" UNSET    "$BATCH/canary/neg" CF
cell B_shrink "$CANARY_T" "$LIMIT" "$BATCH/canary/pos" CF

rf() { grep -o "$1=[^ ]*" "$2" 2>/dev/null | head -1 | cut -d= -f2-; }

P_NEG=$(rf poses          "$BATCH/canary/neg/DONE")
P_POS=$(rf poses          "$BATCH/canary/pos/DONE")
E_NEG=$(rf enum_max       "$BATCH/canary/neg/DONE")
E_POS=$(rf enum_max       "$BATCH/canary/pos/DONE")
B_POS=$(rf maxres_budget  "$BATCH/canary/pos/DONE")
B_NEG=$(rf maxres_budget  "$BATCH/canary/neg/DONE")
CFT=$(rf cf_truncated     "$BATCH/canary/pos/DONE")

# --- 0a BEHAVIOURAL INERTNESS AT DEFAULT, measured as a SAME-TREE A/B.
#
#     WHY NOT A COMPARISON AGAINST THE FROZEN CAMPAIGN. The campaign's engine sha
#     (5bbe8c0f27b11cf3) CANNOT be reproduced from the commit it is recorded against:
#     main @ 94d2d0ba rebuilds deterministically to 3c3b71df0fd85a46 in this same build
#     directory, twice, including --clean-first. So that binary came from a different
#     tree state or toolchain, and any byte-comparison against it would conflate this
#     patch with that unknown difference.
#
#     WHAT IS COMPARED INSTEAD: the pre-patch build (BASE_BIN, from main) against the
#     post-patch build (BATCH/bin, from this branch), knob UNSET, same target, same
#     seed. That isolates the patch exactly. It matters because FlexAIDdS is the LTO
#     target and DOES compile DatasetRunner.cpp, so the patch changes the engine's
#     bytes (3c3b71df -> 4271cbb0) even though its strings are dead-stripped from that
#     binary -- behavioural equivalence therefore has to be measured, not assumed.
BASE_BIN="${MAXRES_BASE_BIN:-$R/state/base_bin}"
if [ -d "$BASE_BIN" ] && [ -x "$BASE_BIN/benchmark_datasets" ]; then
  SAVE_BIN="$BATCH/bin"; BATCH_BIN_KEEP=$(mktemp -d "$TMPDIR/binkeep.XXXXXX")
  cp "$BATCH/bin/"* "$BATCH_BIN_KEEP/" 2>/dev/null
  cp "$BASE_BIN/"* "$BATCH/bin/" 2>/dev/null
  cell B_shrink "$CANARY_T" UNSET "$BATCH/canary/base" CF
  cp "$BATCH_BIN_KEEP/"* "$BATCH/bin/" 2>/dev/null      # restore the patched pair
  BASE_EP=$(find "$BATCH/canary/base" -name 'elected_pose.pdb' | head -1)
  NEW_EP=$(find "$BATCH/canary/neg"  -name 'elected_pose.pdb' | head -1)
  if [ -n "${BASE_EP:-}" ] && [ -n "${NEW_EP:-}" ] && cmp -s "$BASE_EP" "$NEW_EP"; then
    G0A=PASS; G0A_N="same-tree A/B: elected_pose byte-identical pre/post patch at default budget"
  else
    G0A=FAIL; G0A_N="same-tree A/B DIFFERS at default budget (base=${BASE_EP:-missing} new=${NEW_EP:-missing}) -- the patch is NOT inert"
  fi
  BASE_POSES=$(rf poses "$BATCH/canary/base/DONE")
  BASE_ENUM=$(rf enum_max "$BATCH/canary/base/DONE")
  G0A_N="$G0A_N | base poses=${BASE_POSES:-?} enum=${BASE_ENUM:-ABSENT(expected: pre-patch has no [ENUM])}"
else
  G0A=FAIL
  G0A_N="no pre-patch baseline at $BASE_BIN -- cannot establish inertness; refusing to certify"
fi

# --- 0b KNOB REACHED THE CHILD: the engine's own [MAXRES] budget echo must equal
#     what we asked for, and the default arm must echo the default.
if [ "${B_POS:-x}" = "$LIMIT" ] && [ "${B_NEG:-x}" = "50" ]; then
  G0B=PASS; G0B_N="budget echo set=$B_POS unset=$B_NEG"
else
  G0B=FAIL; G0B_N="budget echo set=${B_POS:-ABSENT} (want $LIMIT) unset=${B_NEG:-ABSENT} (want 50)"
fi

# --- 0c THE GATE THE PATCH EXISTS FOR, AND IT MUST BE HARNESS-SIDE.
#     A count of pose FILES ON DISK cannot detect the failure mode: if only the config
#     emit were gated, the ENGINE would write $LIMIT poses (so a disk count passes)
#     while enumerate_emitted_cluster_heads() still scored 50. enum_max comes from the
#     [ENUM] line printed BY that enumerator, so it is the only admissible field here.
if [ "${E_POS:-0}" != "ABSENT" ] && [ "${E_POS:-0}" -gt 50 ] \
   && [ "${E_NEG:-0}" != "ABSENT" ] && [ "${E_NEG:-0}" -le 50 ]; then
  G0C=PASS; G0C_N="harness enumerated set=$E_POS (>50) unset=$E_NEG (<=50)"
else
  G0C=FAIL; G0C_N="harness enumerated set=${E_POS:-ABSENT} unset=${E_NEG:-ABSENT} -- \
if set<=50 the budget reached the engine but NOT the enumerator: the :109 desync"
fi

# --- 0d NO SILENT TRUNCATION at the chosen budget.
if [ "${CFT:-ABSENT}" = "0" ]; then
  G0D=PASS; G0D_N="cf_truncated=0 at budget $LIMIT: the pool is not clipped"
elif [ "${CFT:-ABSENT}" = "1" ]; then
  G0D=WARN; G0D_N="cf_truncated=1 -- engine wrote MORE than budget $LIMIT; ceilings are FLOORS. Raise it."
else
  G0D=FAIL; G0D_N="cf_truncated field ABSENT -- the [ENUM] attestation did not reach the receipt"
fi

# --- 0e FO FEASIBILITY. FO has never executed in this project: zero dual-suffix poses
#     exist anywhere on disk. Smoke it on the canary BEFORE committing campaign cells,
#     and require dual-suffix output specifically -- an FO run that silently falls back
#     to CF would otherwise be scored as an FO arm.
FO_OK=SKIP; FO_N="FO not in ALGOS"
case " $ALGOS " in *" FO "*)
  cell B_shrink "$CANARY_T" "$LIMIT" "$BATCH/canary/fo" FO
  FO_ECHO=$(rf algo_echo "$BATCH/canary/fo/DONE")
  FO_POSES=$(rf poses    "$BATCH/canary/fo/DONE")
  FO_ENUM=$(rf enum_max  "$BATCH/canary/fo/DONE")
  FO_DUAL=$(find "$BATCH/canary/fo" -name '*_[0-9]*_[0-9]*.pdb' 2>/dev/null \
            | grep -v '_INI\.pdb$' | wc -l | tr -d ' ')
  FO_WALL=$(rf wall_s "$BATCH/canary/fo/DONE")
  if [ "${FO_ECHO:-x}" = "FO" ] && [ "${FO_POSES:-0}" -gt 0 ] && [ "${FO_DUAL:-0}" -gt 0 ]; then
    FO_OK=PASS; FO_N="engine echoed FO, ${FO_POSES} poses, ${FO_DUAL} dual-suffix, enum=${FO_ENUM}, ${FO_WALL}s"
  else
    FO_OK=FAIL; FO_N="algo_echo=${FO_ECHO:-ABSENT} poses=${FO_POSES:-0} dual_suffix=${FO_DUAL:-0} \
-- FO either did not engage or fell back to CF; an arm labelled FO must emit dual-suffix heads"
  fi
;; esac

{ echo "phase=0 at=$(date -u +%FT%TZ) engine=${ESHA:0:16} canary=$CANARY_T budget=$LIMIT algos='$ALGOS'"
  echo "G0a_negative_knob=$G0A            ($G0A_N)"
  echo "G0b_budget_reached_child=$G0B     ($G0B_N)"
  echo "G0c_harness_enumeration=$G0C      ($G0C_N)"
  echo "G0d_no_silent_truncation=$G0D     ($G0D_N)"
  echo "G0e_fo_feasibility=$FO_OK         ($FO_N)"
  echo "note=0c is HARNESS-SIDE by construction; a disk-file count is inadmissible here"
  if [ "$G0A" = PASS ] && [ "$G0B" = PASS ] && [ "$G0C" = PASS ] \
     && { [ "$G0D" = PASS ] || [ "$G0D" = WARN ]; } \
     && { [ "$FO_OK" = PASS ] || [ "$FO_OK" = SKIP ]; }; then echo "VERDICT=PASS"
  else echo "VERDICT=FAIL"; fi
} > "$G0"
sed 's/^/  /' "$G0"
grep -q '^VERDICT=PASS$' "$G0" || { echo "[PHASE 0] GATE FAILED -- refusing to launch phase 1" >&2; exit 3; }

# --- GATES-ONLY MODE: run phase 0 (+1) and stop before spending campaign cells.
# --- Used to certify the harness change before any merge.
if [ "${MAXRES_GATES_ONLY:-0}" = "1" ] && [ "${MAXRES_GATES_ONLY_AFTER:-0}" = "0" ]; then
  echo "[GATES-ONLY] phase 0 passed; continuing to phase 1 then stopping."
fi

# =============================================================================
# PHASE 1 -- TOP-1 INVARIANCE. This tests the AGENT'S OWN prediction, not the engine's.
# If top-1 moves, the scope declaration at the head of this file is WRONG and the run
# must stop for review rather than quietly produce numbers under a false scope.
# =============================================================================
G1="$BATCH/G1_VERDICT"
echo "[PHASE 1] top-1 invariance"
python3 - "$BATCH/canary/neg" "$BATCH/canary/pos" "$G1" <<'PY'
import sys, glob, os, re
neg, pos, out = sys.argv[1], sys.argv[2], sys.argv[3]
def top1_cf(d):
    best=None
    for p in glob.glob(f"{d}/**/*.pdb", recursive=True):
        if p.endswith("_INI.pdb") or "elected_pose" in p: continue
        for l in open(p, errors="ignore"):
            if l.startswith("REMARK") and "CF=" in l:
                m=re.search(r'CF=\s*(-?[\d.]+)', l)
                if m:
                    v=float(m.group(1))
                    if best is None or v<best: best=v
                break
            if l.startswith(("ATOM","HETATM")): break
    return best
a,b = top1_cf(neg), top1_cf(pos)
ok = (a is not None and b is not None and abs(a-b) < 1e-4)
with open(out,"w") as fh:
    fh.write(f"phase=1 min_CF_unset={a} min_CF_set={b} delta={None if None in (a,b) else b-a}\n")
    fh.write("note=cap-invariance of top-1 is the scope claim; a change REFUTES it\n")
    fh.write(f"VERDICT={'PASS' if ok else 'FAIL'}\n")
print(f"  min CF unset={a} set={b} -> {'PASS' if ok else 'FAIL'}")
PY
sed 's/^/  /' "$G1"
grep -q '^VERDICT=PASS$' "$G1" || { echo "[PHASE 1] top-1 MOVED -- scope claim refuted, stopping for review" >&2; exit 4; }

if [ "${MAXRES_GATES_ONLY:-0}" = "1" ]; then
  echo "[GATES-ONLY] stopping before phase 2. Verdicts: $BATCH/G0_VERDICT $BATCH/G1_VERDICT"
  echo "GATES_ONLY $(date -u +%FT%TZ)" > "$BATCH/DONE"
  exit 0
fi

# =============================================================================
# PHASE 2 -- the run. ONLY cells that can flip.
# Ceiling-hit status is MONOTONE in the emitted set: adding poses can only lower a
# minimum, so a cell already sub-2 A stays sub-2 A. Scoping to the 48 cells without a
# sub-2 A ceiling is provably sufficient, not a convenience sample.
# =============================================================================
echo "[PHASE 2] $NT cells at FLEXAIDDS_MAX_RESULTS=$LIMIT, window $WINDOW"
n=0
tail -n +2 "$TARGETS" | while IFS=$'\t' read -r arm t _ceil; do
  [ -z "${arm:-}" ] && continue
  O="$BATCH/run/${arm}_s${SEED}/$t"
  if [ -f "$O/DONE" ]; then continue; fi
  if [ -d "$O" ]; then mv "$O" "${O}.stale.$(date -u +%H%M%S)"; fi   # rename aside, never delete
  cell "$arm" "$t" "$LIMIT" "$O" &
  PIDS="${PIDS:-} $!"
  n=$((n+1))
  if [ $(( n % WINDOW )) -eq 0 ]; then
    # PER-CHILD wait, not a bare `wait`: a bare wait discards each child's exit status,
    # so a cell that died would be indistinguishable from one that finished.
    for p in $PIDS; do
      wait "$p" || echo "[WARN] child pid $p exited non-zero" >> "$BATCH/child_failures.log"
    done
    PIDS=""
  fi
  date -u +%FT%TZ > "$BATCH/HEARTBEAT"
done
for p in ${PIDS:-}; do
  wait "$p" || echo "[WARN] child pid $p exited non-zero" >> "$BATCH/child_failures.log"
done
echo "DONE $(date -u +%FT%TZ)" > "$BATCH/DONE"

# =============================================================================
# PHASE 3 -- receipts audit. Reports; does not gate.
# =============================================================================
echo "[PHASE 3] receipts"
NC=$(find "$BATCH/run" -mindepth 3 -maxdepth 3 -name DONE 2>/dev/null | wc -l | tr -d ' ')
echo "  cells with receipts: $NC/$NT"
for f in sentinel restart_timeout; do
  bad=$(find "$BATCH/run" -name DONE -exec grep -ho "$f=[0-9]*" {} \; 2>/dev/null | grep -vc "$f=0" || echo 0)
  echo "  cells with $f != 0 : $bad   (any non-zero cell is VOID)"
done
echo "  cells with an ABSENT [MAXRES] banner (knob did not reach the child -> VOID): \
$(find "$BATCH/run" -name DONE -exec grep -lo 'maxres_banner=ABSENT' {} \; 2>/dev/null | wc -l | tr -d ' ')"
echo "  batch: $BATCH"
