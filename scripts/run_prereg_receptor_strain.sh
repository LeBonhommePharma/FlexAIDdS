#!/usr/bin/env bash
# =============================================================================
# run_prereg_receptor_strain.sh
#
# Staged driver for benchmarks/protocols/PREREG_receptor_strain_2026-09.md.
#
# IT IS NOT EXECUTABLE ON PURPOSE. Run it as:  bash scripts/run_prereg_receptor_strain.sh
# Do not chmod it in the repository.
#
# WHAT IT TESTS (see the pre-registration; do not re-derive it here):
#   whether pricing (1) receptor conformational strain, (2) the per-contact wall
#   ceiling for flexed side chains, and (3) the implicit-solvent reference,
#   removes the 16-target top-1 penalty AND the 38/77 interpenetration failures
#   that receptor flexibility currently costs.
#
# HOW IT IS BUILT
#   * PHASE -1 refuses to start if any docking process is alive. That is the most
#     important line in this file: a campaign launched on top of a running
#     benchmark corrupts both and is not recoverable.
#   * PHASES ARE GATES, NOT STEPS. Each writes VERDICT=PASS|FAIL to a file; the
#     next phase reads that file and refuses to run without PASS. rc=0 is not a
#     gate. A phase that "looks clean" is not a gate.
#   * EVERY GATE IS ASSERTED IN THE EMITTED dock_config.json BEFORE ANY CELL IS
#     SCORED. At least eight past features in this repository were structurally
#     absent because the harness never supplied the precondition while every log
#     looked clean. Phase 2 exists solely to make that impossible here.
#   * NO MACHINE-SPECIFIC ABSOLUTE PATHS. The repository root is derived from the
#     location of this script; everything else comes from FLEXAIDDS_* variables.
#   * The script contains ZERO delete verbs. It cannot remove the wrong path.
#   * It writes NOTHING inside the repository. All output goes under
#     $FLEXAIDDS_RESULTS.
#
# REQUIRED ENVIRONMENT
#   FLEXAIDDS_RESULTS         absolute results root (batch dir is created under it)
#   FLEXAIDDS_CACHE           benchmark cache PARENT (NOT the dataset subdir --
#                             passing the dataset dir is the same defect behind the
#                             incident where this project silently docked an additive)
#   FLEXAIDDS_ORACLE_SITE_DIR site directory for defined-cleft-redock
#   FLEXAIDDS_ARM_BIN         directory holding the PATCHED FlexAIDdS + benchmark_datasets
#   FLEXAIDDS_BASE_BIN        directory holding the PRE-PATCH pair (inertness control).
#                             Build it once from the pre-patch commit into a separate
#                             build dir and copy both binaries there. Without it the
#                             inertness gate CANNOT pass and the campaign is void.
#   FLEXAIDDS_TARGET_LIST     file with the 84 archived target codes, one per line.
#                             Taken VERBATIM from the archive; never re-derived.
# OPTIONAL
#   FLEXAIDDS_REPO            repo root override (default: parent of this script)
#   FLEXAIDDS_PREREG_CANARY   canary target code (default 1LPZ)
#   FLEXAIDDS_PREREG_OMP      omp threads per worker, identical in every arm (default 3)
#   FLEXAIDDS_PREREG_GENS     GA generations (default 3000, the archive value)
#   FLEXAIDDS_PREREG_POP      GA population (default 1000)
#   FLEXAIDDS_PREREG_SEED     seed base (default 12345, frozen by the pre-registration)
#   FLEXAIDDS_PREREG_RESTARTS restarts (default 3)
#   FLEXAIDDS_PREREG_AUTOFLEX flexible side chains (default 5, the archive value)
#   FLEXAIDDS_PREREG_WRITE_FLEXREC  1|0 (default 1) -- pure output, see the disk pilot
#   FLEXAIDDS_PREREG_TIMEOUT  per-complex timeout seconds (default 28800)
#   FLEXAIDDS_PREREG_PHASES   space-separated subset of: 0 1 2 3 4 5 6 (default: all)
#
# Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
# =============================================================================
set -u
set -o pipefail

# =============================================================================
# PHASE -1 -- MACHINE-IDLE GUARD.  THIS RUNS BEFORE ANYTHING ELSE.
#
# A docking campaign is already the single highest-cost thing this machine does.
# Starting a second one silently halves the throughput of both, perturbs every
# wall-clock number, and -- because the engine shares scratch by name -- can make
# one run poison the other's _dockin.sdf. There is no recovery and no way to tell
# afterwards which cells were affected. So: refuse, loudly, and re-check before
# every phase (a human can launch something while this script is running).
# =============================================================================
assert_machine_idle() {
    local where="${1:-startup}" pids ppid_self
    ppid_self="${PPID:-0}"
    pids=$(pgrep -f 'benchmark_datasets|FlexAIDdS|run_maxres\.sh' 2>/dev/null \
           | grep -vx "$$" | grep -vx "$ppid_self" || true)
    if [ -n "${pids:-}" ]; then
        echo "" >&2
        echo "################################################################" >&2
        echo "# REFUSING TO START ($where): docking processes are ALIVE." >&2
        echo "################################################################" >&2
        ps -o pid=,ppid=,etime=,command= -p $(echo "$pids" | tr '\n' ' ') 2>/dev/null >&2
        echo "" >&2
        echo "# This driver will not share the machine with a running benchmark." >&2
        echo "# Wait for it to finish. Do NOT kill it to make room for this." >&2
        echo "################################################################" >&2
        exit 3
    fi
}
assert_machine_idle "startup"

# =============================================================================
# PATHS.  Repo-relative or FLEXAIDDS_*. Never a machine-specific literal.
# =============================================================================
SELF_DIR=$(cd "$(dirname "$0")" && pwd) || exit 2
REPO="${FLEXAIDDS_REPO:-$(cd "$SELF_DIR/.." && pwd)}"
[ -d "$REPO/LIB" ] || { echo "FATAL: $REPO does not look like the FlexAIDdS repo" >&2; exit 2; }

req_abs() {  # $1 varname  $2 value
    case "${2:-}" in
        "")  echo "FATAL: $1 is required and unset" >&2; exit 2 ;;
        /*)  ;;
        *)   echo "FATAL: $1 must be an ABSOLUTE path (got '$2')" >&2; exit 2 ;;
    esac
}
R="${FLEXAIDDS_RESULTS:-}";            req_abs FLEXAIDDS_RESULTS "$R"
CACHE="${FLEXAIDDS_CACHE:-}";          req_abs FLEXAIDDS_CACHE "$CACHE"
SITES="${FLEXAIDDS_ORACLE_SITE_DIR:-}";req_abs FLEXAIDDS_ORACLE_SITE_DIR "$SITES"
ARM_BIN="${FLEXAIDDS_ARM_BIN:-}";      req_abs FLEXAIDDS_ARM_BIN "$ARM_BIN"
BASE_BIN="${FLEXAIDDS_BASE_BIN:-}";    req_abs FLEXAIDDS_BASE_BIN "$BASE_BIN"
TLIST="${FLEXAIDDS_TARGET_LIST:-}";    req_abs FLEXAIDDS_TARGET_LIST "$TLIST"
[ -f "$TLIST" ] || { echo "FATAL: target list not found: $TLIST" >&2; exit 2; }

CANARY="${FLEXAIDDS_PREREG_CANARY:-1LPZ}"
OMP="${FLEXAIDDS_PREREG_OMP:-3}"
GENS="${FLEXAIDDS_PREREG_GENS:-3000}"
POP="${FLEXAIDDS_PREREG_POP:-1000}"
SEED="${FLEXAIDDS_PREREG_SEED:-12345}"
RESTARTS="${FLEXAIDDS_PREREG_RESTARTS:-3}"
AUTOFLEX="${FLEXAIDDS_PREREG_AUTOFLEX:-5}"
WRITE_FLEXREC="${FLEXAIDDS_PREREG_WRITE_FLEXREC:-1}"
TMO="${FLEXAIDDS_PREREG_TIMEOUT:-28800}"
PHASES="${FLEXAIDDS_PREREG_PHASES:-0 1 2 3 4 5 6}"
DATASET=astex_diverse
ARMS="A_rigid B_shrink T1_strain T2_walcap T3_solvref T4_all"
# Archive constants, quoted from the pre-registration. NEVER recomputed here.
ARCH_A=48; ARCH_B=32; ARCH_N=84; ARCH_TOL=3
case "$WRITE_FLEXREC" in 1) WFR_JSON=true ;; *) WFR_JSON=false ;; esac

BATCH="$R/prereg_receptor_strain_$(date -u +%Y%m%d_%H%M%S)"
mkdir -p "$BATCH/bin" "$BATCH/run" "$BATCH/tmp" "$BATCH/verdict" || exit 2
# Absolute per-batch scratch. NEVER inherit TMPDIR: a relative inherited TMPDIR
# has already cost this project a whole arm.
export TMPDIR="$BATCH/tmp"
: > "$TMPDIR/.writable" || { echo "FATAL: TMPDIR not writable: $TMPDIR" >&2; exit 2; }
LOG="$BATCH/driver.log"
say() { echo "$*" | tee -a "$LOG"; }
say "[prereg] batch      = $BATCH"
say "[prereg] repo       = $REPO"
say "[prereg] targets    = $TLIST ($(grep -c '[A-Za-z0-9]' "$TLIST") codes)"
say "[prereg] arms       = $ARMS"
say "[prereg] seed=$SEED restarts=$RESTARTS gens=$GENS pop=$POP omp=$OMP autoflex=$AUTOFLEX"

# ---- verdict plumbing: a phase may not start unless its prerequisites read PASS
verdict_write() {  # $1 name  $2 PASS|FAIL  $3.. notes
    local n=$1 v=$2; shift 2
    { echo "phase=$n at=$(date -u +%FT%TZ)"
      for l in "$@"; do echo "  $l"; done
      echo "VERDICT=$v"
    } > "$BATCH/verdict/$n"
    say "[verdict] $n = $v"
}
verdict_require() {  # $1.. verdict names that must read PASS
    local n
    for n in "$@"; do
        if [ ! -f "$BATCH/verdict/$n" ]; then
            say "FATAL: gate $n has not run; refusing to continue."; exit 4
        fi
        if ! grep -q '^VERDICT=PASS$' "$BATCH/verdict/$n"; then
            say "FATAL: gate $n did not PASS; refusing to continue."
            sed 's/^/    /' "$BATCH/verdict/$n" | tee -a "$LOG"; exit 4
        fi
    done
}
phase_enabled() { case " $PHASES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

# =============================================================================
# ARM DEFINITIONS -- the single place any arm's environment is written down.
# Everything not switched by the arm is IDENTICAL in every arm, on purpose.
# =============================================================================
arm_env() {   # $1 = arm name -> prints "K=V K=V ..." for `env`
    case "$1" in
        A_rigid) : ;;   # rigid: FLEXAIDDS_AUTOFLEX_MAX deliberately UNSET (=0)
        *) printf 'FLEXAIDDS_AUTOFLEX_MAX=%s FLEXAIDDS_AUTOFLEX_METAL_SHRINK=1 ' "$AUTOFLEX" ;;
    esac
    case "$1" in T1_strain|T4_all)  printf 'FLEXAIDDS_RECEPTOR_STRAIN=1 ' ;; esac
    case "$1" in T2_walcap|T4_all)  printf 'FLEXAIDDS_WAL_CAP_MODE=flex ' ;; esac
    case "$1" in T3_solvref|T4_all) printf 'FLEXAIDDS_SOLVATION_REF=crystal ' ;; esac
    # Instruments: ON IN EVERY ARM INCLUDING THE RIGID BASELINE. An instrument
    # that is only on in the treated arm is not an instrument, it is a confound.
    printf 'FLEXAIDDS_CONTACT_PROFILE=1 FLEXAIDDS_WRITE_FLEXED_RECEPTOR=%s ' "$WRITE_FLEXREC"
    # Measurement frame is CRYSTAL in every arm. The flexed frame is a different
    # estimand and is NOT part of this campaign.
    printf 'FLEXAIDDS_PB_RECEPTOR=crystal '
    # Archive invariants.
    printf 'FLEXAIDDS_SCORED_ONLY=1 FLEXAIDDS_NO_SEC=1 '
    printf 'FLEXAIDDS_SEED_BASE=%s FLEXAIDDS_RESTARTS=%s FLEXAIDDS_ORACLE_SITE_DIR=%s ' \
           "$SEED" "$RESTARTS" "$SITES"
}

# Expected dock_config.json contents per arm. Bools and strings are asserted by
# VALUE; the two floating-point keys are asserted by PRESENCE only, because their
# ostream formatting is not something this script can verify without a compiler.
arm_expect() {   # $1 = arm -> one grep -F pattern per line
    local rs=false wcm=legacy sr=dynamic afx=0
    case "$1" in T1_strain|T4_all)  rs=true ;; esac
    case "$1" in T2_walcap|T4_all)  wcm=flex ;; esac
    case "$1" in T3_solvref|T4_all) sr=crystal ;; esac
    case "$1" in A_rigid) afx=0 ;; *) afx=$AUTOFLEX ;; esac
    printf '"receptor_strain": %s\n' "$rs"
    printf '"receptor_strain_temperature_K":\n'
    printf '"wal_cap_mode": "%s"\n' "$wcm"
    printf '"wal_cap_flex":\n'
    printf '"solvation_ref": "%s"\n' "$sr"
    printf '"autoflex_max": %s\n' "$afx"
    printf '"contact_profile": true\n'
    printf '"write_flexed_receptor": %s\n' "$WFR_JSON"
    printf '"pb_receptor": "crystal"\n'
}

# Assert the arm's gates in the CELL'S OWN emitted dock_config.json.
# Returns 0 = all gates present, 1 = at least one missing (cell is VOID).
assert_cell_gates() {   # $1 arm  $2 cell outdir  -> writes $2/GATES
    local arm=$1 O=$2 missing=0 pat ncfg blob
    blob="$O/.dock_config_concat"
    : > "$blob"
    ncfg=0
    while IFS= read -r p; do
        [ -z "$p" ] && continue
        cat "$p" >> "$blob"; ncfg=$((ncfg+1))
    done <<EOF_CFG
$(find "$O" -name dock_config.json 2>/dev/null)
EOF_CFG
    if [ "$ncfg" -eq 0 ]; then
        echo "MISSING dock_config.json under $O" > "$O/GATES"; return 1
    fi
    : > "$O/GATES"
    while IFS= read -r pat; do
        [ -z "$pat" ] && continue
        if grep -qF -- "$pat" "$blob"; then
            echo "ok      $pat" >> "$O/GATES"
        else
            echo "MISSING $pat" >> "$O/GATES"; missing=$((missing+1))
        fi
    done <<EOF_PAT
$(arm_expect "$arm")
EOF_PAT
    echo "n_config_files=$ncfg missing=$missing" >> "$O/GATES"
    [ "$missing" -eq 0 ]
}

# =============================================================================
# CELL RUNNER
# =============================================================================
ESHA=""
cell() {   # $1 arm  $2 target  $3 outdir  [$4 bindir override]
    assert_machine_idle "cell $1/$2"
    local arm=$1 t=$2 O=$3 B=${4:-$BATCH/bin}
    mkdir -p "$O"
    local t0 t1 rc
    t0=$(date +%s)
    ( cd "$B" && env $(arm_env "$arm") TMPDIR="$TMPDIR" \
        ./benchmark_datasets --benchmark "$DATASET" --mode defined-cleft-redock \
          --only-codes "$t" --output "$O" --cache "$CACHE" --clustering CF \
          --threads 1 --omp-threads "$OMP" \
          --ga-population "$POP" --ga-generations "$GENS" \
          --job-timeout-seconds "$TMO" ${ESHA:+--engine-sha256 "$ESHA"} \
          > "$O/run.log" 2> "$O/run.err" )
    rc=$?; t1=$(date +%s)
    local gates=VOID
    if assert_cell_gates "$arm" "$O"; then gates=OK; fi
    local poses csv rmsd sr pb vov keys strain wal solv offrot
    poses=$(find "$O" -name '*.pdb' 2>/dev/null | grep -v '_INI\.pdb$' \
            | grep -v '/flexed_receptor/' | grep -vc 'elected_pose' || echo 0)
    csv=$(find "$O" -name result.csv 2>/dev/null | head -1)
    rmsd=ABSENT; sr=ABSENT; pb=ABSENT; vov=ABSENT; keys=ABSENT
    if [ -n "${csv:-}" ]; then
        rmsd=$(python3 - "$csv" rmsd_to_crystal <<'PY'
import csv,sys
try:
    r=next(csv.DictReader(open(sys.argv[1])))
    print(r.get(sys.argv[2],"ABSENT") or "ABSENT")
except Exception: print("ABSENT")
PY
)
        sr=$(python3 - "$csv" success_rmsd  <<'PY'
import csv,sys
try:
    r=next(csv.DictReader(open(sys.argv[1])))
    print(r.get(sys.argv[2],"ABSENT") or "ABSENT")
except Exception: print("ABSENT")
PY
)
        pb=$(python3 - "$csv" pb_pass <<'PY'
import csv,sys
try:
    r=next(csv.DictReader(open(sys.argv[1])))
    print(r.get(sys.argv[2],"ABSENT") or "ABSENT")
except Exception: print("ABSENT")
PY
)
        vov=$(python3 - "$csv" pb_volume_overlap <<'PY'
import csv,sys
try:
    r=next(csv.DictReader(open(sys.argv[1])))
    print(r.get(sys.argv[2],"ABSENT") or "ABSENT")
except Exception: print("ABSENT")
PY
)
        keys=$(python3 - "$csv" pb_failed_keys <<'PY'
import csv,sys
try:
    r=next(csv.DictReader(open(sys.argv[1])))
    print((r.get(sys.argv[2],"") or "NONE").replace(" ","_"))
except Exception: print("ABSENT")
PY
)
    fi
    # LIVENESS witnesses -- a gate that is set but never fires is a null
    # measurement being reported as a result. Count the events, per cell.
    strain=$(grep -rIha '^REMARK CF\.strain' "$O" 2>/dev/null \
             | awk '{ if ($3+0 != 0) n++ } END { print n+0 }')
    wal=$(grep -rIhac 'WAL_CAP' "$O" 2>/dev/null | awk '{s+=$1} END {print s+0}')
    solv=$(grep -rIhac 'SOLV_REF' "$O" 2>/dev/null | awk '{s+=$1} END {print s+0}')
    offrot=$(grep -rIha 'n_residues_off_input_rotamer' "$O" 2>/dev/null \
             | awk '{ if ($3+0 > 0) n++ } END { print n+0 }')
    printf 'rc=%s arm=%s target=%s gates=%s poses=%s rmsd=%s success_rmsd=%s pb_pass=%s vol_overlap=%s pb_failed=%s strain_nonzero=%s wal_events=%s solv_events=%s offrot_files=%s omp=%s seed=%s engine=%s wall_s=%s at=%s\n' \
      "$rc" "$arm" "$t" "$gates" "${poses:-0}" "${rmsd:-ABSENT}" "${sr:-ABSENT}" \
      "${pb:-ABSENT}" "${vov:-ABSENT}" "${keys:-ABSENT}" "${strain:-0}" "${wal:-0}" \
      "${solv:-0}" "${offrot:-0}" "$OMP" "$SEED" "${ESHA:0:16}" "$((t1-t0))" \
      "$(date -u +%FT%TZ)" > "$O/DONE"
}
rf() { grep -o "$1=[^ ]*" "$2" 2>/dev/null | head -1 | cut -d= -f2-; }

# =============================================================================
# PHASE 0 -- PROVENANCE PREFLIGHT
# =============================================================================
if phase_enabled 0; then
  assert_machine_idle "phase 0"
  say "[PHASE 0] provenance preflight"
  P0=PASS; P0N=""
  for f in FlexAIDdS benchmark_datasets; do
      if [ -x "$ARM_BIN/$f" ]; then cp "$ARM_BIN/$f" "$BATCH/bin/"
      else P0=FAIL; P0N="$P0N missing:$ARM_BIN/$f"; fi
  done
  for f in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def rotobs.lst; do
      src=$(find "$ARM_BIN" -maxdepth 1 -name "$f" 2>/dev/null | head -1)
      [ -n "${src:-}" ] && cp "$src" "$BATCH/bin/"
  done
  [ -x "$BATCH/bin/FlexAIDdS" ] && ESHA=$(shasum -a 256 "$BATCH/bin/FlexAIDdS" | awk '{print $1}')
  # Determinism precondition: the cleft-grid probe-merge fix. Before it, SURFNET
  # merge order permuted the cleftgrid indices that become GA gene 0, so two arms
  # differ for a reason that has nothing to do with any gate under test.
  if git -C "$REPO" merge-base --is-ancestor 8dc88b4e HEAD 2>/dev/null; then
      P0N="$P0N cleft_grid_fix:present"
  else
      P0=FAIL
      P0N="$P0N cleft_grid_fix:ABSENT(8dc88b4e not an ancestor of HEAD -- arms are not comparable)"
  fi
  # All six gates must be known to the flag registry, or FLAGS_DUMP is not a record.
  DUMP=$( (cd "$BATCH/bin" && FLEXAIDDS_FLAGS_DUMP=1 ./benchmark_datasets --benchmark "$DATASET" \
            --list-codes 2>&1) | head -400 )
  # If the dump did not run at all on this code path, that is a WARN: it is a
  # convenience record, not the gate. dock_config.json (phase 2) is the gate, and
  # it is asserted per cell. But if the dump DID run and a variable is absent from
  # it, the registry is genuinely incomplete and that is a FAIL.
  if echo "$DUMP" | grep -q 'FLEXAIDDS_'; then
      for v in FLEXAIDDS_RECEPTOR_STRAIN FLEXAIDDS_WAL_CAP_MODE FLEXAIDDS_SOLVATION_REF \
               FLEXAIDDS_CONTACT_PROFILE FLEXAIDDS_WRITE_FLEXED_RECEPTOR FLEXAIDDS_PB_RECEPTOR; do
          echo "$DUMP" | grep -q "$v" || { P0=FAIL; P0N="$P0N flags_dump_missing:$v"; }
      done
  else
      P0N="$P0N flags_dump:UNAVAILABLE_ON_--list-codes(WARN; dock_config.json in phase 2 is the real gate)"
  fi
  NT=$(grep -c '[A-Za-z0-9]' "$TLIST")
  [ "$NT" -eq "$ARCH_N" ] || { P0=FAIL; P0N="$P0N target_count=$NT want=$ARCH_N"; }
  verdict_write G0_preflight "$P0" "engine=$ESHA" "targets=$NT" "notes:$P0N"
fi

# =============================================================================
# PHASE 1 -- INERTNESS.  Gates unset => byte-identical to the PRE-PATCH binary.
# This is the non-negotiable direction. If it fails, the campaign is VOID (F4).
# =============================================================================
if phase_enabled 1; then
  verdict_require G0_preflight
  assert_machine_idle "phase 1"
  say "[PHASE 1] inertness: all gates UNSET, patched vs pre-patch binary"
  mkdir -p "$BATCH/bin_base"
  for f in FlexAIDdS benchmark_datasets; do
      [ -x "$BASE_BIN/$f" ] && cp "$BASE_BIN/$f" "$BATCH/bin_base/"
  done
  for f in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def rotobs.lst; do
      src=$(find "$BASE_BIN" "$ARM_BIN" -maxdepth 1 -name "$f" 2>/dev/null | head -1)
      [ -n "${src:-}" ] && cp "$src" "$BATCH/bin_base/"
  done
  if [ ! -x "$BATCH/bin_base/benchmark_datasets" ]; then
      verdict_write G1_inertness FAIL \
        "no pre-patch build at FLEXAIDDS_BASE_BIN=$BASE_BIN" \
        "inertness CANNOT be established; per the pre-registration (F4) the campaign is VOID"
  else
      # A bare arm with EVERY gate removed from the environment, on both binaries.
      inert_cell() {  # $1 outdir  $2 bindir
          assert_machine_idle "inertness"
          mkdir -p "$1"
          ( cd "$2" && env -u FLEXAIDDS_RECEPTOR_STRAIN -u FLEXAIDDS_RECEPTOR_STRAIN_T \
                 -u FLEXAIDDS_WAL_CAP_MODE -u FLEXAIDDS_WAL_CAP_FLEX \
                 -u FLEXAIDDS_SOLVATION_REF -u FLEXAIDDS_CONTACT_PROFILE \
                 -u FLEXAIDDS_WRITE_FLEXED_RECEPTOR -u FLEXAIDDS_PB_RECEPTOR \
                 -u FLEXAIDDS_AUTOFLEX_MAX -u FLEXAIDDS_AUTOFLEX_METAL_SHRINK \
                 FLEXAIDDS_SCORED_ONLY=1 FLEXAIDDS_NO_SEC=1 \
                 FLEXAIDDS_SEED_BASE="$SEED" FLEXAIDDS_RESTARTS="$RESTARTS" \
                 FLEXAIDDS_ORACLE_SITE_DIR="$SITES" \
                 FLEXAIDDS_PARALLEL_RESTARTS=0 TMPDIR="$TMPDIR" \
              ./benchmark_datasets --benchmark "$DATASET" --mode defined-cleft-redock \
                --only-codes "$CANARY" --output "$1" --cache "$CACHE" --clustering CF \
                --threads 1 --omp-threads 1 --ga-population "$POP" \
                --ga-generations "$GENS" --job-timeout-seconds "$TMO" \
                > "$1/run.log" 2> "$1/run.err" )
      }
      inert_cell "$BATCH/inert/base" "$BATCH/bin_base"
      inert_cell "$BATCH/inert/post" "$BATCH/bin"
      BEP=$(find "$BATCH/inert/base" -name 'elected_pose.pdb' | head -1)
      NEP=$(find "$BATCH/inert/post" -name 'elected_pose.pdb' | head -1)
      # -r with a DIRECTORY argument, never `grep $(find ...)`: an empty find would
      # otherwise leave grep reading stdin and the driver would hang forever.
      BMD=$(grep -rIh --include='*.pdb' '^REMARK CF\.' "$BATCH/inert/base" 2>/dev/null | shasum -a 256 | awk '{print $1}')
      NMD=$(grep -rIh --include='*.pdb' '^REMARK CF\.' "$BATCH/inert/post" 2>/dev/null | shasum -a 256 | awk '{print $1}')
      SIDE=$(find "$BATCH/inert/post" \( -name '*.cprof.csv' -o -name '*_receptor.pdb' \) 2>/dev/null | wc -l | tr -d ' ')
      NOISE=$(grep -rIhac 'WAL_CAP\|SOLV_REF\|CPROF\|CF\.strain' "$BATCH/inert/post" 2>/dev/null | awk '{s+=$1} END {print s+0}')
      G1=FAIL; G1N="unset"
      if [ -n "${BEP:-}" ] && [ -n "${NEP:-}" ] && cmp -s "$BEP" "$NEP" \
         && [ -n "${BMD:-}" ] && [ "$BMD" = "$NMD" ] \
         && [ "${SIDE:-1}" -eq 0 ] && [ "${NOISE:-1}" -eq 0 ]; then
          G1=PASS; G1N="elected_pose byte-identical; CF-remark digest equal ($BMD); 0 sidecars; 0 gate events"
      else
          G1N="pose_cmp=$( [ -n "${BEP:-}" ] && [ -n "${NEP:-}" ] && cmp -s "$BEP" "$NEP" && echo same || echo DIFFER ) cf_digest_base=${BMD:-ABSENT} cf_digest_post=${NMD:-ABSENT} sidecars=${SIDE:-?} gate_events=${NOISE:-?}"
      fi
      verdict_write G1_inertness "$G1" "canary=$CANARY" "$G1N" \
        "F4: a FAIL here voids the whole campaign; no arm may be reported"
  fi
fi

# =============================================================================
# PHASE 2 -- WIRING.  EVERY GATE MUST APPEAR IN THE EMITTED dock_config.json,
# FOR EVERY ARM, BEFORE ANY CELL IS SCORED. One canary cell per arm.
# =============================================================================
if phase_enabled 2; then
  verdict_require G0_preflight G1_inertness
  assert_machine_idle "phase 2"
  say "[PHASE 2] wiring: dock_config.json must carry every gate, per arm"
  G2=PASS; G2N=""
  for a in $ARMS; do
      O="$BATCH/wiring/$a"
      cell "$a" "$CANARY" "$O"
      if [ "$(rf gates "$O/DONE")" = OK ]; then
          G2N="$G2N $a:ok"
      else
          G2=FAIL
          G2N="$G2N $a:MISSING($(grep -c '^MISSING' "$O/GATES" 2>/dev/null || echo ?))"
          say "  --- $a missing gates:"; grep '^MISSING' "$O/GATES" 2>/dev/null | sed 's/^/      /' | tee -a "$LOG"
      fi
  done
  verdict_write G2_wiring "$G2" "canary=$CANARY" "arms:$G2N" \
    "a gate absent from dock_config.json means the cell cannot be attributed to an arm"
fi

# =============================================================================
# PHASE 3 -- LIVENESS.  A gate that is SET but never FIRES is a null measurement
# being reported as a result. Reuse the phase-2 canary cells; no new docking.
# =============================================================================
if phase_enabled 3; then
  verdict_require G2_wiring
  assert_machine_idle "phase 3"
  say "[PHASE 3] liveness: each treatment must be observed to fire"
  G3=PASS; G3N=""
  chk() {  # $1 arm  $2 field  $3 human label
      local v; v=$(rf "$2" "$BATCH/wiring/$1/DONE")
      if [ "${v:-0}" -gt 0 ] 2>/dev/null; then G3N="$G3N $1/$3=$v"
      else G3=FAIL; G3N="$G3N $1/$3=ZERO(NOT_LIVE)"; fi
  }
  chk T1_strain  strain_nonzero  CF.strain
  chk T4_all     strain_nonzero  CF.strain
  chk T2_walcap  wal_events      WAL_CAP
  chk T4_all     wal_events      WAL_CAP
  chk T3_solvref solv_events     SOLV_REF
  chk T4_all     solv_events     SOLV_REF
  if [ "$WRITE_FLEXREC" = 1 ]; then
      for a in B_shrink T1_strain T2_walcap T3_solvref T4_all; do
          chk "$a" offrot_files n_res_off_rotamer
      done
      # The rigid arm is the NULL control: it must be exactly zero.
      v=$(rf offrot_files "$BATCH/wiring/A_rigid/DONE")
      if [ "${v:-x}" = 0 ]; then G3N="$G3N A_rigid/n_res_off_rotamer=0(null_control_ok)"
      else G3=FAIL; G3N="$G3N A_rigid/n_res_off_rotamer=$v(RIGID_ARM_MOVED_A_SIDECHAIN)"; fi
  fi
  verdict_write G3_liveness "$G3" "$G3N" \
    "F5: a FAIL here means the treatment never reached the model; escalate, do not patch around it"
fi

# =============================================================================
# PHASE 4 -- CONTROL REPLICATION.  A_rigid and B_shrink over the 84 targets.
# If these do not reproduce the archive within +/-3, no comparison against
# 48/84 or 32/84 is licensed and no treated cell may be spent (F6).
# =============================================================================
if phase_enabled 4; then
  verdict_require G2_wiring G3_liveness
  assert_machine_idle "phase 4"
  say "[PHASE 4] control replication (A_rigid, B_shrink) over $ARCH_N targets"
  for a in A_rigid B_shrink; do
      while IFS= read -r t; do
          t=$(echo "$t" | tr -d '[:space:]'); [ -z "$t" ] && continue
          [ -f "$BATCH/run/$a/$t/DONE" ] && continue
          cell "$a" "$t" "$BATCH/run/$a/$t"
      done < "$TLIST"
  done
  count_s1() {  # $1 arm -> "<successes> <complete> <void>"
      local ok=0 n=0 void=0 f sr g
      for f in "$BATCH/run/$1"/*/DONE; do
          [ -f "$f" ] || continue
          n=$((n+1)); g=$(rf gates "$f"); sr=$(rf success_rmsd "$f")
          [ "$g" != OK ] && void=$((void+1))
          [ "$sr" = 1 ] && [ "$g" = OK ] && ok=$((ok+1))
      done
      echo "$ok $n $void"
  }
  set -- $(count_s1 A_rigid); A_OK=$1; A_N=$2; A_V=$3
  set -- $(count_s1 B_shrink); B_OK=$1; B_N=$2; B_V=$3
  G4=PASS; G4N="A_rigid=$A_OK/$A_N(void=$A_V) B_shrink=$B_OK/$B_N(void=$B_V) archive=$ARCH_A/$ARCH_B tol=+-$ARCH_TOL"
  absdiff() { local d=$(( $1 - $2 )); [ $d -lt 0 ] && d=$(( -d )); echo $d; }
  [ "$A_N" -eq "$ARCH_N" ] || { G4=FAIL; G4N="$G4N INCOMPLETE:A_rigid"; }
  [ "$B_N" -eq "$ARCH_N" ] || { G4=FAIL; G4N="$G4N INCOMPLETE:B_shrink"; }
  [ "$A_V" -eq 0 ] && [ "$B_V" -eq 0 ] || { G4=FAIL; G4N="$G4N VOID_CELLS_PRESENT"; }
  [ "$(absdiff "$A_OK" "$ARCH_A")" -le "$ARCH_TOL" ] || { G4=FAIL; G4N="$G4N DRIFT:A_rigid"; }
  [ "$(absdiff "$B_OK" "$ARCH_B")" -le "$ARCH_TOL" ] || { G4=FAIL; G4N="$G4N DRIFT:B_shrink"; }
  verdict_write G4_control "$G4" "$G4N" \
    "NOTE these counts use E2 (the pose the engine ELECTS). The confirmatory endpoint" \
    "is E1 (argmin CF), computed in phase 6. E2 is used here only as the cheap" \
    "drift detector; a PASS here does not license any E1 claim." \
    "F6: a FAIL means the environment is not the archive's; do not spend treated cells"
fi

# =============================================================================
# PHASE 5 -- TREATED ARMS.  Only after the controls replicated.
# =============================================================================
if phase_enabled 5; then
  verdict_require G4_control
  assert_machine_idle "phase 5"
  say "[PHASE 5] treated arms T1..T4 over $ARCH_N targets"
  G5=PASS; G5N=""
  for a in T1_strain T2_walcap T3_solvref T4_all; do
      while IFS= read -r t; do
          t=$(echo "$t" | tr -d '[:space:]'); [ -z "$t" ] && continue
          [ -f "$BATCH/run/$a/$t/DONE" ] && continue
          cell "$a" "$t" "$BATCH/run/$a/$t"
      done < "$TLIST"
      n=$(ls -1 "$BATCH/run/$a"/*/DONE 2>/dev/null | wc -l | tr -d ' ')
      v=$(grep -l 'gates=VOID' "$BATCH/run/$a"/*/DONE 2>/dev/null | wc -l | tr -d ' ')
      G5N="$G5N $a=$n/$ARCH_N(void=$v)"
      { [ "$n" -eq "$ARCH_N" ] && [ "$v" -eq 0 ]; } || G5=FAIL
  done
  verdict_write G5_treated "$G5" "$G5N" \
    "F7: an arm with fewer than $ARCH_N complete, gate-verified cells is INCOMPLETE," \
    "not partially scored. Do not score it."
fi

# =============================================================================
# PHASE 6 -- SCORING.  Endpoints, tiers and tests exactly as pre-registered.
# =============================================================================
if phase_enabled 6; then
  verdict_require G1_inertness G2_wiring G3_liveness G4_control G5_treated
  assert_machine_idle "phase 6"
  say "[PHASE 6] scoring"
  # E1 (argmin CF) needs a per-pose RMSD the harness does not put in result.csv.
  # scripts/backfill_inline_rmsd.py already computes exactly that (cf_top1_rmsd),
  # fail-closed, into a SIDE-BY-SIDE copy. If it is unavailable or leaves any cell
  # unresolved, E1 is DEFERRED and NO confirmatory verdict is issued -- E2 is a
  # different endpoint and may not be substituted (pre-registration section 7).
  BF="$REPO/scripts/backfill_inline_rmsd.py"
  E1_STATUS=DEFERRED
  if [ -f "$BF" ]; then
      RD=""
      for a in $ARMS; do for d in "$BATCH/run/$a"/*; do [ -d "$d" ] && RD="$RD --run-dir $d"; done; done
      if python3 "$BF" $RD --cache "$CACHE" --dataset "$DATASET" \
             > "$BATCH/e1_backfill.log" 2>&1; then
          E1_STATUS=OK
      else
          E1_STATUS="DEFERRED(backfill_failed; see $BATCH/e1_backfill.log)"
      fi
  else
      E1_STATUS="DEFERRED(no $BF)"
  fi
  python3 - "$BATCH" "$ARCH_N" "$E1_STATUS" <<'PYEOF' | tee -a "$LOG" > "$BATCH/RESULTS.txt"
import csv, glob, os, sys
from math import comb
BATCH, N, E1_STATUS = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ARMS = ["A_rigid","B_shrink","T1_strain","T2_walcap","T3_solvref","T4_all"]
WATER = {"minimum_distance_to_waters", "volume_overlap_with_waters"}

def p2(b, c):
    n = b + c
    if n == 0: return 1.0
    k = min(b, c)
    return min(1.0, 2 * sum(comb(n, i) for i in range(k + 1)) / 2 ** n)

def read_done(p):
    d = {}
    for tok in open(p).read().split():
        if "=" in tok:
            k, v = tok.split("=", 1); d[k] = v
    return d

def e1_rmsd(cell):
    for f in glob.glob(os.path.join(cell, "**", "*_backfilled.csv"), recursive=True):
        try:
            r = next(csv.DictReader(open(f)))
        except Exception:
            continue
        v = r.get("cf_top1_rmsd")
        if v not in (None, "", "ABSENT"):
            try: return float(v)
            except ValueError: return None
    return None

data = {}
for a in ARMS:
    cells = {}
    for done in sorted(glob.glob(os.path.join(BATCH, "run", a, "*", "DONE"))):
        cell = os.path.dirname(done); t = os.path.basename(cell)
        d = read_done(done)
        try: rm = float(d.get("rmsd", "nan"))
        except ValueError: rm = float("nan")
        keys = d.get("pb_failed", "NONE")
        keys = set() if keys in ("NONE", "ABSENT", "") else set(keys.split(","))
        cells[t] = {
            "gates": d.get("gates"),
            # GUARD rmsd >= 0 ALWAYS: -1.0 is the shared-_dockin.sdf race sentinel
            # and passes a naive `rmsd < 2` filter.
            "E2": (rm >= 0.0 and rm <= 2.0),
            "E1_rmsd": e1_rmsd(cell),
            "pb_pass": d.get("pb_pass") == "1",
            "failed": keys,
            "vov_fail": "volume_overlap_with_protein" in keys,
        }
    data[a] = cells

def tier(c, which):
    ok = c["E2"] if which == "E2" else (c["E1_rmsd"] is not None and 0.0 <= c["E1_rmsd"] <= 2.0)
    return {
        "V0": ok,
        "V1": ok and c["pb_pass"],
        "V2": ok and not (c["failed"] - WATER),
    }

print("=" * 78)
print("PRE-REGISTERED RESULTS -- benchmarks/protocols/PREREG_receptor_strain_2026-09.md")
print("=" * 78)
print("E1 (argmin CF, PRIMARY) status: %s" % E1_STATUS)
print("")
print("%-11s %6s %8s %8s %8s %8s %8s %8s %6s" %
      ("arm", "cells", "E1/V0", "E1/V1", "E1/V2", "E2/V0", "E2/V1", "E2/V2", "vov_f"))
for a in ARMS:
    cs = data[a]
    if not cs:
        print("%-11s %6s  (no cells)" % (a, 0)); continue
    def cnt(w, t): return sum(1 for c in cs.values() if tier(c, w)[t])
    vov = sum(1 for c in cs.values() if c["vov_fail"])
    e1u = sum(1 for c in cs.values() if c["E1_rmsd"] is None)
    print("%-11s %6d %8s %8s %8s %8d %8d %8d %6d%s" %
          (a, len(cs),
           "n/a" if e1u else cnt("E1", "V0"),
           "n/a" if e1u else cnt("E1", "V1"),
           "n/a" if e1u else cnt("E1", "V2"),
           cnt("E2", "V0"), cnt("E2", "V1"), cnt("E2", "V2"), vov,
           "  (E1 unresolved in %d cells)" % e1u if e1u else ""))

def mcnemar(arm, ref, pred):
    A, B = data.get(arm, {}), data.get(ref, {})
    ts = sorted(set(A) & set(B))
    b = sum(1 for t in ts if pred(A[t]) and not pred(B[t]))
    c = sum(1 for t in ts if pred(B[t]) and not pred(A[t]))
    return b, c, len(ts), p2(b, c)

print("")
print("PAIRED TESTS vs B_shrink  (CONFIRMATORY: T4_all only; T1-T3 are EXPLORATORY)")
for a in ["T1_strain", "T2_walcap", "T3_solvref", "T4_all"]:
    role = "CONFIRMATORY" if a == "T4_all" else "EXPLORATORY "
    b, c, n, p = mcnemar(a, "B_shrink", lambda x: not x["vov_fail"])
    print("  %s %-11s P1 interpenetration : b=%2d c=%2d n=%2d p=%.4g" % (role, a, b, c, n, p))
    b, c, n, p = mcnemar(a, "B_shrink", lambda x: tier(x, "E2")["V0"])
    print("  %s %-11s P2 top-1 (E2)       : b=%2d c=%2d n=%2d p=%.4g" % (role, a, b, c, n, p))
    unres = sum(1 for c in list(data.get(a, {}).values()) + list(data.get("B_shrink", {}).values())
                if c["E1_rmsd"] is None)
    if E1_STATUS == "OK" and unres == 0:
        b, c, n, p = mcnemar(a, "B_shrink", lambda x: tier(x, "E1")["V0"])
        print("  %s %-11s P2 top-1 (E1 PRIM.) : b=%2d c=%2d n=%2d p=%.4g" % (role, a, b, c, n, p))
    else:
        print("  %s %-11s P2 top-1 (E1 PRIM.) : WITHHELD (%d cells unresolved)" % (role, a, unres))

print("")
print("PRE-REGISTERED VERDICT (T4_all, E1, crystal frame)")
# FAIL-CLOSED. An unresolved E1 must never be scored as a MISS -- that is exactly
# how a null measurement gets reported as a refutation. Withhold instead.
def e1_complete(arm):
    cs = data.get(arm, {})
    return bool(cs) and all(c["E1_rmsd"] is not None for c in cs.values())
E1_READY = E1_STATUS == "OK" and e1_complete("T4_all") and e1_complete("B_shrink")
if not E1_READY:
    why = E1_STATUS if E1_STATUS != "OK" else "unresolved in %d T4_all / %d B_shrink cells" % (
        sum(1 for c in data.get("T4_all", {}).values() if c["E1_rmsd"] is None),
        sum(1 for c in data.get("B_shrink", {}).values() if c["E1_rmsd"] is None))
    print("  WITHHELD: E1 is %s. E2 is a DIFFERENT ENDPOINT and may not be" % why)
    print("  substituted for it (pre-registration section 7). Resolve E1, then re-score.")
    print("  An unresolved E1 is NOT a miss and is NOT a refutation.")
elif N != 84:
    print("  WITHHELD: the pre-registered thresholds (>=45/84, <=8/77) are stated for")
    print("  N=84 and 77 PoseBusters cells. This batch has N=%d. Re-derive and" % N)
    print("  re-freeze the pre-registration before any verdict is issued.")
else:
    t4 = data["T4_all"]
    vov = sum(1 for c in t4.values() if c["vov_fail"])
    s1 = sum(1 for c in t4.values() if tier(c, "E1")["V0"])
    b1, c1, _, _ = mcnemar("T4_all", "B_shrink", lambda x: not x["vov_fail"])
    b2, c2, _, _ = mcnemar("T4_all", "B_shrink", lambda x: tier(x, "E1")["V0"])
    P1 = (vov <= 8) and (b1 >= 30) and (c1 <= 2)
    P2 = (s1 >= 45) and (b2 >= 10) and (c2 <= 2)
    print("  P1 interpenetration <=8/77 and b>=30,c<=2 : %s (vov=%d b=%d c=%d)" % (P1, vov, b1, c1))
    print("  P2 top-1 >=45/%d and b>=10,c<=2           : %s (S1=%d b=%d c=%d)" % (N, P2, s1, b2, c2))
    if P1 and P2:
        print("  => CONFIRMED. H1 supported.")
    elif P1 and not P2:
        print("  => F1: REFUTED as a single cause. Two separate problems. Do NOT")
        print("     default-enable the gates and do NOT report this as an accuracy result.")
    elif P2 and not P1:
        print("  => F2: REFUTED. Accuracy moved without the named mechanism. This may NOT")
        print("     be reported as a solvation/strain finding.")
    elif (8 < vov <= 30) or (38 < s1 < 45):
        print("  => INCONCLUSIVE (pre-declared ambiguous zone). No threshold may be moved.")
    else:
        print("  => F3: REFUTED outright.")
print("")
print("Tiers: V0 = RMSD<=2 only (comparable to the archive's 48/84 and 32/84).")
print("       V1 = V0 and PoseBusters-as-run (water checks INCLUDED).")
print("       V2 = V0 and PoseBusters excluding the two crystallographic-water")
print("            checks, which are out of protocol scope for an IMPLICIT-solvent")
print("            model. V2 is reported ALONGSIDE V1, never instead of it.")
PYEOF
  verdict_write G6_scoring PASS "results=$BATCH/RESULTS.txt" "E1_status=$E1_STATUS"
fi

say "[prereg] done. verdicts:"
for f in "$BATCH"/verdict/*; do [ -f "$f" ] && say "  $(basename "$f") = $(grep '^VERDICT=' "$f" | cut -d= -f2)"; done
say "[prereg] batch = $BATCH"
