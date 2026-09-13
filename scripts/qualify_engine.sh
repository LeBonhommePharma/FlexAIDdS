#!/usr/bin/env bash
# =============================================================================
# ENGINE QUALIFICATION GATE  --  a pilot must not launch against an unqualified
# engine. Six checks, each able to FAIL, emitted as a JSON receipt.
#
# WHY THIS EXISTS. On 2026-09-12 a pilot engine was built from branch tip
# 72f70557 and fully instrument-checked (gate reachable, data files staged,
# provenance stamped). It was nonetheless INVALID: commit f227375f, merged to
# main hours later, restructured the wall dispatch, and 72f70557 fails 4 of 17
# cases of the repo's own tests/test_wall_legacy_callsite.py -- in the
# zero-cutoff branch that feeds the pilot's CONTROL arm. Nothing in the process
# asked "has a later commit touched the engine files since this ref?".
# Check 2 is that question. It is the one that would have caught it.
#
# WHY CHECK 6 EXISTS. On 2026-09-13 this gate passed a tree that could not be
# compiled by CI. Commit 4e1043b6 ("wip(dsvib) ... NOT VALIDATED") added calls
# from LIB/ic2cf.cpp into the tENCoM entropy engine; tests/cf_aggregator_stubs.cpp
# did not provide those symbols, so test_cf_aggregator failed to LINK. Every CI
# build job on main went red for two days. This gate -- and every build run by
# hand that week -- configured with BUILD_TESTING=OFF, so the broken target was
# never built and nothing observed the failure. A gate that does not build what
# CI builds is not checking what CI checks.
#
# Usage:  qualify_engine.sh <build_dir> <src_dir> <ref> [receipt.json]
# Exit 0 = QUALIFIED, 1 = DISQUALIFIED, 2 = INDETERMINATE (a check could not run;
# never treat INDETERMINATE as a pass -- an unrunnable check is not a passing one).
set -u
BD=${1:?build dir}; SRC=${2:?source dir}; REF=${3:?ref}; OUT=${4:-$BD/qualification.json}
S=/Users/lp.more/Projects/FlexAIDdS
CONDA=/Users/lp.more/.claude-science/conda/envs/python/bin/python3
export GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_SYSTEM=/dev/null
ENGINE_FILES="LIB/vcfunction.cpp LIB/soft_wall.h LIB/flexaidds_flags.cpp LIB/ic2cf.cpp"
v=0; ind=0; J="{}"
add(){ J=$(printf '%s' "$J" | "$CONDA" -c "
import json,sys
d=json.load(sys.stdin); d['$1']=json.loads('''$2'''); print(json.dumps(d))"); }

# --- 1. the ref exists and is an ancestor of origin/main --------------------
cd "$S" || exit 2
git fetch -q origin 2>/dev/null
if git cat-file -e "${REF}^{commit}" 2>/dev/null; then
  anc=$(git merge-base --is-ancestor "$REF" origin/main 2>/dev/null && echo true || echo false)
else anc=null; fi
[ "$anc" = true ] || { [ "$anc" = null ] && ind=1 || v=1; }
add ancestry "{\"ref\":\"$REF\",\"is_ancestor_of_origin_main\":$anc,\"origin_main\":\"$(git rev-parse --short origin/main)\"}"

# --- 2. THE CHECK THAT WOULD HAVE CAUGHT IT --------------------------------
#     any commit AFTER this ref touching an engine file supersedes the build.
later=$(git rev-list --count "$REF"..origin/main -- $ENGINE_FILES 2>/dev/null || echo -1)
shas=$(git rev-list "$REF"..origin/main -- $ENGINE_FILES 2>/dev/null | head -5 | tr '\n' ' ')
[ "$later" = "-1" ] && ind=1
[ "${later:-0}" -gt 0 ] 2>/dev/null && v=1
add supersession "{\"engine_files\":\"$ENGINE_FILES\",\"later_commits_touching_them\":${later:--1},\"shas\":\"$shas\"}"

# --- 3. the feature is actually in the binary ------------------------------
gs=$(strings "$BD/FlexAIDdS" 2>/dev/null | grep -c 'FLEXAIDDS_WAL_C1' || echo 0)
df=$(ls "$BD"/AMINO.def "$BD"/MC_st0r5.2_6.dat "$BD"/NUCLEOTIDES.def "$BD"/rotobs.lst 2>/dev/null | wc -l | tr -d ' ')
stamp=$(strings "$BD/FlexAIDdS" 2>/dev/null | grep -oE '^[0-9a-f]{8}$' | head -1)
{ [ "$gs" -eq 0 ] || [ "$df" -ne 4 ]; } && v=1
add binary "{\"gate_string\":$gs,\"data_files\":$df,\"stamp\":\"${stamp:-}\",\"engine_sha256\":\"$(shasum -a256 "$BD/FlexAIDdS" 2>/dev/null | cut -d' ' -f1)\",\"runner_sha256\":\"$(shasum -a256 "$BD/benchmark_datasets" 2>/dev/null | cut -d' ' -f1)\"}"

# --- 4. the repo's own dispatch suite, against THIS source tree ------------
T=$SRC/tests/test_wall_legacy_callsite.py
[ -f "$T" ] || cp "$S/tests/test_wall_legacy_callsite.py" "$SRC/tests/" 2>/dev/null
if [ -f "$T" ]; then
  L=$BD.dispatch.log
  ( cd "$SRC" && CXX=c++ "$CONDA" -m pytest tests/test_wall_legacy_callsite.py -q --no-header \
       --rootdir="$SRC" -p no:cacheprovider > "$L" 2>&1 )
  p=$(grep -oE '[0-9]+ passed' "$L" | tail -1 | cut -d' ' -f1); p=${p:-0}
  f=$(grep -oE '[0-9]+ failed' "$L" | tail -1 | cut -d' ' -f1); f=${f:-0}
  e=$(grep -cE 'error in|^ERROR' "$L")
  [ "$e" -gt 0 ] && ind=1
  [ "$f" -gt 0 ] && v=1
  [ "$p" -eq 0 ] && ind=1
else p=0; f=0; e=1; ind=1; fi
add dispatch_suite "{\"passed\":$p,\"failed\":$f,\"errors\":$e,\"source_tree\":\"$SRC\"}"

# --- 5. the gated branch is reached by a real docking run ------------------
G=$BD.reach
if [ -f "$G/off/.setsha" ] && [ -f "$G/on/.setsha" ]; then
  so=$(cat "$G/off/.setsha"); sn=$(cat "$G/on/.setsha"); src=cached
else so=""; sn=""; src=absent; ind=1; fi
[ -n "$so" ] && [ "$so" = "$sn" ] && v=1
add reachability "{\"off_poseset\":\"$so\",\"on_poseset\":\"$sn\",\"diverged\":$([ -n "$so" ] && [ "$so" != "$sn" ] && echo true || echo false),\"source\":\"$src\"}"

# --- 6. THE TREE COMPILES THE WAY CI COMPILES IT ---------------------------
#     BUILD_TESTING=ON is the whole point: test fixtures link against a stub set
#     that production does not use, so a production-clean tree can still fail CI.
#     Set FLEXAIDS_QUALIFY_SKIP_BUILD=1 to skip -- that yields INDETERMINATE, never
#     a pass, because an unrunnable check is not a passing one.
tb=$BD.testbuild
if [ "${FLEXAIDS_QUALIFY_SKIP_BUILD:-0}" = "1" ]; then
  cfg_rc=-1; bld_rc=-1; ct_pass=0; ct_fail=0; ind=1; bnote=skipped_by_env
else
  # USE_SYSTEM_GTEST: FetchContent clones googletest from github at configure
  # time, which fails in a network-restricted sandbox. That is environmental, so
  # it must read INDETERMINATE rather than DISQUALIFIED -- distinguished below by
  # grepping the configure log for the clone failure specifically.
  cmake -S "$SRC" -B "$tb" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
        -DFLEXAIDS_USE_CUDA=OFF -DFLEXAIDS_USE_METAL=OFF -DBUILD_PYTHON_BINDINGS=OFF \
        -DFLEXAIDS_USE_SYSTEM_GTEST=ON > "$tb.cfg.log" 2>&1
  cfg_rc=$?
  if [ "$cfg_rc" -ne 0 ]; then
    if grep -q 'Failed to clone repository' "$tb.cfg.log" 2>/dev/null; then
      bld_rc=-1; ct_pass=0; ct_fail=0; ind=1; bnote=configure_blocked_network
    else
      bld_rc=-1; ct_pass=0; ct_fail=0; v=1; bnote=configure_failed
    fi
  else
    cmake --build "$tb" --parallel 3 > "$tb.build.log" 2>&1
    bld_rc=$?
    if [ "$bld_rc" -ne 0 ]; then
      v=1; ct_pass=0; ct_fail=0; bnote=build_failed
    else
      ctest --test-dir "$tb" --output-on-failure --timeout 240 > "$tb.ctest.log" 2>&1
      ct_pass=$(grep -oE '[0-9]+ tests passed' "$tb.ctest.log" | tail -1 | cut -d' ' -f1)
      ct_pass=${ct_pass:-0}
      ct_fail=$(grep -oE '[0-9]+ tests failed' "$tb.ctest.log" | tail -1 | cut -d' ' -f1)
      ct_fail=${ct_fail:-0}
      if [ "$ct_fail" -gt 0 ]; then v=1; bnote=tests_failed
      elif [ "$ct_pass" -eq 0 ]; then ind=1; bnote=no_tests_ran
      else bnote=ok; fi
    fi
  fi
fi
# the undefined-symbol signature specifically, so the receipt names the 4e1043b6
# failure mode rather than only its exit code
undef=$(grep -ciE 'Undefined symbols|undefined reference' "$tb.build.log" 2>/dev/null || echo 0)
add ci_parity_build "{\"build_testing\":\"ON\",\"configure_rc\":$cfg_rc,\"build_rc\":$bld_rc,\"ctest_passed\":$ct_pass,\"ctest_failed\":$ct_fail,\"undefined_symbol_lines\":$undef,\"note\":\"$bnote\"}"

verdict=QUALIFIED; rc=0
[ "$v" -eq 1 ] && { verdict=DISQUALIFIED; rc=1; }
[ "$v" -eq 0 ] && [ "$ind" -eq 1 ] && { verdict=INDETERMINATE; rc=2; }
add verdict "{\"verdict\":\"$verdict\",\"note\":\"INDETERMINATE is never a pass\"}"
printf '%s' "$J" | "$CONDA" -c 'import json,sys; print(json.dumps(json.load(sys.stdin), indent=1))' > "$OUT"
echo "$verdict"
exit $rc
