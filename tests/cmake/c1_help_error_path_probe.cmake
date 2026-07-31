# C1: error-path leak probe for the shipping FlexAID binary.
#
# Design (board, 2026-07-31):
#   Terminate() throws (LIB/fileio.cpp); main catches after the teardown
#   (LIB/top.cpp). --help takes Terminate(0) and exits 0 having freed nothing.
#   Under Linux ASan+LSan that should report the prologue allocations.
#
# Criterion is INVERTED from a normal leak test:
#   - process path always checked (--help -> usage text; exit 0 when not LSan-aborting)
#   - when EXPECT_LEAK=1 (Linux instrumented builds only): a LeakSanitizer
#     report is a PASS (the finding). Clean under LSan is FAIL — source model
#     wrong, or someone fixed error-path free without updating this probe.
#
# Structural LSan prediction (Bumble, source + reachability — unmeasured until
# Linux CI runs this):
#   Direct:   FA, GB, VC                          (3 roots in main)
#   Indirect: FA->contacts + 5 VC-> arrays        (6, reachable from roots)
#   Total:    9 objects  (Direct + Indirect)
#   Bytes:    DO NOT assert — macOS size classes ≠ glibc; 1,274,880 is mac-only.
#
# A clean LSan run, a total other than 9, or a split other than 3/6 is a finding
# about the source model / ownership graph, not a silent pass.
#
# FORBIDDEN claim: "this test covers top.cpp teardown." It never reaches it.
# Floor-not-ceiling: these nine are the early-exit floor; late Terminate mid-dock
# is unmeasured and is C2's problem (or a future late-error probe).
#
# Inputs ( -D from add_test ):
#   FlexAID_BINARY  absolute path to FlexAID
#   EXPECT_LEAK     ON/OFF — require LSan report with 3/6/9 structure

if(NOT DEFINED FlexAID_BINARY)
  message(FATAL_ERROR "c1_help_error_path_probe: FlexAID_BINARY not set")
endif()
if(NOT EXISTS "${FlexAID_BINARY}")
  message(FATAL_ERROR "c1_help_error_path_probe: binary missing: ${FlexAID_BINARY}")
endif()

if(NOT DEFINED EXPECT_LEAK)
  set(EXPECT_LEAK OFF)
endif()

set(_out_dir "${CMAKE_CURRENT_BINARY_DIR}/c1_help_probe")
file(MAKE_DIRECTORY "${_out_dir}")
set(_stdout "${_out_dir}/help.stdout")
set(_stderr "${_out_dir}/help.stderr")

# Prefer explicit ASAN_OPTIONS from the environment (CI sets detect_leaks=1).
# When EXPECT_LEAK, force detect_leaks so a bare ctest run still probes.
if(EXPECT_LEAK)
  set(ENV{ASAN_OPTIONS} "detect_leaks=1:halt_on_error=0:strict_string_checks=1")
endif()

execute_process(
  COMMAND "${FlexAID_BINARY}" --help
  OUTPUT_FILE "${_stdout}"
  ERROR_FILE  "${_stderr}"
  RESULT_VARIABLE _rc
)

file(READ "${_stdout}" _out_txt)
file(READ "${_stderr}" _err_txt)

# macOS AppleClang: detect_leaks is unsupported and aborts. Re-run without it
# and only enforce the process path.
if(_err_txt MATCHES "detect_leaks is not supported")
  message(STATUS "C1 probe: LSan unsupported on this platform — process path only")
  set(ENV{ASAN_OPTIONS} "detect_leaks=0:halt_on_error=0")
  execute_process(
    COMMAND "${FlexAID_BINARY}" --help
    OUTPUT_FILE "${_stdout}"
    ERROR_FILE  "${_stderr}"
    RESULT_VARIABLE _rc
  )
  file(READ "${_stdout}" _out_txt)
  file(READ "${_stderr}" _err_txt)
  set(EXPECT_LEAK OFF)
endif()

# Usage text must appear on the --help path (stdout). LSan may force a
# non-zero process exit when EXPECT_LEAK — do not require exit 0 first.
if(NOT _out_txt MATCHES "Usage:" AND NOT _out_txt MATCHES "--help")
  message(FATAL_ERROR
    "C1 probe: --help did not print usage text (exit=${_rc})\n"
    "stdout:\n${_out_txt}\nstderr:\n${_err_txt}")
endif()

if(NOT EXPECT_LEAK)
  # Non-instrumented / non-Linux / TSan: process path only — must be clean exit 0.
  # Machine line C1_RESULT=… is what PASS_REGULAR_EXPRESSION keys on so a
  # silent demotion cannot share the green signature of the LSan arm
  # (Honey: "the check that disarms itself").
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR
      "C1_RESULT=FAIL ARM=process_path_only\n"
      "C1 probe: --help expected exit 0, got ${_rc}\n"
      "stdout:\n${_out_txt}\nstderr:\n${_err_txt}")
  endif()
  # Single-line machine receipt first (ctest PASS_REGULAR_EXPRESSION).
  message(STATUS "C1_RESULT=PASS ARM=process_path_only")
  message(STATUS
    "C1 probe ARM=process_path_only: PASS (exit 0, usage present). "
    "LSan inverted arm NOT armed on this config — do not read this green "
    "as 'leaks measured' (TSan / macOS / non-ASan demote here).")
  return()
endif()

# ─── Inverted LSan arm ─────────────────────────────────────────────────────
# Leak report is the expected finding today. LSan typically exits non-zero when
# it reports; that is still a probe PASS if the structure matches prediction.
message(STATUS "C1 probe ARM=lsan_inverted: requiring LeakSanitizer 3 direct / 6 indirect / 9 total")

set(_saw_leak FALSE)
if(_err_txt MATCHES "LeakSanitizer" OR _err_txt MATCHES "detected memory leaks")
  set(_saw_leak TRUE)
endif()

if(NOT _saw_leak)
  # Two clean causes (Bumble, 2026-07-31) — not equally informative:
  #   (1) teardown actually ran / prologue freed  -> falsifies unwind read
  #   (2) stale-root false negative (conservative stack scan keeps FA/GB/VC
  #       "still reachable") -> instrument artefact, not ownership model
  # Only total-9 with wrong split falsifies Bumble's ownership graph.
  # Read the log: nine "still reachable" = stale-root signature; nine absent
  # entirely is not. Do not collapse those into one narrative.
  message(FATAL_ERROR
    "C1_RESULT=FAIL ARM=lsan_inverted\n"
    "C1 probe: Linux ASan+LSan --help reported no definitely-lost leaks.\n"
    "Disambiguate before blaming the ownership model:\n"
    "  - still-reachable count ~9 under LSAN_OPTIONS=report_objects=1 "
    "    -> stale-root false negative (measurability, not wrong model)\n"
    "  - nothing still-reachable and no leaks "
    "    -> Terminate-unwind / free-before-throw / LSan-not-running "
    "    (interesting finding about top.cpp)\n"
    "  - total 9 with split != 3/6 "
    "    -> ownership graph wrong (Bumble's read fails)\n"
    "exit=${_rc}\nstderr:\n${_err_txt}")
endif()

# Sum object counts from Direct/Indirect lines — not header counts alone.
# Real LSan (board samples): "Direct leak of N byte(s) in M object(s)"
# Coalesced stacks produce one header with M>1; counting headers only would
# report 1/1 on a correct 3/6 ownership graph and fail a healthy binary.
# That is the inverted vacuous shape Bumble/Honey called out for "expected 9"
# against the Direct block alone.
#
# Prediction (Bumble, source + LSan reachability — unmeasured until Linux CI):
#   3 direct roots (FA/GB/VC) + 6 indirect (contacts + 5 VC arrays) = 9.
set(_n_direct 0)
string(REGEX MATCHALL
  "Direct leak of [0-9]+ byte\\(s\\) in [0-9]+ object\\(s\\)?"
  _direct_lines "${_err_txt}")
foreach(_line IN LISTS _direct_lines)
  if(_line MATCHES "in ([0-9]+) object")
    math(EXPR _n_direct "${_n_direct} + ${CMAKE_MATCH_1}")
  endif()
endforeach()

set(_n_indirect 0)
string(REGEX MATCHALL
  "Indirect leak of [0-9]+ byte\\(s\\) in [0-9]+ object\\(s\\)?"
  _indirect_lines "${_err_txt}")
foreach(_line IN LISTS _indirect_lines)
  if(_line MATCHES "in ([0-9]+) object")
    math(EXPR _n_indirect "${_n_indirect} + ${CMAKE_MATCH_1}")
  endif()
endforeach()

math(EXPR _n_sum "${_n_direct} + ${_n_indirect}")

# SUMMARY allocation count when present (must reconcile with Direct+Indirect).
set(_summary_n "")
if(_err_txt MATCHES "leaked in ([0-9]+) allocation")
  set(_summary_n "${CMAKE_MATCH_1}")
endif()

# Require a parsed split. SUMMARY-only is not enough — it can pass for the
# wrong ownership graph (9/0, 4/5) and fails the handoff Honey just fixed.
set(_split_ok FALSE)
if(_n_direct EQUAL 3 AND _n_indirect EQUAL 6)
  set(_split_ok TRUE)
endif()

set(_total_ok FALSE)
if(_n_sum EQUAL 9)
  set(_total_ok TRUE)
endif()
if(NOT _summary_n STREQUAL "" AND _summary_n STREQUAL "9")
  set(_total_ok TRUE)
endif()

if(_split_ok AND _total_ok)
  # Cross-check SUMMARY against the split when both are present.
  if(NOT _summary_n STREQUAL "" AND NOT _summary_n STREQUAL "${_n_sum}")
    message(FATAL_ERROR
      "C1_RESULT=FAIL ARM=lsan_inverted\n"
      "C1 probe: SUMMARY allocation count (${_summary_n}) != "
      "Direct+Indirect (${_n_sum}=${_n_direct}+${_n_indirect}).\n"
      "Parser or LSan format disagreement — fix before trusting either number.\n"
      "exit=${_rc}\nstderr:\n${_err_txt}")
  endif()
  # Machine receipt first — ctest PASS_REGULAR_EXPRESSION requires this line.
  message(STATUS
    "C1_RESULT=PASS ARM=lsan_inverted direct=${_n_direct} indirect=${_n_indirect} "
    "summary_allocs=${_summary_n}")
  message(STATUS
    "C1 probe FINDING (expected under current source model): "
    "--help under LSan: ${_n_direct} direct / ${_n_indirect} indirect "
    "(summary_allocs=${_summary_n}, exit=${_rc}). "
    "Teardown skipped via Terminate. NOT production-teardown coverage. "
    "Bytes intentionally not asserted (cross-allocator).")
  return()
endif()

# Total 9 but wrong split: ownership graph is not what the source says.
if(_total_ok AND NOT _split_ok)
  message(FATAL_ERROR
    "C1_RESULT=FAIL ARM=lsan_inverted\n"
    "C1 probe: LSan object total is 9 but Direct/Indirect split is "
    "${_n_direct}/${_n_indirect} (expected 3/6).\n"
    "Reachability/ownership graph differs from FA/GB/VC roots + 6 children. "
    "Do not silence this — it is a real finding about the source model "
    "(9/0 or 4/5 falsifies who owns what).\n"
    "summary_allocs=${_summary_n} exit=${_rc}\nstderr:\n${_err_txt}")
endif()

# Unparseable split (0/0) with only SUMMARY: refuse — that hands off a test
# that cannot distinguish Direct-only from Direct+Indirect.
if(_n_direct EQUAL 0 AND _n_indirect EQUAL 0)
  message(FATAL_ERROR
    "C1_RESULT=FAIL ARM=lsan_inverted\n"
    "C1 probe: LSan reported leaks but Direct/Indirect object counts could "
    "not be parsed (summary_allocs=${_summary_n}).\n"
    "Refusing SUMMARY-only pass — that is the shape that fails on correct "
    "behaviour if the assertion ever reads the Direct line alone.\n"
    "Fix the parser against real LSan output; do not drop the 3/6 assertion.\n"
    "exit=${_rc}\nstderr:\n${_err_txt}")
endif()

# Anything else: wrong total or wrong split.
message(FATAL_ERROR
  "C1_RESULT=FAIL ARM=lsan_inverted\n"
  "C1 probe: LSan reported leaks but not the predicted 9 objects (3 direct / 6 indirect).\n"
  "parsed: direct=${_n_direct} indirect=${_n_indirect} sum=${_n_sum} "
  "summary_allocs=${_summary_n} (expected total 9, split 3/6).\n"
  "Frozen in PLANS/C2_ACCEPTANCE_CRITERION.md: MUST report LSan block; "
  "MUST sum Direct+Indirect to 9; MUST split 3/6; bytes DO NOT compare.\n"
  "exit=${_rc}\nstderr:\n${_err_txt}")
