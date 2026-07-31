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
  # Non-instrumented / non-Linux: process path only — must be clean exit 0.
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR
      "C1 probe: --help expected exit 0, got ${_rc}\n"
      "stdout:\n${_out_txt}\nstderr:\n${_err_txt}")
  endif()
  message(STATUS "C1 probe: process path OK (exit 0, usage present); leak assertion not required on this config")
  return()
endif()

# ─── Inverted LSan arm ─────────────────────────────────────────────────────
# Leak report is the expected finding today. LSan typically exits non-zero when
# it reports; that is still a probe PASS if the structure matches prediction.

set(_saw_leak FALSE)
if(_err_txt MATCHES "LeakSanitizer" OR _err_txt MATCHES "detected memory leaks")
  set(_saw_leak TRUE)
endif()

if(NOT _saw_leak)
  message(FATAL_ERROR
    "C1 probe UNEXPECTED: Linux ASan+LSan run of --help reported no leaks.\n"
    "Either the Terminate-unwind source model is wrong, prologue allocs are freed "
    "before throw, or LSan did not run. Investigate before claiming C1 green "
    "as anything other than this surprising clean result.\n"
    "exit=${_rc}\nstderr:\n${_err_txt}")
endif()

# Count Direct / Indirect report headers (one per leaked object at distinct sites).
# macOS leaks reports a single total of 9; LSan splits reachability — the board
# prediction is 3 direct roots (FA/GB/VC) + 6 indirect (contacts + 5 VC arrays).
string(REGEX MATCHALL "Direct leak of" _direct_hits "${_err_txt}")
string(REGEX MATCHALL "Indirect leak of" _indirect_hits "${_err_txt}")
list(LENGTH _direct_hits _n_direct)
list(LENGTH _indirect_hits _n_indirect)
math(EXPR _n_sum "${_n_direct} + ${_n_indirect}")

# Prefer SUMMARY allocation count when present (allocator-independent object total).
set(_summary_n "")
if(_err_txt MATCHES "leaked in ([0-9]+) allocation")
  set(_summary_n "${CMAKE_MATCH_1}")
endif()

set(_total_ok FALSE)
if(_summary_n STREQUAL "9")
  set(_total_ok TRUE)
elseif(_n_sum EQUAL 9)
  set(_total_ok TRUE)
endif()

set(_split_ok FALSE)
if(_n_direct EQUAL 3 AND _n_indirect EQUAL 6)
  set(_split_ok TRUE)
endif()

if(_total_ok AND _split_ok)
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
    "C1 probe: LSan object total is 9 but Direct/Indirect split is "
    "${_n_direct}/${_n_indirect} (expected 3/6).\n"
    "Reachability/ownership graph differs from FA/GB/VC roots + 6 children. "
    "Do not silence this — it is a real finding about the source model.\n"
    "summary_allocs=${_summary_n} exit=${_rc}\nstderr:\n${_err_txt}")
endif()

# Anything else: wrong total or unparseable.
message(FATAL_ERROR
  "C1 probe: LSan reported leaks but not the predicted 9 objects (3 direct / 6 indirect).\n"
  "parsed: direct=${_n_direct} indirect=${_n_indirect} sum=${_n_sum} "
  "summary_allocs=${_summary_n} (expected total 9, split 3/6).\n"
  "If summary is empty, LSan output format may have changed — fix the parser, "
  "do not drop the structural assertion.\n"
  "exit=${_rc}\nstderr:\n${_err_txt}")
