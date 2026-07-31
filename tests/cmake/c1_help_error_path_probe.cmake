# C1: error-path leak probe for the shipping FlexAID binary.
#
# Design (board, 2026-07-31):
#   Terminate() throws (LIB/fileio.cpp); main catches after the teardown
#   (LIB/top.cpp). --help takes Terminate(0) and exits 0 having freed nothing.
#   Under Linux ASan+LSan that should report the prologue allocations.
#
# Criterion is INVERTED from a normal leak test:
#   - process path always checked (--help -> exit 0, usage text)
#   - when EXPECT_LEAK=1 (Linux instrumented builds only): a LeakSanitizer
#     report is a PASS (the finding). Clean under LSan is FAIL — source model
#     wrong, or someone fixed error-path free without updating this probe.
#
# FORBIDDEN claim: "this test covers top.cpp teardown." It never reaches it.
#
# Inputs ( -D from add_test ):
#   FlexAID_BINARY  absolute path to FlexAID
#   EXPECT_LEAK     ON/OFF — require LSan report

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

# Inverted LSan arm — leak report is the expected finding today.
# LSan typically exits non-zero when it reports; that is still a probe PASS.
set(_saw_leak FALSE)
if(_err_txt MATCHES "LeakSanitizer" OR _err_txt MATCHES "detected memory leaks" OR _err_txt MATCHES "Direct leak")
  set(_saw_leak TRUE)
endif()

if(_saw_leak)
  message(STATUS
    "C1 probe FINDING (expected under current source model): "
    "--help path reports leaks under LSan (teardown skipped via Terminate). "
    "exit=${_rc}. This is NOT production-teardown coverage.")
  return()
endif()

message(FATAL_ERROR
  "C1 probe UNEXPECTED: Linux ASan+LSan run of --help reported no leaks.\n"
  "Either the Terminate-unwind source model is wrong, prologue allocs are freed "
  "before throw, or LSan did not run. Investigate before claiming C1 green "
  "as anything other than this surprising clean result.\n"
  "exit=${_rc}\nstderr:\n${_err_txt}")
