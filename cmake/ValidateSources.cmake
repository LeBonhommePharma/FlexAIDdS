# ValidateSources.cmake
#
# FlexAIDdS — Bulletproof source file guard for CMake builds.
#
# Purpose:
#   Prevent the classic "I added Foo.cpp but forgot to list it in any target"
#   bug that has already caused link failures (GrandPartitionFunction, etc.).
#
# This module provides a single entry point:
#   flexaids_validate_sources([STRICT] [WARN_ONLY] [PYTHON_SCRIPT path])
#
# It executes the Python guard script (scripts/validate_sources.py) at
# configure time and can make configuration fail on orphaned sources.
#
# Recommended usage (near the bottom of the root CMakeLists.txt, after
# all add_subdirectory and target_sources calls):
#
#   include(cmake/ValidateSources.cmake)
#   flexaids_validate_sources()                    # warn only (good for migration)
#   # flexaids_validate_sources(STRICT)            # hard fail (recommended for CI)
#
# The guard also works when doing pure-Python development:
#   python -m pip install -e .
# because setup.py can call the same Python script directly.
#
# Migration:
#   During transition, use WARN_ONLY or provide an ignore file.
#   Once the codebase is clean, switch the top-level call to STRICT.
#
# See also: scripts/validate_sources.py for the actual implementation
# and ignore list management.

cmake_minimum_required(VERSION 3.28)

# ------------------------------------------------------------------------------
# flexaids_validate_sources
# ------------------------------------------------------------------------------
function(flexaids_validate_sources)
    set(options STRICT WARN_ONLY)
    set(oneValueArgs PYTHON_SCRIPT IGNORE_FILE)
    set(multiValueArgs EXTRA_SCAN_DIRS)
    cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Locate the Python guard script
    if(ARGS_PYTHON_SCRIPT)
        set(guard_script "${ARGS_PYTHON_SCRIPT}")
    else()
        set(guard_script "${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_sources.py")
    endif()

    if(NOT EXISTS "${guard_script}")
        message(WARNING "FlexAID source validator not found at ${guard_script} — skipping guard")
        return()
    endif()

    # Find Python
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(NOT Python3_Interpreter_FOUND)
        message(WARNING "Python3 not found — cannot run FlexAID source validator. "
                        "Install Python or set Python3_EXECUTABLE.")
        return()
    endif()

    # Build command line for the guard
    set(guard_cmd
        ${Python3_EXECUTABLE} "${guard_script}"
        --root "${CMAKE_CURRENT_SOURCE_DIR}"
        --cmake-mode
    )

    if(ARGS_IGNORE_FILE)
        list(APPEND guard_cmd --ignore-file "${ARGS_IGNORE_FILE}")
    else()
        # Default ignore file location
        set(default_ignore "${CMAKE_CURRENT_SOURCE_DIR}/build_sources.ignore")
        if(EXISTS "${default_ignore}")
            list(APPEND guard_cmd --ignore-file "${default_ignore}")
        endif()
    endif()

    if(ARGS_EXTRA_SCAN_DIRS)
        foreach(dir ${ARGS_EXTRA_SCAN_DIRS})
            list(APPEND guard_cmd --extra-scan-dir "${dir}")
        endforeach()
    endif()

    # Mode selection
    if(ARGS_STRICT AND NOT ARGS_WARN_ONLY)
        list(APPEND guard_cmd --strict)
        set(failure_message "Orphaned source files detected. Configuration aborted.")
        set(failure_level FATAL_ERROR)
    else()
        list(APPEND guard_cmd --warn-only)
        set(failure_message "Orphaned source files detected (non-fatal during migration).")
        set(failure_level WARNING)
    endif()

    # Run the guard
    message(STATUS "Running FlexAID source validator (scripts/validate_sources.py)...")

    execute_process(
        COMMAND ${guard_cmd}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE guard_result
        OUTPUT_VARIABLE guard_output
        ERROR_VARIABLE  guard_error
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(NOT guard_result EQUAL 0)
        if(guard_output)
            message(${failure_level} "${guard_output}")
        endif()
        if(guard_error)
            message(${failure_level} "${guard_error}")
        endif()
        message(${failure_level} "${failure_message}")
    else()
        if(guard_output)
            # The script prints a nice one-line summary on success
            message(STATUS "${guard_output}")
        else()
            message(STATUS "FlexAID source validator: no orphaned files found.")
        endif()
    endif()
endfunction()
