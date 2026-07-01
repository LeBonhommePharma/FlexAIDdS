# cmake/FlexAIDOptions.cmake
#
# FlexAIDdS — Centralized build options, platform detection, and helper functions.
#
# This module is intended to be included early (right after basic project()
# and language setup, before dependency discovery and add_subdirectory).
#
# Benefits:
#   - Single place for all high-level FLEXAIDS_* toggles
#   - Cleaner root CMakeLists.txt
#   - Easier to document + maintain
#   - Prepares for future components / package config
#
# Usage in root:
#   include(cmake/FlexAIDOptions.cmake)
#
# See also:
#   cmake/FlexAIDDependencies.cmake (planned)
#   cmake/FlexAIDComponents.cmake (planned)
#   cmake/ValidateSources.cmake (already used)

cmake_minimum_required(VERSION 3.28)

# ─────────────────────────────────────────────────────────────────────────────
# Core build options (declared as early as practical)
# ─────────────────────────────────────────────────────────────────────────────

# Code coverage (must be early because it may inject global compile options)
option(FLEXAIDS_ENABLE_COVERAGE "Enable code coverage tracking (requires gcov/lcov)" OFF)

if(FLEXAIDS_ENABLE_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(WARNING "Code coverage is only supported with GCC or Clang")
        set(FLEXAIDS_ENABLE_COVERAGE OFF)
    else()
        message(STATUS "Code coverage enabled — building with -fprofile-arcs -ftest-coverage")
        add_compile_options(-fprofile-arcs -ftest-coverage)
        add_link_options(--coverage)
    endif()
endif()

# GPU / acceleration
option(FLEXAIDS_USE_CUDA    "Enable CUDA GPU evaluation"          OFF)
option(FLEXAIDS_USE_ROCM    "Enable ROCm/HIP GPU evaluation"      OFF)
if(APPLE AND _HAS_OBJCXX)   # _HAS_OBJCXX is set in root before include
    option(FLEXAIDS_USE_METAL   "Enable Metal GPU acceleration"       ON)
else()
    option(FLEXAIDS_USE_METAL   "Enable Metal GPU acceleration"       OFF)
endif()

# SIMD
option(FLEXAIDS_USE_AVX2    "Enable AVX2 SIMD acceleration"       ON)
option(FLEXAIDS_USE_AVX512  "Enable AVX-512 SIMD acceleration"    OFF)
option(FLEXAIDS_USE_SOA_DISTANCES "Route Voronoi hot-path distances through AtomSoA float SoA arrays (C1)" OFF)

# Parallelism & core libs
option(FLEXAIDS_USE_OPENMP  "Enable OpenMP thread parallelism"    ON)
# Eigen3 is a hard requirement — no option to disable
set(FLEXAIDS_USE_EIGEN ON CACHE BOOL "Eigen3 is required" FORCE)
option(FLEXAIDS_USE_256_MATRIX "Enable 256×256 soft contact matrix" ON)
option(FLEXAIDS_USE_MPI     "Enable MPI distributed parallel docking" OFF)

# High-level build toggles
option(BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" OFF)
option(FLEXAIDS_BUILD_CORE  "Build core executables (FlexAID/FlexAIDdS/tENCoM)" ON)

# ─────────────────────────────────────────────────────────────────────────────
# Platform / architecture detection
# ─────────────────────────────────────────────────────────────────────────────

string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" FLEXAIDS_SYSTEM_PROCESSOR_LOWER)
set(FLEXAIDS_IS_ARM FALSE)
set(FLEXAIDS_IS_X86 FALSE)

if(FLEXAIDS_SYSTEM_PROCESSOR_LOWER MATCHES "^(arm64|aarch64)$")
    set(FLEXAIDS_IS_ARM TRUE)
elseif(FLEXAIDS_SYSTEM_PROCESSOR_LOWER MATCHES "^(x86_64|amd64|x64|i[3-6]86)$")
    set(FLEXAIDS_IS_X86 TRUE)
endif()

if(APPLE AND FLEXAIDS_IS_ARM)
    if(FLEXAIDS_USE_AVX2 OR FLEXAIDS_USE_AVX512)
        message(STATUS "Apple Silicon detected — disabling x86 SIMD flags (-mavx2, -mfma, -mavx512*)")
    endif()
    set(FLEXAIDS_USE_AVX2 OFF CACHE BOOL "Enable AVX2 SIMD acceleration" FORCE)
    set(FLEXAIDS_USE_AVX512 OFF CACHE BOOL "Enable AVX-512 SIMD acceleration" FORCE)
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Helper: apply architecture-appropriate SIMD / SoA flags to a target
# ─────────────────────────────────────────────────────────────────────────────
function(flexaids_configure_simd target_name)
    if(MSVC)
        if(FLEXAIDS_IS_X86)
            if(FLEXAIDS_USE_AVX512)
                target_compile_options(${target_name} PRIVATE /arch:AVX512)
                target_compile_definitions(${target_name} PRIVATE FLEXAIDS_USE_AVX512)
                message(STATUS "${target_name}: AVX-512 enabled via /arch:AVX512")
            elseif(FLEXAIDS_USE_AVX2)
                target_compile_options(${target_name} PRIVATE /arch:AVX2)
                target_compile_definitions(${target_name} PRIVATE FLEXAIDS_USE_AVX2)
                message(STATUS "${target_name}: AVX2 enabled via /arch:AVX2")
            endif()
        else()
            message(STATUS "${target_name}: non-x86 MSVC target detected — building without x86 SIMD flags")
        endif()
        return()
    endif()

    if(NOT FLEXAIDS_IS_X86)
        if(FLEXAIDS_IS_ARM)
            # ARM64 (Apple Silicon / aarch64): NEON is part of the baseline ISA,
            # no special -m flag required. Advertise the NEON capability so the
            # SoA distance hot path (C1) can select distance2_1x4_neon (B1).
            target_compile_definitions(${target_name} PRIVATE
                FLEXAIDS_USE_NEON FLEXAIDS_SIMD_NEON)
            message(STATUS "${target_name}: ARM NEON enabled (baseline aarch64 ISA)")
        else()
            message(STATUS "${target_name}: non-x86 target detected (${CMAKE_SYSTEM_PROCESSOR}) — building without SIMD flags")
        endif()
        if(FLEXAIDS_USE_SOA_DISTANCES)
            target_compile_definitions(${target_name} PRIVATE FLEXAIDS_USE_SOA_DISTANCES)
            message(STATUS "${target_name}: SoA distance routing enabled (FLEXAIDS_USE_SOA_DISTANCES)")
        endif()
        return()
    endif()

    if(FLEXAIDS_USE_AVX512)
        check_cxx_compiler_flag("-mavx512f" FLEXAIDS_COMPILER_SUPPORTS_MAVX512F)
        check_cxx_compiler_flag("-mavx512dq" FLEXAIDS_COMPILER_SUPPORTS_MAVX512DQ)
        check_cxx_compiler_flag("-mavx512bw" FLEXAIDS_COMPILER_SUPPORTS_MAVX512BW)
        check_cxx_compiler_flag("-mavx2" FLEXAIDS_COMPILER_SUPPORTS_MAVX2)
        check_cxx_compiler_flag("-mfma" FLEXAIDS_COMPILER_SUPPORTS_MFMA)

        if(FLEXAIDS_COMPILER_SUPPORTS_MAVX512F AND FLEXAIDS_COMPILER_SUPPORTS_MAVX512DQ AND
           FLEXAIDS_COMPILER_SUPPORTS_MAVX512BW AND FLEXAIDS_COMPILER_SUPPORTS_MAVX2 AND
           FLEXAIDS_COMPILER_SUPPORTS_MFMA)
            target_compile_options(${target_name} PRIVATE -mavx512f -mavx512dq -mavx512bw -mavx2 -mfma)
            target_compile_definitions(${target_name} PRIVATE FLEXAIDS_USE_AVX512)
            message(STATUS "${target_name}: AVX-512 enabled")
        else()
            message(WARNING "${target_name}: AVX-512 requested but compiler/target does not support the required flags; building without AVX-512")
        endif()
    elseif(FLEXAIDS_USE_AVX2)
        check_cxx_compiler_flag("-mavx2" FLEXAIDS_COMPILER_SUPPORTS_MAVX2)
        check_cxx_compiler_flag("-mfma" FLEXAIDS_COMPILER_SUPPORTS_MFMA)

        if(FLEXAIDS_COMPILER_SUPPORTS_MAVX2 AND FLEXAIDS_COMPILER_SUPPORTS_MFMA)
            target_compile_options(${target_name} PRIVATE -mavx2 -mfma)
            target_compile_definitions(${target_name} PRIVATE FLEXAIDS_USE_AVX2)
            message(STATUS "${target_name}: AVX2/FMA enabled")
        else()
            message(WARNING "${target_name}: AVX2 requested but compiler/target does not support -mavx2/-mfma; building without x86 SIMD")
        endif()
    endif()

    if(FLEXAIDS_USE_SOA_DISTANCES)
        target_compile_definitions(${target_name} PRIVATE FLEXAIDS_USE_SOA_DISTANCES)
        message(STATUS "${target_name}: SoA distance routing enabled (FLEXAIDS_USE_SOA_DISTANCES)")
    endif()
endfunction()

# ─────────────────────────────────────────────────────────────────────────────
# Small helper used by several tests (kept here for easy reuse)
# ─────────────────────────────────────────────────────────────────────────────
function(flexaids_configure_msvc_test target_name)
    if(MSVC)
        target_compile_definitions(${target_name} PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    endif()
endfunction()

message(STATUS "FlexAIDOptions loaded — options + SIMD helpers ready")
