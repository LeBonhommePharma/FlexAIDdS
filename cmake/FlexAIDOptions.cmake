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
# Default OFF per METHODOLOGY §1 ("any intended behavior change must be opt-in
# behind a flag that defaults OFF, and parity must hold with the flag OFF").
# This path is ranking-affecting, not merely faster: float32 squared distances
# perturb contlist.dist and the near/clash cutoffs, hence Voronoi topology.
# It also could not compile on production x86 until this PR: simd::distance2_1x4
# was defined only in the SSE4.2/NEON/scalar branches of simd_distance.h and was
# missing from the AVX2 and AVX-512 branches. That is now fixed at source (see
# the 4-wide overload added after the ISA branches in simd_distance.h), with
# atom_soa::distance2_1x4 kept as the ISA-portable call site. The reasoning for
# defaulting OFF is unchanged and does not depend on the compile break: enabling
# it would ship a code path no prior CI on those machines has ever executed.
# Turn ON explicitly to benchmark it, and run the Astex-85 A/B before proposing
# it as the default.
option(FLEXAIDS_USE_SOA_DISTANCES "Route Voronoi hot-path distances through AtomSoA float SoA arrays (C1)" OFF)

# ─── Native CPU tuning for flexaid_core (perf vs. portability trade-off) ──
# -mcpu=native (Apple/Clang, arm64) tunes instruction scheduling/selection
# for the exact CPU that runs the build. flexaid_core already compiles with
# -ffast-math; -mcpu=native can additionally shift FMA contraction and
# vectorization choices relative to a generic-aarch64 build, so results are
# NOT guaranteed bit-identical across CPU tunings even though this flag does
# not itself relax IEEE semantics the way -ffast-math does. The resulting
# binary is also NOT portable: it may refuse to run, or silently use a
# codegen path unverified on a different Apple Silicon SKU.
# Default OFF to match BUILD_FLEXAID_FAST and METHODOLOGY §1 (engine SHA
# comparable across machines). Set ON explicitly for local perf work on the
# machine that will also run the binary.
option(FLEXAIDS_MCPU_NATIVE "Compile flexaid_core with -mcpu=native on Apple/Clang arm64 (perf; not portable across machines, may change FP codegen under -ffast-math)" OFF)

# Parallelism & core libs
option(FLEXAIDS_USE_OPENMP  "Enable OpenMP thread parallelism"    ON)
# Eigen3 is a hard requirement — no option to disable
set(FLEXAIDS_USE_EIGEN ON CACHE BOOL "Eigen3 is required" FORCE)
option(FLEXAIDS_USE_256_MATRIX "Enable 256×256 soft contact matrix" ON)
option(FLEXAIDS_USE_MPI     "Enable MPI distributed parallel docking" OFF)

# High-level build toggles
option(BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" OFF)
option(FLEXAIDS_BUILD_CORE  "Build core executables (FlexAID/FlexAIDdS/tENCoM)" ON)

# Grand canonical partition function (Ξ) for competitive binding (P0 foundation).
# Default ON so existing behavior is unchanged. Sources for GrandPartitionFunction,
# TargetServer and MultiSiteGPF are wired to be compiled generally and always
# when BUILD_TESTING=ON (see root CMakeLists.txt force + LIB/CMakeLists.txt).
option(FLEXAIDS_GRAND_CANONICAL "Enable Grand Canonical Partition Function (Ξ) + TargetServer/MultiSiteGPF for competitive binding" ON)
if(FLEXAIDS_GRAND_CANONICAL)
    message(STATUS "FLEXAIDS_GRAND_CANONICAL=ON — GrandPartitionFunction/TargetServer/MultiSiteGPF enabled (forced for BUILD_TESTING)")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Classic FlexAID target: P0 build-level optimization (LTO/native/NDEBUG + PGO)
# ─────────────────────────────────────────────────────────────────────────────
# Historically only the FlexAIDdS target carried the aggressive optimization
# flag block (LTO + INTERPROCEDURAL_OPTIMIZATION + -march=native / -mcpu=native +
# -DNDEBUG + strip); the classic FlexAID target shipped with only -O3 -ffast-math.
# BUILD_FLEXAID_FAST mirrors BUILD_FLEXAIDDS_FAST and lifts that same block onto
# the classic FlexAID target (applied via flexaids_apply_fast_optimization()).
# Default OFF per METHODOLOGY §1. The block includes -march=native, which makes
# the emitted binary a function of the build host's CPU: two correct builds of
# identical source on different runners produce different md5s, so §1 step 1
# ("record both md5s") and the recorded reference md5 at METHODOLOGY.md:131
# cannot be compared across machines. On top of the existing -ffast-math this
# also permits per-host FMA contraction and reassociation, so the numerical
# output moves too — not just the binary. Turn ON explicitly for local perf
# work. (Note: BUILD_FLEXAIDDS_FAST carries the same hazard for the FlexAIDdS
# target and predates this PR; that is a separate fix, not this one's to make.)
option(BUILD_FLEXAID_FAST "Apply the FlexAIDdS LTO + native-arch + NDEBUG optimization block to the classic FlexAID target" OFF)

# Profile-Guided Optimization for the classic FlexAID target. Values:
#   off      — no PGO instrumentation (default)
#   generate — instrumented build (-fprofile-generate); run a representative dock
#              (e.g. 1G9V, fixed FLEXAID_SEED) to emit profile data, then
#              reconfigure with -DFLEXAID_PGO=use and rebuild.
#   use      — optimized build consuming the profile
#              (-fprofile-use [-fprofile-correction on GCC])
set(FLEXAID_PGO "off" CACHE STRING "FlexAID Profile-Guided Optimization mode: off | generate | use")
set_property(CACHE FLEXAID_PGO PROPERTY STRINGS off generate use)

# LTO/native codegen is incompatible with gcov coverage and with ASan/TSan/UBSan
# instrumented builds. Reuse the SAME incompatibility guard the
# BUILD_FLEXAIDDS_FAST block uses (see root CMakeLists.txt) so an instrumented
# build silently drops the fast flags instead of failing to link.
if(FLEXAIDS_ENABLE_COVERAGE
   OR CMAKE_CXX_FLAGS MATCHES "-fsanitize"
   OR CMAKE_C_FLAGS MATCHES "-fsanitize"
   OR CMAKE_EXE_LINKER_FLAGS MATCHES "-fsanitize")
    if(BUILD_FLEXAID_FAST)
        message(STATUS "Disabling BUILD_FLEXAID_FAST (incompatible with coverage/sanitizer flags)")
    endif()
    set(BUILD_FLEXAID_FAST OFF CACHE BOOL "Apply the FlexAIDdS LTO + native-arch + NDEBUG optimization block to the classic FlexAID target" FORCE)
endif()

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
# Helper: lift the FlexAIDdS optimization flag block onto an arbitrary target
# ─────────────────────────────────────────────────────────────────────────────
# Mirrors the flag block applied to the FlexAIDdS target in root CMakeLists.txt
# (the BUILD_FLEXAIDDS_FAST branch): -flto + INTERPROCEDURAL_OPTIMIZATION +
# -march=native / -mcpu=native + -DNDEBUG + link-time strip (-s), or /GL /LTCG on
# MSVC. It ADDS to (does not replace) the target's existing -O3 -ffast-math flags.
#
# No-op unless BUILD_FLEXAID_FAST is ON (the coverage/sanitizer guard above
# force-disables that option so instrumented builds stay linkable).
#
# LTO/OBJECT-library caveat (see root CMakeLists.txt notes): flexaid_core is
# compiled WITHOUT IPO so that non-LTO consumers (tests) can still link its .o
# files. We enable IPO only on the leaf executable and add -flto to its own
# compile/link, exactly as the proven FlexAIDdS target does — mixed LTO/non-LTO
# object linking is handled by the compiler driver.
function(flexaids_apply_fast_optimization target_name)
    if(NOT BUILD_FLEXAID_FAST)
        return()
    endif()
    if(MSVC)
        target_compile_options(${target_name} PRIVATE /GL /DNDEBUG)
        target_link_options(${target_name} PRIVATE /LTCG)
    else()
        target_compile_options(${target_name} PRIVATE -flto -DNDEBUG
            $<$<AND:$<CXX_COMPILER_ID:AppleClang,Clang>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>>:-mcpu=native>
        )
        if(NOT APPLE)
            target_compile_options(${target_name} PRIVATE -march=native)
            target_link_options(${target_name} PRIVATE -flto -s)
        else()
            target_link_options(${target_name} PRIVATE -flto)
        endif()
    endif()

    include(CheckIPOSupported)
    check_ipo_supported(RESULT _flexaid_ipo_ok OUTPUT _flexaid_ipo_err)
    if(_flexaid_ipo_ok)
        set_property(TARGET ${target_name} PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(WARNING "IPO not supported for ${target_name}: ${_flexaid_ipo_err}")
    endif()
    message(STATUS "${target_name}: fast-optimization block enabled (LTO/IPO + native + NDEBUG)")
endfunction()

# ─────────────────────────────────────────────────────────────────────────────
# Helper: apply Profile-Guided Optimization flags (FLEXAID_PGO) to a target
# ─────────────────────────────────────────────────────────────────────────────
# generate → -fprofile-generate ; use → -fprofile-use (+ -fprofile-correction on
# GCC). -fprofile-correction is GCC-only (tolerates inconsistent counters from a
# multi-threaded instrumented run); it is gated to GNU so Clang does not error.
function(flexaids_apply_pgo target_name)
    if(FLEXAID_PGO STREQUAL "off")
        return()
    endif()
    if(MSVC)
        message(WARNING "FLEXAID_PGO=${FLEXAID_PGO} is not wired for MSVC (use /GENPROFILE|/USEPROFILE manually); ignoring")
        return()
    endif()
    if(FLEXAID_PGO STREQUAL "generate")
        target_compile_options(${target_name} PRIVATE -fprofile-generate)
        target_link_options(${target_name} PRIVATE -fprofile-generate)
        message(STATUS "${target_name}: PGO instrumentation build (-fprofile-generate) — "
                       "run a representative dock, then reconfigure with -DFLEXAID_PGO=use")
    elseif(FLEXAID_PGO STREQUAL "use")
        target_compile_options(${target_name} PRIVATE
            -fprofile-use $<$<CXX_COMPILER_ID:GNU>:-fprofile-correction>)
        target_link_options(${target_name} PRIVATE
            -fprofile-use $<$<CXX_COMPILER_ID:GNU>:-fprofile-correction>)
        message(STATUS "${target_name}: PGO optimized build (-fprofile-use)")
    else()
        message(FATAL_ERROR "FLEXAID_PGO must be one of: off | generate | use (got '${FLEXAID_PGO}')")
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
