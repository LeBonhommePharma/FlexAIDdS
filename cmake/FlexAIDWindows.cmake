# cmake/FlexAIDWindows.cmake
#
# MSVC / Windows toolchain defaults and acceleration validation for FlexAIDdS.
# Included from the root CMakeLists.txt after compiler detection.

cmake_minimum_required(VERSION 3.28)

if(NOT MSVC)
    return()
endif()

# Dynamic CRT (/MD) — matches Python extensions, pybind11, and GoogleTest when
# gtest_force_shared_crt=ON.  Per-target /MT flags are intentionally avoided.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL" CACHE STRING "" FORCE)

# MSVC ships OpenMP 2.0 only; default OFF unless the user explicitly enables it.
if(FLEXAIDS_USE_OPENMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "MSVC OpenMP: enabled (legacy v2.0 runtime — no omp simd)")
    else()
        message(WARNING
            "FLEXAIDS_USE_OPENMP=ON on MSVC but OpenMP was not found. "
            "Install LLVM OpenMP via vcpkg or the VS C++ workload.")
    endif()
else()
    message(STATUS "MSVC: OpenMP disabled (legacy v2.0 runtime). "
                   "Pass -DFLEXAIDS_USE_OPENMP=ON to override.")
endif()

# Global definitions used across engine, tools, tests, and Python bindings.
add_compile_definitions(
    _CRT_SECURE_NO_WARNINGS
    _USE_MATH_DEFINES
    NOMINMAX
    WIN32_LEAN_AND_MEAN
)

# Prefer the newest standard MSVC exposes below C++26.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("/std:c++23" FLEXAIDS_MSVC_SUPPORTS_CXX23)
if(FLEXAIDS_MSVC_SUPPORTS_CXX23)
    add_compile_options(/std:c++23)
    set(CMAKE_CXX_STANDARD 23)
    message(STATUS "MSVC ${MSVC_VERSION}: using /std:c++23 (C++26 not yet available)")
else()
    add_compile_options(/std:c++20)
    message(STATUS "MSVC ${MSVC_VERSION}: using /std:c++20 (C++26 not yet available)")
endif()

# ── Configure-time validation of Windows SIMD compiler flags ───────────────
set(FLEXAIDS_MSVC_SUPPORTS_AVX2 FALSE)
set(FLEXAIDS_MSVC_SUPPORTS_AVX512 FALSE)
set(FLEXAIDS_WINDOWS_SIMD_ENABLED "scalar")
if(FLEXAIDS_IS_X86)
    if(FLEXAIDS_USE_AVX512)
        check_cxx_compiler_flag("/arch:AVX512" FLEXAIDS_MSVC_SUPPORTS_AVX512)
        if(FLEXAIDS_MSVC_SUPPORTS_AVX512)
            set(FLEXAIDS_WINDOWS_SIMD_ENABLED "AVX-512 (/arch:AVX512)")
        else()
            message(WARNING "MSVC does not accept /arch:AVX512 on this host — falling back")
            set(FLEXAIDS_USE_AVX512 OFF CACHE BOOL "Enable AVX-512 SIMD acceleration" FORCE)
        endif()
    endif()
    if(FLEXAIDS_USE_AVX2 AND FLEXAIDS_WINDOWS_SIMD_ENABLED STREQUAL "scalar")
        check_cxx_compiler_flag("/arch:AVX2" FLEXAIDS_MSVC_SUPPORTS_AVX2)
        if(FLEXAIDS_MSVC_SUPPORTS_AVX2)
            set(FLEXAIDS_WINDOWS_SIMD_ENABLED "AVX2 (/arch:AVX2)")
        else()
            message(WARNING "MSVC does not accept /arch:AVX2 on this host — scalar fallback")
            set(FLEXAIDS_USE_AVX2 OFF CACHE BOOL "Enable AVX2 SIMD acceleration" FORCE)
        endif()
    endif()
elseif(FLEXAIDS_IS_ARM)
    set(FLEXAIDS_WINDOWS_SIMD_ENABLED "ARM64 NEON (baseline ISA)")
endif()

# CUDA on Windows requires the CUDA Toolkit + MSVC host compiler pairing.
if(FLEXAIDS_USE_CUDA)
    message(STATUS "MSVC + CUDA: host compiler must match VS 2022 toolset used by NVCC")
endif()

# Central helper — common warning level without duplicating CRT or arch flags.
function(flexaids_apply_msvc_target_options target_name)
    if(NOT MSVC)
        return()
    endif()
    get_target_property(_existing_opts ${target_name} COMPILE_OPTIONS)
    if(_existing_opts)
        set(_has_w4 FALSE)
        foreach(_opt IN LISTS _existing_opts)
            if(_opt STREQUAL "/W4" OR _opt STREQUAL "/w4")
                set(_has_w4 TRUE)
            endif()
        endforeach()
        if(NOT _has_w4)
            target_compile_options(${target_name} PRIVATE /W4)
        endif()
    else()
        target_compile_options(${target_name} PRIVATE /W4)
    endif()
endfunction()

# Hot-path helper — fast-math style optimizations safe for docking scoring.
function(flexaids_configure_msvc_hotpath target_name)
    if(NOT MSVC)
        return()
    endif()
    target_compile_options(${target_name} PRIVATE /fp:fast)
endfunction()

# LTO helper for shipping executables (FlexAIDdS, tENCoM).
function(flexaids_configure_msvc_lto target_name)
    if(NOT MSVC)
        return()
    endif()
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err)
    if(_ipo_ok)
        set_property(TARGET ${target_name} PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
        target_compile_options(${target_name} PRIVATE /GL)
        target_link_options(${target_name} PRIVATE /LTCG)
    else()
        message(WARNING "${target_name}: MSVC LTO unavailable — ${_ipo_err}")
    endif()
endfunction()

message(STATUS "── Windows acceleration summary ──")
message(STATUS "  CPU arch     : ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "  SIMD compile : ${FLEXAIDS_WINDOWS_SIMD_ENABLED}")
message(STATUS "  OpenMP       : ${FLEXAIDS_USE_OPENMP}")
message(STATUS "  CUDA         : ${FLEXAIDS_USE_CUDA}")
message(STATUS "  SoA distances: ${FLEXAIDS_USE_SOA_DISTANCES}")
message(STATUS "  Benchmarks   : OFF (UNIX only — use WSL2 for DatasetRunner campaigns)")
message(STATUS "  CRT          : /MD (MultiThreadedDLL)")
message(STATUS "FlexAIDWindows loaded — MSVC helpers active")