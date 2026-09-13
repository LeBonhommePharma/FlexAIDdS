# cmake/FlexAIDDependencies.cmake
#
# FlexAIDdS — External dependency discovery and setup.
#
# This module centralizes:
#   - OpenMP
#   - Eigen3 (vendored / system / pkg-config / FetchContent)
#   - pybind11 (for BUILD_PYTHON_BINDINGS)
#   - CUDA
#   - ROCm / HIP
#   - Metal (macOS)
#   - MPI (must be before LIB/ subdirectory)
#
# Included early, after FlexAIDOptions.cmake and before add_subdirectory(LIB).
#
# Benefits:
#   - Root CMakeLists.txt is dramatically smaller and easier to read
#   - All "find this thing" logic lives in one place
#   - Prepares for proper package config / export in the future
#
# See audit feedback and AGENTS.md for the monolith reduction goal.

cmake_minimum_required(VERSION 3.28)

# ─── OpenMP ───────────────────────────────────────────────────────────────
if(FLEXAIDS_USE_OPENMP)
    if(APPLE)
        set(_flexaids_libomp_prefixes)
        if(DEFINED ENV{HOMEBREW_PREFIX})
            list(APPEND _flexaids_libomp_prefixes "$ENV{HOMEBREW_PREFIX}/opt/libomp")
        endif()
        list(APPEND _flexaids_libomp_prefixes
            "/opt/homebrew/opt/libomp"
            "/usr/local/opt/libomp")

        foreach(_flexaids_libomp_prefix IN LISTS _flexaids_libomp_prefixes)
            if(EXISTS "${_flexaids_libomp_prefix}/include/omp.h" AND
               EXISTS "${_flexaids_libomp_prefix}/lib/libomp.dylib")
                set(OpenMP_CXX_FLAGS
                    "-Xpreprocessor -fopenmp -I${_flexaids_libomp_prefix}/include"
                    CACHE STRING "CXX compiler flags for OpenMP parallelization" FORCE)
                set(OpenMP_CXX_LIB_NAMES
                    "omp"
                    CACHE STRING "CXX compiler libraries for OpenMP parallelization" FORCE)
                set(OpenMP_omp_LIBRARY
                    "${_flexaids_libomp_prefix}/lib/libomp.dylib"
                    CACHE FILEPATH "OpenMP libomp library" FORCE)
                break()
            endif()
        endforeach()
    endif()

    find_package(OpenMP)
    if(NOT OpenMP_CXX_FOUND)
        message(WARNING "FLEXAIDS_USE_OPENMP=ON but OpenMP was not found; CPU fallback will be serial")
    endif()
endif()

# ─── Eigen3 (header-only) ───────────────────────────────────────────────────
#
# Version, hash, resolution order and the divergence gate all live in one
# place.  Nothing about Eigen is decided here — see cmake/FlexAIDEigen.cmake
# for the declared version and the reasoning behind it.
#
# For a build that needs no package manager and no network, initialise the
# vendored submodule:
#   git submodule update --init --recursive LIB/vendor/eigen
# Otherwise the pinned release tarball is downloaded and hash-verified
# automatically, so no manual step is required on any platform.

include(${CMAKE_CURRENT_LIST_DIR}/FlexAIDEigen.cmake)
flexaids_resolve_eigen()

# ─── pybind11 (for Python bindings) ────────────────────────────────────────
if(BUILD_PYTHON_BINDINGS)
    # pybind11's LTO probe (try_compile) inherits CMAKE_CXX_STANDARD from the
    # normal-variable scope.  On MSVC, CMAKE_CXX_STANDARD is already capped to
    # 20 above, but the probe can still fail if try_compile picks up a stale
    # cached value.  Temporarily unset the normal variable so the probe runs
    # against the compiler's default standard, then restore to 20.
    if(MSVC)
        set(_pybind11_saved_cxx_std "${CMAKE_CXX_STANDARD}")  # "20" after MSVC cap
        unset(CMAKE_CXX_STANDARD)   # remove normal variable only; line 35 set no cache entry
    endif()
    find_package(pybind11 CONFIG QUIET)
    if(MSVC AND DEFINED _pybind11_saved_cxx_std)
        set(CMAKE_CXX_STANDARD "${_pybind11_saved_cxx_std}")   # restore as normal variable
        set(CMAKE_CXX_STANDARD_REQUIRED ON)
    endif()
    if(NOT pybind11_FOUND)
        message(WARNING "pybind11 not found. Install with: pip install pybind11[global]")
        set(BUILD_PYTHON_BINDINGS OFF)
    else()
        message(STATUS "pybind11 ${pybind11_VERSION} found — Python bindings enabled")
    endif()
endif()

# ─── CUDA (optional) ─────────────────────────────────────────────────────
if(FLEXAIDS_USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# ─── ROCm / HIP (optional) ───────────────────────────────────────────────
if(FLEXAIDS_USE_ROCM)
    if(CMAKE_VERSION VERSION_LESS "3.21")
        message(FATAL_ERROR "ROCm/HIP support requires CMake >= 3.21")
    endif()

    if(NOT DEFINED ENV{ROCM_PATH})
        set(_rocm_hints /opt/rocm /usr/local/rocm)
    else()
        set(_rocm_hints $ENV{ROCM_PATH})
    endif()
    list(APPEND CMAKE_PREFIX_PATH ${_rocm_hints})

    enable_language(HIP)
    find_package(hip REQUIRED
        HINTS ${_rocm_hints}
        PATH_SUFFIXES lib/cmake/hip hip/lib/cmake/hip)

    set(FLEXAIDS_HIP_ARCHITECTURES "gfx908;gfx90a;gfx942"
        CACHE STRING "Semicolon-separated AMD GPU architectures for HIP")
endif()

# ─── Metal (optional, macOS only) ────────────────────────────────────────
if(FLEXAIDS_USE_METAL)
    if(NOT APPLE)
        message(FATAL_ERROR "Metal acceleration is only supported on macOS")
    endif()
    
    if(NOT _HAS_OBJCXX)
        message(WARNING "Objective-C++ not available — disabling Metal support")
        set(FLEXAIDS_USE_METAL OFF)
    else()
        find_library(METAL_LIBRARY Metal REQUIRED)
        find_library(FOUNDATION_LIBRARY Foundation REQUIRED)
        find_library(METALKIT_LIBRARY MetalKit REQUIRED)
        message(STATUS "Metal acceleration enabled")
    endif()
else()
    if(APPLE)
        message(STATUS "Metal acceleration available but disabled (use -DFLEXAIDS_USE_METAL=ON)")
    endif()
endif()

# ─── MPI detection (must precede add_subdirectory(LIB) so MPI::MPI_CXX
#     is defined when LIB/CMakeLists.txt calls target_link_libraries) ────
if(FLEXAIDS_USE_MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI enabled — distributed parallel docking active")
endif()

message(STATUS "FlexAIDDependencies loaded — external packages resolved")
