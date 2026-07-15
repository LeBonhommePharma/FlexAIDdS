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

# ─── Eigen3 (header-only; vendored at LIB/vendor/eigen) ─────────────────────
# Discovery priority:
#   1. Vendored submodule  — always works, no package manager needed
#   2. System find_package — uses installed libeigen3-dev / brew eigen / etc.
#   3. pkg-config          — additional system fallback
#   4. FetchContent        — downloads from GitLab if nothing else is available
#
# For a zero-dependency build on any platform just run:
#   git clone --recurse-submodules <repo>
# or, if already cloned:
#   git submodule update --init --recursive

set(_eigen_vendor_dir "${CMAKE_SOURCE_DIR}/LIB/vendor/eigen")
if(EXISTS "${_eigen_vendor_dir}/Eigen/Dense")
    message(STATUS "Eigen3: using vendored copy (LIB/vendor/eigen) — no system package needed")
    if(NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
    endif()
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_eigen_vendor_dir}")
    set(Eigen3_FOUND TRUE)
else()
    # Vendored submodule not initialised — try system / network fallbacks
    find_package(Eigen3 3.4 QUIET NO_MODULE)
    if(Eigen3_FOUND)
        message(STATUS "Eigen3 ${Eigen3_VERSION} found via system cmake config")
    else()
        find_package(PkgConfig QUIET)
        if(PKG_CONFIG_FOUND)
            pkg_check_modules(EIGEN3 QUIET eigen3)
        endif()
        if(EIGEN3_FOUND)
            message(STATUS "Eigen3 found via pkg-config")
            if(NOT TARGET Eigen3::Eigen)
                add_library(Eigen3::Eigen INTERFACE IMPORTED)
            endif()
            set_target_properties(Eigen3::Eigen PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIRS}")
            set(Eigen3_FOUND TRUE)
        else()
            # Last resort: download tarball at configure time
            message(STATUS "Eigen3 not found locally — fetching via FetchContent")
            include(FetchContent)
            FetchContent_Declare(eigen3
                URL      https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
                URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
                DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
            FetchContent_Populate(eigen3)
            if(NOT TARGET Eigen3::Eigen)
                add_library(Eigen3::Eigen INTERFACE IMPORTED)
            endif()
            set_target_properties(Eigen3::Eigen PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${eigen3_SOURCE_DIR}")
            set(Eigen3_FOUND TRUE)
            message(STATUS "Eigen3 3.4.0 fetched via FetchContent (headers only)")
        endif()
    endif()
endif()

if(NOT Eigen3_FOUND)
    message(FATAL_ERROR
        "Eigen3 not found. The easiest fix is to initialise the vendored submodule:\n"
        "  git submodule update --init --recursive\n"
        "Alternatively install a system package:\n"
        "  Linux:   sudo apt install libeigen3-dev\n"
        "  macOS:   brew install eigen\n"
        "  Windows: choco install eigen")
endif()

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
