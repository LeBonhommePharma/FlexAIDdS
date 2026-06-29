# Build: cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build

cmake_minimum_required(VERSION 3.28)

# ─── Graceful OBJCXX Fallback ──────────────────────────────────────────────
# Detects OBJCXX compiler availability before requiring it.
# On non-Apple or incomplete toolchain: builds with CXX only, Metal disabled.

set(_ENABLED_LANGUAGES "CXX")

if(APPLE)
    execute_process(
        COMMAND xcrun -f clang++
        OUTPUT_QUIET ERROR_QUIET
        RESULT_VARIABLE _XCRUN_RESULT
    )
    if(_XCRUN_RESULT EQUAL 0)
        list(APPEND _ENABLED_LANGUAGES "OBJCXX")
        set(_HAS_OBJCXX TRUE)
        message(STATUS "Objective-C++ compiler found — Metal support enabled")
    else()
        set(_HAS_OBJCXX FALSE)
        message(WARNING "Objective-C++ compiler not available — Metal support disabled")
        set(FLEXAIDS_USE_METAL OFF CACHE BOOL "Enable Metal GPU acceleration" FORCE)
    endif()
else()
    set(_HAS_OBJCXX FALSE)
endif()

project(FlexAID VERSION 2.0.0 LANGUAGES ${_ENABLED_LANGUAGES})

# ── Embed git provenance in binary ────────────────────────────────────────
execute_process(
    COMMAND git -C "${CMAKE_SOURCE_DIR}" status --porcelain
    OUTPUT_VARIABLE _GIT_STATUS OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
if(_GIT_STATUS STREQUAL "")
    set(FLEXAIDS_GIT_DIRTY "0")
else()
    set(FLEXAIDS_GIT_DIRTY "1")
    message(WARNING "Building from a DIRTY working tree — reproducibility metadata will be flagged.")
endif()

execute_process(
    COMMAND git -C "${CMAKE_SOURCE_DIR}" rev-parse --short HEAD
    OUTPUT_VARIABLE FLEXAIDS_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)

add_compile_definitions(
    FLEXAIDS_GIT_DIRTY=${FLEXAIDS_GIT_DIRTY}
    FLEXAIDS_GIT_COMMIT="${FLEXAIDS_GIT_COMMIT}")

include(CheckCXXCompilerFlag)

# ─── C++ Standard ─────────────────────────────────────────────────────────
# Default to C++26 for full features, but on AppleClang we default to C++23
# for better SDK compatibility (Apple's libc++ in current/future macOS versions
# can have incomplete C++26 support in some math/type_traits headers when
# certain includes are pulled in). Users can force C++26 with
# -DFLEXAIDS_FORCE_CXX26=ON.
set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if(NOT FLEXAIDS_FORCE_CXX26)
        set(CMAKE_CXX_STANDARD 23)
        message(STATUS "AppleClang on macOS: defaulting to C++23 for SDK compatibility. "
                       "Use -DFLEXAIDS_FORCE_CXX26=ON to request C++26 (may hit libc++ template issues).")
    endif()
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0")
        message(FATAL_ERROR
            "FlexAIDdS requires AppleClang ≥ 16 (Xcode 16) for C++23/26.")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "14.0")
        message(FATAL_ERROR
            "FlexAIDdS requires C++26 — found GCC "
            "${CMAKE_CXX_COMPILER_VERSION}. Upgrade to GCC ≥ 14.")
    endif()
elseif(MSVC)
    if(MSVC_VERSION LESS 1940)
        message(FATAL_ERROR
            "FlexAIDdS requires MSVC ≥ 19.40 (Visual Studio 2022 17.10+) for C++26 — "
            "found MSVC ${MSVC_VERSION}.")
    endif()
    # MSVC does not yet implement /std:c++26 even in VS 2022 17.14 (19.44).
    # Cap to C++20 so the Python-bindings smoke test can configure and build on
    # Windows CI.  The full FlexAID engine is intentionally excluded from the
    # Windows CI matrix (see ci.yml) due to deep C++26 gaps in MSVC.
    message(WARNING
        "MSVC ${MSVC_VERSION}: /std:c++26 is not yet available — "
        "building with C++20. Only the _core Python extension will build "
        "correctly; the full FlexAID engine requires GCC ≥ 14 or Clang ≥ 18.")
    set(CMAKE_CXX_STANDARD 20)   # override normal variable — CACHE form leaves it shadowed at 26
endif()

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

# ─── Code coverage configuration ──────────────────────────────────────────────
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

# ─── options ────────────────────────────────────────────────────────────────
option(FLEXAIDS_USE_CUDA    "Enable CUDA GPU evaluation"          OFF)
option(FLEXAIDS_USE_ROCM    "Enable ROCm/HIP GPU evaluation"      OFF)
if(APPLE AND _HAS_OBJCXX)
    option(FLEXAIDS_USE_METAL   "Enable Metal GPU acceleration"       ON)
else()
    option(FLEXAIDS_USE_METAL   "Enable Metal GPU acceleration"       OFF)
endif()
option(FLEXAIDS_USE_AVX2    "Enable AVX2 SIMD acceleration"       ON)
option(FLEXAIDS_USE_AVX512  "Enable AVX-512 SIMD acceleration"    OFF)
option(FLEXAIDS_USE_SOA_DISTANCES "Route Voronoi hot-path distances through AtomSoA float SoA arrays (C1)" ON)
option(FLEXAIDS_USE_OPENMP  "Enable OpenMP thread parallelism"    ON)
# Eigen3 is a hard requirement — no option to disable
set(FLEXAIDS_USE_EIGEN ON)
option(FLEXAIDS_USE_256_MATRIX "Enable 256×256 soft contact matrix" ON)
option(FLEXAIDS_USE_MPI     "Enable MPI distributed parallel docking" OFF)
option(BUILD_PYTHON_BINDINGS "Build Python bindings via pybind11" OFF)
option(FLEXAIDS_BUILD_CORE  "Build core executables (FlexAID/FlexAIDdS/tENCoM)" ON)
option(BUILD_FLEXAIDDS_FAST "Build ultra-fast FlexAIDdS docking executable (LTO + native)" ON)
option(FLEXAIDS_ENABLE_COVERAGE "Enable code coverage tracking (requires gcov/lcov)" OFF)

# ─── platform / SIMD helpers ───────────────────────────────────────────────
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

# ─── OpenMP ───────────────────────────────────────────────────────────────
if(FLEXAIDS_USE_OPENMP)
    find_package(OpenMP)
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

# ─── Core library (LIB/) ────────────────────────────────────────────────
add_subdirectory(LIB)

# ─── executable ─────────────────────────────────────────────────────────
if(FLEXAIDS_BUILD_CORE)
add_executable(FlexAID "")
target_sources(FlexAID PRIVATE
        LIB/top.cpp
        LIB/native_score.cpp
        LIB/fileio.cpp
        LIB/assign_shift.cpp
        LIB/buildic.cpp
        LIB/Vcontacts.cpp
        LIB/ic2cf.cpp
        LIB/vcfunction.cpp
        LIB/create_rebuild_list.cpp
        LIB/read_coor.cpp
        LIB/read_constraints.cpp
        LIB/assign_constraint.cpp
        LIB/modify_pdb.cpp
        LIB/update_constraint.cpp
        LIB/update_optres.cpp
        LIB/update_bonded.cpp
        LIB/read_grid.cpp
        LIB/read_input.cpp
        LIB/read_lig.cpp
        LIB/read_pdb.cpp
        LIB/add2_optimiz_vec.cpp
        LIB/assign_radius.cpp
        LIB/assign_radii.cpp
        LIB/assign_radii_types.cpp
        LIB/assign_types.cpp
        LIB/bondedlist.cpp
        LIB/buildcc.cpp
        LIB/buildlist.cpp
        LIB/geometry.cpp
        LIB/read_conect.cpp
        LIB/residue_conect.cpp
        LIB/wif083.cpp
        LIB/calc_cleftic.cpp
        LIB/write_pdb.cpp
        LIB/gaboom.cpp
        LIB/coarse_init.cpp
        LIB/VibEntropy.cpp
        LIB/GAContext.cpp
        LIB/calc_center.cpp
        LIB/ic_bounds.cpp
        LIB/buildic_point.cpp
        LIB/buildcc_point.cpp
        LIB/read_rmsdst.cpp
        LIB/calc_rmsd_chrom.cpp
        LIB/calc_rmsd.cpp
        LIB/write_rrd.cpp
        LIB/write_sphere.cpp
        LIB/write_grid.cpp
        LIB/write_rrg.cpp
        LIB/read_rotlib.cpp
        LIB/number_of_dihedrals.cpp
        LIB/read_flexscfile.cpp
        LIB/read_rotobs.cpp
        LIB/set_intprob.cpp
        LIB/build_rotamers.cpp
        LIB/receptor_prep.cpp
        LIB/check_clash.cpp
        LIB/read_normalgrid.cpp
        LIB/read_eigen.cpp
        LIB/assign_eigen.cpp
        LIB/print_surfmat.cpp
        LIB/alter_mode.cpp
        LIB/dee_pivot.cpp
        LIB/read_emat.cpp
        LIB/read_spheres.cpp
        LIB/cluster.cpp
        LIB/DensityPeak_Cluster.cpp
        LIB/FastOPTICS_cluster.cpp
        LIB/MinibatchSampler.cpp
        LIB/FOPTICS.cpp
        LIB/fast_optics.cpp
        LIB/BindingMode.cpp
        LIB/rna_structure.cpp
        LIB/shortest_path.cpp
        LIB/assign_shortflex.cpp
        LIB/slice_grid.cpp
        LIB/partition_grid.cpp
        LIB/generate_grid.cpp
        LIB/maps.cpp
        # ── new C++20 modules ──
        LIB/tENCoM/tencm.cpp
        LIB/statmech.cpp
        LIB/encom.cpp             # ← NEW: ENCoM vibrational entropy
        LIB/DiFT/DiFT.cpp
        # ── GIST water displacement + directional H-bond scoring ──
        LIB/GISTEvaluator.cpp
        LIB/HBondEvaluator.cpp
        # ── TurboQuant vector quantization (Zandieh et al. 2025) ──
        LIB/TurboQuant.cpp
        # ── ProcessLigand (Phase 3: standalone ligand preprocessing) ──
        LIB/ProcessLigand/ProcessLigand.cpp
        LIB/ProcessLigand/SmilesParser.cpp
        LIB/ProcessLigand/RingPerception.cpp
        LIB/ProcessLigand/Aromaticity.cpp
        LIB/ProcessLigand/RotatableBonds.cpp
        LIB/ProcessLigand/ValenceChecker.cpp
        LIB/ProcessLigand/SybylTyper.cpp
        LIB/ProcessLigand/FlexAIDWriter.cpp
        LIB/ProcessLigand/ROCmDispatch.cpp
        LIB/ProcessLigand/CoordBuilder.cpp
        LIB/CifReader.cpp
        LIB/MultiModelDock.cpp
        LIB/LibrarySplitter.cpp
        LIB/ReferenceEntropy.cpp
        LIB/ParallelCampaign.cpp
        # ── v1.5 thermodynamic + hardware acceleration modules ──
        LIB/hardware_detect.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        # ── Shannon Thermodynamic GA (entropy-collapse driver) ──
        LIB/shannon_ga.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
        LIB/LigandRingFlex/LigandRingFlex.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/ChiralCenter/ChiralCenterGene.cpp
        # ── CleftDetector + Mol2Reader + SdfReader ──
        LIB/CleftDetector.cpp
        LIB/Mol2Reader.cpp
        LIB/SdfReader.cpp
        # ── CavityDetect: native SURFNET cavity detection ──
        LIB/CavityDetect/CavityDetect.cpp
        # ── PTMAttachment: PTM/glycan modification on target chain ──
        LIB/PTMAttachment/PTMAttachment.cpp
        # ── JSON config parser (zero external dependencies) ──
        LIB/config_parser.cpp
        # ── ThermodynamicEngine: G_bind decomposition (disabled by default) ──
        LIB/ThermodynamicEngine.cpp
        # ── Unified runtime hardware dispatch ──
        LIB/UnifiedHardwareDispatch.cpp
        # ── Parallel grid-decomposed docking ──
        LIB/GridDecomposer.cpp
        LIB/SharedPosePool.cpp
        LIB/ParallelDock.cpp
        # ── CoarseScreen: NRGRank coarse-grained screening ──
        LIB/CoarseScreen.cpp
        LIB/TwoStageScreen.cpp
        # ── TargetServer: target-centric docking architecture ──
        LIB/TargetValidation.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/MultiSiteGPF.cpp
        LIB/TargetKnowledgeBase.cpp
        LIB/TargetServer.cpp
        # ── InStreamClustering: online medoid clustering during GA ──
        LIB/InStreamClustering.cpp
)

target_include_directories(FlexAID PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ChiralCenter
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PTMAttachment
)

# ─── compiler flags ────────────────────────────────────────────────────
if(MSVC)
    target_compile_options(FlexAID PRIVATE /W4 /O2 /MT)
    target_compile_definitions(FlexAID PRIVATE
        _CRT_SECURE_NO_WARNINGS      # suppress MSVC sprintf/strcpy deprecation
        _USE_MATH_DEFINES             # expose M_PI in <cmath>
        NOMINMAX                      # prevent windows.h min/max macros
    )
else()
    target_compile_options(FlexAID PRIVATE -Wall -O3 -ffast-math -fno-finite-math-only)
endif()
flexaids_configure_simd(FlexAID)

# Stage runtime data files next to the FlexAID binary (see note on the
# FlexAIDdS target): the engine resolves MC_st0r5.2_6.dat / AMINO.def /
# NUCLEOTIDES.def relative to the executable's own directory.
add_custom_command(TARGET FlexAID POST_BUILD
    COMMENT "Staging FlexAID runtime data files next to the binary"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/MC_st0r5.2_6.dat"
        "${CMAKE_SOURCE_DIR}/AMINO.def"
        "${CMAKE_SOURCE_DIR}/NUCLEOTIDES.def"
        "$<TARGET_FILE_DIR:FlexAID>")

# ─── Eigen3 link + definition ─────────────────────────────────────────────
if(FLEXAIDS_USE_EIGEN AND Eigen3_FOUND)
    target_link_libraries(FlexAID PRIVATE Eigen3::Eigen)
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_HAS_EIGEN)
elseif(FLEXAIDS_USE_EIGEN AND EIGEN3_FOUND)
    target_include_directories(FlexAID PRIVATE ${EIGEN3_INCLUDE_DIRS})
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_HAS_EIGEN)
endif()

# ─── OpenMP link ────────────────────────────────────────────────────────
if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(FlexAID PRIVATE OpenMP::OpenMP_CXX)
endif()

# ─── CUDA sources ───────────────────────────────────────────────────────
if(FLEXAIDS_USE_CUDA)
    target_sources(FlexAID PRIVATE
        LIB/cuda_eval.cu                          # existing CF batch kernel
        LIB/ShannonThermoStack/shannon_cuda.cu    # Shannon entropy histogram kernel
LIB/tENCoM/tencm_cuda.cu                  # TENCoM contact discovery + Hessian assembly
        LIB/gpu_fast_optics.cu                    # GPU-accelerated FastOPTICS neighbour search
        LIB/TurboQuant.cu                          # TurboQuant batch quantize/dequantize kernels
        # P0 validator: new GPU batch distance kernels
        LIB/gpu_batch_distance.cu
        LIB/gpu_batch_distance.cpp
    )
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_USE_CUDA)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "70;75;80;86;89;90")
        if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.8")
            list(APPEND CMAKE_CUDA_ARCHITECTURES 100 120)
            message(STATUS "CUDA >= 12.8 — Blackwell sm_100 + sm_120 targets enabled")
        elseif(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "12.6")
            list(APPEND CMAKE_CUDA_ARCHITECTURES 100)
            message(STATUS "CUDA >= 12.6 — Blackwell sm_100 target enabled (sm_120 requires CUDA >= 12.8)")
        else()
            message(WARNING "CUDA < 12.6 — Blackwell targets (sm_100, sm_120) disabled")
        endif()
    endif()
    # Host code is C++26 (see CMAKE_CXX_STANDARD above); NVCC 12.x does not accept
    # -std=c++26, so device code stays pinned at C++17 until the CUDA toolkit catches up.
    set_target_properties(FlexAID PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_STANDARD 17)
    target_link_libraries(FlexAID PRIVATE CUDA::cudart CUDA::cudadevrt)
    message(STATUS "CUDA enabled — CF batch + Shannon histogram + TENCoM Hessian kernels")
endif()

if(FLEXAIDS_USE_ROCM)
    target_sources(FlexAID PRIVATE LIB/hip_eval.hip)
    # P0 validator additions
    target_sources(FlexAID PRIVATE LIB/gpu_batch_distance_hip.hip LIB/gpu_fast_optics_hip.hip)
    set_source_files_properties(LIB/hip_eval.hip PROPERTIES LANGUAGE HIP)
    set_property(TARGET FlexAID PROPERTY HIP_ARCHITECTURES ${FLEXAIDS_HIP_ARCHITECTURES})
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_USE_ROCM)
    target_include_directories(FlexAID PRIVATE ${hip_INCLUDE_DIRS})
    target_link_libraries(FlexAID PRIVATE hip::host hip::device)
    message(STATUS "ROCm/HIP enabled — CF batch + runtime detection")
endif()

# ─── Metal sources ──────────────────────────────────────────────────────
if(FLEXAIDS_USE_METAL)
    target_sources(FlexAID PRIVATE LIB/metal_eval.mm)
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_USE_METAL)
    target_compile_options(FlexAID PRIVATE -fno-objc-arc)
    target_link_libraries(FlexAID PRIVATE
        ${METAL_LIBRARY}
        ${FOUNDATION_LIBRARY}
        ${METALKIT_LIBRARY}
    )
    set_source_files_properties(LIB/metal_eval.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # Shannon Metal GPU bridge (Apple Silicon Shannon entropy histogram)
    target_sources(FlexAID PRIVATE LIB/ShannonThermoStack/ShannonMetalBridge.mm)
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_HAS_METAL_SHANNON)
    set_source_files_properties(LIB/ShannonThermoStack/ShannonMetalBridge.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # TENCoM Metal bridge (GPU contact discovery + Hessian assembly)
    target_sources(FlexAID PRIVATE LIB/tENCoM/tencm_metal.mm)
    target_compile_definitions(FlexAID PRIVATE FLEXAIDS_HAS_METAL_TENCM)
    set_source_files_properties(LIB/tENCoM/tencm_metal.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # FastOPTICS Metal kNN (GPU-accelerated neighbor discovery)
    target_sources(FlexAID PRIVATE LIB/gpu_fast_optics_metal.mm)
    set_source_files_properties(LIB/gpu_fast_optics_metal.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # CavityDetect Metal bridge (Task 1.4: GPU probe-sphere generation)
    target_sources(FlexAID PRIVATE LIB/CavityDetect/CavityDetectMetalBridge.mm)
    set_source_files_properties(LIB/CavityDetect/CavityDetectMetalBridge.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # Compile CavityDetect.metal shader → .metallib
    set(CAVITY_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect/CavityDetect.metal")
    set(CAVITY_AIR        "${CMAKE_CURRENT_BINARY_DIR}/CavityDetect.air")
    set(CAVITY_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/CavityDetect.metallib")
    add_custom_command(
        OUTPUT  "${CAVITY_AIR}"
        COMMAND xcrun -sdk macosx metal -c "${CAVITY_METAL_SRC}" -o "${CAVITY_AIR}"
        DEPENDS "${CAVITY_METAL_SRC}"
        COMMENT "Compiling CavityDetect Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${CAVITY_METALLIB}"
        COMMAND xcrun -sdk macosx metallib "${CAVITY_AIR}" -o "${CAVITY_METALLIB}"
        DEPENDS "${CAVITY_AIR}"
        COMMENT "Linking CavityDetect.metallib"
        VERBATIM
    )
    add_custom_target(CavityDetectMetal DEPENDS "${CAVITY_METALLIB}")
    add_dependencies(FlexAID CavityDetectMetal)
    target_compile_definitions(FlexAID PRIVATE
        CAVITY_METALLIB_PATH="${CAVITY_METALLIB}"
    )
    # TurboQuant Metal bridge (TurboQuant.metal batch quantize/dequantize)
    target_sources(FlexAID PRIVATE LIB/TurboQuantMetalBridge.mm)
    set_source_files_properties(LIB/TurboQuantMetalBridge.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # Compile TurboQuant.metal shader → .metallib
    set(TQ_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/TurboQuant.metal")
    set(TQ_AIR        "${CMAKE_CURRENT_BINARY_DIR}/TurboQuant.air")
    set(TQ_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/TurboQuant.metallib")
    add_custom_command(
        OUTPUT  "${TQ_AIR}"
        COMMAND xcrun -sdk macosx metal -c "${TQ_METAL_SRC}" -o "${TQ_AIR}"
        DEPENDS "${TQ_METAL_SRC}"
        COMMENT "Compiling TurboQuant Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${TQ_METALLIB}"
        COMMAND xcrun -sdk macosx metallib "${TQ_AIR}" -o "${TQ_METALLIB}"
        DEPENDS "${TQ_AIR}"
        COMMENT "Linking TurboQuant.metallib"
        VERBATIM
    )
    add_custom_target(TurboQuantMetal DEPENDS "${TQ_METALLIB}")
    add_dependencies(FlexAID TurboQuantMetal)
    target_compile_definitions(FlexAID PRIVATE
        TURBOQUANT_METALLIB_PATH="${TQ_METALLIB}"
    )
    # MetalRMSD bridge (GPU pairwise RMSD for FOPTICS clustering)
    target_sources(FlexAID PRIVATE LIB/MetalRMSDBridge.mm)
    set_source_files_properties(LIB/MetalRMSDBridge.mm PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    # Compile MetalRMSD.metal shader → .metallib
    set(RMSD_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/MetalRMSD.metal")
    set(RMSD_AIR        "${CMAKE_CURRENT_BINARY_DIR}/MetalRMSD.air")
    set(RMSD_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/MetalRMSD.metallib")
    add_custom_command(
        OUTPUT  "${RMSD_AIR}"
        COMMAND xcrun -sdk macosx metal -c "${RMSD_METAL_SRC}" -o "${RMSD_AIR}"
        DEPENDS "${RMSD_METAL_SRC}"
        COMMENT "Compiling MetalRMSD Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${RMSD_METALLIB}"
        COMMAND xcrun -sdk macosx metallib "${RMSD_AIR}" -o "${RMSD_METALLIB}"
        DEPENDS "${RMSD_AIR}"
        COMMENT "Linking MetalRMSD.metallib"
        VERBATIM
    )
    add_custom_target(MetalRMSDMetal DEPENDS "${RMSD_METALLIB}")
    add_dependencies(FlexAID MetalRMSDMetal)
    target_compile_definitions(FlexAID PRIVATE
        METALRMSD_METALLIB_PATH="${RMSD_METALLIB}"
    )
endif()

# ─── Python bindings (pybind11 module) ─────────────────────────────────────────
add_subdirectory(python)

# ─── TeNCoM benchmark (optional standalone binary) ────────────────────────────────
option(ENABLE_TENCOM_BENCHMARK "Build the TeNCoM benchmark binary" OFF)
if(ENABLE_TENCOM_BENCHMARK)
    add_executable(benchmark_tencom
        LIB/tENCoM/benchmark_tencom.cpp
        LIB/tENCoM/tencm.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/encom.cpp
    )
    target_include_directories(benchmark_tencom PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM)
    if(MSVC)
        target_compile_options(benchmark_tencom PRIVATE /O2 /fp:fast)
    else()
        target_compile_options(benchmark_tencom PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(benchmark_tencom)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(benchmark_tencom PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(benchmark_tencom PRIVATE Eigen3::Eigen)

    target_compile_definitions(benchmark_tencom PRIVATE FLEXAIDS_HAS_EIGEN)
    if(MSVC)
        target_compile_definitions(benchmark_tencom PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    endif()
    message(STATUS "TeNCoM benchmark enabled — build target: benchmark_tencom")
endif()

# ─── tENCoM: ultra-fast vibrational entropy tool (HPC-ready) ──────────────────
# Enabled by default. Same LTO + native optimisation as FlexAIDdS.
option(ENABLE_TENCOM_TOOL "Build the tENCoM vibrational entropy differential tool" ON)
if(ENABLE_TENCOM_TOOL)
    add_executable(tENCoM
        LIB/tENCoM/tencom_main.cpp
        LIB/tENCoM/tencom_diff.cpp
        LIB/tENCoM/tencom_output.cpp
        LIB/tENCoM/pdb_calpha.cpp
        LIB/tENCoM/tencm.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/encom.cpp
    )
    target_include_directories(tENCoM PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack)
    if(MSVC)
        target_compile_options(tENCoM PRIVATE /O2 /GL /fp:fast /DNDEBUG /MT)
        target_compile_definitions(tENCoM PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
        target_link_options(tENCoM PRIVATE /LTCG)
    else()
        target_compile_options(tENCoM PRIVATE -O3 -ffast-math -fno-finite-math-only -flto -DNDEBUG)
        if(NOT APPLE)
            target_compile_options(tENCoM PRIVATE -march=native)
            target_link_options(tENCoM PRIVATE -flto -s)
        else()
            target_link_options(tENCoM PRIVATE -flto)
        endif()
    endif()
    flexaids_configure_simd(tENCoM)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(tENCoM PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(tENCoM PRIVATE Eigen3::Eigen)

    target_compile_definitions(tENCoM PRIVATE FLEXAIDS_HAS_EIGEN)
    if(FLEXAIDS_USE_ROCM)
        target_compile_definitions(tENCoM PRIVATE FLEXAIDS_USE_ROCM)
        target_include_directories(tENCoM PRIVATE ${hip_INCLUDE_DIRS})
        target_link_libraries(tENCoM PRIVATE hip::host)
    endif()
    # CUDA GPU acceleration for large proteins (N > 256 residues)
    if(FLEXAIDS_USE_CUDA)
        target_sources(tENCoM PRIVATE LIB/tENCoM/tencm_cuda.cu)
        target_compile_definitions(tENCoM PRIVATE FLEXAIDS_USE_CUDA)
        set_target_properties(tENCoM PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_STANDARD 17)
        target_link_libraries(tENCoM PRIVATE CUDA::cudart CUDA::cudadevrt)
    endif()
    # Metal GPU acceleration (macOS)
    if(FLEXAIDS_USE_METAL)
        target_sources(tENCoM PRIVATE LIB/tENCoM/tencm_metal.mm)
        target_compile_definitions(tENCoM PRIVATE FLEXAIDS_HAS_METAL_TENCM)
        set_source_files_properties(LIB/tENCoM/tencm_metal.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc")
        target_link_libraries(tENCoM PRIVATE
            ${METAL_LIBRARY} ${FOUNDATION_LIBRARY} ${METALKIT_LIBRARY})
    endif()
    message(STATUS "tENCoM tool enabled — build target: tENCoM")

    add_executable(tencom_entropy_diff
        LIB/tENCoM/tencom_entropy_diff.cpp
        LIB/tENCoM/pdb_calpha.cpp
        LIB/tENCoM/tencm.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/encom.cpp
    )
    target_include_directories(tencom_entropy_diff PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM)
    if(MSVC)
        target_compile_options(tencom_entropy_diff PRIVATE /O2 /GL /fp:fast /DNDEBUG /MT)
        target_compile_definitions(tencom_entropy_diff PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
        target_link_options(tencom_entropy_diff PRIVATE /LTCG)
    else()
        target_compile_options(tencom_entropy_diff PRIVATE -O3 -ffast-math -fno-finite-math-only -flto -DNDEBUG)
        if(NOT APPLE)
            target_compile_options(tencom_entropy_diff PRIVATE -march=native)
            target_link_options(tencom_entropy_diff PRIVATE -flto -s)
        else()
            target_link_options(tencom_entropy_diff PRIVATE -flto)
        endif()
    endif()
    flexaids_configure_simd(tencom_entropy_diff)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(tencom_entropy_diff PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(tencom_entropy_diff PRIVATE Eigen3::Eigen)
    target_compile_definitions(tencom_entropy_diff PRIVATE FLEXAIDS_HAS_EIGEN)
    message(STATUS "Legacy tENCoM entropy diff tool enabled — build target: tencom_entropy_diff")
endif()

# ─── VoronoiCFBatch benchmark (Task 2.3) ───────────────────────────────────────
# VoronoiCFBatch.h is header-only; the benchmark binary exercises the std::span
# batch interface and prints a serial-vs-OpenMP speedup table.
option(ENABLE_DIFT_TOOL "Build the DiFT torsional parametrization CLI" ON)
if(ENABLE_DIFT_TOOL)
    add_executable(dift_torsion_fit
        LIB/DiFT/dift_cli.cpp
        LIB/DiFT/DiFT.cpp
    )
    target_include_directories(dift_torsion_fit PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    if(MSVC)
        target_compile_options(dift_torsion_fit PRIVATE /O2)
        target_compile_definitions(dift_torsion_fit PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(dift_torsion_fit PRIVATE -O3)
    endif()
    flexaids_configure_simd(dift_torsion_fit)
    message(STATUS "DiFT torsion fit tool enabled — build target: dift_torsion_fit")
endif()

option(ENABLE_VCFBATCH_BENCHMARK "Build the VoronoiCFBatch benchmark binary" OFF)
if(ENABLE_VCFBATCH_BENCHMARK)
    add_executable(benchmark_vcfbatch LIB/benchmark_vcfbatch.cpp)
    target_include_directories(benchmark_vcfbatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    if(MSVC)
        target_compile_options(benchmark_vcfbatch PRIVATE /O2 /fp:fast)
    else()
        target_compile_options(benchmark_vcfbatch PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(benchmark_vcfbatch)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(benchmark_vcfbatch PRIVATE OpenMP::OpenMP_CXX)
    endif()
    message(STATUS "VoronoiCFBatch benchmark enabled — build target: benchmark_vcfbatch")
endif()

# ─── FlexAIDdS: ultra-fast docking executable ─────────────────────────────
# P0 validator references for GPU/metal shaders (actual add in Metal block above)
#   LIB/gpu_fast_optics_metal.metal LIB/gpu_batch_distance_hip.hip etc.
# Uses flexaid_core OBJECT library (same sources as FlexAID) with aggressive
# LTO and native-arch flags for production deployments.
if(BUILD_FLEXAIDDS_FAST)
    add_executable(FlexAIDdS
        LIB/top.cpp
        $<TARGET_OBJECTS:flexaid_core>
    )
    target_link_libraries(FlexAIDdS PRIVATE flexaid_core)

    if(MSVC)
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_options(FlexAIDdS PRIVATE /O2 /GL /DNDEBUG /MT)
        else()
            target_compile_options(FlexAIDdS PRIVATE /Od /MT)
        endif()
        target_link_options(FlexAIDdS PRIVATE /LTCG)
        target_compile_definitions(FlexAIDdS PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            target_compile_options(FlexAIDdS PRIVATE -O3 -ffast-math -fno-finite-math-only -flto -DNDEBUG
                $<$<AND:$<CXX_COMPILER_ID:AppleClang,Clang>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>>:-mcpu=native>
            )
        else()
            target_compile_options(FlexAIDdS PRIVATE -O2 -ffast-math -fno-finite-math-only -g
                $<$<AND:$<CXX_COMPILER_ID:AppleClang,Clang>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>>:-mcpu=native>
            )
        endif()
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            if(NOT APPLE)
                target_compile_options(FlexAIDdS PRIVATE -march=native)
                target_link_options(FlexAIDdS PRIVATE -flto -s)
            else()
                target_link_options(FlexAIDdS PRIVATE -flto)
            endif()
        elseif(NOT APPLE)
            target_compile_options(FlexAIDdS PRIVATE -march=native)
        endif()
    endif()

    # ─── IPO / LTO via CMake property (supplements the -flto flag above) ────
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_ok OUTPUT ipo_err)
    if(ipo_ok AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        set_property(TARGET FlexAIDdS PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
        set_property(TARGET flexaid_core PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(WARNING "IPO not supported: ${ipo_err}")
    endif()

    flexaids_configure_simd(FlexAIDdS)

    target_link_libraries(FlexAIDdS PRIVATE Eigen3::Eigen)

    # Stage the runtime data files (interaction matrix + atom-type definitions)
    # next to the binary. FlexAID/FlexAIDdS resolve these relative to the
    # executable's own directory (base_path), so without this they are not found
    # unless the binary happens to run from the repo root.
    add_custom_command(TARGET FlexAIDdS POST_BUILD
        COMMENT "Staging FlexAIDdS runtime data files next to the binary"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/MC_st0r5.2_6.dat"
            "${CMAKE_SOURCE_DIR}/AMINO.def"
            "${CMAKE_SOURCE_DIR}/NUCLEOTIDES.def"
            "$<TARGET_FILE_DIR:FlexAIDdS>")

    target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_HAS_EIGEN)

    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(FlexAIDdS PRIVATE OpenMP::OpenMP_CXX)
    endif()

    if(FLEXAIDS_USE_ROCM)
        target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_USE_ROCM)
        target_include_directories(FlexAIDdS PRIVATE ${hip_INCLUDE_DIRS})
        target_link_libraries(FlexAIDdS PRIVATE hip::host hip::device)
        set_property(TARGET FlexAIDdS PROPERTY HIP_ARCHITECTURES ${FLEXAIDS_HIP_ARCHITECTURES})
    endif()

    # ─── CUDA sources for FlexAIDdS (mirrors FlexAID GPU kernel set) ─────────
    if(FLEXAIDS_USE_CUDA)
        target_sources(FlexAIDdS PRIVATE
            LIB/cuda_eval.cu
            LIB/ShannonThermoStack/shannon_cuda.cu
            LIB/tENCoM/tencm_cuda.cu
            LIB/gpu_fast_optics.cu
            LIB/TurboQuant.cu
            LIB/gpu_batch_distance.cu
            LIB/gpu_batch_distance.cpp
        )
        target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_USE_CUDA)
        set_target_properties(FlexAIDdS PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_STANDARD 17)
        target_link_libraries(FlexAIDdS PRIVATE CUDA::cudart CUDA::cudadevrt)
        message(STATUS "FlexAIDdS: CUDA enabled — CF batch + Shannon + TENCoM kernels")
    endif()

    # ─── Metal sources for FlexAIDdS ────────────────────────────────────────
    if(FLEXAIDS_USE_METAL)
        # Metal batch Voronoi evaluation bridge
        target_sources(FlexAIDdS PRIVATE LIB/metal_eval.mm)
        target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_USE_METAL)
        target_compile_options(FlexAIDdS PRIVATE -fno-objc-arc)
        target_link_libraries(FlexAIDdS PRIVATE
            ${METAL_LIBRARY}
            ${FOUNDATION_LIBRARY}
            ${METALKIT_LIBRARY}
        )
        set_source_files_properties(LIB/metal_eval.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        # Shannon Metal GPU bridge (histogram binning on Apple Silicon)
        target_sources(FlexAIDdS PRIVATE LIB/ShannonThermoStack/ShannonMetalBridge.mm)
        target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_HAS_METAL_SHANNON)
        set_source_files_properties(LIB/ShannonThermoStack/ShannonMetalBridge.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        # Shannon Metal shader → .metallib
        set(SHANNON_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack/shannon_metal.metal")
        set(SHANNON_AIR        "${CMAKE_CURRENT_BINARY_DIR}/ShannonEntropy.air")
        set(SHANNON_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/ShannonEntropy.metallib")
        add_custom_command(
            OUTPUT  "${SHANNON_AIR}"
            COMMAND xcrun -sdk macosx metal -c "${SHANNON_METAL_SRC}" -o "${SHANNON_AIR}"
            DEPENDS "${SHANNON_METAL_SRC}"
            COMMENT "Compiling ShannonEntropy Metal shader"
            VERBATIM
        )
        add_custom_command(
            OUTPUT  "${SHANNON_METALLIB}"
            COMMAND xcrun -sdk macosx metallib "${SHANNON_AIR}" -o "${SHANNON_METALLIB}"
            DEPENDS "${SHANNON_AIR}"
            COMMENT "Linking ShannonEntropy.metallib"
            VERBATIM
        )
        add_custom_target(ShannonEntropyMetal DEPENDS "${SHANNON_METALLIB}")
        add_dependencies(FlexAIDdS ShannonEntropyMetal)
        target_compile_definitions(FlexAIDdS PRIVATE
            SHANNON_METALLIB_PATH="${SHANNON_METALLIB}"
        )
        # TENCoM Metal bridge (GPU contact discovery + Hessian assembly)
        target_sources(FlexAIDdS PRIVATE LIB/tENCoM/tencm_metal.mm)
        target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_HAS_METAL_TENCM)
        set_source_files_properties(LIB/tENCoM/tencm_metal.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        # FastOPTICS Metal kNN (GPU-accelerated neighbor discovery)
        target_sources(FlexAIDdS PRIVATE LIB/gpu_fast_optics_metal.mm)
        set_source_files_properties(LIB/gpu_fast_optics_metal.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        # CavityDetect Metal bridge + shader
        target_sources(FlexAIDdS PRIVATE LIB/CavityDetect/CavityDetectMetalBridge.mm)
        set_source_files_properties(LIB/CavityDetect/CavityDetectMetalBridge.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        set(CAVITY_METAL_SRC_DS  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect/CavityDetect.metal")
        set(CAVITY_AIR_DS         "${CMAKE_CURRENT_BINARY_DIR}/CavityDetect_ds.air")
        set(CAVITY_METALLIB_DS    "${CMAKE_CURRENT_BINARY_DIR}/CavityDetect_ds.metallib")
        add_custom_command(
            OUTPUT  "${CAVITY_AIR_DS}"
            COMMAND xcrun -sdk macosx metal -c "${CAVITY_METAL_SRC_DS}" -o "${CAVITY_AIR_DS}"
            DEPENDS "${CAVITY_METAL_SRC_DS}"
            COMMENT "Compiling CavityDetect Metal shader (FlexAIDdS)"
            VERBATIM
        )
        add_custom_command(
            OUTPUT  "${CAVITY_METALLIB_DS}"
            COMMAND xcrun -sdk macosx metallib "${CAVITY_AIR_DS}" -o "${CAVITY_METALLIB_DS}"
            DEPENDS "${CAVITY_AIR_DS}"
            COMMENT "Linking CavityDetect_ds.metallib"
            VERBATIM
        )
        add_custom_target(CavityDetectMetalDS DEPENDS "${CAVITY_METALLIB_DS}")
        add_dependencies(FlexAIDdS CavityDetectMetalDS)
        target_compile_definitions(FlexAIDdS PRIVATE
            CAVITY_METALLIB_PATH="${CAVITY_METALLIB_DS}"
        )
        # TurboQuant Metal bridge + shader
        target_sources(FlexAIDdS PRIVATE LIB/TurboQuantMetalBridge.mm)
        set_source_files_properties(LIB/TurboQuantMetalBridge.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        set(TQ_METAL_SRC_DS  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/TurboQuant.metal")
        set(TQ_AIR_DS         "${CMAKE_CURRENT_BINARY_DIR}/TurboQuant_ds.air")
        set(TQ_METALLIB_DS    "${CMAKE_CURRENT_BINARY_DIR}/TurboQuant_ds.metallib")
        add_custom_command(
            OUTPUT  "${TQ_AIR_DS}"
            COMMAND xcrun -sdk macosx metal -c "${TQ_METAL_SRC_DS}" -o "${TQ_AIR_DS}"
            DEPENDS "${TQ_METAL_SRC_DS}"
            COMMENT "Compiling TurboQuant Metal shader (FlexAIDdS)"
            VERBATIM
        )
        add_custom_command(
            OUTPUT  "${TQ_METALLIB_DS}"
            COMMAND xcrun -sdk macosx metallib "${TQ_AIR_DS}" -o "${TQ_METALLIB_DS}"
            DEPENDS "${TQ_AIR_DS}"
            COMMENT "Linking TurboQuant_ds.metallib"
            VERBATIM
        )
        add_custom_target(TurboQuantMetalDS DEPENDS "${TQ_METALLIB_DS}")
        add_dependencies(FlexAIDdS TurboQuantMetalDS)
        target_compile_definitions(FlexAIDdS PRIVATE
            TURBOQUANT_METALLIB_PATH="${TQ_METALLIB_DS}"
        )
        # MetalRMSD bridge + shader (GPU pairwise RMSD for FOPTICS clustering)
        target_sources(FlexAIDdS PRIVATE LIB/MetalRMSDBridge.mm)
        set_source_files_properties(LIB/MetalRMSDBridge.mm PROPERTIES
            LANGUAGE OBJCXX
            COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        set(RMSD_METAL_SRC_DS  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/MetalRMSD.metal")
        set(RMSD_AIR_DS         "${CMAKE_CURRENT_BINARY_DIR}/MetalRMSD_ds.air")
        set(RMSD_METALLIB_DS    "${CMAKE_CURRENT_BINARY_DIR}/MetalRMSD_ds.metallib")
        add_custom_command(
            OUTPUT  "${RMSD_AIR_DS}"
            COMMAND xcrun -sdk macosx metal -c "${RMSD_METAL_SRC_DS}" -o "${RMSD_AIR_DS}"
            DEPENDS "${RMSD_METAL_SRC_DS}"
            COMMENT "Compiling MetalRMSD Metal shader (FlexAIDdS)"
            VERBATIM
        )
        add_custom_command(
            OUTPUT  "${RMSD_METALLIB_DS}"
            COMMAND xcrun -sdk macosx metallib "${RMSD_AIR_DS}" -o "${RMSD_METALLIB_DS}"
            DEPENDS "${RMSD_AIR_DS}"
            COMMENT "Linking MetalRMSD_ds.metallib"
            VERBATIM
        )
        add_custom_target(MetalRMSDMetalDS DEPENDS "${RMSD_METALLIB_DS}")
        add_dependencies(FlexAIDdS MetalRMSDMetalDS)
        target_compile_definitions(FlexAIDdS PRIVATE
            METALRMSD_METALLIB_PATH="${RMSD_METALLIB_DS}"
        )
        message(STATUS "FlexAIDdS: Metal GPU acceleration enabled (Shannon + CavityDetect + TurboQuant + MetalRMSD)")
    endif()

    message(STATUS "FlexAIDdS fast executable enabled")
endif()

# ─── NATURaL Hammerhead Ribozyme simulation (standalone) ─────────────────────
# Runs the Zhao 2011 RNAP master equation on the 40-nt Schistosoma mansoni
# minimal hammerhead ribozyme in co-transcriptional (RNAP) mode.
# Only depends on RibosomeElongation (pure C++20, no FlexAID engine needed).
option(ENABLE_NATURAL_HAMMERHEAD "Build NATURaL co-transcriptional HHRz simulation" OFF)
if(ENABLE_NATURAL_HAMMERHEAD)
    add_executable(natural_hammerhead
        tools/natural_hammerhead.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
    )
    target_include_directories(natural_hammerhead PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
    )
    if(MSVC)
        target_compile_options(natural_hammerhead PRIVATE /O2 /fp:fast)
        target_compile_definitions(natural_hammerhead PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(natural_hammerhead PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(natural_hammerhead)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(natural_hammerhead PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(natural_hammerhead PRIVATE Eigen3::Eigen)

    target_compile_definitions(natural_hammerhead PRIVATE FLEXAIDS_HAS_EIGEN)
    message(STATUS "NATURaL HHRz simulation enabled — build target: natural_hammerhead")
endif()

# ─── DualAssembly cotranslational docking CLI ────────────────────────────────
# Standalone driver for protofibril × nascent-chain cotranslational docking with the
# Shannon-driven Target/Ligand discriminator. See docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md.
# MVP ships a built-in synthetic GA backend so the full pipeline can be exercised on
# any platform without the heavyweight FlexAID engine.
option(ENABLE_DUAL_ASSEMBLY_TOOL "Build the DualAssembly cotranslational docking CLI" ON)
if(ENABLE_DUAL_ASSEMBLY_TOOL)
    add_executable(dual_assembly
        LIB/NATURaL/dual_assembly_main.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
    )
    target_include_directories(dual_assembly PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
    )
    if(MSVC)
        target_compile_options(dual_assembly PRIVATE /O2 /fp:fast)
        target_compile_definitions(dual_assembly PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(dual_assembly PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(dual_assembly)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(dual_assembly PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(dual_assembly PRIVATE Eigen3::Eigen)
    target_compile_definitions(dual_assembly PRIVATE
        FLEXAIDS_HAS_EIGEN
        FLEXAIDS_OMIT_DUAL_ASSEMBLY_ENGINE)
    message(STATUS "DualAssembly CLI enabled — build target: dual_assembly")
endif()

# ─── Dispatch Benchmark (optional standalone binary) ─────────────────────────
option(ENABLE_DISPATCH_BENCHMARK "Build the hardware dispatch benchmark binary" OFF)
if(ENABLE_DISPATCH_BENCHMARK)
    add_executable(benchmark_dispatch
        LIB/benchmark_dispatch.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp   # provides detect_hardware() + HardwareCapabilities
    )
    target_include_directories(benchmark_dispatch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    if(MSVC)
        target_compile_options(benchmark_dispatch PRIVATE /O2 /fp:fast)
    else()
        target_compile_options(benchmark_dispatch PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(benchmark_dispatch)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(benchmark_dispatch PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(benchmark_dispatch PRIVATE Eigen3::Eigen)

    target_compile_definitions(benchmark_dispatch PRIVATE FLEXAIDS_HAS_EIGEN)
    message(STATUS "Dispatch benchmark enabled — build target: benchmark_dispatch")
endif()

endif() # FLEXAIDS_BUILD_CORE

# ─── Benchmark Dataset Runner ──────────────────────────────────────────────
option(ENABLE_BENCHMARK_DATASETS "Build benchmark dataset runner" ON)
if(ENABLE_BENCHMARK_DATASETS)
    find_package(CURL QUIET)

    add_executable(benchmark_datasets
        LIB/benchmark_datasets.cpp
        LIB/DatasetRunner.cpp
        LIB/AsyncPipeline.cpp
    )
    target_link_libraries(benchmark_datasets PRIVATE flexaid_core)
    target_include_directories(benchmark_datasets PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    if(MSVC)
        target_compile_options(benchmark_datasets PRIVATE /O2 /fp:fast)
        target_compile_definitions(benchmark_datasets PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(benchmark_datasets PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(benchmark_datasets)
    if(CURL_FOUND)
        target_link_libraries(benchmark_datasets PRIVATE CURL::libcurl)
        target_compile_definitions(benchmark_datasets PRIVATE FLEXAIDS_HAS_LIBCURL)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(benchmark_datasets PRIVATE OpenMP::OpenMP_CXX)
    endif()
    # Metal bridge sources — required when flexaid_core was compiled with Metal.
    if(FLEXAIDS_USE_METAL)
        target_sources(benchmark_datasets PRIVATE
            LIB/metal_eval.mm
            LIB/MetalRMSDBridge.mm
            LIB/CavityDetect/CavityDetectMetalBridge.mm
            LIB/ShannonThermoStack/ShannonMetalBridge.mm
            LIB/gpu_fast_optics_metal.mm
        )
        target_compile_definitions(benchmark_datasets PRIVATE
            FLEXAIDS_USE_METAL
            FLEXAIDS_HAS_METAL_SHANNON
        )
        target_compile_options(benchmark_datasets PRIVATE -fno-objc-arc)
        target_link_libraries(benchmark_datasets PRIVATE
            ${METAL_LIBRARY}
            ${FOUNDATION_LIBRARY}
            ${METALKIT_LIBRARY}
        )
        set_source_files_properties(
            LIB/metal_eval.mm
            LIB/MetalRMSDBridge.mm
            LIB/CavityDetect/CavityDetectMetalBridge.mm
            LIB/ShannonThermoStack/ShannonMetalBridge.mm
            LIB/gpu_fast_optics_metal.mm
            PROPERTIES LANGUAGE OBJCXX COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
        )
        add_dependencies(benchmark_datasets
            CavityDetectMetalDS
            ShannonEntropyMetal
            MetalRMSDMetalDS
        )
        target_compile_definitions(benchmark_datasets PRIVATE
            SHANNON_METALLIB_PATH="${SHANNON_METALLIB}"
            METALRMSD_METALLIB_PATH="${RMSD_METALLIB_DS}"
            CAVITY_METALLIB_PATH="${CAVITY_METALLIB_DS}"
            TURBOQUANT_METALLIB_PATH="${TQ_METALLIB_DS}"
        )
    endif()
    message(STATUS "Benchmark dataset runner enabled — build target: benchmark_datasets")
endif()

# ─── FlexAID 2015 Benchmark Convenience Targets ────────────────────────────
# These custom targets invoke the Python DatasetRunner for the FlexAID JCIM
# 2015 validation sequence: Astex Diverse, Astex Non-Native, then HAP2.
# They require Python and the flexaidds package installed.
# Usage: cmake --build build --target flexaid_bench_astex_diverse
#         cmake --build build --target flexaid_bench_astex_non_native
#        cmake --build build --target flexaid_bench_hap2
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_FOUND)
    add_custom_target(flexaid_bench_astex_diverse
        COMMAND ${Python3_EXECUTABLE} -m flexaidds.dataset_runner
                --dataset astex_diverse --tier 2
                --results-dir ${CMAKE_BINARY_DIR}/benchmark_results/astex_diverse
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running Astex Diverse Set benchmark (85 complexes, tier-2)"
        VERBATIM
    )

    add_custom_target(flexaid_bench_astex_non_native
        COMMAND ${Python3_EXECUTABLE} -m flexaidds.dataset_runner
                --dataset astex_nonnative --tier 2
                --results-dir ${CMAKE_BINARY_DIR}/benchmark_results/astex_nonnative
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running Astex Non-Native cross-docking benchmark (65 targets, tier-2)"
        VERBATIM
    )

    add_custom_target(flexaid_bench_hap2
        COMMAND ${Python3_EXECUTABLE} -m flexaidds.dataset_runner
                --dataset hap2 --tier 2
                --results-dir ${CMAKE_BINARY_DIR}/benchmark_results/hap2
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Running HAP2 FlexAID JCIM 2015 benchmark (tier-2)"
        VERBATIM
    )

    message(STATUS "FlexAID 2015 benchmark targets: flexaid_bench_astex_diverse, flexaid_bench_astex_non_native, flexaid_bench_hap2")
else()
    message(STATUS "Python3 not found — FlexAID 2015 benchmark custom targets skipped")
endif()

# ─── Unit Tests (GoogleTest) ───────────────────────────────────────────────
option(BUILD_TESTING "Build unit tests (requires GoogleTest)" OFF)
if(BUILD_TESTING)
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        v1.17.0
        )
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()

    enable_testing()

    # Helper: apply MSVC-specific definitions for test targets
    function(flexaids_configure_msvc_test target_name)
        if(MSVC)
            target_compile_definitions(${target_name} PRIVATE
                _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
        endif()
    endfunction()

    # test_statmech — StatMechEngine unit tests
    add_executable(test_statmech tests/test_statmech.cpp LIB/statmech.cpp LIB/ReportedPoseSelector.cpp
        LIB/TurboQuant.cpp
                   LIB/UnifiedHardwareDispatch.cpp LIB/hardware_detect.cpp)
    target_include_directories(test_statmech PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_statmech PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_statmech PRIVATE /O2)
    else()
        target_compile_options(test_statmech PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_statmech)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_statmech PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_statmech PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_statmech PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_statmech)
    add_test(NAME StatMechTests COMMAND test_statmech)

    add_executable(test_thermo_ledger tests/test_thermo_ledger.cpp LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp LIB/hardware_detect.cpp)
    target_include_directories(test_thermo_ledger PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_thermo_ledger PRIVATE GTest::gtest GTest::gtest_main Eigen3::Eigen)
    if(MSVC)
        target_compile_options(test_thermo_ledger PRIVATE /O2)
    else()
        target_compile_options(test_thermo_ledger PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_thermo_ledger)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_thermo_ledger PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_compile_definitions(test_thermo_ledger PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_thermo_ledger)
    add_test(NAME ThermoLedgerTests COMMAND test_thermo_ledger)

    # test_json_config — JSON parser & config defaults tests (zero dependencies)
    add_executable(test_json_config tests/test_json_config.cpp)
    target_include_directories(test_json_config PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_json_config PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_json_config PRIVATE /O2)
    else()
        target_compile_options(test_json_config PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_json_config)
    add_test(NAME JsonConfigTests COMMAND test_json_config)

    # test_metal_coordination — Metal ion coordination potential tests (header-only)
    add_executable(test_metal_coordination tests/test_metal_coordination.cpp)
    target_include_directories(test_metal_coordination PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_metal_coordination PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_metal_coordination PRIVATE /O2)
    else()
        target_compile_options(test_metal_coordination PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_metal_coordination)
    add_test(NAME MetalCoordinationTests COMMAND test_metal_coordination)

    # test_formal_charges — PDB formal charge assignment tests (header-only)
    add_executable(test_formal_charges tests/test_formal_charges.cpp)
    target_include_directories(test_formal_charges PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_formal_charges PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_formal_charges PRIVATE /O2)
    else()
        target_compile_options(test_formal_charges PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_formal_charges)
    add_test(NAME FormalChargeTests COMMAND test_formal_charges)

    # test_hardware_dispatch — ShannonThermoStack dispatch layer tests
    add_executable(test_hardware_dispatch
        tests/test_hardware_dispatch.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/tENCoM/tencm.cpp)
    target_include_directories(test_hardware_dispatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack)
    target_link_libraries(test_hardware_dispatch PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_hardware_dispatch PRIVATE /O2)
    else()
        target_compile_options(test_hardware_dispatch PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_hardware_dispatch PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_hardware_dispatch PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_hardware_dispatch PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_hardware_dispatch)
    add_test(NAME HardwareDispatchTests COMMAND test_hardware_dispatch)

    # test_hardware_detect_dispatch — unified dispatch layer + hardware detection tests
    add_executable(test_hardware_detect_dispatch
        tests/test_hardware_detect_dispatch.cpp
        LIB/hardware_detect.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp)
    target_include_directories(test_hardware_detect_dispatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_hardware_detect_dispatch PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_hardware_detect_dispatch PRIVATE /O2)
    else()
        target_compile_options(test_hardware_detect_dispatch PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_hardware_detect_dispatch)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_hardware_detect_dispatch PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_hardware_detect_dispatch PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_hardware_detect_dispatch PRIVATE FLEXAIDS_HAS_EIGEN)
    if(FLEXAIDS_USE_CUDA)
        target_compile_definitions(test_hardware_detect_dispatch PRIVATE FLEXAIDS_USE_CUDA)
        target_link_libraries(test_hardware_detect_dispatch PRIVATE CUDA::cudart)
    endif()
    if(FLEXAIDS_USE_METAL)
        target_compile_definitions(test_hardware_detect_dispatch PRIVATE FLEXAIDS_USE_METAL)
    endif()
    if(FLEXAIDS_USE_ROCM)
        target_compile_definitions(test_hardware_detect_dispatch PRIVATE FLEXAIDS_USE_ROCM)
        target_include_directories(test_hardware_detect_dispatch PRIVATE ${hip_INCLUDE_DIRS})
        target_link_libraries(test_hardware_detect_dispatch PRIVATE hip::host)
    endif()
    flexaids_configure_msvc_test(test_hardware_detect_dispatch)
    add_test(NAME HardwareDetectDispatchTests COMMAND test_hardware_detect_dispatch)

    # test_tencom_diff — tENCoM differential engine tests
    add_executable(test_tencom_diff
        tests/test_tencom_diff.cpp
        LIB/tENCoM/tencom_diff.cpp
        LIB/tENCoM/pdb_calpha.cpp
        LIB/tENCoM/tencm.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/encom.cpp
    )
    target_include_directories(test_tencom_diff PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM)
    target_link_libraries(test_tencom_diff PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_tencom_diff PRIVATE /O2)
    else()
        target_compile_options(test_tencom_diff PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_tencom_diff)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_tencom_diff PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_tencom_diff PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_tencom_diff PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_tencom_diff)
    add_test(NAME TENCoMDiffTests COMMAND test_tencom_diff)

    # test_ptm_attachment — PTMAttachment module tests
    add_executable(test_ptm_attachment
        tests/test_ptm_attachment.cpp
        LIB/PTMAttachment/PTMAttachment.cpp
    )
    target_include_directories(test_ptm_attachment PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PTMAttachment
    )
    target_link_libraries(test_ptm_attachment PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ptm_attachment PRIVATE /O2)
    else()
        target_compile_options(test_ptm_attachment PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_ptm_attachment)
    add_test(NAME PTMAttachmentTests COMMAND test_ptm_attachment)

    # test_gaboom — Genetic Algorithm engine unit tests
    add_executable(test_gaboom
        tests/test_gaboom.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/shannon_ga.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_gaboom PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_gaboom PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_gaboom PRIVATE /O2)
    else()
        target_compile_options(test_gaboom PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_gaboom PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_gaboom PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_gaboom PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_gaboom)
    add_test(NAME GABoomTests COMMAND test_gaboom)

    # test_shannon_ga — Shannon Thermodynamic GA entropy-collapse tests
    add_executable(test_shannon_ga
        tests/test_shannon_ga.cpp
        tests/stubs.cpp
        LIB/shannon_ga.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_shannon_ga PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_shannon_ga PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_shannon_ga PRIVATE /O2)
    else()
        target_compile_options(test_shannon_ga PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_shannon_ga PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if(FLEXAIDS_USE_EIGEN AND Eigen3_FOUND)
        target_link_libraries(test_shannon_ga PRIVATE Eigen3::Eigen)
        target_compile_definitions(test_shannon_ga PRIVATE FLEXAIDS_HAS_EIGEN)
    endif()
    flexaids_configure_msvc_test(test_shannon_ga)
    add_test(NAME ShannonGATests COMMAND test_shannon_ga)

    # test_vcontacts — Voronoi contact scoring math unit tests
    add_executable(test_vcontacts
        tests/test_vcontacts.cpp
        tests/stubs.cpp
        LIB/Vcontacts.cpp
        LIB/geometry.cpp
    )
    target_include_directories(test_vcontacts PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_vcontacts PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_vcontacts PRIVATE /O2)
    else()
        target_compile_options(test_vcontacts PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_vcontacts PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_msvc_test(test_vcontacts)
    add_test(NAME VcontactsTests COMMAND test_vcontacts)

    # test_vcontacts_geometry — convex-hull geometry functions (test_point, add_vedge, etc.)
    add_executable(test_vcontacts_geometry
        tests/test_vcontacts_geometry.cpp
        tests/stubs.cpp
        LIB/Vcontacts.cpp
        LIB/geometry.cpp
    )
    target_include_directories(test_vcontacts_geometry PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_vcontacts_geometry PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_vcontacts_geometry PRIVATE /O2)
    else()
        target_compile_options(test_vcontacts_geometry PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_vcontacts_geometry PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_msvc_test(test_vcontacts_geometry)
    add_test(NAME VcontactsGeometryTests COMMAND test_vcontacts_geometry)

    # test_cleft_cavity — CleftDetector + CavityDetect unit tests
    add_executable(test_cleft_cavity
        tests/test_cleft_cavity.cpp
        LIB/CleftDetector.cpp
        LIB/CavityDetect/CavityDetect.cpp
        LIB/geometry.cpp
    )
    target_include_directories(test_cleft_cavity PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_cleft_cavity PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_cleft_cavity PRIVATE /O2)
    else()
        target_compile_options(test_cleft_cavity PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_cleft_cavity PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_msvc_test(test_cleft_cavity)
    add_test(NAME CleftCavityTests COMMAND test_cleft_cavity)

    # test_encom — ENCoM vibrational entropy engine tests
    add_executable(test_encom
        tests/test_encom.cpp
        LIB/encom.cpp
    )
    target_include_directories(test_encom PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_encom PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_encom PRIVATE /O2)
    else()
        target_compile_options(test_encom PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_encom)
    target_link_libraries(test_encom PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_encom PRIVATE FLEXAIDS_HAS_EIGEN)
    add_test(NAME ENCoMTests COMMAND test_encom)

    # test_ring_conformer_library — RingConformerLibrary + ring detection tests
    add_executable(test_ring_conformer_library
        tests/test_ring_conformer_library.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
    )
    target_include_directories(test_ring_conformer_library PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
    )
    target_link_libraries(test_ring_conformer_library PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ring_conformer_library PRIVATE /O2)
    else()
        target_compile_options(test_ring_conformer_library PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_ring_conformer_library)
    add_test(NAME RingConformerLibraryTests COMMAND test_ring_conformer_library)

    # test_sugar_pucker — SugarPucker Cremer-Pople pseudorotation tests
    add_executable(test_sugar_pucker
        tests/test_sugar_pucker.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
    )
    target_include_directories(test_sugar_pucker PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
    )
    target_link_libraries(test_sugar_pucker PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_sugar_pucker PRIVATE /O2)
    else()
        target_compile_options(test_sugar_pucker PRIVATE -O2)
    endif()
    target_link_libraries(test_sugar_pucker PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_sugar_pucker PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_sugar_pucker)
    add_test(NAME SugarPuckerTests COMMAND test_sugar_pucker)

    # test_chiral_center — ChiralCenterGene R/S stereocenter tests
    add_executable(test_chiral_center
        tests/test_chiral_center.cpp
        LIB/ChiralCenter/ChiralCenterGene.cpp
    )
    target_include_directories(test_chiral_center PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ChiralCenter
    )
    target_link_libraries(test_chiral_center PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_chiral_center PRIVATE /O2)
    else()
        target_compile_options(test_chiral_center PRIVATE -O2)
    endif()
    target_link_libraries(test_chiral_center PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_chiral_center PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_chiral_center)
    add_test(NAME ChiralCenterTests COMMAND test_chiral_center)

    # test_cavity_detect — CavityDetect SURFNET cavity detection tests
    add_executable(test_cavity_detect
        tests/test_cavity_detect.cpp
        LIB/CavityDetect/CavityDetect.cpp
    )
    target_include_directories(test_cavity_detect PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_cavity_detect PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_cavity_detect PRIVATE /O2)
    else()
        target_compile_options(test_cavity_detect PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_cavity_detect PRIVATE OpenMP::OpenMP_CXX)
    endif()
flexaids_configure_simd(test_cavity_detect)
    flexaids_configure_msvc_test(test_cavity_detect)
    add_test(NAME CavityDetectTests COMMAND test_cavity_detect)

    # test_mif_grid — MIFGrid molecular interaction field tests
    add_executable(test_mif_grid
        tests/test_mif_grid.cpp
    )
    target_include_directories(test_mif_grid PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_mif_grid PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_mif_grid PRIVATE /O2)
    else()
        target_compile_options(test_mif_grid PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_mif_grid PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_simd(test_mif_grid)
    add_test(NAME MIFGridTests COMMAND test_mif_grid)

    # test_reflig_spectrophore — RefLigSeed + Spectrophore tests
    add_executable(test_reflig_spectrophore
        tests/test_reflig_spectrophore.cpp
    )
    target_include_directories(test_reflig_spectrophore PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_reflig_spectrophore PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_reflig_spectrophore PRIVATE /O2)
    else()
        target_compile_options(test_reflig_spectrophore PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_reflig_spectrophore PRIVATE OpenMP::OpenMP_CXX)
    endif()
    add_test(NAME RefLigSpectrophoreTests COMMAND test_reflig_spectrophore)

    # test_binding_residues — Key binding-site residue identification tests
    add_executable(test_binding_residues
        tests/test_binding_residues.cpp
        LIB/set_intprob.cpp
    )
    target_include_directories(test_binding_residues PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_residues PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_binding_residues PRIVATE /O2)
    else()
        target_compile_options(test_binding_residues PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_binding_residues PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_simd(test_binding_residues)
    add_test(NAME BindingResiduesTests COMMAND test_binding_residues)

    # test_binding_mode_statmech — BindingMode ↔ StatMechEngine integration tests
    add_executable(test_binding_mode_statmech
        tests/test_binding_mode_statmech.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/statmech.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/LigandRingFlex/LigandRingFlex.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
    )
    target_include_directories(test_binding_mode_statmech PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_statmech PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_binding_mode_statmech PRIVATE /O2)
    else()
        target_compile_options(test_binding_mode_statmech PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_binding_mode_statmech PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_binding_mode_statmech PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_binding_mode_statmech PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_binding_mode_statmech)
    add_test(NAME BindingModeStatMechTests COMMAND test_binding_mode_statmech)

    # test_ccbm — Conformer-Coupled Binding Modes (CCBM) unit tests
    add_executable(test_ccbm
        tests/test_ccbm.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_ccbm PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ccbm PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ccbm PRIVATE /O2)
    else()
        target_compile_options(test_ccbm PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_ccbm PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_ccbm PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_ccbm PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_ccbm)
    add_test(NAME CCBMTests COMMAND test_ccbm)

    # test_binding_mode_advanced — BindingMode advanced methods (delta_G, WHAM, representative)
    add_executable(test_binding_mode_advanced
        tests/test_binding_mode_advanced.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_binding_mode_advanced PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_advanced PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_binding_mode_advanced PRIVATE /O2)
    else()
        target_compile_options(test_binding_mode_advanced PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_binding_mode_advanced PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_binding_mode_advanced PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_binding_mode_advanced PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_binding_mode_advanced)
    add_test(NAME BindingModeAdvancedTests COMMAND test_binding_mode_advanced)

    # test_ga_core — GA primitive operations (QuickSort, mutate, crossover, copy)
    add_executable(test_ga_core
        tests/test_ga_core.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/InStreamClustering.cpp
    )
    target_include_directories(test_ga_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ga_core PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ga_core PRIVATE /O2)
    else()
        target_compile_options(test_ga_core PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_ga_core PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_ga_core PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_ga_core PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_ga_core)
    add_test(NAME GACoreTests COMMAND test_ga_core)

    # test_ga_validation — CF scoring functions, eval_span clamping, batch result structure
    add_executable(test_ga_validation
        tests/test_ga_validation.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/InStreamClustering.cpp
    )
    target_include_directories(test_ga_validation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ga_validation PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ga_validation PRIVATE /O2)
    else()
        target_compile_options(test_ga_validation PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_ga_validation PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_ga_validation PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_ga_validation PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_ga_validation)
    add_test(NAME GAValidationTests COMMAND test_ga_validation)

    # test_ga_operators — untested GA helpers (calc_poss, set_bins, validate_dups, etc.)
    add_executable(test_ga_operators
        tests/test_ga_operators.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/calc_rmsp.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/InStreamClustering.cpp
    )
    target_include_directories(test_ga_operators PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ga_operators PRIVATE GTest::gtest GTest::gtest_main Eigen3::Eigen)
    target_compile_definitions(test_ga_operators PRIVATE FLEXAIDS_HAS_EIGEN)
    if(MSVC)
        target_compile_options(test_ga_operators PRIVATE /O2)
    else()
        target_compile_options(test_ga_operators PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_ga_operators PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_msvc_test(test_ga_operators)
    add_test(NAME GAOperatorsTests COMMAND test_ga_operators)

    # test_binding_mode_comparison — Pose/BindingMode comparison operators
    add_executable(test_binding_mode_comparison
        tests/test_binding_mode_comparison.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/calc_rmsp.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/InStreamClustering.cpp
    )
    target_include_directories(test_binding_mode_comparison PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_comparison PRIVATE GTest::gtest GTest::gtest_main Eigen3::Eigen)
    target_compile_definitions(test_binding_mode_comparison PRIVATE FLEXAIDS_HAS_EIGEN)
    if(MSVC)
        target_compile_options(test_binding_mode_comparison PRIVATE /O2)
    else()
        target_compile_options(test_binding_mode_comparison PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_binding_mode_comparison PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_msvc_test(test_binding_mode_comparison)
    add_test(NAME BindingModeComparisonTests COMMAND test_binding_mode_comparison)

    # test_soft_contact_matrix — 256×256 energy matrix I/O, projection, lookup
    if(FLEXAIDS_USE_256_MATRIX)
        add_executable(test_soft_contact_matrix
            tests/test_soft_contact_matrix.cpp
            LIB/ShannonThermoStack/ShannonThermoStack.cpp
            LIB/statmech.cpp
        LIB/TurboQuant.cpp
            LIB/encom.cpp
            LIB/tENCoM/tencm.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
        )
        target_include_directories(test_soft_contact_matrix PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        )
        target_link_libraries(test_soft_contact_matrix PRIVATE GTest::gtest GTest::gtest_main)
        target_compile_definitions(test_soft_contact_matrix PRIVATE FLEXAIDS_USE_256_MATRIX)
        if(MSVC)
            target_compile_options(test_soft_contact_matrix PRIVATE /O2)
        else()
            target_compile_options(test_soft_contact_matrix PRIVATE -O2)
        endif()
        if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
            target_link_libraries(test_soft_contact_matrix PRIVATE OpenMP::OpenMP_CXX)
        endif()
        target_link_libraries(test_soft_contact_matrix PRIVATE Eigen3::Eigen)

        target_compile_definitions(test_soft_contact_matrix PRIVATE FLEXAIDS_HAS_EIGEN)
        add_test(NAME SoftContactMatrixTests COMMAND test_soft_contact_matrix)
    endif()

    add_executable(test_matrix_v2_science tests/test_matrix_v2_science.cpp)
    target_link_libraries(test_matrix_v2_science PRIVATE GTest::gtest GTest::gtest_main)
    add_test(NAME MatrixV2ScienceTests COMMAND test_matrix_v2_science
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

    # test_compare_energy_matrices — compare 256×256 vs 40×40 energy matrices
    if(FLEXAIDS_USE_256_MATRIX)
        add_executable(test_compare_energy_matrices
            tests/test_compare_energy_matrices.cpp
        )
        target_include_directories(test_compare_energy_matrices PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        )
        target_link_libraries(test_compare_energy_matrices PRIVATE GTest::gtest GTest::gtest_main)
        target_compile_definitions(test_compare_energy_matrices PRIVATE FLEXAIDS_USE_256_MATRIX)
        if(MSVC)
            target_compile_options(test_compare_energy_matrices PRIVATE /O2)
        else()
            target_compile_options(test_compare_energy_matrices PRIVATE -O2)
        endif()
        if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
            target_link_libraries(test_compare_energy_matrices PRIVATE OpenMP::OpenMP_CXX)
        endif()
        target_link_libraries(test_compare_energy_matrices PRIVATE Eigen3::Eigen)

        target_compile_definitions(test_compare_energy_matrices PRIVATE FLEXAIDS_HAS_EIGEN)
        add_test(NAME CompareEnergyMatricesTests COMMAND test_compare_energy_matrices)
    endif()

    # test_ion_handling — assign_radii / assign_types HETATM + ion coverage
    #                     + tENCoM ion contact surface (TmContSct) tests
    add_executable(test_ion_handling
        tests/test_ion_handling.cpp
        LIB/assign_radii.cpp
        LIB/assign_types.cpp
        LIB/fileio.cpp
        LIB/tENCoM/tencm.cpp)
    target_include_directories(test_ion_handling PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_ion_handling PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ion_handling PRIVATE /O2)
    else()
        target_compile_options(test_ion_handling PRIVATE -O2)
    endif()
    target_link_libraries(test_ion_handling PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_ion_handling PRIVATE FLEXAIDS_HAS_EIGEN)
    add_test(NAME IonHandlingTests COMMAND test_ion_handling)

    # test_entropy_ga — SMFREE fitness model, Boltzmann blending, entropy fields
    add_executable(test_entropy_ga
        tests/test_entropy_ga.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_entropy_ga PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_entropy_ga PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_entropy_ga PRIVATE /O2)
    else()
        target_compile_options(test_entropy_ga PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_entropy_ga PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_entropy_ga PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_entropy_ga PRIVATE FLEXAIDS_HAS_EIGEN)
    add_test(NAME EntropyGATests COMMAND test_entropy_ga)

    # test_cube_grid — GridKey reproducibility, diagonal formula, neighbor probing
    add_executable(test_cube_grid
        tests/test_cube_grid.cpp
        LIB/maps.cpp
    )
    target_include_directories(test_cube_grid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_cube_grid PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_cube_grid PRIVATE /O2)
    else()
        target_compile_options(test_cube_grid PRIVATE -O2)
    endif()
    add_test(NAME CubeGridTests COMMAND test_cube_grid)

    # test_parallel_dock — GridDecomposer, SharedPosePool, StatMechEngine merge tests
    add_executable(test_parallel_dock
        tests/test_parallel_dock.cpp
        LIB/GridDecomposer.cpp
        LIB/SharedPosePool.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_parallel_dock PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_parallel_dock PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_parallel_dock PRIVATE /O2)
    else()
        target_compile_options(test_parallel_dock PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_parallel_dock PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_parallel_dock PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_parallel_dock PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_parallel_dock)
    add_test(NAME ParallelDockTests COMMAND test_parallel_dock)

    # test_knowledge_pool — TargetKnowledgeBase + SharedPosePool + MPIMergeable
    add_executable(test_knowledge_pool
        tests/test_knowledge_pool.cpp
        LIB/TargetKnowledgeBase.cpp
        LIB/SharedPosePool.cpp
    )
    target_include_directories(test_knowledge_pool PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_knowledge_pool PRIVATE GTest::gtest GTest::gtest_main)
    target_compile_options(test_knowledge_pool PRIVATE -std=c++26)
    if(MSVC)
        target_compile_options(test_knowledge_pool PRIVATE /O2)
    else()
        target_compile_options(test_knowledge_pool PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_knowledge_pool)
    add_test(NAME KnowledgePoolTests COMMAND test_knowledge_pool)

    # ── Mol2Reader / SdfReader tests ──────────────────────────────────────
    add_executable(test_mol2_sdf_reader
        tests/test_mol2_sdf_reader.cpp
        tests/stubs.cpp
        LIB/Mol2Reader.cpp
        LIB/SdfReader.cpp
    )
    target_include_directories(test_mol2_sdf_reader PRIVATE LIB)
    target_link_libraries(test_mol2_sdf_reader PRIVATE GTest::gtest_main)
    target_compile_options(test_mol2_sdf_reader PRIVATE -std=c++26)
    flexaids_configure_simd(test_mol2_sdf_reader)
    flexaids_configure_msvc_test(test_mol2_sdf_reader)
    add_test(NAME Mol2SdfReaderTests COMMAND test_mol2_sdf_reader)

    # test_turboquant — TurboQuant quantisation library unit tests
    add_executable(test_turboquant
        tests/test_turboquant.cpp
        LIB/TurboQuant.cpp
    )
    target_include_directories(test_turboquant PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_turboquant PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_turboquant PRIVATE /O2)
    else()
        target_compile_options(test_turboquant PRIVATE -O2)
    endif()
    target_link_libraries(test_turboquant PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_turboquant PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_turboquant)
    add_test(NAME TurboQuantTests COMMAND test_turboquant)

    # test_dataset_runner — benchmark dataset runner tests
    add_executable(test_dataset_runner
        tests/test_dataset_runner.cpp
        LIB/DatasetRunner.cpp
        LIB/AsyncPipeline.cpp
    )
    target_include_directories(test_dataset_runner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_dataset_runner PRIVATE flexaid_core GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_dataset_runner PRIVATE /O2)
        target_compile_definitions(test_dataset_runner PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(test_dataset_runner PRIVATE -O2)
    endif()
    if(CURL_FOUND)
        target_link_libraries(test_dataset_runner PRIVATE CURL::libcurl)
        target_compile_definitions(test_dataset_runner PRIVATE FLEXAIDS_HAS_LIBCURL)
    endif()
    flexaids_configure_msvc_test(test_dataset_runner)
    add_test(NAME DatasetRunnerTests COMMAND test_dataset_runner)

    # test_calc_rmsd_chrom_bounds — OOB snapshot genes must not crash coord cache
    add_executable(test_calc_rmsd_chrom_bounds
        tests/test_calc_rmsd_chrom_bounds.cpp
    )
    target_include_directories(test_calc_rmsd_chrom_bounds PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_calc_rmsd_chrom_bounds PRIVATE flexaid_core GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_calc_rmsd_chrom_bounds PRIVATE /O2)
    else()
        target_compile_options(test_calc_rmsd_chrom_bounds PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_calc_rmsd_chrom_bounds)
    add_test(NAME CalcRmsdChromBoundsTests COMMAND test_calc_rmsd_chrom_bounds)

    # test_cofactor_blacklist — ligand-extraction cofactor blacklist (1HNN SAH/SKF)
    add_executable(test_cofactor_blacklist
        tests/test_cofactor_blacklist.cpp
        LIB/DatasetRunner.cpp
        LIB/AsyncPipeline.cpp
    )
    target_include_directories(test_cofactor_blacklist PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_cofactor_blacklist PRIVATE flexaid_core GTest::gtest GTest::gtest_main)
    target_compile_definitions(test_cofactor_blacklist PRIVATE
        FLEXAIDS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    if(MSVC)
        target_compile_options(test_cofactor_blacklist PRIVATE /O2)
        target_compile_definitions(test_cofactor_blacklist PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(test_cofactor_blacklist PRIVATE -O2)
    endif()
    if(CURL_FOUND)
        target_link_libraries(test_cofactor_blacklist PRIVATE CURL::libcurl)
        target_compile_definitions(test_cofactor_blacklist PRIVATE FLEXAIDS_HAS_LIBCURL)
    endif()
    flexaids_configure_msvc_test(test_cofactor_blacklist)
    add_test(NAME CofactorBlacklistTests COMMAND test_cofactor_blacklist)

    # test_binding_mode_vibrational — Phase 3 ENCoM vibrational entropy integration
    add_executable(test_binding_mode_vibrational
        tests/test_binding_mode_vibrational.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/geometry.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
        LIB/BindingMode.cpp
        LIB/fast_optics.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_binding_mode_vibrational PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_vibrational PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_binding_mode_vibrational PRIVATE /O2)
    else()
        target_compile_options(test_binding_mode_vibrational PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_binding_mode_vibrational PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_binding_mode_vibrational PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_binding_mode_vibrational PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_binding_mode_vibrational)
    add_test(NAME BindingModeVibrationalTests COMMAND test_binding_mode_vibrational)

    # test_fast_optics — FastOPTICS clustering algorithm
    add_executable(test_fast_optics
        tests/test_fast_optics.cpp
        LIB/fast_optics.cpp
    )
    target_include_directories(test_fast_optics PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_fast_optics PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_fast_optics PRIVATE /O2)
    else()
        target_compile_options(test_fast_optics PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_fast_optics)
    add_test(NAME FastOPTICSTests COMMAND test_fast_optics)

    # test_tencom_entropy_diff — tENCoM entropy differential
    add_executable(test_tencom_entropy_diff
        tests/test_tencom_entropy_diff.cpp
        LIB/tENCoM/tencom_diff.cpp
        LIB/tENCoM/pdb_calpha.cpp
        LIB/tENCoM/tencm.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/encom.cpp
    )
    target_include_directories(test_tencom_entropy_diff PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM)
    target_link_libraries(test_tencom_entropy_diff PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_tencom_entropy_diff PRIVATE /O2)
    else()
        target_compile_options(test_tencom_entropy_diff PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_tencom_entropy_diff)
    target_link_libraries(test_tencom_entropy_diff PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_tencom_entropy_diff PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_tencom_entropy_diff)
    add_test(NAME TencomEntropyDiffTests COMMAND test_tencom_entropy_diff)

    # ── TargetValidation tests ──────────────────────────────────────────
    add_executable(test_target_validation
        tests/test_target_validation.cpp
        LIB/TargetValidation.cpp
    )
    target_include_directories(test_target_validation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_target_validation PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_target_validation PRIVATE /O2)
    else()
        target_compile_options(test_target_validation PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_target_validation)
    add_test(NAME TargetValidationTests COMMAND test_target_validation)

    # ── GrandPartitionFunction tests ──────────────────────────────────────
    add_executable(test_grand_partition
        tests/test_grand_partition.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/statmech.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_grand_partition PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_grand_partition PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_grand_partition PRIVATE /O2)
    else()
        target_compile_options(test_grand_partition PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_grand_partition)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_grand_partition PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_grand_partition PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_grand_partition PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_grand_partition)
    add_test(NAME GrandPartitionTests COMMAND test_grand_partition)

    # ── TargetServer integration tests ────────────────────────────────────
    add_executable(test_target_server
        tests/test_target_server.cpp
        LIB/TargetServer.cpp
        LIB/TargetValidation.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/TargetKnowledgeBase.cpp
        LIB/statmech.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_target_server PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_target_server PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_target_server PRIVATE /O2)
    else()
        target_compile_options(test_target_server PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_target_server)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_target_server PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_target_server PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_target_server PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_target_server)
    add_test(NAME TargetServerTests COMMAND test_target_server)

    # ── MultiSiteGPF tests ──────────────────────────────────────────────
    add_executable(test_multi_site_gpf
        tests/test_multi_site_gpf.cpp
        LIB/MultiSiteGPF.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/statmech.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_multi_site_gpf PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_multi_site_gpf PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_multi_site_gpf PRIVATE /O2)
    else()
        target_compile_options(test_multi_site_gpf PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_multi_site_gpf)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_multi_site_gpf PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_multi_site_gpf PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_multi_site_gpf PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_multi_site_gpf)
    add_test(NAME MultiSiteGPFTests COMMAND test_multi_site_gpf)

    # test_gist_hbond — GIST evaluator + directional H-bond scoring
    add_executable(test_gist_hbond tests/test_gist_hbond.cpp
                   LIB/GISTEvaluator.cpp LIB/HBondEvaluator.cpp)
    target_include_directories(test_gist_hbond PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_gist_hbond PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_gist_hbond PRIVATE /O2)
    else()
        target_compile_options(test_gist_hbond PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_gist_hbond)
    add_test(NAME GISTHBondTests COMMAND test_gist_hbond)

    # test_buffer_safety — bounded legacy config parsing regressions
    add_executable(test_buffer_safety
        tests/test_buffer_safety.cpp
    )
    target_include_directories(test_buffer_safety PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_buffer_safety PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_buffer_safety PRIVATE /O2)
    else()
        target_compile_options(test_buffer_safety PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_buffer_safety)
    add_test(NAME BufferSafetyTests COMMAND test_buffer_safety)

    message(STATUS "Unit tests enabled — build targets: test_statmech, test_hardware_dispatch, test_tencom_diff, test_target_validation, test_grand_partition, test_target_server, test_multi_site_gpf")
    # test_hbond_potential — Angular-dependent H-bond potential unit tests
    add_executable(test_hbond_potential
        tests/test_hbond_potential.cpp
    )
    target_include_directories(test_hbond_potential PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_hbond_potential PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_hbond_potential PRIVATE /O2)
    else()
        target_compile_options(test_hbond_potential PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_hbond_potential)
    add_test(NAME HBondPotentialTests COMMAND test_hbond_potential)

    # test_gist_grid — GIST desolvation grid unit tests
    add_executable(test_gist_grid
        tests/test_gist_grid.cpp
    )
    target_include_directories(test_gist_grid PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_gist_grid PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_gist_grid PRIVATE /O2)
    else()
        target_compile_options(test_gist_grid PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_gist_grid)
    add_test(NAME GISTGridTests COMMAND test_gist_grid)

    # test_ga_diversity — GA population diversity monitoring unit tests
    add_executable(test_ga_diversity
        tests/test_ga_diversity.cpp
    )
    target_include_directories(test_ga_diversity PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_ga_diversity PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ga_diversity PRIVATE /O2)
    else()
        target_compile_options(test_ga_diversity PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_ga_diversity)
    add_test(NAME GADiversityTests COMMAND test_ga_diversity)

    # test_soa_conversion — AoS↔SoA atom data conversion unit tests
    add_executable(test_soa_conversion
        tests/test_soa_conversion.cpp
    )
    target_include_directories(test_soa_conversion PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_soa_conversion PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_soa_conversion PRIVATE /O2)
    else()
        target_compile_options(test_soa_conversion PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_soa_conversion)
    flexaids_configure_msvc_test(test_soa_conversion)
    add_test(NAME SoAConversionTests COMMAND test_soa_conversion)

    # test_dift — Discrete Fourier Transform torsional parametrization
    add_executable(test_dift
        tests/test_dift.cpp
        LIB/DiFT/DiFT.cpp
    )
    target_include_directories(test_dift PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/DiFT
    )
    target_link_libraries(test_dift PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_dift PRIVATE /Od /Zi)
    else()
        target_compile_options(test_dift PRIVATE -O0 -g3)
    endif()
    flexaids_configure_msvc_test(test_dift)
    add_test(NAME DiFTTests COMMAND test_dift)

    # test_distributed_backend — Distributed backend abstraction unit tests
    add_executable(test_distributed_backend
        tests/test_distributed_backend.cpp
    )
    target_include_directories(test_distributed_backend PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_distributed_backend PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_distributed_backend PRIVATE /O2)
    else()
        target_compile_options(test_distributed_backend PRIVATE -O2)
    endif()
    flexaids_configure_msvc_test(test_distributed_backend)
    add_test(NAME DistributedBackendTests COMMAND test_distributed_backend)

    # test_ga_context — GAContext re-entrant state + GPUContextPool singleton
    add_executable(test_ga_context
        tests/test_ga_context.cpp
        LIB/GAContext.cpp
    )
    target_include_directories(test_ga_context PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
    )
    target_link_libraries(test_ga_context PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ga_context PRIVATE /O2)
    else()
        target_compile_options(test_ga_context PRIVATE -O2)
    endif()
    target_link_libraries(test_ga_context PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_ga_context PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_ga_context)
    add_test(NAME GAContextTests COMMAND test_ga_context)

    if(NOT TARGET test_natural)
        add_executable(test_natural
            tests/test_natural.cpp
            tests/stubs.cpp
            LIB/NATURaL/RibosomeElongation.cpp
            LIB/NATURaL/TransloconInsertion.cpp
            LIB/NATURaL/NucleationDetector.cpp
            LIB/NATURaL/NATURaLDualAssembly.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
            LIB/Vcontacts.cpp
            LIB/geometry.cpp
        )
        target_include_directories(test_natural PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        )
        target_link_libraries(test_natural PRIVATE GTest::gtest GTest::gtest_main)
        if(MSVC)
            target_compile_options(test_natural PRIVATE /O2)
            target_compile_definitions(test_natural PRIVATE
                _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
        else()
            target_compile_options(test_natural PRIVATE -O2)
        endif()
        if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
            target_link_libraries(test_natural PRIVATE OpenMP::OpenMP_CXX)
        endif()
        target_link_libraries(test_natural PRIVATE Eigen3::Eigen)
        target_compile_definitions(test_natural PRIVATE FLEXAIDS_HAS_EIGEN)
        flexaids_configure_simd(test_natural)
        flexaids_configure_msvc_test(test_natural)
        add_test(NAME NATURaLTests COMMAND test_natural)
    endif()

    # test_nascent_chain_scheduler — cotranslational DualAssembly suite
    add_executable(test_nascent_chain_scheduler
        tests/test_nascent_chain_scheduler.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/tENCoM/tencm.cpp
        LIB/encom.cpp
    )
    target_include_directories(test_nascent_chain_scheduler PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
    )
    target_link_libraries(test_nascent_chain_scheduler PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_nascent_chain_scheduler PRIVATE /O2)
        target_compile_definitions(test_nascent_chain_scheduler PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(test_nascent_chain_scheduler PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_nascent_chain_scheduler PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_nascent_chain_scheduler PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_nascent_chain_scheduler PRIVATE
        FLEXAIDS_HAS_EIGEN
        FLEXAIDS_OMIT_DUAL_ASSEMBLY_ENGINE)
    flexaids_configure_simd(test_nascent_chain_scheduler)
    flexaids_configure_msvc_test(test_nascent_chain_scheduler)
    add_test(NAME NascentChainSchedulerTests COMMAND test_nascent_chain_scheduler)

    if(TARGET dual_assembly)
        add_test(NAME DualAssemblyRealGAFailsClosed
            COMMAND dual_assembly
                --real-ga
                --target-pdb fake_protofibril.pdb
                --sequence AAAAA
                --no-sim-c)
        set_tests_properties(DualAssemblyRealGAFailsClosed PROPERTIES
            WILL_FAIL TRUE)
    endif()

    message(STATUS "Unit tests enabled — build targets: test_statmech, test_hardware_dispatch, test_tencom_diff, test_target_validation, test_grand_partition, test_target_server, test_ga_context, test_nascent_chain_scheduler")
endif()

# ─── Swift Bridge (macOS only) ────────────────────────────────────────────
option(BUILD_SWIFT_BRIDGE "Build Swift bridge static library (macOS only)" OFF)
if(BUILD_SWIFT_BRIDGE)
    if(NOT APPLE)
        message(FATAL_ERROR "BUILD_SWIFT_BRIDGE requires macOS (Apple platform)")
    endif()

    enable_language(OBJCXX)

    find_library(FOUNDATION_LIBRARY Foundation REQUIRED)

    add_library(FlexAIDCore STATIC
        swift/Sources/FlexAIDCore/FXStatMechEngine.mm
        swift/Sources/FlexAIDCore/FXENCoMEngine.mm
        swift/Sources/FlexAIDCore/FXBoltzmannLUT.mm
        swift/Sources/FlexAIDCore/FXGA.mm
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/encom.cpp
    )
    target_include_directories(FlexAIDCore PUBLIC
        swift/Sources/FlexAIDCore/include
        LIB
    )
    set_target_properties(FlexAIDCore PROPERTIES
        CXX_STANDARD ${CMAKE_CXX_STANDARD}   # follows global (capped to 20 on MSVC)
        CXX_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE ON
    )
    set_source_files_properties(
        swift/Sources/FlexAIDCore/FXStatMechEngine.mm
        swift/Sources/FlexAIDCore/FXENCoMEngine.mm
        swift/Sources/FlexAIDCore/FXBoltzmannLUT.mm
        swift/Sources/FlexAIDCore/FXGA.mm
        PROPERTIES LANGUAGE OBJCXX
    )
    target_compile_options(FlexAIDCore PRIVATE -std=c++26 -O2)
    target_link_libraries(FlexAIDCore PUBLIC ${FOUNDATION_LIBRARY})

    message(STATUS "Swift bridge enabled — build target: FlexAIDCore (static library)")

    # test_unified_dispatch — runtime HardwareDispatcher tests
    add_executable(test_unified_dispatch
        tests/test_unified_dispatch.cpp
        LIB/UnifiedHardwareDispatch.cpp
    )
    target_include_directories(test_unified_dispatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_unified_dispatch PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_unified_dispatch PRIVATE /O2)
    else()
        target_compile_options(test_unified_dispatch PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_unified_dispatch)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_unified_dispatch PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_unified_dispatch PRIVATE Eigen3::Eigen)

    target_compile_definitions(test_unified_dispatch PRIVATE FLEXAIDS_HAS_EIGEN)
    add_test(NAME UnifiedDispatchTests COMMAND test_unified_dispatch)

    # test_coarse_screen — NRGRank CoarseScreen + TwoStageScreen unit tests
    add_executable(test_coarse_screen
        tests/test_coarse_screen.cpp
        LIB/CoarseScreen.cpp
        LIB/TwoStageScreen.cpp
    )
    target_include_directories(test_coarse_screen PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_coarse_screen PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_coarse_screen PRIVATE /O2)
    else()
        target_compile_options(test_coarse_screen PRIVATE -O2)
    endif()
    flexaids_configure_simd(test_coarse_screen)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_coarse_screen PRIVATE OpenMP::OpenMP_CXX)
    endif()
    flexaids_configure_msvc_test(test_coarse_screen)
    add_test(NAME CoarseScreenTests COMMAND test_coarse_screen)

    # test_process_ligand — ProcessLigand pipeline: SmilesParser, RingPerception,
    # Aromaticity, RotatableBonds, ValenceChecker, SybylTyper
    add_executable(test_process_ligand
        tests/test_process_ligand.cpp
        LIB/ProcessLigand/SmilesParser.cpp
        LIB/ProcessLigand/RingPerception.cpp
        LIB/ProcessLigand/Aromaticity.cpp
        LIB/ProcessLigand/RotatableBonds.cpp
        LIB/ProcessLigand/ValenceChecker.cpp
        LIB/ProcessLigand/SybylTyper.cpp
        LIB/ProcessLigand/CoordBuilder.cpp
        LIB/ProcessLigand/FlexAIDWriter.cpp
        LIB/ProcessLigand/ProcessLigand.cpp
    )
    target_include_directories(test_process_ligand PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ProcessLigand
    )
    target_link_libraries(test_process_ligand PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_process_ligand PRIVATE /Od /Zi)
        target_compile_definitions(test_process_ligand PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(test_process_ligand PRIVATE -O0 -g3)
    endif()
    target_link_libraries(test_process_ligand PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_process_ligand PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_process_ligand)
    add_test(NAME ProcessLigandTests COMMAND test_process_ligand)

    # test_ligand_ring_flex — LigandRingFlex unified ring flexibility interface
    add_executable(test_ligand_ring_flex
        tests/test_ligand_ring_flex.cpp
        LIB/LigandRingFlex/LigandRingFlex.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
    )
    target_include_directories(test_ligand_ring_flex PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
    )
    target_link_libraries(test_ligand_ring_flex PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ligand_ring_flex PRIVATE /Od /Zi)
    else()
        target_compile_options(test_ligand_ring_flex PRIVATE -O0 -g3)
    endif()
    target_link_libraries(test_ligand_ring_flex PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_ligand_ring_flex PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_ligand_ring_flex)
    add_test(NAME LigandRingFlexTests COMMAND test_ligand_ring_flex)

    # test_natural — RibosomeElongation, TransloconInsertion, NucleationSeedDetector
    if(NOT TARGET test_natural)
        add_executable(test_natural
            tests/test_natural.cpp
            LIB/NATURaL/RibosomeElongation.cpp
            LIB/NATURaL/TransloconInsertion.cpp
            LIB/NATURaL/NucleationDetector.cpp
            LIB/NATURaL/NATURaLDualAssembly.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
            LIB/Vcontacts.cpp
            LIB/geometry.cpp
        )
        target_include_directories(test_natural PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        )
        target_link_libraries(test_natural PRIVATE GTest::gtest GTest::gtest_main)
        if(MSVC)
            target_compile_options(test_natural PRIVATE /O2)
            target_compile_definitions(test_natural PRIVATE
                _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
        else()
            target_compile_options(test_natural PRIVATE -O2)
        endif()
        if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
            target_link_libraries(test_natural PRIVATE OpenMP::OpenMP_CXX)
        endif()
        target_link_libraries(test_natural PRIVATE Eigen3::Eigen)
        target_compile_definitions(test_natural PRIVATE FLEXAIDS_HAS_EIGEN)
        flexaids_configure_simd(test_natural)
        flexaids_configure_msvc_test(test_natural)
        add_test(NAME NATURaLTests COMMAND test_natural)
    endif()

    # test_docking_pipeline — end-to-end integration: ProcessLigand → StatMech → BindingMode
    add_executable(test_docking_pipeline
        tests/test_docking_pipeline.cpp
        tests/stubs.cpp
        LIB/ProcessLigand/SmilesParser.cpp
        LIB/ProcessLigand/RingPerception.cpp
        LIB/ProcessLigand/Aromaticity.cpp
        LIB/ProcessLigand/RotatableBonds.cpp
        LIB/ProcessLigand/ValenceChecker.cpp
        LIB/ProcessLigand/SybylTyper.cpp
        LIB/ProcessLigand/CoordBuilder.cpp
        LIB/ProcessLigand/FlexAIDWriter.cpp
        LIB/ProcessLigand/ProcessLigand.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp
        LIB/gaboom.cpp
        LIB/VibEntropy.cpp
        LIB/InStreamClustering.cpp
        LIB/GAContext.cpp
        LIB/Vcontacts.cpp
        LIB/BindingMode.cpp
        LIB/encom.cpp
        LIB/geometry.cpp
        LIB/fast_optics.cpp
        LIB/tENCoM/tencm.cpp
        LIB/ShannonThermoStack/ShannonThermoStack.cpp
        LIB/CleftDetector.cpp
        LIB/NATURaL/NATURaLDualAssembly.cpp
        LIB/NATURaL/RibosomeElongation.cpp
        LIB/NATURaL/TransloconInsertion.cpp
        LIB/NATURaL/NucleationDetector.cpp
        LIB/NATURaL/NascentChainScheduler.cpp
        LIB/NATURaL/FibrilGrowthOracle.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/NATURaL/DualAssemblyRunner.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_docking_pipeline PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ProcessLigand
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_docking_pipeline PRIVATE GTest::gtest GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_docking_pipeline PRIVATE /O2)
        target_compile_definitions(test_docking_pipeline PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(test_docking_pipeline PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_docking_pipeline PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_docking_pipeline PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_docking_pipeline PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_simd(test_docking_pipeline)
    flexaids_configure_msvc_test(test_docking_pipeline)
    add_test(NAME DockingPipelineTests COMMAND test_docking_pipeline)
    set_tests_properties(DockingPipelineTests PROPERTIES LABELS "integration")

endif()

# =============================================================================
# ── D2: Profile-Guided Optimization (PGO) ────────────────────────────────────
# Two-stage workflow:
#   Stage 1: cmake -DFLEXAIDS_PGO_GENERATE=ON ... && cmake --build && <training run>
#   Stage 2: cmake -DFLEXAIDS_PGO_USE=/path/to/default.profdata ... && cmake --build
# Clang only (Apple Silicon default). GCC uses -fprofile-arcs/-fprofile-use.
# =============================================================================
option(FLEXAIDS_PGO_GENERATE "Stage 1: instrument binary for PGO profiling" OFF)
set(FLEXAIDS_PGO_USE "" CACHE STRING "Stage 2: path to .profdata (clang) or dir (gcc) for PGO use")
if(FLEXAIDS_PGO_GENERATE)
    message(STATUS "[PGO] Stage 1: building instrumented binary (FLEXAIDS_PGO_GENERATE=ON)")
    target_compile_options(FlexAIDdS PRIVATE -fprofile-generate)
    target_link_options(FlexAIDdS PRIVATE -fprofile-generate)
elseif(FLEXAIDS_PGO_USE)
    message(STATUS "[PGO] Stage 2: building with profile data: ${FLEXAIDS_PGO_USE}")
    target_compile_options(FlexAIDdS PRIVATE -fprofile-use=${FLEXAIDS_PGO_USE} -fprofile-correction)
    target_link_options(FlexAIDdS PRIVATE -fprofile-use=${FLEXAIDS_PGO_USE})
endif()

# FlexAID source validator guard (MUST be the last thing in this file)
# =============================================================================
#
# This runs after every add_subdirectory / target_sources / add_test etc.
# so it sees the complete declared source set.
#
# P0 integration: Supports FLEXAIDS_STRICT_SOURCE_VALIDATION for CI.
# Developers get the safe WARN_ONLY default.
# =============================================================================
include(cmake/ValidateSources.cmake)

if(FLEXAIDS_STRICT_SOURCE_VALIDATION)
    flexaids_validate_sources(STRICT)
else()
    flexaids_validate_sources(WARN_ONLY)
endif()
