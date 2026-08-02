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

# ─── Centralized options + platform + SIMD helpers ─────────────────────────
# Extracted to cmake/FlexAIDOptions.cmake to fight monolith growth.
# All FLEXAIDS_* toggles, architecture detection, flexaids_configure_simd(),
# and coverage global flags now live in one place.
include(cmake/FlexAIDOptions.cmake)

# External dependencies (Eigen, OpenMP, CUDA, HIP, Metal, MPI, pybind11, etc.)
# Extracted to fight CMake monolith. Must precede add_subdirectory(LIB).
include(cmake/FlexAIDDependencies.cmake)

if(APPLE AND FLEXAIDS_USE_METAL)
    if(NOT FLEXAIDS_METALC)
        execute_process(
            COMMAND xcrun -f metal
            OUTPUT_VARIABLE _FLEXAIDS_XCRUN_METALC
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE _FLEXAIDS_XCRUN_METALC_RESULT
        )
        if(_FLEXAIDS_XCRUN_METALC_RESULT EQUAL 0 AND _FLEXAIDS_XCRUN_METALC)
            set(FLEXAIDS_METALC "${_FLEXAIDS_XCRUN_METALC}" CACHE FILEPATH "Metal shader compiler")
        endif()
    endif()

    if(FLEXAIDS_METALC)
        get_filename_component(_FLEXAIDS_METAL_TOOL_DIR "${FLEXAIDS_METALC}" DIRECTORY)
    endif()

    if(NOT FLEXAIDS_METALLIB AND _FLEXAIDS_METAL_TOOL_DIR)
        find_program(FLEXAIDS_METALLIB
            NAMES metallib
            HINTS "${_FLEXAIDS_METAL_TOOL_DIR}"
            NO_DEFAULT_PATH
        )
    endif()

    if(NOT FLEXAIDS_METALC OR NOT EXISTS "${FLEXAIDS_METALC}")
        message(FATAL_ERROR "FLEXAIDS_USE_METAL=ON but no executable Metal compiler was found. Run: xcodebuild -downloadComponent MetalToolchain")
    endif()
    if(NOT FLEXAIDS_METALLIB OR NOT EXISTS "${FLEXAIDS_METALLIB}")
        message(FATAL_ERROR "FLEXAIDS_USE_METAL=ON but metallib was not found beside ${FLEXAIDS_METALC}. Run: xcodebuild -downloadComponent MetalToolchain")
    endif()

    set(FLEXAIDS_METAL_MODULE_CACHE
        "${CMAKE_BINARY_DIR}/metal_module_cache"
        CACHE PATH "Run-local Metal compiler module cache"
    )
    file(MAKE_DIRECTORY "${FLEXAIDS_METAL_MODULE_CACHE}")
    message(STATUS "Metal compiler: ${FLEXAIDS_METALC}")
    message(STATUS "Metal library linker: ${FLEXAIDS_METALLIB}")
    message(STATUS "Metal module cache: ${FLEXAIDS_METAL_MODULE_CACHE}")
endif()

# Optional WebGPU (Dawn / wgpu-native) backend detection — additive, skips
# gracefully when the dependency isn't present on this system.
include(cmake/FlexAIDWebGPU.cmake)

# ─── PoseBusters CLI auto-detection (optional) ─────────────────────────────
# The upstream PoseBusters `bust` CLI is an OPTIONAL validation bridge used by
# LIB/PoseBust/BustCli.cpp. Everything builds and runs without it; when it is
# absent the bridge simply reports backend "bust_cli_missing" and the native
# Apache-2.0 NativePoseQC path is used instead.
#
# At runtime the binary resolves `bust` in this order (BustCli.cpp):
#   1. $FLEXAIDDS_POSEBUSTERS_BIN
#   2. PATH
#   3. $FLEXAIDDS_ROOT/.venv-posebusters/bin/bust
#   4. ./.venv-posebusters/bin/bust or ../.venv-posebusters/bin/bust
#   5. the compile-time default baked in below, if configure found one
#
# Manual install (if the detection below reports NOT found):
#   python3 -m venv .venv-posebusters
#   .venv-posebusters/bin/pip install posebusters
# then either re-run cmake, or export the path explicitly:
#   export FLEXAIDDS_POSEBUSTERS_BIN=/abs/path/to/bust
find_program(POSEBUSTERS_BIN
    NAMES bust
    HINTS "${CMAKE_CURRENT_SOURCE_DIR}/.venv-posebusters/bin"
    PATHS ENV PATH
)
if(POSEBUSTERS_BIN)
    message(STATUS "PoseBusters CLI found: ${POSEBUSTERS_BIN}")
    add_compile_definitions(FLEXAIDDS_POSEBUSTERS_BIN_DEFAULT="${POSEBUSTERS_BIN}")
else()
    message(STATUS "PoseBusters CLI not found (optional) — upstream bust bridge "
                   "will report 'bust_cli_missing'. To enable: "
                   "python3 -m venv .venv-posebusters && "
                   ".venv-posebusters/bin/pip install posebusters, then re-run cmake.")
endif()

# ─── Core library (LIB/) ────────────────────────────────────────────────
add_subdirectory(LIB)

# Component targets (FlexAID::Core, ::VoronoiCF, ::StatMech, etc.)
# See cmake/FlexAIDComponents.cmake for the audit-driven split direction.
include(cmake/FlexAIDComponents.cmake)

# Professional helpers (test boilerplate reduction, common patterns)
include(cmake/FlexAIDHelpers.cmake)

# ─── executable ─────────────────────────────────────────────────────────
if(FLEXAIDS_BUILD_CORE)
add_executable(FlexAID "")
# Slimmed: use central OBJECT library instead of duplicating ~150 source files.
# Only the dedicated entry points + the core objects.
target_sources(FlexAID PRIVATE LIB/top.cpp)
target_link_libraries(FlexAID PRIVATE flexaid_core)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(FlexAID PRIVATE m)
endif()

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

# ─── P0 build-level wins: lift the FlexAIDdS LTO/native/NDEBUG block onto FlexAID ─
# The classic FlexAID target historically shipped with only -O3 -ffast-math (above)
# while the FlexAIDdS target carried an aggressive LTO + INTERPROCEDURAL_OPTIMIZATION
# + -march=native/-mcpu=native + -DNDEBUG + strip block. flexaids_apply_fast_optimization
# lifts that same block onto FlexAID, guarded by BUILD_FLEXAID_FAST (default ON, but
# force-disabled under coverage/sanitizer builds by the SAME guard used for
# BUILD_FLEXAIDDS_FAST — see cmake/FlexAIDOptions.cmake). FLEXAID_PGO adds optional
# profile-generate/profile-use instrumentation. Both are no-ops when disabled.
flexaids_apply_fast_optimization(FlexAID)
flexaids_apply_pgo(FlexAID)

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

# ─── Metal sources (attach to flexaid_core once) ─────────────────────────
# hardware_detect.cpp, UnifiedHardwareDispatch.cpp, gaboom.cpp, FOPTICS.cpp and
# other flexaid_core TUs reference metal_eval / metal_rmsd when FLEXAIDS_USE_METAL
# is on. Implementations MUST live on flexaid_core (or a library it pulls in),
# not only on individual executables. Attaching .mm files solely to FlexAID /
# FlexAIDdS left undefined symbols under LTO / multi-target links (Homebrew
# --with-metal). All consumers of flexaid_core inherit the OBJCXX objects and
# Metal frameworks via PUBLIC linkage below.
if(FLEXAIDS_USE_METAL)
    set(_FLEXAIDS_METAL_BRIDGES
        LIB/metal_eval.mm
        LIB/MetalRMSDBridge.mm
        LIB/ShannonThermoStack/ShannonMetalBridge.mm
        LIB/tENCoM/tencm_metal.mm
        LIB/gpu_fast_optics_metal.mm
        LIB/CavityDetect/CavityDetectMetalBridge.mm
        LIB/TurboQuantMetalBridge.mm
    )
    target_sources(flexaid_core PRIVATE ${_FLEXAIDS_METAL_BRIDGES})
    target_compile_options(flexaid_core PRIVATE -fno-objc-arc)
    target_link_libraries(flexaid_core PUBLIC
        ${METAL_LIBRARY}
        ${FOUNDATION_LIBRARY}
        ${METALKIT_LIBRARY}
    )
    set_source_files_properties(${_FLEXAIDS_METAL_BRIDGES} PROPERTIES
        LANGUAGE OBJCXX
        COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    target_compile_definitions(flexaid_core PUBLIC
        FLEXAIDS_USE_METAL
        FLEXAIDS_HAS_METAL_SHANNON
        FLEXAIDS_HAS_METAL_TENCM
    )
endif()

# ─── WebGPU sources (attach to flexaid_core once, additive/optional) ───────
# Sibling of the Metal block above: same pattern (opt-in cache option,
# graceful skip when the dependency is absent — see cmake/FlexAIDWebGPU.cmake),
# targeting Dawn or wgpu-native via the vendor-neutral webgpu.h API so the
# same src/backends/webgpu/webgpu_eval.cpp compiles against either.
if(FLEXAIDS_USE_WEBGPU)
    target_sources(flexaid_core PRIVATE src/backends/webgpu/webgpu_eval.cpp)
    target_include_directories(flexaid_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/backends/webgpu
        ${FLEXAIDS_WEBGPU_INCLUDE_DIRS}
    )
    target_link_libraries(flexaid_core PUBLIC ${FLEXAIDS_WEBGPU_LIBRARIES})
    target_compile_definitions(flexaid_core PUBLIC FLEXAIDS_USE_WEBGPU)

    # Stage the WGSL shader next to the FlexAID/FlexAIDdS binaries at build
    # time, mirroring the runtime-data staging pattern used for
    # MC_st0r5.2_6.dat etc. — webgpu_eval.cpp loads it via
    # FLEXAIDDS_WEBGPU_SHADER_DIR or the binary's own working directory.
    add_custom_target(flexaid_webgpu_shaders ALL
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/webgpu_shaders"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/src/backends/webgpu/cf_eval.wgsl"
            "${CMAKE_BINARY_DIR}/webgpu_shaders/cf_eval.wgsl"
        COMMENT "Staging WebGPU CF shader (cf_eval.wgsl)"
    )
    add_dependencies(flexaid_core flexaid_webgpu_shaders)
endif()

if(FLEXAIDS_USE_METAL)
    # ── Compile .metal shaders → .metallib (single shared set) ────────────
    set(CAVITY_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect/CavityDetect.metal")
    set(CAVITY_AIR        "${CMAKE_CURRENT_BINARY_DIR}/CavityDetect.air")
    set(CAVITY_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/CavityDetect.metallib")
    add_custom_command(
        OUTPUT  "${CAVITY_AIR}"
        COMMAND "${FLEXAIDS_METALC}" "-fmodules-cache-path=${FLEXAIDS_METAL_MODULE_CACHE}" -c "${CAVITY_METAL_SRC}" -o "${CAVITY_AIR}"
        DEPENDS "${CAVITY_METAL_SRC}"
        COMMENT "Compiling CavityDetect Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${CAVITY_METALLIB}"
        COMMAND "${FLEXAIDS_METALLIB}" "${CAVITY_AIR}" -o "${CAVITY_METALLIB}"
        DEPENDS "${CAVITY_AIR}"
        COMMENT "Linking CavityDetect.metallib"
        VERBATIM
    )
    add_custom_target(CavityDetectMetal DEPENDS "${CAVITY_METALLIB}")

    set(TQ_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/TurboQuant.metal")
    set(TQ_AIR        "${CMAKE_CURRENT_BINARY_DIR}/TurboQuant.air")
    set(TQ_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/TurboQuant.metallib")
    add_custom_command(
        OUTPUT  "${TQ_AIR}"
        COMMAND "${FLEXAIDS_METALC}" "-fmodules-cache-path=${FLEXAIDS_METAL_MODULE_CACHE}" -c "${TQ_METAL_SRC}" -o "${TQ_AIR}"
        DEPENDS "${TQ_METAL_SRC}"
        COMMENT "Compiling TurboQuant Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${TQ_METALLIB}"
        COMMAND "${FLEXAIDS_METALLIB}" "${TQ_AIR}" -o "${TQ_METALLIB}"
        DEPENDS "${TQ_AIR}"
        COMMENT "Linking TurboQuant.metallib"
        VERBATIM
    )
    add_custom_target(TurboQuantMetal DEPENDS "${TQ_METALLIB}")

    set(RMSD_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/MetalRMSD.metal")
    set(RMSD_AIR        "${CMAKE_CURRENT_BINARY_DIR}/MetalRMSD.air")
    set(RMSD_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/MetalRMSD.metallib")
    add_custom_command(
        OUTPUT  "${RMSD_AIR}"
        COMMAND "${FLEXAIDS_METALC}" "-fmodules-cache-path=${FLEXAIDS_METAL_MODULE_CACHE}" -c "${RMSD_METAL_SRC}" -o "${RMSD_AIR}"
        DEPENDS "${RMSD_METAL_SRC}"
        COMMENT "Compiling MetalRMSD Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${RMSD_METALLIB}"
        COMMAND "${FLEXAIDS_METALLIB}" "${RMSD_AIR}" -o "${RMSD_METALLIB}"
        DEPENDS "${RMSD_AIR}"
        COMMENT "Linking MetalRMSD.metallib"
        VERBATIM
    )
    add_custom_target(MetalRMSDMetal DEPENDS "${RMSD_METALLIB}")

    set(SHANNON_METAL_SRC  "${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack/shannon_metal.metal")
    set(SHANNON_AIR        "${CMAKE_CURRENT_BINARY_DIR}/ShannonEntropy.air")
    set(SHANNON_METALLIB   "${CMAKE_CURRENT_BINARY_DIR}/ShannonEntropy.metallib")
    add_custom_command(
        OUTPUT  "${SHANNON_AIR}"
        COMMAND "${FLEXAIDS_METALC}" "-fmodules-cache-path=${FLEXAIDS_METAL_MODULE_CACHE}" -c "${SHANNON_METAL_SRC}" -o "${SHANNON_AIR}"
        DEPENDS "${SHANNON_METAL_SRC}"
        COMMENT "Compiling ShannonEntropy Metal shader"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  "${SHANNON_METALLIB}"
        COMMAND "${FLEXAIDS_METALLIB}" "${SHANNON_AIR}" -o "${SHANNON_METALLIB}"
        DEPENDS "${SHANNON_AIR}"
        COMMENT "Linking ShannonEntropy.metallib"
        VERBATIM
    )
    add_custom_target(ShannonEntropyMetal DEPENDS "${SHANNON_METALLIB}")

    # Path macros for bridges that load metallibs at runtime.
    # PUBLIC so every TU that includes the bridges sees the same paths.
    target_compile_definitions(flexaid_core PUBLIC
        CAVITY_METALLIB_PATH="${CAVITY_METALLIB}"
        TURBOQUANT_METALLIB_PATH="${TQ_METALLIB}"
        METALRMSD_METALLIB_PATH="${RMSD_METALLIB}"
        SHANNON_METALLIB_PATH="${SHANNON_METALLIB}"
    )
    add_dependencies(flexaid_core
        CavityDetectMetal
        TurboQuantMetal
        MetalRMSDMetal
        ShannonEntropyMetal
    )
    # Convenience aliases for older custom-target names referenced by helpers.
    if(NOT TARGET CavityDetectMetalDS)
        add_custom_target(CavityDetectMetalDS DEPENDS CavityDetectMetal)
    endif()
    if(NOT TARGET MetalRMSDMetalDS)
        add_custom_target(MetalRMSDMetalDS DEPENDS MetalRMSDMetal)
    endif()
    if(NOT TARGET TurboQuantMetalDS)
        add_custom_target(TurboQuantMetalDS DEPENDS TurboQuantMetal)
    endif()
    # Back-compat path vars used by benchmark_datasets / attach helpers.
    set(CAVITY_METALLIB_DS "${CAVITY_METALLIB}")
    set(RMSD_METALLIB_DS "${RMSD_METALLIB}")
    set(TQ_METALLIB_DS "${TQ_METALLIB}")
    message(STATUS "Metal bridges attached to flexaid_core (all executables inherit them)")
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

    add_executable(ligand_tencom_pose
        LIB/ligand_tencom_pose.cpp
        LIB/tENCoM/tencm.cpp
        LIB/VibEntropy.cpp
        LIB/PoseBust/Loaders.cpp
        LIB/TurboQuant.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(ligand_tencom_pose PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PoseBust)
    target_link_libraries(ligand_tencom_pose PRIVATE Eigen3::Eigen)
    target_compile_definitions(ligand_tencom_pose PRIVATE FLEXAIDS_HAS_EIGEN)
    if(MSVC)
        target_compile_options(ligand_tencom_pose PRIVATE /O2)
        target_compile_definitions(ligand_tencom_pose PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(ligand_tencom_pose PRIVATE -O3)
    endif()
    flexaids_configure_simd(ligand_tencom_pose)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(ligand_tencom_pose PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if(FLEXAIDS_USE_METAL)
        target_sources(ligand_tencom_pose PRIVATE LIB/tENCoM/tencm_metal.mm)
        target_compile_definitions(ligand_tencom_pose PRIVATE FLEXAIDS_HAS_METAL_TENCM)
        set_source_files_properties(LIB/tENCoM/tencm_metal.mm PROPERTIES
            LANGUAGE OBJCXX COMPILE_FLAGS "-std=c++26 -fno-objc-arc")
        target_link_libraries(ligand_tencom_pose PRIVATE
            ${METAL_LIBRARY} ${FOUNDATION_LIBRARY} ${METALKIT_LIBRARY})
    endif()
    message(STATUS "Ligand tENCoM/Eigen pose validator enabled — build target: ligand_tencom_pose")
endif()

# ─── P5: NRGRank two-stage screening driver (flexaid_screen) ───────────────────
# Standalone CLI that runs the rigid NRGRank coarse screen (CoarseScreen /
# TwoStageScreen) as a pre-filter in front of the expensive GA. Reuses the same
# classes exercised by test_coarse_screen; no coupling to the GA/top.cpp engine.
option(ENABLE_TWO_STAGE_SCREEN "Build the NRGRank two-stage screening CLI (flexaid_screen)" ON)
if(ENABLE_TWO_STAGE_SCREEN
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LIB/two_stage_screen_main.cpp"
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LIB/CoarseScreen.cpp"
   AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LIB/TwoStageScreen.cpp")
    add_executable(flexaid_screen
        LIB/two_stage_screen_main.cpp
        LIB/CoarseScreen.cpp
        LIB/TwoStageScreen.cpp
    )
    target_include_directories(flexaid_screen PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    set_target_properties(flexaid_screen PROPERTIES
        CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(flexaid_screen PRIVATE /O2)
        target_compile_definitions(flexaid_screen PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(flexaid_screen PRIVATE -O3)
    endif()
    flexaids_configure_simd(flexaid_screen)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(flexaid_screen PRIVATE OpenMP::OpenMP_CXX)
    endif()
    message(STATUS "Two-stage screening CLI enabled — build target: flexaid_screen")
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

# ─── ProcessLigand CLI tool ─────────────────────────────────────────────────
# Standalone ligand preprocessing CLI (SMILES/SDF/MOL2 -> FlexAID input files).
# Uses the ProcessLigand module (now wired as first-class target, no longer orphan).
option(ENABLE_PROCESS_LIGAND "Build the ProcessLigand standalone CLI tool (flexaidds-process-ligand)" OFF)
if(ENABLE_PROCESS_LIGAND)
    add_executable(flexaidds_process_ligand
        LIB/ProcessLigand/main.cpp
    )
    target_link_libraries(flexaidds_process_ligand PRIVATE flexaid_core)
    target_include_directories(flexaidds_process_ligand PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ProcessLigand
    )
    if(MSVC)
        target_compile_options(flexaidds_process_ligand PRIVATE /O2)
        target_compile_definitions(flexaidds_process_ligand PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(flexaidds_process_ligand PRIVATE -O3)
    endif()
    message(STATUS "ProcessLigand CLI tool enabled — build target: flexaidds_process_ligand")
endif()

# ─── FlexAIDdS: ultra-fast docking executable ─────────────────────────────
# P0 validator references for GPU/metal shaders (actual add in Metal block above)
#   LIB/gpu_fast_optics_metal.metal LIB/gpu_batch_distance_hip.hip etc.
# Uses flexaid_core OBJECT library (same sources as FlexAID) with aggressive
# LTO and native-arch flags for production deployments.
option(BUILD_FLEXAIDDS_FAST "Build ultra-fast FlexAIDdS docking executable (LTO + native)" ON)
# LTO is incompatible with gcov coverage and with ASan/TSan/UBSan instrumented builds.
if(FLEXAIDS_ENABLE_COVERAGE
   OR CMAKE_CXX_FLAGS MATCHES "-fsanitize"
   OR CMAKE_C_FLAGS MATCHES "-fsanitize"
   OR CMAKE_EXE_LINKER_FLAGS MATCHES "-fsanitize")
    if(BUILD_FLEXAIDDS_FAST)
        message(STATUS "Disabling BUILD_FLEXAIDDS_FAST (incompatible with coverage/sanitizer flags)")
    endif()
    set(BUILD_FLEXAIDDS_FAST OFF CACHE BOOL "Build ultra-fast FlexAIDdS docking executable (LTO + native)" FORCE)
endif()
if(BUILD_FLEXAIDDS_FAST)
    # --redock <PDBid> reuses DatasetRunner cognate prep (download + ligand SDF +
    # apo strip). Same supporting sources as test_dataset_runner / benchmark_datasets.
    set(FLEXAIDDS_REDOCK_SOURCES
        LIB/DatasetRunner.cpp
        LIB/DatasetRunnerStats.cpp
        LIB/DatasetRunnerProvenance.cpp
        LIB/AsyncPipeline.cpp
        LIB/PoseBust/Loaders.cpp
        LIB/PoseBust/ChecksChemistry.cpp
        LIB/PoseBust/ChecksGeometry.cpp
        LIB/PoseBust/ChecksProtein.cpp
        LIB/PoseBust/Engine.cpp
        LIB/PoseBust/BustCli.cpp
    )
    add_executable(FlexAIDdS LIB/top.cpp ${FLEXAIDDS_REDOCK_SOURCES})
    target_link_libraries(FlexAIDdS PRIVATE flexaid_core)
    target_include_directories(FlexAIDdS PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PoseBust
    )
    # Gate DatasetRunner-backed --redock path in top.cpp (shared with FlexAID).
    target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDDS_ENABLE_REDOCK=1)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(FlexAIDdS PRIVATE m)
    endif()
    if(CURL_FOUND)
        target_link_libraries(FlexAIDdS PRIVATE CURL::libcurl)
        target_compile_definitions(FlexAIDdS PRIVATE FLEXAIDS_HAS_LIBCURL)
    endif()

    if(MSVC)
        target_compile_options(FlexAIDdS PRIVATE /O2 /GL /DNDEBUG /MT)
        target_link_options(FlexAIDdS PRIVATE /LTCG)
        target_compile_definitions(FlexAIDdS PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(FlexAIDdS PRIVATE -O3 -ffast-math -fno-finite-math-only -flto -DNDEBUG
            $<$<AND:$<CXX_COMPILER_ID:AppleClang,Clang>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>>:-mcpu=native>
        )
        if(NOT APPLE)
            target_compile_options(FlexAIDdS PRIVATE -march=native)
            target_link_options(FlexAIDdS PRIVATE -flto -s)
        else()
            target_link_options(FlexAIDdS PRIVATE -flto)
        endif()
    endif()

    # ─── IPO / LTO via CMake property (supplements the -flto flag above) ────
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_ok OUTPUT ipo_err)
    if(ipo_ok)
        # IPO on flexaid_core only — FlexAIDdS links with -flto. Do not enable IPO on
        # the OBJECT library itself: other targets (FlexAID, tests) link those .o files
        # without -flto and clang then reports "file not recognized" on LTO bitcode.
        set_property(TARGET FlexAIDdS PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
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

    # Metal OBJCXX bridges + frameworks come from flexaid_core when
    # FLEXAIDS_USE_METAL=ON (see "Metal sources (attach to flexaid_core once)").
    if(FLEXAIDS_USE_METAL)
        message(STATUS "FlexAIDdS: Metal GPU acceleration via flexaid_core")
    endif()

    message(STATUS "FlexAIDdS fast executable enabled (with --redock RCSB cognate prep)")
endif()

# Attach Metal dispatch runtime to targets that compile hardware_detect /
# UnifiedHardwareDispatch themselves WITHOUT linking flexaid_core.
# Targets that already link flexaid_core already get metal_eval.mm etc.
function(flexaids_attach_metal_dispatch_runtime target_name)
    if(NOT FLEXAIDS_USE_METAL)
        return()
    endif()
    # Avoid duplicate-symbol if the target already pulls flexaid_core objects.
    get_target_property(_flexaids_metal_link_libs ${target_name} LINK_LIBRARIES)
    if(_flexaids_metal_link_libs MATCHES "flexaid_core")
        return()
    endif()
    target_sources(${target_name} PRIVATE
        LIB/metal_eval.mm
        LIB/ShannonThermoStack/ShannonMetalBridge.mm
    )
    target_compile_definitions(${target_name} PRIVATE
        FLEXAIDS_USE_METAL
        FLEXAIDS_HAS_METAL_SHANNON
    )
    if(DEFINED SHANNON_METALLIB)
        target_compile_definitions(${target_name} PRIVATE
            SHANNON_METALLIB_PATH="${SHANNON_METALLIB}"
        )
    endif()
    target_compile_options(${target_name} PRIVATE -fno-objc-arc)
    target_link_libraries(${target_name} PRIVATE
        ${METAL_LIBRARY}
        ${FOUNDATION_LIBRARY}
        ${METALKIT_LIBRARY}
    )
    set_source_files_properties(
        LIB/metal_eval.mm
        LIB/ShannonThermoStack/ShannonMetalBridge.mm
        PROPERTIES LANGUAGE OBJCXX COMPILE_FLAGS "-std=c++26 -fno-objc-arc"
    )
    if(TARGET ShannonEntropyMetal)
        add_dependencies(${target_name} ShannonEntropyMetal)
    endif()
endfunction()

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
    flexaids_attach_metal_dispatch_runtime(benchmark_dispatch)
    message(STATUS "Dispatch benchmark enabled — build target: benchmark_dispatch")
endif()

endif() # FLEXAIDS_BUILD_CORE

# ─── Native CavityDetector CLI ─────────────────────────────────────────────
option(ENABLE_CAVITY_DETECT_CLI "Build native FlexAIDdS cavity detector CLI" ON)
if(ENABLE_CAVITY_DETECT_CLI)
    add_executable(cavity_detect_cli
        LIB/cavity_detect_cli.cpp
    )
    target_link_libraries(cavity_detect_cli PRIVATE flexaid_core)
    target_include_directories(cavity_detect_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    if(MSVC)
        target_compile_options(cavity_detect_cli PRIVATE /O2 /fp:fast)
        target_compile_definitions(cavity_detect_cli PRIVATE
            _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
    else()
        target_compile_options(cavity_detect_cli PRIVATE -O3 -ffast-math -fno-finite-math-only)
    endif()
    flexaids_configure_simd(cavity_detect_cli)
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(cavity_detect_cli PRIVATE OpenMP::OpenMP_CXX)
    endif()
    # Metal CavityDetect bridge + frameworks come from flexaid_core.
    if(FLEXAIDS_USE_METAL AND TARGET CavityDetectMetal)
        add_dependencies(cavity_detect_cli CavityDetectMetal)
    endif()
    message(STATUS "Native CavityDetector CLI enabled -- build target: cavity_detect_cli")
endif()

# ─── Benchmark Dataset Runner ──────────────────────────────────────────────
option(ENABLE_BENCHMARK_DATASETS "Build benchmark dataset runner" ON)
if(ENABLE_BENCHMARK_DATASETS)
    find_package(CURL QUIET)

    set(FLEXAIDDS_POSEBUST_SOURCES
        LIB/PoseBust/Loaders.cpp
        LIB/PoseBust/ChecksChemistry.cpp
        LIB/PoseBust/ChecksGeometry.cpp
        LIB/PoseBust/ChecksProtein.cpp
        LIB/PoseBust/Engine.cpp
        LIB/PoseBust/BustCli.cpp
    )

    add_executable(benchmark_datasets
        LIB/benchmark_datasets.cpp
        LIB/DatasetRunner.cpp
        LIB/FleetRunner.cpp
        LIB/FleetRunner.h
        LIB/DatasetRunnerStats.cpp
        LIB/AsyncPipeline.cpp
        # NativePoseQC diagnostic plus upstream PoseBusters CLI bridge
        ${FLEXAIDDS_POSEBUST_SOURCES}
        LIB/PoseBust/BustCli.h
        LIB/PoseBust/ChecksChemistry.h
        LIB/PoseBust/ChecksGeometry.h
        LIB/PoseBust/ChecksProtein.h
        LIB/PoseBust/Engine.h
        LIB/PoseBust/Loaders.h
        LIB/PoseBust/Types.h
    )
    target_link_libraries(benchmark_datasets PRIVATE flexaid_core)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(benchmark_datasets PRIVATE m)
    endif()
    target_include_directories(benchmark_datasets PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PoseBust
    )
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
    # Metal bridges/frameworks are on flexaid_core when FLEXAIDS_USE_METAL=ON.
    if(FLEXAIDS_USE_METAL)
        add_dependencies(benchmark_datasets
            CavityDetectMetal
            ShannonEntropyMetal
            MetalRMSDMetal
            TurboQuantMetal
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
    # P0 wiring: ensure GPF sources (GrandPartitionFunction, TargetServer, MultiSiteGPF)
    # are ALWAYS compiled for test coverage, regardless of user option.
    # This + default=ON in FlexAIDOptions + conditional in LIB/CMakeLists.txt
    # guarantees the modules are present for ctest -R 'Grand|Target|...'
    set(FLEXAIDS_GRAND_CANONICAL ON CACHE BOOL "Force GPF modules for BUILD_TESTING (P0 gate: always compile Grand/Target/MultiSite)" FORCE)
    message(STATUS "BUILD_TESTING=ON → FLEXAIDS_GRAND_CANONICAL forced ON to guarantee GrandPartitionFunction/TargetServer/MultiSiteGPF compilation and tests")

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

    # flexaids_configure_msvc_test is provided by cmake/FlexAIDOptions.cmake
    # (single definition, available for all test targets below).

    # test_statmech — StatMechEngine unit tests
    flexaids_add_unit_test(test_statmech
        SOURCES
            tests/test_statmech.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME StatMechTests
    )

    # test_thermo_gate — thermodynamic impossibility gate (header-only, no deps)
    flexaids_add_unit_test(test_thermo_gate
        SOURCES tests/test_thermo_gate.cpp
        TEST_NAME ThermoGateTests
    )

    # test_thermo_ledger — Thermo ledger tests
    flexaids_add_unit_test(test_thermo_ledger
        SOURCES
            tests/test_thermo_ledger.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME ThermoLedgerTests
    )

    # test_json_config — JSON parser & config defaults tests (zero dependencies)
    flexaids_add_unit_test(test_json_config
        SOURCES tests/test_json_config.cpp
        TEST_NAME JsonConfigTests
    )

    # test_protocol_config — typed ProtocolConfig env adapter + JSON round-trip
    
    flexaids_add_unit_test(test_fleet_runner
        SOURCES
            tests/test_fleet_runner.cpp
            LIB/FleetRunner.cpp
            LIB/FleetRunner.h
        TEST_NAME FleetRunnerTests
    )

flexaids_add_unit_test(test_protocol_config
        SOURCES
            tests/test_protocol_config.cpp
            LIB/ProtocolConfig.cpp
            LIB/RunReceipt.cpp
        TEST_NAME ProtocolConfigTests
    )

    # test_posebust — NativePoseQC clean-room + BustCli (real Astex 1G9V artifacts)
    flexaids_add_unit_test(test_posebust
        SOURCES
            tests/test_posebust.cpp
            LIB/PoseBust/Loaders.cpp
            LIB/PoseBust/ChecksChemistry.cpp
            LIB/PoseBust/ChecksGeometry.cpp
            LIB/PoseBust/ChecksProtein.cpp
            LIB/PoseBust/Engine.cpp
            LIB/PoseBust/BustCli.cpp
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PoseBust
        TEST_NAME PoseBustTests
    )

    # test_soft_wall — overlap-based soft-core clash potential (header-only)
    flexaids_add_unit_test(test_soft_wall
        SOURCES tests/test_soft_wall.cpp
        TEST_NAME SoftWallTests
    )

    # test_memetic_gate — Wave 3.4 FLEXAIDDS_MEMETIC + WALL_PILOT_PASS (header-only)
    flexaids_add_unit_test(test_memetic_gate
        SOURCES tests/test_memetic_gate.cpp
        TEST_NAME MemeticGateTests
    )

    # test_niche_distance — G4.2 Cartesian vs gene-space niche (header-only)
    flexaids_add_unit_test(test_niche_distance
        SOURCES tests/test_niche_distance.cpp
        TEST_NAME NicheDistanceTests
    )

    # test_new_search_arch — S4 A+B phenotype_unique + basin_reinject (header-only)
    flexaids_add_unit_test(test_new_search_arch
        SOURCES tests/test_new_search_arch.cpp
        TEST_NAME NewSearchArchTests
    )

    # test_metal_coordination — Metal ion coordination potential tests (header-only)
    flexaids_add_unit_test(test_metal_coordination
        SOURCES tests/test_metal_coordination.cpp
        TEST_NAME MetalCoordinationTests
    )

    # test_cf_aggregator — pins the REAL get_cf_evalue() term set (links ic2cf.cpp,
    # not the reduced stub in tests/stubs.cpp)
    flexaids_add_unit_test(test_cf_aggregator
        SOURCES
            tests/test_cf_aggregator.cpp
            tests/cf_aggregator_stubs.cpp
            LIB/ic2cf.cpp
        TEST_NAME CFAggregatorTests
    )

    # test_pb_vdw_parity — pins get_element(type) vs atoms[].element as PoseBusters
    # vdW radius sources (links write_pdb.cpp for the real get_element())
    flexaids_add_unit_test(test_pb_vdw_parity
        SOURCES
            tests/test_pb_vdw_parity.cpp
            LIB/write_pdb.cpp
            LIB/fileio.cpp
        TEST_NAME PbVdwParityTests
    )

    # test_formal_charges — PDB formal charge assignment tests (header-only)
    flexaids_add_unit_test(test_formal_charges
        SOURCES tests/test_formal_charges.cpp
        TEST_NAME FormalChargeTests
    )

    # test_hardware_dispatch — ShannonThermoStack dispatch layer tests
    flexaids_add_unit_test(test_hardware_dispatch
        SOURCES
            tests/test_hardware_dispatch.cpp
            LIB/ShannonThermoStack/ShannonThermoStack.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
            LIB/tENCoM/tencm.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME HardwareDispatchTests
    )

    # test_hardware_detect_dispatch — unified dispatch layer + hardware detection tests
    add_executable(test_hardware_detect_dispatch
        tests/test_hardware_detect_dispatch.cpp
        LIB/hardware_detect.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/statmech.cpp
        LIB/TurboQuant.cpp)
    target_include_directories(test_hardware_detect_dispatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_hardware_detect_dispatch PRIVATE GTest::gtest_main)
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
    flexaids_attach_metal_dispatch_runtime(test_hardware_detect_dispatch)
    if(FLEXAIDS_USE_ROCM)
        target_compile_definitions(test_hardware_detect_dispatch PRIVATE FLEXAIDS_USE_ROCM)
        target_include_directories(test_hardware_detect_dispatch PRIVATE ${hip_INCLUDE_DIRS})
        target_link_libraries(test_hardware_detect_dispatch PRIVATE hip::host)
    endif()
    flexaids_configure_msvc_test(test_hardware_detect_dispatch)
    add_test(NAME HardwareDetectDispatchTests COMMAND test_hardware_detect_dispatch)

    # test_unified_dispatch — runtime HardwareDispatcher tests
    add_executable(test_unified_dispatch
        tests/test_unified_dispatch.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_unified_dispatch PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_unified_dispatch PRIVATE GTest::gtest_main)
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
    flexaids_attach_metal_dispatch_runtime(test_unified_dispatch)
    flexaids_configure_msvc_test(test_unified_dispatch)
    add_test(NAME UnifiedDispatchTests COMMAND test_unified_dispatch)

    # test_tencom_diff — tENCoM differential engine tests
    flexaids_add_unit_test(test_tencom_diff
        SOURCES
            tests/test_tencom_diff.cpp
            LIB/tENCoM/tencom_diff.cpp
            LIB/tENCoM/pdb_calpha.cpp
            LIB/tENCoM/tencm.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
            LIB/encom.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME TENCoMDiffTests
    )

    # test_ptm_attachment — PTMAttachment module tests
    flexaids_add_unit_test(test_ptm_attachment
        SOURCES
            tests/test_ptm_attachment.cpp
            LIB/PTMAttachment/PTMAttachment.cpp
        TEST_NAME PTMAttachmentTests
    )

    # test_gaboom — Genetic Algorithm engine unit tests
    add_executable(test_gaboom
        tests/test_gaboom.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
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
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_gaboom PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
    )
    target_link_libraries(test_gaboom PRIVATE GTest::gtest_main)
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
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
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_shannon_ga PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
    )
    target_link_libraries(test_shannon_ga PRIVATE GTest::gtest_main)
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
    target_link_libraries(test_vcontacts PRIVATE GTest::gtest_main)
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
    target_link_libraries(test_vcontacts_geometry PRIVATE GTest::gtest_main)
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
    target_link_libraries(test_cleft_cavity PRIVATE GTest::gtest_main)
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
    flexaids_add_unit_test(test_encom
        SOURCES
            tests/test_encom.cpp
            LIB/encom.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        TEST_NAME ENCoMTests
    )

    # test_ring_conformer_library — RingConformerLibrary + ring detection tests
    flexaids_add_unit_test(test_ring_conformer_library
        SOURCES
            tests/test_ring_conformer_library.cpp
            LIB/LigandRingFlex/RingConformerLibrary.cpp
        TEST_NAME RingConformerLibraryTests
    )

    # test_sugar_pucker — SugarPucker Cremer-Pople pseudorotation tests
    flexaids_add_unit_test(test_sugar_pucker
        SOURCES
            tests/test_sugar_pucker.cpp
            LIB/LigandRingFlex/SugarPucker.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        TEST_NAME SugarPuckerTests
    )

    # test_chiral_center — ChiralCenterGene R/S stereocenter tests
    flexaids_add_unit_test(test_chiral_center
        SOURCES
            tests/test_chiral_center.cpp
            LIB/ChiralCenter/ChiralCenterGene.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        TEST_NAME ChiralCenterTests
    )

    # test_cavity_detect — CavityDetect SURFNET cavity detection tests
    flexaids_add_unit_test(test_cavity_detect
        SOURCES
            tests/test_cavity_detect.cpp
            LIB/CavityDetect/CavityDetect.cpp
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME CavityDetectTests
    )

    # test_read_spheres — read_spheres HETATM sphere PDB ingest
    flexaids_add_unit_test(test_read_spheres
        SOURCES
            tests/test_read_spheres.cpp
            LIB/read_spheres.cpp
            LIB/fileio.cpp
        TEST_NAME read_spheres
    )

    # test_ensemble_pipeline — 4-layer reproducibility helpers (header-only)
    flexaids_add_unit_test(test_ensemble_pipeline
        SOURCES
            tests/test_ensemble_pipeline.cpp
        TEST_NAME EnsemblePipelineTests
    )

    # test_mif_grid — MIFGrid molecular interaction field tests
    flexaids_add_unit_test(test_mif_grid
        SOURCES tests/test_mif_grid.cpp
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME MIFGridTests
    )

    # test_reflig_spectrophore — RefLigSeed + Spectrophore tests
    flexaids_add_unit_test(test_reflig_spectrophore
        SOURCES tests/test_reflig_spectrophore.cpp
        LINK_OPENMP
        TEST_NAME RefLigSpectrophoreTests
    )

    # test_binding_residues — Key binding-site residue identification tests
    flexaids_add_unit_test(test_binding_residues
        SOURCES
            tests/test_binding_residues.cpp
            LIB/set_intprob.cpp
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME BindingResiduesTests
    )

    # test_binding_mode_statmech — BindingMode ↔ StatMechEngine integration tests
    add_executable(test_binding_mode_statmech
        tests/test_binding_mode_statmech.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_binding_mode_statmech PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_statmech PRIVATE GTest::gtest_main)
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

    # test_classic_entropy_ranking — classic FlexAID soft-β / ACF rank-0 (rollback: force_cf_rank_emission)
    add_executable(test_classic_entropy_ranking
        tests/test_classic_entropy_ranking.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/ProtocolConfig.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
    target_include_directories(test_classic_entropy_ranking PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_classic_entropy_ranking PRIVATE GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_classic_entropy_ranking PRIVATE /O2)
    else()
        target_compile_options(test_classic_entropy_ranking PRIVATE -O2)
    endif()
    if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(test_classic_entropy_ranking PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_link_libraries(test_classic_entropy_ranking PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_classic_entropy_ranking PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_classic_entropy_ranking)
    add_test(NAME ClassicEntropyRankingTests COMMAND test_classic_entropy_ranking)

    # test_ccbm — Conformer-Coupled Binding Modes (CCBM) unit tests
    add_executable(test_ccbm
        tests/test_ccbm.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_ccbm PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ccbm PRIVATE GTest::gtest_main)
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
        LIB/ProtocolConfig.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
    target_link_libraries(test_binding_mode_advanced PRIVATE GTest::gtest_main)
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_ga_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ga_core PRIVATE GTest::gtest_main)
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_ga_validation PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ga_validation PRIVATE GTest::gtest_main)
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_ga_operators PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_ga_operators PRIVATE GTest::gtest_main Eigen3::Eigen)
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_binding_mode_comparison PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_comparison PRIVATE GTest::gtest_main Eigen3::Eigen)
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
        target_link_libraries(test_soft_contact_matrix PRIVATE GTest::gtest_main)
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

    # test_compare_energy_matrices — compare 256×256 vs 40×40 energy matrices
    if(FLEXAIDS_USE_256_MATRIX)
        add_executable(test_compare_energy_matrices
            tests/test_compare_energy_matrices.cpp
        )
        target_include_directories(test_compare_energy_matrices PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB
            ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        )
        target_link_libraries(test_compare_energy_matrices PRIVATE GTest::gtest_main)
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
    flexaids_add_unit_test(test_ion_handling
        SOURCES
            tests/test_ion_handling.cpp
            LIB/assign_radii.cpp
            LIB/assign_types.cpp
            LIB/fileio.cpp
            LIB/tENCoM/tencm.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        TEST_NAME IonHandlingTests
    )

    # test_entropy_ga — SMFREE fitness model, Boltzmann blending, entropy fields
    add_executable(test_entropy_ga
        tests/test_entropy_ga.cpp
        tests/stubs.cpp
        LIB/gaboom.cpp
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_entropy_ga PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_entropy_ga PRIVATE GTest::gtest_main)
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
    flexaids_add_unit_test(test_cube_grid
        SOURCES
            tests/test_cube_grid.cpp
            LIB/maps.cpp
        TEST_NAME CubeGridTests
    )

    # test_parallel_dock — GridDecomposer, SharedPosePool, StatMechEngine merge tests
    flexaids_add_unit_test(test_parallel_dock
        SOURCES
            tests/test_parallel_dock.cpp
            LIB/GridDecomposer.cpp
            LIB/SharedPosePool.cpp
            LIB/statmech.cpp
            LIB/TurboQuant.cpp
            LIB/UnifiedHardwareDispatch.cpp
            LIB/hardware_detect.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        CONFIGURE_SIMD
        LINK_OPENMP
        TEST_NAME ParallelDockTests
    )

    # test_knowledge_pool — TargetKnowledgeBase + SharedPosePool + MPIMergeable
    flexaids_add_unit_test(test_knowledge_pool
        SOURCES
            tests/test_knowledge_pool.cpp
            LIB/TargetKnowledgeBase.cpp
            LIB/SharedPosePool.cpp
        COMPILE_OPTIONS -std=c++26
        TEST_NAME KnowledgePoolTests
    )

    # ── Mol2Reader / SdfReader tests ──────────────────────────────────────
    flexaids_add_unit_test(test_mol2_sdf_reader
        SOURCES
            tests/test_mol2_sdf_reader.cpp
            tests/stubs.cpp
            LIB/Mol2Reader.cpp
            LIB/SdfReader.cpp
            LIB/buildcc.cpp
            LIB/buildic.cpp
            LIB/buildlist.cpp
            LIB/read_pdb.cpp
            LIB/read_coor.cpp
            LIB/read_conect.cpp
            LIB/assign_radii_types.cpp
            LIB/geometry.cpp
            LIB/LigandRingFlex/SugarPucker.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN FLEXAIDS_READER_REAL_GEOMETRY
        COMPILE_OPTIONS -std=c++26
        CONFIGURE_SIMD
        GTEST_MAIN_ONLY
        TEST_NAME Mol2SdfReaderTests
    )

    # ── top.cpp helpers (seam A: logic lifted out of the main() TU) ───────
    flexaids_add_unit_test(test_top_helpers
        SOURCES
            tests/test_top_helpers.cpp
            LIB/top_helpers.cpp
        COMPILE_OPTIONS -std=c++26
        GTEST_MAIN_ONLY
        TEST_NAME TopHelpersTests
    )

    # ── read_lig latm inclusive range (pilot8 missing last HETTYP atom) ───
    # Uses real fileio.cpp (OpenFile_B / Terminate throws) — no stubs OpenFile.
    flexaids_add_unit_test(test_read_lig_latm
        SOURCES
            tests/test_read_lig_latm.cpp
            LIB/read_lig.cpp
            LIB/fileio.cpp
            LIB/read_conect.cpp
            LIB/buildlist.cpp
            LIB/assign_shift.cpp
            LIB/bondedlist.cpp
            LIB/update_bonded.cpp
            LIB/shortest_path.cpp
            LIB/assign_shortflex.cpp
            LIB/free_resid.cpp
            LIB/buildcc.cpp
            LIB/geometry.cpp
        COMPILE_OPTIONS -std=c++26
        GTEST_MAIN_ONLY
        TEST_NAME ReadLigLatmTests
    )

    # test_turboquant — TurboQuant quantisation library unit tests
    flexaids_add_unit_test(test_turboquant
        SOURCES
            tests/test_turboquant.cpp
            LIB/TurboQuant.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        TEST_NAME TurboQuantTests
    )

    # test_dataset_runner — benchmark dataset runner tests
    add_executable(test_dataset_runner
        tests/test_dataset_runner.cpp
        LIB/DatasetRunner.cpp
        LIB/DatasetRunnerStats.cpp
        LIB/DatasetRunnerProvenance.cpp
        LIB/AsyncPipeline.cpp
        ${FLEXAIDDS_POSEBUST_SOURCES}
    )
    target_include_directories(test_dataset_runner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_dataset_runner PRIVATE flexaid_core GTest::gtest_main)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(test_dataset_runner PRIVATE m)
    endif()
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

    # test_cofactor_blacklist — ligand-extraction cofactor blacklist (1HNN SAH/SKF)
    add_executable(test_cofactor_blacklist
        tests/test_cofactor_blacklist.cpp
        LIB/DatasetRunner.cpp
        LIB/DatasetRunnerStats.cpp
        LIB/AsyncPipeline.cpp
        ${FLEXAIDDS_POSEBUST_SOURCES}
    )
    target_include_directories(test_cofactor_blacklist PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_cofactor_blacklist PRIVATE flexaid_core GTest::gtest_main)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(test_cofactor_blacklist PRIVATE m)
    endif()
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
        LIB/ProtocolConfig.cpp
    )
    target_include_directories(test_binding_mode_vibrational PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
    )
    target_link_libraries(test_binding_mode_vibrational PRIVATE GTest::gtest_main)
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
    flexaids_add_unit_test(test_fast_optics
        SOURCES
            tests/test_fast_optics.cpp
            LIB/fast_optics.cpp
        TEST_NAME FastOPTICSTests
    )

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
    target_link_libraries(test_tencom_entropy_diff PRIVATE GTest::gtest_main)
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
    flexaids_add_unit_test(test_target_validation
        SOURCES
            tests/test_target_validation.cpp
            LIB/TargetValidation.cpp
        TEST_NAME TargetValidationTests
    )

    # ── GrandPartitionFunction tests ──────────────────────────────────────
    add_executable(test_grand_partition
        tests/test_grand_partition.cpp
        LIB/GrandPartitionFunction.cpp
        LIB/MultiSiteGPF.cpp
        LIB/statmech.cpp
        LIB/UnifiedHardwareDispatch.cpp
        LIB/hardware_detect.cpp
    )
    target_include_directories(test_grand_partition PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_grand_partition PRIVATE GTest::gtest_main)
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
    target_link_libraries(test_target_server PRIVATE GTest::gtest_main)
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
    target_link_libraries(test_multi_site_gpf PRIVATE GTest::gtest_main)
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
    flexaids_add_unit_test(test_gist_hbond
        SOURCES
            tests/test_gist_hbond.cpp
            LIB/GISTEvaluator.cpp
            LIB/HBondEvaluator.cpp
        TEST_NAME GISTHBondTests
    )

    # test_buffer_safety — bounded legacy config parsing regressions
    flexaids_add_unit_test(test_buffer_safety
        SOURCES tests/test_buffer_safety.cpp
        TEST_NAME BufferSafetyTests
    )

    # test_hungarian_rmsd_bounds — post-cluster RMSDST Hungarian SEGV regression
    # (cad_only_no_pose_pdb: atoms[k+l] walked past ligand into NULL coor_ref)
    flexaids_add_unit_test(test_hungarian_rmsd_bounds
        SOURCES
            tests/test_hungarian_rmsd_bounds.cpp
            tests/hungarian_rmsd_stubs.cpp
            LIB/calc_rmsd.cpp
            LIB/geometry.cpp
        TEST_NAME HungarianRmsdBoundsTests
    )

    # test_shell_exec — path-safe shell quoting / argv exec helpers
    flexaids_add_unit_test(test_shell_exec
        SOURCES tests/test_shell_exec.cpp
        TEST_NAME ShellExecTests
    )

    message(STATUS "Unit tests enabled — build targets: test_statmech, test_hardware_dispatch, test_tencom_diff, test_target_validation, test_grand_partition, test_target_server, test_multi_site_gpf")
    # test_hbond_potential — Angular-dependent H-bond potential unit tests
    flexaids_add_unit_test(test_hbond_potential
        SOURCES tests/test_hbond_potential.cpp
        TEST_NAME HBondPotentialTests
    )

    # test_gist_grid — GIST desolvation grid unit tests
    flexaids_add_unit_test(test_gist_grid
        SOURCES tests/test_gist_grid.cpp
        TEST_NAME GISTGridTests
    )

    # test_ga_diversity — GA population diversity monitoring unit tests
    flexaids_add_unit_test(test_ga_diversity
        SOURCES tests/test_ga_diversity.cpp
        TEST_NAME GADiversityTests
    )

    # test_soa_conversion — AoS↔SoA atom data conversion unit tests
    flexaids_add_unit_test(test_soa_conversion
        SOURCES tests/test_soa_conversion.cpp
        CONFIGURE_SIMD
        TEST_NAME SoAConversionTests
    )

    # test_dift — Discrete Fourier Transform torsional parametrization
    add_executable(test_dift
        tests/test_dift.cpp
        LIB/DiFT/DiFT.cpp
    )
    target_include_directories(test_dift PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/DiFT
    )
    target_link_libraries(test_dift PRIVATE GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_dift PRIVATE /Od /Zi)
    else()
        target_compile_options(test_dift PRIVATE -O0 -g3)
    endif()
    flexaids_configure_msvc_test(test_dift)
    add_test(NAME DiFTTests COMMAND test_dift)

    # test_distributed_backend — Distributed backend abstraction unit tests
    flexaids_add_unit_test(test_distributed_backend
        SOURCES tests/test_distributed_backend.cpp
        TEST_NAME DistributedBackendTests
    )

    # test_ga_context — GAContext re-entrant state + GPUContextPool singleton
    flexaids_add_unit_test(test_ga_context
        SOURCES
            tests/test_ga_context.cpp
            LIB/GAContext.cpp
        LINK_LIBRARIES Eigen3::Eigen
        DEFINES FLEXAIDS_HAS_EIGEN
        TEST_NAME GAContextTests
    )

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
        target_link_libraries(test_natural PRIVATE GTest::gtest_main)
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
    target_link_libraries(test_nascent_chain_scheduler PRIVATE GTest::gtest_main)
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

    # test_process_ligand — ProcessLigand pipeline (must live under BUILD_TESTING,
    # not BUILD_SWIFT_BRIDGE, or the target is never generated when Swift is OFF).
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
    target_link_libraries(test_process_ligand PRIVATE GTest::gtest_main)
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

    # ── CMA-ES mock unit tests (opt-in search backend; no flexaid_core) ───────
    # Links only LIB/cmaes_search.cpp + tests/cmaes_mock_seams_stub.cpp + GTest.
    # See also tests/CMakeLists_cmaes.fragment for the standalone registration recipe.
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/LIB/cmaes_search.cpp"
       AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cmaes_search.cpp")
        include(GoogleTest)
        set(_cmaes_test_sources
            tests/test_cmaes_search.cpp
            LIB/cmaes_search.cpp
        )
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/cmaes_mock_seams_stub.cpp")
            list(APPEND _cmaes_test_sources tests/cmaes_mock_seams_stub.cpp)
        endif()
        add_executable(test_cmaes_search ${_cmaes_test_sources})
        target_include_directories(test_cmaes_search PRIVATE ${CMAKE_SOURCE_DIR}/LIB)
        target_link_libraries(test_cmaes_search PRIVATE GTest::gtest_main)
        set_target_properties(test_cmaes_search PROPERTIES
            CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
        if(NOT MSVC)
            target_compile_options(test_cmaes_search PRIVATE -O2)
        endif()
        gtest_discover_tests(test_cmaes_search
            TEST_PREFIX "CmaesSearch."
            DISCOVERY_MODE PRE_TEST
            PROPERTIES LABELS "cmaes;unit;mock" TIMEOUT 120)
        add_test(NAME CmaesSearchMock COMMAND test_cmaes_search)
        set_tests_properties(CmaesSearchMock PROPERTIES LABELS "cmaes;unit;mock" TIMEOUT 120)
        message(STATUS "Registered test_cmaes_search (mock CMA-ES unit tests)")
        unset(_cmaes_test_sources)
    endif()

    # test_coarse_screen — NRGRank CoarseScreen + TwoStageScreen unit tests
    add_executable(test_coarse_screen
        tests/test_coarse_screen.cpp
        LIB/CoarseScreen.cpp
        LIB/TwoStageScreen.cpp
    )
    target_include_directories(test_coarse_screen PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB)
    target_link_libraries(test_coarse_screen PRIVATE GTest::gtest_main)
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

    # ─── C1: FlexAID --help error-path leak probe ─────────────────────────
    # Shipping binary only. --help takes Terminate(0) and skips the entire
    # teardown (throw from fileio.cpp, catch after free block in top.cpp).
    # Process path is always checked. On Linux *Address*Sanitizer builds the
    # criterion is inverted: an LSan report is PASS (the finding). Structural
    # prediction: 9 objects = 3 DIRECT (FA/GB/VC) + 6 INDIRECT (contacts + 5
    # VC arrays). Bytes are NOT asserted (macOS size classes ≠ glibc).
    # Does NOT cover teardown — only C2 can.
    # Gate MUST match the sanitizer, not the flag family: -fsanitize=thread
    # contains "fsanitize" but ships no LSan; a bare "fsanitize" match arms
    # EXPECT_LEAK under TSan and the inverted criterion fails a clean TSan
    # run (Honey, 2026-07-31). Match address specifically so ASan+UBSan still
    # arms and TSan demotes to process-path-only.
    # See tests/cmake/c1_help_error_path_probe.cmake.
    if(TARGET FlexAID)
        string(CONCAT _c1_san_flags
            "${CMAKE_CXX_FLAGS}" "${CMAKE_C_FLAGS}" "${CMAKE_EXE_LINKER_FLAGS}")
        set(_c1_expect_leak OFF)
        # Match -fsanitize=address and -fsanitize=undefined,address etc.
        # Do NOT match -fsanitize=thread (TSan has no LSan).
        if(UNIX AND NOT APPLE AND _c1_san_flags MATCHES "fsanitize=[a-z,]*address")
            set(_c1_expect_leak ON)
        endif()
        add_test(
            NAME FlexAID_C1_HelpErrorPathProbe
            COMMAND ${CMAKE_COMMAND}
                -DFlexAID_BINARY=$<TARGET_FILE:FlexAID>
                -DEXPECT_LEAK=${_c1_expect_leak}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/c1_help_error_path_probe.cmake
        )
        set_tests_properties(FlexAID_C1_HelpErrorPathProbe PROPERTIES
            TIMEOUT 60
            LABELS "c1;lsan_probe;top_cpp"
        )
        if(_c1_expect_leak)
            set_property(TEST FlexAID_C1_HelpErrorPathProbe APPEND PROPERTY
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=1:halt_on_error=0:strict_string_checks=1")
            # PASS_REGULAR_EXPRESSION: green is not enough — output must prove
            # which arm ran. If the gate silently demotes, this string is
            # absent and ctest fails (Honey: check that disarms itself).
            # Note: when PASS_REGULAR_EXPRESSION is set, CMake ignores the
            # process exit code; FAIL lines must not match the PASS pattern
            # (they use C1_RESULT=FAIL).
            set_property(TEST FlexAID_C1_HelpErrorPathProbe APPEND PROPERTY
                PASS_REGULAR_EXPRESSION "C1_RESULT=PASS ARM=lsan_inverted")
            message(STATUS
                "C1 --help probe ARM=lsan_inverted "
                "(Linux + -fsanitize=*address* — an LSan report is required; "
                "the direct/indirect split is asserted only when "
                "C1_EXPECT_DIRECT/C1_EXPECT_INDIRECT are set, otherwise "
                "measured and reported)")
        else()
            set_property(TEST FlexAID_C1_HelpErrorPathProbe APPEND PROPERTY
                PASS_REGULAR_EXPRESSION "C1_RESULT=PASS ARM=process_path_only")
            message(STATUS
                "C1 --help probe ARM=process_path_only "
                "(not Linux+address; TSan/macOS/Release demote here — "
                "ctest green on this arm does NOT mean LSan ran)")
        endif()
    endif()


    message(STATUS "Unit tests enabled — build targets: test_statmech, test_hardware_dispatch, test_tencom_diff, test_target_validation, test_grand_partition, test_target_server, test_ga_context, test_nascent_chain_scheduler, test_process_ligand, test_cmaes_search, test_coarse_screen, FlexAID_C1_HelpErrorPathProbe")
endif()

# ─── Swift Bridge (macOS only) ────────────────────────────────────────────
# EXPERIMENTAL / STUB: Swift bindings and prototypes are NOT core components.
# Default OFF. See .gitmodules for recommended submodule management.
# Swift/TS code evolves separately; do not depend on in core C++ builds.
option(BUILD_SWIFT_BRIDGE "Build Swift bridge static library (macOS only) [experimental stub, not core]" OFF)
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

    if(BUILD_TESTING)
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
    target_link_libraries(test_ligand_ring_flex PRIVATE GTest::gtest_main)
    if(MSVC)
        target_compile_options(test_ligand_ring_flex PRIVATE /Od /Zi)
    else()
        target_compile_options(test_ligand_ring_flex PRIVATE -O0 -g3)
    endif()
    target_link_libraries(test_ligand_ring_flex PRIVATE Eigen3::Eigen)
    target_compile_definitions(test_ligand_ring_flex PRIVATE FLEXAIDS_HAS_EIGEN)
    flexaids_configure_msvc_test(test_ligand_ring_flex)
    add_test(NAME LigandRingFlexTests COMMAND test_ligand_ring_flex)
    endif() # BUILD_TESTING inside SWIFT bridge block

    # test_natural — RibosomeElongation, TransloconInsertion, NucleationSeedDetector
    # Guarded: an earlier block already defines this target with the same
    # NATURaLTests name; without the guard a fresh configure fails with a
    # duplicate-target error.
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
    target_link_libraries(test_natural PRIVATE GTest::gtest_main)
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
        LIB/coarse_init.cpp
        LIB/ThermodynamicEngine.cpp
        LIB/LigandRingFlex/RingConformerLibrary.cpp
        LIB/LigandRingFlex/SugarPucker.cpp
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
    target_link_libraries(test_docking_pipeline PRIVATE GTest::gtest_main)
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

# ─── probe_cf: frozen-pose CF scoring diagnostic ─────────────────────────────
# Standalone driver that scores a receptor + frozen ligand pose through the
# engine's own score_native_pose()/vcfunction() path (no GA) and emits a
# decomposed CF breakdown as JSON. It drives the FlexAIDdS binary in
# FLEXAIDDS_SCORE_NATIVE + FLEXAIDDS_NATIVE_ONLY mode (see tools/probe_cf.cpp),
# so it needs no link against flexaid_core and cannot diverge from the engine.
option(ENABLE_PROBE_CF "Build the probe_cf frozen-pose CF diagnostic" ON)
# probe_cf uses POSIX-only <unistd.h>/getpid/popen — gate off on MSVC.
if(ENABLE_PROBE_CF AND NOT MSVC)
    add_executable(probe_cf tools/probe_cf.cpp)
    set_target_properties(probe_cf PROPERTIES
        CXX_STANDARD ${CMAKE_CXX_STANDARD} CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_definitions(probe_cf PRIVATE _CRT_SECURE_NO_WARNINGS)
    endif()
    message(STATUS "probe_cf frozen-pose CF diagnostic enabled — build target: probe_cf")
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
