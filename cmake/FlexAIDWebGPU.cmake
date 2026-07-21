# FlexAIDWebGPU.cmake — optional WebGPU (Dawn or wgpu-native) backend detection.
#
# Additive and best-effort: if neither Dawn nor wgpu-native is found on the
# system, FLEXAIDS_USE_WEBGPU is forced OFF and the target is skipped
# gracefully (no FATAL_ERROR), matching the CPU-default/GPU-additive
# contract for src/backends/.
#
# Detection order:
#   1. User-supplied FLEXAIDS_WEBGPU_ROOT (CMAKE_PREFIX_PATH hint)
#   2. Dawn config package (find_package(Dawn CONFIG))
#   3. wgpu-native prebuilt (libwgpu_native + webgpu.h via pkg-config or WGPU_NATIVE_DIR)

option(FLEXAIDS_USE_WEBGPU "Enable experimental WebGPU compute backend (Dawn or wgpu-native)" OFF)

if(FLEXAIDS_USE_WEBGPU)
    set(_flexaids_webgpu_found FALSE)
    set(FLEXAIDS_WEBGPU_LIBRARIES "")
    set(FLEXAIDS_WEBGPU_INCLUDE_DIRS "")

    if(FLEXAIDS_WEBGPU_ROOT)
        list(APPEND CMAKE_PREFIX_PATH "${FLEXAIDS_WEBGPU_ROOT}")
    endif()

    find_package(Dawn CONFIG QUIET)
    if(Dawn_FOUND OR TARGET dawn::webgpu_dawn)
        set(_flexaids_webgpu_found TRUE)
        list(APPEND FLEXAIDS_WEBGPU_LIBRARIES dawn::webgpu_dawn)
        message(STATUS "WebGPU backend: found Dawn (dawn::webgpu_dawn)")
    endif()

    if(NOT _flexaids_webgpu_found)
        find_library(WGPU_NATIVE_LIBRARY
            NAMES wgpu_native wgpu
            HINTS "${FLEXAIDS_WEBGPU_ROOT}" ENV WGPU_NATIVE_DIR
            PATH_SUFFIXES lib release debug
        )
        find_path(WGPU_NATIVE_INCLUDE_DIR
            NAMES webgpu/webgpu.h webgpu.h
            HINTS "${FLEXAIDS_WEBGPU_ROOT}" ENV WGPU_NATIVE_DIR
            PATH_SUFFIXES include
        )
        if(WGPU_NATIVE_LIBRARY AND WGPU_NATIVE_INCLUDE_DIR)
            set(_flexaids_webgpu_found TRUE)
            list(APPEND FLEXAIDS_WEBGPU_LIBRARIES "${WGPU_NATIVE_LIBRARY}")
            list(APPEND FLEXAIDS_WEBGPU_INCLUDE_DIRS "${WGPU_NATIVE_INCLUDE_DIR}")
            message(STATUS "WebGPU backend: found wgpu-native (${WGPU_NATIVE_LIBRARY})")
        endif()
    endif()

    if(NOT _flexaids_webgpu_found)
        message(STATUS
            "FLEXAIDS_USE_WEBGPU=ON but neither Dawn nor wgpu-native was found — "
            "skipping WebGPU backend (set FLEXAIDS_WEBGPU_ROOT or WGPU_NATIVE_DIR "
            "to point at an installed copy). Falling back to FLEXAIDS_USE_WEBGPU=OFF.")
        set(FLEXAIDS_USE_WEBGPU OFF CACHE BOOL "Enable experimental WebGPU compute backend" FORCE)
    endif()
endif()
