# cmake/FlexAIDHelpers.cmake
#
# FlexAIDdS — Professional CMake helper functions.
#
# This centralizes repeated boilerplate for a clean, maintainable build system.
# Reduces duplication in root CMakeLists.txt (monolith reduction + professionalism).

cmake_minimum_required(VERSION 3.28)

# ------------------------------------------------------------------------------
# flexaids_add_unit_test
#
# Adds a GoogleTest unit test with standard FlexAIDdS conventions:
#   - Standard include paths (LIB + common submodules)
#   - GTest linking (gtest_main by default; the imported target carries gtest)
#   - Common compile flags (MSVC vs others; -O2 default)
#   - Optional: SIMD configuration, OpenMP link, custom compile opts
#   - MSVC test helper + explicit DEFINES
#   - Registers the test (custom TEST_NAME optional)
#
# Usage:
#   flexaids_add_unit_test(my_test
#       SOURCES
#           tests/test_my.cpp
#           LIB/some.cpp
#       LINK_LIBRARIES
#           Eigen3::Eigen
#           flexaid_core
#       DEFINES
#           FOO=1
#       INCLUDES
#           ${CMAKE_CURRENT_SOURCE_DIR}/extra
#       COMPILE_OPTIONS
#           -std=c++26
#       CONFIGURE_SIMD
#       LINK_OPENMP
#   )
#
# Flags (options):
#   CONFIGURE_SIMD     - call flexaids_configure_simd
#   LINK_OPENMP        - conditionally link OpenMP::OpenMP_CXX if enabled
#   GTEST_MAIN_ONLY    - link only GTest::gtest_main (omit gtest)
#   NO_DEFAULT_COMPILE_OPTS - skip the built-in /O2 or -O2
#
# Minimal:
#   flexaids_add_unit_test(test_foo tests/test_foo.cpp)
# ------------------------------------------------------------------------------
function(flexaids_add_unit_test name)
    set(options CONFIGURE_SIMD LINK_OPENMP GTEST_MAIN_ONLY NO_DEFAULT_COMPILE_OPTS)
    set(oneValueArgs TEST_NAME)
    set(multiValueArgs SOURCES LINK_LIBRARIES DEFINES INCLUDES COMPILE_OPTIONS)
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT ARG_SOURCES)
        # Allow simple call: flexaids_add_unit_test(name source1 source2 ...)
        set(ARG_SOURCES ${ARGN})
    endif()

    add_executable(${name} ${ARG_SOURCES})

    target_include_directories(${name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/tENCoM
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ShannonThermoStack
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/LigandRingFlex
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/ChiralCenter
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PTMAttachment
        ${ARG_INCLUDES}
    )

    if(ARG_GTEST_MAIN_ONLY)
        target_link_libraries(${name} PRIVATE
            GTest::gtest_main
            ${ARG_LINK_LIBRARIES}
        )
    else()
        target_link_libraries(${name} PRIVATE
            GTest::gtest_main
            ${ARG_LINK_LIBRARIES}
        )
    endif()

    if(NOT ARG_NO_DEFAULT_COMPILE_OPTS)
        if(MSVC)
            target_compile_options(${name} PRIVATE /O2)
            target_compile_definitions(${name} PRIVATE
                _CRT_SECURE_NO_WARNINGS
                _USE_MATH_DEFINES
                NOMINMAX
            )
        else()
            target_compile_options(${name} PRIVATE -O2)
        endif()
    endif()

    if(ARG_COMPILE_OPTIONS)
        target_compile_options(${name} PRIVATE ${ARG_COMPILE_OPTIONS})
    endif()

    flexaids_configure_msvc_test(${name})

    if(ARG_CONFIGURE_SIMD)
        flexaids_configure_simd(${name})
    endif()

    if(ARG_LINK_OPENMP)
        if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
            target_link_libraries(${name} PRIVATE OpenMP::OpenMP_CXX)
        endif()
    endif()

    foreach(def ${ARG_DEFINES})
        target_compile_definitions(${name} PRIVATE ${def})
    endforeach()

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(${name} PRIVATE m)
    endif()

    if(ARG_TEST_NAME)
        add_test(NAME ${ARG_TEST_NAME} COMMAND ${name})
    else()
        add_test(NAME ${name} COMMAND ${name})
    endif()
endfunction()

message(STATUS "FlexAIDHelpers loaded — test and build helpers available")
