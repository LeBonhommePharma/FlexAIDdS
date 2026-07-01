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
#   - Standard include paths
#   - GTest linking
#   - Common compile flags (MSVC vs others)
#   - SIMD / MSVC helper configuration
#   - Registers the test
#
# Usage:
#   flexaids_add_unit_test(my_test
#       SOURCES
#           tests/test_my.cpp
#           LIB/some.cpp
#       LINK_LIBRARIES
#           Eigen3::Eigen
#           OpenMP::OpenMP_CXX
#       DEFINES
#           FOO=1
#       INCLUDES
#           ${CMAKE_CURRENT_SOURCE_DIR}/extra
#   )
#
# Minimal:
#   flexaids_add_unit_test(test_foo tests/test_foo.cpp)
# ------------------------------------------------------------------------------
function(flexaids_add_unit_test name)
    set(options "")
    set(oneValueArgs TEST_NAME)
    set(multiValueArgs SOURCES LINK_LIBRARIES DEFINES INCLUDES)
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
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/NATURaL
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/CavityDetect
        ${CMAKE_CURRENT_SOURCE_DIR}/LIB/PTMAttachment
        ${ARG_INCLUDES}
    )

    target_link_libraries(${name} PRIVATE
        GTest::gtest
        GTest::gtest_main
        ${ARG_LINK_LIBRARIES}
    )

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

    flexaids_configure_msvc_test(${name})

    foreach(def ${ARG_DEFINES})
        target_compile_definitions(${name} PRIVATE ${def})
    endforeach()

    if(ARG_TEST_NAME)
        add_test(NAME ${ARG_TEST_NAME} COMMAND ${name})
    else()
        add_test(NAME ${name} COMMAND ${name})
    endif()
endfunction()

message(STATUS "FlexAIDHelpers loaded — test and build helpers available")
