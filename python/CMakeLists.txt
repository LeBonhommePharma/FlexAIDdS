# python/CMakeLists.txt — pybind11 Python bindings for FlexAIDdS
#
# Builds the _core Python extension module exposing StatMechEngine,
# ENCoMEngine, TurboQuant, and related types to Python.

if(NOT BUILD_PYTHON_BINDINGS OR NOT pybind11_FOUND)
    return()
endif()

pybind11_add_module(_core
    bindings/core_bindings.cpp
    ${CMAKE_SOURCE_DIR}/LIB/statmech.cpp
    ${CMAKE_SOURCE_DIR}/LIB/TurboQuant.cpp
    ${CMAKE_SOURCE_DIR}/LIB/UnifiedHardwareDispatch.cpp
    ${CMAKE_SOURCE_DIR}/LIB/hardware_detect.cpp
    ${CMAKE_SOURCE_DIR}/LIB/encom.cpp
    ${CMAKE_SOURCE_DIR}/LIB/ShannonThermoStack/ShannonThermoStack.cpp
    ${CMAKE_SOURCE_DIR}/LIB/tENCoM/tencm.cpp
    ${CMAKE_SOURCE_DIR}/LIB/DiFT/DiFT.cpp
    ${CMAKE_SOURCE_DIR}/LIB/fast_optics.cpp
    ${CMAKE_SOURCE_DIR}/LIB/GrandPartitionFunction.cpp
    ${CMAKE_SOURCE_DIR}/LIB/GrandCanonicalEngine.cpp
)
target_compile_definitions(_core PRIVATE FLEXAIDS_ENABLE_MUVT=1)

# Conditionally add 256x256 matrix bindings
if(FLEXAIDS_USE_256_MATRIX)
    target_sources(_core PRIVATE bindings/bindings_matrix.cpp)
    target_compile_definitions(_core PRIVATE FLEXAIDS_USE_256_MATRIX)
endif()

target_include_directories(_core PRIVATE
    ${CMAKE_SOURCE_DIR}/LIB
    ${CMAKE_SOURCE_DIR}/LIB/ShannonThermoStack
    ${CMAKE_SOURCE_DIR}/LIB/tENCoM
)

target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})

if(MSVC)
    target_compile_definitions(_core PRIVATE
        _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES NOMINMAX)
endif()

# Link Eigen3 (required)
target_link_libraries(_core PRIVATE Eigen3::Eigen)
target_compile_definitions(_core PRIVATE FLEXAIDS_HAS_EIGEN)

# OpenMP
if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(_core PRIVATE OpenMP::OpenMP_CXX)
endif()

# SIMD
flexaids_configure_simd(_core)

message(STATUS "Python bindings enabled — build target: _core.so")
message(STATUS "Install with: pip install -e python/")
