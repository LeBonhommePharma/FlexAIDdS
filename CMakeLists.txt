# Updated CMakeLists.txt with Bonhomme Fleet support

# Add FleetRunner
add_library(FleetRunner LIB/FleetRunner.cpp)
target_link_libraries(FleetRunner PUBLIC CoreLib Metal::Metal)

# New benchmark target
target_sources(flexaidds PRIVATE tools/metal_microbench.cpp)

# Tests
if(BUILD_TESTING)
    add_executable(test_fleet_runner tests/cpp/test_fleet_runner.cpp)
    target_link_libraries(test_fleet_runner PRIVATE FleetRunner gtest)
endif()

# Ensure all new files compiled with C++26, Metal flags, etc.
# Full integration for DatasetRunner and fleet mode
