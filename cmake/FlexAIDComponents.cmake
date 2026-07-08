# cmake/FlexAIDComponents.cmake
#
# FlexAIDdS — Logical component targets (modern CMake best practice).
#
# This provides stable, namespaced targets that:
#   - Hide raw names like "flexaid_core"
#   - Allow granular linking (e.g. a tool only needs VoronoiCF + Eigen)
#   - Prepare for future extraction of smaller static/object libs
#   - Make the audit-recommended "split" visible in CMake code
#
# Include *after* add_subdirectory(LIB) because we create ALIASes and
# INTERFACE targets that may depend on flexaid_core or its deps.
#
# Example usage in new or refactored code:
#   target_link_libraries(my_tool PRIVATE FlexAID::Core)
#   target_link_libraries(vcf_bench PRIVATE FlexAID::VoronoiCF)
#
# See AGENTS.md / audit feedback for monolith reduction goals.

cmake_minimum_required(VERSION 3.28)

# ─────────────────────────────────────────────────────────────────────────────
# FlexAID::Core — the main engine (OBJECT library today for build speed)
# ─────────────────────────────────────────────────────────────────────────────
if(TARGET flexaid_core AND NOT TARGET FlexAID::Core)
    add_library(FlexAID::Core ALIAS flexaid_core)
    message(STATUS "FlexAIDComponents: FlexAID::Core alias created")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# FlexAID::VoronoiCF — header-oriented contact function batch interface
# (as suggested in audit for component split example)
# Pattern: create plain target then ALIAS with :: namespace (avoids reserved name issues).
# ─────────────────────────────────────────────────────────────────────────────
add_library(FlexAID_VoronoiCF INTERFACE)
target_include_directories(FlexAID_VoronoiCF INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB
)
target_sources(FlexAID_VoronoiCF INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/VoronoiCFBatch_SoA.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LIB/VoronoiCFBatch.h
)
target_link_libraries(FlexAID_VoronoiCF INTERFACE
    Eigen3::Eigen
)
if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(FlexAID_VoronoiCF INTERFACE OpenMP::OpenMP_CXX)
endif()
target_compile_definitions(FlexAID_VoronoiCF INTERFACE FLEXAIDS_HAS_VORONOI_CF)

add_library(FlexAID::VoronoiCF ALIAS FlexAID_VoronoiCF)
message(STATUS "FlexAIDComponents: FlexAID::VoronoiCF INTERFACE target ready")

# ─────────────────────────────────────────────────────────────────────────────
# FlexAID::StatMech (placeholder for future finer split)
# ─────────────────────────────────────────────────────────────────────────────
add_library(FlexAID_StatMech INTERFACE)
target_link_libraries(FlexAID_StatMech INTERFACE FlexAID::Core)
target_compile_definitions(FlexAID_StatMech INTERFACE FLEXAIDS_HAS_STATMECH)
add_library(FlexAID::StatMech ALIAS FlexAID_StatMech)
message(STATUS "FlexAIDComponents: FlexAID::StatMech (thin alias to Core) ready")

# ─────────────────────────────────────────────────────────────────────────────
# FlexAID::ENCoM (vibrational entropy component hook)
# ─────────────────────────────────────────────────────────────────────────────
add_library(FlexAID_ENCoM INTERFACE)
target_link_libraries(FlexAID_ENCoM INTERFACE FlexAID::Core)
target_compile_definitions(FlexAID_ENCoM INTERFACE FLEXAIDS_HAS_ENCOM)
add_library(FlexAID::ENCoM ALIAS FlexAID_ENCoM)
message(STATUS "FlexAIDComponents: FlexAID::ENCoM ready")

# Future components can follow the same pattern:
#  - tENCoM
#  - ShannonThermo
#  - GrandCanonical (see FLEXAIDS_GRAND_CANONICAL option + TargetServer/MultiSiteGPF)
#  - etc.
# Note: Grand canonical sources are already in flexaid_core when FLEXAIDS_GRAND_CANONICAL=ON (P0).

# Export a variable so root / other modules know components were loaded
set(FLEXAIDS_COMPONENTS_LOADED TRUE)
