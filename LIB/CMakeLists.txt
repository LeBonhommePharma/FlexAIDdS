# LIB/CMakeLists.txt — FlexAID core library
#
# Defines the flexaid_core OBJECT library containing all core C++ sources.
# GPU-specific sources (CUDA, Metal, ROCm) are added conditionally.

add_library(flexaid_core OBJECT
    top.cpp
    fileio.cpp
    assign_shift.cpp
    buildic.cpp
    Vcontacts.cpp
    ic2cf.cpp
    vcfunction.cpp
    create_rebuild_list.cpp
    read_coor.cpp
    read_constraints.cpp
    assign_constraint.cpp
    modify_pdb.cpp
    update_constraint.cpp
    update_optres.cpp
    update_bonded.cpp
    read_grid.cpp
    read_input.cpp
    read_lig.cpp
    read_pdb.cpp
    add2_optimiz_vec.cpp
    assign_radius.cpp
    assign_radii.cpp
    assign_radii_types.cpp
    assign_types.cpp
    bondedlist.cpp
    buildcc.cpp
    buildlist.cpp
    geometry.cpp
    read_conect.cpp
    residue_conect.cpp
    wif083.cpp
    calc_cleftic.cpp
    write_pdb.cpp
    gaboom.cpp
    GAContext.cpp
    calc_center.cpp
    ic_bounds.cpp
    buildic_point.cpp
    buildcc_point.cpp
    read_rmsdst.cpp
    calc_rmsd_chrom.cpp
    calc_rmsd.cpp
    write_rrd.cpp
    write_sphere.cpp
    write_grid.cpp
    write_rrg.cpp
    read_rotlib.cpp
    number_of_dihedrals.cpp
    read_flexscfile.cpp
    read_rotobs.cpp
    set_intprob.cpp
    build_rotamers.cpp
    check_clash.cpp
    read_normalgrid.cpp
    read_eigen.cpp
    assign_eigen.cpp
    print_surfmat.cpp
    alter_mode.cpp
    dee_pivot.cpp
    read_emat.cpp
    read_spheres.cpp
    cluster.cpp
    DensityPeak_Cluster.cpp
    FastOPTICS_cluster.cpp
    FOPTICS.cpp
    fast_optics.cpp
    BindingMode.cpp
    rna_structure.cpp
    shortest_path.cpp
    assign_shortflex.cpp
    slice_grid.cpp
    partition_grid.cpp
    generate_grid.cpp
    maps.cpp
    # ── New C++20 modules ──
    tENCoM/tencm.cpp
    statmech.cpp
    encom.cpp
    # ── GIST water displacement + directional H-bond scoring ──
    GISTEvaluator.cpp
    HBondEvaluator.cpp
    # ── TurboQuant vector quantization ──
    TurboQuant.cpp
    # ── ProcessLigand ──
    ProcessLigand/ProcessLigand.cpp
    ProcessLigand/SmilesParser.cpp
    ProcessLigand/RingPerception.cpp
    ProcessLigand/Aromaticity.cpp
    ProcessLigand/RotatableBonds.cpp
    ProcessLigand/ValenceChecker.cpp
    ProcessLigand/SybylTyper.cpp
    ProcessLigand/FlexAIDWriter.cpp
    ProcessLigand/ROCmDispatch.cpp
    ProcessLigand/CoordBuilder.cpp
    CifReader.cpp
    MultiModelDock.cpp
    LibrarySplitter.cpp
    ReferenceEntropy.cpp
    ParallelCampaign.cpp
    # ── Hardware acceleration modules ──
    hardware_detect.cpp
    UnifiedHardwareDispatch.cpp
    ShannonThermoStack/ShannonThermoStack.cpp
    LigandRingFlex/RingConformerLibrary.cpp
    LigandRingFlex/SugarPucker.cpp
    LigandRingFlex/LigandRingFlex.cpp
    NATURaL/NATURaLDualAssembly.cpp
    NATURaL/RibosomeElongation.cpp
    NATURaL/TransloconInsertion.cpp
    NATURaL/NucleationDetector.cpp
    ChiralCenter/ChiralCenterGene.cpp
    # ── CleftDetector + Mol2Reader + SdfReader ──
    CleftDetector.cpp
    Mol2Reader.cpp
    SdfReader.cpp
    # ── CavityDetect ──
    CavityDetect/CavityDetect.cpp
    # ── PTMAttachment ──
    PTMAttachment/PTMAttachment.cpp
    # ── JSON config parser ──
    config_parser.cpp
    # ── Unified runtime hardware dispatch ──
    UnifiedHardwareDispatch.cpp
    # ── Parallel grid-decomposed docking ──
    GridDecomposer.cpp
    SharedPosePool.cpp
    ParallelDock.cpp
    # ── CoarseScreen ──
    CoarseScreen.cpp
    TwoStageScreen.cpp
    # ── TargetServer ──
    TargetValidation.cpp
    GrandPartitionFunction.cpp
    TargetKnowledgeBase.cpp
    TargetServer.cpp
    # ── MPI transport (stubs when MPI is off) ──
    MPITransport.cpp
)

target_include_directories(flexaid_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/tENCoM
    ${CMAKE_CURRENT_SOURCE_DIR}/ShannonThermoStack
    ${CMAKE_CURRENT_SOURCE_DIR}/LigandRingFlex
    ${CMAKE_CURRENT_SOURCE_DIR}/NATURaL
    ${CMAKE_CURRENT_SOURCE_DIR}/ChiralCenter
    ${CMAKE_CURRENT_SOURCE_DIR}/CavityDetect
    ${CMAKE_CURRENT_SOURCE_DIR}/PTMAttachment
)

set_target_properties(flexaid_core PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

# ─── Compiler flags ──────────────────────────────────────────────────────
if(MSVC)
    target_compile_options(flexaid_core PRIVATE /W4 /O2 /MT)
    target_compile_definitions(flexaid_core PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _USE_MATH_DEFINES
        NOMINMAX
    )
else()
    target_compile_options(flexaid_core PRIVATE -Wall -O3 -ffast-math)
endif()

# ─── Required dependencies ───────────────────────────────────────────────
target_link_libraries(flexaid_core PUBLIC Eigen3::Eigen)
target_compile_definitions(flexaid_core PUBLIC FLEXAIDS_HAS_EIGEN)

# ─── OpenMP (optional) ──────────────────────────────────────────────────
if(FLEXAIDS_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(flexaid_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# ─── MPI (optional) ─────────────────────────────────────────────────────
if(FLEXAIDS_USE_MPI)
    target_link_libraries(flexaid_core PUBLIC MPI::MPI_CXX)
    target_compile_definitions(flexaid_core PUBLIC FLEXAIDS_USE_MPI)
endif()

# ─── SIMD configuration ─────────────────────────────────────────────────
flexaids_configure_simd(flexaid_core)
