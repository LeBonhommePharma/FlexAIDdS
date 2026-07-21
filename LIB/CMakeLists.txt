# LIB/CMakeLists.txt — FlexAID core library
#
# Defines the flexaid_core OBJECT library containing all core C++ sources.
# GPU-specific sources (CUDA, Metal, ROCm) are added conditionally.
#
# Source list is in a variable to make maintenance easier and support
# future component splits (see audit monolith reduction).

set(FLEXAID_CORE_SOURCES
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
# FLEXAIDDS_CMAES_CMAKE_BEGIN
    cmaes_search.cpp
# FLEXAIDDS_CMAES_CMAKE_END
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
    receptor_prep.cpp
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
    MinibatchSampler.cpp
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
    # ── Level-3 H(ω) ligand vibrational-mode Shannon entropy (diagnostic) ──
    VibEntropy.cpp
    # ── DiFT torsional Fourier parametrization ──
    DiFT/DiFT.cpp
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
    # ── Shannon Thermodynamic GA (entropy-collapse driver) ──
    shannon_ga.cpp
    LigandRingFlex/RingConformerLibrary.cpp
    LigandRingFlex/SugarPucker.cpp
    LigandRingFlex/LigandRingFlex.cpp
    NATURaL/NATURaLDualAssembly.cpp
    NATURaL/RibosomeElongation.cpp
    NATURaL/TransloconInsertion.cpp
    NATURaL/NucleationDetector.cpp
    NATURaL/NascentChainScheduler.cpp
    NATURaL/FibrilGrowthOracle.cpp
    NATURaL/DualAssemblyRunner.cpp
    ChiralCenter/ChiralCenterGene.cpp
    # ── CleftDetector + Mol2Reader + SdfReader ──
    CleftDetector.cpp
    DirectLigandIC.h
    Mol2Reader.cpp
    SdfReader.cpp
    # ── CavityDetect ──
    CavityDetect/CavityDetect.cpp
    # ── PTMAttachment ──
    PTMAttachment/PTMAttachment.cpp
    # ── JSON config parser ──
    config_parser.cpp
    # ── Typed protocol knobs (env adapter; audit getenv consolidation) ──
    ProtocolConfig.cpp
    # ── Campaign RUN_RECEIPT.json provenance (DatasetRunner + C0 schema) ──
    RunReceipt.cpp
    # ── ThermodynamicEngine: G_bind decomposition (disabled by default) ──
    ThermodynamicEngine.cpp
    # ── Coarse pocket-scan gen-0 seeding ──
    coarse_init.cpp
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
    MultiSiteGPF.cpp
    TargetKnowledgeBase.cpp
    TargetServer.cpp
    # ── MPI transport (stubs when MPI is off) ──
    MPITransport.cpp
    # ── InStreamClustering: online medoid clustering during GA ──
    InStreamClustering.cpp
    # ── P0 source validator completion: previously orphaned core units ──
    BinarySnapshot.cpp
    build_close.cpp
    cffunction.cpp
    cpu_eval.cpp
    cpu_fast_optics.cpp
    dee_print.cpp
    direct_input.cpp
    spfunction.cpp
    tENCoM/tencm_cpu_fallback.cpp
    # ── Native-pose CF scoring diagnostic (DatasetRunner IdealPopulation) ──
    native_score.cpp
)

add_library(flexaid_core OBJECT ${FLEXAID_CORE_SOURCES})

# Export the source list for use by parent scope / other modules if needed
# (e.g. for documentation or future finer-grained components).
set(FLEXAID_CORE_SOURCES "${FLEXAID_CORE_SOURCES}" PARENT_SCOPE)

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
    CXX_STANDARD ${CMAKE_CXX_STANDARD}   # follows global (capped to 20 on MSVC)
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
    # -fno-finite-math-only: keep -ffast-math optimizations but preserve
    # infinity()/NaN semantics used as sentinels in fast_optics, BindingMode,
    # SharedPosePool, GrandPartitionFunction, UnifiedHardwareDispatch, etc.
    target_compile_options(flexaid_core PRIVATE
        -Wall -O3 -ffast-math -fno-finite-math-only)

    # -mcpu=native (opt-in, FLEXAIDS_MCPU_NATIVE — see cmake/FlexAIDOptions.cmake):
    # tunes flexaid_core codegen for the exact build machine's CPU. Scoped to
    # this OBJECT library only. Not portable across machines; may shift FP
    # codegen (FMA contraction, vector width) under the -ffast-math already
    # enabled above, so a build with this ON is not guaranteed bit-identical
    # to one with it OFF. Default ON here matches local perf-tuning use;
    # portable/reproducible builds should configure with
    # -DFLEXAIDS_MCPU_NATIVE=OFF.
    if(FLEXAIDS_MCPU_NATIVE)
        target_compile_options(flexaid_core PRIVATE
            $<$<AND:$<CXX_COMPILER_ID:AppleClang,Clang>,$<STREQUAL:${CMAKE_SYSTEM_PROCESSOR},arm64>>:-mcpu=native>)
        message(STATUS "flexaid_core: -mcpu=native enabled (FLEXAIDS_MCPU_NATIVE=ON) — build is tuned for THIS machine and is not portable; use -DFLEXAIDS_MCPU_NATIVE=OFF for portable/reproducible builds")
    endif()
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

# libm symbols must resolve on Linux, especially under LTO where the linker
# does not auto-pull libm.
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(flexaid_core PUBLIC m)
endif()

# ─── Metal (optional) ───────────────────────────────────────────────────
# Propagate FLEXAIDS_USE_METAL so hardware_detect.cpp / UnifiedHardwareDispatch
# / gaboom / FOPTICS compile their Metal call sites. The OBJCXX bridges
# (metal_eval.mm, MetalRMSDBridge.mm, …) and Metal frameworks are attached to
# flexaid_core in the root CMakeLists.txt after this subdirectory returns —
# see "Metal sources (attach to flexaid_core once)". Do NOT leave Metal
# implementations only on executables: that breaks multi-target and LTO links
# (Homebrew --with-metal undefined-symbol failures).
if(FLEXAIDS_USE_METAL)
    target_compile_definitions(flexaid_core PUBLIC FLEXAIDS_USE_METAL)
endif()

# ─── SIMD configuration ─────────────────────────────────────────────────
flexaids_configure_simd(flexaid_core)

# NOTE: FlexAID::Core alias is created in cmake/FlexAIDComponents.cmake (after subdir)
# to avoid duplicate target definition.
