// ParallelDock.h — Orchestrator for massively parallel grid-decomposed docking
//
// Combines octree spatial decomposition of the cube grid with independent
// GA instances per region. Results are aggregated via StatMechEngine to
// compute the global partition function and thermodynamic properties.
#pragma once

#include "flexaid.h"
#include "gaboom.h"
#include "GAContext.h"
#include "statmech.h"
#include "GridDecomposer.h"
#include "SharedPosePool.h"
#include <vector>
#include <functional>

struct ParallelDockConfig {
    int target_regions       = 128;   // number of spatial regions
    int min_points_per_region = 50;   // merge regions smaller than this
    int pose_pool_size       = 256;   // shared pool capacity
    int exchange_interval    = 10;    // generations between pool reads
    int seed_from_pool_count = 5;     // how many pool poses to inject per exchange
    bool use_mpi             = false; // true for distributed (MPI), false for thread-based
};

struct RegionResult {
    int region_id;
    statmech::Thermodynamics local_thermo;
    std::vector<double> energies;
    std::vector<double> multiplicities;  // double to match StatMechEngine API (post C-1)
    double best_energy;
    float  best_coor[3];
    int    num_snapshots;

    RegionResult() : region_id(-1), best_energy(1e30), best_coor{0,0,0}, num_snapshots(0) {}
};

class ParallelDockManager {
public:
    ParallelDockManager(
        FA_Global* FA, GB_Global* GB, VC_Global* VC,
        atom* atoms, resid* residue,
        gridpoint* cleftgrid,
        const ParallelDockConfig& config
    );

    // Phase 1: Decompose grid into octree regions
    void decompose();

    // Phase 2: Run all GA instances
    //   - MPI mode: distributed across ranks (call from all ranks)
    //   - Thread mode: OpenMP parallel over regions on single machine
    void run(cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));

    // Phase 3: Aggregate results into global partition function
    statmech::Thermodynamics aggregate() const;

    // Access per-region results
    const std::vector<RegionResult>& region_results() const { return results_; }

    // Build merged StatMechEngine from all regions
    statmech::StatMechEngine get_global_engine() const;

    // Access regions (for inspection/visualization)
    const std::vector<GridRegion>& regions() const { return regions_; }

private:
    FA_Global* FA_;
    GB_Global* GB_;
    VC_Global* VC_;
    atom* atoms_;
    resid* residue_;
    gridpoint* cleftgrid_;
    ParallelDockConfig config_;

    std::vector<GridRegion> regions_;
    std::vector<RegionResult> results_;
    SharedPosePool pool_;

    // Run a single region's GA and return its result
    RegionResult run_region(
        const GridRegion& region,
        unsigned int rng_seed,
        cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*)
    );

    // Create deep copies of mutable state for a region
    struct RegionWorkspace {
        FA_Global fa;
        GB_Global gb;
        VC_Global vc;
        std::vector<atom> atoms_copy;
        std::vector<resid> residue_copy;
        GAContext ga_ctx;  // per-region GA state for re-entrant execution
    };
    RegionWorkspace create_workspace() const;
};
