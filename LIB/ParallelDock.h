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

    // ── Best chromosome storage for downstream PDB output ──────────────────
    // Deep copy of the best chromosome's genes array (indexed 0..num_genes-1).
    // After run(), get_best_chromosome() uses this to populate chrom[0] in
    // top.cpp so that the standard clustering/output path produces a real PDB.
    std::vector<gene> best_genes;       // genes[0..num_genes-1] of best pose
    chromosome        best_chrom_meta;  // evalue, cf, fitnes, etc. (genes ptr = NULL)
    int               best_local_grd_idx;  // local subgrid index (1-based) for gene[0]

    // Grid index mapping: best_local_grd_idx k → global index region_grid_indices[k-1]
    // Copied from GridRegion::grid_indices at run_region() time.
    std::vector<int>  region_grid_indices;

    RegionResult() : region_id(-1), best_energy(1e30), best_coor{0,0,0},
                     num_snapshots(0), best_local_grd_idx(-1) {
        std::memset(&best_chrom_meta, 0, sizeof(best_chrom_meta));
    }
};

class ParallelDockManager {
public:
    // parent_gene_lim: the genlim array from top.cpp covering genes 0..parent_num_genes-1.
    //   genes[1..N-1] (flexible bonds, rotations) are COPIED into each region's workspace
    //   so that GA() runs with correct gene bounds.  gene[0] is overridden per subgrid.
    //   Pass nullptr to fall back to default [0,1] limits (produces wrong poses — only for
    //   rigid-body docking where gene[0] is the only meaningful gene).
    ParallelDockManager(
        FA_Global* FA, GB_Global* GB, VC_Global* VC,
        atom* atoms, resid* residue,
        gridpoint* cleftgrid,
        const ParallelDockConfig& config,
        genlim* parent_gene_lim = nullptr,
        int     parent_num_genes = 0
    );

    // Phase 1: Decompose grid into octree regions
    void decompose();

    // Phase 2: Run all GA instances
    //   - MPI mode: distributed across ranks (call from all ranks)
    //   - Thread mode: OpenMP parallel over regions on single machine
    void run(cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));

    // Phase 3: Aggregate results into global partition function
    statmech::Thermodynamics aggregate() const;

    // ── Best chromosome extraction ─────────────────────────────────────────
    // Finds the region with lowest best_energy and copies its best chromosome
    // into out_chrom (caller must have pre-allocated out_chrom.genes with at
    // least num_genes entries).  Also writes the GLOBAL cleftgrid index to
    // out_chrom.genes[0].to_ic so that ic2cf() operates on the correct point.
    //
    // Returns true on success, false if no valid chromosome was found.
    // Called by top.cpp after run() to feed chrom[0] into the standard
    // clustering/write_rrd output path.
    bool get_best_chromosome(chromosome& out_chrom, int& out_global_grd_idx) const;

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

    // Parent gene limits are accepted by the constructor for compatibility
    // with callers that still pass them, but region workspaces now inherit
    // their mutable state from the copied globals.
    [[maybe_unused]] genlim* parent_gene_lim_;
    [[maybe_unused]] int     parent_num_genes_;

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
