// ParallelDock.cpp — Orchestrator for parallel grid-decomposed docking
#include "ParallelDock.h"
#include "fileio.h"

#include <cstdio>
#include <cstring>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef FLEXAIDS_USE_MPI
#include "MPITransport.h"
#endif

// ============================================================================
// Construction
// ============================================================================

ParallelDockManager::ParallelDockManager(
    FA_Global* FA, GB_Global* GB, VC_Global* VC,
    atom* atoms, resid* residue,
    gridpoint* cleftgrid,
    const ParallelDockConfig& config)
    : FA_(FA), GB_(GB), VC_(VC),
      atoms_(atoms), residue_(residue),
      cleftgrid_(cleftgrid),
      config_(config),
      pool_(config.pose_pool_size)
{}

// ============================================================================
// Phase 1: Grid decomposition
// ============================================================================

void ParallelDockManager::decompose() {
    regions_ = GridDecomposer::decompose_octree(
        cleftgrid_,
        FA_->num_grd,
        config_.target_regions,
        config_.min_points_per_region
    );

    printf("ParallelDock: decomposed into %d regions\n", (int)regions_.size());
}

// ============================================================================
// Phase 2: Run parallel GA instances
// ============================================================================

ParallelDockManager::RegionWorkspace ParallelDockManager::create_workspace() const {
    RegionWorkspace ws;

    // Shallow copy globals
    ws.fa = *FA_;
    ws.gb = *GB_;
    ws.vc = *VC_;

    // Deep copy mutable arrays
    ws.atoms_copy.assign(atoms_, atoms_ + FA_->atm_cnt);
    ws.residue_copy.assign(residue_, residue_ + FA_->res_cnt + 1);

    return ws;
}

void ParallelDockManager::run(
    cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*))
{
    if (regions_.empty()) {
        fprintf(stderr, "ParallelDock: no regions — call decompose() first\n");
        return;
    }

    int n_regions = (int)regions_.size();
    results_.resize(n_regions);

    // Seed generator for per-region RNG
    std::mt19937 seed_gen(42);

#ifdef FLEXAIDS_USE_MPI
    // MPI distributed mode: each rank processes a subset of regions
    int rank = MPITransport::rank();
    int world = MPITransport::world_size();

    // Round-robin assignment
    for (int r = rank; r < n_regions; r += world) {
        unsigned int seed = seed_gen() + r;
        results_[r] = run_region(regions_[r], seed, target);

        // Publish best to shared pool
        SharedPose sp;
        sp.energy = results_[r].best_energy;
        std::memcpy(sp.grid_coor, results_[r].best_coor, 3 * sizeof(float));
        sp.source_region = r;
        pool_.publish(sp);
    }

    // Gather all results to rank 0
    // (simplified: each rank sends its results)
    auto all_results = MPITransport::gather_results(results_, n_regions);
    if (rank == 0) results_ = std::move(all_results);

#else
    // Thread-based mode: OpenMP parallel over regions
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for (int r = 0; r < n_regions; r++) {
        unsigned int seed;
        #pragma omp critical
        { seed = seed_gen() + r; }

        results_[r] = run_region(regions_[r], seed, target);

        // Publish best to shared pool (thread-safe)
        SharedPose sp;
        sp.energy = results_[r].best_energy;
        std::memcpy(sp.grid_coor, results_[r].best_coor, 3 * sizeof(float));
        sp.source_region = r;
        pool_.publish(sp);
    }
#endif

    printf("ParallelDock: completed %d region GA runs\n", n_regions);
}

// ============================================================================
// Single region GA execution
// ============================================================================

RegionResult ParallelDockManager::run_region(
    const GridRegion& region,
    unsigned int rng_seed,
    cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*))
{
    RegionResult result;
    result.region_id = region.region_id;

    // Extract subgrid for this region
    int sub_num_grd;
    gridpoint* subgrid = GridDecomposer::extract_subgrid(
        cleftgrid_, region, sub_num_grd);
    if (!subgrid) {
        fprintf(stderr, "ParallelDock: failed to allocate subgrid for region %d\n",
                region.region_id);
        return result;
    }

    // Create per-region workspace (deep copy of mutable state)
    RegionWorkspace ws = create_workspace();

    // Override grid parameters for this region
    ws.fa.num_grd = sub_num_grd;

    // Allocate chromosomes and gene limits for this region's GA
    int num_genes = ws.gb.num_genes;
    int num_chrom = ws.gb.num_chrom;

    chromosome* chrom = (chromosome*)calloc(num_chrom * 2, sizeof(chromosome));
    chromosome* chrom_snapshot = (chromosome*)calloc(
        num_chrom * ws.gb.max_generations, sizeof(chromosome));

    if (!chrom || !chrom_snapshot) {
        free(subgrid);
        free(chrom);
        free(chrom_snapshot);
        return result;
    }

    // Allocate genes for each chromosome
    for (int i = 0; i < num_chrom * 2; i++) {
        chrom[i].genes = (gene*)calloc(num_genes, sizeof(gene));
    }
    for (int i = 0; i < num_chrom * ws.gb.max_generations; i++) {
        chrom_snapshot[i].genes = (gene*)calloc(num_genes, sizeof(gene));
    }

    genlim* gene_lim = (genlim*)calloc(num_genes, sizeof(genlim));

    // Set gene[0] limits to match subgrid size
    gene_lim[0].min = 1.0;
    gene_lim[0].max = (double)(sub_num_grd - 1);
    gene_lim[0].del = gene_lim[0].max - gene_lim[0].min;
    set_bins(&gene_lim[0]);

    // Copy remaining gene limits from original
    // (flexible bonds, rotamers, etc. are the same across regions)
    // gene_lim[1..num_genes-1] would need to be initialized from
    // the parent GA configuration — this is a simplified version

    int memchrom = num_chrom * 2;
    char gainpfile[256] = "";

    // Run the GA on this region's subgrid with per-region context
    GA(&ws.fa, &ws.gb, &ws.vc,
       &chrom, &chrom_snapshot,
       &gene_lim,
       ws.atoms_copy.data(),
       ws.residue_copy.data(),
       &subgrid,
       gainpfile,
       &memchrom,
       target,
       &ws.ga_ctx);

    // Collect results: snapshot energies for partition function
    int n_snap = ws.gb.num_chrom;  // snapshot from last generation
    result.best_energy = 1e30;

    statmech::StatMechEngine regional_engine(ws.fa.temperature);
    for (int i = 0; i < n_snap; i++) {
        double e = chrom[i].evalue;
        regional_engine.add_sample(e);
        result.energies.push_back(e);
        result.multiplicities.push_back(1);
        if (e < result.best_energy) {
            result.best_energy = e;
            int grd_idx = (int)chrom[i].genes[0].to_ic;
            if (grd_idx >= 0 && grd_idx < sub_num_grd) {
                std::memcpy(result.best_coor, subgrid[grd_idx].coor, 3 * sizeof(float));
            }
        }
    }
    result.num_snapshots = n_snap;
    result.local_thermo = regional_engine.compute();

    // Cleanup
    for (int i = 0; i < num_chrom * 2; i++) free(chrom[i].genes);
    for (int i = 0; i < num_chrom * ws.gb.max_generations; i++)
        free(chrom_snapshot[i].genes);
    free(chrom);
    free(chrom_snapshot);
    free(gene_lim);
    free(subgrid);

    return result;
}

// ============================================================================
// Phase 3: Aggregate results
// ============================================================================

statmech::StatMechEngine ParallelDockManager::get_global_engine() const {
    statmech::StatMechEngine global(FA_->temperature);

    for (const auto& r : results_) {
        if (r.energies.empty()) continue;
        global.merge_samples(
            std::span<const double>(r.energies),
            std::span<const double>(r.multiplicities)
        );
    }

    return global;
}

statmech::Thermodynamics ParallelDockManager::aggregate() const {
    auto engine = get_global_engine();

    if (engine.size() == 0) {
        statmech::Thermodynamics td{};
        td.temperature = FA_->temperature;
        return td;
    }

    auto td = engine.compute();

    printf("ParallelDock aggregate: %zu total samples across %d regions\n",
           engine.size(), (int)results_.size());
    printf("  F = %.4f kcal/mol, <E> = %.4f, S = %.6f kcal/mol/K\n",
           td.free_energy, td.mean_energy, td.entropy);

    return td;
}
