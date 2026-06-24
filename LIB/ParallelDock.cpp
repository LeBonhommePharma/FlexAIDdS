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
    const ParallelDockConfig& config,
    genlim* parent_gene_lim,
    int     parent_num_genes)
    : FA_(FA), GB_(GB), VC_(VC),
      atoms_(atoms), residue_(residue),
      cleftgrid_(cleftgrid),
      config_(config),
      parent_gene_lim_(parent_gene_lim),
      parent_num_genes_(parent_num_genes),
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

    // Store grid index mapping for later chromosome remapping
    result.region_grid_indices = region.grid_indices;

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

    // ── gene[0]: translation — limited to this sub-region's grid ───────────
    gene_lim[0].min = 1.0;
    gene_lim[0].max = (double)(sub_num_grd - 1);
    gene_lim[0].del = gene_lim[0].max - gene_lim[0].min;
    set_bins(&gene_lim[0]);

    // ── genes[1..N-1]: rotation + flexible bonds — copy from parent ─────────
    // Previously these were left uninitialised (the "simplified version" TODO).
    // Uninitialised gene limits cause GA() to search the full [0,1] cube for
    // dihedral angles and produce physically nonsensical conformations.
    if (parent_gene_lim_ && parent_num_genes_ >= num_genes) {
        for (int g = 1; g < num_genes; g++) {
            gene_lim[g] = parent_gene_lim_[g];
        }
    } else {
        // Fallback when parent limits aren't provided (rigid-body docking only).
        // Rotational genes are typically [-π, π] in IC space → map to [0, 1].
        for (int g = 1; g < num_genes; g++) {
            gene_lim[g].min = 0.0;
            gene_lim[g].max = 1.0;
            gene_lim[g].del = 1.0;
            set_bins(&gene_lim[g]);
        }
        if (num_genes > 1) {
            fprintf(stderr,
                "ParallelDock WARNING: region %d using fallback gene limits "
                "for genes 1..%d — pass parent_gene_lim for correct flex poses\n",
                region.region_id, num_genes - 1);
        }
    }

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

    // Collect results: snapshot energies for partition function + best chromosome
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

            // Record grid position
            int grd_idx = (int)chrom[i].genes[0].to_ic;
            result.best_local_grd_idx = grd_idx;
            if (grd_idx >= 0 && grd_idx < sub_num_grd) {
                std::memcpy(result.best_coor, subgrid[grd_idx].coor, 3 * sizeof(float));
            }

            // ── Deep-copy the full best chromosome ──────────────────────────
            // This is the fix: store genes[0..num_genes-1] so that top.cpp
            // can populate chrom[0] and feed it into the standard
            // clustering / write_rrd output pipeline.
            result.best_genes.assign(chrom[i].genes, chrom[i].genes + num_genes);

            // Copy all chromosome metadata (evalue, cf, fitnes, status,
            // boltzmann_weight, free_energy, ring flex arrays).
            // Leave .genes = nullptr — we own the data in result.best_genes.
            result.best_chrom_meta = chrom[i];
            result.best_chrom_meta.genes = nullptr;
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
// Best chromosome extraction — called by top.cpp after run()
// ============================================================================

bool ParallelDockManager::get_best_chromosome(chromosome& out_chrom,
                                               int& out_global_grd_idx) const
{
    // Find region with lowest best_energy
    const RegionResult* best_region = nullptr;
    for (const auto& rr : results_) {
        if (rr.best_genes.empty()) continue;
        if (!best_region || rr.best_energy < best_region->best_energy)
            best_region = &rr;
    }

    if (!best_region || best_region->best_genes.empty()) {
        fprintf(stderr, "ParallelDock: get_best_chromosome: no valid region result\n");
        out_global_grd_idx = -1;
        return false;
    }

    int num_genes = (int)best_region->best_genes.size();

    // Validate that out_chrom.genes is allocated (caller's responsibility)
    if (!out_chrom.genes) {
        fprintf(stderr, "ParallelDock: get_best_chromosome: out_chrom.genes is NULL\n");
        out_global_grd_idx = -1;
        return false;
    }

    // Copy chromosome metadata
    const chromosome& meta = best_region->best_chrom_meta;
    out_chrom.evalue          = meta.evalue;
    out_chrom.app_evalue      = meta.app_evalue;
    out_chrom.cf              = meta.cf;
    out_chrom.fitnes          = meta.fitnes;
    out_chrom.status          = 'n';   // force "evaluated" so clustering uses it
    out_chrom.boltzmann_weight = meta.boltzmann_weight;
    out_chrom.free_energy     = meta.free_energy;
    std::memcpy(out_chrom.ring_phases, meta.ring_phases, sizeof(out_chrom.ring_phases));
    std::memcpy(out_chrom.ring_six,    meta.ring_six,    sizeof(out_chrom.ring_six));

    // Copy genes[0..num_genes-1]
    for (int g = 0; g < num_genes; g++) {
        out_chrom.genes[g] = best_region->best_genes[g];
    }

    // ── Remap gene[0]: local subgrid index (1-based) → global cleftgrid index ──
    // GridDecomposer::extract_subgrid() places original cleftgrid[0] at subgrid[0]
    // and region.grid_indices[k-1] at subgrid[k] for k = 1..sub_num_grd-1.
    int local_idx = best_region->best_local_grd_idx;
    out_global_grd_idx = -1;

    if (local_idx >= 1 &&
        local_idx <= (int)best_region->region_grid_indices.size()) {
        out_global_grd_idx = best_region->region_grid_indices[local_idx - 1];
    } else {
        fprintf(stderr,
            "ParallelDock: get_best_chromosome: local_idx %d out of range "
            "(region %d has %d grid points)\n",
            local_idx, best_region->region_id,
            (int)best_region->region_grid_indices.size());
        return false;
    }

    // Apply the global index to gene[0] so ic2cf() uses the correct grid point
    out_chrom.genes[0].to_ic = (double)out_global_grd_idx;

    printf("ParallelDock: best pose — region %d, local_idx %d → global_idx %d, "
           "evalue=%.3f\n",
           best_region->region_id, local_idx, out_global_grd_idx,
           out_chrom.evalue);

    return true;
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
