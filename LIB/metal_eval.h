// metal_eval.h — Metal GPU evaluation API for batched chromosome scoring
//
// Mirrors cuda_eval.cuh but targets Apple GPU via Metal compute pipelines.
// Enabled with -DFLEXAIDS_USE_METAL.
//
// Usage:
//   MetalEvalCtx* ctx = metal_eval_init(...);   // once, before GA loop
//   metal_eval_batch(ctx, ...);                  // every generation
//   metal_eval_shutdown(ctx);                    // once, after GA loop
#pragma once

#ifdef FLEXAIDS_USE_METAL

#include <cstddef>

// Number of samples used to pre-sample each energy-matrix density curve.
// Must match the value in metal_eval.mm and the host-side sampling loop.
static constexpr int METAL_EMAT_SAMPLES = 128;

// Opaque handle to all Metal device-resident state.
struct MetalEvalCtx;

// Runtime availability probe. Compile-time Metal support does not guarantee
// that the current process can actually acquire a default Metal device.
bool metal_eval_runtime_available();

// Rich hardware capabilities (P2). Used by launcher for resource-aware decisions.
struct MetalCapabilities {
    bool   available;
    char   device_name[128];
    size_t unified_memory_bytes;   // 0 if unknown
    size_t max_buffer_length;
    int    gpu_core_estimate;      // rough (0 if unknown)
}; 

// Fills the struct. Safe to call even if Metal is not compiled in.
void metal_eval_get_capabilities(MetalCapabilities* out);

// Allocate Metal device buffers and compile the compute shader.
//   n_atoms        – total atom count
//   n_types        – number of atom types (energy_matrix dimension)
//   max_pop        – maximum population size (upper bound on batch)
//   lig_first      – 0-based index of first ligand atom
//   lig_last       – 0-based index of last ligand atom
//   perm           – van-der-Waals permeability (FA->permeability)
//   soft wall    – configure via metal_eval_set_soft_wall (soft_wall.h parity)
//   h_atom_xyz     – host atom coordinates   [n_atoms × 3, float]
//   h_atom_type    – host atom type array    [n_atoms, int, 0-based]
//   h_atom_radius  – host atom radii         [n_atoms, float]
//   h_emat_sampled – pre-sampled energy matrix
//                    [n_types × n_types × n_emat_samples, float]
//   n_emat_samples – samples per type-pair curve (must equal METAL_EMAT_SAMPLES)
MetalEvalCtx* metal_eval_init(int   n_atoms,
                               int   n_types,
                               int   max_pop,
                               int   lig_first,
                               int   lig_last,
                               float perm,
                               const float* h_atom_xyz,
                               const int*   h_atom_type,
                               const float* h_atom_radius,
                               const float* h_emat_sampled,
                               int   n_emat_samples);

// Evaluate a batch of chromosomes on the Metal GPU.
//   ctx        – context from metal_eval_init
//   pop_size   – number of chromosomes to evaluate this call
//   n_genes    – genes per chromosome
//   h_genes    – host gene array [pop_size × n_genes, double]
//   h_com_out  – host output: complementarity CF   [pop_size, double]
//   h_wal_out  – host output: wall/clash energy     [pop_size, double]
//   h_sas_out  – host output: solvent-accessible    [pop_size, double]
void metal_eval_batch(MetalEvalCtx* ctx,
                      int           pop_size,
                      int           n_genes,
                      const double* h_genes,
                      double*       h_com_out,
                      double*       h_wal_out,
                      double*       h_sas_out);

// Free all Metal device resources.
void metal_eval_shutdown(MetalEvalCtx* ctx);

// Soft-core wall config (default cutoff=0.40, k_wal=50). Parity with soft_wall.h.
void metal_eval_set_soft_wall(MetalEvalCtx* ctx, float soft_wall_cutoff, float k_wal);


// ─── Multi-complex batched evaluation (GPU utilization maximization) ──────────
//
// Dispatches N independent chromosome populations in a single Metal command buffer.
// Each entry i corresponds to one docking complex (worker); the GPU evaluates
// n_complex × pop_size chromosomes in one dispatch instead of N serial calls.
//
// Thread mapping: thread k → complex_id = k / pop_size, chrom_id = k % pop_size
//
// Requirements:
//   - All complexes must share the same n_atoms, n_types, and pop_size.
//   - The single shared context (ctx) must have been init'd with max_pop = n_complex × pop_size.
//   - h_genes[i]   points to the gene array for complex i  [pop_size × n_genes, double]
//   - h_com_out[i] / h_wal_out[i] / h_sas_out[i] are per-complex output buffers [pop_size]
//
// All buffers are host-side; the implementation copies in, dispatches, and copies out.
//
// Per-complex entry for multi-complex batch dispatch.
// Each entry carries its own atom data so complexes with different
// receptors/ligands can be batched into a single GPU kernel launch.
// When n_complex == 1 the fast path delegates to metal_eval_batch().
struct MetalMultiBatchEntry {
    // Gene input
    const double* h_genes;      // [pop_size × n_genes]
    // Per-complex atom data (different receptor+ligand each entry)
    const float*  h_atom_xyz;   // [n_atoms × 3]
    const int*    h_atom_type;  // [n_atoms]
    const float*  h_atom_radius;// [n_atoms]
    int           n_atoms;
    int           lig_first;
    int           lig_last;
    int           n_types;
    float         perm;
    // Output
    double*       h_com_out;    // [pop_size]
    double*       h_wal_out;    // [pop_size]
    double*       h_sas_out;    // [pop_size]
};

void metal_eval_batch_multi(MetalEvalCtx*              ctx,
                             int                        n_complex,
                             int                        pop_size,
                             int                        n_genes,
                             const MetalMultiBatchEntry* entries);

#endif  // FLEXAIDS_USE_METAL
