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

// Allocate Metal device buffers and compile the compute shader.
//   n_atoms        – total atom count
//   n_types        – number of atom types (energy_matrix dimension)
//   max_pop        – maximum population size (upper bound on batch)
//   lig_first      – 0-based index of first ligand atom
//   lig_last       – 0-based index of last ligand atom
//   perm           – van-der-Waals permeability (FA->permeability)
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

#endif  // FLEXAIDS_USE_METAL
