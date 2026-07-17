// cuda_eval.cuh — CUDA kernel interface for batched chromosome evaluation
//
// When FLEXAIDS_USE_CUDA is defined (and nvcc compiles this TU),
// FlexAIDdS offloads the per-chromosome complementarity-function (CF)
// evaluation to the GPU.
//
// Full-fidelity scoring (vs. the previous single-scalar energy table):
//   – Energy matrix sampled at N_CUDA_EMAT_SAMPLES points per type-pair
//     and interpolated on-device (gpu_get_yval).
//   – WAL (clash) term: soft-core k·o² or legacy capped r⁻¹² (soft_wall.h parity).
//   – SAS contribution: per-ligand-atom remaining exposed surface area,
//     accumulated via shared-memory atomic-add across protein contacts.
//   – COM term: energy-matrix yval × normalised contact area.
//
// Context lifetime (persistent across generate() calls):
//   cuda_eval_init()      – once at GA startup; uploads constant atom data
//   cuda_eval_batch()     – every generation (gene upload + kernel + readback)
//   cuda_eval_shutdown()  – once at GA teardown
//
// The caller (gaboom.cpp) maintains a static CudaEvalCtx* and re-initialises
// only when atom count or type count changes.
#pragma once

#ifdef FLEXAIDS_USE_CUDA

#include <cstddef>

// Samples per type-pair energy curve.  Must match N_EMAT_SAMPLES in cuda_eval.cu.
static constexpr int CUDA_EMAT_SAMPLES = 128;

// Opaque handle to all device-resident data.
struct CudaEvalCtx;

// Allocate device memory and upload constant atom data.
//   n_atoms        – total atom count
//   n_types        – number of atom types (energy_matrix dimension)
//   max_pop        – maximum population size (upper bound on any batch)
//   lig_first      – 0-based index of first ligand atom
//   lig_last       – 0-based index of last ligand atom
//   perm           – van-der-Waals permeability (FA->permeability)
//   h_atom_xyz     – host atom coordinates   [n_atoms × 3, float]
//   h_atom_type    – host atom type array    [n_atoms, int, 0-based]
//   h_atom_radius  – host atom radii         [n_atoms, float]
//   h_emat_sampled – pre-sampled energy-matrix density functions
//                    [n_types × n_types × CUDA_EMAT_SAMPLES, float]
CudaEvalCtx* cuda_eval_init(int   n_atoms,
                             int   n_types,
                             int   max_pop,
                             int   max_genes,
                             int   lig_first,
                             int   lig_last,
                             float perm,
                             const float* h_atom_xyz,
                             const int*   h_atom_type,
                             const float* h_atom_radius,
                             const float* h_emat_sampled);

// Evaluate a batch of chromosomes on the GPU.
//   ctx        – context from cuda_eval_init
//   pop_size   – number of chromosomes to evaluate this call
//   n_genes    – genes per chromosome
//   h_genes    – host gene array [pop_size × n_genes, double]
//   h_com_out  – host output: complementarity CF   [pop_size, double]
//   h_wal_out  – host output: wall/clash energy     [pop_size, double]
//   h_sas_out  – host output: solvent-accessible    [pop_size, double]
void cuda_eval_batch(CudaEvalCtx*  ctx,
                     int           pop_size,
                     int           n_genes,
                     const double* h_genes,
                     double*       h_com_out,
                     double*       h_wal_out,
                     double*       h_sas_out);

// Free all device memory.
void cuda_eval_shutdown(CudaEvalCtx* ctx);

// Soft-core wall config (default cutoff=0.40, k_wal=50). Parity with soft_wall.h.
void cuda_eval_set_soft_wall(CudaEvalCtx* ctx, float soft_wall_cutoff, float k_wal);


#endif  // FLEXAIDS_USE_CUDA
