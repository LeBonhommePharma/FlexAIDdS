// cpu_eval.h — CPU fallback for batched chromosome fitness evaluation
//
// Provides the same API as cuda_eval.cuh / metal_eval.h but runs entirely
// on the CPU using OpenMP thread parallelism.  Used when no GPU backend
// (CUDA, ROCm, Metal) is available.
//
// Scoring pipeline (identical to GPU kernels):
//   1. Decode translation genes (tx, ty, tz) from the gene vector.
//   2. For each ligand-protein atom pair:
//        a. Compute inter-atomic distance r.
//        b. Normalised contact area via linear switching function.
//        c. COM: energy-matrix lookup via linear interpolation.
//        d. WAL: soft_wall_fitness_energy_f (soft-core k·o² or legacy capped r⁻¹²).
//        e. SAS: subtract contact area from per-ligand-atom SAS counter.
//   3. SAS energy contribution from remaining exposed surface area.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#pragma once

#include <cstddef>

// Number of samples per type-pair energy curve (matches CUDA/Metal).
static constexpr int CPU_EMAT_SAMPLES = 128;

// Opaque handle to CPU evaluation context.
struct CpuEvalCtx;

// Allocate CPU evaluation context.
//   n_atoms        – total atom count
//   n_types        – number of atom types (energy_matrix dimension)
//   max_pop        – maximum population size
//   max_genes      – maximum genes per chromosome
//   lig_first      – 0-based index of first ligand atom
//   lig_last       – 0-based index of last ligand atom
//   perm           – van-der-Waals permeability
//   h_atom_xyz     – atom coordinates   [n_atoms × 3, float]
//   h_atom_type    – atom type array    [n_atoms, int, 0-based]
//   h_atom_radius  – atom radii         [n_atoms, float]
//   h_emat_sampled – pre-sampled energy matrix
//                    [n_types × n_types × CPU_EMAT_SAMPLES, float]
CpuEvalCtx* cpu_eval_init(int   n_atoms,
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

// Evaluate a batch of chromosomes on the CPU.
//   ctx        – context from cpu_eval_init
//   pop_size   – number of chromosomes to evaluate
//   n_genes    – genes per chromosome
//   h_genes    – gene array [pop_size × n_genes, double]
//   h_com_out  – output: complementarity CF   [pop_size, double]
//   h_wal_out  – output: wall/clash energy     [pop_size, double]
//   h_sas_out  – output: solvent-accessible    [pop_size, double]
void cpu_eval_batch(CpuEvalCtx*  ctx,
                    int           pop_size,
                    int           n_genes,
                    const double* h_genes,
                    double*       h_com_out,
                    double*       h_wal_out,
                    double*       h_sas_out);

// Free CPU evaluation context.
void cpu_eval_shutdown(CpuEvalCtx* ctx);

// Configure soft-core wall (default soft_wall_cutoff=0.40, k_wal=50).
// Must match FA->soft_wall_cutoff / FA->k_wal_stiff for parity with vcfunction.
void cpu_eval_set_soft_wall(CpuEvalCtx* ctx, float soft_wall_cutoff, float k_wal);

