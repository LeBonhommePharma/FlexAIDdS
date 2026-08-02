// cuda_eval.cuh — CUDA kernel interface for batched chromosome evaluation
//
// When FLEXAIDS_USE_CUDA is defined (and nvcc compiles this TU),
// FlexAIDdS offloads the per-chromosome complementarity-function (CF)
// evaluation to the GPU.
//
// FULL-FIDELITY scoring (P4).  The GPU CF now reproduces the CPU
// ic2cf/vcfunction CF *assembly* within a drift tolerance sufficient to
// preserve pose ranking (NOT bit-for-bit).  The GPU fills every scoring
// channel that feeds get_cf_evalue():
//
//     com  + wal + sas + con + elec + hbond + gist_desolv + pb_clash
//
//   – com : energy-matrix yval × normalised contact area, optionally
//           distance-weighted (exp(-r/r0)) to match the VCT dist-weight term.
//   – wal : KWALL × (r⁻¹² − (perm·rAB)⁻¹²) repulsive wall (unchanged).
//   – sas : per-ligand-atom residual solvent-accessible surface scored
//           against the solvent column (T-1) of the energy matrix.
//   – con : covalent-constraint Gaussian restraint (KDIST baseline minus
//           KDIST·GetValueFromGaussian(dist,bond_len,max_dist)).
//   – elec: distance-dependent-dielectric Coulomb over contacting pairs.
//   – hbond: distance-Gaussian × representative angle term, gated by
//            donor/acceptor complementarity or salt bridge (see error note).
//   – gist_desolv: per-ligand-atom trilinear GIST grid desolvation.
//   – pb_clash: all-pairs ligand↔receptor PoseBusters vdW-overlap penalty
//               plus the pocket-presence centroid penalty.
//
// ── Contact-area model (drift-tolerant; the recommended plan route (b)) ──────
// The branchy analytic Voronoi (voronoi_poly2/calc_areas) is deliberately NOT
// ported.  Instead each ligand↔receptor pair contributes a normalised contact
// area from a C0 linear switching function:
//
//     rel_area(r) = 1                              for r <= rA+rB
//                 = 1 - (r-(rA+rB))/(2·Rw)         for rA+rB < r < rA+rB+2·Rw
//                 = 0                              for r >= rA+rB+2·Rw
//
// This maps one pair → one thread with no divergent polyhedron clipping, keeps
// the receptor + energy matrix resident on-device, and batches the whole
// population per generation.  Expected error vs. the analytic Voronoi area is a
// smooth per-contact bias (the linear ramp over-/under-estimates the true
// spherical-cap overlap by up to ~15–20 % on individual faces) that is highly
// correlated across poses, so it largely cancels in pose *ranking*; the drift
// policy (ranking-preserving, not bit-exact) is what makes it admissible.
// Terms that depend on true absolute geometry (angular H-bond directionality,
// metal-CN, vct-entropy, tENCoM) are NOT reproduced here — see the gaboom.cpp
// divergence guard.
//
// Context lifetime (persistent across generate() calls):
//   cuda_eval_init()      – once at GA startup; uploads rigid atom data
//   cuda_eval_set_extra() – once; uploads pb-vdw / charge / hbond-flag /
//                           constraint / GIST device arrays (may pass nullptr)
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

// ── Shared GPU-CF POD types ──────────────────────────────────────────────────
// Defined here AND (identically) in metal_eval.h.  gaboom.cpp may include both
// headers under their respective backend guards, so a shared include-guard macro
// prevents a duplicate-definition error while keeping every symbol inside the
// P4 track's assigned files (no new header introduced).
#ifndef FLEXAIDS_GPU_CF_TYPES_DEFINED
#define FLEXAIDS_GPU_CF_TYPES_DEFINED

// Per-batch scalar parameters (weights/config; fixed per dock, cheap to pass).
struct GpuCfParams {
    float dw_r0;             // VCT distance-weight r0 (<=0 => distance weighting OFF)
    int   elec_on;           // 1 => compute Coulomb elec term
    float dielectric;        // FA->dielectric
    float hbond_weight;      // 0 => H-bond term OFF (else standard H-bond weight)
    float hbond_salt_weight; // salt-bridge weight
    float hbond_opt_dist;    // optimal donor..acceptor distance
    float hbond_sigma_dist;  // distance Gaussian sigma
    float hbond_angle_repr;  // representative angle term in [0,1] (drift model)
    float pb_clash_weight;   // 0 => pb_clash term OFF
    float pb_clash_ratio;    // PoseBusters clash ratio
    float pb_clash_exponent; // severity exponent p
    float pb_pocket_weight;  // 0 => pocket-presence penalty OFF
    float pb_pocket_radius;  // pocket radius (A)
    float kdist;             // KDIST constant for covalent constraints
    float sas_weight;        // FA->sas_weight applied to the sas channel
    int   solvent_flat;      // 1 => FA->solventterm active (sas = solventterm*SAS)
    float solventterm;       // FA->solventterm (only used when solvent_flat)
};

// Host output channels, one double array of length pop_size each.
struct GpuCfResults {
    double* com;
    double* wal;
    double* sas;
    double* con;
    double* elec;
    double* hbond;
    double* gist_desolv;
    double* pb_clash;
};

// Rigid, once-per-dock device inputs beyond xyz/type/radius/emat.  Any pointer
// may be nullptr and the corresponding term is then disabled on-device.
struct GpuCfExtraStatic {
    const float* atom_pbvdw;   // [n_atoms] PoseBusters element vdW radius
    const float* atom_charge;  // [n_atoms] partial (or RESP) charge
    const int*   atom_hflags;  // [n_atoms] bit0 donor, bit1 acceptor, bit2 heavy(non-H)

    int          n_cons;       // number of covalent (type-1) constraints
    const int*   cons_i;       // [n_cons] 0-based atom index
    const int*   cons_j;       // [n_cons] 0-based atom index
    const float* cons_bondlen; // [n_cons]
    const float* cons_maxdist; // [n_cons]

    // GIST grid (gist_nx==0 => GIST disabled).  Layout: i*ny*nz + j*nz + k.
    int          gist_nx, gist_ny, gist_nz;
    float        gist_origin[3];
    float        gist_inv_delta[3]; // 1/delta[k][k]
    float        gist_weight;
    const float* gist_data;    // [nx*ny*nz]
};

#endif // FLEXAIDS_GPU_CF_TYPES_DEFINED

// Opaque handle to all device-resident data.
struct CudaEvalCtx;

// Allocate device memory and upload constant atom data.
//   n_atoms        – total atom count
//   n_types        – number of atom types (energy_matrix dimension)
//   max_pop        – maximum population size (upper bound on any batch)
//   max_genes      – maximum genes per chromosome (buffer bound)
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

// Upload the optional full-fidelity static arrays (pb-vdw / charge / hbond
// flags / constraints / GIST).  Call once after cuda_eval_init.  Passing
// extra==nullptr (or leaving individual pointers null) disables the matching
// terms; the base com/wal/sas path is unaffected.
void cuda_eval_set_extra(CudaEvalCtx* ctx, const GpuCfExtraStatic* extra);

// Evaluate a batch of chromosomes on the GPU (full-fidelity CF).
//   ctx      – context from cuda_eval_init
//   pop_size – number of chromosomes to evaluate this call
//   n_genes  – genes per chromosome
//   h_genes  – host gene array [pop_size × n_genes, double]
//   params   – per-batch scalar weights/config (must be non-null)
//   out      – host output channels; each array has length >= pop_size
void cuda_eval_batch(CudaEvalCtx*        ctx,
                     int                 pop_size,
                     int                 n_genes,
                     const double*       h_genes,
                     const GpuCfParams*  params,
                     const GpuCfResults* out);

// Free all device memory.
void cuda_eval_shutdown(CudaEvalCtx* ctx);

#endif  // FLEXAIDS_USE_CUDA
