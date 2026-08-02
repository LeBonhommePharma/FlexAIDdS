// metal_eval.h — Metal GPU evaluation API for batched chromosome scoring
//
// Mirrors cuda_eval.cuh but targets Apple GPU via Metal compute pipelines.
// Enabled with -DFLEXAIDS_USE_METAL.
//
// FULL-FIDELITY scoring (P4): the Metal CF now fills the same eight scoring
// channels as the CUDA path — com, wal, sas, con, elec, hbond, gist_desolv,
// pb_clash — using the identical drift-tolerant linear-switching contact-area
// model (see cuda_eval.cuh for the model definition and error bound).
//
// Usage:
//   MetalEvalCtx* ctx = metal_eval_init(...);   // once, before GA loop
//   metal_eval_set_extra(ctx, &extra);          // once (optional full-fidelity)
//   metal_eval_batch(ctx, ...);                  // every generation
//   metal_eval_shutdown(ctx);                    // once, after GA loop
#pragma once

#ifdef FLEXAIDS_USE_METAL

#include <cstddef>

// Number of samples used to pre-sample each energy-matrix density curve.
// Must match the value in metal_eval.mm and the host-side sampling loop.
static constexpr int METAL_EMAT_SAMPLES = 128;

// ── Shared GPU-CF POD types ──────────────────────────────────────────────────
// Identical to the block in cuda_eval.cuh; guarded so including both headers in
// gaboom.cpp does not redefine the structs.
#ifndef FLEXAIDS_GPU_CF_TYPES_DEFINED
#define FLEXAIDS_GPU_CF_TYPES_DEFINED

struct GpuCfParams {
    float dw_r0;
    int   elec_on;
    float dielectric;
    float hbond_weight;
    float hbond_salt_weight;
    float hbond_opt_dist;
    float hbond_sigma_dist;
    float hbond_angle_repr;
    float pb_clash_weight;
    float pb_clash_ratio;
    float pb_clash_exponent;
    float pb_pocket_weight;
    float pb_pocket_radius;
    float kdist;
    float sas_weight;
    int   solvent_flat;
    float solventterm;
};

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

struct GpuCfExtraStatic {
    const float* atom_pbvdw;
    const float* atom_charge;
    const int*   atom_hflags;

    int          n_cons;
    const int*   cons_i;
    const int*   cons_j;
    const float* cons_bondlen;
    const float* cons_maxdist;

    int          gist_nx, gist_ny, gist_nz;
    float        gist_origin[3];
    float        gist_inv_delta[3];
    float        gist_weight;
    const float* gist_data;
};

#endif // FLEXAIDS_GPU_CF_TYPES_DEFINED

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

// Upload the optional full-fidelity static arrays (pb-vdw / charge / hbond
// flags / constraints / GIST). Call once after metal_eval_init. Passing
// extra==nullptr (or leaving individual pointers null) disables the matching
// terms; the base com/wal/sas path is unaffected.
void metal_eval_set_extra(MetalEvalCtx* ctx, const GpuCfExtraStatic* extra);

// Evaluate a batch of chromosomes on the Metal GPU (full-fidelity CF).
//   ctx      – context from metal_eval_init
//   pop_size – number of chromosomes to evaluate this call
//   n_genes  – genes per chromosome
//   h_genes  – host gene array [pop_size × n_genes, double]
//   params   – per-batch scalar weights/config (must be non-null)
//   out      – host output channels; each array has length >= pop_size
void metal_eval_batch(MetalEvalCtx*       ctx,
                      int                 pop_size,
                      int                 n_genes,
                      const double*       h_genes,
                      const GpuCfParams*  params,
                      const GpuCfResults* out);

// Free all Metal device resources.
void metal_eval_shutdown(MetalEvalCtx* ctx);

// Stash the per-batch CF params into the context. Used by the multi-complex
// path (whose public signature is fixed by GPUContextPool.h and cannot take a
// params argument) and as a fallback source of params for internal delegations.
// Call once per dock before dispatch; params are a per-dock constant.
void metal_eval_set_params(MetalEvalCtx* ctx, const GpuCfParams* params);

// ─── Multi-complex batched evaluation (GPU utilization maximization) ──────────
//
// Dispatches N independent chromosome populations in a single Metal command
// buffer for many-ligand single-dispatch screening. Thread mapping:
// thread k → complex_id = k / pop_size, chrom_id = k % pop_size.
//
// Fidelity note (P4): the multi path is a screening pre-filter and computes the
// com (distance-weighted) / wal / sas channels only, using the params stashed
// via metal_eval_set_params(). The con / elec / hbond / gist_desolv / pb_clash
// channels are NOT computed here (those need per-complex extra arrays that the
// GPUContextPool queue does not carry) and remain whatever the caller left them.
// Single-complex GA runs (batch_n<=1) go through metal_eval_batch above, which
// is full-fidelity. This signature and MetalMultiBatchEntry are held fixed
// because GPUContextPool.h stores and dispatches them.
//
// Requirements:
//   - All complexes share the same n_genes and pop_size.
//   - The shared context must be init'd with max_pop = n_complex × pop_size.
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

void metal_eval_batch_multi(MetalEvalCtx*               ctx,
                             int                         n_complex,
                             int                         pop_size,
                             int                         n_genes,
                             const MetalMultiBatchEntry* entries);

#endif  // FLEXAIDS_USE_METAL
