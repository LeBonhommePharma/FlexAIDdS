// webgpu_eval.h — WebGPU GPU evaluation API for batched chromosome scoring
//
// Cross-platform sibling of LIB/metal_eval.h, targeting Dawn or wgpu-native
// via the common webgpu.h header. Enabled with -DFLEXAIDS_USE_WEBGPU.
// Same experimental/reduced scoring model as the Metal path (see the parity
// note in cf_eval.wgsl) — selected only via --backend webgpu.
//
// Usage:
//   WebGPUEvalCtx* ctx = webgpu_eval_init(...);   // once, before GA loop
//   webgpu_eval_batch(ctx, ...);                   // every generation
//   webgpu_eval_shutdown(ctx);                     // once, after GA loop
#pragma once

#ifdef FLEXAIDS_USE_WEBGPU

#include <cstddef>

static constexpr int WEBGPU_EMAT_SAMPLES = 128;

// Opaque handle to all WebGPU device-resident state.
struct WebGPUEvalCtx;

// Runtime availability probe. Compile-time WebGPU support does not guarantee
// that the current process can actually acquire an adapter/device (Dawn or
// wgpu-native must be able to find a Vulkan/Metal/D3D12 backend at runtime).
bool webgpu_eval_runtime_available();

struct WebGPUCapabilities {
    bool available;
    char adapter_name[128];
    char backend_type[32];   // "Vulkan" / "Metal" / "D3D12" / "OpenGL"
};

// Fills the struct. Safe to call even if WebGPU is not compiled in.
void webgpu_eval_get_capabilities(WebGPUCapabilities* out);

// Allocate WebGPU device buffers and compile the compute pipeline.
// Same parameter contract as metal_eval_init() (see LIB/metal_eval.h).
WebGPUEvalCtx* webgpu_eval_init(int   n_atoms,
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

// Evaluate a batch of chromosomes on the WebGPU device.
//   ctx        – context from webgpu_eval_init
//   pop_size   – number of chromosomes to evaluate this call
//   n_genes    – genes per chromosome
//   h_genes    – host gene array [pop_size × n_genes, double]
//   h_com_out  – host output: complementarity CF   [pop_size, double]
//   h_wal_out  – host output: wall/clash energy     [pop_size, double]
//   h_sas_out  – host output: solvent-accessible    [pop_size, double]
void webgpu_eval_batch(WebGPUEvalCtx* ctx,
                        int            pop_size,
                        int            n_genes,
                        const double*  h_genes,
                        double*        h_com_out,
                        double*        h_wal_out,
                        double*        h_sas_out);

// Free all WebGPU device resources.
void webgpu_eval_shutdown(WebGPUEvalCtx* ctx);

#endif // FLEXAIDS_USE_WEBGPU
