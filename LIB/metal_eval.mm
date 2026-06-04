// metal_eval.mm — Metal GPU batched chromosome evaluation
//
// Implements the same full-fidelity CF scoring as cuda_eval.cu but using
// Apple's Metal compute API.  The MSL kernel is compiled at runtime from
// an embedded string; no separate .metal compilation step is needed.
//
// Scoring pipeline (per chromosome, one threadgroup per chromosome):
//   1. Decode translation genes (tx, ty, tz) from the gene vector.
//   2. For each ligand-protein atom pair:
//        a. Compute inter-atomic distance r.
//        b. Approximate contact area (linear switching 0→1 as r→rA+rB).
//        c. Look up energy value via linear interpolation in the
//           pre-sampled density-function table.
//        d. Accumulate COM contribution and WAL (clash) energy.
//        e. Subtract contact area from per-ligand-atom SAS counter
//           using a CAS-loop float atomic in threadgroup memory.
//   3. Compute SAS energy contribution for each ligand atom using the
//      remaining exposed surface and the solvent column of the energy matrix.
//   4. Reduce COM, WAL, SAS across the threadgroup and write outputs.

#ifdef FLEXAIDS_USE_METAL

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>

#include "metal_eval.h"
#include "flexaid_exception.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>
#include <string>

// ─── MSL kernel source (embedded) ────────────────────────────────────────────
static const char* kMSLSource = R"MSL(
#include <metal_stdlib>
#include <metal_atomic>
using namespace metal;

#define N_EMAT_SAMPLES 128

// GPU-side linear interpolation into the pre-sampled energy-matrix table.
static float gpu_get_yval(device const float* emat_sampled,
                           int t1, int t2, int T, float rel_area)
{
    int base = (t1 * T + t2) * N_EMAT_SAMPLES;
    rel_area = clamp(rel_area, 0.0f, 1.0f);
    float kf  = rel_area * float(N_EMAT_SAMPLES - 1);
    int   k0  = int(kf);
    int   k1  = min(k0 + 1, N_EMAT_SAMPLES - 1);
    float frac = kf - float(k0);
    return emat_sampled[base + k0] * (1.0f - frac)
         + emat_sampled[base + k1] * frac;
}

// CAS-loop float atomic subtract in threadgroup memory (Metal 2-compatible).
static void tg_atomic_sub_float(threadgroup float* ptr, float val)
{
    threadgroup atomic_uint* ap = (threadgroup atomic_uint*)ptr;
    uint old_bits, new_bits;
    do {
        old_bits = atomic_load_explicit(ap, memory_order_relaxed);
        float new_val = as_type<float>(old_bits) - val;
        new_bits = as_type<uint>(new_val);
    } while (!atomic_compare_exchange_weak_explicit(
                ap, &old_bits, new_bits,
                memory_order_relaxed, memory_order_relaxed));
}

// Params packed into one buffer for convenience.
struct EvalParams {
    int   N;           // total atom count
    int   T;           // atom type count
    int   n_genes;
    int   lig_first;
    int   lig_last;
    float perm;
    int   pad0;
    int   pad1;
};

kernel void kernel_eval_cf_full(
    device const float*    atom_xyz        [[ buffer(0) ]],
    device const int*      atom_type       [[ buffer(1) ]],
    device const float*    atom_radius     [[ buffer(2) ]],
    device const float*    emat_sampled    [[ buffer(3) ]],
    device const float*    genes_f         [[ buffer(4) ]],  // float cast of double genes
    device float*          cf_com_out      [[ buffer(5) ]],
    device float*          cf_wal_out      [[ buffer(6) ]],
    device float*          cf_sas_out      [[ buffer(7) ]],
    constant EvalParams&   p               [[ buffer(8) ]],
    threadgroup float*     lig_sas         [[ threadgroup(0) ]],
    uint tid                               [[ thread_position_in_threadgroup ]],
    uint chrom_id                          [[ threadgroup_position_in_grid ]],
    uint blockDim                          [[ threads_per_threadgroup ]])
{
    const int n_lig   = p.lig_last - p.lig_first + 1;
    const int n_pro   = p.N - n_lig;
    const int n_pairs = n_lig * n_pro;

    // Load translation from genes (first 3 genes).
    const int gbase = int(chrom_id) * p.n_genes;
    const float tx = genes_f[gbase + 0];
    const float ty = genes_f[gbase + 1];
    const float tz = genes_f[gbase + 2];

    // Initialise per-ligand SAS to full surface area.
    for (int la = int(tid); la < n_lig && la < 256; la += int(blockDim)) {
        float ra  = atom_radius[p.lig_first + la];
        float rwa = ra + 1.4f;
        lig_sas[la] = 4.0f * M_PI_F * rwa * rwa;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float local_com = 0.0f, local_wal = 0.0f;

    for (int pr = int(tid); pr < n_pairs; pr += int(blockDim)) {
        const int li      = pr / n_pro;
        const int pro_rel = pr % n_pro;
        const int ai      = p.lig_first + li;
        const int aj      = (pro_rel < p.lig_first) ? pro_rel : (pro_rel + n_lig);

        const float lx = atom_xyz[ai * 3 + 0] + tx;
        const float ly = atom_xyz[ai * 3 + 1] + ty;
        const float lz = atom_xyz[ai * 3 + 2] + tz;
        const float dx = lx - atom_xyz[aj * 3 + 0];
        const float dy = ly - atom_xyz[aj * 3 + 1];
        const float dz = lz - atom_xyz[aj * 3 + 2];
        const float r  = sqrt(dx*dx + dy*dy + dz*dz + 1e-10f);

        const float rA    = atom_radius[ai];
        const float rB    = atom_radius[aj];
        const float rsum  = rA + rB;
        const float rwa_A = rA + 1.4f;
        const float surf_A = 4.0f * M_PI_F * rwa_A * rwa_A;
        const float outer_r = rsum + 2.8f;  // rA + rB + 2*Rw

        // Normalised contact area (0..1), linear switching.
        float rel_area = 0.0f;
        if      (r < rsum)    rel_area = 1.0f;
        else if (r < outer_r) rel_area = 1.0f - (r - rsum) / (outer_r - rsum);

        // Subtract from ligand-atom SAS using CAS float atomic.
        if (rel_area > 0.0f && li < 256) {
            tg_atomic_sub_float(&lig_sas[li], rel_area * surf_A);
        }

        // Complementarity energy (sampled energy matrix lookup).
        const int ti  = atom_type[ai];
        const int tj  = atom_type[aj];
        const float yval = gpu_get_yval(emat_sampled, ti, tj, p.T, rel_area);
        local_com += yval * rel_area;

        // WAL: repulsive wall energy when r < perm * (rA+rB).
        const float clash_r = p.perm * rsum;
        if (r < clash_r && r > 0.0f) {
            const float inv_r12  = 1.0f / pow(r,       12.0f);
            const float inv_cr12 = 1.0f / pow(clash_r, 12.0f);
            local_wal += 1.0e6f * (inv_r12 - inv_cr12);
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // SIMD reduction for COM and WAL.
    local_com = simd_sum(local_com);
    local_wal = simd_sum(local_wal);

    // SAS contribution: each ligand atom's remaining exposed area.
    float local_sas = 0.0f;
    for (int la = int(tid); la < n_lig && la < 256; la += int(blockDim)) {
        const float sas_rem  = max(0.0f, lig_sas[la]);
        const float rwa_la   = atom_radius[p.lig_first + la] + 1.4f;
        const float surf_la  = 4.0f * M_PI_F * rwa_la * rwa_la;
        const float sas_norm = sas_rem / surf_la;
        const int   ti_la    = atom_type[p.lig_first + la];
        const float yval_sas = gpu_get_yval(emat_sampled, ti_la, p.T - 1, p.T, sas_norm);
        local_sas += yval_sas * sas_norm;
    }
    local_sas = simd_sum(local_sas);

    if (tid == 0) {
        cf_com_out[chrom_id] = local_com;
        cf_wal_out[chrom_id] = local_wal;
        cf_sas_out[chrom_id] = local_sas;
    }
}
)MSL";

// ─── context structure ────────────────────────────────────────────────────────
struct MetalEvalCtx {
    id<MTLDevice>              device;
    id<MTLCommandQueue>        queue;
    id<MTLComputePipelineState> pipeline;

    id<MTLBuffer> buf_atom_xyz;
    id<MTLBuffer> buf_atom_type;
    id<MTLBuffer> buf_atom_radius;
    id<MTLBuffer> buf_emat_sampled;
    id<MTLBuffer> buf_genes_f;
    id<MTLBuffer> buf_com_out;
    id<MTLBuffer> buf_wal_out;
    id<MTLBuffer> buf_sas_out;

    int n_atoms;
    int n_types;
    int max_pop;
    int max_genes;
    int lig_first;
    int lig_last;
    float perm;
};

// ─── host API ────────────────────────────────────────────────────────────────

bool metal_eval_runtime_available()
{
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    return device != nil;
}

void metal_eval_get_capabilities(MetalCapabilities* out)
{
    if (!out) return;

    memset(out, 0, sizeof(*out));

    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
        fprintf(stderr, "[DEBUG] MTLCreateSystemDefaultDevice() returned NULL\n");
        out->available = false;
        return;
    }

    fprintf(stderr, "[DEBUG] Metal device acquired: %s\n", device.name.UTF8String ? device.name.UTF8String : "unknown");
    out->available = true;
    strncpy(out->device_name, device.name.UTF8String ? device.name.UTF8String : "Apple GPU", sizeof(out->device_name)-1);

#if defined(__MAC_OS_X_VERSION_MIN_REQUIRED) && __MAC_OS_X_VERSION_MIN_REQUIRED >= 101500
    if (@available(macOS 10.15, *)) {
        out->unified_memory_bytes = device.hasUnifiedMemory ? device.maxBufferLength : 0;
    }
#endif
    out->max_buffer_length = device.maxBufferLength;

    // Rough core estimate for M-series (best effort)
    NSString* name = device.name;
    if ([name containsString:@"M3 Pro"]) out->gpu_core_estimate = 18;
    else if ([name containsString:@"M3 Max"]) out->gpu_core_estimate = 30;
    else if ([name containsString:@"M3"]) out->gpu_core_estimate = 10;
    else out->gpu_core_estimate = 8; // conservative default
}

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
                               int   n_emat_samples)
{
    MetalEvalCtx* ctx = new MetalEvalCtx();
    ctx->n_atoms   = n_atoms;
    ctx->n_types   = n_types;
    ctx->max_pop   = max_pop;
    ctx->lig_first = lig_first;
    ctx->lig_last  = lig_last;
    ctx->perm      = perm;

    // Device & queue.
    ctx->device = MTLCreateSystemDefaultDevice();
    if (!ctx->device) {
        fprintf(stderr, "metal_eval: no Metal device found\n");
        delete ctx;
        return nullptr;
    }
    ctx->queue = [ctx->device newCommandQueue];

    // Compile kernel.
    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kMSLSource];
    id<MTLLibrary> lib = [ctx->device newLibraryWithSource:src
                                                   options:nil
                                                     error:&err];
    if (!lib) {
        fprintf(stderr, "metal_eval: shader compile error: %s\n",
                [[err localizedDescription] UTF8String]);
        delete ctx;
        return nullptr;
    }
    id<MTLFunction> fn = [lib newFunctionWithName:@"kernel_eval_cf_full"];
    ctx->pipeline = [ctx->device newComputePipelineStateWithFunction:fn error:&err];
    if (!ctx->pipeline) {
        fprintf(stderr, "metal_eval: pipeline error: %s\n",
                [[err localizedDescription] UTF8String]);
        delete ctx;
        return nullptr;
    }

    // Allocate constant device buffers (uploaded once).
    auto newBuf = [&](const void* data, size_t bytes) -> id<MTLBuffer> {
        return [ctx->device newBufferWithBytes:data
                                       length:bytes
                                      options:MTLResourceStorageModeShared];
    };

    ctx->buf_atom_xyz    = newBuf(h_atom_xyz,    (size_t)n_atoms * 3 * sizeof(float));
    ctx->buf_atom_type   = newBuf(h_atom_type,   (size_t)n_atoms     * sizeof(int));
    ctx->buf_atom_radius = newBuf(h_atom_radius, (size_t)n_atoms     * sizeof(float));
    ctx->buf_emat_sampled= newBuf(h_emat_sampled,
                                  (size_t)n_types * n_types * n_emat_samples * sizeof(float));

    // Mutable per-batch buffers.
    // Use max 256 genes as upper bound; actual n_genes validated in batch call.
    ctx->max_genes = 256;
    ctx->buf_genes_f = [ctx->device newBufferWithLength:(size_t)max_pop * ctx->max_genes * sizeof(float)
                                                options:MTLResourceStorageModeShared];
    ctx->buf_com_out = [ctx->device newBufferWithLength:(size_t)max_pop * sizeof(float)
                                                options:MTLResourceStorageModeShared];
    ctx->buf_wal_out = [ctx->device newBufferWithLength:(size_t)max_pop * sizeof(float)
                                                options:MTLResourceStorageModeShared];
    ctx->buf_sas_out = [ctx->device newBufferWithLength:(size_t)max_pop * sizeof(float)
                                                options:MTLResourceStorageModeShared];

    return ctx;
}

void metal_eval_batch(MetalEvalCtx* ctx,
                      int           pop_size,
                      int           n_genes,
                      const double* h_genes,
                      double*       h_com_out,
                      double*       h_wal_out,
                      double*       h_sas_out)
{
    // Validate against allocated buffer sizes — throw on overflow.
    if (pop_size > ctx->max_pop) {
        throw FlexAIDException("metal_eval: pop_size " + std::to_string(pop_size) +
            " exceeds max_pop " + std::to_string(ctx->max_pop));
    }
    if (n_genes > ctx->max_genes) {
        throw FlexAIDException("metal_eval: n_genes " + std::to_string(n_genes) +
            " exceeds max_genes " + std::to_string(ctx->max_genes));
    }

    // Convert double genes to float for the GPU.
    // Pack with n_genes stride to match kernel indexing (genes_f[chrom_id * n_genes + g]).
    float* genes_f = (float*)[ctx->buf_genes_f contents];
    for (int c = 0; c < pop_size; ++c)
        for (int g = 0; g < n_genes; ++g)
            genes_f[c * n_genes + g] = (float)h_genes[c * n_genes + g];

    // Build EvalParams.
    struct EvalParams {
        int N, T, n_genes, lig_first, lig_last;
        float perm;
        int pad0, pad1;
    };
    EvalParams ep = { ctx->n_atoms, ctx->n_types, n_genes,
                      ctx->lig_first, ctx->lig_last, ctx->perm, 0, 0 };

    id<MTLBuffer> buf_params = [ctx->device
        newBufferWithBytes:&ep
                   length:sizeof(ep)
                  options:MTLResourceStorageModeShared];

    // Encode and dispatch.
    id<MTLCommandBuffer>       cb  = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:ctx->pipeline];
    [enc setBuffer:ctx->buf_atom_xyz    offset:0 atIndex:0];
    [enc setBuffer:ctx->buf_atom_type   offset:0 atIndex:1];
    [enc setBuffer:ctx->buf_atom_radius offset:0 atIndex:2];
    [enc setBuffer:ctx->buf_emat_sampled offset:0 atIndex:3];
    [enc setBuffer:ctx->buf_genes_f     offset:0 atIndex:4];
    [enc setBuffer:ctx->buf_com_out     offset:0 atIndex:5];
    [enc setBuffer:ctx->buf_wal_out     offset:0 atIndex:6];
    [enc setBuffer:ctx->buf_sas_out     offset:0 atIndex:7];
    [enc setBuffer:buf_params           offset:0 atIndex:8];

    // Threadgroup shared memory: 256 floats for per-ligand SAS.
    [enc setThreadgroupMemoryLength:256 * sizeof(float) atIndex:0];

    NSUInteger threadsPerGroup = 256;
    MTLSize    gridSize        = { (NSUInteger)pop_size, 1, 1 };
    MTLSize    groupSize       = { threadsPerGroup, 1, 1 };
    [enc dispatchThreadgroups:gridSize threadsPerThreadgroup:groupSize];
    [enc endEncoding];

    [cb commit];
    [cb waitUntilCompleted];

    if (cb.status == MTLCommandBufferStatusError) {
        NSString* desc = cb.error ? cb.error.localizedDescription : @"unknown";
        throw FlexAIDException("metal_eval: GPU command buffer error: " +
            std::string([desc UTF8String]));
    }

    // Copy results back (float → double).
    const float* com_f = (const float*)[ctx->buf_com_out contents];
    const float* wal_f = (const float*)[ctx->buf_wal_out contents];
    const float* sas_f = (const float*)[ctx->buf_sas_out contents];
    for (int c = 0; c < pop_size; ++c) {
        h_com_out[c] = (double)com_f[c];
        h_wal_out[c] = (double)wal_f[c];
        h_sas_out[c] = (double)sas_f[c];
    }
}

void metal_eval_shutdown(MetalEvalCtx* ctx)
{
    if (!ctx) return;
    // ARC-managed objects are released automatically.
    delete ctx;
}

#endif  // FLEXAIDS_USE_METAL
