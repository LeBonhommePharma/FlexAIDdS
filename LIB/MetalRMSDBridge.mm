// MetalRMSDBridge.mm — Objective-C++ bridge to Metal GPU pairwise RMSD kernel
//
// Persistent device/pipeline/queue caching eliminates per-call init overhead.
// Handles N=0, N=1 (trivial), N<=4096 (single dispatch), and N>4096 (tiled).
//
// Apache-2.0 (C) 2026 Le Bonhomme Pharma

#include "MetalRMSDBridge.h"

#ifdef FLEXAIDS_USE_METAL

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <mutex>
#include <string>

namespace metal_rmsd {

// ─── Persistent Metal context (singleton, thread-safe) ──────────────────────

struct MetalRMSDContext {
    id<MTLDevice>               device;
    id<MTLCommandQueue>         queue;
    id<MTLComputePipelineState> pairwisePipeline;
    bool                        valid;
    std::string                 deviceInfo;
};

static MetalRMSDContext& get_context() {
    static MetalRMSDContext ctx{};
    static std::once_flag flag;

    std::call_once(flag, [&] {
        ctx.valid = false;

        ctx.device = MTLCreateSystemDefaultDevice();
        if (!ctx.device) return;

        ctx.queue = [ctx.device newCommandQueue];
        if (!ctx.queue) return;

        // Load the MetalRMSD.metallib (path injected by CMake at build time)
        NSError* err = nil;
        NSString* metallib_path = @METALRMSD_METALLIB_PATH;
        NSURL* metallib_url = [NSURL fileURLWithPath:metallib_path];
        id<MTLLibrary> lib = [ctx.device newLibraryWithURL:metallib_url error:&err];

        if (!lib) {
            // Fallback: try default library (shader compiled into app bundle)
            lib = [ctx.device newDefaultLibrary];
        }
        if (!lib) return;

        id<MTLFunction> fn = [lib newFunctionWithName:@"pairwise_rmsd_squared"];
        if (!fn) return;

        ctx.pairwisePipeline = [ctx.device newComputePipelineStateWithFunction:fn error:&err];
        if (!ctx.pairwisePipeline) return;

        ctx.valid = true;
        ctx.deviceInfo = std::string([[ctx.device name] UTF8String]);
    });

    return ctx;
}

// ─── CPU fallback pairwise RMSD ─────────────────────────────────────────────

static void cpu_pairwise_rmsd(const float* coords, int n_conf, int n_atoms,
                               std::vector<float>& dist_matrix)
{
    dist_matrix.resize(static_cast<size_t>(n_conf) * n_conf, 0.0f);
    int stride = 3 * n_atoms;

    for (int i = 0; i < n_conf; ++i) {
        dist_matrix[static_cast<size_t>(i) * n_conf + i] = 0.0f; // diagonal
        for (int j = i + 1; j < n_conf; ++j) {
            float sum_sq = 0.0f;
            const float* ci = coords + static_cast<ptrdiff_t>(i) * stride;
            const float* cj = coords + static_cast<ptrdiff_t>(j) * stride;
            for (int k = 0; k < stride; ++k) {
                float diff = ci[k] - cj[k];
                sum_sq += diff * diff;
            }
            float rmsd = sqrtf(sum_sq / static_cast<float>(n_atoms));
            dist_matrix[static_cast<size_t>(i) * n_conf + j] = rmsd;
            dist_matrix[static_cast<size_t>(j) * n_conf + i] = rmsd;
        }
    }
}

// ─── GPU pairwise RMSD (single dispatch, N <= 4096) ─────────────────────────

static bool gpu_pairwise_rmsd_single(const float* coords, int n_conf, int n_atoms,
                                      std::vector<float>& dist_matrix)
{
    auto& ctx = get_context();
    if (!ctx.valid || !ctx.pairwisePipeline) return false;

    NSUInteger N = static_cast<NSUInteger>(n_conf);
    NSUInteger M3 = static_cast<NSUInteger>(n_atoms) * 3;
    uint64_t total_pairs = static_cast<uint64_t>(N) * (N - 1) / 2;

    // Allocate output for upper-triangular raw sums
    NSUInteger coord_bytes = N * M3 * sizeof(float);
    NSUInteger output_bytes = static_cast<NSUInteger>(total_pairs) * sizeof(float);

    id<MTLBuffer> coordBuf = [ctx.device newBufferWithBytes:coords
                                                     length:coord_bytes
                                                    options:MTLResourceStorageModeShared];
    id<MTLBuffer> outBuf = [ctx.device newBufferWithLength:output_bytes
                                                   options:MTLResourceStorageModeShared];

    if (!coordBuf || !outBuf) return false;

    // Dispatch
    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    [enc setComputePipelineState:ctx.pairwisePipeline];
    [enc setBuffer:coordBuf offset:0 atIndex:0];
    [enc setBuffer:outBuf   offset:0 atIndex:1];
    [enc setBytes:&N       length:sizeof(NSUInteger) atIndex:2];
    [enc setBytes:&M3      length:sizeof(NSUInteger) atIndex:3];

    NSUInteger tpg = ctx.pairwisePipeline.maxTotalThreadsPerThreadgroup;
    MTLSize grid     = MTLSizeMake(static_cast<NSUInteger>(total_pairs), 1, 1);
    MTLSize threadgrp = MTLSizeMake(tpg, 1, 1);

    [enc dispatchThreads:grid threadsPerThreadgroup:threadgrp];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    if (cmd.status == MTLCommandBufferStatusError) {
        NSString* desc = cmd.error ? cmd.error.localizedDescription : @"unknown";
        fprintf(stderr, "[MetalRMSD] GPU command buffer error: %s\n",
                [desc UTF8String]);
        return false;
    }

    // Read back upper-triangular results and expand to full N x N matrix
    const float* raw = static_cast<const float*>(outBuf.contents);
    float inv_n_atoms = 1.0f / static_cast<float>(n_atoms);

    dist_matrix.assign(static_cast<size_t>(N) * N, 0.0f);

    uint64_t idx = 0;
    for (NSUInteger i = 0; i < N; ++i) {
        for (NSUInteger j = i + 1; j < N; ++j) {
            float rmsd_sq = raw[idx] * inv_n_atoms; // divide by n_atoms for RMSD^2
            float rmsd = sqrtf(rmsd_sq);
            dist_matrix[static_cast<size_t>(i) * N + j] = rmsd;
            dist_matrix[static_cast<size_t>(j) * N + i] = rmsd;
            ++idx;
        }
    }

    return true;
}

// ─── GPU pairwise RMSD (tiled, N > 4096) ────────────────────────────────────
//
// Tiles over row ranges to keep GPU memory bounded.  Each tile computes
// pairs for rows [row_start, row_end) against all j > i.  Results are
// written directly into the correct positions of the full N x N matrix.

static bool gpu_pairwise_rmsd_tiled(const float* coords, int n_conf, int n_atoms,
                                     std::vector<float>& dist_matrix)
{
    auto& ctx = get_context();
    if (!ctx.valid) return false;

    NSUInteger N = static_cast<NSUInteger>(n_conf);
    NSUInteger M3 = static_cast<NSUInteger>(n_atoms) * 3;
    NSUInteger coord_bytes = N * M3 * sizeof(float);

    // Upload coords once (shared across all tiles)
    id<MTLBuffer> coordBuf = [ctx.device newBufferWithBytes:coords
                                                     length:coord_bytes
                                                    options:MTLResourceStorageModeShared];
    if (!coordBuf) return false;

    dist_matrix.assign(static_cast<size_t>(N) * N, 0.0f);

    // Tile parameters: each tile processes tile_size rows
    // Memory per tile: ~tile_size * N * 4 bytes (output) + coord_bytes
    // For N=16384, tile_size=512 → ~32 MB output + ~150 MB coords = ~182 MB
    const NSUInteger tile_size = 512;
    float inv_n_atoms = 1.0f / static_cast<float>(n_atoms);

    // We need the tiled kernel pipeline
    NSError* err = nil;
    NSString* metallib_path = @METALRMSD_METALLIB_PATH;
    NSURL* metallib_url = [NSURL fileURLWithPath:metallib_path];
    id<MTLLibrary> lib = [ctx.device newLibraryWithURL:metallib_url error:&err];
    if (!lib) {
        lib = [ctx.device newDefaultLibrary];
    }
    if (!lib) return false;

    id<MTLFunction> tiledFn = [lib newFunctionWithName:@"pairwise_rmsd_squared_tiled"];
    if (!tiledFn) {
        // No tiled kernel — fall back to CPU for safety
        return false;
    }

    id<MTLComputePipelineState> tiledPipeline =
        [ctx.device newComputePipelineStateWithFunction:tiledFn error:&err];
    if (!tiledPipeline) return false;

    for (NSUInteger row_start = 0; row_start < N; row_start += tile_size) {
        NSUInteger row_end = std::min(row_start + tile_size, N);

        // Count pairs in this tile
        uint64_t tile_pairs = 0;
        for (NSUInteger i = row_start; i < row_end; ++i) {
            tile_pairs += (N - 1 - i);
        }

        NSUInteger output_bytes = static_cast<NSUInteger>(tile_pairs) * sizeof(float);
        id<MTLBuffer> outBuf = [ctx.device newBufferWithLength:output_bytes
                                                       options:MTLResourceStorageModeShared];
        if (!outBuf) return false;

        id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        [enc setComputePipelineState:tiledPipeline];
        [enc setBuffer:coordBuf offset:0 atIndex:0];
        [enc setBuffer:outBuf   offset:0 atIndex:1];
        [enc setBytes:&N         length:sizeof(NSUInteger) atIndex:2];
        [enc setBytes:&M3        length:sizeof(NSUInteger) atIndex:3];

        uint32_t rs = static_cast<uint32_t>(row_start);
        uint32_t re = static_cast<uint32_t>(row_end);
        uint64_t zero_offset = 0;
        [enc setBytes:&rs          length:sizeof(uint32_t) atIndex:4];
        [enc setBytes:&re          length:sizeof(uint32_t) atIndex:5];
        [enc setBytes:&zero_offset length:sizeof(uint64_t)    atIndex:6];

        NSUInteger tpg = tiledPipeline.maxTotalThreadsPerThreadgroup;
        MTLSize grid     = MTLSizeMake(static_cast<NSUInteger>(tile_pairs), 1, 1);
        MTLSize threadgrp = MTLSizeMake(tpg, 1, 1);

        [enc dispatchThreads:grid threadsPerThreadgroup:threadgrp];
        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];

        if (cmd.status == MTLCommandBufferStatusError) {
            return false;
        }

        // Expand tile results into full matrix
        const float* raw = static_cast<const float*>(outBuf.contents);
        uint64_t idx = 0;
        for (NSUInteger i = row_start; i < row_end; ++i) {
            for (NSUInteger j = i + 1; j < N; ++j) {
                float rmsd = sqrtf(raw[idx] * inv_n_atoms);
                dist_matrix[static_cast<size_t>(i) * N + j] = rmsd;
                dist_matrix[static_cast<size_t>(j) * N + i] = rmsd;
                ++idx;
            }
        }
    }

    return true;
}

// ─── Public API ──────────────────────────────────────────────────────────────

bool compute_pairwise_rmsd_metal(const float* coords, int n_conf, int n_atoms,
                                  std::vector<float>& dist_matrix)
{
    // Edge cases: N=0 or N=1
    if (n_conf <= 0) {
        dist_matrix.clear();
        return true; // trivially handled
    }
    if (n_conf == 1) {
        dist_matrix.assign(1, 0.0f);
        return true;
    }

    // Check Metal availability
    auto& ctx = get_context();
    if (!ctx.valid) {
        cpu_pairwise_rmsd(coords, n_conf, n_atoms, dist_matrix);
        return false;
    }

    // Choose dispatch strategy based on N
    bool gpu_ok;
    if (n_conf <= 4096) {
        gpu_ok = gpu_pairwise_rmsd_single(coords, n_conf, n_atoms, dist_matrix);
    } else {
        gpu_ok = gpu_pairwise_rmsd_tiled(coords, n_conf, n_atoms, dist_matrix);
    }

    if (!gpu_ok) {
        cpu_pairwise_rmsd(coords, n_conf, n_atoms, dist_matrix);
        return false;
    }

    return true;
}

bool is_metal_rmsd_available() {
    return get_context().valid;
}

const char* metal_rmsd_device_info() {
    auto& ctx = get_context();
    if (!ctx.valid) return "Metal RMSD unavailable";
    return ctx.deviceInfo.c_str();
}

} // namespace metal_rmsd

#else // FLEXAIDS_USE_METAL not defined

#include "MetalRMSDBridge.h"
#include <cmath>

namespace metal_rmsd {

static void cpu_pairwise_rmsd(const float* coords, int n_conf, int n_atoms,
                               std::vector<float>& dist_matrix)
{
    dist_matrix.resize(static_cast<size_t>(n_conf) * n_conf, 0.0f);
    int stride = 3 * n_atoms;

    for (int i = 0; i < n_conf; ++i) {
        dist_matrix[static_cast<size_t>(i) * n_conf + i] = 0.0f;
        for (int j = i + 1; j < n_conf; ++j) {
            float sum_sq = 0.0f;
            const float* ci = coords + static_cast<ptrdiff_t>(i) * stride;
            const float* cj = coords + static_cast<ptrdiff_t>(j) * stride;
            for (int k = 0; k < stride; ++k) {
                float diff = ci[k] - cj[k];
                sum_sq += diff * diff;
            }
            float rmsd = sqrtf(sum_sq / static_cast<float>(n_atoms));
            dist_matrix[static_cast<size_t>(i) * n_conf + j] = rmsd;
            dist_matrix[static_cast<size_t>(j) * n_conf + i] = rmsd;
        }
    }
}

bool compute_pairwise_rmsd_metal(const float* coords, int n_conf, int n_atoms,
                                  std::vector<float>& dist_matrix)
{
    if (n_conf <= 0) {
        dist_matrix.clear();
        return false;
    }
    if (n_conf == 1) {
        dist_matrix.assign(1, 0.0f);
        return false;
    }
    cpu_pairwise_rmsd(coords, n_conf, n_atoms, dist_matrix);
    return false;
}

bool is_metal_rmsd_available() { return false; }

const char* metal_rmsd_device_info() { return "Metal unavailable (not macOS or not enabled)"; }

} // namespace metal_rmsd

#endif // FLEXAIDS_USE_METAL
