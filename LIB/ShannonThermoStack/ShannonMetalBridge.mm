// ShannonMetalBridge.mm — Objective-C++ bridge to Metal GPU kernels
//
// Compiled only on APPLE targets with FLEXAIDS_HAS_METAL_SHANNON defined.
// Persistent device/pipeline/queue caching eliminates per-call init overhead.
// Dispatches: histogram, Boltzmann weights, parallel sum, log-sum-exp.
//
// NOTE: Metal Shading Language does not support double. All GPU buffers use
// float (FP32). This bridge converts double↔float at the host boundary.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include "ShannonMetalBridge.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <mutex>
#include <string>

namespace ShannonMetalBridge {

// ─── Persistent Metal context (singleton, thread-safe init) ─────────────────

struct MetalContext {
    id<MTLDevice>               device;
    id<MTLCommandQueue>         queue;
    id<MTLLibrary>              library;
    // Cached pipelines
    id<MTLComputePipelineState> histogramPipeline;
    id<MTLComputePipelineState> boltzmannPipeline;
    id<MTLComputePipelineState> sumReducePipeline;
    id<MTLComputePipelineState> logSumExpPipeline;
    // P2: reusable device buffers (grow-only pool)
    id<MTLBuffer>               pooled_energy_buf;
    id<MTLBuffer>               pooled_weight_buf;
    id<MTLBuffer>               pooled_bin_buf;
    id<MTLBuffer>               pooled_exp_buf;
    id<MTLBuffer>               pooled_partial_buf;
    NSUInteger                  pooled_energy_cap  = 0;
    NSUInteger                  pooled_weight_cap  = 0;
    NSUInteger                  pooled_bin_cap     = 0;
    NSUInteger                  pooled_exp_cap     = 0;
    NSUInteger                  pooled_partial_cap = 0;
    bool                        valid;
    std::string                 deviceInfo;
};

static id<MTLBuffer> ensure_buffer(MetalContext& ctx,
                                   id<MTLBuffer>* slot,
                                   NSUInteger* cap,
                                   NSUInteger need_bytes)
{
    if (*cap < need_bytes) {
        *slot = [ctx.device newBufferWithLength:need_bytes
                                         options:MTLResourceStorageModeShared];
        *cap = need_bytes;
    }
    return *slot;
}

static MetalContext& get_context() {
    static MetalContext ctx{};
    static std::once_flag flag;

    std::call_once(flag, [&] {
        ctx.valid = false;

        ctx.device = MTLCreateSystemDefaultDevice();
        if (!ctx.device) return;

        ctx.queue = [ctx.device newCommandQueue];
        if (!ctx.queue) return;

        // Load library from default bundle or compiled metallib
        NSError* err = nil;
        ctx.library = [ctx.device newDefaultLibrary];
        if (!ctx.library) return;

        // Build pipelines for each kernel
        auto makePipeline = [&](NSString* name) -> id<MTLComputePipelineState> {
            id<MTLFunction> fn = [ctx.library newFunctionWithName:name];
            if (!fn) return nil;
            return [ctx.device newComputePipelineStateWithFunction:fn error:&err];
        };

        ctx.histogramPipeline  = makePipeline(@"shannon_histogram");
        ctx.boltzmannPipeline  = makePipeline(@"boltzmann_weights_batch");
        ctx.sumReducePipeline  = makePipeline(@"parallel_sum_reduce");
        ctx.logSumExpPipeline  = makePipeline(@"log_sum_exp_shifted");

        ctx.valid = (ctx.histogramPipeline != nil);

        ctx.deviceInfo = std::string([[ctx.device name] UTF8String]);
    });

    return ctx;
}

// ─── CPU fallback Shannon entropy from bin counts ───────────────────────────

static double cpu_shannon_from_bins(const std::vector<int>& bins) {
    int total = 0;
    for (int c : bins) total += c;
    if (total == 0) return 0.0;

    double H = 0.0;
    for (int c : bins) {
        if (c > 0) {
            double p = static_cast<double>(c) / total;
            H -= p * std::log(p);
        }
    }
    return H;
}

// ─── CPU fallback histogram ─────────────────────────────────────────────────

static double cpu_shannon_fallback(const std::vector<double>& energies, int num_bins) {
    double min_v = *std::min_element(energies.begin(), energies.end());
    double max_v = *std::max_element(energies.begin(), energies.end());
    if (max_v - min_v < 1e-12) return 0.0;
    double bw = (max_v - min_v) / num_bins + 1e-10;
    std::vector<int> bins(num_bins, 0);
    for (double e : energies) {
        int b = std::min(std::max((int)((e - min_v) / bw), 0), num_bins - 1);
        bins[b]++;
    }
    return cpu_shannon_from_bins(bins);
}

// ─── double → float conversion helper ───────────────────────────────────────

static std::vector<float> to_float(const std::vector<double>& src) {
    std::vector<float> dst(src.size());
    for (size_t i = 0; i < src.size(); ++i)
        dst[i] = static_cast<float>(src[i]);
    return dst;
}

// ─── Shannon entropy (GPU) ──────────────────────────────────────────────────

double compute_shannon_entropy_metal(const std::vector<double>& energies,
                                     int num_bins)
{
    if (energies.empty()) return 0.0;
    if (num_bins <= 0) num_bins = 20;

    auto& ctx = get_context();
    if (!ctx.valid || !ctx.histogramPipeline) {
        return cpu_shannon_fallback(energies, num_bins);
    }

    NSUInteger n = energies.size();
    float min_v = static_cast<float>(*std::min_element(energies.begin(), energies.end()));
    float max_v = static_cast<float>(*std::max_element(energies.begin(), energies.end()));
    if (max_v - min_v < 1e-12f) return 0.0;
    float bw = (max_v - min_v) / num_bins + 1e-10f;

    // Convert double → float for GPU
    std::vector<float> energies_f = to_float(energies);

    const NSUInteger energy_bytes = n * sizeof(float);
    const NSUInteger bin_bytes = static_cast<NSUInteger>(num_bins) * sizeof(int);
    id<MTLBuffer> energy_buf = ensure_buffer(ctx, &ctx.pooled_energy_buf,
                                               &ctx.pooled_energy_cap, energy_bytes);
    memcpy(energy_buf.contents, energies_f.data(), energy_bytes);

    id<MTLBuffer> bin_buf = ensure_buffer(ctx, &ctx.pooled_bin_buf,
                                          &ctx.pooled_bin_cap, bin_bytes);
    memset(bin_buf.contents, 0, bin_bytes);

    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    [enc setComputePipelineState:ctx.histogramPipeline];
    [enc setBuffer:energy_buf offset:0 atIndex:0];
    [enc setBuffer:bin_buf    offset:0 atIndex:1];
    [enc setBytes:&n          length:sizeof(NSUInteger) atIndex:2];
    [enc setBytes:&num_bins   length:sizeof(int)        atIndex:3];
    [enc setBytes:&min_v      length:sizeof(float)      atIndex:4];
    [enc setBytes:&bw         length:sizeof(float)      atIndex:5];

    MTLSize tpg = MTLSizeMake(256, 1, 1);
    MTLSize ng  = MTLSizeMake((n + 255) / 256, 1, 1);
    [enc dispatchThreadgroups:ng threadsPerThreadgroup:tpg];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    int* bin_data = static_cast<int*>(bin_buf.contents);
    std::vector<int> bins(bin_data, bin_data + num_bins);
    return cpu_shannon_from_bins(bins);
}

// ─── Boltzmann weights (GPU) ────────────────────────────────────────────────

std::vector<double> compute_boltzmann_weights_metal(
    const std::vector<double>& energies,
    double beta,
    double& sum_w,
    double& E_min)
{
    const NSUInteger n = energies.size();
    sum_w = 0.0;
    E_min = 0.0;

    if (energies.empty()) return {};

    auto& ctx = get_context();

    // Pre-compute E_min on CPU (single pass)
    E_min = *std::min_element(energies.begin(), energies.end());
    float neg_beta_f = static_cast<float>(-beta);
    float E_min_f = static_cast<float>(E_min);

    if (!ctx.valid || !ctx.boltzmannPipeline) {
        // CPU fallback
        std::vector<double> weights(n);
        for (NSUInteger i = 0; i < n; ++i) {
            weights[i] = std::exp(-beta * (energies[i] - E_min));
            sum_w += weights[i];
        }
        return weights;
    }

    // Convert double → float for GPU
    std::vector<float> energies_f = to_float(energies);

    // GPU path
    const NSUInteger n_bytes = n * sizeof(float);
    id<MTLBuffer> energy_buf = ensure_buffer(ctx, &ctx.pooled_energy_buf,
                                               &ctx.pooled_energy_cap, n_bytes);
    memcpy(energy_buf.contents, energies_f.data(), n_bytes);

    id<MTLBuffer> weight_buf = ensure_buffer(ctx, &ctx.pooled_weight_buf,
                                             &ctx.pooled_weight_cap, n_bytes);

    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    [enc setComputePipelineState:ctx.boltzmannPipeline];
    [enc setBuffer:energy_buf  offset:0 atIndex:0];
    [enc setBuffer:weight_buf  offset:0 atIndex:1];
    uint32_t n32 = static_cast<uint32_t>(n);
    [enc setBytes:&n32       length:sizeof(uint32_t) atIndex:2];
    [enc setBytes:&neg_beta_f length:sizeof(float)    atIndex:3];
    [enc setBytes:&E_min_f   length:sizeof(float)     atIndex:4];

    MTLSize tpg = MTLSizeMake(256, 1, 1);
    MTLSize ng  = MTLSizeMake((n + 255) / 256, 1, 1);
    [enc dispatchThreadgroups:ng threadsPerThreadgroup:tpg];
    [enc endEncoding];

    // If sum reduction pipeline is available, use GPU for sum too
    if (ctx.sumReducePipeline && n > 1024) {
        NSUInteger numGroups = (n + 255) / 256;
        const NSUInteger partial_bytes = numGroups * sizeof(float);
        id<MTLBuffer> partial_buf = ensure_buffer(ctx, &ctx.pooled_partial_buf,
                                                    &ctx.pooled_partial_cap, partial_bytes);

        id<MTLComputeCommandEncoder> enc2 = [cmd computeCommandEncoder];
        [enc2 setComputePipelineState:ctx.sumReducePipeline];
        [enc2 setBuffer:weight_buf   offset:0 atIndex:0];
        [enc2 setBuffer:partial_buf  offset:0 atIndex:1];
        [enc2 setBytes:&n32          length:sizeof(uint32_t) atIndex:2];
        [enc2 dispatchThreadgroups:MTLSizeMake(numGroups, 1, 1) threadsPerThreadgroup:tpg];
        [enc2 endEncoding];

        [cmd commit];
        [cmd waitUntilCompleted];

        // Final sum on CPU from partial sums (small array, float→double)
        float* partials = static_cast<float*>(partial_buf.contents);
        sum_w = 0.0;
        for (NSUInteger i = 0; i < numGroups; ++i)
            sum_w += static_cast<double>(partials[i]);
    } else {
        [cmd commit];
        [cmd waitUntilCompleted];

        // Sum on CPU (float→double)
        float* w = static_cast<float*>(weight_buf.contents);
        sum_w = 0.0;
        for (NSUInteger i = 0; i < n; ++i)
            sum_w += static_cast<double>(w[i]);
    }

    // Copy results (float→double)
    float* w = static_cast<float*>(weight_buf.contents);
    std::vector<double> weights(n);
    for (NSUInteger i = 0; i < n; ++i)
        weights[i] = static_cast<double>(w[i]);
    return weights;
}

// ─── Log-sum-exp (GPU) ─────────────────────────────────────────────────────

double log_sum_exp_metal(const std::vector<double>& values) {
    if (values.empty())
        return -std::numeric_limits<double>::infinity();

    const NSUInteger n = values.size();
    double x_max = *std::max_element(values.begin(), values.end());
    if (!std::isfinite(x_max)) return x_max;

    auto& ctx = get_context();
    if (!ctx.valid || !ctx.logSumExpPipeline) {
        // CPU fallback
        double sum = 0.0;
        for (double v : values)
            sum += std::exp(v - x_max);
        return x_max + std::log(sum);
    }

    // Convert double → float for GPU
    std::vector<float> values_f = to_float(values);
    float x_max_f = static_cast<float>(x_max);

    // GPU: compute exp(x - x_max) then sum
    const NSUInteger n_bytes = n * sizeof(float);
    id<MTLBuffer> val_buf = ensure_buffer(ctx, &ctx.pooled_energy_buf,
                                          &ctx.pooled_energy_cap, n_bytes);
    memcpy(val_buf.contents, values_f.data(), n_bytes);

    id<MTLBuffer> exp_buf = ensure_buffer(ctx, &ctx.pooled_exp_buf,
                                          &ctx.pooled_exp_cap, n_bytes);

    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    uint32_t n32 = static_cast<uint32_t>(n);
    [enc setComputePipelineState:ctx.logSumExpPipeline];
    [enc setBuffer:val_buf offset:0 atIndex:0];
    [enc setBuffer:exp_buf offset:0 atIndex:1];
    [enc setBytes:&n32     length:sizeof(uint32_t) atIndex:2];
    [enc setBytes:&x_max_f length:sizeof(float)    atIndex:3];

    MTLSize tpg = MTLSizeMake(256, 1, 1);
    MTLSize ng  = MTLSizeMake((n + 255) / 256, 1, 1);
    [enc dispatchThreadgroups:ng threadsPerThreadgroup:tpg];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // CPU sum of exp-shifted values (float→double)
    float* exp_data = static_cast<float*>(exp_buf.contents);
    double sum = 0.0;
    for (NSUInteger i = 0; i < n; ++i)
        sum += static_cast<double>(exp_data[i]);

    return x_max + std::log(sum);
}

// ─── Utility ────────────────────────────────────────────────────────────────

bool is_metal_available() {
    return get_context().valid;
}

std::string metal_device_info() {
    auto& ctx = get_context();
    if (!ctx.valid) return "Metal unavailable";
    return ctx.deviceInfo;
}

} // namespace ShannonMetalBridge
