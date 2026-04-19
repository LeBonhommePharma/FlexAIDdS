// ShannonMetalBridge.mm — Objective-C++ bridge to Metal GPU kernels
//
// Compiled only on APPLE targets with FLEXAIDS_HAS_METAL_SHANNON defined.
// Persistent device/pipeline/queue caching eliminates per-call init overhead.
// Dispatches: histogram, Boltzmann weights, parallel sum, log-sum-exp.
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
    bool                        valid;
    std::string                 deviceInfo;
};

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
    double min_v = *std::min_element(energies.begin(), energies.end());
    double max_v = *std::max_element(energies.begin(), energies.end());
    if (max_v - min_v < 1e-12) return 0.0;
    double bw = (max_v - min_v) / num_bins + 1e-10;

    // Buffers (shared memory — zero-copy on Apple Silicon)
    id<MTLBuffer> energy_buf = [ctx.device newBufferWithBytes:energies.data()
                                                       length:n * sizeof(double)
                                                      options:MTLResourceStorageModeShared];
    id<MTLBuffer> bin_buf = [ctx.device newBufferWithLength:num_bins * sizeof(int)
                                                    options:MTLResourceStorageModeShared];
    memset(bin_buf.contents, 0, num_bins * sizeof(int));

    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    [enc setComputePipelineState:ctx.histogramPipeline];
    [enc setBuffer:energy_buf offset:0 atIndex:0];
    [enc setBuffer:bin_buf    offset:0 atIndex:1];
    [enc setBytes:&n          length:sizeof(NSUInteger) atIndex:2];
    [enc setBytes:&num_bins   length:sizeof(int)        atIndex:3];
    [enc setBytes:&min_v      length:sizeof(double)     atIndex:4];
    [enc setBytes:&bw         length:sizeof(double)     atIndex:5];

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
    double neg_beta = -beta;

    if (!ctx.valid || !ctx.boltzmannPipeline) {
        // CPU fallback
        std::vector<double> weights(n);
        for (NSUInteger i = 0; i < n; ++i) {
            weights[i] = std::exp(-beta * (energies[i] - E_min));
            sum_w += weights[i];
        }
        return weights;
    }

    // GPU path
    id<MTLBuffer> energy_buf = [ctx.device newBufferWithBytes:energies.data()
                                                       length:n * sizeof(double)
                                                      options:MTLResourceStorageModeShared];
    id<MTLBuffer> weight_buf = [ctx.device newBufferWithLength:n * sizeof(double)
                                                       options:MTLResourceStorageModeShared];

    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    [enc setComputePipelineState:ctx.boltzmannPipeline];
    [enc setBuffer:energy_buf  offset:0 atIndex:0];
    [enc setBuffer:weight_buf  offset:0 atIndex:1];
    uint32_t n32 = static_cast<uint32_t>(n);
    [enc setBytes:&n32       length:sizeof(uint32_t) atIndex:2];
    [enc setBytes:&neg_beta  length:sizeof(double)   atIndex:3];
    [enc setBytes:&E_min     length:sizeof(double)   atIndex:4];

    MTLSize tpg = MTLSizeMake(256, 1, 1);
    MTLSize ng  = MTLSizeMake((n + 255) / 256, 1, 1);
    [enc dispatchThreadgroups:ng threadsPerThreadgroup:tpg];
    [enc endEncoding];

    // If sum reduction pipeline is available, use GPU for sum too
    if (ctx.sumReducePipeline && n > 1024) {
        NSUInteger numGroups = (n + 255) / 256;
        id<MTLBuffer> partial_buf = [ctx.device newBufferWithLength:numGroups * sizeof(double)
                                                            options:MTLResourceStorageModeShared];

        id<MTLComputeCommandEncoder> enc2 = [cmd computeCommandEncoder];
        [enc2 setComputePipelineState:ctx.sumReducePipeline];
        [enc2 setBuffer:weight_buf   offset:0 atIndex:0];
        [enc2 setBuffer:partial_buf  offset:0 atIndex:1];
        [enc2 setBytes:&n32          length:sizeof(uint32_t) atIndex:2];
        [enc2 dispatchThreadgroups:MTLSizeMake(numGroups, 1, 1) threadsPerThreadgroup:tpg];
        [enc2 endEncoding];

        [cmd commit];
        [cmd waitUntilCompleted];

        // Final sum on CPU from partial sums (small array)
        double* partials = static_cast<double*>(partial_buf.contents);
        sum_w = 0.0;
        for (NSUInteger i = 0; i < numGroups; ++i)
            sum_w += partials[i];
    } else {
        [cmd commit];
        [cmd waitUntilCompleted];

        // Sum on CPU
        double* w = static_cast<double*>(weight_buf.contents);
        sum_w = 0.0;
        for (NSUInteger i = 0; i < n; ++i)
            sum_w += w[i];
    }

    // Copy results
    double* w = static_cast<double*>(weight_buf.contents);
    return std::vector<double>(w, w + n);
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

    // GPU: compute exp(x - x_max) then sum
    id<MTLBuffer> val_buf = [ctx.device newBufferWithBytes:values.data()
                                                    length:n * sizeof(double)
                                                   options:MTLResourceStorageModeShared];
    id<MTLBuffer> exp_buf = [ctx.device newBufferWithLength:n * sizeof(double)
                                                    options:MTLResourceStorageModeShared];

    id<MTLCommandBuffer> cmd = [ctx.queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    uint32_t n32 = static_cast<uint32_t>(n);
    [enc setComputePipelineState:ctx.logSumExpPipeline];
    [enc setBuffer:val_buf offset:0 atIndex:0];
    [enc setBuffer:exp_buf offset:0 atIndex:1];
    [enc setBytes:&n32     length:sizeof(uint32_t) atIndex:2];
    [enc setBytes:&x_max   length:sizeof(double)   atIndex:3];

    MTLSize tpg = MTLSizeMake(256, 1, 1);
    MTLSize ng  = MTLSizeMake((n + 255) / 256, 1, 1);
    [enc dispatchThreadgroups:ng threadsPerThreadgroup:tpg];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // CPU sum of exp-shifted values
    double* exp_data = static_cast<double*>(exp_buf.contents);
    double sum = 0.0;
    for (NSUInteger i = 0; i < n; ++i)
        sum += exp_data[i];

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
