// gpu_fast_optics_metal.mm — Metal GPU kNN bridge for FastOPTICS
//
// Objective-C++ bridge that creates Metal compute pipeline, uploads point
// data, launches the gpuFastOPTICSKernel, and downloads results.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#ifdef FLEXAIDS_USE_METAL

#include "gpu_fast_optics_metal.h"
#include "FOPTICS.h"
#import <Metal/Metal.h>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>

void metal_foptics_knn(
    const std::vector<std::pair<chromosome*, std::vector<float>>>& points,
    int k, int nDim,
    std::vector<std::vector<int>>& out_neighbors,
    std::vector<std::vector<float>>& out_distances)
{
    int N = static_cast<int>(points.size());
    if (N == 0 || k <= 0) return;

    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            fprintf(stderr, "[metal_foptics_knn] No Metal device available\n");
            return;
        }

        id<MTLCommandQueue> queue = [device newCommandQueue];
        if (!queue) {
            fprintf(stderr, "[metal_foptics_knn] Failed to create command queue\n");
            return;
        }

        // Locate pre-compiled metallib
        NSBundle* bundle = [NSBundle mainBundle];
        NSString* metallibPath = [bundle pathForResource:@"gpu_fast_optics_metal"
                                                 ofType:@"metallib"];

        if (!metallibPath) {
            NSString* exePath = [bundle executablePath];
            NSString* exeDir = [exePath stringByDeletingLastPathComponent];
            metallibPath = [exeDir stringByAppendingPathComponent:
                            @"gpu_fast_optics_metal.metallib"];
        }

        if (!metallibPath ||
            ![[NSFileManager defaultManager] fileExistsAtPath:metallibPath]) {
            fprintf(stderr, "[metal_foptics_knn] gpu_fast_optics_metal.metallib not found\n");
            return;
        }

        NSError* error = nil;
        NSURL* metallibURL = [NSURL fileURLWithPath:metallibPath];
        id<MTLLibrary> library = [device newLibraryWithURL:metallibURL
                                                     error:&error];
        if (!library) {
            fprintf(stderr, "[metal_foptics_knn] Failed to load metallib: %s\n",
                    error ? [[error localizedDescription] UTF8String] : "unknown");
            return;
        }

        id<MTLFunction> kernelFunc = [library newFunctionWithName:@"gpuFastOPTICSKernel"];
        if (!kernelFunc) {
            fprintf(stderr, "[metal_foptics_knn] Kernel function not found\n");
            return;
        }

        id<MTLComputePipelineState> pipeline =
            [device newComputePipelineStateWithFunction:kernelFunc error:&error];
        if (!pipeline) {
            fprintf(stderr, "[metal_foptics_knn] Failed to create pipeline: %s\n",
                    error ? [[error localizedDescription] UTF8String] : "unknown");
            return;
        }

        // Flatten points into row-major array
        std::vector<float> h_points(N * nDim);
        for (int i = 0; i < N; ++i) {
            const auto& coords = points[i].second;
            for (int d = 0; d < nDim; ++d) {
                h_points[i * nDim + d] =
                    (d < static_cast<int>(coords.size())) ? coords[d] : 0.0f;
            }
        }

        // Create Metal buffers (Shared mode for unified memory on Apple Silicon)
        MTLResourceOptions bufOpts = MTLResourceStorageModeShared;

        id<MTLBuffer> d_points = [device newBufferWithBytes:h_points.data()
                                                     length:static_cast<NSUInteger>(N * nDim * sizeof(float))
                                                    options:bufOpts];

        id<MTLBuffer> d_knn_idx = [device newBufferWithLength:static_cast<NSUInteger>(N * k * sizeof(int))
                                                       options:bufOpts];

        id<MTLBuffer> d_knn_dist = [device newBufferWithLength:static_cast<NSUInteger>(N * k * sizeof(float))
                                                        options:bufOpts];

        // Scalar constant buffers for N, D, k
        id<MTLBuffer> d_N = [device newBufferWithBytes:&N
                                                length:sizeof(int)
                                               options:bufOpts];
        id<MTLBuffer> d_D = [device newBufferWithBytes:&nDim
                                                length:sizeof(int)
                                               options:bufOpts];
        id<MTLBuffer> d_k = [device newBufferWithBytes:&k
                                                length:sizeof(int)
                                               options:bufOpts];

        if (!d_points || !d_knn_idx || !d_knn_dist) {
            fprintf(stderr, "[metal_foptics_knn] Failed to allocate Metal buffers\n");
            return;
        }

        // Create command buffer and encoder
        id<MTLCommandBuffer> cmdBuffer = [queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [cmdBuffer computeCommandEncoder];

        [encoder setComputePipelineState:pipeline];
        [encoder setBuffer:d_points   offset:0 atIndex:0];
        [encoder setBuffer:d_knn_idx  offset:0 atIndex:1];
        [encoder setBuffer:d_knn_dist offset:0 atIndex:2];
        [encoder setBuffer:d_N        offset:0 atIndex:3];
        [encoder setBuffer:d_D        offset:0 atIndex:4];
        [encoder setBuffer:d_k        offset:0 atIndex:5];

        // Dispatch: N threadgroups, 256 threads per group
        MTLSize threadgroupSize = MTLSizeMake(256, 1, 1);
        MTLSize gridSize = MTLSizeMake(static_cast<NSUInteger>(N), 1, 1);
        [encoder dispatchThreadgroups:gridSize
                threadsPerThreadgroup:threadgroupSize];

        [encoder endEncoding];
        [cmdBuffer commit];
        [cmdBuffer waitUntilCompleted];

        if (cmdBuffer.status != MTLCommandBufferStatusCompleted) {
            fprintf(stderr, "[metal_foptics_knn] Command buffer failed\n");
            return;
        }

        // Download results (Shared mode: direct pointer access)
        int*   h_knn_idx  = static_cast<int*>(d_knn_idx.contents);
        float* h_knn_dist = static_cast<float*>(d_knn_dist.contents);

        out_neighbors.resize(N);
        out_distances.resize(N);
        for (int i = 0; i < N; ++i) {
            out_neighbors[i].clear();
            out_distances[i].clear();
            int base = i * k;
            for (int j = 0; j < k; ++j) {
                if (h_knn_idx[base + j] >= 0) {
                    out_neighbors[i].push_back(h_knn_idx[base + j]);
                    out_distances[i].push_back(h_knn_dist[base + j]);
                }
            }
        }
    }
}

#endif // FLEXAIDS_USE_METAL
