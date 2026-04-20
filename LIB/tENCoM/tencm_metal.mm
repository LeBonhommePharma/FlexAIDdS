// tencm_metal.mm — Objective-C++ Metal bridge for TENCoM GPU kernels
//
// Implements the C++ interface declared in tencm_metal.h using Apple Metal
// compute pipelines for contact discovery and Hessian assembly.
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma / NRGlab

#include "tencm_metal.h"

#ifdef FLEXAIDS_HAS_METAL_TENCM

#import <Metal/Metal.h>
#include <cstdio>
#include <cstring>

namespace tencm { namespace metal {

static id<MTLDevice>        s_device   = nil;
static id<MTLCommandQueue>  s_queue    = nil;
static id<MTLComputePipelineState> s_contact_pipeline  = nil;
static id<MTLComputePipelineState> s_hessian_pipeline  = nil;
static bool s_initialized = false;
static bool s_available   = false;

bool init() {
    if (s_initialized) return s_available;

    s_initialized = true;
    s_available   = false;

    s_device = MTLCreateSystemDefaultDevice();
    if (!s_device) {
        std::fprintf(stderr, "[tencm_metal] No Metal device found\n");
        return false;
    }

    s_queue = [s_device newCommandQueue];
    if (!s_queue) {
        std::fprintf(stderr, "[tencm_metal] Failed to create command queue\n");
        return false;
    }

    NSError* error = nil;
    NSBundle* bundle = [NSBundle mainBundle];
    NSString* metallibPath = [bundle pathForResource:@"tencm_kernels"
                                              ofType:@"metallib"];

    if (metallibPath) {
        NSURL* metallibURL = [NSURL fileURLWithPath:metallibPath];
        id<MTLLibrary> library = [s_device newLibraryWithURL:metallibURL
                                                        error:&error];
        if (library) {
            MTLComputePipelineDescriptor* cpd =
                [[MTLComputePipelineDescriptor alloc] init];

            cpd.computeFunction = [library newFunctionWithName:@"tencm_contacts"];
            s_contact_pipeline = [s_device newComputePipelineStateWithDescriptor:cpd
                                                                          options:0
                                                                       reflection:nil
                                                                            error:&error];

            cpd.computeFunction = [library newFunctionWithName:@"tencm_hessian"];
            s_hessian_pipeline = [s_device newComputePipelineStateWithDescriptor:cpd
                                                                          options:0
                                                                       reflection:nil
                                                                            error:&error];
        } else {
            std::fprintf(stderr, "[tencm_metal] metallib load failed: %s\n",
                         error.localizedDescription.UTF8String);
        }
    } else {
        std::fprintf(stderr, "[tencm_metal] tencm_kernels.metallib not found — "
                     "GPU kernels unavailable, using CPU fallback\n");
    }

    s_available = true;
    std::fprintf(stderr, "[tencm_metal] initialized on %s\n",
                 s_device.name.UTF8String);
    return true;
}

void shutdown() {
    s_hessian_pipeline = nil;
    s_contact_pipeline = nil;
    s_queue  = nil;
    s_device = nil;
    s_initialized = false;
    s_available   = false;
}

bool is_available() {
    return s_available;
}

int build_contacts_gpu(const float* ca_xyz, int N,
                       float cutoff, float k0,
                       std::vector<int>&   contacts_ij,
                       std::vector<float>& contacts_k,
                       std::vector<float>& contacts_r0) {
    if (!s_available || N < GPU_THRESHOLD || !s_contact_pipeline) return 0;

    @autoreleasepool {
        int n3 = N * 3;
        size_t xyz_bytes = n3 * sizeof(float);

        id<MTLBuffer> xyz_buf = [s_device newBufferWithBytes:ca_xyz
                                                      length:xyz_bytes
                                                     options:MTLResourceStorageModeShared];

        int max_contacts = N * (N - 1) / 2;
        size_t ij_bytes = max_contacts * 2 * sizeof(int);
        size_t ck_bytes = max_contacts * sizeof(float);

        id<MTLBuffer> ij_buf = [s_device newBufferWithLength:ij_bytes
                                                      options:MTLResourceStorageModeShared];
        id<MTLBuffer> ck_buf = [s_device newBufferWithLength:ck_bytes
                                                      options:MTLResourceStorageModeShared];
        id<MTLBuffer> r0_buf = [s_device newBufferWithLength:ck_bytes
                                                      options:MTLResourceStorageModeShared];

        id<MTLBuffer> count_buf = [s_device newBufferWithLength:sizeof(int)
                                                         options:MTLResourceStorageModeShared];
        memset(count_buf.contents, 0, sizeof(int));

        struct Params { int N; float cutoff; float k0; };
        Params p{N, cutoff, k0};
        id<MTLBuffer> params_buf = [s_device newBufferWithBytes:&p
                                                          length:sizeof(p)
                                                         options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        [enc setComputePipelineState:s_contact_pipeline];
        [enc setBuffer:xyz_buf     offset:0 atIndex:0];
        [enc setBuffer:ij_buf      offset:0 atIndex:1];
        [enc setBuffer:ck_buf      offset:0 atIndex:2];
        [enc setBuffer:r0_buf      offset:0 atIndex:3];
        [enc setBuffer:count_buf   offset:0 atIndex:4];
        [enc setBuffer:params_buf  offset:0 atIndex:5];

        MTLSize grid   = MTLSizeMake(N, N, 1);
        MTLSize thread = MTLSizeMake(16, 16, 1);
        [enc dispatchThreads:grid threadsPerThreadgroup:thread];
        [enc endEncoding];

        [cmd commit];
        [cmd waitUntilCompleted];

        int count = *static_cast<const int*>(count_buf.contents);
        if (count > 0) {
            const int*   ij_data = static_cast<const int*>(ij_buf.contents);
            const float* ck_data = static_cast<const float*>(ck_buf.contents);
            const float* r0_data = static_cast<const float*>(r0_buf.contents);

            contacts_ij.assign(ij_data, ij_data + count * 2);
            contacts_k.assign(ck_data, ck_data + count);
            contacts_r0.assign(r0_data, r0_data + count);
        }
        return count;
    }
}

void assemble_hessian_gpu(const float* ca_xyz, int N,
                          const int* contacts_ij,
                          const float* contacts_k,
                          int M, int C,
                          double* H_out) {
    if (!s_available || !s_hessian_pipeline) return;

    memset(H_out, 0, M * M * sizeof(double));

    @autoreleasepool {
        size_t n3 = N * 3;
        id<MTLBuffer> xyz_buf = [s_device newBufferWithBytes:ca_xyz
                                                      length:n3 * sizeof(float)
                                                     options:MTLResourceStorageModeShared];
        id<MTLBuffer> ij_buf = [s_device newBufferWithBytes:contacts_ij
                                                      length:C * 2 * sizeof(int)
                                                     options:MTLResourceStorageModeShared];
        id<MTLBuffer> ck_buf = [s_device newBufferWithBytes:contacts_k
                                                      length:C * sizeof(float)
                                                     options:MTLResourceStorageModeShared];
        id<MTLBuffer> H_buf = [s_device newBufferWithLength:M * M * sizeof(double)
                                                     options:MTLResourceStorageModeShared];

        struct Params { int N; int M; int C; };
        Params p{N, M, C};
        id<MTLBuffer> params_buf = [s_device newBufferWithBytes:&p
                                                          length:sizeof(p)
                                                         options:MTLResourceStorageModeShared];

        id<MTLCommandBuffer> cmd = [s_queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        [enc setComputePipelineState:s_hessian_pipeline];
        [enc setBuffer:xyz_buf     offset:0 atIndex:0];
        [enc setBuffer:ij_buf      offset:0 atIndex:1];
        [enc setBuffer:ck_buf      offset:0 atIndex:2];
        [enc setBuffer:H_buf       offset:0 atIndex:3];
        [enc setBuffer:params_buf  offset:0 atIndex:4];

        MTLSize grid   = MTLSizeMake(M, M, 1);
        MTLSize thread = MTLSizeMake(16, 16, 1);
        [enc dispatchThreads:grid threadsPerThreadgroup:thread];
        [enc endEncoding];

        [cmd commit];
        [cmd waitUntilCompleted];

        memcpy(H_out, H_buf.contents, M * M * sizeof(double));
    }
}

}}  // namespace tencm::metal

#endif  // FLEXAIDS_HAS_METAL_TENCM
