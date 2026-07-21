// webgpu_eval.cpp — WebGPU batched chromosome evaluation (Dawn / wgpu-native)
//
// Cross-platform sibling of LIB/metal_eval.mm. Uses the vendor-neutral
// webgpu.h C API implemented by both Dawn and wgpu-native, so this file
// compiles unmodified against either backend (CMake picks whichever was
// detected — see cmake/FlexAIDWebGPU.cmake).
//
// Experimental approximate GPU chromosome evaluation — see the parity note
// in cf_eval.wgsl. Not used for production GA fitness; opt-in via
// --backend webgpu.

#ifdef FLEXAIDS_USE_WEBGPU

#include "webgpu_eval.h"
#include <webgpu/webgpu.h>

#include <cstdio>
#include <cstring>
#include <vector>
#include <atomic>
#include <fstream>
#include <sstream>
#include <string>

namespace {

std::string read_shader_source() {
    // Shipped next to the binary by CMake (see FLEXAIDS_WEBGPU_SHADER_DIR).
    const char* dir = std::getenv("FLEXAIDDS_WEBGPU_SHADER_DIR");
    std::string path = dir ? (std::string(dir) + "/cf_eval.wgsl")
                            : "cf_eval.wgsl";
    std::ifstream f(path);
    if (!f) return {};
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

} // namespace

struct WebGPUEvalCtx {
    WGPUInstance device_instance = nullptr;
    WGPUAdapter  adapter         = nullptr;
    WGPUDevice   device          = nullptr;
    WGPUQueue    queue           = nullptr;
    WGPUComputePipeline pipeline = nullptr;

    WGPUBuffer buf_atom_xyz    = nullptr;
    WGPUBuffer buf_atom_type   = nullptr;
    WGPUBuffer buf_atom_radius = nullptr;
    WGPUBuffer buf_emat        = nullptr;
    WGPUBuffer buf_genes       = nullptr;
    WGPUBuffer buf_com_out     = nullptr;
    WGPUBuffer buf_wal_out     = nullptr;
    WGPUBuffer buf_sas_out     = nullptr;
    WGPUBuffer buf_params      = nullptr;
    WGPUBuffer buf_readback[3] = {nullptr, nullptr, nullptr};

    WGPUBindGroup bind_group = nullptr;

    int n_atoms = 0, n_types = 0, max_pop = 0;
    int lig_first = 0, lig_last = 0;
    float perm = 0.0f;
};

static bool g_probed = false;
static bool g_available = false;

bool webgpu_eval_runtime_available() {
    if (g_probed) return g_available;
    g_probed = true;
    WGPUInstanceDescriptor desc{};
    WGPUInstance inst = wgpuCreateInstance(&desc);
    g_available = (inst != nullptr);
    if (inst) wgpuInstanceRelease(inst);
    return g_available;
}

void webgpu_eval_get_capabilities(WebGPUCapabilities* out) {
    if (!out) return;
    std::memset(out, 0, sizeof(*out));
    out->available = webgpu_eval_runtime_available();
    std::snprintf(out->adapter_name, sizeof(out->adapter_name), "%s",
                  out->available ? "WebGPU adapter (Dawn/wgpu-native)" : "unavailable");
}

static WGPUBuffer make_buffer(WGPUDevice dev, uint64_t size, WGPUBufferUsage usage,
                               const void* init_data) {
    WGPUBufferDescriptor bd{};
    bd.size  = size;
    bd.usage = usage;
    bd.mappedAtCreation = init_data != nullptr;
    WGPUBuffer buf = wgpuDeviceCreateBuffer(dev, &bd);
    if (init_data && buf) {
        void* mapped = wgpuBufferGetMappedRange(buf, 0, size);
        std::memcpy(mapped, init_data, size);
        wgpuBufferUnmap(buf);
    }
    return buf;
}

WebGPUEvalCtx* webgpu_eval_init(int n_atoms, int n_types, int max_pop,
                                 int lig_first, int lig_last, float perm,
                                 const float* h_atom_xyz, const int* h_atom_type,
                                 const float* h_atom_radius,
                                 const float* h_emat_sampled, int n_emat_samples) {
    if (n_emat_samples != WEBGPU_EMAT_SAMPLES) {
        std::fprintf(stderr, "webgpu_eval_init: n_emat_samples mismatch (%d != %d)\n",
                     n_emat_samples, WEBGPU_EMAT_SAMPLES);
        return nullptr;
    }
    if (!webgpu_eval_runtime_available()) return nullptr;

    auto* ctx = new WebGPUEvalCtx();
    ctx->n_atoms = n_atoms; ctx->n_types = n_types; ctx->max_pop = max_pop;
    ctx->lig_first = lig_first; ctx->lig_last = lig_last; ctx->perm = perm;

    WGPUInstanceDescriptor idesc{};
    ctx->device_instance = wgpuCreateInstance(&idesc);

    // Synchronous adapter/device request (Dawn/wgpu-native both support the
    // blocking C-callback pattern under CMAKE-detected FLEXAIDS_USE_WEBGPU).
    struct AdapterResult { WGPUAdapter adapter = nullptr; bool done = false; } ar;
    WGPURequestAdapterOptions opts{};
    wgpuInstanceRequestAdapter(
        ctx->device_instance, &opts,
        [](WGPURequestAdapterStatus status, WGPUAdapter adapter, const char*, void* userdata) {
            auto* r = static_cast<AdapterResult*>(userdata);
            if (status == WGPURequestAdapterStatus_Success) r->adapter = adapter;
            r->done = true;
        },
        &ar);
    if (!ar.adapter) { delete ctx; return nullptr; }
    ctx->adapter = ar.adapter;

    struct DeviceResult { WGPUDevice device = nullptr; bool done = false; } dr;
    WGPUDeviceDescriptor ddesc{};
    wgpuAdapterRequestDevice(
        ctx->adapter, &ddesc,
        [](WGPURequestDeviceStatus status, WGPUDevice device, const char*, void* userdata) {
            auto* r = static_cast<DeviceResult*>(userdata);
            if (status == WGPURequestDeviceStatus_Success) r->device = device;
            r->done = true;
        },
        &dr);
    if (!dr.device) { delete ctx; return nullptr; }
    ctx->device = dr.device;
    ctx->queue  = wgpuDeviceGetQueue(ctx->device);

    std::string src = read_shader_source();
    if (src.empty()) {
        std::fprintf(stderr, "webgpu_eval_init: cf_eval.wgsl not found "
                     "(set FLEXAIDDS_WEBGPU_SHADER_DIR)\n");
        delete ctx; return nullptr;
    }
    WGPUShaderModuleWGSLDescriptor wgsl_desc{};
    wgsl_desc.chain.sType = WGPUSType_ShaderModuleWGSLDescriptor;
    wgsl_desc.code = src.c_str();
    WGPUShaderModuleDescriptor smod_desc{};
    smod_desc.nextInChain = reinterpret_cast<WGPUChainedStruct*>(&wgsl_desc);
    WGPUShaderModule shader = wgpuDeviceCreateShaderModule(ctx->device, &smod_desc);

    WGPUComputePipelineDescriptor pdesc{};
    pdesc.compute.module     = shader;
    pdesc.compute.entryPoint = "kernel_eval_cf_full";
    ctx->pipeline = wgpuDeviceCreateComputePipeline(ctx->device, &pdesc);
    wgpuShaderModuleRelease(shader);

    const uint64_t xyz_bytes    = static_cast<uint64_t>(n_atoms) * 3 * sizeof(float);
    const uint64_t type_bytes   = static_cast<uint64_t>(n_atoms) * sizeof(int32_t);
    const uint64_t radius_bytes = static_cast<uint64_t>(n_atoms) * sizeof(float);
    const uint64_t emat_bytes   = static_cast<uint64_t>(n_types) * n_types * n_emat_samples * sizeof(float);

    ctx->buf_atom_xyz    = make_buffer(ctx->device, xyz_bytes,
                                        WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, h_atom_xyz);
    ctx->buf_atom_type   = make_buffer(ctx->device, type_bytes,
                                        WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, h_atom_type);
    ctx->buf_atom_radius = make_buffer(ctx->device, radius_bytes,
                                        WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, h_atom_radius);
    ctx->buf_emat        = make_buffer(ctx->device, emat_bytes,
                                        WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, h_emat_sampled);

    const uint64_t genes_bytes = static_cast<uint64_t>(max_pop) * 32 * sizeof(float); // upper-bound gene stride
    ctx->buf_genes   = make_buffer(ctx->device, genes_bytes, WGPUBufferUsage_Storage | WGPUBufferUsage_CopyDst, nullptr);
    ctx->buf_com_out = make_buffer(ctx->device, max_pop * sizeof(float), WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, nullptr);
    ctx->buf_wal_out = make_buffer(ctx->device, max_pop * sizeof(float), WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, nullptr);
    ctx->buf_sas_out = make_buffer(ctx->device, max_pop * sizeof(float), WGPUBufferUsage_Storage | WGPUBufferUsage_CopySrc, nullptr);
    ctx->buf_params  = make_buffer(ctx->device, 32, WGPUBufferUsage_Uniform | WGPUBufferUsage_CopyDst, nullptr);

    for (auto& b : ctx->buf_readback)
        b = make_buffer(ctx->device, max_pop * sizeof(float), WGPUBufferUsage_MapRead | WGPUBufferUsage_CopyDst, nullptr);

    WGPUBindGroupLayout layout = wgpuComputePipelineGetBindGroupLayout(ctx->pipeline, 0);
    WGPUBindGroupEntry entries[9] = {};
    WGPUBuffer bufs[9] = { ctx->buf_atom_xyz, ctx->buf_atom_type, ctx->buf_atom_radius,
                           ctx->buf_emat, ctx->buf_genes, ctx->buf_com_out,
                           ctx->buf_wal_out, ctx->buf_sas_out, ctx->buf_params };
    for (int i = 0; i < 9; ++i) {
        entries[i].binding = static_cast<uint32_t>(i);
        entries[i].buffer  = bufs[i];
        entries[i].size    = WGPU_WHOLE_SIZE;
    }
    WGPUBindGroupDescriptor bgdesc{};
    bgdesc.layout     = layout;
    bgdesc.entryCount = 9;
    bgdesc.entries    = entries;
    ctx->bind_group = wgpuDeviceCreateBindGroup(ctx->device, &bgdesc);
    wgpuBindGroupLayoutRelease(layout);

    return ctx;
}

void webgpu_eval_batch(WebGPUEvalCtx* ctx, int pop_size, int n_genes,
                        const double* h_genes, double* h_com_out,
                        double* h_wal_out, double* h_sas_out) {
    if (!ctx) return;

    std::vector<float> genes_f(static_cast<size_t>(pop_size) * n_genes);
    for (size_t i = 0; i < genes_f.size(); ++i) genes_f[i] = static_cast<float>(h_genes[i]);
    wgpuQueueWriteBuffer(ctx->queue, ctx->buf_genes, 0, genes_f.data(), genes_f.size() * sizeof(float));

    struct { int32_t N, T, n_genes, lig_first, lig_last; float perm; int32_t pad0, pad1; } params{
        ctx->n_atoms, ctx->n_types, n_genes, ctx->lig_first, ctx->lig_last, ctx->perm, 0, 0
    };
    wgpuQueueWriteBuffer(ctx->queue, ctx->buf_params, 0, &params, sizeof(params));

    WGPUCommandEncoderDescriptor encdesc{};
    WGPUCommandEncoder enc = wgpuDeviceCreateCommandEncoder(ctx->device, &encdesc);

    WGPUComputePassDescriptor passdesc{};
    WGPUComputePassEncoder pass = wgpuCommandEncoderBeginComputePass(enc, &passdesc);
    wgpuComputePassEncoderSetPipeline(pass, ctx->pipeline);
    wgpuComputePassEncoderSetBindGroup(pass, 0, ctx->bind_group, 0, nullptr);
    wgpuComputePassEncoderDispatchWorkgroups(pass, static_cast<uint32_t>(pop_size), 1, 1);
    wgpuComputePassEncoderEnd(pass);
    wgpuComputePassEncoderRelease(pass);

    const uint64_t out_bytes = static_cast<uint64_t>(pop_size) * sizeof(float);
    WGPUBuffer src[3] = { ctx->buf_com_out, ctx->buf_wal_out, ctx->buf_sas_out };
    for (int i = 0; i < 3; ++i)
        wgpuCommandEncoderCopyBufferToBuffer(enc, src[i], 0, ctx->buf_readback[i], 0, out_bytes);

    WGPUCommandBufferDescriptor cbdesc{};
    WGPUCommandBuffer cmd = wgpuCommandEncoderFinish(enc, &cbdesc);
    wgpuCommandEncoderRelease(enc);
    wgpuQueueSubmit(ctx->queue, 1, &cmd);
    wgpuCommandBufferRelease(cmd);

    double* outs[3] = { h_com_out, h_wal_out, h_sas_out };
    for (int i = 0; i < 3; ++i) {
        struct MapResult { std::atomic<bool> done{false}; } mr;
        wgpuBufferMapAsync(
            ctx->buf_readback[i], WGPUMapMode_Read, 0, out_bytes,
            [](WGPUBufferMapAsyncStatus, void* userdata) {
                static_cast<MapResult*>(userdata)->done.store(true);
            },
            &mr);
        while (!mr.done.load()) { /* Dawn/wgpu-native drive the device via their own poll loop */ }
        const float* mapped = static_cast<const float*>(
            wgpuBufferGetConstMappedRange(ctx->buf_readback[i], 0, out_bytes));
        for (int p = 0; p < pop_size; ++p) outs[i][p] = static_cast<double>(mapped[p]);
        wgpuBufferUnmap(ctx->buf_readback[i]);
    }
}

void webgpu_eval_shutdown(WebGPUEvalCtx* ctx) {
    if (!ctx) return;
    for (auto b : { ctx->buf_atom_xyz, ctx->buf_atom_type, ctx->buf_atom_radius, ctx->buf_emat,
                     ctx->buf_genes, ctx->buf_com_out, ctx->buf_wal_out, ctx->buf_sas_out, ctx->buf_params,
                     ctx->buf_readback[0], ctx->buf_readback[1], ctx->buf_readback[2] })
        if (b) wgpuBufferRelease(b);
    if (ctx->bind_group) wgpuBindGroupRelease(ctx->bind_group);
    if (ctx->pipeline)   wgpuComputePipelineRelease(ctx->pipeline);
    if (ctx->queue)      wgpuQueueRelease(ctx->queue);
    if (ctx->device)     wgpuDeviceRelease(ctx->device);
    if (ctx->adapter)    wgpuAdapterRelease(ctx->adapter);
    if (ctx->device_instance) wgpuInstanceRelease(ctx->device_instance);
    delete ctx;
}

#endif // FLEXAIDS_USE_WEBGPU
