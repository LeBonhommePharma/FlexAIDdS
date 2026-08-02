// metal_eval.mm — Full-fidelity Metal GPU batched chromosome evaluation
//
// P4: brings the single-complex Metal CF up to parity with the CPU
// ic2cf/vcfunction assembly within the ranking-preserving drift tolerance. The
// single-complex MSL kernel now fills the eight scoring channels that feed
// get_cf_evalue():
//
//     com, wal, sas, con, elec, hbond, gist_desolv, pb_clash
//
// Contact-area model: the same C0 linear-switching approximation documented in
// cuda_eval.cuh (NOT the branchy analytic Voronoi). Terms that require true
// absolute geometry (angular H-bond directionality, metal-CN, vct-entropy,
// tENCoM) are not reproduced on-device — the gaboom.cpp divergence guard warns
// when their weights are non-zero.
//
// The multi-complex screening path keeps its GPUContextPool-fixed public API and
// remains a com(dist-weighted)/wal/sas pre-filter (see metal_eval.h).
//
// Reduction correctness: the previous experimental kernels reduced with
// simd_sum() only, which sums within a single 32-lane SIMD group while the
// threadgroup runs 256 threads — silently dropping 7/8 of the per-pair
// contributions. Both kernels below use an explicit power-of-two threadgroup
// tree reduction (tg_reduce_sum / tg_reduce_min) that includes every thread.

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

// ─── host/MSL shared kernel-parameter POD (single-complex) ───────────────────
// Layout MUST match the MSL `KP` struct below (all 4-byte scalars, tight).
struct MetalKP {
    int   N, T, n_genes, lig_first, lig_last;
    float perm;
    int   has_pbvdw, has_charge, has_hflags;
    int   n_cons;
    int   gist_nx, gist_ny, gist_nz;
    float gox, goy, goz, gidx, gidy, gidz, gist_weight;
    float dw_r0;
    int   elec_on;
    float dielectric;
    float hbond_weight, hbond_salt_weight, hbond_opt_dist, hbond_sigma_dist, hbond_angle_repr;
    float pb_clash_weight, pb_clash_ratio, pb_clash_exponent, pb_pocket_weight, pb_pocket_radius;
    float kdist, sas_weight;
    int   solvent_flat;
    float solventterm;
};

// Multi-complex shared params (screening pre-filter: com/wal/sas only).
struct MetalMultiParams {
    int   pop_size, n_genes, T;
    float perm;
    float dw_r0;
    float sas_weight;
    int   solvent_flat;
    float solventterm;
};

static constexpr int METAL_CF_CHANNELS = 8;

// ─── MSL kernel source (embedded) ────────────────────────────────────────────
static const char* kMSLSource = R"MSL(
#include <metal_stdlib>
#include <metal_atomic>
using namespace metal;

#define N_EMAT_SAMPLES 128
#define TG_THREADS 256
#define KWALL_F     1.0e6f
#define KCOULOMB_F  332.0637f
#define HBOND_PAIR_MIN -2.0f
#define Q_SALT 0.30f

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

// GetValueFromGaussian, base clamped >=0 (drift-safe).
static float gpu_gaussian(float x, float mx, float zero)
{
    float dz = zero - mx;
    float denom = dz * dz;
    if (denom < 1e-12f) return 0.0f;
    float base = (-(x - zero) * (x - (2.0f * mx - zero))) / denom;
    if (base <= 0.0f) return 0.0f;
    return pow(base, 50.0f);
}

// Atomic float subtract on a threadgroup counter (Metal-2 compatible CAS loop).
static void tg_atomic_sub_float(threadgroup float* ptr, float val)
{
    threadgroup atomic_uint* ap = (threadgroup atomic_uint*)ptr;
    uint old_bits = atomic_load_explicit(ap, memory_order_relaxed);
    uint new_bits;
    do {
        float nv = as_type<float>(old_bits) - val;
        new_bits = as_type<uint>(nv);
    } while (!atomic_compare_exchange_weak_explicit(
                ap, &old_bits, new_bits,
                memory_order_relaxed, memory_order_relaxed));
}

// Threadgroup tree reductions (TG_THREADS must be a power of two).
static float tg_reduce_sum(threadgroup float* s, uint tid, float v)
{
    s[tid] = v;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint stride = TG_THREADS / 2; stride > 0; stride >>= 1) {
        if (tid < stride) s[tid] += s[tid + stride];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    float r = s[0];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return r;
}
static float tg_reduce_min(threadgroup float* s, uint tid, float v)
{
    s[tid] = v;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint stride = TG_THREADS / 2; stride > 0; stride >>= 1) {
        if (tid < stride) s[tid] = min(s[tid], s[tid + stride]);
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    float r = s[0];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return r;
}

// Single-complex kernel parameters (matches host MetalKP).
struct KP {
    int   N, T, n_genes, lig_first, lig_last;
    float perm;
    int   has_pbvdw, has_charge, has_hflags;
    int   n_cons;
    int   gist_nx, gist_ny, gist_nz;
    float gox, goy, goz, gidx, gidy, gidz, gist_weight;
    float dw_r0;
    int   elec_on;
    float dielectric;
    float hbond_weight, hbond_salt_weight, hbond_opt_dist, hbond_sigma_dist, hbond_angle_repr;
    float pb_clash_weight, pb_clash_ratio, pb_clash_exponent, pb_pocket_weight, pb_pocket_radius;
    float kdist, sas_weight;
    int   solvent_flat;
    float solventterm;
};

// ─── Single-complex full-fidelity kernel ──────────────────────────────────────
kernel void kernel_eval_cf_full(
    device const float*    atom_xyz        [[ buffer(0) ]],
    device const int*      atom_type       [[ buffer(1) ]],
    device const float*    atom_radius     [[ buffer(2) ]],
    device const float*    emat_sampled    [[ buffer(3) ]],
    device const float*    genes_f         [[ buffer(4) ]],
    device float*          out             [[ buffer(5) ]],   // [CH * max_pop]
    constant KP&           p               [[ buffer(6) ]],
    device const float*    atom_pbvdw      [[ buffer(7) ]],
    device const float*    atom_charge     [[ buffer(8) ]],
    device const int*      atom_hflags     [[ buffer(9) ]],
    device const int*      cons_i          [[ buffer(10) ]],
    device const int*      cons_j          [[ buffer(11) ]],
    device const float*    cons_bondlen    [[ buffer(12) ]],
    device const float*    cons_maxdist    [[ buffer(13) ]],
    device const float*    gist_data       [[ buffer(14) ]],
    constant int&          max_pop         [[ buffer(15) ]],
    threadgroup float*     lig_sas         [[ threadgroup(0) ]],
    threadgroup float*     red             [[ threadgroup(1) ]],
    uint tid                               [[ thread_position_in_threadgroup ]],
    uint chrom_id                          [[ threadgroup_position_in_grid ]],
    uint blockDim                          [[ threads_per_threadgroup ]])
{
    const int n_lig   = p.lig_last - p.lig_first + 1;
    const int n_pro   = p.N - n_lig;
    const int n_pairs = n_lig * n_pro;

    const int gbase = int(chrom_id) * p.n_genes;
    const float tx = genes_f[gbase + 0];
    const float ty = genes_f[gbase + 1];
    const float tz = genes_f[gbase + 2];

    for (int la = int(tid); la < n_lig && la < 256; la += int(blockDim)) {
        float ra  = atom_radius[p.lig_first + la];
        float rwa = ra + 1.4f;
        lig_sas[la] = 4.0f * M_PI_F * rwa * rwa;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    const bool do_dw    = (p.dw_r0 > 0.0f);
    const float inv_r0  = do_dw ? (1.0f / p.dw_r0) : 0.0f;
    const bool do_elec  = (p.elec_on != 0) && (p.has_charge != 0);
    const bool do_hbond = (p.hbond_weight != 0.0f || p.hbond_salt_weight != 0.0f)
                          && (p.has_hflags != 0) && (p.has_charge != 0);
    const bool do_pb    = (p.pb_clash_weight > 0.0f) && (p.has_pbvdw != 0);

    float local_com = 0.0f, local_wal = 0.0f;
    float local_elec = 0.0f, local_hbond = 0.0f, local_pb = 0.0f;

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
        const float outer_r = rsum + 2.8f;

        float rel_area = 0.0f;
        if      (r < rsum)    rel_area = 1.0f;
        else if (r < outer_r) rel_area = 1.0f - (r - rsum) / (outer_r - rsum);
        const bool in_contact = (rel_area > 0.0f);

        if (in_contact && li < 256)
            tg_atomic_sub_float(&lig_sas[li], rel_area * surf_A);

        const int ti  = atom_type[ai];
        const int tj  = atom_type[aj];
        const float yval = gpu_get_yval(emat_sampled, ti, tj, p.T, rel_area);
        float com_pair = yval * rel_area;
        if (do_dw && in_contact) com_pair *= exp(-r * inv_r0);
        local_com += com_pair;

        if (do_elec && in_contact && r > 0.5f) {
            const float qA = atom_charge[ai];
            const float qB = atom_charge[aj];
            if (qA != 0.0f && qB != 0.0f)
                local_elec += KCOULOMB_F * qA * qB / (p.dielectric * r * r);
        }

        if (do_hbond && in_contact) {
            const int fa = atom_hflags[ai];
            const int fb = atom_hflags[aj];
            const bool a_to_b = ((fa & 0x1) != 0) && ((fb & 0x2) != 0);
            const bool b_to_a = ((fb & 0x1) != 0) && ((fa & 0x2) != 0);
            const float qa = atom_charge[ai];
            const float qb = atom_charge[aj];
            const bool salt = (qa <= -Q_SALT && qb >= Q_SALT) ||
                              (qb <= -Q_SALT && qa >= Q_SALT);
            if (a_to_b || b_to_a || salt) {
                const float dd = (r - p.hbond_opt_dist) / p.hbond_sigma_dist;
                const float E_dist = exp(-0.5f * dd * dd);
                const float w = salt ? p.hbond_salt_weight : p.hbond_weight;
                float E = w * E_dist * p.hbond_angle_repr;
                if (E < HBOND_PAIR_MIN) E = HBOND_PAIR_MIN;
                local_hbond += E;
            }
        }

        const float clash_r = p.perm * rsum;
        if (r < clash_r && r > 0.0f) {
            const float r2 = r*r; const float r4 = r2*r2; const float r6 = r4*r2;
            const float inv_r12 = 1.0f / (r6 * r6);
            const float cr2 = clash_r*clash_r; const float cr4 = cr2*cr2; const float cr6 = cr4*cr2;
            const float inv_cr12 = 1.0f / (cr6 * cr6);
            local_wal += KWALL_F * (inv_r12 - inv_cr12);
        }

        if (do_pb) {
            const float cr_pb = p.pb_clash_ratio * (atom_pbvdw[ai] + atom_pbvdw[aj]);
            const float o = cr_pb - r;
            if (o > 0.0f && r > 1.0e-6f)
                local_pb += pow(o, p.pb_clash_exponent);
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float com   = tg_reduce_sum(red, tid, local_com);
    float wal   = tg_reduce_sum(red, tid, local_wal);
    float elec  = tg_reduce_sum(red, tid, local_elec);
    float hbond = tg_reduce_sum(red, tid, local_hbond);
    float pb    = tg_reduce_sum(red, tid, local_pb);

    const bool do_gist = (p.gist_nx > 0);
    float local_sas = 0.0f, local_gist = 0.0f;
    float local_cx = 0.0f, local_cy = 0.0f, local_cz = 0.0f;
    for (int la = int(tid); la < n_lig && la < 256; la += int(blockDim)) {
        const int aidx = p.lig_first + la;
        const float sas_rem = max(0.0f, lig_sas[la]);
        const float rwa_la  = atom_radius[aidx] + 1.4f;
        const float surf_la = 4.0f * M_PI_F * rwa_la * rwa_la;
        const float sas_norm = sas_rem / surf_la;
        float contribution;
        if (p.solvent_flat != 0) {
            contribution = p.solventterm * sas_rem;
        } else {
            const int ti_la = atom_type[aidx];
            const float yval_sas = gpu_get_yval(emat_sampled, ti_la, p.T - 1, p.T, sas_norm);
            contribution = yval_sas * sas_norm;
        }
        local_sas += p.sas_weight * contribution;

        const float lx = atom_xyz[aidx*3+0] + tx;
        const float ly = atom_xyz[aidx*3+1] + ty;
        const float lz = atom_xyz[aidx*3+2] + tz;
        local_cx += lx; local_cy += ly; local_cz += lz;

        if (do_gist) {
            float fx = (lx - p.gox) * p.gidx;
            float fy = (ly - p.goy) * p.gidy;
            float fz = (lz - p.goz) * p.gidz;
            if (fx >= 0.0f && fx < float(p.gist_nx - 1) &&
                fy >= 0.0f && fy < float(p.gist_ny - 1) &&
                fz >= 0.0f && fz < float(p.gist_nz - 1)) {
                int ix = int(fx), iy = int(fy), iz = int(fz);
                float ddx = fx - ix, ddy = fy - iy, ddz = fz - iz;
                int ny = p.gist_ny, nz = p.gist_nz;
                #define GIDX(i,j,k) (((i)*ny + (j))*nz + (k))
                float c000 = gist_data[GIDX(ix,   iy,   iz  )];
                float c001 = gist_data[GIDX(ix,   iy,   iz+1)];
                float c010 = gist_data[GIDX(ix,   iy+1, iz  )];
                float c011 = gist_data[GIDX(ix,   iy+1, iz+1)];
                float c100 = gist_data[GIDX(ix+1, iy,   iz  )];
                float c101 = gist_data[GIDX(ix+1, iy,   iz+1)];
                float c110 = gist_data[GIDX(ix+1, iy+1, iz  )];
                float c111 = gist_data[GIDX(ix+1, iy+1, iz+1)];
                #undef GIDX
                float c00 = c000*(1.0f-ddz) + c001*ddz;
                float c01 = c010*(1.0f-ddz) + c011*ddz;
                float c10 = c100*(1.0f-ddz) + c101*ddz;
                float c11 = c110*(1.0f-ddz) + c111*ddz;
                float c0  = c00*(1.0f-ddy) + c01*ddy;
                float c1  = c10*(1.0f-ddy) + c11*ddy;
                local_gist += p.gist_weight * (c0*(1.0f-ddx) + c1*ddx);
            }
        }
    }
    float sas    = tg_reduce_sum(red, tid, local_sas);
    float gist   = tg_reduce_sum(red, tid, local_gist);
    float cx_sum = tg_reduce_sum(red, tid, local_cx);
    float cy_sum = tg_reduce_sum(red, tid, local_cy);
    float cz_sum = tg_reduce_sum(red, tid, local_cz);

    if (p.pb_pocket_weight > 0.0f && p.has_hflags != 0 && n_lig > 0) {
        const float inv_nl = 1.0f / float(n_lig);
        const float gx = cx_sum * inv_nl;
        const float gy = cy_sum * inv_nl;
        const float gz = cz_sum * inv_nl;
        float local_min = 3.4e38f;
        for (int a = int(tid); a < p.N; a += int(blockDim)) {
            if (a >= p.lig_first && a <= p.lig_last) continue;
            if ((atom_hflags[a] & 0x4) == 0) continue;
            const float ddx = gx - atom_xyz[a*3+0];
            const float ddy = gy - atom_xyz[a*3+1];
            const float ddz = gz - atom_xyz[a*3+2];
            local_min = min(local_min, ddx*ddx + ddy*ddy + ddz*ddz);
        }
        float best_d2 = tg_reduce_min(red, tid, local_min);
        if (best_d2 < 3.0e38f) {
            float d = sqrt(best_d2);
            float over = d - p.pb_pocket_radius;
            if (over > 0.0f) pb += p.pb_pocket_weight * over * over;
        }
    }

    float con = 0.0f;
    if (tid == 0 && p.n_cons > 0) {
        for (int c = 0; c < p.n_cons; ++c) {
            int i = cons_i[c], j = cons_j[c];
            float ix = atom_xyz[i*3+0], iy = atom_xyz[i*3+1], iz = atom_xyz[i*3+2];
            float jx = atom_xyz[j*3+0], jy = atom_xyz[j*3+1], jz = atom_xyz[j*3+2];
            if (i >= p.lig_first && i <= p.lig_last) { ix += tx; iy += ty; iz += tz; }
            if (j >= p.lig_first && j <= p.lig_last) { jx += tx; jy += ty; jz += tz; }
            float ddx = ix-jx, ddy = iy-jy, ddz = iz-jz;
            float d = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
            con += p.kdist - p.kdist * gpu_gaussian(d, cons_bondlen[c], cons_maxdist[c]);
        }
    }

    if (tid == 0) {
        int mp = max_pop;
        out[0*mp + int(chrom_id)] = com;
        out[1*mp + int(chrom_id)] = wal;
        out[2*mp + int(chrom_id)] = sas;
        out[3*mp + int(chrom_id)] = con;
        out[4*mp + int(chrom_id)] = elec;
        out[5*mp + int(chrom_id)] = hbond;
        out[6*mp + int(chrom_id)] = gist;
        out[7*mp + int(chrom_id)] = pb;
    }
}

// ─── Multi-complex screening kernel (com dist-weighted / wal / sas only) ──────
struct MultiParams {
    int   pop_size, n_genes, T;
    float perm;
    float dw_r0;
    float sas_weight;
    int   solvent_flat;
    float solventterm;
};

struct ComplexDesc {
    int atom_offset, n_atoms, lig_first, lig_last;
    int gene_offset, result_offset, pad0, pad1;
};

kernel void kernel_eval_cf_multi(
    device const float*        atom_xyz_all    [[ buffer(0) ]],
    device const int*          atom_type_all   [[ buffer(1) ]],
    device const float*        atom_radius_all [[ buffer(2) ]],
    device const float*        emat_sampled    [[ buffer(3) ]],
    device const float*        genes_f_all     [[ buffer(4) ]],
    device float*              cf_com_out      [[ buffer(5) ]],
    device float*              cf_wal_out      [[ buffer(6) ]],
    device float*              cf_sas_out      [[ buffer(7) ]],
    constant MultiParams&      mp              [[ buffer(8) ]],
    device const ComplexDesc*  descs           [[ buffer(9) ]],
    threadgroup float*         lig_sas         [[ threadgroup(0) ]],
    threadgroup float*         red             [[ threadgroup(1) ]],
    uint tid                                   [[ thread_position_in_threadgroup ]],
    uint global_chrom_id                       [[ threadgroup_position_in_grid ]],
    uint blockDim                              [[ threads_per_threadgroup ]])
{
    const int ci  = int(global_chrom_id) / mp.pop_size;
    const int chi = int(global_chrom_id) % mp.pop_size;
    device const ComplexDesc& d = descs[ci];

    const int n_lig   = d.lig_last - d.lig_first + 1;
    const int n_pro   = d.n_atoms - n_lig;
    const int n_pairs = n_lig * n_pro;

    const int gbase = d.gene_offset + chi * mp.n_genes;
    const float tx = genes_f_all[gbase + 0];
    const float ty = genes_f_all[gbase + 1];
    const float tz = genes_f_all[gbase + 2];

    for (int la = int(tid); la < n_lig && la < 256; la += int(blockDim)) {
        float ra  = atom_radius_all[d.atom_offset + d.lig_first + la];
        float rwa = ra + 1.4f;
        lig_sas[la] = 4.0f * M_PI_F * rwa * rwa;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    const bool do_dw   = (mp.dw_r0 > 0.0f);
    const float inv_r0 = do_dw ? (1.0f / mp.dw_r0) : 0.0f;

    float local_com = 0.0f, local_wal = 0.0f;

    for (int pr = int(tid); pr < n_pairs; pr += int(blockDim)) {
        const int li      = pr / n_pro;
        const int pro_rel = pr % n_pro;
        const int ai      = d.atom_offset + d.lig_first + li;
        const int pro_loc = (pro_rel < d.lig_first) ? pro_rel : (pro_rel + n_lig);
        const int aj      = d.atom_offset + pro_loc;

        const float lx = atom_xyz_all[ai*3+0] + tx;
        const float ly = atom_xyz_all[ai*3+1] + ty;
        const float lz = atom_xyz_all[ai*3+2] + tz;
        const float dx = lx - atom_xyz_all[aj*3+0];
        const float dy = ly - atom_xyz_all[aj*3+1];
        const float dz = lz - atom_xyz_all[aj*3+2];
        const float r  = sqrt(dx*dx + dy*dy + dz*dz + 1e-10f);

        const float rA = atom_radius_all[ai];
        const float rB = atom_radius_all[aj];
        const float rsum = rA + rB;
        const float rwa_A = rA + 1.4f;
        const float surf_A = 4.0f * M_PI_F * rwa_A * rwa_A;
        const float outer_r = rsum + 2.8f;

        float rel_area = 0.0f;
        if      (r < rsum)    rel_area = 1.0f;
        else if (r < outer_r) rel_area = 1.0f - (r - rsum) / (outer_r - rsum);
        const bool in_contact = (rel_area > 0.0f);

        if (in_contact && li < 256)
            tg_atomic_sub_float(&lig_sas[li], rel_area * surf_A);

        const int ti = atom_type_all[ai];
        const int tj = atom_type_all[aj];
        const float yval = gpu_get_yval(emat_sampled, ti, tj, mp.T, rel_area);
        float com_pair = yval * rel_area;
        if (do_dw && in_contact) com_pair *= exp(-r * inv_r0);
        local_com += com_pair;

        const float clash_r = mp.perm * rsum;
        if (r < clash_r && r > 0.0f) {
            const float r2=r*r; const float r4=r2*r2; const float r6=r4*r2;
            const float inv_r12 = 1.0f/(r6*r6);
            const float cr2=clash_r*clash_r; const float cr4=cr2*cr2; const float cr6=cr4*cr2;
            const float inv_cr12 = 1.0f/(cr6*cr6);
            local_wal += KWALL_F * (inv_r12 - inv_cr12);
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float com = tg_reduce_sum(red, tid, local_com);
    float wal = tg_reduce_sum(red, tid, local_wal);

    float local_sas = 0.0f;
    for (int la = int(tid); la < n_lig && la < 256; la += int(blockDim)) {
        const int aidx = d.atom_offset + d.lig_first + la;
        const float sas_rem = max(0.0f, lig_sas[la]);
        const float rwa_la  = atom_radius_all[aidx] + 1.4f;
        const float surf_la = 4.0f * M_PI_F * rwa_la * rwa_la;
        const float sas_norm = sas_rem / surf_la;
        float contribution;
        if (mp.solvent_flat != 0) {
            contribution = mp.solventterm * sas_rem;
        } else {
            const int ti_la = atom_type_all[aidx];
            const float yval_sas = gpu_get_yval(emat_sampled, ti_la, mp.T - 1, mp.T, sas_norm);
            contribution = yval_sas * sas_norm;
        }
        local_sas += mp.sas_weight * contribution;
    }
    float sas = tg_reduce_sum(red, tid, local_sas);

    if (tid == 0) {
        const int out = d.result_offset + chi;
        cf_com_out[out] = com;
        cf_wal_out[out] = wal;
        cf_sas_out[out] = sas;
    }
}
)MSL";

// ─── context structure ────────────────────────────────────────────────────────
struct MetalEvalCtx {
    id<MTLDevice>              device;
    id<MTLCommandQueue>        queue;
    id<MTLComputePipelineState> pipeline;        // single-complex kernel
    id<MTLComputePipelineState> pipeline_multi;  // multi-complex kernel

    id<MTLBuffer> buf_atom_xyz;
    id<MTLBuffer> buf_atom_type;
    id<MTLBuffer> buf_atom_radius;
    id<MTLBuffer> buf_emat_sampled;
    id<MTLBuffer> buf_genes_f;
    id<MTLBuffer> buf_out;          // single merged output [CH * max_pop]

    // Extra full-fidelity static buffers (always allocated; dummy when absent).
    id<MTLBuffer> buf_pbvdw;
    id<MTLBuffer> buf_charge;
    id<MTLBuffer> buf_hflags;
    id<MTLBuffer> buf_cons_i;
    id<MTLBuffer> buf_cons_j;
    id<MTLBuffer> buf_cons_bondlen;
    id<MTLBuffer> buf_cons_maxdist;
    id<MTLBuffer> buf_gist;

    int   has_pbvdw, has_charge, has_hflags;
    int   n_cons;
    int   gist_nx, gist_ny, gist_nz;
    float gist_origin[3];
    float gist_inv_delta[3];
    float gist_weight;

    GpuCfParams cur_params;   // stashed by metal_eval_set_params()
    bool        have_params;

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
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return device != nil;
    }
}

void metal_eval_get_capabilities(MetalCapabilities* out)
{
    if (!out) return;
    memset(out, 0, sizeof(*out));
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) { out->available = false; return; }
        out->available = true;
        const char* name_cstr = device.name.UTF8String ? device.name.UTF8String : "Apple GPU";
        strncpy(out->device_name, name_cstr, sizeof(out->device_name) - 1);
#if defined(__MAC_OS_X_VERSION_MIN_REQUIRED) && __MAC_OS_X_VERSION_MIN_REQUIRED >= 101500
        if (@available(macOS 10.15, *)) {
            out->unified_memory_bytes = device.hasUnifiedMemory ? device.maxBufferLength : 0;
        }
#endif
        out->max_buffer_length = device.maxBufferLength;
        NSString* name = device.name;
        if ([name containsString:@"M3 Pro"])       out->gpu_core_estimate = 18;
        else if ([name containsString:@"M3 Max"])  out->gpu_core_estimate = 30;
        else if ([name containsString:@"M3"])      out->gpu_core_estimate = 10;
        else                                        out->gpu_core_estimate = 8;
    }
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
    memset(ctx, 0, sizeof(*ctx));
    ctx->n_atoms   = n_atoms;
    ctx->n_types   = n_types;
    ctx->max_pop   = max_pop;
    ctx->lig_first = lig_first;
    ctx->lig_last  = lig_last;
    ctx->perm      = perm;

    ctx->device = MTLCreateSystemDefaultDevice();
    if (!ctx->device) {
        fprintf(stderr, "metal_eval: no Metal device found\n");
        delete ctx;
        return nullptr;
    }
    ctx->queue = [ctx->device newCommandQueue];

    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:kMSLSource];
    id<MTLLibrary> lib = [ctx->device newLibraryWithSource:src options:nil error:&err];
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
    id<MTLFunction> fn_multi = [lib newFunctionWithName:@"kernel_eval_cf_multi"];
    if (fn_multi) {
        ctx->pipeline_multi = [ctx->device newComputePipelineStateWithFunction:fn_multi error:&err];
        if (!ctx->pipeline_multi)
            fprintf(stderr, "metal_eval: pipeline_multi compile warning: %s\n",
                    err ? [[err localizedDescription] UTF8String] : "unknown");
    }

    auto newBuf = [&](const void* data, size_t bytes) -> id<MTLBuffer> {
        return [ctx->device newBufferWithBytes:data length:bytes
                                       options:MTLResourceStorageModeShared];
    };
    auto newLen = [&](size_t bytes) -> id<MTLBuffer> {
        return [ctx->device newBufferWithLength:(bytes ? bytes : sizeof(float))
                                       options:MTLResourceStorageModeShared];
    };

    ctx->buf_atom_xyz    = newBuf(h_atom_xyz,    (size_t)n_atoms * 3 * sizeof(float));
    ctx->buf_atom_type   = newBuf(h_atom_type,   (size_t)n_atoms     * sizeof(int));
    ctx->buf_atom_radius = newBuf(h_atom_radius, (size_t)n_atoms     * sizeof(float));
    ctx->buf_emat_sampled= newBuf(h_emat_sampled,
                                  (size_t)n_types * n_types * n_emat_samples * sizeof(float));

    ctx->max_genes = 256;
    ctx->buf_genes_f = newLen((size_t)max_pop * ctx->max_genes * sizeof(float));
    ctx->buf_out     = newLen((size_t)METAL_CF_CHANNELS * max_pop * sizeof(float));

    // Dummy extra buffers by default (disabled terms never index them).
    ctx->buf_pbvdw        = newLen(0);
    ctx->buf_charge       = newLen(0);
    ctx->buf_hflags       = newLen(0);
    ctx->buf_cons_i       = newLen(0);
    ctx->buf_cons_j       = newLen(0);
    ctx->buf_cons_bondlen = newLen(0);
    ctx->buf_cons_maxdist = newLen(0);
    ctx->buf_gist         = newLen(0);

    return ctx;
}

void metal_eval_set_extra(MetalEvalCtx* ctx, const GpuCfExtraStatic* extra)
{
    if (!ctx || !extra) return;
    const int N = ctx->n_atoms;
    auto newBuf = [&](const void* data, size_t bytes) -> id<MTLBuffer> {
        return [ctx->device newBufferWithBytes:data length:bytes
                                       options:MTLResourceStorageModeShared];
    };

    if (extra->atom_pbvdw) {
        ctx->buf_pbvdw  = newBuf(extra->atom_pbvdw,  (size_t)N * sizeof(float));
        ctx->has_pbvdw  = 1;
    }
    if (extra->atom_charge) {
        ctx->buf_charge = newBuf(extra->atom_charge, (size_t)N * sizeof(float));
        ctx->has_charge = 1;
    }
    if (extra->atom_hflags) {
        ctx->buf_hflags = newBuf(extra->atom_hflags, (size_t)N * sizeof(int));
        ctx->has_hflags = 1;
    }

    ctx->n_cons = extra->n_cons;
    if (extra->n_cons > 0 && extra->cons_i) {
        ctx->buf_cons_i       = newBuf(extra->cons_i,       (size_t)extra->n_cons * sizeof(int));
        ctx->buf_cons_j       = newBuf(extra->cons_j,       (size_t)extra->n_cons * sizeof(int));
        ctx->buf_cons_bondlen = newBuf(extra->cons_bondlen, (size_t)extra->n_cons * sizeof(float));
        ctx->buf_cons_maxdist = newBuf(extra->cons_maxdist, (size_t)extra->n_cons * sizeof(float));
    }

    ctx->gist_nx = extra->gist_nx;
    ctx->gist_ny = extra->gist_ny;
    ctx->gist_nz = extra->gist_nz;
    for (int k = 0; k < 3; ++k) {
        ctx->gist_origin[k]    = extra->gist_origin[k];
        ctx->gist_inv_delta[k] = extra->gist_inv_delta[k];
    }
    ctx->gist_weight = extra->gist_weight;
    if (extra->gist_nx > 0 && extra->gist_data) {
        size_t vox = (size_t)extra->gist_nx * extra->gist_ny * extra->gist_nz * sizeof(float);
        ctx->buf_gist = newBuf(extra->gist_data, vox);
    }
}

void metal_eval_set_params(MetalEvalCtx* ctx, const GpuCfParams* params)
{
    if (!ctx || !params) return;
    ctx->cur_params  = *params;
    ctx->have_params = true;
}

static void fill_kp(MetalKP& kp, const MetalEvalCtx* ctx,
                    int n_genes, const GpuCfParams& P)
{
    kp.N = ctx->n_atoms; kp.T = ctx->n_types; kp.n_genes = n_genes;
    kp.lig_first = ctx->lig_first; kp.lig_last = ctx->lig_last;
    kp.perm = ctx->perm;
    kp.has_pbvdw = ctx->has_pbvdw; kp.has_charge = ctx->has_charge; kp.has_hflags = ctx->has_hflags;
    kp.n_cons = ctx->n_cons;
    kp.gist_nx = ctx->gist_nx; kp.gist_ny = ctx->gist_ny; kp.gist_nz = ctx->gist_nz;
    kp.gox = ctx->gist_origin[0]; kp.goy = ctx->gist_origin[1]; kp.goz = ctx->gist_origin[2];
    kp.gidx = ctx->gist_inv_delta[0]; kp.gidy = ctx->gist_inv_delta[1]; kp.gidz = ctx->gist_inv_delta[2];
    kp.gist_weight = ctx->gist_weight;
    kp.dw_r0 = P.dw_r0; kp.elec_on = P.elec_on; kp.dielectric = P.dielectric;
    kp.hbond_weight = P.hbond_weight; kp.hbond_salt_weight = P.hbond_salt_weight;
    kp.hbond_opt_dist = P.hbond_opt_dist; kp.hbond_sigma_dist = P.hbond_sigma_dist;
    kp.hbond_angle_repr = P.hbond_angle_repr;
    kp.pb_clash_weight = P.pb_clash_weight; kp.pb_clash_ratio = P.pb_clash_ratio;
    kp.pb_clash_exponent = P.pb_clash_exponent;
    kp.pb_pocket_weight = P.pb_pocket_weight; kp.pb_pocket_radius = P.pb_pocket_radius;
    kp.kdist = P.kdist; kp.sas_weight = P.sas_weight;
    kp.solvent_flat = P.solvent_flat; kp.solventterm = P.solventterm;
}

void metal_eval_batch(MetalEvalCtx*       ctx,
                      int                 pop_size,
                      int                 n_genes,
                      const double*       h_genes,
                      const GpuCfParams*  params,
                      const GpuCfResults* out)
{
    if (pop_size > ctx->max_pop) {
        throw FlexAIDException("metal_eval: pop_size " + std::to_string(pop_size) +
            " exceeds max_pop " + std::to_string(ctx->max_pop));
    }
    if (n_genes > ctx->max_genes) {
        throw FlexAIDException("metal_eval: n_genes " + std::to_string(n_genes) +
            " exceeds max_genes " + std::to_string(ctx->max_genes));
    }

    float* genes_f = (float*)[ctx->buf_genes_f contents];
    for (int c = 0; c < pop_size; ++c)
        for (int g = 0; g < n_genes; ++g)
            genes_f[c * n_genes + g] = (float)h_genes[c * n_genes + g];

    MetalKP kp; fill_kp(kp, ctx, n_genes, *params);
    id<MTLBuffer> buf_kp = [ctx->device newBufferWithBytes:&kp length:sizeof(kp)
                                                   options:MTLResourceStorageModeShared];
    int mp_val = ctx->max_pop;
    id<MTLBuffer> buf_mp = [ctx->device newBufferWithBytes:&mp_val length:sizeof(int)
                                                   options:MTLResourceStorageModeShared];

    id<MTLCommandBuffer>         cb  = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:ctx->pipeline];
    [enc setBuffer:ctx->buf_atom_xyz     offset:0 atIndex:0];
    [enc setBuffer:ctx->buf_atom_type    offset:0 atIndex:1];
    [enc setBuffer:ctx->buf_atom_radius  offset:0 atIndex:2];
    [enc setBuffer:ctx->buf_emat_sampled offset:0 atIndex:3];
    [enc setBuffer:ctx->buf_genes_f      offset:0 atIndex:4];
    [enc setBuffer:ctx->buf_out          offset:0 atIndex:5];
    [enc setBuffer:buf_kp                offset:0 atIndex:6];
    [enc setBuffer:ctx->buf_pbvdw        offset:0 atIndex:7];
    [enc setBuffer:ctx->buf_charge       offset:0 atIndex:8];
    [enc setBuffer:ctx->buf_hflags       offset:0 atIndex:9];
    [enc setBuffer:ctx->buf_cons_i       offset:0 atIndex:10];
    [enc setBuffer:ctx->buf_cons_j       offset:0 atIndex:11];
    [enc setBuffer:ctx->buf_cons_bondlen offset:0 atIndex:12];
    [enc setBuffer:ctx->buf_cons_maxdist offset:0 atIndex:13];
    [enc setBuffer:ctx->buf_gist         offset:0 atIndex:14];
    [enc setBuffer:buf_mp                offset:0 atIndex:15];
    [enc setThreadgroupMemoryLength:256 * sizeof(float) atIndex:0]; // lig_sas
    [enc setThreadgroupMemoryLength:256 * sizeof(float) atIndex:1]; // reduction scratch

    MTLSize gridSize  = { (NSUInteger)pop_size, 1, 1 };
    MTLSize groupSize = { 256, 1, 1 };
    [enc dispatchThreadgroups:gridSize threadsPerThreadgroup:groupSize];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];

    if (cb.status == MTLCommandBufferStatusError) {
        NSString* desc = cb.error ? cb.error.localizedDescription : @"unknown";
        throw FlexAIDException("metal_eval: GPU command buffer error: " +
            std::string([desc UTF8String]));
    }

    const float* o = (const float*)[ctx->buf_out contents];
    const int mp = ctx->max_pop;
    auto ch = [&](int c) { return o + (size_t)c * mp; };
    for (int c = 0; c < pop_size; ++c) {
        if (out->com)         out->com[c]         = (double)ch(0)[c];
        if (out->wal)         out->wal[c]         = (double)ch(1)[c];
        if (out->sas)         out->sas[c]         = (double)ch(2)[c];
        if (out->con)         out->con[c]         = (double)ch(3)[c];
        if (out->elec)        out->elec[c]        = (double)ch(4)[c];
        if (out->hbond)       out->hbond[c]       = (double)ch(5)[c];
        if (out->gist_desolv) out->gist_desolv[c] = (double)ch(6)[c];
        if (out->pb_clash)    out->pb_clash[c]    = (double)ch(7)[c];
    }
}

void metal_eval_shutdown(MetalEvalCtx* ctx)
{
    if (!ctx) return;
    // ARC-managed objects are released automatically.
    delete ctx;
}

// ─── Multi-complex batched evaluation (screening pre-filter) ─────────────────
void metal_eval_batch_multi(MetalEvalCtx*               ctx,
                             int                         n_complex,
                             int                         pop_size,
                             int                         n_genes,
                             const MetalMultiBatchEntry* entries)
{
    if (n_complex <= 0) return;

    // Params: use whatever was stashed (default-zero = all extra terms off).
    GpuCfParams P;
    if (ctx->have_params) P = ctx->cur_params;
    else { memset(&P, 0, sizeof(P)); P.sas_weight = 1.0f; }

    // Fast path / no-multi-kernel fallback: delegate to the full-fidelity single
    // batch. (Only com/wal/sas outputs exist on the entry, so the extra channels
    // are dropped — consistent with the multi contract.)
    if (n_complex == 1 || !ctx->pipeline_multi) {
        for (int ci = 0; ci < n_complex; ++ci) {
            GpuCfResults r; memset(&r, 0, sizeof(r));
            r.com = entries[ci].h_com_out;
            r.wal = entries[ci].h_wal_out;
            r.sas = entries[ci].h_sas_out;
            metal_eval_batch(ctx, pop_size, n_genes, entries[ci].h_genes, &P, &r);
        }
        return;
    }

    const int total_pop = n_complex * pop_size;

    int total_atoms = 0;
    for (int ci = 0; ci < n_complex; ++ci) total_atoms += entries[ci].n_atoms;

    std::vector<float> xyz_all (total_atoms * 3);
    std::vector<int>   type_all(total_atoms);
    std::vector<float> rad_all (total_atoms);

    struct ComplexDescHost {
        int atom_offset, n_atoms, lig_first, lig_last;
        int gene_offset, result_offset, pad0, pad1;
    };
    std::vector<ComplexDescHost> descs(n_complex);

    int atom_off = 0;
    for (int ci = 0; ci < n_complex; ++ci) {
        const int na = entries[ci].n_atoms;
        memcpy(xyz_all.data()  + atom_off * 3, entries[ci].h_atom_xyz,    na * 3 * sizeof(float));
        memcpy(type_all.data() + atom_off,     entries[ci].h_atom_type,   na     * sizeof(int));
        memcpy(rad_all.data()  + atom_off,     entries[ci].h_atom_radius, na     * sizeof(float));
        descs[ci] = { atom_off, na, entries[ci].lig_first, entries[ci].lig_last,
                      ci * pop_size * n_genes, ci * pop_size, 0, 0 };
        atom_off += na;
    }

    std::vector<float> genes_all(total_pop * n_genes);
    for (int ci = 0; ci < n_complex; ++ci) {
        float* dst = genes_all.data() + ci * pop_size * n_genes;
        const double* src = entries[ci].h_genes;
        for (int c = 0; c < pop_size; ++c)
            for (int g = 0; g < n_genes; ++g)
                dst[c * n_genes + g] = (float)src[c * n_genes + g];
    }

    auto mk = [&](const void* d, size_t n) {
        return [ctx->device newBufferWithBytes:d length:(n ? n : sizeof(float))
                                      options:MTLResourceStorageModeShared];
    };
    id<MTLBuffer> buf_xyz   = mk(xyz_all.data(),  xyz_all.size()  * sizeof(float));
    id<MTLBuffer> buf_type  = mk(type_all.data(), type_all.size() * sizeof(int));
    id<MTLBuffer> buf_rad   = mk(rad_all.data(),  rad_all.size()  * sizeof(float));
    id<MTLBuffer> buf_genes = mk(genes_all.data(),genes_all.size()* sizeof(float));
    id<MTLBuffer> buf_com   = [ctx->device newBufferWithLength:total_pop*sizeof(float)
                                                       options:MTLResourceStorageModeShared];
    id<MTLBuffer> buf_wal   = [ctx->device newBufferWithLength:total_pop*sizeof(float)
                                                       options:MTLResourceStorageModeShared];
    id<MTLBuffer> buf_sas   = [ctx->device newBufferWithLength:total_pop*sizeof(float)
                                                       options:MTLResourceStorageModeShared];

    MetalMultiParams mp;
    mp.pop_size = pop_size; mp.n_genes = n_genes; mp.T = entries[0].n_types;
    mp.perm = entries[0].perm;
    mp.dw_r0 = P.dw_r0; mp.sas_weight = P.sas_weight;
    mp.solvent_flat = P.solvent_flat; mp.solventterm = P.solventterm;
    id<MTLBuffer> buf_mp    = mk(&mp, sizeof(mp));
    id<MTLBuffer> buf_descs = mk(descs.data(), descs.size() * sizeof(ComplexDescHost));

    id<MTLCommandBuffer>         cb  = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:ctx->pipeline_multi];
    [enc setBuffer:buf_xyz               offset:0 atIndex:0];
    [enc setBuffer:buf_type              offset:0 atIndex:1];
    [enc setBuffer:buf_rad               offset:0 atIndex:2];
    [enc setBuffer:ctx->buf_emat_sampled offset:0 atIndex:3];
    [enc setBuffer:buf_genes             offset:0 atIndex:4];
    [enc setBuffer:buf_com               offset:0 atIndex:5];
    [enc setBuffer:buf_wal               offset:0 atIndex:6];
    [enc setBuffer:buf_sas               offset:0 atIndex:7];
    [enc setBuffer:buf_mp                offset:0 atIndex:8];
    [enc setBuffer:buf_descs             offset:0 atIndex:9];
    [enc setThreadgroupMemoryLength:256 * sizeof(float) atIndex:0];
    [enc setThreadgroupMemoryLength:256 * sizeof(float) atIndex:1];

    MTLSize gridSize  = { (NSUInteger)total_pop, 1, 1 };
    MTLSize groupSize = { 256, 1, 1 };
    [enc dispatchThreadgroups:gridSize threadsPerThreadgroup:groupSize];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];

    if (cb.status == MTLCommandBufferStatusError) {
        NSString* desc = cb.error ? cb.error.localizedDescription : @"unknown";
        throw FlexAIDException("metal_eval_batch_multi: GPU error: " +
            std::string([desc UTF8String]));
    }

    const float* com_f = (const float*)[buf_com contents];
    const float* wal_f = (const float*)[buf_wal contents];
    const float* sas_f = (const float*)[buf_sas contents];
    for (int ci = 0; ci < n_complex; ++ci) {
        const int base = ci * pop_size;
        for (int c = 0; c < pop_size; ++c) {
            entries[ci].h_com_out[c] = (double)com_f[base + c];
            entries[ci].h_wal_out[c] = (double)wal_f[base + c];
            entries[ci].h_sas_out[c] = (double)sas_f[base + c];
        }
    }
}

#endif  // FLEXAIDS_USE_METAL
