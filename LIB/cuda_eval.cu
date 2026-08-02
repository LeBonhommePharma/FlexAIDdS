// cuda_eval.cu — Full-fidelity CUDA batched chromosome evaluation kernels
//
// Architecture:
//   Grid:  pop_size threadblocks  (one chromosome per block)
//   Block: 256 threads            (cooperative reduction over atom pairs)
//
// Each block evaluates one chromosome and produces the eight CF scoring
// channels that feed get_cf_evalue() on the host:
//     com, wal, sas, con, elec, hbond, gist_desolv, pb_clash
//
// See cuda_eval.cuh for the drift-tolerant contact-area model and the list of
// terms deliberately NOT reproduced on-device (angular H-bond directionality,
// metal-CN, vct-entropy, tENCoM), which the gaboom.cpp divergence guard warns
// about when their weights are non-zero.
//
// Performance characteristics (unchanged from the com/wal/sas baseline):
//   – Receptor atoms + energy matrix stay resident on-device (persistent ctx).
//   – Whole population batched per generation (one launch, pop_size blocks).
//   – CUDA stream for async H↔D transfers; pinned host staging buffers.
//   – Multiplication chains replace powf() in the r⁻¹² / severity terms.

#ifdef FLEXAIDS_USE_CUDA

#include "cuda_eval.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// ─── error-checking macro ─────────────────────────────────────────────────────
#include "flexaid_exception.h"
#include <string>
#define CUDA_CHECK(call) do {                                             \
    cudaError_t _e = (call);                                              \
    if (_e != cudaSuccess) {                                              \
        throw FlexAIDException(std::string("CUDA error at ") +           \
            __FILE__ + ":" + std::to_string(__LINE__) + "  " +           \
            cudaGetErrorString(_e));                                      \
    }                                                                     \
} while (0)

// ─── constants ────────────────────────────────────────────────────────────────
static constexpr int   BLOCK_SIZE     = 256;
static constexpr int   NWARPS         = BLOCK_SIZE / 32;
static constexpr int   N_EMAT_SAMPLES = 128;   // must match CUDA_EMAT_SAMPLES in .cuh
static constexpr float Rw             = 1.4f;  // water probe radius (Å)
static constexpr float KWALL_F        = 1.0e6f;
static constexpr float KCOULOMB_F     = 332.0637f;  // matches KCOULOMB in vcfunction.cpp
static constexpr float HBOND_PAIR_MIN = -2.0f;      // per-pair H-bond floor
static constexpr float Q_SALT         = 0.30f;      // salt-bridge charge threshold

// Number of CF channels reduced per block (com,wal,sas,con,elec,hbond,gist,pb).
static constexpr int   N_CF_CHANNELS  = 8;

// Maximum ligand atoms handled in shared-memory SAS/GIST accumulators.
// Ligands with more atoms fall back to zero contribution for the overflow atoms.
static constexpr int MAX_LIG_SAS = 512;

// ─── context ─────────────────────────────────────────────────────────────────
struct CudaEvalCtx {
    // Base rigid data
    float*  d_atom_xyz;      // [n_atoms × 3]
    int*    d_atom_type;     // [n_atoms]
    float*  d_atom_radius;   // [n_atoms]
    float*  d_emat_sampled;  // [n_types × n_types × N_EMAT_SAMPLES]
    double* d_genes;         // [max_pop × max_genes]

    // Extra full-fidelity static data (nullptr / 0 when unused)
    float*  d_atom_pbvdw;    // [n_atoms]
    float*  d_atom_charge;   // [n_atoms]
    int*    d_atom_hflags;   // [n_atoms]
    int     n_cons;
    int*    d_cons_i;        // [n_cons]
    int*    d_cons_j;        // [n_cons]
    float*  d_cons_bondlen;  // [n_cons]
    float*  d_cons_maxdist;  // [n_cons]
    int     gist_nx, gist_ny, gist_nz;
    float   gist_origin[3];
    float   gist_inv_delta[3];
    float   gist_weight;
    float*  d_gist_data;     // [nx*ny*nz]

    // Merged device output: N_CF_CHANNELS × max_pop doubles
    double* d_cf_out;

    // Pinned host buffers for async DMA
    double* h_genes_pinned;  // [max_pop × max_genes]
    double* h_cf_pinned;     // [N_CF_CHANNELS × max_pop]

    cudaStream_t stream;

    int   n_atoms;
    int   n_types;
    int   max_pop;
    int   max_genes;
    int   lig_first;
    int   lig_last;
    float perm;
};

// ─── device helpers ──────────────────────────────────────────────────────────

// Interpolated energy-matrix lookup (matches host get_yval on the sampled curve).
__device__ __forceinline__ float gpu_get_yval(
        const float* __restrict__ emat_sampled,
        int t1, int t2, int T, float rel_area)
{
    int   base = (t1 * T + t2) * N_EMAT_SAMPLES;
    rel_area   = fmaxf(0.0f, fminf(1.0f, rel_area));
    float kf   = rel_area * (N_EMAT_SAMPLES - 1.0f);
    int   k0   = (int)kf;
    int   k1   = min(k0 + 1, N_EMAT_SAMPLES - 1);
    float frac = kf - (float)k0;
    return emat_sampled[base + k0] * (1.0f - frac)
         + emat_sampled[base + k1] * frac;
}

// GetValueFromGaussian (geometry.cpp): pow(-(x-zero)*(x-(2*max-zero))/(zero-max)^2, 50).
// Base is clamped to >=0 to avoid NaN from a negative base at a non-integer-treated
// exponent (drift-safe; the CPU form is ~0 well outside the well anyway).
__device__ __forceinline__ float gpu_gaussian(float x, float mx, float zero)
{
    float dz   = zero - mx;
    float denom = dz * dz;
    if (denom < 1e-12f) return 0.0f;
    float base = (-(x - zero) * (x - (2.0f * mx - zero))) / denom;
    if (base <= 0.0f) return 0.0f;
    return powf(base, 50.0f);
}

// Block-wide sum reduction; result broadcast to all threads via smem[0].
__device__ __forceinline__ float block_reduce_sum(float v, float* smem, int tid)
{
    for (int o = warpSize / 2; o > 0; o >>= 1)
        v += __shfl_down_sync(0xFFFFFFFF, v, o);
    int lane = tid & 31, wid = tid >> 5;
    if (lane == 0) smem[wid] = v;
    __syncthreads();
    if (wid == 0) {
        float x = (lane < NWARPS) ? smem[lane] : 0.0f;
        for (int o = NWARPS / 2; o > 0; o >>= 1)
            x += __shfl_down_sync(0xFFFFFFFF, x, o);
        if (lane == 0) smem[0] = x;
    }
    __syncthreads();
    float r = smem[0];
    __syncthreads();
    return r;
}

// Block-wide min reduction; result broadcast to all threads via smem[0].
__device__ __forceinline__ float block_reduce_min(float v, float* smem, int tid)
{
    for (int o = warpSize / 2; o > 0; o >>= 1)
        v = fminf(v, __shfl_down_sync(0xFFFFFFFF, v, o));
    int lane = tid & 31, wid = tid >> 5;
    if (lane == 0) smem[wid] = v;
    __syncthreads();
    if (wid == 0) {
        float x = (lane < NWARPS) ? smem[lane] : 3.4e38f;
        for (int o = NWARPS / 2; o > 0; o >>= 1)
            x = fminf(x, __shfl_down_sync(0xFFFFFFFF, x, o));
        if (lane == 0) smem[0] = x;
    }
    __syncthreads();
    float r = smem[0];
    __syncthreads();
    return r;
}

// ─── full-fidelity CF kernel ──────────────────────────────────────────────────
__global__ void kernel_eval_cf_full(
    const float*  __restrict__ atom_xyz,       // [N × 3]
    const int*    __restrict__ atom_type,      // [N]
    const float*  __restrict__ atom_radius,    // [N]
    const float*  __restrict__ emat_sampled,   // [T × T × N_EMAT_SAMPLES]
    const double* __restrict__ genes,          // [pop × n_genes]
    // extra static (may be null)
    const float*  __restrict__ atom_pbvdw,     // [N] or null
    const float*  __restrict__ atom_charge,    // [N] or null
    const int*    __restrict__ atom_hflags,    // [N] or null
    int n_cons,
    const int*    __restrict__ cons_i,
    const int*    __restrict__ cons_j,
    const float*  __restrict__ cons_bondlen,
    const float*  __restrict__ cons_maxdist,
    int gist_nx, int gist_ny, int gist_nz,
    float gox, float goy, float goz,
    float gidx, float gidy, float gidz,
    float gist_weight,
    const float* __restrict__ gist_data,       // [nx*ny*nz] or null
    // params
    GpuCfParams P,
    // outputs (N_CF_CHANNELS × pop, contiguous per channel)
    double* __restrict__ out,
    int max_pop,
    int N, int T, int n_genes,
    int lig_first, int lig_last, float perm)
{
    const int chrom_id = blockIdx.x;
    const int tid      = threadIdx.x;

    // Translation from genes (first three genes encode tx, ty, tz).
    __shared__ float tx, ty, tz;
    if (tid == 0) {
        tx = (float)genes[chrom_id * n_genes + 0];
        ty = (float)genes[chrom_id * n_genes + 1];
        tz = (float)genes[chrom_id * n_genes + 2];
    }

    // Per-ligand SAS accumulator (initialised to 4π(rA+Rw)²).
    __shared__ float lig_sas[MAX_LIG_SAS];
    const int n_lig = lig_last - lig_first + 1;
    const int n_pro = N - n_lig;
    for (int la = tid; la < n_lig && la < MAX_LIG_SAS; la += BLOCK_SIZE) {
        float rwa   = atom_radius[lig_first + la] + Rw;
        lig_sas[la] = 4.0f * 3.141592653589793f * rwa * rwa;
    }
    __syncthreads();

    const bool do_dw    = (P.dw_r0 > 0.0f);
    const float inv_r0  = do_dw ? (1.0f / P.dw_r0) : 0.0f;
    const bool do_elec  = (P.elec_on != 0) && (atom_charge != nullptr);
    const bool do_hbond = (P.hbond_weight != 0.0f || P.hbond_salt_weight != 0.0f) &&
                          (atom_hflags != nullptr) && (atom_charge != nullptr);
    const bool do_pb    = (P.pb_clash_weight > 0.0f) && (atom_pbvdw != nullptr);

    // ── pair loop ────────────────────────────────────────────────────────────
    const int n_pairs = n_lig * n_pro;
    float local_com = 0.0f, local_wal = 0.0f;
    float local_elec = 0.0f, local_hbond = 0.0f, local_pb = 0.0f;

    for (int pr = tid; pr < n_pairs; pr += BLOCK_SIZE) {
        const int li      = pr / n_pro;
        const int pro_rel = pr % n_pro;
        const int ai      = lig_first + li;
        const int aj      = (pro_rel < lig_first) ? pro_rel : (pro_rel + n_lig);

        // Ligand atom position with translation applied.
        const float lx = atom_xyz[ai * 3 + 0] + tx;
        const float ly = atom_xyz[ai * 3 + 1] + ty;
        const float lz = atom_xyz[ai * 3 + 2] + tz;

        const float dx = lx - atom_xyz[aj * 3 + 0];
        const float dy = ly - atom_xyz[aj * 3 + 1];
        const float dz = lz - atom_xyz[aj * 3 + 2];
        const float r  = sqrtf(dx*dx + dy*dy + dz*dz + 1e-10f);

        const float rA    = atom_radius[ai];
        const float rB    = atom_radius[aj];
        const float rsum  = rA + rB;
        const float rwa_A = rA + Rw;
        const float surf_A = 4.0f * 3.141592653589793f * rwa_A * rwa_A;
        const float outer_r = rsum + 2.0f * Rw;  // rA + rB + 2·Rw

        // Normalised contact area: linear from 1 at r=rsum to 0 at r=outer_r.
        float rel_area = 0.0f;
        if      (r < rsum)    rel_area = 1.0f;
        else if (r < outer_r) rel_area = 1.0f - (r - rsum) / (outer_r - rsum);

        const bool in_contact = (rel_area > 0.0f);

        // Subtract contact area from this ligand atom's SAS counter.
        if (in_contact && li < MAX_LIG_SAS) {
            atomicAdd(&lig_sas[li], -rel_area * surf_A);
        }

        // COM: energy-matrix interpolation scaled by normalised contact area,
        // optionally distance-weighted (VCT dist-weight term).
        const int   ti   = atom_type[ai];
        const int   tj   = atom_type[aj];
        const float yval = gpu_get_yval(emat_sampled, ti, tj, T, rel_area);
        float com_pair   = yval * rel_area;
        if (do_dw && in_contact) com_pair *= __expf(-r * inv_r0);
        local_com += com_pair;

        // ELEC: distance-dependent-dielectric Coulomb, contacting pairs only.
        if (do_elec && in_contact && r > 0.5f) {
            const float qA = atom_charge[ai];
            const float qB = atom_charge[aj];
            if (qA != 0.0f && qB != 0.0f)
                local_elec += KCOULOMB_F * qA * qB / (P.dielectric * r * r);
        }

        // HBOND: distance-Gaussian × representative angle term (drift model).
        if (do_hbond && in_contact) {
            const int  fa = atom_hflags[ai];
            const int  fb = atom_hflags[aj];
            const bool donor_a    = (fa & 0x1) != 0;
            const bool acceptor_a = (fa & 0x2) != 0;
            const bool donor_b    = (fb & 0x1) != 0;
            const bool acceptor_b = (fb & 0x2) != 0;
            const float qa = atom_charge[ai];
            const float qb = atom_charge[aj];
            const bool a_to_b = donor_a && acceptor_b;
            const bool b_to_a = donor_b && acceptor_a;
            const bool salt   = (qa <= -Q_SALT && qb >= Q_SALT) ||
                                (qb <= -Q_SALT && qa >= Q_SALT);
            if (a_to_b || b_to_a || salt) {
                const float dd = (r - P.hbond_opt_dist) / P.hbond_sigma_dist;
                const float E_dist = __expf(-0.5f * dd * dd);
                const float w = salt ? P.hbond_salt_weight : P.hbond_weight;
                float E = w * E_dist * P.hbond_angle_repr;
                if (E < HBOND_PAIR_MIN) E = HBOND_PAIR_MIN;
                local_hbond += E;
            }
        }

        // WAL: repulsive wall energy when r < perm × (rA+rB).
        const float clash_r = perm * rsum;
        if (r < clash_r) {
            const float r2  = r * r;
            const float r4  = r2 * r2;
            const float r6  = r4 * r2;
            const float inv_r12  = 1.0f / (r6 * r6);
            const float cr2 = clash_r * clash_r;
            const float cr4 = cr2 * cr2;
            const float cr6 = cr4 * cr2;
            const float inv_cr12 = 1.0f / (cr6 * cr6);
            local_wal += KWALL_F * (inv_r12 - inv_cr12);
        }

        // PB_CLASH: PoseBusters vdW-overlap severity (independent of Voronoi).
        if (do_pb) {
            const float cr_pb = P.pb_clash_ratio * (atom_pbvdw[ai] + atom_pbvdw[aj]);
            const float o     = cr_pb - r;
            if (o > 0.0f && r > 1.0e-6f)
                local_pb += powf(o, P.pb_clash_exponent);
        }
    }
    // Pair loop done; ensure all atomicAdds to lig_sas are visible.
    __syncthreads();

    // ── block reductions (single shared scratch, reused per channel) ──────────
    __shared__ float red[NWARPS];
    float com   = block_reduce_sum(local_com,   red, tid);
    float wal   = block_reduce_sum(local_wal,   red, tid);
    float elec  = block_reduce_sum(local_elec,  red, tid);
    float hbond = block_reduce_sum(local_hbond, red, tid);
    float pb    = block_reduce_sum(local_pb,    red, tid);

    // ── SAS + GIST contribution (per ligand atom) ────────────────────────────
    const bool do_gist = (gist_nx > 0) && (gist_data != nullptr);
    float local_sas = 0.0f, local_gist = 0.0f;
    // Centroid accumulation for the pocket penalty (ligand atoms).
    float local_cx = 0.0f, local_cy = 0.0f, local_cz = 0.0f;
    if (n_lig <= MAX_LIG_SAS) {
        for (int la = tid; la < n_lig; la += BLOCK_SIZE) {
            const int   aidx     = lig_first + la;
            const float sas_rem  = fmaxf(0.0f, lig_sas[la]);
            const float rwa_la   = atom_radius[aidx] + Rw;
            const float surf_la  = 4.0f * 3.141592653589793f * rwa_la * rwa_la;
            const float sas_norm = sas_rem / surf_la;
            float contribution;
            if (P.solvent_flat) {
                contribution = P.solventterm * sas_rem;
            } else {
                const int ti_la = atom_type[aidx];
                const float yval_sas = gpu_get_yval(emat_sampled, ti_la, T - 1, T, sas_norm);
                contribution = yval_sas * sas_norm;   // normalised-area branch
            }
            local_sas += P.sas_weight * contribution;

            // Translated ligand coords for GIST + centroid.
            const float lx = atom_xyz[aidx * 3 + 0] + tx;
            const float ly = atom_xyz[aidx * 3 + 1] + ty;
            const float lz = atom_xyz[aidx * 3 + 2] + tz;
            local_cx += lx; local_cy += ly; local_cz += lz;

            if (do_gist) {
                float fx = (lx - gox) * gidx;
                float fy = (ly - goy) * gidy;
                float fz = (lz - goz) * gidz;
                if (fx >= 0.0f && fx < (gist_nx - 1) &&
                    fy >= 0.0f && fy < (gist_ny - 1) &&
                    fz >= 0.0f && fz < (gist_nz - 1)) {
                    int ix = (int)fx, iy = (int)fy, iz = (int)fz;
                    float ddx = fx - ix, ddy = fy - iy, ddz = fz - iz;
                    #define GIDX(i,j,k) (((i)*gist_ny + (j))*gist_nz + (k))
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
                    local_gist += gist_weight * (c0*(1.0f-ddx) + c1*ddx);
                }
            }
        }
    }
    float sas       = block_reduce_sum(local_sas,  red, tid);
    float gist      = block_reduce_sum(local_gist, red, tid);
    float cx_sum    = block_reduce_sum(local_cx,   red, tid);
    float cy_sum    = block_reduce_sum(local_cy,   red, tid);
    float cz_sum    = block_reduce_sum(local_cz,   red, tid);

    // ── pocket-presence penalty (added into pb channel) ──────────────────────
    if (P.pb_pocket_weight > 0.0f && atom_hflags != nullptr && n_lig > 0) {
        const float inv_nl = 1.0f / (float)n_lig;
        const float gx = cx_sum * inv_nl;
        const float gy = cy_sum * inv_nl;
        const float gz = cz_sum * inv_nl;
        // Nearest receptor HEAVY atom to the ligand centroid (brute force).
        float local_min = 3.4e38f;
        for (int a = tid; a < N; a += BLOCK_SIZE) {
            if (a >= lig_first && a <= lig_last) continue;   // skip ligand
            if ((atom_hflags[a] & 0x4) == 0)      continue;  // heavy only
            const float ddx = gx - atom_xyz[a*3+0];
            const float ddy = gy - atom_xyz[a*3+1];
            const float ddz = gz - atom_xyz[a*3+2];
            const float d2 = ddx*ddx + ddy*ddy + ddz*ddz;
            local_min = fminf(local_min, d2);
        }
        float best_d2 = block_reduce_min(local_min, red, tid);
        if (best_d2 < 3.0e38f) {
            const float d = sqrtf(best_d2);
            const float over = d - P.pb_pocket_radius;
            if (over > 0.0f) pb += P.pb_pocket_weight * over * over;
        }
    }

    // ── covalent-constraint term (con), computed by thread 0 over the list ────
    float con = 0.0f;
    if (tid == 0 && n_cons > 0 && cons_i != nullptr) {
        for (int c = 0; c < n_cons; ++c) {
            const int i = cons_i[c];
            const int j = cons_j[c];
            float ix = atom_xyz[i*3+0], iy = atom_xyz[i*3+1], iz = atom_xyz[i*3+2];
            float jx = atom_xyz[j*3+0], jy = atom_xyz[j*3+1], jz = atom_xyz[j*3+2];
            if (i >= lig_first && i <= lig_last) { ix += tx; iy += ty; iz += tz; }
            if (j >= lig_first && j <= lig_last) { jx += tx; jy += ty; jz += tz; }
            const float ddx = ix - jx, ddy = iy - jy, ddz = iz - jz;
            const float d = sqrtf(ddx*ddx + ddy*ddy + ddz*ddz);
            // Baseline KDIST (one ligand-side atom per covalent constraint) minus
            // the Gaussian restraint reward — see .cuh model note.
            con += P.kdist - P.kdist * gpu_gaussian(d, cons_bondlen[c], cons_maxdist[c]);
        }
    }

    // ── write outputs (thread 0) ─────────────────────────────────────────────
    if (tid == 0) {
        double* com_o  = out + 0 * max_pop;
        double* wal_o  = out + 1 * max_pop;
        double* sas_o  = out + 2 * max_pop;
        double* con_o  = out + 3 * max_pop;
        double* elec_o = out + 4 * max_pop;
        double* hb_o   = out + 5 * max_pop;
        double* gist_o = out + 6 * max_pop;
        double* pb_o   = out + 7 * max_pop;
        com_o[chrom_id]  = (double)com;
        wal_o[chrom_id]  = (double)wal;
        sas_o[chrom_id]  = (double)sas;
        con_o[chrom_id]  = (double)con;
        elec_o[chrom_id] = (double)elec;
        hb_o[chrom_id]   = (double)hbond;
        gist_o[chrom_id] = (double)gist;
        pb_o[chrom_id]   = (double)pb;
    }
}

// ─── host API ────────────────────────────────────────────────────────────────

CudaEvalCtx* cuda_eval_init(int   n_atoms,
                             int   n_types,
                             int   max_pop,
                             int   max_genes,
                             int   lig_first,
                             int   lig_last,
                             float perm,
                             const float* h_atom_xyz,
                             const int*   h_atom_type,
                             const float* h_atom_radius,
                             const float* h_emat_sampled)
{
    CudaEvalCtx* ctx = new CudaEvalCtx;
    std::memset(ctx, 0, sizeof(*ctx));   // null all extra/gist pointers by default
    ctx->n_atoms   = n_atoms;
    ctx->n_types   = n_types;
    ctx->max_pop   = max_pop;
    ctx->max_genes = max_genes;
    ctx->lig_first = lig_first;
    ctx->lig_last  = lig_last;
    ctx->perm      = perm;

    const size_t xyz_bytes  = (size_t)n_atoms * 3          * sizeof(float);
    const size_t type_bytes = (size_t)n_atoms               * sizeof(int);
    const size_t rad_bytes  = (size_t)n_atoms               * sizeof(float);
    const size_t em_bytes   = (size_t)n_types * n_types * N_EMAT_SAMPLES * sizeof(float);
    const size_t gene_bytes = (size_t)max_pop * max_genes   * sizeof(double);
    const size_t cf_bytes   = (size_t)max_pop               * sizeof(double);

    // Device allocations
    CUDA_CHECK(cudaMalloc(&ctx->d_atom_xyz,     xyz_bytes));
    CUDA_CHECK(cudaMalloc(&ctx->d_atom_type,    type_bytes));
    CUDA_CHECK(cudaMalloc(&ctx->d_atom_radius,  rad_bytes));
    CUDA_CHECK(cudaMalloc(&ctx->d_emat_sampled, em_bytes));
    CUDA_CHECK(cudaMalloc(&ctx->d_genes,        gene_bytes));
    // Merged output: N_CF_CHANNELS channels contiguous for a single DMA readback
    CUDA_CHECK(cudaMalloc(&ctx->d_cf_out,       N_CF_CHANNELS * cf_bytes));

    // Pinned host memory for async transfers
    CUDA_CHECK(cudaMallocHost(&ctx->h_genes_pinned, gene_bytes));
    CUDA_CHECK(cudaMallocHost(&ctx->h_cf_pinned,    N_CF_CHANNELS * cf_bytes));

    // Create dedicated stream for async operations
    CUDA_CHECK(cudaStreamCreate(&ctx->stream));

    // Upload constant atom data (blocking — only done once at init)
    CUDA_CHECK(cudaMemcpy(ctx->d_atom_xyz,     h_atom_xyz,     xyz_bytes,  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_atom_type,    h_atom_type,    type_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_atom_radius,  h_atom_radius,  rad_bytes,  cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_emat_sampled, h_emat_sampled, em_bytes,   cudaMemcpyHostToDevice));

    return ctx;
}

void cuda_eval_set_extra(CudaEvalCtx* ctx, const GpuCfExtraStatic* extra)
{
    if (!ctx || !extra) return;
    const int   N = ctx->n_atoms;
    const size_t fN = (size_t)N * sizeof(float);
    const size_t iN = (size_t)N * sizeof(int);

    auto up_f = [&](float** dev, const float* host) {
        if (!host) return;
        if (!*dev) CUDA_CHECK(cudaMalloc(dev, fN));
        CUDA_CHECK(cudaMemcpy(*dev, host, fN, cudaMemcpyHostToDevice));
    };
    up_f(&ctx->d_atom_pbvdw,  extra->atom_pbvdw);
    up_f(&ctx->d_atom_charge, extra->atom_charge);
    if (extra->atom_hflags) {
        if (!ctx->d_atom_hflags) CUDA_CHECK(cudaMalloc(&ctx->d_atom_hflags, iN));
        CUDA_CHECK(cudaMemcpy(ctx->d_atom_hflags, extra->atom_hflags, iN, cudaMemcpyHostToDevice));
    }

    // Constraints
    ctx->n_cons = extra->n_cons;
    if (extra->n_cons > 0 && extra->cons_i) {
        const size_t ci = (size_t)extra->n_cons * sizeof(int);
        const size_t cf = (size_t)extra->n_cons * sizeof(float);
        CUDA_CHECK(cudaMalloc(&ctx->d_cons_i,       ci));
        CUDA_CHECK(cudaMalloc(&ctx->d_cons_j,       ci));
        CUDA_CHECK(cudaMalloc(&ctx->d_cons_bondlen, cf));
        CUDA_CHECK(cudaMalloc(&ctx->d_cons_maxdist, cf));
        CUDA_CHECK(cudaMemcpy(ctx->d_cons_i,       extra->cons_i,       ci, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ctx->d_cons_j,       extra->cons_j,       ci, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ctx->d_cons_bondlen, extra->cons_bondlen, cf, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(ctx->d_cons_maxdist, extra->cons_maxdist, cf, cudaMemcpyHostToDevice));
    }

    // GIST grid
    ctx->gist_nx = extra->gist_nx;
    ctx->gist_ny = extra->gist_ny;
    ctx->gist_nz = extra->gist_nz;
    for (int k = 0; k < 3; ++k) {
        ctx->gist_origin[k]    = extra->gist_origin[k];
        ctx->gist_inv_delta[k] = extra->gist_inv_delta[k];
    }
    ctx->gist_weight = extra->gist_weight;
    if (extra->gist_nx > 0 && extra->gist_data) {
        const size_t vox = (size_t)extra->gist_nx * extra->gist_ny * extra->gist_nz * sizeof(float);
        CUDA_CHECK(cudaMalloc(&ctx->d_gist_data, vox));
        CUDA_CHECK(cudaMemcpy(ctx->d_gist_data, extra->gist_data, vox, cudaMemcpyHostToDevice));
    }
}

void cuda_eval_batch(CudaEvalCtx*        ctx,
                     int                 pop_size,
                     int                 n_genes,
                     const double*       h_genes,
                     const GpuCfParams*  params,
                     const GpuCfResults* out)
{
    // Warn if ligand exceeds shared-memory SAS capacity.
    {
        int n_lig = ctx->lig_last - ctx->lig_first + 1;
        if (n_lig > MAX_LIG_SAS) {
            fprintf(stderr, "cuda_eval: n_lig=%d exceeds MAX_LIG_SAS=%d, "
                    "SAS/GIST contribution will be zero for atoms beyond limit\n",
                    n_lig, MAX_LIG_SAS);
        }
    }

    // Validate against allocated buffer sizes — throw on overflow.
    if (pop_size > ctx->max_pop) {
        throw FlexAIDException("cuda_eval: pop_size " + std::to_string(pop_size) +
            " exceeds max_pop " + std::to_string(ctx->max_pop));
    }
    if (n_genes > ctx->max_genes) {
        throw FlexAIDException("cuda_eval: n_genes " + std::to_string(n_genes) +
            " exceeds max_genes " + std::to_string(ctx->max_genes));
    }

    const size_t gene_bytes = (size_t)pop_size * n_genes * sizeof(double);
    const size_t cf_bytes   = (size_t)pop_size            * sizeof(double);

    // Copy genes into pinned buffer, then async upload to device
    memcpy(ctx->h_genes_pinned, h_genes, gene_bytes);
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_genes, ctx->h_genes_pinned, gene_bytes,
                               cudaMemcpyHostToDevice, ctx->stream));

    GpuCfParams P = *params;

    kernel_eval_cf_full<<<pop_size, BLOCK_SIZE, 0, ctx->stream>>>(
        ctx->d_atom_xyz,
        ctx->d_atom_type,
        ctx->d_atom_radius,
        ctx->d_emat_sampled,
        ctx->d_genes,
        ctx->d_atom_pbvdw,
        ctx->d_atom_charge,
        ctx->d_atom_hflags,
        ctx->n_cons,
        ctx->d_cons_i, ctx->d_cons_j, ctx->d_cons_bondlen, ctx->d_cons_maxdist,
        ctx->gist_nx, ctx->gist_ny, ctx->gist_nz,
        ctx->gist_origin[0], ctx->gist_origin[1], ctx->gist_origin[2],
        ctx->gist_inv_delta[0], ctx->gist_inv_delta[1], ctx->gist_inv_delta[2],
        ctx->gist_weight,
        ctx->d_gist_data,
        P,
        ctx->d_cf_out,
        ctx->max_pop,
        ctx->n_atoms,
        ctx->n_types,
        n_genes,
        ctx->lig_first,
        ctx->lig_last,
        ctx->perm);

    CUDA_CHECK(cudaGetLastError());

    // Single merged readback: all channels → pinned host buffer
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_cf_pinned, ctx->d_cf_out,
                               N_CF_CHANNELS * cf_bytes,
                               cudaMemcpyDeviceToHost, ctx->stream));

    // Wait for all stream operations to complete
    CUDA_CHECK(cudaStreamSynchronize(ctx->stream));

    // Scatter merged results to caller's output arrays (per-channel stride = max_pop).
    const double* base = ctx->h_cf_pinned;
    auto ch = [&](int c) { return base + (size_t)c * ctx->max_pop; };
    if (out->com)         memcpy(out->com,         ch(0), cf_bytes);
    if (out->wal)         memcpy(out->wal,         ch(1), cf_bytes);
    if (out->sas)         memcpy(out->sas,         ch(2), cf_bytes);
    if (out->con)         memcpy(out->con,         ch(3), cf_bytes);
    if (out->elec)        memcpy(out->elec,        ch(4), cf_bytes);
    if (out->hbond)       memcpy(out->hbond,       ch(5), cf_bytes);
    if (out->gist_desolv) memcpy(out->gist_desolv, ch(6), cf_bytes);
    if (out->pb_clash)    memcpy(out->pb_clash,    ch(7), cf_bytes);
}

void cuda_eval_shutdown(CudaEvalCtx* ctx)
{
    if (!ctx) return;
    cudaFree(ctx->d_atom_xyz);
    cudaFree(ctx->d_atom_type);
    cudaFree(ctx->d_atom_radius);
    cudaFree(ctx->d_emat_sampled);
    cudaFree(ctx->d_genes);
    cudaFree(ctx->d_cf_out);
    cudaFree(ctx->d_atom_pbvdw);
    cudaFree(ctx->d_atom_charge);
    cudaFree(ctx->d_atom_hflags);
    cudaFree(ctx->d_cons_i);
    cudaFree(ctx->d_cons_j);
    cudaFree(ctx->d_cons_bondlen);
    cudaFree(ctx->d_cons_maxdist);
    cudaFree(ctx->d_gist_data);
    cudaFreeHost(ctx->h_genes_pinned);
    cudaFreeHost(ctx->h_cf_pinned);
    cudaStreamDestroy(ctx->stream);
    delete ctx;
}

#endif  // FLEXAIDS_USE_CUDA
