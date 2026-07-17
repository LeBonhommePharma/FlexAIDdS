// cpu_eval.cpp — CPU fallback for batched chromosome fitness evaluation
//
// Equivalent to cuda_eval.cu / metal_eval.mm but runs on the CPU with
// OpenMP parallelism.  Each chromosome is evaluated independently.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include "cpu_eval.h"
#include "soft_wall.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#  include <omp.h>
#endif

// ─── constants (match GPU kernels exactly) ──────────────────────────────────
static constexpr int   N_EMAT_SAMPLES = 128;
static constexpr float Rw             = 1.4f;   // water probe radius (Å)
static constexpr float PI_F           = 3.141592653589793f;
static constexpr int   MAX_LIG_SAS    = 512;

// ─── context ────────────────────────────────────────────────────────────────
struct CpuEvalCtx {
    std::vector<float> atom_xyz;       // [n_atoms × 3]
    std::vector<int>   atom_type;      // [n_atoms]
    std::vector<float> atom_radius;    // [n_atoms]
    std::vector<float> emat_sampled;   // [n_types × n_types × N_EMAT_SAMPLES]

    int   n_atoms;
    int   n_types;
    int   max_pop;
    int   max_genes;
    int   lig_first;
    int   lig_last;
    float perm;
    float soft_wall_cutoff = 0.40f;
    float k_wal = 50.0f;
};

// ─── helper: interpolated energy-matrix lookup (matches GPU) ────────────────
static inline float cpu_get_yval(const float* emat_sampled,
                                 int t1, int t2, int T, float rel_area)
{
    int   base = (t1 * T + t2) * N_EMAT_SAMPLES;
    rel_area   = std::max(0.0f, std::min(1.0f, rel_area));
    float kf   = rel_area * (N_EMAT_SAMPLES - 1.0f);
    int   k0   = static_cast<int>(kf);
    int   k1   = std::min(k0 + 1, N_EMAT_SAMPLES - 1);
    float frac = kf - static_cast<float>(k0);
    return emat_sampled[base + k0] * (1.0f - frac)
         + emat_sampled[base + k1] * frac;
}

// ─── single-chromosome evaluation ───────────────────────────────────────────
static void eval_one_chromosome(
    const CpuEvalCtx* ctx,
    int chrom_id,
    int n_genes,
    const double* genes,
    double& com_out,
    double& wal_out,
    double& sas_out)
{
    const float* atom_xyz    = ctx->atom_xyz.data();
    const int*   atom_type   = ctx->atom_type.data();
    const float* atom_radius = ctx->atom_radius.data();
    const float* emat        = ctx->emat_sampled.data();
    const int    N           = ctx->n_atoms;
    const int    T           = ctx->n_types;
    const int    lig_first   = ctx->lig_first;
    const int    lig_last    = ctx->lig_last;
    const float  perm        = ctx->perm;

    // Translation from genes (first three genes encode tx, ty, tz)
    const float tx = static_cast<float>(genes[chrom_id * n_genes + 0]);
    const float ty = static_cast<float>(genes[chrom_id * n_genes + 1]);
    const float tz = static_cast<float>(genes[chrom_id * n_genes + 2]);

    const int n_lig = lig_last - lig_first + 1;
    const int n_pro = N - n_lig;

    // Per-ligand SAS accumulator initialised to full surface area
    float lig_sas[MAX_LIG_SAS];
    const int n_lig_capped = std::min(n_lig, MAX_LIG_SAS);
    for (int la = 0; la < n_lig_capped; ++la) {
        float rwa = atom_radius[lig_first + la] + Rw;
        lig_sas[la] = 4.0f * PI_F * rwa * rwa;
    }

    float local_com = 0.0f;
    float local_wal = 0.0f;

    // Pair loop: all ligand-protein atom pairs
    const int n_pairs = n_lig * n_pro;
    for (int pr = 0; pr < n_pairs; ++pr) {
        const int li      = pr / n_pro;
        const int pro_rel = pr % n_pro;
        const int ai      = lig_first + li;
        const int aj      = (pro_rel < lig_first) ? pro_rel : (pro_rel + n_lig);

        // Ligand atom position with translation applied
        const float lx = atom_xyz[ai * 3 + 0] + tx;
        const float ly = atom_xyz[ai * 3 + 1] + ty;
        const float lz = atom_xyz[ai * 3 + 2] + tz;

        const float dx = lx - atom_xyz[aj * 3 + 0];
        const float dy = ly - atom_xyz[aj * 3 + 1];
        const float dz = lz - atom_xyz[aj * 3 + 2];
        const float r  = std::sqrt(dx*dx + dy*dy + dz*dz + 1e-10f);

        const float rA    = atom_radius[ai];
        const float rB    = atom_radius[aj];
        const float rsum  = rA + rB;
        const float rwa_A = rA + Rw;
        const float surf_A = 4.0f * PI_F * rwa_A * rwa_A;
        const float outer_r = rsum + 2.0f * Rw;

        // Normalised contact area: linear from 1 at r=rsum to 0 at r=outer_r
        float rel_area = 0.0f;
        if      (r < rsum)    rel_area = 1.0f;
        else if (r < outer_r) rel_area = 1.0f - (r - rsum) / (outer_r - rsum);

        // Subtract contact area from this ligand atom's SAS counter
        if (rel_area > 0.0f && li < MAX_LIG_SAS) {
            lig_sas[li] -= rel_area * surf_A;
        }

        // COM: energy-matrix interpolation scaled by normalised contact area
        const int   ti   = atom_type[ai];
        const int   tj   = atom_type[aj];
        const float yval = cpu_get_yval(emat, ti, tj, T, rel_area);
        local_com += yval * rel_area;

        // WAL: soft-core k·o² or legacy capped r^-12 (parity with vcfunction)
        const float clash_r = perm * rsum;
        if (r < clash_r) {
            local_wal += soft_wall_fitness_energy_f(
                r, clash_r, ctx->soft_wall_cutoff, ctx->k_wal);
        }
    }

    // SAS contribution: remaining exposed area per ligand atom
    float local_sas = 0.0f;
    if (n_lig <= MAX_LIG_SAS) {
        for (int la = 0; la < n_lig; ++la) {
            const float sas_rem  = std::max(0.0f, lig_sas[la]);
            const float rwa_la   = atom_radius[lig_first + la] + Rw;
            const float surf_la  = 4.0f * PI_F * rwa_la * rwa_la;
            const float sas_norm = sas_rem / surf_la;
            const int   ti_la    = atom_type[lig_first + la];
            const float yval_sas = cpu_get_yval(emat, ti_la, T - 1, T, sas_norm);
            local_sas += yval_sas * sas_norm;
        }
    }

    com_out = static_cast<double>(local_com);
    wal_out = static_cast<double>(local_wal);
    sas_out = static_cast<double>(local_sas);
}

// ─── host API ───────────────────────────────────────────────────────────────

CpuEvalCtx* cpu_eval_init(int   n_atoms,
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
    auto* ctx = new CpuEvalCtx;
    ctx->n_atoms   = n_atoms;
    ctx->n_types   = n_types;
    ctx->max_pop   = max_pop;
    ctx->max_genes = max_genes;
    ctx->lig_first = lig_first;
    ctx->lig_last  = lig_last;
    ctx->perm      = perm;
    ctx->soft_wall_cutoff = 0.40f;
    ctx->k_wal            = 50.0f;

    ctx->atom_xyz.assign(h_atom_xyz, h_atom_xyz + n_atoms * 3);
    ctx->atom_type.assign(h_atom_type, h_atom_type + n_atoms);
    ctx->atom_radius.assign(h_atom_radius, h_atom_radius + n_atoms);

    const int emat_size = n_types * n_types * N_EMAT_SAMPLES;
    ctx->emat_sampled.assign(h_emat_sampled, h_emat_sampled + emat_size);

    return ctx;
}

void cpu_eval_batch(CpuEvalCtx*  ctx,
                    int           pop_size,
                    int           n_genes,
                    const double* h_genes,
                    double*       h_com_out,
                    double*       h_wal_out,
                    double*       h_sas_out)
{
    if (!ctx) return;

    if (pop_size > ctx->max_pop) {
        fprintf(stderr, "cpu_eval: pop_size %d exceeds max_pop %d\n",
                pop_size, ctx->max_pop);
        return;
    }
    if (n_genes > ctx->max_genes) {
        fprintf(stderr, "cpu_eval: n_genes %d exceeds max_genes %d\n",
                n_genes, ctx->max_genes);
        return;
    }

    int n_lig = ctx->lig_last - ctx->lig_first + 1;
    if (n_lig > MAX_LIG_SAS) {
        fprintf(stderr, "cpu_eval: n_lig=%d exceeds MAX_LIG_SAS=%d, "
                "SAS contribution will be zero for atoms beyond limit\n",
                n_lig, MAX_LIG_SAS);
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 4)
#endif
    for (int c = 0; c < pop_size; ++c) {
        eval_one_chromosome(ctx, c, n_genes, h_genes,
                            h_com_out[c], h_wal_out[c], h_sas_out[c]);
    }
}

void cpu_eval_set_soft_wall(CpuEvalCtx* ctx, float soft_wall_cutoff, float k_wal)
{
    if (!ctx) return;
    ctx->soft_wall_cutoff = soft_wall_cutoff;
    ctx->k_wal = (k_wal > 0.0f) ? k_wal : 50.0f;
}

void cpu_eval_shutdown(CpuEvalCtx* ctx)
{
    delete ctx;
}
