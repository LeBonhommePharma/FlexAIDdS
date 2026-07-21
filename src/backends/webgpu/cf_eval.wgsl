// cf_eval.wgsl — WebGPU/WGSL port of the Metal single-complex CF kernel
// (LIB/metal_eval.mm :: kernel_eval_cf_full).
//
// Experimental approximate GPU chromosome evaluation, mirroring the same
// reduced scoring model as the Metal path: only the first three genes are
// treated as a Cartesian translation (tx, ty, tz) applied to the ligand
// atoms; COM (complementarity), WAL (clash) and SAS terms are accumulated
// per ligand-protein atom pair. This does NOT reproduce the full CPU
// ic2cf/vcfunction path (rotations, torsions, Voronoi tessellation) — see
// the parity note in metal_eval.mm. Selected only via --backend webgpu.
//
// Workgroup layout: one workgroup per chromosome, WORKGROUP_SIZE threads.

const N_EMAT_SAMPLES: u32 = 128u;
const WORKGROUP_SIZE: u32 = 256u;
const PI: f32 = 3.14159265358979323846;

struct EvalParams {
    N:         i32,
    T:         i32,
    n_genes:   i32,
    lig_first: i32,
    lig_last:  i32,
    perm:      f32,
    pad0:      i32,
    pad1:      i32,
};

@group(0) @binding(0) var<storage, read>       atom_xyz:     array<f32>;
@group(0) @binding(1) var<storage, read>       atom_type:    array<i32>;
@group(0) @binding(2) var<storage, read>       atom_radius:  array<f32>;
@group(0) @binding(3) var<storage, read>       emat_sampled: array<f32>;
@group(0) @binding(4) var<storage, read>       genes_f:      array<f32>;
@group(0) @binding(5) var<storage, read_write>  cf_com_out:   array<f32>;
@group(0) @binding(6) var<storage, read_write>  cf_wal_out:   array<f32>;
@group(0) @binding(7) var<storage, read_write>  cf_sas_out:   array<f32>;
@group(0) @binding(8) var<uniform>             p:            EvalParams;

// Bitcast-CAS float add into a u32 storage atomic (WGSL has no native
// atomic<f32>), same trick as the Metal tg_atomic_sub_float CAS loop.
var<workgroup> lig_sas_bits: array<atomic<u32>, 256u>;
var<workgroup> reduce_com: array<f32, WORKGROUP_SIZE>;
var<workgroup> reduce_wal: array<f32, WORKGROUP_SIZE>;

fn atomic_add_f32(idx: u32, val: f32) {
    loop {
        let old_bits = atomicLoad(&lig_sas_bits[idx]);
        let new_val  = bitcast<f32>(old_bits) + val;
        let new_bits = bitcast<u32>(new_val);
        let r = atomicCompareExchangeWeak(&lig_sas_bits[idx], old_bits, new_bits);
        if (r.exchanged) { break; }
    }
}

fn gpu_get_yval(t1: i32, t2: i32, T: i32, rel_area_in: f32) -> f32 {
    let base = (t1 * T + t2) * i32(N_EMAT_SAMPLES);
    let rel_area = clamp(rel_area_in, 0.0, 1.0);
    let kf = rel_area * f32(N_EMAT_SAMPLES - 1u);
    let k0 = i32(kf);
    let k1 = min(k0 + 1, i32(N_EMAT_SAMPLES) - 1);
    let frac = kf - f32(k0);
    return emat_sampled[base + k0] * (1.0 - frac) + emat_sampled[base + k1] * frac;
}

@compute @workgroup_size(WORKGROUP_SIZE)
fn kernel_eval_cf_full(
    @builtin(local_invocation_id) local_id: vec3<u32>,
    @builtin(workgroup_id) wg_id: vec3<u32>)
{
    let tid      = local_id.x;
    let chrom_id = wg_id.x;

    let n_lig   = p.lig_last - p.lig_first + 1;
    let n_pro   = p.N - n_lig;
    let n_pairs = n_lig * n_pro;

    let gbase = i32(chrom_id) * p.n_genes;
    let tx = genes_f[gbase + 0];
    let ty = genes_f[gbase + 1];
    let tz = genes_f[gbase + 2];

    var la: i32 = i32(tid);
    while (la < n_lig && la < 256) {
        let ra  = atom_radius[p.lig_first + la];
        let rwa = ra + 1.4;
        atomicStore(&lig_sas_bits[u32(la)], bitcast<u32>(4.0 * PI * rwa * rwa));
        la += i32(WORKGROUP_SIZE);
    }
    workgroupBarrier();

    var local_com: f32 = 0.0;
    var local_wal: f32 = 0.0;

    var pr: i32 = i32(tid);
    while (pr < n_pairs) {
        let li      = pr / n_pro;
        let pro_rel = pr % n_pro;
        let ai      = p.lig_first + li;
        var aj: i32;
        if (pro_rel < p.lig_first) { aj = pro_rel; } else { aj = pro_rel + n_lig; }

        let lx = atom_xyz[ai * 3 + 0] + tx;
        let ly = atom_xyz[ai * 3 + 1] + ty;
        let lz = atom_xyz[ai * 3 + 2] + tz;
        let dx = lx - atom_xyz[aj * 3 + 0];
        let dy = ly - atom_xyz[aj * 3 + 1];
        let dz = lz - atom_xyz[aj * 3 + 2];
        let r  = sqrt(dx * dx + dy * dy + dz * dz + 1e-10);

        let rA = atom_radius[ai];
        let rB = atom_radius[aj];
        let rsum = rA + rB;
        let rwa_A = rA + 1.4;
        let surf_A = 4.0 * PI * rwa_A * rwa_A;
        let outer_r = rsum + 2.8;

        var rel_area: f32 = 0.0;
        if (r < rsum) {
            rel_area = 1.0;
        } else if (r < outer_r) {
            rel_area = 1.0 - (r - rsum) / (outer_r - rsum);
        }

        if (rel_area > 0.0 && li < 256) {
            atomic_add_f32(u32(li), -rel_area * surf_A);
        }

        let ti = atom_type[ai];
        let tj = atom_type[aj];
        let yval = gpu_get_yval(ti, tj, p.T, rel_area);
        local_com += yval * rel_area;

        let clash_r = p.perm * rsum;
        if (r < clash_r && r > 0.0) {
            let inv_r12  = 1.0 / pow(r, 12.0);
            let inv_cr12 = 1.0 / pow(clash_r, 12.0);
            local_wal += 1.0e6 * (inv_r12 - inv_cr12);
        }

        pr += i32(WORKGROUP_SIZE);
    }

    reduce_com[tid] = local_com;
    reduce_wal[tid] = local_wal;
    workgroupBarrier();

    // Tree reduction across the workgroup.
    var stride: u32 = WORKGROUP_SIZE / 2u;
    loop {
        if (stride == 0u) { break; }
        if (tid < stride) {
            reduce_com[tid] += reduce_com[tid + stride];
            reduce_wal[tid] += reduce_wal[tid + stride];
        }
        workgroupBarrier();
        stride = stride / 2u;
    }
    if (tid == 0u) {
        cf_com_out[chrom_id] = reduce_com[0];
        cf_wal_out[chrom_id] = reduce_wal[0];
    }
    workgroupBarrier();

    var local_sas: f32 = 0.0;
    la = i32(tid);
    while (la < n_lig && la < 256) {
        let sas_rem  = max(0.0, bitcast<f32>(atomicLoad(&lig_sas_bits[u32(la)])));
        let rwa_la   = atom_radius[p.lig_first + la] + 1.4;
        let surf_la  = 4.0 * PI * rwa_la * rwa_la;
        let sas_norm = sas_rem / surf_la;
        let ti_la    = atom_type[p.lig_first + la];
        let yval_sas = gpu_get_yval(ti_la, p.T - 1, p.T, sas_norm);
        local_sas += yval_sas * sas_norm;
        la += i32(WORKGROUP_SIZE);
    }
    reduce_com[tid] = local_sas;
    workgroupBarrier();
    stride = WORKGROUP_SIZE / 2u;
    loop {
        if (stride == 0u) { break; }
        if (tid < stride) {
            reduce_com[tid] += reduce_com[tid + stride];
        }
        workgroupBarrier();
        stride = stride / 2u;
    }

    if (tid == 0u) {
        cf_sas_out[chrom_id] = reduce_com[0];
    }
}
