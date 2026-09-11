#include "gaboom.h"
#include "fileio.h"
#include "LigandRingFlex/LigandRingFlex.h"   // Phase 2: ring pucker apply
#include "VibEntropy.h"
#include "tENCoM/tencm.h"
#include "encom.h"      // encom::ENCoMEngine::compute_vibrational_entropy (classical HO S_vib)

#include <mutex>
#include <array>
#include <algorithm>   // std::find / std::min / std::max for the flexbond DOF set

#ifdef _OPENMP
#  include <omp.h>
#endif

// A3 perf: file-scope struct definitions so thread_local scratch can hold them.
struct IC2CFSavedAtom      { int idx; atom value; };
struct IC2CFSavedResRot    { int idx; int  rot;   };

// Per-thread reusable buffers — eliminates per-eval heap allocs under OMP.
struct IC2CFScratch {
	std::vector<IC2CFSavedAtom>           saved_atoms;      // reserve 64
	std::vector<IC2CFSavedResRot>         saved_res_rots;   // reserve 16
	std::vector<char>                     saved_res_seen;   // sized to res_cnt+1
	std::vector<std::pair<int,int>>       intraclashes;
	// FLEXAIDDS_DSVIB per-pose receptor field gather. Per-thread so the OMP
	// scoring path allocates once, not once per pose; untouched when the gate
	// is off (the term returns before reaching them).
	std::vector<std::array<float,3>>      dsvib_field;
	std::vector<char>                     dsvib_picked;
};
// One IC2CFScratch per OMP thread (no synchronisation needed — each thread
// has exclusive access to its own scratch throughout ic2cf execution).
static thread_local IC2CFScratch tl_ic2cf_scratch;

// Per-pose tENCoM vibrational Shannon entropy H(ω) in nats (single-rep H_pop).
// Skipped when tencom_weight == 0 so existing benchmarks pay no ANM cost.
static double compute_ligand_h_rep(const FA_Global* FA, const atom* atoms) {
	if (!FA || FA->tencom_weight <= 0.0f) return 0.0;
	const int lig_start = (FA->resligand && FA->resligand->fatm)
	                      ? FA->resligand->fatm[0] : -1;
	const int lig_end_incl = (FA->resligand && FA->resligand->latm)
	                      ? FA->resligand->latm[0] : -1;
	if (lig_start < 0 || lig_end_incl < lig_start) return 0.0;

	tencm::TorsionalENM lig_enm;
	lig_enm.build_from_ligand(atoms, lig_start, lig_end_incl + 1);
	if (!lig_enm.is_built()) return 0.0;

	// Rigid-body modes REMOVED, not merely sign-filtered. The former
	// `nm.eigenvalue > 0.0` kept whichever of the six landed positive from
	// round-off, which set the lower log-frequency bin edge and moved H_pop by
	// ~56% (BU72: 1.896 -> 2.950) with a survivor count unstable under 1e-6 A
	// jitter. This value feeds tencom_weight * cf.h_rep, so it must be
	// reproducible.
	const std::vector<double> eigs = lig_enm.vibrational_eigenvalues();
	if (eigs.empty()) return 0.0;

	const std::vector<std::vector<double>> single = { eigs };
	return vibentropy::compute_vib_entropy_collapse(single).H_pop;
}

// ═══ FLEXAIDDS_DSVIB — thermodynamic ligand -T*dS_vib inside the objective ═══
//
// WHAT IS DIFFERENT FROM compute_ligand_h_rep ABOVE, since both sit in the same
// file and both call tENCoM:
//
//   h_rep        32-bin SHANNON entropy of the log-frequency spectrum of a
//                CARTESIAN 3N ligand ANM. Dimensionless (nats). Enters the CF as
//                tencom_weight * h_rep, i.e. only through a FITTED weight, which
//                makes it a heuristic ranking feature. Undefined over an empty
//                spectrum and identically 0 over a single mode, so it cannot be
//                evaluated at all on the 4 zero-DOF and 5 one-DOF ligands of the
//                Astex roster.
//   this term    classical-harmonic S_vib over the ligand's TORSIONAL (dihedral)
//                DOFs -- a sum over modes, so 0 modes gives exactly 0 and a rigid
//                ligand is correct rather than undefined. In kcal/mol/K, so
//                -T*dS_vib is in kcal/mol and enters dG with NO fitted weight.
//
// THE CYCLE, AND WHY ONLY TWO TERMS ARE COMPUTED. dS_vib(bind) =
// S(complex) - S(apo) - S(ligand). No builder in tencm.cpp assembles a joint
// protein+ligand Hessian, so S(complex) is not directly computable. It IS
// computable in the rigid-receptor limit these arms run in (autoflex_max = 0):
// with every receptor internal coordinate frozen, the receptor contributes no
// mode to either the complex or the apo spectrum, so S(complex) - S(apo) is
// EXACTLY the ligand's torsional entropy in the static receptor field. The apo
// term is not dropped, it is identically absent; the emitted columns say so.
//
// MODE-COUNT CONSERVATION IS CHECKED, NOT ASSUMED. The absolute S_vib of
// encom.cpp carries an UNCALIBRATED eigenvalue->omega scale, so each term alone
// is model-scale only. In the difference the whole calibration-bearing part is
// M*kB*(1 + ln(kB T / (hbar c))) and cancels -- but only if both states have the
// same M surviving modes. vibrational_eigenvalues() can drop a non-finite or
// non-positive eigenvalue in one state and not the other, which would leave a
// calibration-dependent residue masquerading as physics. A mismatch is refused.
static void compute_ligand_dsvib(const FA_Global* FA, const atom* atoms, cfstr* cf)
{
	cf->dsvib_status = 0;
	if (!FA || FA->dsvib_mode <= 0) return;

	const int lig_start    = (FA->resligand && FA->resligand->fatm)
	                         ? FA->resligand->fatm[0] : -1;
	const int lig_end_incl = (FA->resligand && FA->resligand->latm)
	                         ? FA->resligand->latm[0] : -1;
	if (lig_start < 0 || lig_end_incl < lig_start) { cf->dsvib_status = 4; return; }

	// Temperature from the RUN, never a literal. 0 means the configuration never
	// established one, and a thermodynamic term with an invented temperature is
	// worse than no term at all.
	const double T_K = static_cast<double>(FA->dsvib_T_K);
	if (!(T_K > 0.0)) { cf->dsvib_status = 5; return; }

	// ── fixed receptor field nodes, built ONCE ──────────────────────────────
	// The receptor is rigid for the whole run (autoflex_max = 0), so its heavy
	// atoms and the cell grid over them are run-invariant. Rebuilding per pose
	// would dominate the cost; caching makes the per-pose work proportional to
	// the LIGAND size, not the receptor's.
	struct FieldGrid {
		std::vector<std::array<float,3>> pts;
		std::vector<int>   starts;      // CSR offsets, size nx*ny*nz + 1
		std::vector<int>   items;
		float lo[3]{0.f,0.f,0.f};
		float cell = 9.0f;
		int   nx = 0, ny = 0, nz = 0;
		int   n_receptor_heavy = 0;
	};
	static FieldGrid   grid;
	static std::once_flag grid_once;
	std::call_once(grid_once, [&]() {
		grid.cell = tencm::DEFAULT_RC;     // one rule: the ENM contact cutoff
		auto is_h = [](const atom& a) noexcept {
			const char* e = a.element;
			while (*e == ' ') ++e;
			return e[0] == 'H' && (e[1] == '\0' || e[1] == ' ');
		};
		const int n_end = (FA->atm_cnt_real > 0) ? FA->atm_cnt_real : FA->atm_cnt;
		for (int ai = 1; ai <= n_end; ++ai) {
			if (ai >= lig_start && ai <= lig_end_incl) continue;   // ligand is not its own field
			if (is_h(atoms[ai])) continue;
			grid.pts.push_back({ atoms[ai].coor[0], atoms[ai].coor[1], atoms[ai].coor[2] });
		}
		grid.n_receptor_heavy = static_cast<int>(grid.pts.size());
		if (grid.pts.empty()) return;

		float hi[3];
		for (int d = 0; d < 3; ++d) { grid.lo[d] = grid.pts[0][d]; hi[d] = grid.pts[0][d]; }
		for (const auto& p : grid.pts)
			for (int d = 0; d < 3; ++d) {
				if (p[d] < grid.lo[d]) grid.lo[d] = p[d];
				if (p[d] > hi[d])      hi[d]      = p[d];
			}
		grid.nx = std::max(1, static_cast<int>((hi[0] - grid.lo[0]) / grid.cell) + 1);
		grid.ny = std::max(1, static_cast<int>((hi[1] - grid.lo[1]) / grid.cell) + 1);
		grid.nz = std::max(1, static_cast<int>((hi[2] - grid.lo[2]) / grid.cell) + 1);

		const std::size_t ncell = static_cast<std::size_t>(grid.nx) * grid.ny * grid.nz;
		std::vector<int> counts(ncell + 1, 0);
		auto cell_of = [&](const std::array<float,3>& p) {
			int ix = static_cast<int>((p[0] - grid.lo[0]) / grid.cell);
			int iy = static_cast<int>((p[1] - grid.lo[1]) / grid.cell);
			int iz = static_cast<int>((p[2] - grid.lo[2]) / grid.cell);
			ix = std::min(std::max(ix, 0), grid.nx - 1);
			iy = std::min(std::max(iy, 0), grid.ny - 1);
			iz = std::min(std::max(iz, 0), grid.nz - 1);
			return (static_cast<std::size_t>(iz) * grid.ny + iy) * grid.nx + ix;
		};
		for (const auto& p : grid.pts) ++counts[cell_of(p) + 1];
		for (std::size_t c = 1; c <= ncell; ++c) counts[c] += counts[c - 1];
		grid.starts = counts;
		grid.items.assign(grid.pts.size(), 0);
		std::vector<int> fill(ncell, 0);
		for (std::size_t i = 0; i < grid.pts.size(); ++i) {
			const std::size_t c = cell_of(grid.pts[i]);
			grid.items[static_cast<std::size_t>(grid.starts[c]) + fill[c]] = static_cast<int>(i);
			++fill[c];
		}
		fprintf(stderr,
		        "[DSVIB-INIT] receptor_heavy_nodes=%d grid=%dx%dx%d cell=%.2fA "
		        "cutoff=%.2fA T=%.2fK basis=torsional units=kcal/mol/K\n",
		        grid.n_receptor_heavy, grid.nx, grid.ny, grid.nz,
		        grid.cell, static_cast<double>(tencm::DEFAULT_RC), T_K);
	});
	if (grid.pts.empty()) { cf->dsvib_status = 4; return; }

	// ── per-pose: gather receptor heavy atoms within cutoff of the ligand ───
	const float rc = tencm::DEFAULT_RC;
	std::vector<std::array<float,3>>& field = tl_ic2cf_scratch.dsvib_field;
	std::vector<char>& picked = tl_ic2cf_scratch.dsvib_picked;
	field.clear();
	picked.assign(grid.pts.size(), 0);
	auto is_h2 = [](const atom& a) noexcept {
		const char* e = a.element;
		while (*e == ' ') ++e;
		return e[0] == 'H' && (e[1] == '\0' || e[1] == ' ');
	};
	const float rc2 = rc * rc;
	for (int ai = lig_start; ai <= lig_end_incl; ++ai) {
		if (is_h2(atoms[ai])) continue;
		const float* L = atoms[ai].coor;
		int ix = static_cast<int>((L[0] - grid.lo[0]) / grid.cell);
		int iy = static_cast<int>((L[1] - grid.lo[1]) / grid.cell);
		int iz = static_cast<int>((L[2] - grid.lo[2]) / grid.cell);
		for (int dz = -1; dz <= 1; ++dz)
		for (int dy = -1; dy <= 1; ++dy)
		for (int dx = -1; dx <= 1; ++dx) {
			const int jx = ix + dx, jy = iy + dy, jz = iz + dz;
			if (jx < 0 || jy < 0 || jz < 0 || jx >= grid.nx || jy >= grid.ny || jz >= grid.nz)
				continue;
			const std::size_t c =
				(static_cast<std::size_t>(jz) * grid.ny + jy) * grid.nx + jx;
			for (int t = grid.starts[c]; t < grid.starts[c + 1]; ++t) {
				const int pi = grid.items[static_cast<std::size_t>(t)];
				if (picked[static_cast<std::size_t>(pi)]) continue;
				const auto& P = grid.pts[static_cast<std::size_t>(pi)];
				const float ddx = P[0] - L[0], ddy = P[1] - L[1], ddz = P[2] - L[2];
				if (ddx*ddx + ddy*ddy + ddz*ddz > rc2) continue;
				picked[static_cast<std::size_t>(pi)] = 1;
				field.push_back(P);
			}
		}
	}

	// ── DOF SET: THE ENGINE'S OWN FLEXIBLE BONDS, NOT A RE-PERCEPTION ────────
	// tencm's atom[] overloads perceive the DOF set with a PROXY rule: every
	// distinct rec[1]--rec[0] pair carried by a recs=='m' ligand atom. recs=='m'
	// flags every reconstruction-frame atom, not only the atoms whose dihedral is
	// actually a search variable, so the proxy OVER-COUNTS. Measured on 1JD0:
	// the proxy gives 3 where the engine's own FA->nflexbonds is 1. An entropy
	// computed over 3 torsions when the search moves 1 is not the quantity this
	// term claims to be, and the difference is not a rounding error -- it is a
	// different coordinate set.
	//
	// The engine's authoritative DOF set is the contiguous flexbond gene range
	// map_par[fb0 .. fb0+nflexbonds-1]: only add2_optimiz_vec.cpp:242 increments
	// nflexbonds, and it is the one site that sets a real .bnd and uses
	// delta_flexible. That is the same delimiter gaboom.cpp and the ligtopo
	// sidecar (top.cpp) already use. Building the topology from it makes
	// n_torsion_dofs() == FA->nflexbonds true BY CONSTRUCTION rather than by
	// agreement, so there is ONE perception in this process instead of two that
	// have to be cross-checked afterwards.
	const int fb0 = FA->map_par_flexbond_first_index;
	const int nfb = (fb0 >= 0) ? FA->nflexbonds : 0;
	cf->dsvib_nflex = nfb;

	auto is_h3 = [](const atom& a) noexcept {
		const char* e = a.element;
		while (*e == ' ') ++e;
		return e[0] == 'H' && (e[1] == '\0' || e[1] == ' ');
	};

	// Heavy-atom node array over the ligand span. Same rule and same ordering as
	// tencm's own perception, so coordinates and topology index the same nodes.
	const std::size_t lspan = static_cast<std::size_t>(lig_end_incl + 1 - lig_start);
	std::vector<std::array<float,3>> lxyz;
	std::vector<int> node_of(lspan, -1);
	lxyz.reserve(lspan);
	for (int ai = lig_start; ai <= lig_end_incl; ++ai) {
		if (is_h3(atoms[ai])) continue;
		node_of[static_cast<std::size_t>(ai - lig_start)] = static_cast<int>(lxyz.size());
		lxyz.push_back({ atoms[ai].coor[0], atoms[ai].coor[1], atoms[ai].coor[2] });
	}
	const int Na = static_cast<int>(lxyz.size());
	if (Na < 3) { cf->dsvib_status = 4; return; }

	tencm::TorsionalENM::LigandTopology topo;
	topo.adjacency.assign(static_cast<std::size_t>(Na), {});
	for (int ai = lig_start; ai <= lig_end_incl; ++ai) {
		const int na = node_of[static_cast<std::size_t>(ai - lig_start)];
		if (na < 0) continue;
		for (int t = 1; t <= atoms[ai].bond[0] && t < 7; ++t) {
			const int nb = atoms[ai].bond[t];
			if (nb < lig_start || nb > lig_end_incl) continue;
			const int nn = node_of[static_cast<std::size_t>(nb - lig_start)];
			if (nn >= 0) topo.adjacency[static_cast<std::size_t>(na)].push_back(nn);
		}
	}

	// Rotatable bonds straight off the engine's gene range. The typ/bnd test is
	// the ligtopo sidecar's own secondary discriminator: a real flexible bond
	// carries typ==2 and bnd>=0, while a rigid-body placement pseudo-dihedral
	// carries bnd=-1. If either test fires the contiguous-range assumption has
	// broken, and the honest response is to refuse rather than assemble a DOF set
	// that is neither the engine's nor the proxy's.
	int n_range_bad = 0, n_unmapped = 0;
	for (int p = fb0; fb0 >= 0 && p < fb0 + nfb; ++p) {
		if (FA->map_par[p].typ != 2 || FA->map_par[p].bnd < 0) { ++n_range_bad; continue; }
		const int a  = FA->map_par[p].atm;
		const int b0 = atoms[a].rec[0];
		const int b1 = atoms[a].rec[1];
		if (b0 < lig_start || b0 > lig_end_incl ||
		    b1 < lig_start || b1 > lig_end_incl) { ++n_unmapped; continue; }
		const int n0 = node_of[static_cast<std::size_t>(b0 - lig_start)];
		const int n1 = node_of[static_cast<std::size_t>(b1 - lig_start)];
		if (n0 < 0 || n1 < 0) { ++n_unmapped; continue; }
		const std::pair<int,int> key{ std::min(n0,n1), std::max(n0,n1) };
		if (std::find(topo.rot_bonds.begin(), topo.rot_bonds.end(), key)
		    == topo.rot_bonds.end()) topo.rot_bonds.push_back(key);
	}
	if (n_range_bad > 0 || n_unmapped > 0) {
		cf->dsvib_ndof   = static_cast<int>(topo.rot_bonds.size());
		cf->dsvib_status = 6;
		return;
	}

	// ── the two spectra: SAME DOFs, SAME spring, SAME cutoff, one difference ──
	// Explicit-topology overloads, so neither state re-perceives anything.
	tencm::TorsionalENM enm_free, enm_field;
	enm_free.build_from_ligand_torsional(lxyz.data(), Na, topo, rc);
	enm_field.build_from_ligand_torsional_in_field(
		lxyz.data(), Na, topo,
		field.empty() ? nullptr : field.data(),
		static_cast<int>(field.size()), rc);

	cf->dsvib_ndof = enm_free.n_torsion_dofs();

	// The by-construction equality, asserted in the output instead of a comment.
	// Both states must also agree with each other: they are handed the same topo
	// object, so a difference here would mean the builder mutated it.
	if (cf->dsvib_ndof != nfb ||
	    enm_field.n_torsion_dofs() != cf->dsvib_ndof) {
		cf->dsvib_status = 7;
		return;
	}

	// Fewer than 2 rotatable bonds: the builder refuses (the 32-bin Shannon
	// functional is undefined there), but the THERMODYNAMIC answer is not
	// undefined -- a rigid ligand has no internal torsional entropy to lose, so
	// S_vib = 0 exactly in both states and dS_vib = 0 exactly. Report that as a
	// computed zero with its own status, not as a failure and not as a silent 0.
	if (!enm_free.is_built() || !enm_field.is_built()) {
		if (cf->dsvib_ndof < 2) {
			cf->dsvib_s_free = 0.0; cf->dsvib_s_field = 0.0;
			cf->dsvib_ds = 0.0;     cf->minus_T_dsvib = 0.0;
			cf->dsvib_status = 2;
		} else {
			cf->dsvib_status = 4;
		}
		return;
	}

	const std::vector<double> ev_free  = enm_free.vibrational_eigenvalues();
	const std::vector<double> ev_field = enm_field.vibrational_eigenvalues();
	if (ev_free.size() != ev_field.size() || ev_free.empty()) {
		// The calibration prefactor only cancels when both spectra carry the same
		// number of modes. It does not here, so the difference would be
		// contaminated by the uncalibrated scale. Refuse; do not emit a number.
		cf->dsvib_status = 3;
		return;
	}

	auto to_modes = [](const std::vector<double>& ev) {
		std::vector<encom::NormalMode> m(ev.size());
		for (std::size_t i = 0; i < ev.size(); ++i) {
			m[i].index      = static_cast<int>(i) + 1;
			m[i].eigenvalue = ev[i];
			m[i].frequency  = std::sqrt(ev[i]);
		}
		return m;
	};
	// eigenvalue_cutoff = 0.0, NOT the 1e-6 default: a second, different filter
	// here could drop a mode in one state and not the other and re-break the
	// mode-count conservation that was just checked.
	const encom::VibrationalEntropy S_free =
		encom::ENCoMEngine::compute_vibrational_entropy(to_modes(ev_free),  T_K, 0.0);
	const encom::VibrationalEntropy S_field =
		encom::ENCoMEngine::compute_vibrational_entropy(to_modes(ev_field), T_K, 0.0);
	if (S_free.n_modes != S_field.n_modes) { cf->dsvib_status = 3; return; }

	cf->dsvib_s_free  = S_free.S_vib_kcal_mol_K;
	cf->dsvib_s_field = S_field.S_vib_kcal_mol_K;
	cf->dsvib_ds      = cf->dsvib_s_field - cf->dsvib_s_free;
	cf->minus_T_dsvib = -T_K * cf->dsvib_ds;
	if (!std::isfinite(cf->minus_T_dsvib)) {
		cf->dsvib_s_free = 0.0; cf->dsvib_s_field = 0.0;
		cf->dsvib_ds = 0.0;     cf->minus_T_dsvib = 0.0;
		cf->dsvib_status = 4;
		return;
	}
	cf->dsvib_status = 1;
}

/******************************************************************************
 * SUBROUTINE ic2cf gets a vector with internal coordinates rebuilds the 
 * cartesian coordinates and calculates the complementarity function. Its 
 * input vector has the list of ic's that are to be optimized and a global
 * vector contains the information of what kind of variable each item in icv
 * is and to which residue it belongs.
 *****************************************************************************/

//THE PROCEDURE SHOULD RECEIVE A 2ND SET OF GENES THAT ENCODES FOR THE ROTAMER DISTRIBUTION IN THE BPK

cfstr ic2cf(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue,
			gridpoint* cleftgrid,int npar, double* icv)
{
  
	// static int nbranch = 0;
	
	int i,j,k;
	int cat;    /* atom number constrained to the one considered */

	// Value-initialize every CF field (never leave elec/gist_desolv/etc. garbage).
	cfstr cf{};
	
	int rclash=0;

	psFlexDEE_Node psFlexDEENode;
	int dee_val;
	int rotflag;
  
	unsigned int grd_idx;
	unsigned int rot_idx;
	//float threshold=10.0;
	//int nflxchk=-1;
	int normalmode=-1;
	int deelig_list[100];
	
	//int rigid_clash=0;
	//float rand=0.0;
	//float min_dis=10.0;
  
	// copy values from icv into respective srtructure atom ic fields 
	// andcompute the ic of a constrained atom prior to reconstruction
	
	//for(i=0;i<npar;i++){printf("[%8.3f]",icv[i]);}printf("\n");
	//PAUSE;

	//printf("NEW INDIVIDUAL=");
	for(i=0;i<npar;i++){
		//printf("[%8.3f]",icv[i]);

		if(FA->map_par[i].typ==-1) { //by index
			
			grd_idx = (uint)icv[i];
			atoms[FA->map_par[i].atm].dis = cleftgrid[grd_idx].dis;
			atoms[FA->map_par[i].atm].ang = cleftgrid[grd_idx].ang;
			atoms[FA->map_par[i].atm].dih = cleftgrid[grd_idx].dih;
			
		}else if(FA->map_par[i].typ==0)  {
			atoms[FA->map_par[i].atm].dis = (float)icv[i];
			
		}else if(FA->map_par[i].typ==1)  {
			atoms[FA->map_par[i].atm].ang = (float)icv[i];
			
		}else if(FA->map_par[i].typ==2)  {
			atoms[FA->map_par[i].atm].dih = (float)icv[i];
			
			j=FA->map_par[i].atm;
			cat=atoms[j].rec[3];
			if(cat != 0){
				while(cat != FA->map_par[i].atm){
					atoms[cat].dih=atoms[j].dih + atoms[cat].shift; 
					j=cat;
					cat=atoms[j].rec[3];
				}
			}
			
		}else if(FA->map_par[i].typ==3) { //by index
			grd_idx = (uint)icv[i];
			//printf("icv(index): %d\n", grd_idx);
			//PAUSE;
      
			// serves as flag , but also as grid index
			normalmode=grd_idx;
      
		}else if(FA->map_par[i].typ==4)  {
			// WAS: rot_idx = (uint)(icv[i]+0.5f);  -- (uint) of a NEGATIVE float is
			// undefined behaviour and yielded a huge index -> SIGSEGV on fatm[].
			rot_idx = (uint)rot_gene_index(icv[i],
			          &residue[atoms[FA->map_par[i].atm].ofres], "ic2cf");
      
			residue[atoms[FA->map_par[i].atm].ofres].rot=(int)rot_idx;
      
			/*
			  printf("residue[%d].rot[%d] - fatm=%d - latm=%d\n",
			  residue[atoms[FA->map_par[i].atm].ofres].number,
			  residue[atoms[FA->map_par[i].atm].ofres].rot,
			  residue[atoms[FA->map_par[i].atm].ofres].fatm[rot_idx],
			  residue[atoms[FA->map_par[i].atm].ofres].latm[rot_idx]);
			*/
      
		}
    
	}
	//printf("HERE\n");
	//PAUSE;
  
	// do not alter default (ini) protein conf.
	if(normalmode > -1){
		alter_mode(atoms,residue,FA->normal_grid[normalmode],FA->res_cnt,FA->normal_modes);
	}

	// Save FA->ori and the mutable atom/residue state before any modification.
	// ic2cf() is called repeatedly on shared atom/residue buffers in serial
	// paths, so every evaluation must leave the caller's baseline unchanged.
	float ori_save[3] = {FA->ori[0], FA->ori[1], FA->ori[2]};
	// A3: use thread-local pre-allocated scratch (no heap alloc per eval)
	IC2CFScratch& scr = tl_ic2cf_scratch;
	scr.saved_atoms.clear();
	scr.saved_res_rots.clear();
	if(scr.saved_atoms.capacity() < 64)    scr.saved_atoms.reserve(64);
	if(scr.saved_res_rots.capacity() < 16) scr.saved_res_rots.reserve(16);
	// ensure saved_res_seen is large enough and zeroed for this call
	const int res_cnt_p1 = FA->res_cnt + 1;
	if((int)scr.saved_res_seen.size() < res_cnt_p1)
		scr.saved_res_seen.assign(res_cnt_p1, 0);
	else
		std::fill(scr.saved_res_seen.begin(), scr.saved_res_seen.begin() + res_cnt_p1, 0);
	std::vector<IC2CFSavedAtom>&   saved_atoms    = scr.saved_atoms;
	std::vector<IC2CFSavedResRot>& saved_res_rots = scr.saved_res_rots;
	std::vector<char>&             saved_res_seen  = scr.saved_res_seen;
	for (int r = 0; r < FA->nors; ++r)
		for (int m = 0; m < FA->nmov[r]; ++m) {
			int ai = FA->mov[r][m];
			saved_atoms.push_back(IC2CFSavedAtom{ai, atoms[ai]});
		}
	for (i = 0; i < npar; ++i) {
		if (FA->map_par[i].typ == 4) {
			int ri = atoms[FA->map_par[i].atm].ofres;
			if (ri >= 0 && ri <= FA->res_cnt && !saved_res_seen[ri]) {
				saved_res_rots.push_back(IC2CFSavedResRot{ri, residue[ri].rot});
				saved_res_seen[ri] = 1;
			}
		}
	}

	// ── Ring pucker apply (LigandRingFlex Phase 2) ───────────────────────────
	// Snap each furanose ring's internal dihedrals (.dih) to the Cremer-Pople
	// pucker phase carried by the current chromosome (loaded into
	// FA->ring_cur_phases before this call). Must run AFTER the icv→.dih copy
	// loop above and BEFORE buildcc(), so the reconstructed Cartesian coords
	// reflect the pucker. Gated OFF by default; ring bonds are excluded from
	// map_par, so this never perturbs the standard torsional genes.
	if (FA->ring_flex_active && FA->ring_flex_template &&
	    FA->ring_flex_template->has_rings() && FA->ring_n_sugars > 0) {
		const ligand_ring_flex::RingFlexGenes& tmpl = *FA->ring_flex_template;
		std::vector<float> phases(FA->ring_cur_phases,
		                          FA->ring_cur_phases + FA->ring_n_sugars);
		sugar_pucker::apply_sugar_puckers(
			atoms, tmpl.sugar_ring_indices, phases, tmpl.sugar_types);
	}

	// Shared restore for any error path that mutated FA->ori / atoms / rotamers.
	// Must run before return so the next evaluation cannot inherit contamination.
	auto restore_ic2cf_baseline = [&]() {
		FA->ori[0] = ori_save[0];
		FA->ori[1] = ori_save[1];
		FA->ori[2] = ori_save[2];
		for (const auto& sa : saved_atoms) {
			atoms[sa.idx] = sa.value;
		}
		for (const auto& sr : saved_res_rots) {
			residue[sr.idx].rot = sr.rot;
		}
	};

	/* rebuild cartesian coordinates of optimized residues*/
	for(i=0;i<FA->nors;i++){ //number of optimized residues
		if (!buildcc(FA,atoms,FA->nmov[i],FA->mov[i])) {
			// Explicit reconstruction failure: restore baseline, do not score
			// with half-mutated FA / atom state (serial eval contamination).
			restore_ic2cf_baseline();
			cfstr cf_bad{};
			cf_bad.wal = 1.0e12;
			cf_bad.rclash = 1;
			return cf_bad;
		}
	}

	// Out-of-bounds penalty: if any moved atom lands >200Å beyond the protein
	// bounding box, the ligand has escaped the grid.  Restore FA->ori and the
	// moved-atom coordinates so serial callers do not inherit corrupted state,
	// then return maximum penalty — the chromosome stays in the population but
	// ranks last, and the GA evolves away from it naturally.
	{
		const float margin = 200.0f;
		bool oob = false;
		for (int r = 0; r < FA->nors && !oob; ++r) {
			for (int m = 0; m < FA->nmov[r] && !oob; ++m) {
				int ai = FA->mov[r][m];
				for (int j = 0; j < 3; ++j) {
					if (atoms[ai].coor[j] < FA->globalmin[j] - margin ||
					    atoms[ai].coor[j] > FA->globalmax[j] + margin) {
						oob = true;
						break;
					}
				}
			}
		}
		if (oob) {
			restore_ic2cf_baseline();
			cfstr cf_oob{};
			cf_oob.com = 99999.0;
			return cf_oob;
		}
	}

	// A3: reuse thread-local intraclashes buffer
	scr.intraclashes.clear();
	std::vector<std::pair<int,int>>& intraclashes = scr.intraclashes;
	bool error;
	double penalty = vcfunction(FA,VC,atoms,residue,intraclashes,&error);
	if(error){
		// Fail-closed: restore FA ori, moved atoms, and residue rotamers so a
		// subsequent evaluation cannot see a contaminated pose from this call.
		restore_ic2cf_baseline();
		// wal=penalty is the only nonzero energy field; rclash=1 marks clash.
		cfstr cf_clash{};
		cf_clash.wal    = penalty;
		cf_clash.rclash = 1;
		return cf_clash;
	}
	
	// cf already value-initialized; re-zero explicit energy accumulators.
	cf = cfstr{};
	cf.com = 0.0;
	cf.wal = 0.0;
	cf.sas = 0.0;
	cf.con = 0.0;
	cf.elec = 0.0;
	cf.hbond = 0.0;
	cf.gist_desolv = 0.0;
	cf.metal_coord = 0.0;
	cf.pb_clash = 0.0;
	cf.h_rep = 0.0;
	cf.entropy = 0.0;
	cf.rclash = 0;
    
	for(i=0;i<FA->num_optres;i++){
    
		resid* res = &residue[FA->optres[i].rnum];
		
		// flexible side-chain optimization
		if ( !FA->optres[i].type ) {
  
			if ( FA->optres[i].cf.rclash == 1 ) { 

				/*
				  printf("%s %c %d is clashing\n",
				  residue[FA->optres[i].rnum].name,
				  residue[FA->optres[i].rnum].chn,
				  residue[FA->optres[i].rnum].number);
				*/
	
				rclash = 1; 
				
			}
      
		}else{
			
			//int fatm = res->fatm[0];
			if(FA->deelig_flex){
				std::vector< std::pair<int,int> >::iterator it;
				for(it=intraclashes.begin(); it!=intraclashes.end(); ++it)
				{
					for(k=1; k<=res->fdih; k++){
						deelig_list[k] = -1000;
					}
					
					// flex bonds list
					int fbindex = 0;
					int* fblist = res->shortflex[it->first][it->second];
					
					//printf("between[%d][%d]\n", atoms[it->first+fatm].number, atoms[it->second+fatm].number);
					//cout << fblist[fbindex] << endl;
					/*
					printf("fblist = [");
					while(fblist[fbindex] != -1){
						printf("%d,", fblist[fbindex]);
						fbindex++;
					}
					printf("]\n");
					*/
					
					fbindex = 0;
					while(fblist[fbindex] != -1){
						if(atoms[res->bond[fblist[fbindex]]].par != NULL){
							deelig_list[fblist[fbindex]] =
								(int)(atoms[res->bond[fblist[fbindex]]].dih + 0.5);
						}
						fbindex++;
					}
										
					struct deelig_node_struct* node = FA->deelig_root_node;
				[[maybe_unused]] bool add = false;
					
						for(k=1; k<=res->fdih; k++){
						std::map<int, struct deelig_node_struct*>::iterator it;
						it = node->childs.find(deelig_list[k]);
						
						if(it == node->childs.end()){
							struct deelig_node_struct* deelig_child_node = new struct deelig_node_struct;
							
							//if(k==1) cout << "new node added " << deelig_list[k] << endl;
							node->childs[deelig_list[k]] = deelig_child_node;
							
							deelig_child_node->parent = node;
							node = deelig_child_node;
							add = true;
						}else{
							node = it->second;
						}
					}
					
					/*
					if(add) { 
						printf("deelig list = [");
						for(k=1; k<=res->fdih; k++){
							printf("%d,", deelig_list[k]);
						}
						printf("]\n");
						
						nbranch++; cout << "total branches " << nbranch << endl;
					}
					*/
				}
			}
		}
		
		/*
		  printf("optres[%2d].cf  .wal = %.3f\n               .com = %.3f\n               .sas = %.3f\n               .con = %.3f\n",
		  i,FA->optres[i].cf.wal,FA->optres[i].cf.com,FA->optres[i].cf.sas,FA->optres[i].cf.con);
		*/
        
		//sum += (FA->optres[i].cf.com - FA->optres[i].cf.wal + FA->optres[i].cf.sas - FA->optres[i].cf.con);
    
		cf.com += FA->optres[i].cf.com;
		cf.wal += FA->optres[i].cf.wal;
		cf.sas += FA->optres[i].cf.sas;
		cf.con += FA->optres[i].cf.con;
		// Production aggregator must include every get_cf_evalue() term.
		// Historical bug: elec and gist_desolv were zeroed above then never
		// summed from optres, silently dropping electrostatics / GIST desolv.
		cf.elec += FA->optres[i].cf.elec;
		cf.gist_desolv += FA->optres[i].cf.gist_desolv;
		cf.metal_coord += FA->optres[i].cf.metal_coord;
		cf.hbond += FA->optres[i].cf.hbond;
		cf.entropy += FA->optres[i].cf.entropy;
		cf.pb_clash += FA->optres[i].cf.pb_clash;

	}

  
	// add rotamer list to dee list
	// When running inside an OpenMP parallel region each thread has its own
	// copy of FA (thread-local), so DEE updates are serialised via a critical
	// section and written back to FA directly (the thread-local FA shares the
	// psFlexDEENode pointer with the master; the critical section prevents
	// concurrent linked-list corruption).
#ifdef _OPENMP
	if (FA->useflexdee > 0 && rclash && !omp_in_parallel()) {
#else
	if (FA->useflexdee > 0 && rclash) {
#endif
    
		NEW( psFlexDEENode, sFlexDEE_Node );

		psFlexDEENode->rotlist = (int*)malloc(FA->nflxsc_real*sizeof(int));
    
    
		// fill rotamer list
		k=0;
		rotflag=0;
    
		for(j=0;j<FA->nflxsc;j++){
      
			if(residue[FA->flex_res[j].inum].trot > 0   &&
			   FA->flex_res[j].cflag != 0){
	
				psFlexDEENode->rotlist[k++] = residue[FA->flex_res[j].inum].rot;
	
				if(residue[FA->flex_res[j].inum].rot != 0) { rotflag=1; }
			}
      
		}
    
    
		// do not add initial conformation to DEE list
		if ( rotflag ) {
      
			/*
			  printf("\n-----------------\nCreating new node...\n");
			  printf("DEE to add = ");for(k=0;k<FA->nflxsc_real;k++){printf("%3d",psFlexDEENode->rotlist[k]);}printf("\n");
	
			  dee_print(FA->psFlexDEENode,FA->nflxsc_real);
	
			  //getchar();
			  */
      
    
			if( FA->psFlexDEENode ) {
	
				//FA->psFlexDEENode = FA->psFlexDEENode->last;
	
				while ( FA->psFlexDEENode->next != NULL ) {
					FA->psFlexDEENode = FA->psFlexDEENode->next;
				}
	
				dee_val = dee_pivot(psFlexDEENode,&FA->psFlexDEENode,1,FA->FlexDEE_Nodes,(int)((FA->FlexDEE_Nodes+1)/2),FA->FlexDEE_Nodes,FA->nflxsc_real);
	
				if ( dee_val == 1 ) {
	  
					if ( FA->psFlexDEENode->next == NULL ) {
	    
						psFlexDEENode->next = NULL;
						psFlexDEENode->prev = FA->psFlexDEENode;
						FA->psFlexDEENode->next = psFlexDEENode;
	    
						psFlexDEENode->first = FA->psFlexDEENode;
	    
						dee_last(FA->psFlexDEENode,psFlexDEENode);
	    
					} else {
	    
						psFlexDEENode->first = FA->psFlexDEENode->first;
						psFlexDEENode->last = FA->psFlexDEENode->last;
	    
						psFlexDEENode->next = FA->psFlexDEENode->next;
						psFlexDEENode->prev = FA->psFlexDEENode;
						FA->psFlexDEENode->next = psFlexDEENode; 
						psFlexDEENode->next->prev = psFlexDEENode;
	    
					}
	  
					FA->FlexDEE_Nodes++;
	  
				} else if ( dee_val == -1 ) {
	  
					if ( FA->psFlexDEENode->prev == NULL ) {
	    
						psFlexDEENode->prev = NULL;
						psFlexDEENode->next = FA->psFlexDEENode;
						FA->psFlexDEENode->prev = psFlexDEENode;
	    
						psFlexDEENode->last = FA->psFlexDEENode;
	    
						dee_first(FA->psFlexDEENode,psFlexDEENode);
	    
					} else {
	    
						psFlexDEENode->first = FA->psFlexDEENode->first;
						psFlexDEENode->last = FA->psFlexDEENode->last;
	    
						psFlexDEENode->prev = FA->psFlexDEENode->prev;
						psFlexDEENode->next = FA->psFlexDEENode;
						FA->psFlexDEENode->prev = psFlexDEENode; 
						psFlexDEENode->prev->next = psFlexDEENode;
	    
					}
	  
					FA->FlexDEE_Nodes++;
	  
				} else {
	  
					FREE(psFlexDEENode);
	  
				}
	
			} else {
	
				FA->psFlexDEENode = psFlexDEENode;
	
				FA->psFlexDEENode->next = NULL;
				FA->psFlexDEENode->prev = NULL;
	
				FA->psFlexDEENode->first = FA->psFlexDEENode;
				FA->psFlexDEENode->last = FA->psFlexDEENode;
	
				FA->FlexDEE_Nodes++;
	
			}
      
		}
    
	}

	// Restore FA->ori only — NOT atoms[] — on the normal scoring exit.
	//
	// FA->ori is global scratch that drifts cumulatively if not reset; always
	// restore it so the next ic2cf call uses the correct receptor-centre frame.
	//
	// atoms[] and residue[].rot are intentionally NOT restored here.
	// The GA parallel path works on thread-private copies that are discarded
	// after scoring, so atom state persistence doesn't matter.
	// The serial output path (cluster.cpp:240→302, top.cpp:1242→1287) calls
	// ic2cf precisely to populate atoms[] with the final docked Cartesian pose,
	// then immediately passes atoms[] to write_pdb.  Restoring atoms on normal
	// exit would silently write the pre-call (initial) structure — every cluster
	// output PDB would be a copy of the start conformation, not the docked pose.
	// Atom restore belongs only on the OOB/penalty early-exit path above, where
	// we must leave a clean baseline for the next chromosome evaluation.
	FA->ori[0] = ori_save[0];
	FA->ori[1] = ori_save[1];
	FA->ori[2] = ori_save[2];

	cf.h_rep = compute_ligand_h_rep(FA, atoms);
	compute_ligand_dsvib(FA, atoms, &cf);

	return cf;

}

#ifdef _WIN32
double get_apparent_cf_evalue(cfstr* cf) {
#else
	double get_apparent_cf_evalue(cfstr* cf) {
#endif
		return cf->com + cf->wal + cf->sas + cf->elec + cf->hbond + cf->gist_desolv + cf->metal_coord + cf->entropy + cf->pb_clash;
	}

#ifdef _WIN32
	double get_cf_evalue(cfstr* cf, FA_Global* FA) {
#else
		double get_cf_evalue(cfstr* cf, FA_Global* FA) {
#endif
			double total = cf->com + cf->wal + cf->sas + cf->con + cf->elec
			             + cf->hbond + cf->gist_desolv + cf->metal_coord + cf->entropy + cf->pb_clash;
			if (FA && FA->tencom_weight > 0.0f) {
				total += static_cast<double>(FA->tencom_weight) * cf->h_rep;
			}
			// -T*dS_vib enters with UNIT coefficient. There is deliberately no
			// weight field to multiply here: the term is already a free energy in
			// kcal/mol, and a fitted weight is exactly what makes the h_rep line
			// above a heuristic instead of physics. Any weight sweep is a
			// sensitivity diagnostic run from outside, never a fit inside the CF.
			// Guarded on status == 1 so a refused computation (degenerate spectrum,
			// mode-count mismatch, missing temperature) contributes nothing rather
			// than a zero that reads like a measured zero.
			if (FA && FA->dsvib_mode > 0 && cf->dsvib_status == 1) {
				total += cf->minus_T_dsvib;
			}
			return total;
		}
