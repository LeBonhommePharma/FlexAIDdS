#include "Vcontacts.h"
#include "soft_wall.h"
#include "GISTEvaluator.h"
#include "HBondEvaluator.h"
#include "hbond_potential.h"
#include "metal_coordination.h"
#include "GISTGrid.h"
#include "ProtocolConfig.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <unordered_map>

#define DEBUG_LEVEL 0

// Coulomb constant: 332.0637 kcal·Å/(mol·e²)
#define KCOULOMB 332.0637

// ── File-scope magic-static env flags (hoisted from vcfunction body) ─────────
// Each reads the environment exactly once (thread-safe C++11 static init) so
// std::getenv() is never called inside the hot per-contact loop.
static const bool softcore_wal =
    (std::getenv("FLEXAIDDS_SOFTCORE_WAL") != nullptr);
static const double softcore_floor_frac = [](){
    const char* s = std::getenv("FLEXAIDDS_SOFTCORE_FLOOR");
    double v = 0.5;   // default: hard wall below 50 % of contact radius
    if(s){ double p = strtod(s, nullptr); if(p > 0.0 && p < 0.7) v = p; }
    return v;
}();
static const bool dist_weight_con =
    (std::getenv("FLEXAIDDS_DIST_WEIGHT_CON") != nullptr);
static const float con_r0 = []() {
    const char* s = std::getenv("FLEXAIDDS_CON_R0");
    float v = 3.5f;
    if(s){ float p = strtof(s, nullptr); if(p > 0.0f){ v = p; } }
    return v;
}();
static const bool no_sas = (std::getenv("FLEXAIDDS_NO_SAS") != nullptr);

// FLEXAIDDS_POLAR_DESOLV_WEIGHT (float >= 0, default 0.0 = OFF): per-ligand-atom
// polar-burial desolvation penalty, folded into the cfs->sas channel (which
// ic2cf already sums into the CF total, so no new struct field / no ic2cf edit
// is needed). Scoped fix for the com-over-burial CF-inversion (1G9V): a polar
// ligand atom that gets buried on binding WITHOUT forming a compensating
// H-bond / salt-bridge / metal coordination pays a desolvation cost. Well
// satisfied polar atoms (s_i -> 1) pay ~0, so deep, H-bonded natives are
// protected. Default 0.0 -> the block is skipped entirely (zero overhead,
// bit-identical). Read once (magic static, thread-safe) — never in the loop.
//   penalty_i = w * P(type_i) * b_i * (1 - s_i)   (>= 0, added to cfs->sas)
//     b_i = 1 - SAS/surfA          buried fraction, clamped to [0,1]
//     s_i = min(1, (|E_hb_i|+|E_elec_i|+|E_mc_i|) / E_ref),  E_ref = 2.5
static const double polar_desolv_weight = [](){
    const char* s = std::getenv("FLEXAIDDS_POLAR_DESOLV_WEIGHT");
    double v = 0.0;
    if(s){ double p = strtod(s, nullptr); if(p > 0.0) v = p; }
    return v;
}();
// Reference satisfaction energy (kcal-scale) that saturates s_i to 1.
static const double polar_desolv_eref = 2.5;
// P(type_i): polar propensity by VCT row (atoms[].type).
//   0.0  C rows 1-5, S rows 17-21, halogens 23-26, DUMMY 39 (and any other)
//   1.0  neutral polar N/O: rows 6-8, 10-14, 16
//   1.75 charged: N.4 (row 9), O.co2 (row 15), P.3 (row 22)
static inline double polar_desolv_ptype(int t){
    switch(t){
        case 9:  case 15: case 22:
            return 1.75;
        case 6:  case 7:  case 8:
        case 10: case 11: case 12: case 13: case 14:
        case 16:
            return 1.0;
        default:
            return 0.0;
    }
}


// ── pb_clash receptor grid cache (Route A hoist) ─────────────────────────────
// The receptor is rigid across all CF evals in one dock session. The cell-list
// grid over receptor atoms (atoms where optres==NULL) is therefore loop-invariant
// and can be built once per dock rather than once per eval.
//
// Invalidation key: (atm_cnt, pb_clash_ratio).  atm_cnt changes between docks;
// pb_clash_ratio is a config value that only changes at dock-init time.
// Within a single dock the receptor atom positions are fixed, so a matching key
// guarantees a valid grid.
//
// Always active whenever pb_clash_weight>0.  This was previously gated behind
// FLEXAIDDS_PB_CLASH_GRID_HOIST (default OFF) as a caution during initial
// testing; the hoist is verified bit-identical, and leaving it opt-in meant a
// config that enabled pb_clash but forgot the env var silently reproduced the
// v134 timeout.  CF math is unchanged — only grid *construction* is skipped on
// evals 2…N.

struct PBClashGridCache {
    // Invalidation key
    int    atm_cnt   = -1;
    double ratio     = -1.0;
    // Cached grid data
    double rmin[3]   = {};
    int    dim[3]    = {};
    double cell      = 0.0;
    double inv_cell  = 0.0;
    std::vector<int>    rec_idx;                        // rigid receptor atom indices
    std::vector<double> rec_vdw;                        // precomputed vdW per rec atom (parallel to rec_idx)
    std::vector<char>   rec_metal;                      // 1 if rec atom is a coordinating metal (parallel to rec_idx)
    std::vector<char>   rec_heavy;                      // 1 if rec atom is a heavy (non-H) atom (parallel to rec_idx)
    std::unordered_map<int, std::vector<int>> grid;     // cell-list: flat cell key → rec_idx positions
    bool valid = false;

    bool matches(int a, double r) const { return valid && (a == atm_cnt) && (r == ratio); }
};
static thread_local PBClashGridCache pb_cache;

// FLEXAIDDS_PB_METAL_CARVEOUT (truthy, default OFF): exclude ligand/receptor
// pairs in which either partner is a coordinating metal (Zn, Fe, Mg, Ca, Mn,
// Co, Ni, Cu) from the pb_clash penalty. Genuine coordination bonds sit at
// ~1.9-2.3 A, well inside pb_clash_ratio*(vdw_i+vdw_j), so without this they
// are scored as hard clashes and fight the cf.metal_coord Morse reward.
// Default OFF → the pair set is unchanged and CF is bit-identical.
static const bool pb_metal_carveout =
    (std::getenv("FLEXAIDDS_PB_METAL_CARVEOUT") != nullptr &&
     std::getenv("FLEXAIDDS_PB_METAL_CARVEOUT")[0] != '0');

// FLEXAIDDS_PB_VDW_CACHED (truthy, default OFF): read the precomputed
// atoms[].pb_vdw_radius (populated once in update_optres() from the structure
// file's element column) instead of re-running posebusters_vdw_radius() on
// get_element(atoms[].type) per ligand atom per eval. Vcontacts.cpp already
// reads the cached field; this brings vcfunction.cpp onto the same source.
//
// NOT a pure perf swap — the two sources are NOT fully interchangeable, which is
// why this is flag-gated rather than applied unconditionally. The cause is the
// READER tables, not the PoseBusters table: `element` records the true element,
// while `type` is a VCT slot, and the readers deliberately map some elements
// onto a different element's slot or onto DUMMY. type -> element via
// get_element() therefore loses the true element for exactly those atoms.
// tests/test_pb_vdw_parity.cpp encodes the reader tables and pins the result:
//   * every element mapped onto its OWN row: identical on both paths (no-op).
//   * iodine: readers assign the BROMINE row (type 25) on purpose, so the type
//     path gives 1.90 while the element path gives the true iodine 2.10.
//   * Na / K: no row in the 40-type table, so they land in DUMMY and the type
//     path falls back to the NRG radius 1.70, against 2.40 / 2.80 by element.
//   * hydrogen: lands in DUMMY too, but the readers give H an NRG radius of
//     exactly 1.20, which IS the PoseBusters H radius — so H does NOT diverge.
// Enabling this widens the clash cutoff for I/Na/K toward the chemically correct
// radius, but it is a scoring change and must be benchmarked before it becomes a
// default. Default OFF → bit-identical.
static const bool pb_vdw_cached =
    (std::getenv("FLEXAIDDS_PB_VDW_CACHED") != nullptr &&
     std::getenv("FLEXAIDDS_PB_VDW_CACHED")[0] != '0');

// FLEXAIDDS_PB_POCKET_DEBUG: emit the per-eval [PB_POCKET] trace line.
static const bool pb_pocket_debug =
    (std::getenv("FLEXAIDDS_PB_POCKET_DEBUG") != nullptr);

// FLEXAIDDS_WAL_COERCIVE: remove WAL_CONTACT_CAP ceiling on the soft-core
// fitness wall so deep clashes can overcome unbounded CF.com overpacking.
// Default OFF → bit-identical to current behaviour.
static const bool wal_coercive =
    (std::getenv("FLEXAIDDS_WAL_COERCIVE") != nullptr &&
     std::getenv("FLEXAIDDS_WAL_COERCIVE")[0] != '0');
// FLEXAIDDS_WAL_STIFF: override k_wal (default WAL_CONTACT_CAP=50) for
// stiffness sweeps at benchmark time.  0 or unset → existing default.
static const double wal_stiff = [](){
    const char* s = std::getenv("FLEXAIDDS_WAL_STIFF");
    if (!s) return 0.0;
    double v = strtod(s, nullptr);
    return (v > 0.0) ? v : 0.0;
}();


// FLEXAIDDS_CONTACTS_EPOCH: O(1) clear of FA->contacts (avoids a
// MAX_ATOM_NUMBER-sized memset every eval; see below). Default off —
// legacy memset-every-call behavior is unchanged unless this is set.
static const bool contacts_epoch_mode =
    (std::getenv("FLEXAIDDS_CONTACTS_EPOCH") != nullptr);

double vcfunction(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue, std::vector<std::pair<int,int> > & intraclashes, bool* error)
{
	int    rnum=0;
	int    type=1;
	
	// reset all values pointed
	// contacts[] is a per-atom "already visited this eval" flag, read via
	// truthy check and written as 1 (see the two use sites below). In the
	// default (legacy) path we still memset it to 0 every call — bit-identical
	// to before. Under FLEXAIDDS_CONTACTS_EPOCH we instead stamp visited slots
	// with a monotonically increasing per-eval epoch and treat any slot whose
	// stamp != current epoch as unset, making the "clear" O(1) instead of an
	// O(MAX_ATOM_NUMBER) memset. contacts_epoch lives in FA_Global (copied
	// per-thread alongside the contacts pointer in gaboom.cpp/VoronoiCFBatch.h),
	// so each thread's counter tracks only its own private contacts buffer.
	if(contacts_epoch_mode){
		// Wraparound guard: an int epoch could in principle repeat a stale
		// stamp after ~2^31 evals (never reached in one restart, but this
		// keeps arbitrarily long campaigns correct). Fall back to one real
		// clear and restart the counter from 1.
		if(FA->contacts_epoch >= INT_MAX - 1){
			memset(FA->contacts,0,MAX_ATOM_NUMBER*sizeof(int));
			FA->contacts_epoch = 0;
		}
		++FA->contacts_epoch;
	} else {
		memset(FA->contacts,0,MAX_ATOM_NUMBER*sizeof(int));
	}
	memset(FA->contributions,0,FA->ntypes*FA->ntypes*sizeof(float));
	
	// reset CF values
	for(int j=0; j<FA->num_optres; ++j){
		FA->optres[j].cf.rclash=0;
		FA->optres[j].cf.wal=0.0;
		FA->optres[j].cf.com=0.0;
		FA->optres[j].cf.con=0.0;
		FA->optres[j].cf.elec=0.0;
		FA->optres[j].cf.gist=0.0;
		FA->optres[j].cf.hbond=0.0;
		FA->optres[j].cf.gist_desolv=0.0;
		FA->optres[j].cf.metal_coord=0.0;
		FA->optres[j].cf.entropy=0.0;
		FA->optres[j].cf.pb_clash=0.0;
		FA->optres[j].cf.sas=0.0;
		FA->optres[j].cf.totsas=0.0;
	}
	
	double permea = (double)FA->permeability;
	double dee_clash = (double)FA->dee_clash;

	// Lever 2: per-optres VCT contact tally for intensive CF.com normalization.
	// Only allocated/used when the flag is set so the default path is unchanged.
	std::vector<long> vct_ncon;
	if(FA->vct_normalize_contacts){ vct_ncon.assign(FA->num_optres, 0); }
	
	// allocate
	//float  matrix[FA->ntypes*FA->ntypes];
	/*
	// empty surface matrix values
	for(i=0;i<FA->ntypes;++i){
	for(j=0;j<FA->ntypes;++j){
	matrix[i][j]=0.0;
	}
	}
	*/  
    
	//printf("=============NEW INDIVIDUAL==============\n");
    
	double clash_value;
	*error = false;
	int rv = Vcontacts(FA,atoms,residue,VC,&clash_value,false);
	if(rv < 0){
		*error = true;
		
		for(int i=0;i<FA->atm_cnt_real;i++){
			if(VC->Calc[i].score){
				VC->Calc[i].atom = NULL;
			}
		}
		
		if(!FA->vindex){ free(VC->box); }
		
		if(rv == -1){
			FA->skipped++;
			return(POLYHEDRON_PENALTY);
		}else if(rv == -2){
			FA->clashed++;
			return(clash_value);
		}
	}
	
	for(int i=0; i<FA->atm_cnt_real; ++i) {
		if(VC->Calc[i].atom == NULL) continue;
		
		cfstr* cfs = NULL;
#if DEBUG_LEVEL > 0
		cfstr cfs_atom;
#endif

		// atom from which contacts are calculated
		int atomzero = FA->num_atm[VC->Calc[i].atom->number];
		
		// number of constraints for atomzero
		int nconszero = atoms[atomzero].ncons;
		
		if(atoms[atomzero].optres != NULL)
		{
			// the residue optimizable
			rnum = atoms[atomzero].optres->rnum;
			type = atoms[atomzero].optres->type;
			cfs = &atoms[atomzero].optres->cf;
		}
		else
		{
			continue;
		}
		
		//printf("-------------------------------\nAtom[%4d]=%s in Residue[%4d]\n",VC->Calc[i].atom->number,atoms[FA->num_atm[VC->Calc[i].atom->number]].name,VC->Calc[i].residue->number);
		
		
#if DEBUG_LEVEL > 0
		double com_atm=0.0;
		double Ewall_atm=0.0;
#endif
		
		double radA  = (double)VC->Calc[i].atom->radius;
		double radoA = radA + Rw;
		
		double SAS = 4.0*PI*radoA*radoA;
       		double surfA = SAS;
		
		if(FA->useacs && atoms[atomzero].acs < 0.0){
			// accessible contact surface with solvent/atom		
			// ACS = Total surface area - surface areas of bonded contacts (atoms with a bond/angle between them)
			atoms[atomzero].acs = surfA;
		
			int currindex = VC->ca_index[i];
			
			while(currindex != -1) {
				int atomcont = FA->num_atm[VC->Calc[VC->ca_rec[currindex].atom].atom->number];
				int intramolecular = atoms[atomcont].ofres == atoms[atomzero].ofres;
				
				// get first atom of residue
				int fatm = residue[rnum].fatm[0];
				
				if(intramolecular && residue[rnum].bonded != NULL &&
				   residue[rnum].bonded[atomcont-fatm][atomzero-fatm] >= 0)
				{
					atoms[atomzero].acs -= VC->ca_rec[currindex].area;
				}
				
				currindex = VC->ca_rec[currindex].prev;
			}
			
			if(atoms[atomzero].acs < 0.0){ atoms[atomzero].acs = 0.0f; }
			//printf("after ACS=%.3f\n", ACS);
		}
		
#if DEBUG_LEVEL > 0
		cfs_atom.sas = 0.0;
#endif

		if(atoms[atomzero].ncons > 0){

			for(int j=0;j<atoms[atomzero].ncons;j++){
	
				/*
				  double radC = atoms[atomzero].number==atoms[atomzero].cons[j]->anum1 ? 
				  (double)atoms[FA->num_atm[atoms[atomzero].cons[j]->anum2]].radius:
				  (double)atoms[FA->num_atm[atoms[atomzero].cons[j]->anum1]].radius;
				*/
				
				// maximum penalty value (starting penalty)
				// default value if atoms are not interacting
				//cfs->con += KANGLE*(radA+radC+2.0*Rw);

				if(atoms[atomzero].cons[j]->type == 1){
					cfs->con += KDIST;
				}
				
				//printf("constraint for atom[%d]: %.3f\n", atoms[atomzero].number,cfs->con);
			}
      
			/*
			  printf("default constraint value: %.3f\n",cfs->con);
			  getchar();
			*/
		}


#if DEBUG_LEVEL > 0
		printf("==================================================================================\n");
		printf("ATOM :: RES C RNUM  ANUM  T  RAD ::   COMPL  (W)  DIST   AREA ::     CF.COM     CF.WAL\n");
		printf("----------------------------------------------------------------------------------\n");
		printf("ATOM :: %3s %c %4d %5d %2d %4.2f\n",
		       residue[atoms[atomzero].ofres].name,
		       residue[atoms[atomzero].ofres].chn,
		       residue[atoms[atomzero].ofres].number,
		       atoms[atomzero].number,
		       atoms[atomzero].type,
		       atoms[atomzero].radius);
#endif

	[[maybe_unused]] int contnum = 0;  // number of contacts (excluding bloops away atoms)
		int metal_cn_count = 0;  // coordination number counter for metal atoms
		// Per-ligand-atom satisfaction accumulators for the polar-burial
		// desolvation penalty (FLEXAIDDS_POLAR_DESOLV_WEIGHT). Signed sums over
		// this atom's contacts; |.| taken when forming s_i near cfs->sas below.
		// Only touched when the term is enabled (weight != 0) — zero overhead OFF.
		double pd_hb = 0.0, pd_elec = 0.0, pd_mc = 0.0;
		int currindex = VC->ca_index[i];
		
		while(currindex != -1) {
			
			double radB  = (double)VC->Calc[VC->ca_rec[currindex].atom].atom->radius;
			// double radoB = radB + Rw;
			// double surfB = 4.0*PI*radoB*radoB;
			
			double rAB   = radA+radB;
			
			//double complementarity = 0.0;
			double area = VC->ca_rec[currindex].area;
			
			struct energy_matrix* energy_matrix = &FA->energy_matrix[(VC->Calc[i].atom->type-1)*FA->ntypes +
										 (VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)];

			//double yval = get_yval(energy_matrix,area/((surfA+surfB)/2.0));
			// always use normalized areas in density functions
			double yval = get_yval(energy_matrix,area/surfA);
			
			SAS -= area;
			
			// number of contacts counter
			contnum++;

#if DEBUG_LEVEL > 0
			cfs_atom.com  =  0.0;
			cfs_atom.wal  =  0.0;
#endif
			// atom in contact with atom zero
			int atomcont = FA->num_atm[VC->Calc[VC->ca_rec[currindex].atom].atom->number];
			
			int intramolecular = 0;
			int intraresidue = 0;
			if(atoms[atomcont].ofres == atoms[atomzero].ofres){
				intraresidue = 1;
				intramolecular = 1;
				
			}else if(residue[atoms[atomcont].ofres].type == 0 &&
				 residue[atoms[atomzero].ofres].type == 0){

				intramolecular = 1;
			}
			
			// get first atom of residue
			int fatm = residue[rnum].fatm[0];
		
			// is contact atom bonded to atom zero
			// if YES, skip contact atom
			if(intraresidue && residue[rnum].bonded != NULL)
			{
				// always skip atoms forming a bond or angle with each other
				if(residue[rnum].bonded[atomcont-fatm][atomzero-fatm] >= 0)
				{

#if DEBUG_LEVEL > 2
					printf("    (B) %3s %c %4d %5d %2d %4.2f :: %7.4f (%s) %6.2f %6.2f :: %10.3f %10.3f\n",
					       VC->Calc[VC->ca_rec[currindex].atom].residue->name,
					       VC->Calc[VC->ca_rec[currindex].atom].residue->chn,
					       VC->Calc[VC->ca_rec[currindex].atom].residue->number,
					       VC->Calc[VC->ca_rec[currindex].atom].atom->number,
					       VC->Calc[VC->ca_rec[currindex].atom].atom->type,
					       VC->Calc[VC->ca_rec[currindex].atom].atom->radius,
					       
					       yval, energy_matrix->weight ? "Y": "N",
					       VC->ca_rec[currindex].dist,
					       VC->ca_rec[currindex].area,
					       
					       cfs_atom.com, cfs_atom.wal);

#endif
					currindex = VC->ca_rec[currindex].prev;
					continue;	  
				}
			}
			
			if(contacts_epoch_mode
			   ? (FA->contacts[VC->Calc[VC->ca_rec[currindex].atom].atom->number] == FA->contacts_epoch)
			   : (FA->contacts[VC->Calc[VC->ca_rec[currindex].atom].atom->number] != 0)){
				//printf("%d already calculated\n",VC->Calc[VC->ca_rec[currindex].atom].atom->number );
				currindex = VC->ca_rec[currindex].prev;
				continue;
			}
			
			// covalently bonded flag
			bool covalent = false;
			
			// number of constraints for contact atom
			int nconscont = atoms[atomcont].ncons;
      			double dist_opt = 0.0;
			
			// do contacting atoms have the same constraint
			constraint* cons = NULL;
			if(nconszero > 0 && nconscont > 0){
				for(int j=0;j<nconszero;j++){
					for(int k=0;k<nconscont;k++){
						if(atoms[atomzero].cons[j]->id == atoms[atomcont].cons[k]->id){
							cons = &FA->constraints[atoms[atomzero].cons[j]->id];
							break;
						}
					}
					if(cons != NULL){break;}
				}
				
				if(cons != NULL){
					
					if(cons->type == 1){
						
						covalent = true;
						dist_opt = cons->bond_len;
						
						cfs->con -= KDIST * GetValueFromGaussian(VC->ca_rec[currindex].dist,dist_opt,cons->max_dist);						
					}
				}
	
				//printf("constraint[%d] applies.\n",cons->id);
			}
			
			//	coorB = VC->Calc[VC->ca_rec[currindex].atom].coor;

			// CHECK IF CLASH
			float clash_distance = permea*rAB;
			if(covalent && permea*dist_opt < permea*rAB){
				clash_distance = permea*dist_opt;
			}
                        
			if (VC->ca_rec[currindex].dist < clash_distance){
				
				// Fast multiplication chain for r⁻¹² (replaces slow pow() calls)
			double d  = VC->ca_rec[currindex].dist;
			double d2 = d * d; double d4 = d2 * d2; double d6 = d4 * d2;
			double inv_d12 = 1.0 / (d6 * d6);
			double cr = permea * rAB;
			double cr2 = cr * cr; double cr4 = cr2 * cr2; double cr6 = cr4 * cr2;
			double inv_cr12 = 1.0 / (cr6 * cr6);
			double Ewall = KWALL * (inv_d12 - inv_cr12);

			// ── Softcore WAL (FLEXAIDDS_SOFTCORE_WAL, default OFF) ─────────────
			// The bare r^-12 wall spikes too steeply just inside the contact
			// radius cr = permea*rAB (= r_min, where Ewall=0). A near-native pose
			// with 0.2-0.4 Å crystal coordinate error can then score WORSE than a
			// decoy: the spike dominates the complementarity term. When enabled,
			// for r < r_softcrit = 0.7*cr we replace r^-12 with a C1-continuous
			// downward-opening parabola that matches BOTH the value V(r_softcrit)
			// and the slope V'(r_softcrit) of the r^-12 form at the transition
			// (no value jump, no kink) and levels off to a finite plateau as
			// r->0 (parabola maximum sits at r=0). This is in ADDITION to the
			// WAL_CONTACT_CAP overflow guard below, which still clamps fitness.
			// softcore_wal / softcore_floor_frac are file-scope statics (hoisted).
			if(softcore_wal){
				const double r_softcrit = 0.7 * cr;
				// ── Hard-floor gate (FLEXAIDDS_SOFTCORE_FLOOR, default 0.5) ────────
				// The parabola plateau can sit below WAL_CONTACT_CAP=50 for deeply
				// buried large-atom pairs, letting catastrophically buried poses score
				// nearly free.  For d < r_hardfloor we revert to the raw r^-12 wall
				// so the existing per-contact cap of 50 applies.  The parabola only
				// operates in the moderate-penetration band [r_hardfloor, r_softcrit).
				// softcore_floor_frac is a file-scope static (hoisted).
				const double r_hardfloor = softcore_floor_frac * cr;
				if(d < r_softcrit){
					if(d >= r_hardfloor){
						// Moderate penetration: C1-continuous parabola (existing logic)
						double rsc2 = r_softcrit * r_softcrit;
						double rsc4 = rsc2 * rsc2; double rsc6 = rsc4 * rsc2;
						double inv_rsc12 = 1.0 / (rsc6 * rsc6);
						double V_sc = KWALL * (inv_rsc12 - inv_cr12);
						double absVp = 12.0 * KWALL * (inv_rsc12 / r_softcrit);
						double u = r_softcrit - d;             // >= 0 inside softcore
						Ewall = V_sc + absVp * u
						        - (absVp / (2.0 * r_softcrit)) * u * u;
					}
					// else: d < r_hardfloor — leave Ewall as hard r^-12; WAL_CONTACT_CAP=50 applies below
				}
			}

			// ── Per-contact wall ceiling (root-cause fix #2) ──────────────────
			// The raw r^-12 wall energy is unbounded as d->0. A single near-clash
			// can fire +586 (1TT1) or +2578 (1M2Z), dwarfing the complementarity
			// term (com ~ -180..-316) and saturating fitness across the entire GA
			// population -> selection gradient collapses -> clone population ->
			// 0 clustering modes. Bounding each contact's contribution to the
			// fitness keeps the penalty differentiable (one bad contact can at
			// most cancel one favorable com pair) so com-driven selection survives.
			// The wall still penalizes clashes by count; it just no longer explodes.
			// NOTE: the *raw* Ewall (not the capped value) is deliberately used for
			// the DEE-elimination threshold below, so intramolecular-clash pruning
			// is unaffected by the cap.
			// Overlap-based soft-core (v43) or legacy capped r^-12 — shared with
			// Vcontacts pre-filter via soft_wall.h.  DEE/intramolecular checks
			// below still use the raw r^-12 Ewall.
			double Ewall_fitness;
			if (FA->soft_wall_cutoff > 0.0f) {
				Ewall_fitness = soft_wall_fitness_energy(
				    d, cr, FA->soft_wall_cutoff, wal_coercive, wal_stiff);
			} else {
				Ewall_fitness = (Ewall > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall;
			}

			cfs->wal += Ewall_fitness;
			// NOTE: PoseBust clash penalty is NOT computed here. The Vcontacts
			// loop enumerates only Voronoi SURFACE-contact pairs; deep
			// interpenetration pairs are geometrically occluded and never
			// visited, undercounting the true clash ~23x (validated on 1SJ0).
			// The PoseBust term is instead computed by an independent all-pairs
			// cell-list scan once per evaluation (see pb_clash scan after the
			// main loop), matching PoseBusters check_intermolecular_distance.

#if DEBUG_LEVEL > 0
				Ewall_atm += Ewall;
				cfs_atom.wal += Ewall;
#endif
				// ligand intramolecular clash exceeds threshold
				// add an entry in the dee elimination
				if(intramolecular && type == 1 && Ewall > DEE_WALL_THRESHOLD){
					intraclashes.push_back(std::pair<int,int>(atomzero-fatm,atomcont-fatm));
				}
				
				// Treat everything as rigid
				if ( VC->ca_rec[currindex].dist <= dee_clash*rAB ) { cfs->rclash=1; }

			}
			
			
			if( !covalent ){

				if(FA->intramolecular || !intramolecular) {

					double contribution = 0.0;
					if(energy_matrix->weight){
						if(FA->normalize_area){
							contribution = yval*area/surfA;
						}else{
							contribution = yval*area;
						}
					}else{
						contribution = yval;
					}

					if(FA->useacs){
						//printf("USE ACS\n");
						//printf("default contribution=%.3f\n", contribution);
						contribution = contribution * atoms[atomzero].acs/surfA * FA->acsweight;
						//printf("after contribution=%.3f\n", contribution);
					}

					// VCT distance-weighted contacts. Multiply the matrix
					// complementarity by exp(-r/r0) so distal contacts contribute
					// far less than proximal ones — breaking the VCT degeneracy
					// where an off-native pose's distal-arm contacts tie the
					// crystal pose's tight ones. Applied here (before the H-bond /
					// per-type bookkeeping) so every downstream use of contribution
					// sees the weighted value. r_onset = 0 so w(0) = 1.0.
					//
					// Two independent activation paths, both default OFF as a pair
					// and never meant to compound:
					//   (a) env gate FLEXAIDDS_DIST_WEIGHT_CON — r0 from
					//       FLEXAIDDS_CON_R0 (float Å, default 3.5), expf decay.
					//       Same magic-static pattern as the softcore-WAL gate
					//       above (read once, thread-safe, no getenv in the loop).
					//   (b) P9 runner knob FA->vct_dist_weight_r0 > 0 — exp decay
					//       (legacy; config default 4.0 Å).
					// When the env gate is set it takes precedence over (b) so the
					// matrix complementarity is decayed exactly once. r0 <= 0 on
					// the runner knob disables that path (legacy behaviour).
					// dist_weight_con / con_r0 are file-scope statics (hoisted).
					if(dist_weight_con){
						contribution *= expf(-(float)VC->ca_rec[currindex].dist / con_r0);
					}else if(FA->vct_dist_weight_r0 > 0.0){
						double w_r = exp(-VC->ca_rec[currindex].dist /
						                 FA->vct_dist_weight_r0);
						contribution *= w_r;
					}

					// Directional H-bond angular correction (v52: REMOVED).
					// v51 attempted to convert the old sign-error BONUS (contribution<0
					// *(hb_mult-1)<0 * hbond_weight<0 = negative = more attractive for
					// bad geometry) into an explicit frustration penalty via
					// 0.75*lost_attraction. This over-corrected: the VCT grid (0.375Å)
					// means all crystal-pose H-bonds have some angular imperfection, so
					// the penalty accumulated across 20-40 contacts and wiped out the
					// native CF entirely (cf_native(1JD0) collapsed from ~-50 to +0.5).
					// The Gaussian H-bond bonus in compute_hbond_energy() (below) handles
					// directionality correctly and is left intact. This block is a no-op.

					cfs->com += contribution;

					// Lever 2: count this com-contributing contact for the
					// current optimizable residue (intensive normalization below).
					if(FA->vct_normalize_contacts){
						long oi = (long)(atoms[atomzero].optres - FA->optres);
						if(oi >= 0 && oi < (long)FA->num_optres){ vct_ncon[oi]++; }
					}

					// Coulomb electrostatic term (distance-dependent dielectric)
					// Uses RESP charges when available, otherwise standard partial charges.
					{
						double qA = atoms[atomzero].has_resp
						            ? (double)atoms[atomzero].resp_charge
						            : (double)atoms[atomzero].charge;
						double qB = atoms[atomcont].has_resp
						            ? (double)atoms[atomcont].resp_charge
						            : (double)atoms[atomcont].charge;
						if(FA->use_elec && qA != 0.0 && qB != 0.0){
							double dist = VC->ca_rec[currindex].dist;
							if(dist > 0.5){ // avoid singularity
								// E_elec = (332.0637 * qA * qB) / (eps * r)
								// distance-dependent dielectric: eps = dielectric * r
								double E_elec = KCOULOMB * qA * qB / (FA->dielectric * dist * dist);
								cfs->elec += E_elec;
								if(polar_desolv_weight != 0.0){ pd_elec += E_elec; }
							}
						}
					}

					// Angular-dependent hydrogen bond potential (Gaussian bell)
					// v58: search off / rank on — GA skips hbond unless use_hbond_search;
					// post-GA rank re-score sets hbond_rank_rescore with use_hbond_rank.
					if (FA->use_hbond &&
					    ((FA->use_hbond_search && !FA->hbond_rank_rescore) ||
					     (FA->use_hbond_rank && FA->hbond_rank_rescore))) {
						double dist = VC->ca_rec[currindex].dist;
						double E_hb = hbond::compute_hbond_energy(
							atoms, atomzero, atomcont, dist,
							FA->hbond_optimal_dist, FA->hbond_optimal_angle,
							FA->hbond_sigma_dist, FA->hbond_sigma_angle,
							FA->hbond_weight, FA->hbond_salt_bridge_weight);
						cfs->hbond += E_hb;
						if(polar_desolv_weight != 0.0){ pd_hb += E_hb; }
					}

					// Metal ion coordination potential (Gaussian well)
					if (FA->use_metal_coord) {
						double dist = VC->ca_rec[currindex].dist;
						// Charge-aware weight: reduce non-electrostatic component
						// when Coulomb is also active for this pair
						double mc_weight = FA->metal_coord_weight;
						if (FA->use_elec) {
							float mqA = atoms[atomzero].has_resp
							            ? atoms[atomzero].resp_charge
							            : atoms[atomzero].charge;
							float mqB = atoms[atomcont].has_resp
							            ? atoms[atomcont].resp_charge
							            : atoms[atomcont].charge;
							if (mqA != 0.0f && mqB != 0.0f)
								mc_weight *= 0.3;
						}
						double E_mc = metal_coord::compute_metal_coord_energy(
							atoms, atomzero, atomcont, dist,
							mc_weight, FA->metal_coord_sigma);
						cfs->metal_coord += E_mc;
						if(polar_desolv_weight != 0.0){ pd_mc += E_mc; }
						// Track coordination number for CN penalty
						if (E_mc != 0.0) {
							int ta = atoms[atomzero].type;
							int tb = atoms[atomcont].type;
							bool center_is_metal = metal_coord::is_metal_type(ta);
							// Only count if atomzero is the metal (avoid double-counting)
							if (center_is_metal) {
								auto aff = metal_coord::get_donor_affinity(ta, tb);
								if (aff && metal_coord::is_coordinating(dist, aff->ideal_dist))
									metal_cn_count++;
							}
						}
					}

#if DEBUG_LEVEL > 0
					cfs_atom.com += contribution;
#endif

					FA->contributions[(VC->Calc[i].atom->type-1)*FA->ntypes+(VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)] += contribution;
					if((VC->Calc[i].atom->type-1) != (VC->Calc[VC->ca_rec[currindex].atom].atom->type-1))
						FA->contributions[(VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)*FA->ntypes+(VC->Calc[i].atom->type-1)] += contribution;
					
				}
				/*
				  else{
				  printf("skipped intramolecular contact: %d %d\n",
				  atoms[atomzero].number, atoms[atomcont].number); //VC->Calc[VC->ca_rec[currindex].atom].atom->number);
				  } 
				*/
			}
			
			/*
			// generate surface area matrix
			matrix[VC->Calc[i].atom->type][VC->Calc[VC->ca_rec[currindex].atom].atom->type] += area;
			
			if(VC->Calc[i].atom->type != VC->Calc[VC->ca_rec[currindex].atom].atom->type){
			matrix[VC->Calc[VC->ca_rec[currindex].atom].atom->type][VC->Calc[i].atom->type] += area;
			}
			*/
					
#if DEBUG_LEVEL > 0
			printf("        %3s %c %4d %5d %2d %4.2f :: %7.4f (%s) %6.2f %6.2f :: %10.3f %10.3f\n",
			       VC->Calc[VC->ca_rec[currindex].atom].residue->name,
			       VC->Calc[VC->ca_rec[currindex].atom].residue->chn,
			       VC->Calc[VC->ca_rec[currindex].atom].residue->number,
			       VC->Calc[VC->ca_rec[currindex].atom].atom->number,
			       VC->Calc[VC->ca_rec[currindex].atom].atom->type,
			       VC->Calc[VC->ca_rec[currindex].atom].atom->radius,
						
			       yval, energy_matrix->weight ? "Y": "N",
			       VC->ca_rec[currindex].dist,
			       VC->ca_rec[currindex].area,
			       
			       cfs_atom.com, cfs_atom.wal);
			
#endif

			// skip to next contact
			currindex = VC->ca_rec[currindex].prev;
		}
		
		//    printf("Atom[%d]=%d has %d contacts\n",VC->Calc[i].,VC->Calc[i].atom->number,contnum);
		//    printf("Atom[%d] COM=[%8.2f]\tWAL=[%8.2f]\n",VC->Calc[i].atomnum,com_atm,Ewall_atm);

		// Metal coordination number penalty: applied after all pairwise
		// contacts are processed. Penalizes deviations from ideal CN.
		if (FA->use_metal_coord && metal_cn_count > 0) {
			const auto* mp = metal_coord::get_metal_params(atoms[atomzero].type);
			if (mp) {
				double E_cn = metal_coord::cn_penalty(
					metal_cn_count, mp->ideal_cn, FA->metal_coord_cn_weight);
				cfs->metal_coord += E_cn;
			}
		}

		if(SAS < 0.0){ SAS = 0.0; }
		cfs->totsas += SAS;
		
		double contribution = 0.0;
		if(FA->solventterm){
			contribution = (double)FA->solventterm * SAS;
			//printf("SP: multiply ST=%.3f with SAS.area=%.3f\n", (double)FA->solventterm, SAS);
		} else {
			struct energy_matrix* energy_matrix = &FA->energy_matrix[(VC->Calc[i].atom->type-1)*FA->ntypes +
										 (FA->ntypes-1)];
			//printf("type1: %d\ttype2: %d\n", energy_matrix->type1, energy_matrix->type2);
			
			double yval = get_yval(energy_matrix,SAS/surfA);
			
			if(energy_matrix->weight){
				if(FA->normalize_area){
					contribution = yval * SAS / surfA;
				}else{
					contribution = yval * SAS;
				}
				//printf("Weight: multiply yval=%.3f by SAS.area=%.3f\n", yval, SAS);
			}
			else {
				contribution = yval;
				/*
				  if(VC->Calc[i].type == 3){
				  printf("Density: add yval=%.3f for norm.SAS=%.3f for atom %d\n",
				  yval, SAS/surfA, VC->Calc[i].atomnum);
				  }
				*/
			}
		}
		
		if(FA->useacs){
			contribution = contribution * atoms[atomzero].acs/surfA * FA->acsweight;
		}

		// ── Entropy-ablation hook: FLEXAIDDS_NO_SAS ──────────────────────────
		// The SAS channel (residual solvent-accessible surface scored against the
		// solvent pseudo-type via the MC_st0r5.2_6 density function) is the
		// engine's implicit-desolvation / hydration-shell proxy: burying surface
		// on binding pays a per-atom desolvation cost. Setting FLEXAIDDS_NO_SAS
		// zeroes cfs->sas for every atom, ablating the solvation channel so its
		// contribution to pose ranking can be isolated. Read once (magic static,
		// thread-safe) to avoid a getenv() in the inner per-contact loop.
		// no_sas is a file-scope static (hoisted).
		if(no_sas){ contribution = 0.0; }

		cfs->sas += FA->sas_weight * contribution;

		// ── Polar-burial desolvation penalty (FLEXAIDDS_POLAR_DESOLV_WEIGHT) ──
		// Scoped fix for the com-over-burial CF-inversion. A polar LIGAND atom
		// (type==1 optres) that gets buried without a compensating H-bond /
		// salt-bridge / metal contact pays a desolvation cost, added (>=0) to
		// the SAS channel. H-bond-satisfied polar atoms have s_i -> 1 so the
		// penalty vanishes, protecting deep well-satisfied natives. Skipped
		// entirely when the weight is 0 (default) -> zero behaviour change.
		if(polar_desolv_weight != 0.0 && type == 1){
			double Ptype = polar_desolv_ptype(atoms[atomzero].type);
			if(Ptype > 0.0){
				double b_i = 1.0 - (SAS / surfA);          // buried fraction
				if(b_i < 0.0) b_i = 0.0; else if(b_i > 1.0) b_i = 1.0;
				double sat = (std::fabs(pd_hb) + std::fabs(pd_elec) +
				              std::fabs(pd_mc)) / polar_desolv_eref;
				double s_i = (sat < 1.0) ? sat : 1.0;      // satisfaction, [0,1]
				double penalty = polar_desolv_weight * Ptype * b_i * (1.0 - s_i);
				cfs->sas += penalty;                       // positive = cost
			}
		}

		FA->contributions[(VC->Calc[i].atom->type-1)*FA->ntypes + (FA->ntypes-1)] += contribution;
		FA->contributions[(FA->ntypes-1)*FA->ntypes + (VC->Calc[i].atom->type-1)] += contribution;
		
		FA->contacts[VC->Calc[i].atom->number] =
			contacts_epoch_mode ? FA->contacts_epoch : 1;

		// GIST desolvation: accumulate grid-based water displacement energy
		if (FA->use_gist && FA->gist_evaluator != NULL) {
			const auto* grid = static_cast<const gist::GISTGrid*>(FA->gist_evaluator);
			double E_gist = FA->gist_weight *
				grid->desolvation_energy(
					atoms[atomzero].coor[0],
					atoms[atomzero].coor[1],
					atoms[atomzero].coor[2]);
			cfs->gist_desolv += E_gist;
		}

#if DEBUG_LEVEL > 1
		printf("CF.SAS is %.3f for %d contacts with contribution %.3f\n",
		       SAS, contnum, contribution);
#endif

		(void)contnum;
	}


	// Lever 2: convert CF.com from an extensive (sum-of-contacts) score to an
	// intensive one. The raw com rewards sheer contact count, so a ligand jammed
	// into a tight groove can out-score the native H-bond-driven pose purely by
	// burying more surface. Dividing by the contact count scores the *mean*
	// per-contact complementarity instead; the VCT_NREF rescale keeps the term's
	// magnitude comparable to the SAS/wall channels (a bare 1/N collapses com
	// into the noise floor and lets the orientation-independent SAS baseline win).
	if(FA->vct_normalize_contacts){
		constexpr double VCT_NREF = 100.0;
		for(int j=0; j<FA->num_optres; ++j){
			if(vct_ncon[j] > 0){
				FA->optres[j].cf.com *= VCT_NREF / (double)vct_ncon[j];
			}
		}
	}

	// ── P3: soft lower clamp on the CF.com channel (FLEXAIDDS_COM_FLOOR) ────────
	// Enabler, NOT a standalone accuracy fix (see commit message). The favorable
	// (negative) com term is unbounded below: an overpacked non-native pose can
	// drive com → −∞ and swamp every attractive term (H-bond, metal, elec), so a
	// downstream orientation-aware rescorer (P1) can never out-vote it. This
	// installs a *soft* floor at −F: a strictly monotone, bounded-below squashing
	// of com so the runaway tail is capped while pose order by com is preserved.
	//
	// Transition scale S is fixed at 5 (not F). Using F as the softplus width made
	// the floor perturb healthy com values far above −F (e.g. F=130 rewrote com=0
	// to ~+41). With S=5 the map is near-identity for com ≳ −120 and hard-floors
	// the over-burial tail below ~−140 when F=130.
	//
	//   softfloor(x) = −F + S·softplus((x + F)/S),  softplus(z)=max(z,0)+log1p(e^−|z|)
	//     x ≫ −F  ⇒ softfloor(x) → x      (near-identity; ranking untouched)
	//     x → −∞  ⇒ softfloor(x) → −F     (bounded; no term can swamp the sum)
	//     softfloor′ ∈ (0,1]              (monotone ⇒ rank-preserving)
	//
	// Env-gated, DEFAULT-OFF: unset or F≤0 ⇒ skipped entirely ⇒ bit-identical.
	// FLEXAIDDS_COM_FLOOR=F sets the floor magnitude (kcal/mol-equivalent CF units).
	if(const char* com_floor_env = std::getenv("FLEXAIDDS_COM_FLOOR")){
		const double F = std::atof(com_floor_env);
		if(F > 0.0){
			// Floor at −F; transition scale S << F so near-identity for com >> −F.
			// With S=5: near-identity for com > −120, hard floor for com < −140.
			const double S = 5.0;
			for(int j=0; j<FA->num_optres; ++j){
				double& com = FA->optres[j].cf.com;
				const double z  = (com + F) / S;
				const double az = (z < 0.0) ? -z : z;
				const double softplus = (z > 0.0 ? z : 0.0) + std::log1p(std::exp(-az));
				com = -F + S * softplus;
			}
		}
	}

	// Shannon entropy of contact-type distribution as VCT false-minimum penalty.
	// H = -Σ p_ij log₂(p_ij) over pairwise atom-type contribution matrix.
	// False minima have HIGH entropy (many diffuse type-pair contacts);
	// native sites have LOW entropy (few dominant, specific complementary pairs).
	// Penalty: CF += λ·H demotes false minima relative to native binding mode.
	// contributions[] is distance-weighted (exp(-r/r0) applied in contact loop).
	// Only negative contributions (favorable contacts) define the distribution.
	if (FA->vct_entropy_weight > 0.0) {
		double total_abs = 0.0;
		const int ntypes2 = FA->ntypes * FA->ntypes;
		for (int k = 0; k < ntypes2; ++k) {
			if (FA->contributions[k] < 0.0f)
				total_abs += (double)(-FA->contributions[k]);
		}
		double H = 0.0;
		if (total_abs > 0.0) {
			for (int k = 0; k < ntypes2; ++k) {
				if (FA->contributions[k] < 0.0f) {
					double p = (double)(-FA->contributions[k]) / total_abs;
					if (p > 0.0) H -= p * std::log2(p);
				}
			}
		}
		double entropy_penalty = FA->vct_entropy_weight * H;
		for (int j = 0; j < FA->num_optres; ++j) {
			if (FA->optres[j].type == 1) {  // ligand optres only
				FA->optres[j].cf.entropy += entropy_penalty;
				break;
			}
		}
	}

	// Penalize Freesurf.

	//printf("FREESURF(SAS)=[%8.2f]\n",SAStot);
	//printf("FINAL SUM COM=[%8.2f]\tWAL=[%8.2f]\n",com,Ewall);
  
	//print_surfmat(matrix,"surf.mat");


/*
  #if DEBUG_LEVEL > 0
  printf("\n");
  printf("CF.sum = %.3f\n", cfs->com + cfs->sas + cfs->wal);
  printf("CF.com = %.3f\n", cfs->com);
  printf("CF.sas = %.3f\n", cfs->sas);
  printf("CF.wal = %.3f\n", cfs->wal);
  getchar(); 
  #endif
*/
	//getchar();

	// GIST water displacement scoring (applied once per evaluation)
	if(FA->use_gist && FA->gist_evaluator != NULL){
		const GISTEvaluator* gist =
			static_cast<const GISTEvaluator*>(FA->gist_evaluator);
		double gist_score = gist->score_ligand(atoms, FA);
		// Distribute GIST score to first ligand optres
		for(int j=0; j<FA->num_optres; ++j){
			if(FA->optres[j].type == 1){
				FA->optres[j].cf.gist += gist_score;
				break;
			}
		}
	}

	// ── PoseBust physical-realism clash penalty — all-pairs cell-list scan ──
	// Computed ONCE per pose evaluation (like GIST above), NOT inside the
	// Vcontacts contact loop. The Vcontacts loop enumerates only Voronoi
	// surface-contact pairs, so deep interpenetration (a ligand atom buried
	// inside receptor atoms) is occluded and undercounts the clash ~23x
	// (validated on 1SJ0: engine raw 3.01 vs all-pairs 70.79). This scan
	// visits EVERY ligand-atom / receptor-atom pair within the max PB vdW-sum
	// cutoff via a uniform spatial grid (cell-list), matching PoseBusters'
	// check_intermolecular_distance semantics: a pair clashes when the
	// interatomic distance d < pb_clash_ratio*(vdw_i+vdw_j) (element vdW radii
	// from posebusters_vdw_radius(), soft_wall.h). Penalty is uncapped and
	// severity-scaled: weight * sum(max(0, cr_pb - d)^p), so a physically
	// impossible interpenetration always out-weighs the unbounded CF.com.
	// OFF unless pb_clash_weight>0.
	//
	// The same receptor cell list also backs the pocket-presence penalty (P1), so
	// the block runs when either term is enabled and each term is gated
	// individually inside. OFF unless pb_clash_weight>0 or pb_pocket_weight>0.
	if (FA->pb_clash_weight > 0.0 || flexaids::pb_pocket_enabled(FA->pb_pocket_weight)) {
		// Partition atoms into rigid-receptor vs movable-ligand (optres != NULL).
		// Ligand atoms are few (tens); receptor atoms are many (thousands) — so
		// build the grid over RECEPTOR atoms and scan each ligand atom's cell nbhd.
		const double PB_MAX_VDW = 2.10;               // max element vdW (I) in soft_wall.h table
		// cell and inv_cell are derived from pb_clash_ratio, which is set once at dock init.
		// They are the same across all evals in a dock session.

		// ── Route A: hoist receptor grid construction out of the eval loop ────
		// The receptor cell-list grid is built once per dock (keyed on atm_cnt +
		// pb_clash_ratio) and reused across all evals.  The CF math is identical;
		// only the grid *construction* moves out of the eval loop.  This is
		// unconditional whenever pb_clash_weight>0 — see the cache comment above.
		if (!pb_cache.matches(FA->atm_cnt, FA->pb_clash_ratio)) {
			// Invalidate and rebuild.
			pb_cache.valid   = false;
			pb_cache.atm_cnt = FA->atm_cnt;
			pb_cache.ratio   = FA->pb_clash_ratio;
			pb_cache.cell    = FA->pb_clash_ratio * 2.0 * PB_MAX_VDW;
			pb_cache.inv_cell= 1.0 / pb_cache.cell;

			pb_cache.rec_idx.clear();
			pb_cache.rec_vdw.clear();
			pb_cache.rec_metal.clear();
			pb_cache.rec_heavy.clear();
			double rmax[3] = {-1e30,-1e30,-1e30};
			pb_cache.rmin[0]=pb_cache.rmin[1]=pb_cache.rmin[2]=1e30;

			for (int ai = 1; ai <= FA->atm_cnt; ++ai) {
				if (atoms[ai].optres == NULL) {
					pb_cache.rec_idx.push_back(ai);
					for (int k=0;k<3;++k){
						double c=atoms[ai].coor[k];
						if(c<pb_cache.rmin[k])pb_cache.rmin[k]=c;
						if(c>rmax[k])rmax[k]=c;
					}
					// Precompute receptor vdW radius (avoids repeated get_element calls in inner loop).
					// Both sides of a pair must draw from the same source, so the cached-field
					// switch applies here too — see the pb_vdw_cached comment above.
					const char* re_raw = get_element(atoms[ai].type);
					while (*re_raw == ' ') ++re_raw;
					pb_cache.rec_vdw.push_back(
						pb_vdw_cached ? atoms[ai].pb_vdw_radius
						              : posebusters_vdw_radius(re_raw, atoms[ai].radius));
					pb_cache.rec_metal.push_back(
						pb_metal_carveout && posebusters_is_coordinating_metal(re_raw) ? 1 : 0);
					// Heavy = not hydrogen. Type 39 ("Du") is FlexAID's dummy/hydrogen
					// tag; "H" covers readers that type hydrogens explicitly.
					pb_cache.rec_heavy.push_back(
						(std::strcmp(re_raw, "H") == 0 || std::strcmp(re_raw, "Du") == 0) ? 0 : 1);
				}
			}
			for (int k=0;k<3;++k){
				pb_cache.dim[k] = std::max(1, (int)((rmax[k]-pb_cache.rmin[k])*pb_cache.inv_cell) + 1);
			}
			pb_cache.grid.clear();
			pb_cache.grid.reserve(pb_cache.rec_idx.size());
			auto cidx_build = [&](int cx,int cy,int cz){
				return (cx*pb_cache.dim[1]+cy)*pb_cache.dim[2]+cz;
			};
			for (int i = 0; i < (int)pb_cache.rec_idx.size(); ++i) {
				int ri = pb_cache.rec_idx[i];
				int cx=(int)((atoms[ri].coor[0]-pb_cache.rmin[0])*pb_cache.inv_cell);
				int cy=(int)((atoms[ri].coor[1]-pb_cache.rmin[1])*pb_cache.inv_cell);
				int cz=(int)((atoms[ri].coor[2]-pb_cache.rmin[2])*pb_cache.inv_cell);
				if(cx<0)cx=0; if(cx>=pb_cache.dim[0])cx=pb_cache.dim[0]-1;
				if(cy<0)cy=0; if(cy>=pb_cache.dim[1])cy=pb_cache.dim[1]-1;
				if(cz<0)cz=0; if(cz>=pb_cache.dim[2])cz=pb_cache.dim[2]-1;
				pb_cache.grid[cidx_build(cx,cy,cz)].push_back(i); // store position in rec_idx, not atom index
			}
			pb_cache.valid = true;
		}

		// Determine whether we're using the cached path or the legacy per-eval path.
		const bool use_cache = pb_cache.valid;

		// Local grid variables: either point at the cache or build fresh (legacy).
		double l_rmin[3]  = { 1e30, 1e30, 1e30};
		double l_rmax[3]  = {-1e30,-1e30,-1e30};
		int    l_dim[3]   = {};
		double l_cell     = FA->pb_clash_ratio * 2.0 * PB_MAX_VDW;
		double l_inv_cell = 1.0 / l_cell;
		std::vector<int> lig_idx;
		lig_idx.reserve(64);
		// Legacy: also collect rec_idx and build grid locally when not using cache.
		std::vector<int>    l_rec_idx;
		std::unordered_map<int,std::vector<int>> l_grid;

		if (!use_cache) {
			// Legacy path: collect receptor atoms, bounding box, and grid (same as before).
			l_rec_idx.reserve(4096);
			for (int ai = 1; ai <= FA->atm_cnt; ++ai) {
				if (atoms[ai].optres != NULL) {
					lig_idx.push_back(ai);
				} else {
					l_rec_idx.push_back(ai);
					for (int k=0;k<3;++k){ double c=atoms[ai].coor[k]; if(c<l_rmin[k])l_rmin[k]=c; if(c>l_rmax[k])l_rmax[k]=c; }
				}
			}
			for (int k=0;k<3;++k){ l_dim[k] = std::max(1, (int)((l_rmax[k]-l_rmin[k])*l_inv_cell) + 1); }
			auto cidx = [&](int cx,int cy,int cz){ return (cx*l_dim[1]+cy)*l_dim[2]+cz; };
			l_grid.reserve(l_rec_idx.size());
			for (int ri : l_rec_idx) {
				int cx=(int)((atoms[ri].coor[0]-l_rmin[0])*l_inv_cell);
				int cy=(int)((atoms[ri].coor[1]-l_rmin[1])*l_inv_cell);
				int cz=(int)((atoms[ri].coor[2]-l_rmin[2])*l_inv_cell);
				if(cx<0)cx=0; if(cx>=l_dim[0])cx=l_dim[0]-1;
				if(cy<0)cy=0; if(cy>=l_dim[1])cy=l_dim[1]-1;
				if(cz<0)cz=0; if(cz>=l_dim[2])cz=l_dim[2]-1;
				l_grid[cidx(cx,cy,cz)].push_back(ri);
			}
		} else {
			// Cached path: only need to collect ligand atoms.
			for (int ai = 1; ai <= FA->atm_cnt; ++ai) {
				if (atoms[ai].optres != NULL) lig_idx.push_back(ai);
			}
		}

		// Select which grid/metadata to use for the scan.
		const double* rmin       = use_cache ? pb_cache.rmin      : l_rmin;
		const int*    dim        = use_cache ? pb_cache.dim        : l_dim;
		const double  inv_cell_s = use_cache ? pb_cache.inv_cell   : l_inv_cell;
		const std::unordered_map<int,std::vector<int>>& grid_ref =
			use_cache ? pb_cache.grid : l_grid;

		if (!lig_idx.empty() && (use_cache ? !pb_cache.rec_idx.empty() : !l_rec_idx.empty())) {
			auto cidx = [&](int cx,int cy,int cz){ return (cx*dim[1]+cy)*dim[2]+cz; };
			if (FA->pb_clash_weight > 0.0) {
			double pb_pen = 0.0;
			long pb_nclash = 0;
			long pb_ncarved = 0;
			for (int li : lig_idx) {
				const double lx=atoms[li].coor[0], ly=atoms[li].coor[1], lz=atoms[li].coor[2];
				// The element string is only needed for the table lookup and the metal
				// test; when the cached radius is in use and the carve-out is off,
				// skip the ~26-branch string compare entirely.
				const char* le_raw = "";
				if (!pb_vdw_cached || pb_metal_carveout) {
					le_raw = get_element(atoms[li].type);
					while (*le_raw == ' ') ++le_raw;
				}
				const double vdw_l = pb_vdw_cached
					? atoms[li].pb_vdw_radius
					: posebusters_vdw_radius(le_raw, atoms[li].radius);
				// Metal carve-out, ligand side: when the ligand atom is itself a
				// coordinating metal every pair it forms is carved out.
				const bool lig_is_metal =
					pb_metal_carveout && posebusters_is_coordinating_metal(le_raw);
				int cx=(int)((lx-rmin[0])*inv_cell_s);
				int cy=(int)((ly-rmin[1])*inv_cell_s);
				int cz=(int)((lz-rmin[2])*inv_cell_s);
				for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz){
					int nx=cx+dx,ny=cy+dy,nz=cz+dz;
					if(nx<0||ny<0||nz<0||nx>=dim[0]||ny>=dim[1]||nz>=dim[2]) continue;
					auto it=grid_ref.find(cidx(nx,ny,nz));
					if(it==grid_ref.end()) continue;
					for(int pos : it->second){
						double vdw_r;
						int ri;
						bool rec_is_metal;
						if (use_cache) {
							// pos is an index into rec_idx / rec_vdw / rec_metal
							ri    = pb_cache.rec_idx[pos];
							vdw_r = pb_cache.rec_vdw[pos];
							rec_is_metal = (pb_cache.rec_metal[pos] != 0);
						} else {
							// pos is the atom index directly (legacy path)
							ri = pos;
							const char* re_raw = get_element(atoms[ri].type);
							while (*re_raw == ' ') ++re_raw;
							vdw_r = pb_vdw_cached
								? atoms[ri].pb_vdw_radius
								: posebusters_vdw_radius(re_raw, atoms[ri].radius);
							rec_is_metal =
								pb_metal_carveout && posebusters_is_coordinating_metal(re_raw);
						}
						// Metal carve-out: skip the pair if EITHER partner is a
						// coordinating metal. No-op when the carve-out is off, since
						// both flags are then hard-false.
						if (lig_is_metal || rec_is_metal) { ++pb_ncarved; continue; }
						const double ddx=lx-atoms[ri].coor[0];
						const double ddy=ly-atoms[ri].coor[1];
						const double ddz=lz-atoms[ri].coor[2];
						const double d=std::sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
						if(d < 1.0e-6) continue;
						const double cr_pb=FA->pb_clash_ratio*(vdw_l+vdw_r);
						const double o=cr_pb-d;
						if(o>0.0){ pb_pen += std::pow(o, FA->pb_clash_exponent); ++pb_nclash; }
					}
				}
			}
			if (std::getenv("FLEXAIDDS_PB_CLASH_DEBUG")) {
				fprintf(stderr, "[PB_CLASH_DEBUG] lig_atoms=%zu rec_atoms=%zu nclash=%ld ncarved=%ld raw_pen=%.4f cached=%d metal_carveout=%d\n",
					lig_idx.size(),
					use_cache ? pb_cache.rec_idx.size() : l_rec_idx.size(),
					pb_nclash, pb_ncarved, pb_pen, (int)use_cache,
					(int)pb_metal_carveout);
			}
			pb_pen *= FA->pb_clash_weight;
			// Distribute to the first ligand optres (same convention as GIST).
			for(int j=0;j<FA->num_optres;++j){
				if(FA->optres[j].type == 1){ FA->optres[j].cf.pb_clash += pb_pen; break; }
			}
			}  // end pb_clash_weight > 0.0

			// ── P1: pocket-presence soft penalty ──────────────────────────────
			// Catches the "ligand left the pocket" false minimum — the dominant
			// Astex failure class here (1SG0 elected a pose whose centroid sat
			// 29.64 A from the true binding site). PoseBusters reports this
			// post-dock as protein-ligand_maximum_distance; steering on it during
			// the GA search is what stops the population settling there at all.
			//
			// Measure: distance from the ligand centroid to the nearest receptor
			// HEAVY atom. Penalty: a soft quadratic ramp that is exactly zero
			// while the centroid is within pb_pocket_radius, so shallow and
			// surface binders are untouched and only real drift is charged.
			//
			// The nearest-atom search walks the already-hoisted cell list outward
			// in Chebyshev shells and stops as soon as no unscanned cell can hold
			// anything closer ((R-1)*cell >= d_best). For a ligand sitting in the
			// pocket that terminates at shell 0 or 1 — O(1). A drifted pose costs
			// more shells, but those are exactly the poses being penalized and are
			// rare once the term is steering the population.
			if (flexaids::pb_pocket_enabled(FA->pb_pocket_weight)) {
				double sx = 0.0, sy = 0.0, sz = 0.0;
				for (int li : lig_idx) {
					sx += atoms[li].coor[0];
					sy += atoms[li].coor[1];
					sz += atoms[li].coor[2];
				}
				const double n_lig = (double)lig_idx.size();
				const double gx = sx / n_lig, gy = sy / n_lig, gz = sz / n_lig;

				const int ccx = (int)((gx-rmin[0])*inv_cell_s);
				const int ccy = (int)((gy-rmin[1])*inv_cell_s);
				const int ccz = (int)((gz-rmin[2])*inv_cell_s);

				// Enough shells to reach every cell from any centroid, in or out of box.
				const int max_shell = dim[0] + dim[1] + dim[2] +
					std::abs(ccx) + std::abs(ccy) + std::abs(ccz) + 2;
				const double cell_s = use_cache ? pb_cache.cell : l_cell;

				double best_d2 = 1e30;
				for (int R = 0; R <= max_shell; ++R) {
					// No cell at Chebyshev radius R or beyond can beat best_d2 once
					// the inner boundary of this shell is already further away.
					if (best_d2 < 1e29) {
						const double reach = (double)(R-1) * cell_s;
						if (reach > 0.0 && reach*reach >= best_d2) break;
					}
					for(int dx=-R;dx<=R;++dx)for(int dy=-R;dy<=R;++dy)for(int dz=-R;dz<=R;++dz){
						// Shell only — the interior was scanned at radius R-1.
						if (R > 0 && std::abs(dx) != R && std::abs(dy) != R && std::abs(dz) != R)
							continue;
						const int nx=ccx+dx, ny=ccy+dy, nz=ccz+dz;
						if(nx<0||ny<0||nz<0||nx>=dim[0]||ny>=dim[1]||nz>=dim[2]) continue;
						auto it=grid_ref.find(cidx(nx,ny,nz));
						if(it==grid_ref.end()) continue;
						for(int pos : it->second){
							int ri;
							bool heavy;
							if (use_cache) {
								ri    = pb_cache.rec_idx[pos];
								heavy = (pb_cache.rec_heavy[pos] != 0);
							} else {
								ri = pos;
								const char* re_raw = get_element(atoms[ri].type);
								while (*re_raw == ' ') ++re_raw;
								heavy = !(std::strcmp(re_raw,"H")==0 || std::strcmp(re_raw,"Du")==0);
							}
							if (!heavy) continue;
							const double ddx=gx-atoms[ri].coor[0];
							const double ddy=gy-atoms[ri].coor[1];
							const double ddz=gz-atoms[ri].coor[2];
							const double d2=ddx*ddx+ddy*ddy+ddz*ddz;
							if(d2 < best_d2) best_d2 = d2;
						}
					}
				}

				double pocket_pen = 0.0;
				double centroid_d = 0.0;
				if (best_d2 < 1e29) {
					centroid_d = std::sqrt(best_d2);
					const double over = centroid_d - FA->pb_pocket_radius;
					if (over > 0.0) pocket_pen = FA->pb_pocket_weight * over * over;
				}

				if (pb_pocket_debug) {
					fprintf(stderr, "[PB_POCKET] centroid_dist=%.2f pocket_radius=%.2f penalty=%.2f\n",
						centroid_d, FA->pb_pocket_radius, pocket_pen);
				}

				if (pocket_pen != 0.0) {
					// Same attribution convention as pb_clash / GIST: first ligand optres.
					for(int j=0;j<FA->num_optres;++j){
						if(FA->optres[j].type == 1){ FA->optres[j].cf.pb_clash += pocket_pen; break; }
					}
				}
			}
		}
	}

	for(int i=0;i<FA->atm_cnt_real;i++){
		if(VC->Calc[i].score){
			VC->Calc[i].atom = NULL;
		}
	}

	if(!FA->vindex){ free(VC->box); }
	
	return(0.0);
  
}

/* A2 perf: flat-array piecewise-linear interpolation replaces linked-list walk.
   Boundary logic matches the original:
     ra < fx[0]            → 0 (no left-bound data)
     ra in [fx[i], fx[i+1]) → linear interpolation using precomputed slope
     ra >= fx[n-1]         → fy[n-1] (no right-bound data; clamp to last y)
   Linear scan over n (typically 5–15 breakpoints) is cache-friendlier than
   branching binary search at this scale. */
double get_yval(struct energy_matrix* em, double relative_area)
{
	if(!em->energy_values) return 0.0;
	// single-scalar weight case (no linked-list walk needed either way)
	if(em->weight) return (double)em->energy_values->y;
	// flat-array path
	const int n = em->flat_n;
	if(n == 0) return 0.0;
	const float ra = (float)relative_area;
	const float* fx = em->flat_x;
	const float* fy = em->flat_y;
	if(ra < fx[0]) return 0.0;                   // below first breakpoint
	if(ra >= fx[n-1]) return (double)fy[n-1];    // at or beyond last breakpoint
	int i = 0;
	while(i < n-2 && ra >= fx[i+1]) ++i;         // find segment (n small → linear scan)
	return (double)(fy[i] + em->flat_slope[i] * (ra - fx[i]));
}
