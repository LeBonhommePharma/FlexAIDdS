#include "Vcontacts.h"
#include "GISTEvaluator.h"
#include "HBondEvaluator.h"
#include "hbond_potential.h"
#include "metal_coordination.h"
#include "GISTGrid.h"
#include <cmath>
#include <cstdlib>
#include <vector>

#define DEBUG_LEVEL 0

// Coulomb constant: 332.0637 kcal·Å/(mol·e²)
#define KCOULOMB 332.0637

double vcfunction(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue, std::vector<std::pair<int,int> > & intraclashes, bool* error)
{
	int    rnum=0;
	int    type=1;
	
	// reset all values pointed
	memset(FA->contacts,0,MAX_ATOM_NUMBER*sizeof(int));
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
			
			if(FA->contacts[VC->Calc[VC->ca_rec[currindex].atom].atom->number]){
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
			// Read once (magic static, thread-safe) — no getenv in the hot loop.
			static const bool softcore_wal =
				(std::getenv("FLEXAIDDS_SOFTCORE_WAL") != nullptr);
			if(softcore_wal){
				const double r_softcrit = 0.7 * cr;
				// ── Hard-floor gate (FLEXAIDDS_SOFTCORE_FLOOR, default 0.5) ────────
				// The parabola plateau can sit below WAL_CONTACT_CAP=50 for deeply
				// buried large-atom pairs, letting catastrophically buried poses score
				// nearly free.  For d < r_hardfloor we revert to the raw r^-12 wall
				// so the existing per-contact cap of 50 applies.  The parabola only
				// operates in the moderate-penetration band [r_hardfloor, r_softcrit).
				static const double softcore_floor_frac = [](){
					const char* s = std::getenv("FLEXAIDDS_SOFTCORE_FLOOR");
					double v = 0.5;   // default: hard wall below 50 % of contact radius
					if(s){ double p = strtod(s, nullptr); if(p > 0.0 && p < 0.7) v = p; }
					return v;
				}();
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
			constexpr double WAL_CONTACT_CAP = 50.0;
			double Ewall_fitness = (Ewall > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall;

			cfs->wal += Ewall_fitness;

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
					static const bool dist_weight_con =
						(std::getenv("FLEXAIDDS_DIST_WEIGHT_CON") != nullptr);
					static const float con_r0 = []() {
						const char* s = std::getenv("FLEXAIDDS_CON_R0");
						float v = 3.5f;
						if(s){ float p = strtof(s, nullptr); if(p > 0.0f){ v = p; } }
						return v;
					}();
					if(dist_weight_con){
						contribution *= expf(-(float)VC->ca_rec[currindex].dist / con_r0);
					}else if(FA->vct_dist_weight_r0 > 0.0){
						double w_r = exp(-VC->ca_rec[currindex].dist /
						                 FA->vct_dist_weight_r0);
						contribution *= w_r;
					}

					// Directional H-bond angular correction:
					// Scale complementarity by angular multiplier for H-bond pairs
					if(FA->use_hbond){
						double hb_mult = hbond::evaluate_contact(
							&atoms[atomzero], &atoms[atomcont],
							atoms, VC->ca_rec[currindex].dist);
						if(hb_mult < 1.0){
							double hb_correction = contribution * (hb_mult - 1.0) * FA->hbond_weight;
							cfs->hbond += hb_correction;
						}
					}

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
							}
						}
					}

					// Angular-dependent hydrogen bond potential (Gaussian bell)
					if (FA->use_hbond) {
						double dist = VC->ca_rec[currindex].dist;
						double E_hb = hbond::compute_hbond_energy(
							atoms, atomzero, atomcont, dist,
							FA->hbond_optimal_dist, FA->hbond_optimal_angle,
							FA->hbond_sigma_dist, FA->hbond_sigma_angle,
							FA->hbond_weight, FA->hbond_salt_bridge_weight);
						cfs->hbond += E_hb;
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
		static const bool no_sas = (std::getenv("FLEXAIDDS_NO_SAS") != nullptr);
		if(no_sas){ contribution = 0.0; }

		cfs->sas += contribution;

		FA->contributions[(VC->Calc[i].atom->type-1)*FA->ntypes + (FA->ntypes-1)] += contribution;
		FA->contributions[(FA->ntypes-1)*FA->ntypes + (VC->Calc[i].atom->type-1)] += contribution;
		
		FA->contacts[VC->Calc[i].atom->number] = 1;

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

	for(int i=0;i<FA->atm_cnt_real;i++){
		if(VC->Calc[i].score){
			VC->Calc[i].atom = NULL;
		}
	}

	if(!FA->vindex){ free(VC->box); }
	
	return(0.0);
  
}

double get_yval(struct energy_matrix* energy_matrix, double relative_area)
{
	double yval = 0.0;
	if(energy_matrix->energy_values == NULL) return 0.0;
	
	// a single value in matrix (weighted by area)
	if(energy_matrix->weight)
		yval = energy_matrix->energy_values->y;
	else { // density function
		struct energy_values* xyval = energy_matrix->energy_values;

		while(xyval->next_value != NULL && relative_area > xyval->next_value->x){
			/*
			  printf("x=%.3f next_value.x=%.3f next_value.y=%.3f\n",
			  xyval->x, xyval->next_value->x, xyval->next_value->y);
			*/
			xyval = xyval->next_value;
		}
		
		if(xyval->x > relative_area){
			// no left bound data
			yval = 0.0;
		}else if(xyval->next_value == NULL){
			// no right bound data
			yval = xyval->y;
		}else{
			yval = xyval->y + 
				( relative_area - xyval->x ) / (xyval->next_value->x - xyval->x ) *
				( xyval->next_value->y - xyval->y );
		}
		
		/*
		  if(energy_matrix->type2 == 40){
		  printf("stopped at x=%.3f with y=%.3f\n", xyval->x, xyval->y);
		  if(xyval->next_value != NULL){
		  printf("next is x=%.3f with y=%.3f\n", xyval->next_value->x, xyval->next_value->y);
		  }
		  printf("prob func. yval=%.3f for relative_area %.3f for [%d][%d]\n", yval, relative_area,
		  energy_matrix->type1, energy_matrix->type2);
		  printf("calculated y=%.3f\n", yval);
		  getchar();
		  }
		*/
	}

	return yval;
}
