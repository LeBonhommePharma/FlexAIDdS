#include "flexaid.h"
#include "fileio.h"

#define N_NEW_PAR  5

/***************************************************************************** 
 * SUBROUTINE add2_optimiz_vec builds the vector with the values of the ic's
 * that are going to be optimized for one or more residues/ligands.
 *****************************************************************************/
void add2_optimiz_vec(FA_Global* FA,atom* atoms,resid* residue,int val[], char chain, const char* extras){
	int i;
	int at;
	//int rot;
	
	// val[0]=residue number
	// val[1]=0 -> optimization of global position; 
	// val[1]=n -> optimization of dihedral bond number n
		
	i=1;
	while(residue[i].number != val[0] || residue[i].chn != chain){
		i++;
		if(i == FA->res_cnt) break;
	}
	at=i;
  
	buildic(FA,atoms,residue,at);

// PRODUCTION GUARD 1 of 2: RESERVE BEFORE ANY POINTER IS TAKEN
//
// Every branch below stores raw pointers into the reallocatable arrays:
//     atoms[...].par              = &FA->map_par[FA->npar];
//     FA->map_par_sidechain_first = &FA->map_par[FA->npar];
//     FA->map_par_sidechain_last  = &FA->map_par[FA->npar];
// realloc_par() relocates map_par, so a realloc anywhere in this call leaves
// every pointer taken earlier in it out of range -- including the ligand
// branches, when the first side-chain gene triggers the growth. Reserving up
// front makes the arrays immovable for the rest of the call; it is a no-op when
// capacity already suffices, so correct paths are unchanged.
//
// Defence in depth only. The reservation that MATTERS is the global one in
// top.cpp / read_input.cpp before the FIRST add2_optimiz_vec call.
//
// MECHANISM CORRECTED 2026-09-05 (external code audit). What is measured and what
// is not:
//   MEASURED  add2_optimiz_vec() stores raw pointers into FA->map_par and
//             realloc_par() relocates it, so a growth mid-setup leaves earlier
//             captures outside the live array. validate_map_par_pointers()
//             reports this as [PAR-PTR]; 11 of 18 on 1R55, also 4 of 15 and
//             3 of 16, across 21 logs on disk.
//   MEASURED  exit 139 on 1R55 / 1R58 / 1HQ2 with as few as one flexible side
//             chain at population >= 400; population 100 and the rigid arm clean.
//   REFUTED   that the out-of-range map_par pointer is what faults. GREPPED AT
//             HEAD: atoms[].par is assigned 8 times and DEREFERENCED ZERO times
//             in LIB; map_par_sidechain_first/last likewise, their only
//             non-assignment read being inside validate_map_par_pointers itself.
//             All three are WRITE-ONLY, so a stale one cannot fault.
//   ACTUAL    the crash was ASan-localised to atoms[].optres, dereferenced at
//             vcfunction.cpp:796-798 and GISTEvaluator.cpp:206, against an array
//             relocated by build_rotamers.cpp:308 (see vcfunction.cpp:748).
//             reserve_optres() is the guard for that array.
// So this reservation buys uniform map_par lifetime and a clean detector, not
// crash safety. Full discussion at top.cpp:2661.
	// MEASURED: exit 139 with as few as ONE flexible side chain at population
	// >= 400 (1R55/1R58/1HQ2); population 100 and the rigid arm clean.
	// Reserving up front makes the arrays immovable for the rest of the call.
	// No-op when capacity already suffices, so correct paths are unchanged.
	// Defence in depth only. The reservation that MATTERS is the global one in
	// top.cpp / read_input.cpp before the FIRST add2_optimiz_vec call: this
	// function is called up to four times per run (extras "" x3, then "SC",
	// "NM"), and growing the array on a LATER call invalidates the pointers the
	// EARLIER calls captured at lines 128-207. Measured: reserving here first
	// grew MIN_PAR 6 -> 86 on the "SC" call and left 11 of 18 ligand pointers
	// dangling. With the global reservation in place this is a no-op.
	reserve_par(FA, FA->npar + FA->nflxsc + FA->normal_modes + 8);
	const optmap* par_base_on_entry = FA->map_par;

	if(strcmp(extras,"SC") == 0){

		for(i=0;i<FA->nflxsc;i++){
            
			if(residue[FA->flex_res[i].inum].trot > 0){
				
				//printf("new par flex sc: %s %c %d\n", residue[FA->flex_res[i].inum].name,residue[FA->flex_res[i].inum].chn, residue[FA->flex_res[i].inum].number);
				//printf("npar: %d - MIN_PAR: %d\n", FA->npar,FA->MIN_PAR);
                
				if(FA->npar==FA->MIN_PAR){ realloc_par(FA,&FA->MIN_PAR); }
				
				FA->map_par[FA->npar].atm = residue[FA->flex_res[i].inum].fatm[0];
				FA->map_par[FA->npar].typ = 4;
				FA->opt_par[FA->npar] = 0.0;
				FA->del_opt_par[FA->npar] = FA->delta_index;
				FA->min_opt_par[FA->npar] = 0.0;
				FA->max_opt_par[FA->npar] = (double)residue[FA->flex_res[i].inum].trot;
				FA->map_opt_par[FA->npar] = 1;

				atoms[residue[FA->flex_res[i].inum].fatm[0]].par = &FA->map_par[FA->npar];
				
				printf("npar=%d map_par[%d].typ=%d map_par[%d].atm=%d opt_par[%d]=%f\n",FA->npar,	   
				       FA->npar,FA->map_par[FA->npar].typ,
				       FA->npar,atoms[FA->map_par[FA->npar].atm].number,
				       FA->npar,FA->opt_par[FA->npar]);
				//	  printf("min=%f max=%f del=%f\n",FA->min_opt_par[FA->npar],FA->max_opt_par[FA->npar],FA->del_opt_par[FA->npar]);

				if(FA->map_par_sidechain_first == NULL){
					FA->map_par_sidechain_first = &FA->map_par[FA->npar];
					FA->map_par_sidechain_first_index = FA->npar;
				}
				FA->map_par_sidechain_last = &FA->map_par[FA->npar];
				
				FA->npar++;
                
			}
            
		}
        
		// PRODUCTION GUARD 2: prove the reservation held and every captured pointer
		// still lands inside the live array. Loud failure, never silent corruption.
		if(FA->map_par != par_base_on_entry){
			fprintf(stderr,"[PAR-PTR] add2_optimiz_vec(SC): map_par RELOCATED during "
			        "the call despite reservation (npar=%d MIN_PAR=%d)\n",
			        FA->npar, FA->MIN_PAR);
			const char* rf = getenv("FLEXAIDDS_PAR_PTR_FATAL");
			if(rf && rf[0] == '1') Terminate(2);
		}
		validate_map_par_pointers(FA, atoms, FA->atm_cnt, "add2_optimiz_vec(SC)");

	}else if(strcmp(extras,"NM") == 0){
                
		if (FA->normal_modes > 0){
            
			if(FA->npar==FA->MIN_PAR){ realloc_par(FA,&FA->MIN_PAR); }
            
			FA->map_par[FA->npar].typ = 3;
			FA->opt_par[FA->npar] = 0.0;
			FA->del_opt_par[FA->npar] = FA->delta_index;
			FA->min_opt_par[FA->npar] = FA->normalindex_min;
			FA->max_opt_par[FA->npar] = FA->normalindex_max;
			FA->map_opt_par[FA->npar] = 1;
            
			printf("npar=%d map_par[%d].typ=%d map_par[%d].atm=%d opt_par[%d]=%f\n",FA->npar,	   
			       FA->npar,FA->map_par[FA->npar].typ,
			       FA->npar,atoms[FA->map_par[FA->npar].atm].number,
			       FA->npar,FA->opt_par[FA->npar]);
			//      printf("min=%f max=%f del=%f\n",FA->min_opt_par[FA->npar],FA->max_opt_par[FA->npar],FA->del_opt_par[FA->npar]);
            
			FA->npar++;
		}
        
	}else if(val[1] == -1){    // (1 degree of freedom of translation)

		if(FA->npar==FA->MIN_PAR){ realloc_par(FA,&FA->MIN_PAR); }

		FA->map_par[FA->npar].typ = -1; //anchor point in space (3 degrees of freedom of translation)
		FA->map_par[FA->npar].atm = residue[at].gpa[0];
		FA->opt_par[FA->npar] = 0.0; //sets default position to sphere index 0
		FA->del_opt_par[FA->npar] = FA->delta_index;
		FA->min_opt_par[FA->npar] = FA->index_min;
		FA->max_opt_par[FA->npar] = FA->index_max;
		FA->map_opt_par[FA->npar] = 1;
        
		atoms[residue[at].gpa[0]].par = &FA->map_par[FA->npar];
		
		printf("npar=%d map_par[%d].typ=%d map_par[%d].atm=%d opt_par[%d]=%f\n",FA->npar,
		       FA->npar,FA->map_par[FA->npar].typ,
		       FA->npar,atoms[FA->map_par[FA->npar].atm].number,
		       FA->npar,FA->opt_par[FA->npar]);
		//PAUSE;
        
		FA->npar++;
		
		FA->translational = 1;
		
	}else if(val[1] == 0){          // (3 degrees of freedom of rotation)
    
		if(FA->npar==FA->MIN_PAR){ realloc_par(FA,&FA->MIN_PAR); }
    
		for (i=1;i<=3;i++) {
			FA->map_par[FA->npar].bnd=0;

			if (i==1) {
				FA->map_par[FA->npar].typ = 1;  //ang
				FA->map_par[FA->npar].atm = residue[at].gpa[1];
				FA->opt_par[FA->npar] = atoms[FA->map_par[FA->npar].atm].ang;
				FA->del_opt_par[FA->npar] = FA->delta_angle;
				FA->min_opt_par[FA->npar] = -180.0;
				FA->max_opt_par[FA->npar] = 180.0;
				FA->map_opt_par[FA->npar] = 0;

				atoms[residue[at].gpa[1]].par = &FA->map_par[FA->npar];
				
			}else if (i==2) {
				FA->map_par[FA->npar].typ = 2;  //dih
				FA->map_par[FA->npar].bnd = -1;
				FA->map_par[FA->npar].atm=residue[at].gpa[1];
				FA->opt_par[FA->npar]=atoms[FA->map_par[FA->npar].atm].dih;
				FA->del_opt_par[FA->npar] = FA->delta_dihedral;
				FA->min_opt_par[FA->npar] = -180.0;
				FA->max_opt_par[FA->npar] = 180.0;
				FA->map_opt_par[FA->npar] = 0;
				
				atoms[residue[at].gpa[1]].par = &FA->map_par[FA->npar];

			}else if (i==3) {
				FA->map_par[FA->npar].typ = 2;
				FA->map_par[FA->npar].atm=residue[at].gpa[2];
				FA->opt_par[FA->npar]=atoms[FA->map_par[FA->npar].atm].dih;
				FA->del_opt_par[FA->npar] = FA->delta_dihedral;
				FA->min_opt_par[FA->npar] = -180.0;
				FA->max_opt_par[FA->npar] = 180.0;
				FA->map_opt_par[FA->npar] = 0;
				
				atoms[residue[at].gpa[2]].par = &FA->map_par[FA->npar];

			}

			printf("npar=%d map_par[%d].typ=%d map_par[%d].atm=%d opt_par[%d]=%f\n",FA->npar,
			       FA->npar,FA->map_par[FA->npar].typ,
			       FA->npar,atoms[FA->map_par[FA->npar].atm].number,
			       FA->npar,FA->opt_par[FA->npar]);
			//PAUSE;
            
			FA->npar++;
		}
  
	}else if(val[1] > 0){   // val[1] > 0 (dihedrals rotation)
        
		if(FA->npar==FA->MIN_PAR){ realloc_par(FA,&FA->MIN_PAR); }
		
		/* dihedral angle optimization */
		//FA->intramolecular = 1;
		FA->map_par[FA->npar].typ = 2;
		FA->map_par[FA->npar].bnd = val[1];
		FA->map_par[FA->npar].atm = residue[at].bond[val[1]];
		FA->opt_par[FA->npar] = atoms[FA->map_par[FA->npar].atm].dih;
		FA->del_opt_par[FA->npar] = FA->delta_flexible;
		FA->min_opt_par[FA->npar] = -180.0;
		FA->max_opt_par[FA->npar] = 180.0;
		FA->map_opt_par[FA->npar] = 0;

		atoms[residue[at].bond[val[1]]].par = &FA->map_par[FA->npar];
		
		FA->nflexbonds++;

		printf("npar=%d map_par[%d].typ=%d map_par[%d].atm=%d opt_par[%d]=%f\n",FA->npar,
		       FA->npar,FA->map_par[FA->npar].typ,
		       FA->npar,atoms[FA->map_par[FA->npar].atm].number,
		       FA->npar,FA->opt_par[FA->npar]);
		//PAUSE;

		if(FA->map_par_flexbond_first == NULL){
			FA->map_par_flexbond_first = &FA->map_par[FA->npar];
			FA->map_par_flexbond_first_index = FA->npar;
		}
		FA->map_par_flexbond_last = &FA->map_par[FA->npar];
        
		FA->npar++;
	}
	
	return;
}


// Capacity reservation. Added during the SIGSEGV investigation; it was not the
// fix (see MECHANISM CORRECTED below).
//
// add2_optimiz_vec() stores raw pointers into FA->map_par (atoms[].par,
// map_par_sidechain_first, map_par_sidechain_last) and realloc_par() relocates
// the array, so a mid-loop realloc leaves earlier captures out of range; the
// flexible-side-chain loop takes one per residue, so with >=2 flexible residues a
// realloc can fire between two captures. reserve_par() grows to a capacity known
// up front so no realloc occurs while pointers are being taken, and is a no-op
// when capacity already suffices, so it cannot change a correct path.
//
// MECHANISM CORRECTED 2026-09-05 (external code audit). What is measured and what
// is not:
//   MEASURED  add2_optimiz_vec() stores raw pointers into FA->map_par and
//             realloc_par() relocates it, so a growth mid-setup leaves earlier
//             captures outside the live array. validate_map_par_pointers()
//             reports this as [PAR-PTR]; 11 of 18 on 1R55, also 4 of 15 and
//             3 of 16, across 21 logs on disk.
//   MEASURED  exit 139 on 1R55 / 1R58 / 1HQ2 with as few as one flexible side
//             chain at population >= 400; population 100 and the rigid arm clean.
//   REFUTED   that the out-of-range map_par pointer is what faults. GREPPED AT
//             HEAD: atoms[].par is assigned 8 times and DEREFERENCED ZERO times
//             in LIB; map_par_sidechain_first/last likewise, their only
//             non-assignment read being inside validate_map_par_pointers itself.
//             All three are WRITE-ONLY, so a stale one cannot fault.
//   ACTUAL    the crash was ASan-localised to atoms[].optres, dereferenced at
//             vcfunction.cpp:796-798 and GISTEvaluator.cpp:206, against an array
//             relocated by build_rotamers.cpp:308 (see vcfunction.cpp:748).
//             reserve_optres() is the guard for that array.
// So this reservation buys uniform map_par lifetime and a clean detector, not
// crash safety. Full discussion at top.cpp:2661.
//
// reserve_par() grows the arrays to a capacity known up front, so no realloc can
// occur while pointers are being taken. It is a no-op when capacity suffices, so
// it cannot change behaviour on any path that was already correct.
// ── PRODUCTION GUARD 2 of 2: VALIDATE, do not assume ────────────────────────
// Asserts every non-null atoms[].par points INSIDE the live map_par array.
// A violation means a relocation happened while pointers were live. That is a
// real invariant break and worth reporting -- but it is NOT the crash cause:
// atoms[].par is write-only (dereferenced 0x in LIB), so a stale one cannot
// fault. The SIGSEGV was atoms[].optres (vcfunction.cpp:796-798). See the
// MECHANISM CORRECTED block above. Report loudly anyway: an invariant that
// breaks silently is how the optres defect stayed hidden.
// FLEXAIDDS_PAR_PTR_WARN=1 downgrades to a warning for bisecting.
int validate_map_par_pointers(FA_Global* FA, atom* atoms, int atm_cnt,
                              const char* where)
{
	if(FA == NULL || FA->map_par == NULL || atoms == NULL) return 0;
	const optmap* lo = FA->map_par;
	const optmap* hi = FA->map_par + FA->npar;   // one past the last live gene
	int bad = 0, checked = 0;
	for(int i = 1; i <= atm_cnt; i++){
		const optmap* p = atoms[i].par;
		if(p == NULL) continue;
		checked++;
		if(p < lo || p >= hi) bad++;
	}
	const optmap* sc[2] = { FA->map_par_sidechain_first, FA->map_par_sidechain_last };
	for(int k = 0; k < 2; k++){
		if(sc[k] == NULL) continue;
		checked++;
		if(sc[k] < lo || sc[k] >= hi) bad++;
	}
	if(bad > 0){
		fprintf(stderr,
			"[PAR-PTR] %s: %d of %d captured map_par pointer(s) lie OUTSIDE the live "
			"array [%p,%p) (npar=%d MIN_PAR=%d). A realloc relocated map_par while "
			"pointers were live -- this is the dangling-pointer defect, not a data issue.\n",
			(where ? where : "?"), bad, checked, (const void*)lo, (const void*)hi,
			FA->npar, FA->MIN_PAR);
		// WARN by default: this condition predates the guard and is reachable in
		// rigid runs, so a fatal default would halt every production run rather
		// than surface the defect. Opt in to aborting when bisecting.
		const char* fatal = getenv("FLEXAIDDS_PAR_PTR_FATAL");
		if(fatal && fatal[0] == '1') Terminate(2);
	}
	return bad;
}

// Reserve FA->optres capacity so the array cannot relocate while atoms[].optres
// pointers are live. Same defect class as map_par; see update_optres.cpp.
// FLEXAIDDS_PAR_RESERVE=0 also disables this, for baseline measurement.
void reserve_optres(FA_Global* FA, int need){
	const char* rsv = getenv("FLEXAIDDS_PAR_RESERVE");
	if(rsv && rsv[0] == '0') return;
	if(need < 1) need = 1;
	if(FA->MIN_OPTRES >= need) return;
	OptRes* grown = (OptRes*)realloc(FA->optres, (size_t)need * sizeof(OptRes));
	if(!grown){
		fprintf(stderr,"ERROR: reserve_optres could not allocate %d entries\n", need);
		Terminate(2);
	}
	memset(&grown[FA->MIN_OPTRES], 0,
	       (size_t)(need - FA->MIN_OPTRES) * sizeof(OptRes));
	FA->optres = grown;
	FA->MIN_OPTRES = need;
}

void reserve_par(FA_Global* FA, int need){
	// FLEXAIDDS_PAR_RESERVE=0 disables the reservation, restoring STOCK growth
	// behaviour (MIN_PAR=6 + fixed steps) so the validator can measure BASELINE
	// exposure. Without this switch, any violation count is confounded by the
	// reservation's own relocation and cannot establish whether unpatched runs
	// were affected.
	const char* rsv = getenv("FLEXAIDDS_PAR_RESERVE");
	if(rsv && rsv[0] == '0') return;
	int guard = 0;
	while(FA->MIN_PAR <= need){
		realloc_par(FA,&FA->MIN_PAR);
		if(++guard > 1000){
			fprintf(stderr,"ERROR: reserve_par runaway (need=%d MIN_PAR=%d)\n",
			        need, FA->MIN_PAR);
			Terminate(2);
		}
	}
}

void realloc_par(FA_Global* FA, int* MIN_PAR){
	
	*MIN_PAR += N_NEW_PAR;

	FA->map_par = (optmap*)realloc(FA->map_par,(*MIN_PAR)*sizeof(optmap));
	FA->opt_par = (double*)realloc(FA->opt_par,(*MIN_PAR)*sizeof(double));
	FA->del_opt_par = (double*)realloc(FA->del_opt_par,(*MIN_PAR)*sizeof(double));
	FA->min_opt_par = (double*)realloc(FA->min_opt_par,(*MIN_PAR)*sizeof(double));
	FA->max_opt_par = (double*)realloc(FA->max_opt_par,(*MIN_PAR)*sizeof(double));
	FA->map_opt_par = (int*)realloc(FA->map_opt_par,(*MIN_PAR)*sizeof(int));
	
	if(!FA->map_par || !FA->opt_par ||
	   !FA->del_opt_par || !FA->min_opt_par || 
	   !FA->max_opt_par || !FA->map_opt_par){
		
		fprintf(stderr,"ERROR: memory allocation error for opt_par\n");
		Terminate(2);
	}
	
	memset(&FA->map_par[(*MIN_PAR)-N_NEW_PAR],0,N_NEW_PAR*sizeof(optmap));
	memset(&FA->opt_par[(*MIN_PAR)-N_NEW_PAR],0,N_NEW_PAR*sizeof(double));
	memset(&FA->del_opt_par[(*MIN_PAR)-N_NEW_PAR],0,N_NEW_PAR*sizeof(double));
	memset(&FA->min_opt_par[(*MIN_PAR)-N_NEW_PAR],0,N_NEW_PAR*sizeof(double));
	memset(&FA->max_opt_par[(*MIN_PAR)-N_NEW_PAR],0,N_NEW_PAR*sizeof(double));
	memset(&FA->map_opt_par[(*MIN_PAR)-N_NEW_PAR],0,N_NEW_PAR*sizeof(int));
      	
}
