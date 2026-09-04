#include <cstdlib>
#include <cstdio>
#include "flexaid.h"
#include "fileio.h"
#include "soft_wall.h"  // posebusters_vdw_radius() for the pb_vdw_radius cache

/*****************************************************/
/******* this function updates the atoms str  ********/
/******* to assign the pointer pointing at    ********/
/******* the residue to optimize during dock  ********/
/*****************************************************/

void update_optres(atom* atoms, resid* residue, int atm_cnt, OptRes* optres_ptr,int num_optres)
{
    
	int i,j;

	// ── pb-vdw-precompute ──────────────────────────────────────────────
	// Cache each atom's PoseBusters vdW radius ONCE so the pb_clash hot loop
	// in get_contlist4() (Vcontacts.cpp) reads a double instead of running a
	// ~26-branch element string-compare per atom, twice per pair, for every
	// one of the ~2M CF evals per restart.
	//
	// Correctness: the cached value is byte-for-byte the same call the hot
	// loop makes today — posebusters_vdw_radius(atoms[j].element, atoms[j].radius)
	// — so the clash penalty is bit-identical. This runs on the master atoms[]
	// array on every load path (top.cpp / read_input.cpp / direct_input.cpp all
	// call update_optres() here, AFTER assign_radii_types() and build_rotamers()
	// have finalized both .element and .radius), covering receptor, cofactor,
	// ion, rotamer and ligand atoms. The GA per-thread atom copies are struct-
	// copied from this master array every evaluation (gaboom.cpp tl_atoms), so
	// the cache propagates into every thread with no extra work and no races.
	for(j=1;j<=atm_cnt;j++){
		atoms[j].pb_vdw_radius =
			posebusters_vdw_radius(atoms[j].element, atoms[j].radius);
	}

	// CLEAR FIRST. Without this pass an atom keeps a pointer into a PREVIOUS
	// FA->optres allocation forever, because the loop below only assigns on a
	// match and never resets a non-match. FA->optres is realloc'd per flexible
	// residue (build_rotamers.cpp:308), so a retained pointer is a freed one.
	// This function is the single owner of the atoms[].optres mapping and runs
	// AFTER build_rotamers in both entry paths (top.cpp:2695, read_input.cpp:707),
	// so rebuilding the mapping from scratch here is both correct and sufficient.
	for(j=1;j<=atm_cnt;j++){ atoms[j].optres = NULL; }

	for(i=0;i<num_optres;i++){
        
		for(j=1;j<=atm_cnt;j++){
			
			if(atoms[j].ofres == optres_ptr[i].rnum)
			{
				if(residue[atoms[j].ofres].type != 0 || !atoms[j].isbb)
				{
					atoms[j].optres = &optres_ptr[i];
					//printf("atom %d has optres\n", atoms[j].number);
				}
			}
            
		}
        
	}

	// ── DIAGNOSTIC, opt-in via FLEXAIDDS_OPTRES_DIAG. Default output unchanged.
	//
	// WHY THIS EXISTS. FLEXAIDDS_SCORED_ONLY writes only atoms carrying an optres
	// back-pointer (write_pdb.cpp:44,150). MEASURED: a flexible cell writes EXACTLY
	// the same 36 ligand atoms as the rigid one, so the flexed side chains are not
	// getting a pointer. Three candidate causes were ruled out by reading source --
	// isbb propagates through the struct copy at build_rotamers.cpp:157; the
	// condition above DOES admit side-chain atoms (protein type==0 AND !isbb); and
	// reserve_optres() already runs BEFORE build_rotamers on all three load paths,
	// so the pointer-invalidation route is closed. The remaining candidate is that
	// the flexible residues never became OptRes entries at all --
	// build_rotamers.cpp:352 is the only incrementer of num_optres -- which this
	// print settles outright: num_optres==1 means ligand only.
	// PERMANENT GUARD, always on. An OptRes entry with rnum==0 at index>0 is
	// unmapped: no atom has ofres==0, so the match above cannot fire and that
	// residue is silently absent from every per-residue consumer
	// (ic2cf.cpp:299, cluster.cpp:732, BindingMode.cpp:921/1091, FOPTICS.cpp:407,
	// DensityPeak_Cluster.cpp:614, top.cpp:2886) AND from FLEXAIDDS_SCORED_ONLY
	// output. That is the failure this project keeps paying for: a capability
	// fully present in the engine whose precondition is never supplied, with every
	// log, config and exit code looking clean. It must never be silent again.
	{
		int unmapped = 0;
		for(i=1;i<num_optres;i++){ if(optres_ptr[i].rnum == 0) ++unmapped; }
		if(unmapped > 0){
			fprintf(stderr,"[OPTRES] WARNING: %d of %d optimizable-residue entries have "
			        "rnum=0 and are UNMAPPED; their side chains will be absent from the "
			        "per-residue CF decomposition and from scored-only output. See "
			        "build_rotamers.cpp (slot index) -- this is a defect, not a warning "
			        "to ignore.\n", unmapped, num_optres);
			fflush(stderr);
		}
	}

	{
		const char* diag = getenv("FLEXAIDDS_OPTRES_DIAG");
		if(diag != NULL && diag[0] == '1'){
			int with = 0, prot = 0, lig = 0;
			for(j=1;j<=atm_cnt;j++){
				if(atoms[j].optres == NULL) continue;
				++with;
				if(residue[atoms[j].ofres].type == 0) ++prot; else ++lig;
			}
			printf("[OPTRES-DIAG] num_optres=%d atoms_with_optres=%d "
			       "(protein=%d ligand=%d) atm_cnt=%d\n",
			       num_optres, with, prot, lig, atm_cnt);
			for(i=0;i<num_optres;i++){
				printf("[OPTRES-DIAG]   optres[%d] rnum=%d restype=%d\n",
				       i, optres_ptr[i].rnum, residue[optres_ptr[i].rnum].type);
			}
			fflush(stdout);
		}
	}
    
	return;
    
}
