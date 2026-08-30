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
    
	return;
    
}
