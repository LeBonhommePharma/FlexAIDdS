#include "flexaid.h"

// Sole owner of the per-residue teardown order.
//
// Before this file the sequence existed twice, byte-identical: the production
// copy in top.cpp's cleanup loop and the test copy in tests/cleanup_fa. Nothing
// tied them together, so a change to one silently diverged from the other --
// and the production copy is the one no test target compiles, which is exactly
// the pair that must not drift. Both now call this.
//
// Order is load-bearing. bonded / shortpath / shortflex are three sibling
// matrices dimensioned by the same natm, and natm is recovered from the
// residue's own fatm/latm -- so those two must still be readable when the trio
// is released. Freeing fatm/latm first loses the dimension and the free walks
// the wrong length.

void free_resid(resid* residue)
{
	if(residue == nullptr){ return; }

	// The natm derivation is guarded on the trio rather than run
	// unconditionally: read_pdb.cpp:42,51,52,55,56 leaves slot 0 with all three
	// NULL and its fatm[0] never written by anything, so computing natm there
	// would read uninitialised memory for a result the free_* helpers ignore.
	// latm[0] IS written -- read_coor.cpp:299 stores through residue[res_cnt-1]
	// and res_cnt is 1 on the first residue, which is why the read_pdb.cpp:44-45
	// allocation must stay. Guarding here is what lets callers loop from 0 with
	// no special case for that slot.
	if(residue->fatm != nullptr && residue->latm != nullptr &&
	   (residue->bonded    != nullptr ||
	    residue->shortpath != nullptr ||
	    residue->shortflex != nullptr)){

		const int natm = residue->latm[0]-residue->fatm[0]+1;

		free_bonded(residue, natm);
		free_shortpath(residue, natm);
		free_shortflex(residue, natm);
	}

	if(residue->gpa  != nullptr){ free(residue->gpa);  residue->gpa  = nullptr; }
	if(residue->fatm != nullptr){ free(residue->fatm); residue->fatm = nullptr; }
	if(residue->latm != nullptr){ free(residue->latm); residue->latm = nullptr; }
	if(residue->bond != nullptr){ free(residue->bond); residue->bond = nullptr; }
}
