#include "flexaid.h"
#include "fileio.h"

// the function updates the matrix of bonded atoms
// the matrix starts at [0,0] with residue[ires].fatm
// works/available for rotamers also

// Releases residue->bonded allocated by update_bonded().
// tot must be the same value passed to the matching call.

void free_bonded(resid* residue, int tot)
{
  if(residue == nullptr || residue->bonded == nullptr){ return; }

  for(int i=0; i<tot; i++)
    {
      free(residue->bonded[i]);
    }
  free(residue->bonded);
  residue->bonded = nullptr;
}

void update_bonded(resid* residue, int tot, int nlist, int* list, int* nbr)
{

  int i,j;             // dumb counters
  int fatm,latm;       // first-last atom of residue i

  /**************************************************/
  /*******     allocate memory for matrix ***********/
  /**************************************************/
  
  // The guard covers both allocation levels, so a second call on the same
  // residue reuses rows sized by the FIRST call's tot. That is safe only
  // because tot is the same expression at every call site --
  // latm[0]-fatm[0]+1, rot-independent -- at read_lig.cpp:464,
  // build_rotamers.cpp:337 (the genuine re-entry, on protein residues),
  // Mol2Reader.cpp:422 and SdfReader.cpp:752. If a caller ever passes a
  // larger tot, the loops below walk past the previous allocation.
  if(residue->bonded == NULL)
    {

      residue->bonded = (int**)malloc(tot*sizeof(int*));
      
      if(residue->bonded == NULL)
	{
	  fprintf(stderr,"ERROR: Could not allocate memory for residue->bonded\n");
	  Terminate(2);
	}
      else
	{
	  for(i=0;i<tot;i++){

	    residue->bonded[i] = (int*)malloc(tot*sizeof(int));

	    if(residue->bonded[i] == NULL)
	      {
		fprintf(stderr,"ERROR: Could not allocate memory for residue->bonded[%d]\n",i);
		Terminate(2);
	      }	    
	    
	  }
	    
      }
     
      // fill matrix with -1 (nothing found up to nloops)
      for(i=0;i<tot;i++)
	for(j=0;j<tot;j++)
	  residue->bonded[i][j] = -1;
	
      /*
	printf("new bonded residue!\n");
	getchar();
      */
    }


  /**************************************************/
  /*******     update content of matrix   ***********/
  /**************************************************/

  /*
    printf("updating...\nmatrix so far...\n");
    for(i=0;i<tot;++i){for(j=0;j<tot;++j){printf("[%d-%d]=%d\n",i,j,residue->bonded[i][j]);}}
  */

  //  printf("will write matrix to %p\n",residue->bonded);

  // assign first atom to start from 0
  fatm = residue->fatm[0];
  latm = residue->latm[0];
  
  //printf("fatm:%d\tlatm:%d\n",fatm,latm);

  // loop through neighbours list
  for(i=0;i<nlist;i++){
    
    /*
    printf("writing %d against %d with value %d\n",list[0]-fatm,list[i]-fatm,nbr[i]);
    printf("internal atom numbers list[0]=%d\tlist[%d]=%d\tfatm=%d\n",
    	   list[0],i,list[i],fatm);
    */
    
    if(list[i]-fatm >= 0)
      {
	if(list[i]-fatm < (latm-fatm+1))
	  {
	    residue->bonded[list[0]-fatm][list[i]-fatm] = nbr[i];
	  }
      }
  }


  return;

}
