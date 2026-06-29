#include "gaboom.h"

#include <cstring>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

// Per-thread scratch for calc_rmsd_chrom: shallow copies of atoms[] and
// residue[] so parallel coord-cache fills never mutate shared GA state.
struct CalcRmsdChromScratch {
	std::vector<atom>   atoms;
	std::vector<resid>  residue;
	std::vector<double> opt_par;
};
static thread_local CalcRmsdChromScratch tl_calc_rmsd_chrom_scratch;

/******************************************************************************
 * SUBROUTINE calc_rmsd_chrom calculates the rmsd between any two chromosomes  
 * present in the population.
 *****************************************************************************/
float calc_rmsd_chrom(FA_Global* FA, GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, int chrom_a, int chrom_b,
                      float* coor_a_dest, float* coor_b_dest, bool calc_rmsd){

	float rmsd_chrom=0.0f;
	int i = 0,j = 0,k = 0,l = 0,m = 0;
	int cat;
	int rot;

	// Thread-local copies: callers may invoke this in parallel (coord cache).
	const int atm_n = FA->atm_cnt;
	const int res_n = FA->res_cnt + 1;
	CalcRmsdChromScratch& scr = tl_calc_rmsd_chrom_scratch;
	if ((int)scr.atoms.size() < atm_n)    scr.atoms.resize(static_cast<std::size_t>(atm_n));
	if ((int)scr.residue.size() < res_n)  scr.residue.resize(static_cast<std::size_t>(res_n));
	if ((int)scr.opt_par.size() < npar)   scr.opt_par.resize(static_cast<std::size_t>(npar));
	if (atm_n > 0)
		std::memcpy(scr.atoms.data(), atoms, static_cast<std::size_t>(atm_n) * sizeof(atom));
	if (res_n > 0)
		std::memcpy(scr.residue.data(), residue, static_cast<std::size_t>(res_n) * sizeof(resid));
	atom*  work_atoms    = scr.atoms.data();
	resid* work_residue  = scr.residue.data();
	double* opt_par      = scr.opt_par.data();
    
    float coor_a[MAX_ATM_HET*3];
    float coor_b[MAX_ATM_HET*3];
    
    if(coor_a_dest == NULL){
        coor_a_dest = coor_a;
    }
    
    if(coor_b_dest == NULL){
        coor_b_dest = coor_b;
    }
    
	uint grd_idx;
	int normalmode=-1;
	int rot_idx=0;

	for(k=0;k<2;k++)
	{
		j=chrom_a;
		if(k==1){j=chrom_b;}
        
		for(i=0;i<npar;i++){ opt_par[i] = chrom[j].genes[i].to_ic; }
  
		/*
		  printf("%2d (",j);
		  for(l=0;l<GB->num_genes;l++) printf("%12.6f ",opt_par[l]);
		  printf(") ");
		  printf(" value=%11.6f fitnes=%11.6f\n",chrom[j].evalue,chrom[j].fitnes);
		*/
    
    
		for(i=0;i<npar;i++)
		{
			//printf("[%8.3f]",opt_par[i]);
      
			if(FA->map_par[i].typ==-1) 
			{ //by index
				grd_idx = (uint)opt_par[i];
				//printf("opt_par(index): %d\n", grd_idx);
				//PAUSE;
				work_atoms[FA->map_par[i].atm].dis = cleftgrid[grd_idx].dis;
				work_atoms[FA->map_par[i].atm].ang = cleftgrid[grd_idx].ang;
				work_atoms[FA->map_par[i].atm].dih = cleftgrid[grd_idx].dih;
	
			}
			else if(FA->map_par[i].typ==0) 
			{
				work_atoms[FA->map_par[i].atm].dis = (float)opt_par[i];
			}
			else if(FA->map_par[i].typ==1) 
			{
				work_atoms[FA->map_par[i].atm].ang = (float)opt_par[i];
			}
			else if(FA->map_par[i].typ==2)
			{
				work_atoms[FA->map_par[i].atm].dih = (float)opt_par[i];
	
				j=FA->map_par[i].atm;
				cat=work_atoms[j].rec[3];
				if(cat != 0){
					while(cat != FA->map_par[i].atm){
						work_atoms[cat].dih=work_atoms[j].dih + work_atoms[cat].shift; 
						j=cat;
						cat=work_atoms[j].rec[3];
					}
				}
			}else if(FA->map_par[i].typ==3)
			{ //by index
				grd_idx = (uint)opt_par[i];
				//printf("opt_par(index): %d\n", grd_idx);
				//PAUSE;
	
				// serves as flag , but also as grid index
				normalmode=grd_idx;
	
			}else if(FA->map_par[i].typ==4)
			{
				rot_idx = (int)(opt_par[i]+0.5);
	
				work_residue[work_atoms[FA->map_par[i].atm].ofres].rot=rot_idx;
	
				/*
				  printf("residue[%d].rot[%d] - fatm=%d - latm=%d\n",
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].number,
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].rot,
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].fatm[rot_idx],
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].latm[rot_idx]);
				*/
			}
      
		}

		if(normalmode > -1)
			alter_mode(work_atoms,work_residue,FA->normal_grid[normalmode],FA->res_cnt,FA->normal_modes);
  
		/* rebuild cartesian coordinates of optimized residues*/
		for(i=0;i<FA->nors;i++){
			/*printf("nors=%d opt_res[%d]=%d nmov[%d]=%d\n",
			  i,i,FA->opt_res[i],i,FA->nmov[i]);
			  for(j=0;j<FA->nmov[i];j++){
			  printf("mov[%d][%d]=%d\n",i,j,FA->mov[i][j]);
			  }
			  PAUSE;*/
			buildcc(FA,work_atoms,FA->nmov[i],FA->mov[i]);
		}

		// residue that is optimized geometrically (ligand)
		l=work_atoms[FA->map_par[0].atm].ofres;

		m=0;
		rot=work_residue[l].rot;
		for(i=work_residue[l].fatm[rot];i<=work_residue[l].latm[rot];i++){
			//printf("i:%d %f %f %f\n",i,work_atoms[i].coor[0],
			//     work_atoms[i].coor[1],
			//     work_atoms[i].coor[2]);
			for(j=0;j<3;j++)
			{
				if(k==0) coor_a_dest[m*3+j]=work_atoms[i].coor[j];
				if(k==1) coor_b_dest[m*3+j]=work_atoms[i].coor[j];
			}
			m++;
		}    
	}
  
    if(calc_rmsd){
        for(i=0;i<m;i++)
            rmsd_chrom += sqrdist(&coor_a_dest[i*3],&coor_b_dest[i*3]);
        
        rmsd_chrom = sqrt(rmsd_chrom/((float)m));
    }
	//printf("RMSD=%f\n",rmsd_chrom);
	//PAUSE;

	return rmsd_chrom;
}