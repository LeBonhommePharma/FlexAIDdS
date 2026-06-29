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

static bool calc_rmsd_atom_ok(int idx, int atm_n) {
	return idx >= 0 && idx < atm_n;
}

static bool calc_rmsd_res_ok(int idx, int res_n) {
	return idx >= 0 && idx < res_n;
}

/******************************************************************************
 * SUBROUTINE calc_rmsd_chrom calculates the rmsd between any two chromosomes  
 * present in the population.
 *****************************************************************************/
float calc_rmsd_chrom(FA_Global* FA, GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, int chrom_a, int chrom_b,
                      float* coor_a_dest, float* coor_b_dest, bool calc_rmsd, int* out_n_atoms){

	float rmsd_chrom=0.0f;
	int i = 0,k = 0,l = 0,m = 0;
	int chrom_idx = 0;
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
		normalmode=-1;
		chrom_idx=chrom_a;
		if(k==1){chrom_idx=chrom_b;}
        
		// Mirror eval_chromosome(): clamp ICs to gene_lim before applying them.
		// Snapshot chromosomes can carry raw crossover values; unclamped grid
		// indices dereference cleftgrid[] out of bounds → intermittent SIGSEGV
		// in the post-GA clustering coord-cache fill.
		for(i=0;i<npar;i++){
			double ic = chrom[chrom_idx].genes[i].to_ic;
			if(gene_lim != nullptr){
				if(ic > gene_lim[i].max) ic = gene_lim[i].max;
				else if(ic < gene_lim[i].min) ic = gene_lim[i].min;
			}
			opt_par[i] = ic;
		}
  
		/*
		  printf("%2d (",chrom_idx);
		  for(l=0;l<GB->num_genes;l++) printf("%12.6f ",opt_par[l]);
		  printf(") ");
		  printf(" value=%11.6f fitnes=%11.6f\n",chrom[chrom_idx].evalue,chrom[chrom_idx].fitnes);
		*/
    
    
		for(i=0;i<npar;i++)
		{
			const int atm_i = FA->map_par[i].atm;
			if(!calc_rmsd_atom_ok(atm_i, atm_n))
				continue;

			//printf("[%8.3f]",opt_par[i]);
      
			if(FA->map_par[i].typ==-1) 
			{ //by index
				grd_idx = (uint)opt_par[i];
				if(FA->num_grd > 0 && grd_idx >= (uint)FA->num_grd)
					grd_idx = (uint)FA->num_grd - 1;
				//printf("opt_par(index): %d\n", grd_idx);
				//PAUSE;
				work_atoms[atm_i].dis = cleftgrid[grd_idx].dis;
				work_atoms[atm_i].ang = cleftgrid[grd_idx].ang;
				work_atoms[atm_i].dih = cleftgrid[grd_idx].dih;
	
			}
			else if(FA->map_par[i].typ==0) 
			{
				work_atoms[atm_i].dis = (float)opt_par[i];
			}
			else if(FA->map_par[i].typ==1) 
			{
				work_atoms[atm_i].ang = (float)opt_par[i];
			}
			else if(FA->map_par[i].typ==2)
			{
				const int root_atm = atm_i;
				work_atoms[root_atm].dih = (float)opt_par[i];
	
				int atm_j = root_atm;
				cat = work_atoms[atm_j].rec[3];
				if(cat != 0 && cat != root_atm){
					int steps = 0;
					while(cat != root_atm && steps < atm_n){
						if(!calc_rmsd_atom_ok(cat, atm_n))
							break;
						work_atoms[cat].dih = work_atoms[atm_j].dih + work_atoms[cat].shift;
						atm_j = cat;
						cat = work_atoms[atm_j].rec[3];
						++steps;
					}
				}
			}else if(FA->map_par[i].typ==3)
			{ //by index
				grd_idx = (uint)opt_par[i];
				if(FA->normal_modes > 0 && grd_idx >= (uint)FA->normal_modes)
					grd_idx = (uint)FA->normal_modes - 1;
				//printf("opt_par(index): %d\n", grd_idx);
				//PAUSE;
	
				// serves as flag , but also as grid index
				normalmode=(int)grd_idx;
	
			}else if(FA->map_par[i].typ==4)
			{
				rot_idx = (int)(opt_par[i]+0.5);
				const int res_i = work_atoms[atm_i].ofres;
				if(calc_rmsd_res_ok(res_i, res_n)){
					const int trot  = work_residue[res_i].trot;
					if(trot > 0){
						if(rot_idx < 0) rot_idx = 0;
						else if(rot_idx >= trot) rot_idx = trot - 1;
					}
					work_residue[res_i].rot = rot_idx;
				}
	
				/*
				  printf("residue[%d].rot[%d] - fatm=%d - latm=%d\n",
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].number,
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].rot,
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].fatm[rot_idx],
				  work_residue[work_atoms[FA->map_par[i].atm].ofres].latm[rot_idx]);
				*/
			}
      
		}

		if(normalmode > -1 && FA->normal_modes > 0 && normalmode < FA->normal_modes)
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
		m=0;
		const int lig_atm = FA->map_par[0].atm;
		if(!calc_rmsd_atom_ok(lig_atm, atm_n))
			continue;
		l = work_atoms[lig_atm].ofres;
		if(!calc_rmsd_res_ok(l, res_n))
			continue;

		rot=work_residue[l].rot;
		if(work_residue[l].trot > 0){
			if(rot < 0) rot = 0;
			else if(rot >= work_residue[l].trot) rot = work_residue[l].trot - 1;
		}
		const int max_coord_atoms = MAX_ATM_HET;
		int fatm = work_residue[l].fatm[rot];
		int latm = work_residue[l].latm[rot];
		if(!calc_rmsd_atom_ok(fatm, atm_n) || !calc_rmsd_atom_ok(latm, atm_n) || fatm > latm)
			continue;
		for(i=fatm;i<=latm;i++){
			if(!calc_rmsd_atom_ok(i, atm_n))
				break;
			//printf("i:%d %f %f %f\n",i,work_atoms[i].coor[0],
			//     work_atoms[i].coor[1],
			//     work_atoms[i].coor[2]);
			if(m >= max_coord_atoms) break;
			for(int d=0;d<3;d++)
			{
				if(k==0) coor_a_dest[m*3+d]=work_atoms[i].coor[d];
				if(k==1) coor_b_dest[m*3+d]=work_atoms[i].coor[d];
			}
			m++;
		}
		if(k == 0 && out_n_atoms != nullptr)
			*out_n_atoms = m;
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