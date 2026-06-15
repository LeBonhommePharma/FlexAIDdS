#include "gaboom.h"
#include "fileio.h"
#include "simd_distance.h"
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

void cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	bool Hungarian = false;
	int i,j;
	cfstr cf;                                /* complementarity function value */
	resid *res_ptr = NULL;
	cfstr* cf_ptr = NULL;

	FILE* outfile_ptr = NULL;

	// will 
	double partition_function = 0.0;

	float rmsd = 0.0f;
	int num_of_results = FA->max_results;
	int num_of_clusters = 0;
	int n_unclus = 0;

	char sufix[24];
	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];

	// Clustering Variable Definitions
	int* Clus_GAPOP = NULL;
	float* Clus_RMSDT = NULL;
	double* Clus_ACF = NULL;
	double* Clus_TCF = NULL;
	int* Clus_TOP = NULL;
	int* Clus_FRE = NULL;

    ////////////////////////////////
    //// memory allocation for /////
    //// for clustering chrom  /////
    ////////////////////////////////
	Clus_GAPOP = (int*)malloc(num_chrom*sizeof(int));		// Population (Cluster)	
	Clus_RMSDT = (float*)malloc(num_chrom*sizeof(float));	// RMSD
	Clus_ACF = (double*)malloc(num_chrom*sizeof(double));	// Apparent CF
	Clus_TCF = (double*)malloc(num_chrom*sizeof(double));	// Total CF
	Clus_TOP = (int*)malloc(num_chrom*sizeof(int));			// Best Chromosome in Cluster
	Clus_FRE = (int*)malloc(num_chrom*sizeof(int));			// Frequency
      
	if(!Clus_GAPOP || !Clus_RMSDT || !Clus_ACF ||
	   !Clus_TCF   || !Clus_TOP   || !Clus_FRE)   
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}
      
    ////////////////////////////////
    //////       END         ///////
    ////////////////////////////////
  
    /******************************************************************/
  
    //-------------------------------------------------------
    // fixed center clustering of chrmosome population around highest ranking
    // solutions with rmsd_threshold angstrons threshold
    // Clus_GAPOP[i]=j assigns for each chromosome i to which cluster it belongs
    // as described by j, the chromosome index of the "cluster head", i.e.,that with the
    // higest CF value.
	n_unclus=num_chrom;
	num_of_clusters=0;
	
	// CLustering Variable Initialization and partition_function calculation
	for(j=0;j<num_chrom;++j)
	{
		Clus_GAPOP[j]=-1;
		Clus_ACF[j]=0.0;
		Clus_TCF[j]=0.0;
		Clus_TOP[j]=0;
		Clus_FRE[j]=0;
		if(FA->temperature){
			partition_function += exp((-1.0) * FA->beta * chrom[j].app_evalue);
		}
	}
    //printf("n_unclus=%d\n",n_unclus);
    //PAUSE;
	
	// Guard: single-pose fast path — PF is undefined (or underflows to 0) for
	// n≤1.  When the GA collapses to 1 chromosome, skip Boltzmann-weighted
	// clustering entirely: there is exactly one cluster of size 1, P=1, and
	// the entropy term T·p·log(p) = 0.  Pre-populate arrays and set n_unclus=0
	// so the while-loop below is skipped; the emit loop runs once for chrom[0].
	if (num_chrom <= 1) {
		n_unclus        = 0;
		num_of_clusters = num_of_results = (num_chrom == 1) ? 1 : 0;
		if (num_chrom == 1) {
			Clus_GAPOP[0]=0; Clus_RMSDT[0]=0.0f;
			Clus_ACF[0]=Clus_TCF[0]=chrom[0].app_evalue;
			Clus_TOP[0]=0;   Clus_FRE[0]=1;
		}
	} else if (FA->temperature && partition_function == 0.0) {
		fprintf(stderr,"ERROR: The Partition Function is NULL in the clustering step.\n");
		Terminate(2);
	}
	
	// ── Pre-compute Cartesian coordinates for all chromosomes ──────────
	// Eliminates redundant buildcc() calls in the O(N²) clustering loop.
	// DensityPeak_Cluster already uses this pattern (Chrom[i].Coord).
	const int nAtoms_clus = residue[atoms[FA->map_par[0].atm].ofres].latm[0]
	                      - residue[atoms[FA->map_par[0].atm].ofres].fatm[0] + 1;
	const int coord_stride = nAtoms_clus * 3;
	std::vector<float> coord_cache((size_t)num_chrom * coord_stride);

	for(int c = 0; c < num_chrom; ++c)
	{
		if(c + 1 < num_chrom) {
			calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,
			                GB->num_genes, c, c+1,
			                &coord_cache[c * coord_stride],
			                &coord_cache[(c+1) * coord_stride], false);
		} else {
			calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,
			                GB->num_genes, c, c,
			                &coord_cache[c * coord_stride], NULL, false);
		}
	}

	// Clustering part — uses cached coordinates + SIMD distance
	while(n_unclus > 0)
	{
		for(j=0;j<num_chrom;++j){if(Clus_GAPOP[j]==-1){break;}}
		//printf("at chromosome j=%d with app_evalue=%.3f\n", j, chrom[j].app_evalue);
        if(j < 0 || j >= num_chrom) {break;} // this break intends to break the while() loop if j is invalid (which should not occur but still does)
        Clus_GAPOP[j]=j;
		Clus_RMSDT[j]=0.0;
		n_unclus--;
		if(FA->temperature > 0)
		{
			double Pj = exp((-1.0) * FA->beta * chrom[j].app_evalue) / partition_function;
			Clus_ACF[num_of_clusters] = (double)( ( Pj * chrom[j].app_evalue) + (FA->temperature * Pj * log(Pj)) );
			Clus_TCF[num_of_clusters] = (double)( ( Pj * chrom[j].app_evalue) + (FA->temperature * Pj * log(Pj)) );
		}
		else
		{
			Clus_TCF[num_of_clusters] = chrom[j].app_evalue;
			Clus_ACF[num_of_clusters] = chrom[j].app_evalue;
		}
		Clus_TOP[num_of_clusters]=j;
		Clus_FRE[num_of_clusters]++;

		const float* coor_j = &coord_cache[j * coord_stride];

		// OpenMP parallelised inner RMSD loop over unclustered chromosomes.
		// Each thread computes RMSD from cached coordinates using SIMD.
		// Cluster assignments use a critical section (rare — only on match).
		#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for(i=j+1;i<num_chrom;++i)
		{
			if(Clus_GAPOP[i]==-1)
			{
				const float* coor_i = &coord_cache[i * coord_stride];
				float d = flexaids::sum_sq_distances_f(coor_i, coor_j, coord_stride);
				float loc_rmsd = sqrtf(d / (float)nAtoms_clus);

				if(loc_rmsd <= FA->cluster_rmsd)
				{
					#ifdef _OPENMP
					#pragma omp critical
					#endif
					{
						Clus_GAPOP[i]=j;
						Clus_RMSDT[i]=loc_rmsd;
						n_unclus--;
						if(FA->temperature){
							double Pi = exp((-1.0) * FA->beta * chrom[i].app_evalue) / partition_function;
							Clus_ACF[num_of_clusters] += (double)( (Pi * chrom[i].app_evalue) + (FA->temperature * Pi * log(Pi)) );
						}else{
							Clus_ACF[num_of_clusters] += chrom[i].app_evalue;
						}
						Clus_FRE[num_of_clusters]++;
					}
				}
			}
		}
		num_of_clusters++;

		// quit storing clusters up to N max results
		if(num_of_clusters == num_of_results){break;}
	}

	// Cap num_of_results to actual cluster count before sorting — arrays are
	// malloc'd to num_chrom elements; sorting past num_of_clusters overflows
	// when num_chrom < FA->max_results (e.g. entropy-collapsed run with 2 poses).
	if(num_of_clusters < num_of_results){num_of_results=num_of_clusters;}

	if(FA->temperature)
	{
		// Reordering the clusters properly by lowest ACF values first (after considering cluster's entropy !)
		QuickSort_Clusters(Clus_TOP, Clus_FRE, Clus_TCF, Clus_ACF, Clus_GAPOP, 0, num_of_results-1);
	}

	// ── Emit rank-0 = lowest-CF representative ───────────────────────────────
	// The ACF ordering above is the entropy-augmented cluster free energy, which
	// can rank a densely-populated cluster at WORSE representative CF ahead of a
	// sparse cluster at BETTER CF (e.g. 1KZK: emitted _0.pdb = -390.45 while a
	// -420.62 representative existed in another cluster). DatasetRunner measures
	// both the RMSD and best_score on the lowest "REMARK CF=" pose, and that
	// REMARK is chrom[Clus_TOP[j]].evalue (see emission below). Re-sort the
	// clusters by representative evalue ascending (most negative first) so the
	// emitted _0.pdb is always the lowest-CF cluster representative the GA found.
	// Selection sort over num_of_results (≤ max_results, small) reusing
	// swap_clusters to keep all five parallel arrays consistent.
	for(int a=0; a<num_of_results-1; ++a)
	{
		int best_idx = a;
		for(int b=a+1; b<num_of_results; ++b)
		{
			if(chrom[Clus_TOP[b]].evalue < chrom[Clus_TOP[best_idx]].evalue)
				best_idx = b;
		}
		if(best_idx != a)
			swap_clusters(&Clus_TOP[a], &Clus_FRE[a], &Clus_TCF[a], &Clus_ACF[a], &Clus_GAPOP[a],
			              &Clus_TOP[best_idx], &Clus_FRE[best_idx], &Clus_TCF[best_idx], &Clus_ACF[best_idx], &Clus_GAPOP[best_idx]);
	}

	// print cluster information
	snprintf(tmp_end_strfile, MAX_PATH__, "%s.cad", end_strfile);
	if (FA->htpmode == false)
	{
		if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr))
		{
			Terminate(6);
		}
		else
		{
			for(i=0;i<num_of_clusters;++i)
			{
				fprintf(outfile_ptr,"Cluster %d: TOP=%d TCF=%f ACF=%f freq=%d\n",i,
					Clus_TOP[i],Clus_TCF[i],
					Clus_ACF[i], Clus_FRE[i]);
			}
			if(num_of_clusters > 1)
			{
				fprintf(outfile_ptr,"RMSD between clusters\n");
				for(i=0;i<num_of_clusters;++i)
				{
					for(j=i+1;j<num_of_clusters;++j)
					{
						rmsd=calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,Clus_TOP[i],Clus_TOP[j], NULL, NULL, true);
						fprintf(outfile_ptr,"rmsd(%d,%d)=%f\n",i,j,rmsd);
					}
				}
			}
		}
		CloseFile_B(&outfile_ptr,"w");
	}
        //num_of_results=1;
      
	printf("num_of_clusters=%d num_of_results=%d\n",num_of_clusters,num_of_results);
	
        // output results, 10% of the number of chromosomes or 
        // the number of clusters, the smallest.
      
	for(j=0;j<num_of_results;++j)
	{
		// get parameters of fittest individual in population
		// after optimization -> best docking candidate
    
		// cf=chrom[Clus_TOP[j]].app_evalue;

		for(int k=0; k<GB->num_genes; ++k)
		{
			FA->opt_par[k] = chrom[Clus_TOP[j]].genes[k].to_ic;
		}

		// Ring pucker (LigandRingFlex Phase 2): load the winning chromosome's
		// furanose phases so the emitted pose's Cartesian coords match the
		// puckered conformation the GA scored. No-op unless ring flex is active.
		if (FA->ring_flex_active) {
			for (int s = 0; s < FA->ring_n_sugars && s < MAX_RING_FLEX; ++s)
				FA->ring_cur_phases[s] = chrom[Clus_TOP[j]].ring_phases[s];
		}

		// Rebuild atom coordinates for PDB output and populate the per-optres
		// CF breakdown as a side effect. The returned re-score is NOT used for
		// the REMARK CF line: at emission Vcontacts can early-return a raw
		// uncapped clash penalty (up to ~1e15) that bypasses the per-contact
		// wall cap and blows up the reported CF. Report the stored chromosome
		// evalue instead — that is exactly what the GA optimized.
		cf=ic2cf(FA,VC,atoms,residue,cleftgrid,GB->num_genes,FA->opt_par);

		size_t remark_len = 0;
		remark[0] = '\0';
		safe_remark_cat(remark, "REMARK optimized structure\n", &remark_len);

		snprintf(tmpremark, MAX_REMARK, "REMARK CF=%8.5f\n",chrom[Clus_TOP[j]].evalue);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.app=%8.5f\n",chrom[Clus_TOP[j]].app_evalue);
		safe_remark_cat(remark, tmpremark, &remark_len);

		for(i=0;i<FA->num_optres;++i)
		{
	  
			res_ptr = &residue[FA->optres[i].rnum];
			cf_ptr = &FA->optres[i].cf;
	  
			snprintf(tmpremark, MAX_REMARK, "REMARK optimizable residue %s %c %d\n",
				res_ptr->name,res_ptr->chn,res_ptr->number);
			safe_remark_cat(remark, tmpremark, &remark_len);

			snprintf(tmpremark, MAX_REMARK, "REMARK CF.com=%8.5f\n",cf_ptr->com);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.sas=%8.5f\n",cf_ptr->sas);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.wal=%8.5f\n",cf_ptr->wal);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.con=%8.5f\n",cf_ptr->con);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.gist=%8.5f\n",cf_ptr->gist);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.hbond=%8.5f\n",cf_ptr->hbond);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.metal=%8.5f\n",cf_ptr->metal_coord);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.elec=%8.5f\n",cf_ptr->elec);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK Residue has an overall SAS of %.3f\n",cf_ptr->totsas);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}

		snprintf(tmpremark, MAX_REMARK, "REMARK Cluster %d: Rank (top):%d Average CF:%8.5f Frequency:%d\n",
			j,Clus_TOP[j],Clus_ACF[j],Clus_FRE[j]);
		safe_remark_cat(remark, tmpremark, &remark_len);
		for(i=0;i<FA->npar;++i)
		{
			snprintf(tmpremark, MAX_REMARK, "REMARK [%8.3f]\n",FA->opt_par[i]);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}
		//snprintf(tmpremark, MAX_REMARK, "REMARK seed=%ld\n",FA->seed_ini);
		safe_remark_cat(remark, tmpremark, &remark_len);
		if(FA->refstructure == 1){
			snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n",
				calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian));
			safe_remark_cat(remark, tmpremark, &remark_len);
			Hungarian = true;
			snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n",
				calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian));
			safe_remark_cat(remark, tmpremark, &remark_len);
		}
		snprintf(tmpremark, MAX_REMARK, "REMARK inputs: %s & %s\n",dockinp,gainp);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(sufix, sizeof(sufix), "_%d.pdb",j);
		snprintf(tmp_end_strfile, MAX_PATH__, "%s%s", end_strfile, sufix);
		//printf("filename=<%s>\n",tmp_end_strfile);
		//PAUSE;
		write_pdb(FA,atoms,residue,tmp_end_strfile,remark);

		// ── Write member-CF sidecar (.mcf) for Boltzmann Z+H cluster selection ──
		// DatasetRunner reads these to compute Z_cluster = sum exp(-CF_i/kT)
		// and Shannon entropy H = -sum p_i*log(p_i) for composite cluster scoring.
		// Format: one app_evalue per line; head chromosome first, then members.
		// Uses Clus_GAPOP[k] == Clus_TOP[j] to reconstruct all cluster members.
		{
			std::string mcf_path(tmp_end_strfile);
			if (mcf_path.size() > 4 &&
			    mcf_path.substr(mcf_path.size() - 4) == ".pdb")
				mcf_path = mcf_path.substr(0, mcf_path.size() - 4) + ".mcf";
			FILE* mf = fopen(mcf_path.c_str(), "w");
			if (mf) {
				// Head chromosome (guaranteed first so member_cfs[0] == head CF)
				fprintf(mf, "%.6f\n", chrom[Clus_TOP[j]].app_evalue);
				// Members: scan all chromosomes assigned to this cluster head
				for (int k = 0; k < num_chrom; ++k) {
					if (k != Clus_TOP[j] && Clus_GAPOP[k] == Clus_TOP[j])
						fprintf(mf, "%.6f\n", chrom[k].app_evalue);
				}
				fclose(mf);
			}
		}
	}
      
        // print the RMSD between each chrom. and the reference structure if there is one.
	if(FA->refstructure == 1){write_rrd(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,Clus_GAPOP,Clus_RMSDT,end_strfile); }


	// Clusters memory de-allocation
	if(Clus_GAPOP != NULL) free(Clus_GAPOP);
	if(Clus_RMSDT != NULL) free(Clus_RMSDT);
	if(Clus_ACF   != NULL) free(Clus_ACF);
	if(Clus_TCF   != NULL) free(Clus_TCF);
	if(Clus_TOP   != NULL) free(Clus_TOP);
	if(Clus_FRE   != NULL) free(Clus_FRE);
	
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*                  QuickSort functions for Clusters                   */
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void QuickSort_Clusters(int* TOP, int* FRE, double* TCF, double* ACF, int* GAPOP, int beg, int end)
{
	int l,r,p;
	double piv;

	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		piv = ACF[p];
		
		while(1)
		{
			while( (l<=r) && QS_ASC(ACF[l],piv) <= 0 ) ++l;
			while( (l<=r) && QS_ASC(ACF[r],piv)  > 0 ) --r;
			
			if (l > r) break;
			
			swap_clusters(&TOP[l], &FRE[l], &TCF[l], &ACF[l], &GAPOP[l],&TOP[r], &FRE[r], &TCF[r], &ACF[r], &GAPOP[r]);
			
			if (p == r) p=l;
			++l;--r;
		}
		swap_clusters(&TOP[p], &FRE[p], &TCF[p], &ACF[p], &GAPOP[p],&TOP[r], &FRE[r], &TCF[r], &ACF[r], &GAPOP[r]);
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_Clusters(TOP, FRE, TCF, ACF, GAPOP, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_Clusters(TOP, FRE, TCF, ACF, GAPOP, l, end);
			end = r;
		}
	}
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*                   Swap Function for Clusters                        */
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_clusters(int* TOPx, int* FREx, double* TCFx, double* ACFx, int* GAPOPx, int* TOPy, int* FREy, double* TCFy, double* ACFy, int* GAPOPy)
{
	int TOPt, FREt, GAPOPt;
	double TCFt, ACFt;
	TOPt = *TOPx; *TOPx = *TOPy; *TOPy = TOPt;
	FREt = *FREx; *FREx = *FREy; *FREy = FREt;
	TCFt = *TCFx; *TCFx = *TCFy; *TCFy = TCFt;
	ACFt = *ACFx; *ACFx = *ACFy; *ACFy = ACFt;
	GAPOPt = *GAPOPx; *GAPOPx = *GAPOPy; *GAPOPy = GAPOPt;
}
