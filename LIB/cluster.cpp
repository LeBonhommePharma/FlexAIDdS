#include "gaboom.h"
#include "fileio.h"
#include "simd_distance.h"
#include "statmech.h"
#include "SoftBetaFreeEnergy.h"
#include "EnvFlags.h"
#include "ClusterRepMode.h"
#include "TargetServer.h"
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

void cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, target::TargetServer* ts, const std::string& ligand_name)
{
	bool Hungarian = false;
	int i,j;
	cfstr cf;                                /* complementarity function value */
	resid *res_ptr = NULL;
	cfstr* cf_ptr = NULL;

	FILE* outfile_ptr = NULL;

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
	
	// Clustering variable initialization.
	for(j=0;j<num_chrom;++j)
	{
		Clus_GAPOP[j]=-1;
		Clus_ACF[j]=0.0;
		Clus_TCF[j]=0.0;
		Clus_TOP[j]=0;
		Clus_FRE[j]=0;
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
						Clus_FRE[num_of_clusters]++;
					}
				}
			}
		}
		// Basin score: cluster-local soft-β free energy G̃ = H̃ − T·S̃ ≡ ACF.
		// Shared identity with BindingMode classic ranking and DatasetRunner S1
		// (LIB/SoftBetaFreeEnergy.h). Soft-β: β = 1/T (FA->beta), not 1/(k_B T).
		//
		// Election uses the DUPLICATE-INVARIANT strict variant. The legacy
		// soft_beta::acf() is marked "Diagnostic only — prefer
		// free_energy_strict() for claim re-ranking" by its own header
		// (SoftBetaFreeEnergy.h:135-136), yet it was driving the shipped
		// election. Because G̃ = Emin − T·ln Z with Z summed over MEMBERS,
		// multiplicity alone lowers G̃: at T=300 every population doubling buys
		// ~300·ln2 ≈ 208 CF units, which exceeds the entire CF spread on 72/85
		// Astex targets — i.e. the election degenerates into a popularity
		// contest (agrees with largest-cluster 82/85 but lowest-CF only 29/85).
		// free_energy_strict(UniqueGeometry) collapses exact-CF duplicates
		// before the same free energy, so cloned/re-emitted members can no
		// longer inflate a basin. Same units, same T; ranking only changes
		// where multiplicity was the deciding term.
		//
		// Set FLEXAIDDS_ELECT_LEGACY_ACF=1 to restore the legacy diagnostic
		// path bit-identically (A/B control against pre-fix baselines).
		Clus_TCF[num_of_clusters] = chrom[j].app_evalue;
		Clus_ACF[num_of_clusters] = chrom[j].app_evalue;
		if (FA->temperature > 0 && FA->beta > 0.0) {
			// Same env parser and same semantics as BindingMode.cpp, so the
			// two election paths cannot disagree about which arm is active.
			static const bool legacy_acf =
				flexaids::env_bool("FLEXAIDDS_ELECT_LEGACY_ACF");
			std::vector<double> member_energies;
			member_energies.reserve(static_cast<size_t>(num_chrom));
			for (int k = 0; k < num_chrom; ++k) {
				if (Clus_GAPOP[k] == j && std::isfinite(chrom[k].app_evalue))
					member_energies.push_back(chrom[k].app_evalue);
			}
			const double T_soft = static_cast<double>(FA->temperature);
			Clus_ACF[num_of_clusters] =
				legacy_acf
					? flexaids::soft_beta::acf(member_energies, T_soft)
					: flexaids::soft_beta::free_energy_strict(
					      member_energies, T_soft,
					      flexaids::soft_beta::StrictRerankMode::UniqueGeometry).G;
		}
		num_of_clusters++;

		// quit storing clusters up to N max results
		if(num_of_clusters == num_of_results){break;}
	}

	// Cap num_of_results to actual cluster count before sorting — arrays are
	// malloc'd to num_chrom elements; sorting past num_of_clusters overflows
	// when num_chrom < FA->max_results (e.g. entropy-collapsed run with 2 poses).
	if(num_of_clusters < num_of_results){num_of_results=num_of_clusters;}

	// ── Within-cluster representative election (P2) ────────────────────────────
	// FLEXAIDDS_CLUSTER_REP gate (see ClusterRepMode.h). This SUPERSEDES the
	// default-ON, Boltzmann-CF-weighted medoid committed at HEAD 3e674479c, which
	// violated two hard constraints: (1) it was DEFAULT-ON, so the unset default
	// no longer reproduced prior behavior; (2) it re-injected the CF signal
	// (Boltzmann weights) into a selector whose entire premise is CF-independence
	// (within-target Spearman(CF,RMSD) ≈ 0 — CF is orthogonal to pose correctness).
	//
	//   unset|lowcf → greedy lowest-CF head kept as-is (DEFAULT; bit-identical)
	//   medoid      → pure UNWEIGHTED geometric medoid (the fix; CF plays no role):
	//                    medoid = argmin_m  Σ_{n≠m} ‖x_m − x_n‖²      (≥3 members)
	//   bmedoid     → Boltzmann-CF-weighted medoid (HEAD variant; ablation only):
	//                    w_n = exp(−β·(E_n − E_min))/Z                 (≥2 members)
	//   center      → n/a for CF/leader backend (no density center) ⇒ treated as lowcf
	//
	// Reuses coord_cache (no geometry rebuilds). Never touches Clus_ACF — the
	// between-cluster ranking is representative-independent, so this only changes
	// WHICH member is emitted per cluster, not cluster order.
	const flexaids::ClusterRepMode rep_mode = flexaids::cluster_rep_mode();
	// Provenance keyed by chromosome index (survives the QuickSort_Clusters
	// permutation below): rep_shift_src[new_head] = old_head when a pick moved.
	std::vector<int> rep_shift_src(static_cast<size_t>(num_chrom), -1);

	if (rep_mode == flexaids::ClusterRepMode::MEDOID ||
	    rep_mode == flexaids::ClusterRepMode::BMEDOID)
	{
		const bool boltzmann  = (rep_mode == flexaids::ClusterRepMode::BMEDOID);
		// bmedoid weights need a temperature/β; the pure medoid is geometry-only.
		const bool weights_ok = !boltzmann || (FA->temperature > 0 && FA->beta > 0.0);
		// medoid: ≥3 members (a 2-member medoid is degenerate). bmedoid: ≥2 (HEAD).
		const int  min_members = boltzmann ? 2 : 3;

		if (weights_ok) {
			int n_refined = 0;
			for (int cl = 0; cl < num_of_clusters; ++cl) {
				if (Clus_FRE[cl] < min_members) continue;

				const int old_head = Clus_TOP[cl];

				// Members: Clus_GAPOP[k] == old_head for the head and every
				// chromosome absorbed into this cluster.
				std::vector<int> members;
				members.reserve(static_cast<size_t>(Clus_FRE[cl]));
				for (int k = 0; k < num_chrom; ++k) {
					if (Clus_GAPOP[k] == old_head)
						members.push_back(k);
				}
				if (static_cast<int>(members.size()) < min_members) continue;

				// Weights: uniform (pure medoid) or Boltzmann (bmedoid ablation).
				std::vector<double> weights(members.size(), 1.0);
				if (boltzmann) {
					double E_min = std::numeric_limits<double>::infinity();
					for (int k : members) {
						const double e = static_cast<double>(chrom[k].evalue);
						if (std::isfinite(e)) E_min = std::min(E_min, e);
					}
					if (!std::isfinite(E_min)) continue;
					double Z = 0.0;
					for (size_t mi = 0; mi < members.size(); ++mi) {
						const double e = static_cast<double>(chrom[members[mi]].evalue);
						const double w = std::isfinite(e)
						               ? std::exp(-FA->beta * (e - E_min)) : 0.0;
						weights[mi] = w;
						Z += w;
					}
					if (Z <= 0.0) continue;
					for (double& w : weights) w /= Z;
				}

				// argmin_candidate Σ_n w_n · sqRMSD(candidate, n). Squared units
				// (skip sqrt — monotone, argmin identical). Strict < keeps the
				// lowest array index on ties (deterministic, seed-independent).
				double best_cost   = std::numeric_limits<double>::max();
				int    best_member = old_head;
				for (size_t mi = 0; mi < members.size(); ++mi) {
					const float* cand = &coord_cache[members[mi] * coord_stride];
					double cost = 0.0;
					for (size_t ni = 0; ni < members.size(); ++ni) {
						if (ni == mi) continue;
						const float* other = &coord_cache[members[ni] * coord_stride];
						const float sq = flexaids::sum_sq_distances_f(cand, other, coord_stride);
						cost += weights[ni] * static_cast<double>(sq);
					}
					if (cost < best_cost) {
						best_cost   = cost;
						best_member = members[mi];
					}
				}

				if (best_member == old_head) continue;  // head already optimal

				// Remap Clus_GAPOP old_head→best_member (covers old_head itself);
				// update head tables. Leave Clus_ACF untouched (ranking invariant).
				for (int k = 0; k < num_chrom; ++k) {
					if (Clus_GAPOP[k] == old_head)
						Clus_GAPOP[k] = best_member;
				}
				Clus_GAPOP[best_member] = best_member;  // new head is self-referential
				Clus_TOP[cl] = best_member;
				Clus_TCF[cl] = chrom[best_member].app_evalue;
				rep_shift_src[best_member] = old_head;   // provenance for REMARK

				const double cost_angstrom2 = best_cost / static_cast<double>(nAtoms_clus);
				fprintf(stdout,
				        "[MEDOID_REFINE] mode=%s cluster %d: head %d→%d "
				        "(CF %.4f→%.4f, freq=%d, wRMSD²=%.4f Å²)\n",
				        flexaids::cluster_rep_mode_name(rep_mode),
				        cl, old_head, best_member,
				        chrom[old_head].evalue, chrom[best_member].evalue,
				        Clus_FRE[cl], cost_angstrom2);
				++n_refined;
			}
			if (n_refined > 0)
				fprintf(stdout,
				        "[MEDOID_REFINE] %d/%d clusters refined (mode=%s)\n",
				        n_refined, num_of_clusters,
				        flexaids::cluster_rep_mode_name(rep_mode));
		}
	}

	if(FA->temperature)
	{
		// Reordering the clusters properly by lowest ACF values first (after considering cluster's entropy !)
		// Classic FlexAID contract: this ACF order IS emission order when T>0.
		QuickSort_Clusters(Clus_TOP, Clus_FRE, Clus_TCF, Clus_ACF,
		                   0, num_of_results-1);
	}

	// ── Rank-0 emission policy ───────────────────────────────────────────────
	// Classic FlexAID (default when T>0 && !force_cf_rank_emission): keep ACF
	// order so dense entropy-favored basins can beat sparse lowest-CF clusters.
	// P3b rollback (force_cf_rank_emission or T==0): re-sort by representative
	// evalue so _0.pdb is always lowest CF (commit cd9004d behavior).
	// Single gate — flip FA->force_cf_rank_emission to restore old product path.
	const bool classic_entropy_emit =
		(FA->temperature > 0) && !FA->force_cf_rank_emission;
	if (!classic_entropy_emit)
	{
		for(int a=0; a<num_of_results-1; ++a)
		{
			int best_idx = a;
			for(int b=a+1; b<num_of_results; ++b)
			{
				if(chrom[Clus_TOP[b]].evalue < chrom[Clus_TOP[best_idx]].evalue)
					best_idx = b;
			}
			if(best_idx != a)
					swap_clusters(&Clus_TOP[a], &Clus_FRE[a], &Clus_TCF[a], &Clus_ACF[a],
					              &Clus_TOP[best_idx], &Clus_FRE[best_idx],
					              &Clus_TCF[best_idx], &Clus_ACF[best_idx]);
		}
	}
	else if (num_of_results > 0)
	{
		// Debuggable proof that rank-0 is ACF, not necessarily min CF.
		double best_cf = chrom[Clus_TOP[0]].evalue;
		for (int a = 1; a < num_of_results; ++a) {
			if (chrom[Clus_TOP[a]].evalue < best_cf)
				best_cf = chrom[Clus_TOP[a]].evalue;
		}
		fprintf(stdout,
			"[ENTROPY_RANK] classic FlexAID: rank-0 by ACF (ACF=%.4f CF=%.4f freq=%d); "
			"best_CF_among_emitted=%.4f (CF re-sort off; set force_cf_rank_emission to restore P3b)\n",
			Clus_ACF[0], chrom[Clus_TOP[0]].evalue, Clus_FRE[0], best_cf);
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
	fflush(stdout);
	
        // output results, 10% of the number of chromosomes or 
        // the number of clusters, the smallest.
      
	for(j=0;j<num_of_results;++j)
	{
		printf("emitting ranked pose %d/%d (TOP chrom=%d)...\n",
		       j + 1, num_of_results, Clus_TOP[j]);
		fflush(stdout);
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

		// Rebuild atom coordinates for PDB output and score that exact geometry.
		// A clash penalty is evidence about the emitted pose, never a reason to
		// substitute a stale search score.
		cf=ic2cf(FA,VC,atoms,residue,cleftgrid,GB->num_genes,FA->opt_par);
		const double emitted_cf = get_cf_evalue(&cf, FA);
		const double score_delta = std::abs(emitted_cf - chrom[Clus_TOP[j]].evalue);
		const bool score_pose_consistent = std::isfinite(emitted_cf) &&
		                                   score_delta <= 1e-4;
		if (!score_pose_consistent) {
			fprintf(stderr,
			        "WARNING: cluster %d stored CF=%.8f emitted-pose CF=%.8f "
			        "delta=%.8f\n",
			        j, chrom[Clus_TOP[j]].evalue, emitted_cf, score_delta);
		}

		size_t remark_len = 0;
		remark[0] = '\0';
		safe_remark_cat(remark, "REMARK optimized structure\n", &remark_len);

		snprintf(tmpremark, MAX_REMARK, "REMARK CF=%8.5f\n", emitted_cf);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.search=%8.5f\n",
		         chrom[Clus_TOP[j]].evalue);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.pose_score_delta=%.8f\n",
		         score_delta);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.pose_score_consistent=%s\n",
		         score_pose_consistent ? "true" : "false");
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.app=%8.5f\n",chrom[Clus_TOP[j]].app_evalue);
		safe_remark_cat(remark, tmpremark, &remark_len);

		// P2 provenance — emitted ONLY for non-default modes so the default
		// (lowcf) PDB stays byte-identical to pre-medoid HEAD (acceptance gate #1).
		if (rep_mode != flexaids::ClusterRepMode::LOWCF) {
			snprintf(tmpremark, MAX_REMARK, "REMARK cluster_rep_mode=%s\n",
			         flexaids::cluster_rep_mode_name(rep_mode));
			safe_remark_cat(remark, tmpremark, &remark_len);
			const int head_idx = Clus_TOP[j];
			if (head_idx >= 0 && head_idx < num_chrom && rep_shift_src[head_idx] >= 0) {
				snprintf(tmpremark, MAX_REMARK,
				         "REMARK cluster_rep_shifted=1 head_cf=%8.5f rep_cf=%8.5f\n",
				         chrom[rep_shift_src[head_idx]].evalue, chrom[head_idx].evalue);
				safe_remark_cat(remark, tmpremark, &remark_len);
			}
		}

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
		// CF-proxy ledger from cluster-member optimizer records (display/plugin
		// only; does not change GA ranking or cluster selection).
		{
			const double T = (FA->temperature > 0)
				? static_cast<double>(FA->temperature) : 300.0;
			statmech::StatMechEngine engine(
				T, statmech::make_contact_function_optimizer_provenance());
			engine.add_sample(static_cast<double>(chrom[Clus_TOP[j]].app_evalue));
			for (int k = 0; k < num_chrom; ++k) {
				if (k != Clus_TOP[j] && Clus_GAPOP[k] == Clus_TOP[j])
					engine.add_sample(static_cast<double>(chrom[k].app_evalue));
			}
			const statmech::Thermodynamics td = engine.compute();
			// Explicit metadata prevents legacy numeric keys from being promoted to
			// physical thermodynamic claims. Ranking remains soft_beta_G = ACF.
			snprintf(tmpremark, MAX_REMARK, "REMARK thermo_schema_version = 2\n");
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK thermo_claim_validity = proxy_only\n");
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK thermo_energy_domain = cf_arbitrary_units\n");
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK thermo_ensemble_measure = optimizer_samples\n");
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK thermo_reference_state = bound_only\n");
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK proxy_free_energy = %.6f\n", td.free_energy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			// Deprecated compatibility key; see thermo_claim_validity above.
			snprintf(tmpremark, MAX_REMARK, "REMARK free_energy = %.6f\n", td.free_energy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK soft_beta_G = %.6f\n", Clus_ACF[j]);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK enthalpy = %.6f\n", td.mean_energy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK entropy = %.8f\n", td.entropy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK heat_capacity = %.8f\n", td.heat_capacity);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK temperature = %.2f\n", T);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK binding_mode = %d\n", j);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK pose_rank = 1\n");
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK frequency = %d\n", Clus_FRE[j]);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}
		for(i=0;i<FA->npar;++i)
		{
			snprintf(tmpremark, MAX_REMARK, "REMARK [%8.3f]\n",FA->opt_par[i]);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}
		//snprintf(tmpremark, MAX_REMARK, "REMARK seed=%ld\n",FA->seed_ini);
		if(FA->refstructure == 1){
			const double rmsd_raw = calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian);
			snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n",
				rmsd_raw);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK rmsd_raw = %.5f\n", rmsd_raw);
			safe_remark_cat(remark, tmpremark, &remark_len);
			Hungarian = true;
			const double rmsd_sym = calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian);
			snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n",
				rmsd_sym);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK rmsd_sym = %.5f\n", rmsd_sym);
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

	// ── IP-5: opt-in full-population audit dump (.pop.tsv) ──────────────────────
	// Sibling to the fixed-width .rrd, adding the columns the .rrd lacks: CF
	// components (com/wal) and the elected-representative join (pose_id / is_elected)
	// so selection experiments become measurable ("was the near-native population
	// pose the one elected?"). Gated on refstructure==1 (native available for RMSD)
	// AND FLEXAIDDS_DUMP_POP=1 so DEFAULT artifacts (.pdb/.cad/.mcf/.rrd) are
	// byte-unchanged. The per-chromosome ic2cf re-score here is audit-only (never
	// on the benchmark hot path). Runs last: only frees follow, so mutating
	// FA->opt_par / FA->optres is safe.
	if (FA->refstructure == 1) {
		const char* dump_env = std::getenv("FLEXAIDDS_DUMP_POP");
		if (dump_env && std::atoi(dump_env) != 0) {
			std::vector<int> elected_rank(static_cast<size_t>(num_chrom), -1);
			for (int r = 0; r < num_of_results; ++r) {
				const int t = Clus_TOP[r];
				if (t >= 0 && t < num_chrom) elected_rank[t] = r;
			}
			char pop_path[MAX_PATH__];
			snprintf(pop_path, MAX_PATH__, "%s.pop.tsv", end_strfile);
			FILE* pf = fopen(pop_path, "w");
			if (pf) {
				fprintf(pf, "idx\tcluster\trmsd_to_head\trmsd_raw\trmsd_sym\t"
				            "cf_total\tcf_com\tcf_wal\tpose_id\tis_elected\n");
				bool Hung;
				for (int c = 0; c < num_chrom; ++c) {
					for (int g = 0; g < GB->num_genes; ++g)
						FA->opt_par[g] = chrom[c].genes[g].to_ic;
					Hung = false;
					const double rr = calc_rmsd(FA,atoms,residue,cleftgrid,
					                            FA->npar,FA->opt_par,Hung);
					Hung = true;
					const double rs = calc_rmsd(FA,atoms,residue,cleftgrid,
					                            FA->npar,FA->opt_par,Hung);
					cfstr cc = ic2cf(FA,VC,atoms,residue,cleftgrid,
					                 GB->num_genes,FA->opt_par);
					(void)cc;  // side effect: populates FA->optres[].cf components
					double com_sum = 0.0, wal_sum = 0.0;
					for (int r = 0; r < FA->num_optres; ++r) {
						com_sum += FA->optres[r].cf.com;
						wal_sum += FA->optres[r].cf.wal;
					}
					fprintf(pf, "%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\n",
					        c, Clus_GAPOP[c], Clus_RMSDT[c], rr, rs,
					        chrom[c].evalue, com_sum, wal_sum,
					        elected_rank[c], elected_rank[c] >= 0 ? 1 : 0);
				}
				fclose(pf);
				fprintf(stdout, "[DUMP_POP] wrote %s (%d population rows)\n",
				        pop_path, num_chrom);
			} else {
				fprintf(stderr, "WARNING: [DUMP_POP] could not open %s for write\n",
				        pop_path);
			}
		}
	}


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
void QuickSort_Clusters(int* TOP, int* FRE, double* TCF, double* ACF,
	                    int beg, int end)
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
			
			swap_clusters(&TOP[l], &FRE[l], &TCF[l], &ACF[l],
			              &TOP[r], &FRE[r], &TCF[r], &ACF[r]);
			
			if (p == r) p=l;
			++l;--r;
		}
		swap_clusters(&TOP[p], &FRE[p], &TCF[p], &ACF[p],
		              &TOP[r], &FRE[r], &TCF[r], &ACF[r]);
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_Clusters(TOP, FRE, TCF, ACF, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_Clusters(TOP, FRE, TCF, ACF, l, end);
			end = r;
		}
	}
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*                   Swap Function for Clusters                        */
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_clusters(int* TOPx, int* FREx, double* TCFx, double* ACFx,
	               int* TOPy, int* FREy, double* TCFy, double* ACFy)
{
	int TOPt, FREt;
	double TCFt, ACFt;
	TOPt = *TOPx; *TOPx = *TOPy; *TOPy = TOPt;
	FREt = *FREx; *FREx = *FREy; *FREy = FREt;
	TCFt = *TCFx; *TCFx = *TCFy; *TCFy = TCFt;
	ACFt = *ACFx; *ACFx = *ACFy; *ACFy = ACFt;
}
