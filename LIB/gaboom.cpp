#include "gaboom.h"
#include "Vcontacts.h"
#include "fileio.h"
#include "flexaid_exception.h"
#include "ga_constants.h"
#include "hardware_dispatch.h"
#include "MIFGrid.h"
#include "CavityDetect/SpatialGrid.h"

#include <random>
#include <functional>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <array>
#include <memory>
#include <span>
#include <unordered_set>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Eigen/Dense>
#include <cmath>

#ifdef FLEXAIDS_USE_CUDA
#include "cuda_eval.cuh"
#endif

#ifdef FLEXAIDS_USE_METAL
#include "metal_eval.h"
#endif

#include "statmech.h"
#include "tENCoM/tencm.h"
#include "ShannonThermoStack/ShannonThermoStack.h"
#include "ga_diversity.h"
#include "TurboQuant.h"
#include "GAContext.h"
#include "GPUContextPool.h"
#include "fast_optics.hpp"
#include "NATURaL/NATURaLDualAssembly.h"

// in milliseconds
# define SLEEP GA_SLEEP_MS

#ifdef _WIN32
# include <windows.h>
#else
# include <unistd.h>
#endif


/// ═══ CCBM: Add receptor conformer strain energy to chromosome evalue ═══
/// When multi-model mode is ON, the model gene selects the receptor conformer
/// and the strain energy of that conformer is added to the CF-based evalue.
/// This makes the GA search the joint (ligand_pose, receptor_conformer) space.
static inline void ccbm_inject_strain(FA_Global* FA, chromosome& chrom, const genlim* gene_lim) {
    if (!FA->multi_model || FA->n_models <= 1 || FA->model_gene_index < 0) return;
    int mg = FA->model_gene_index;
    // Decode discrete model index from the gene value (round to nearest integer)
    int model_idx = static_cast<int>(std::round(chrom.genes[mg].to_ic));
    // Clamp to valid range
    if (model_idx < 0) model_idx = 0;
    if (model_idx >= FA->n_models) model_idx = FA->n_models - 1;
    // Snap gene IC value to exact integer for discrete gene
    chrom.genes[mg].to_ic = static_cast<double>(model_idx);
    // Add strain energy
    double strain = FA->model_strain[model_idx];
    chrom.evalue += strain;
    chrom.app_evalue += strain;
}

// Forward declarations for functions defined later in this file
int reproduce(FA_Global* FA,GB_Global* GB,VC_Global* VC, chromosome* chrom, const genlim* gene_lim,
               atom* atoms,resid* residue,gridpoint* cleftgrid,char* repmodel,
               double mutprob, double crossprob, int print,
               std::function<int32_t()> & dice,
               std::unordered_map<size_t, int> & duplicates,
               cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
               GAContext& ctx);

void calculate_fitness(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,
                       atom* atoms,resid* residue,gridpoint* cleftgrid,char method[], int pop_size, int print,
                       cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
                       GAContext& ctx);

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int GA(FA_Global* FA, GB_Global* GB,VC_Global* VC,chromosome** chrom,chromosome** chrom_snapshot,
       genlim** gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,char gainpfile[],
       int* memchrom, cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
       GAContext* ctx){

	// Create a stack-local context if none was provided (backward compat)
	GAContext local_ctx;
	if (!ctx) ctx = &local_ctx;

	int i;
	int print=0;

	//char tmp_rrgfile[MAX_PATH__];
	//int rrg_flag;
	//int rrg_skip=100;

	char outfile[MAX_PATH__];
	int n_chrom_snapshot=0;
	char gridfile[MAX_PATH__];
	char gridfilename[MAX_PATH__];

	int geninterval=GA_DEFAULT_GEN_INTERVAL;
	int popszpartition=GA_DEFAULT_POP_PARTITION;

	int  state=0;
	char PAUSEFILE[MAX_PATH__];
	char ABORTFILE[MAX_PATH__];
	char STOPFILE[MAX_PATH__];

	const int INTERVAL = GA_STATE_CHECK_INTERVAL; // sleep interval between checking file state

	*memchrom=0; //num chrom allocated in memory

	// for generation random doubles from [0,1[ (mutation crossover operators)
#ifdef _WIN32
	snprintf(PAUSEFILE,MAX_PATH__,"%s\\.pause",FA->state_path);
	snprintf(ABORTFILE,MAX_PATH__,"%s\\.abort",FA->state_path);
	snprintf(STOPFILE,MAX_PATH__,"%s\\.stop",FA->state_path);
#else
	snprintf(PAUSEFILE,MAX_PATH__,"%s/.pause",FA->state_path);
	snprintf(ABORTFILE,MAX_PATH__,"%s/.abort",FA->state_path);
	snprintf(STOPFILE,MAX_PATH__,"%s/.stop",FA->state_path);
#endif

	GB->num_genes=FA->npar;
	if(GB->num_genes == 0){
		fprintf(stderr,"ERROR: no parameters to optimize.\n");
		Terminate(1);
	}

	// ═══ CCBM: Add discrete gene for receptor model index when multi-model is ON ═══
	if (FA->multi_model && FA->n_models > 1) {
		FA->model_gene_index = GB->num_genes;  // last gene is the model selector
		GB->num_genes++;  // add one gene for receptor conformer selection
		printf("CCBM: multi-model mode enabled with %d conformers, model_gene_index=%d\n",
		       FA->n_models, FA->model_gene_index);
	} else {
		FA->model_gene_index = -1;  // no model gene
	}

	printf("num_genes=%d\n",GB->num_genes);

	printf("file in GA is <%s>\n",gainpfile);

	if (gainpfile[0] != '\0') {
		//GB->rrg_skip=0;
		GB->adaptive_ga=0;
		GB->num_print=GA_DEFAULT_NUM_PRINT;
		GB->print_int=GA_DEFAULT_PRINT_INT;
		GB->seed = GA_DEFAULT_SEED;

	// Entropy convergence defaults (opt-in)
	GB->entropy_convergence    = 0;
	GB->entropy_check_interval = GA_DEFAULT_ENTROPY_CHECK_INTERVAL;
	GB->entropy_window         = GA_DEFAULT_ENTROPY_WINDOW;
	GB->entropy_rel_threshold  = GA_DEFAULT_ENTROPY_REL_THRESHOLD;

	printf("file in GA is <%s>\n",gainpfile);

		read_gainputs(FA,GB,&geninterval,&popszpartition,gainpfile);
	} else {
		printf("No GA input file — using pre-configured parameters\n");
	}
	unsigned int tt;
	if (GB->seed==0)
	{
		tt = static_cast<unsigned int>(time(0));
	}
	else
	{
		tt = GB->seed;
	}
	//tt = (unsigned)1;
	printf("srand=%u\n", tt);
	srand(tt);
	std::mt19937 rng(tt);

	std::uniform_int_distribution<int32_t> one_to_max_int32( 0, MAX_RANDOM_VALUE );
	std::function<int32_t()> dice = [&]() { return one_to_max_int32(rng); };

	(*gene_lim) = (genlim*)malloc(GB->num_genes*sizeof(genlim));
	if(!(*gene_lim)){
		fprintf(stderr,"ERROR: memory allocation error for gene_lim.\n");
		Terminate(2);
	}

	long int at = 0;

	if(strcmp(GB->pop_init_method,"RANDOM") == 0){
		set_gene_lim(FA, GB, (*gene_lim));
		// ═══ CCBM: Set gene limits for model selection gene ═══
		if (FA->multi_model && FA->n_models > 1 && FA->model_gene_index >= 0) {
			int mg = FA->model_gene_index;
			(*gene_lim)[mg].min = 0.0;
			(*gene_lim)[mg].max = static_cast<double>(FA->n_models - 1);
			(*gene_lim)[mg].del = 1.0;  // discrete steps
			(*gene_lim)[mg].map = -1;   // no mapping
			printf("CCBM: model gene %d: min=0 max=%d delta=1 (discrete)\n",
			       mg, FA->n_models - 1);
		}
		set_bins((*gene_lim),GB->num_genes);

	}else if(strcmp(GB->pop_init_method,"IPFILE") == 0){
		at = read_pop_init_file(FA, GB, (*gene_lim), GB->pop_init_file);
		if(!at){
			fprintf(stderr,"ERROR: Unknown format for pop init file.\n");
			Terminate(10);
		}
	}

	if(GB->print_int < 0){ GB->print_int = 1; }

	//if(GB->rrg_skip > 0){ rrg_skip = GB->rrg_skip; }

	if(GB->num_print > GB->num_chrom){ GB->num_print = GB->num_chrom; }

	if(popszpartition > GB->num_chrom){ popszpartition = GB->num_chrom; }

	if(FA->opt_grid){
		printf("will partition grid every %d generations considering %d individuals\n",
		       geninterval, popszpartition);
	}

	validate_dups(GB, (*gene_lim), GB->num_genes);

	(*memchrom) = GB->num_chrom;
	if(strcmp(GB->rep_model,"STEADY")==0){
		(*memchrom) += GB->ssnum;
	}else if(strcmp(GB->rep_model,"BOOM")==0){
		(*memchrom) += (int)(GB->pbfrac*(double)GB->num_chrom);
	}

	//printf("memchrom=%d\n",(*memchrom));
	//printf("num_genes=%d\n",GB->num_genes);

	// *** chrom
	(*chrom) = (chromosome*)malloc((*memchrom)*sizeof(chromosome));
	if(!(*chrom)){
		fprintf(stderr,"ERROR: memory allocation error for chrom.\n");
		Terminate(2);
	}

	for(i=0;i<(*memchrom);++i)
	{
		(*chrom)[i].genes = (gene*)malloc(GB->num_genes*sizeof(gene));

		if(!(*chrom)[i].genes){
			fprintf(stderr,"ERROR: memory allocation error for chrom[%d].genes.\n",i);
			Terminate(2);
		}

		(*chrom)[i].app_evalue = 0.0;
		(*chrom)[i].evalue = 0.0;
		(*chrom)[i].fitnes = 0.0;
		(*chrom)[i].boltzmann_weight = 0.0;
		(*chrom)[i].free_energy = 0.0;
		(*chrom)[i].status = ' ';
	}

	// *** chrom_snapshot
	(*chrom_snapshot) = (chromosome*)malloc((GB->num_chrom*GB->max_generations)*sizeof(chromosome));
	if(!(*chrom_snapshot))
	{
		fprintf(stderr,"ERROR: memory allocation error for chrom_snapshot.\n");
		Terminate(2);
	}

	for(i=0;i<(GB->num_chrom*GB->max_generations);++i)
	{
		(*chrom_snapshot)[i].genes = (gene*)malloc(GB->num_genes*sizeof(gene));

		if(!(*chrom_snapshot)[i].genes){
			fprintf(stderr,"ERROR: memory allocation error for chrom_snapshot[%d].genes.\n",i);
			Terminate(2);
		}

		(*chrom_snapshot)[i].app_evalue = 0.0;
		(*chrom_snapshot)[i].evalue = 0.0;
		(*chrom_snapshot)[i].fitnes = 0.0;
		(*chrom_snapshot)[i].boltzmann_weight = 0.0;
		(*chrom_snapshot)[i].free_energy = 0.0;
		(*chrom_snapshot)[i].status = ' ';
		//printf("chrom_snapshot[%d] allocated at address %p!\n", i, &(*chrom_snapshot)[i]);
	}

	printf("alpha %lf peaks %lf scale %lf\n",GB->alpha,GB->peaks,GB->scale);
	GB->sig_share=0.0;

	for(i=0;i<GB->num_genes;i++)
	{
		//printf("max=%6.3f min=%6.3f del=%6.3f\n",(*gene_lim)[i].max,(*gene_lim)[i].min,(*gene_lim)[i].del);
		//PAUSE;
		GB->sig_share += ((*gene_lim)[i].max-(*gene_lim)[i].min)*((*gene_lim)[i].max-(*gene_lim)[i].min);
	}
	GB->sig_share = sqrt(GB->sig_share/(double)GB->num_genes)/(2.0*pow(GB->peaks,(1.0/(double)GB->num_genes)));
	GB->sig_share /= GB->scale;
	printf("SIGMA_SHARE=%f\n",GB->sig_share);
	fflush(stdout);

	// for(i=0;i<GB->num_genes;i++) {
	//printf("GENE(%d)=[%8.3f,%8.3f,%8.3f,%d]\n",
	//	   i,(*gene_lim)[i].min,(*gene_lim)[i].max,(*gene_lim)[i].del);
	//PAUSE;

	std::unordered_map<size_t, int> duplicates;

	populate_chromosomes(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),
			     GB->pop_init_method,target,GB->pop_init_file,at,0,print,dice,duplicates);
	//}

	//print_pop((*chrom),(*gene_lim),GB->num_chrom,GB->num_genes);

	/*
	  for(i=0;i<GB->num_genes;i++){
	  printf("%d %f %f %f\n",i,GB->min_opt_par[i],GB->max_opt_par[i],GB->del_opt_par[i]);
	  PAUSE;
	  }
	*/

	int save_num_chrom = (int)(GB->num_chrom*SAVE_CHROM_FRACTION);
	int nrejected = 0;

	// Entropy convergence tracking
	std::vector<double> entropy_history;
	bool entropy_converged = false;
	if (GB->entropy_convergence) {
		entropy_history.reserve(GB->max_generations / GB->entropy_check_interval + 1);
	}

	////////////////////////////////
	// Stagnation detection: terminate GA when best fitness stops improving
	const int STAGNATION_WINDOW = 100;   // check every N generations
	const int STAGNATION_LIMIT  = 300;   // break after this many stagnant windows
	double prev_best_fitness = -1e30;
	int    stagnation_count  = 0;
	bool   ga_stagnant = false;

	////// Genetic Algorithm ///////
	////////////////////////////////
	for(i=0;i<GB->max_generations;i++)
	{
		///////////////////////////////////////////////////

		state=check_state(PAUSEFILE,ABORTFILE,STOPFILE,INTERVAL);

		if(state == -1){
			return(state);
		}else if(state == 1){
			break;
		}

		////////////////////////////////

		////////////////////////////////

		////////////////////////////////

		//printf("chrom_snapshot[%d] at address %p\n", i*GB->num_chrom, chrom_snapshot[i*GB->num_chrom]);
		if (	FA->opt_grid                    &&     // if a OPTGRD line was specified
		    	((i+1) % geninterval) == 0      &&     // is factor of
		    	(i+1) != GB->max_generations 	)      // discard the last generation
		{

			//need to sort in decreasing order of energy
			QuickSort((*chrom),0,GB->num_chrom-1,true);

			//printf("Partionning grid...(%d)\n",FA->popszpartition);
			partition_grid(FA,(*chrom),(*gene_lim),atoms,residue,cleftgrid,popszpartition,1);

			if(FA->output_range){
#ifdef _WIN32
				snprintf(gridfile,MAX_PATH__,"%s\\grid.%d.prt.pdb",FA->temp_path,i+1);
#else
				snprintf(gridfile,MAX_PATH__,"%s/grid.%d.prt.pdb",FA->temp_path,i+1);
#endif
				write_grid(FA,(*cleftgrid),gridfile);
			}

			slice_grid(FA,(*gene_lim),atoms,residue,cleftgrid);

			if(FA->output_range){
#ifdef _WIN32
				snprintf(gridfile,MAX_PATH__,"%s\\grid.%d.slc.pdb",FA->temp_path,i+1);
#else
				snprintf(gridfile,MAX_PATH__,"%s/grid.%d.slc.pdb",FA->temp_path,i+1);
#endif
				write_grid(FA,(*cleftgrid),gridfile);
			}

			// Recompute MIF for adapted grid
			if (FA->mif_enabled || FA->grid_prio_percent < 100.0f) {
				std::vector<atom> protein_atoms(atoms, atoms + FA->atm_cnt_real);
				cavity_detect::SpatialGrid sg;
				sg.build(protein_atoms);
				auto mif = mif::compute_mif(*cleftgrid, FA->num_grd,
				                             atoms, FA->atm_cnt_real, sg);
				free(FA->mif_energies); free(FA->mif_sorted); free(FA->mif_cdf);
				FA->mif_count = static_cast<int>(mif.sorted_indices.size());
				FA->mif_energies = static_cast<float*>(
				    malloc(mif.energies.size() * sizeof(float)));
				FA->mif_sorted = static_cast<int*>(
				    malloc(mif.sorted_indices.size() * sizeof(int)));
				std::copy_n(mif.energies.data(), mif.energies.size(), FA->mif_energies);
				std::copy_n(mif.sorted_indices.data(), mif.sorted_indices.size(), FA->mif_sorted);
				mif::build_sampling_cdf(mif, FA->mif_temperature);
				FA->mif_cdf = static_cast<double*>(
				    malloc(mif.cdf.size() * sizeof(double)));
				std::copy_n(mif.cdf.data(), mif.cdf.size(), FA->mif_cdf);
			}

			validate_dups(GB, (*gene_lim), GB->num_genes);

			//repopulate unselected individuals
			populate_chromosomes(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),
					     GB->pop_init_method,target,GB->pop_init_file,at,popszpartition,print,dice,duplicates);
		}

		print = ( (i+1) % GB->print_int == 0 ) ? 1 : 0;
		//if(print) { printf("Generation: %5d\n",i+1); }

		//print_par(chrom,gene_lim,20,GB->num_genes);
		//PAUSE;

		/*
		  rrg_flag=0;
		  if((i/rrg_skip)*rrg_skip == i) rrg_flag=1;
		  if((rrg_flag==1) && (GB->outgen==1)){
		  if(FA->refstructure == 1){
		  snprintf(tmp_rrgfile,MAX_PATH__,"%s_%d.rrg",FA->rrgfile,i);
		  //printf("%s\n",tmp_rrgfile);
		  //PAUSE;
		  write_rrg(FA,GB,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),tmp_rrgfile);
		  }
		  }
		*/


		//before reproducing for an extra generation, evaluate if population has converged.
		//before calculating get avg and max fitness of the whole pop.
		fitness_stats(GB,(*chrom),GB->num_chrom);

		//printf("------fitness stats-------\navg=%8.3f\tmax=%8.3f\n",GB->fit_avg,GB->fit_max);
        //getchar();

		// Stagnation detection: check if best fitness has plateaued
		if ((i + 1) % STAGNATION_WINDOW == 0 && i > 0) {
			if (std::abs(GB->fit_max - prev_best_fitness) < 1e-6) {
				stagnation_count += STAGNATION_WINDOW;
				if (stagnation_count >= STAGNATION_LIMIT) {
					printf("GA terminated early: fitness stagnant for %d generations (best=%.4f)\n", stagnation_count, GB->fit_max);
					ga_stagnant = true;
					break;
				}
			} else {
				stagnation_count = 0;
			}
			prev_best_fitness = GB->fit_max;
		}

		// Entropy convergence check (opt-in via ENTRCNVG config keyword)
		if (GB->entropy_convergence &&
		    ((i + 1) % GB->entropy_check_interval == 0)) {
			std::vector<double> pop_energies(GB->num_chrom);
			for (int c = 0; c < GB->num_chrom; ++c)
				pop_energies[c] = (*chrom)[c].evalue;
			double H = shannon_thermo::compute_shannon_entropy(
				pop_energies, shannon_thermo::DEFAULT_HIST_BINS);
			entropy_history.push_back(H);

			if (shannon_thermo::detect_entropy_plateau(
			        entropy_history, GB->entropy_window,
			        GB->entropy_rel_threshold)) {
				printf("Entropy convergence at generation %d "
				       "(H=%.4f nats, stable for %d checks)\n",
				       i + 1, H, GB->entropy_window);
				entropy_converged = true;
				break;
			}
		}

		// Diversity monitoring: detect and mitigate premature entropy collapse
		if (GB->diversity_monitoring &&
		    ((i + 1) % GB->diversity_check_interval == 0)) {
			auto dm = ga_diversity::compute_diversity(
				*chrom, GB->num_chrom, GB->num_genes, *gene_lim,
				GB->diversity_collapse_threshold);
			if (dm.collapse_detected && (i + 1) < GB->max_generations / 2) {
				// Only trigger catastrophic mutation in the first half of generations
				ga_diversity::catastrophic_mutation(
					*chrom, GB->num_chrom, GB->num_genes, *gene_lim,
					GB->catastrophic_mutation_fraction, rng);
				GB->catastrophic_mutation_count++;
				fprintf(stderr, "[DIVERSITY] Catastrophic mutation #%d at gen %d "
				        "(H_allele=%.3f, min_gene=%.3f)\n",
				        GB->catastrophic_mutation_count, i + 1,
				        dm.allele_entropy, dm.min_gene_entropy);
			}
		}

		nrejected = reproduce(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),
				      GB->rep_model,GB->mut_rate,GB->cross_rate,print,dice,duplicates,target,*ctx);

		save_snapshot(&(*chrom_snapshot)[i*GB->num_chrom],(*chrom),save_num_chrom,GB->num_genes);
		n_chrom_snapshot += save_num_chrom;


		if(strcmp(GB->fitness_model,"PSHARE")==0){
			QuickSort((*chrom),0,GB->num_chrom-1,false);

			if(print){
				printf("best by fitnes\n");
				print_par((*chrom),(*gene_lim),GB->num_print,GB->num_genes, stdout);
			}
		}

	}

	printf("%d ligand conformers rejected\n", nrejected);
	if (entropy_converged)
		printf("GA terminated early by entropy convergence\n");
		if (ga_stagnant)
			printf("GA terminated early by fitness stagnation\n");

	QuickSort((*chrom),0,GB->num_chrom-1,true);

	snprintf(outfile,MAX_PATH__,"%s_par.res",FA->rrgfile);
	if (FA->htpmode == false) {write_par((*chrom),(*gene_lim),i+1,outfile,GB->num_chrom,GB->num_genes);}

	printf("sorting chrom_snapshot\n");
	//quicksort_app_evalue((*chrom_snapshot),0,n_chrom_snapshot-1);
	QuickSort((*chrom_snapshot),0,n_chrom_snapshot-1,true);

	/*
	  printf("Save snapshot == END ==\n");
	  print_par((*chrom_snapshot),(*gene_lim),n_chrom_snapshot,GB->num_genes);
	*/

	printf("removing duplicates\n");
	n_chrom_snapshot = remove_dups((*chrom_snapshot),n_chrom_snapshot,GB->num_genes);

	/*
		printf("Save snapshot == END ==\n");
		print_par((*chrom_snapshot),(*gene_lim),n_chrom_snapshot,GB->num_genes);
	*/

	// Thermodynamic analysis of the final conformational ensemble
	if(n_chrom_snapshot > 0) {
		double T_K = (FA->temperature > 0) ? static_cast<double>(FA->temperature) : GA_DEFAULT_TEMPERATURE_K;
		statmech::StatMechEngine sme(T_K);
		for(int s = 0; s < n_chrom_snapshot; ++s)
			sme.add_sample((*chrom_snapshot)[s].evalue);

		// Optional super-cluster pre-filtering for faster Shannon entropy collapse
		if (FA->use_super_cluster && n_chrom_snapshot > 4) {
			std::vector<fast_optics::Point> energy_pts(n_chrom_snapshot);
			for (int s = 0; s < n_chrom_snapshot; ++s)
				energy_pts[s].coords = { (*chrom_snapshot)[s].evalue };

			fast_optics::FastOPTICS foptics(energy_pts, std::max(GA_FOPTICS_MIN_POINTS, n_chrom_snapshot / GA_FOPTICS_DIVISOR));
			auto sc_indices = foptics.extractSuperCluster(fast_optics::ClusterMode::SUPER_CLUSTER_ONLY);

			if (!sc_indices.empty() && sc_indices.size() < static_cast<size_t>(n_chrom_snapshot)) {
				statmech::StatMechEngine sme_filtered(T_K);
				for (size_t idx : sc_indices)
					sme_filtered.add_sample((*chrom_snapshot)[idx].evalue);

				printf("--- SuperCluster pre-filter: %zu / %d poses selected ---\n",
				       sc_indices.size(), n_chrom_snapshot);
				sme = sme_filtered;
			}
		}

		statmech::Thermodynamics td = sme.compute();
		printf("--- Thermodynamics (T = %.1f K, N = %d conformers) ---\n",
		       td.temperature, n_chrom_snapshot);
		printf("  Helmholtz free energy  F  = %10.4f kcal/mol\n", td.free_energy);
		printf("  Mean energy          <E>  = %10.4f kcal/mol\n", td.mean_energy);
		printf("  Energy std dev        σ_E = %10.4f kcal/mol\n", td.std_energy);
		printf("  Heat capacity         C_v = %10.4f kcal/(mol·K)\n", td.heat_capacity);
		printf("  Entropy (conf)        S   = %10.6f kcal/(mol·K)\n", td.entropy);

		// ── Phase 2.5: TurboQuant ensemble compression ──────────────
		// Quantize the conformational ensemble energy vectors using TurboQuant
		// (Zandieh et al. 2025, arXiv:2504.19874) for near-optimal distortion.
		// This compresses the population energy representations while preserving
		// inner product structure needed for Boltzmann-weighted Shannon entropy.
		//
		// When TQENS (use_tqens) is enabled, we use QuantizedEnsemble with a
		// multi-dimensional energy descriptor (com, wal, sas, elec) per conformer
		// and compute the approximate partition function via unbiased inner-product
		// preserving TurboQuantProd quantization.  We compare to exact StatMechEngine
		// Boltzmann weights and log the empirical bias and max weight error.
		//
		// TurboQuant MSE bound: D_mse ≤ sqrt(3π/2) · 1/4^b ≈ 2.7/4^b
		// At b=3 bits/coordinate: D_mse ≈ 0.03 (97% fidelity)
		if (n_chrom_snapshot > GA_TQENS_MIN_SNAPSHOTS) {
			constexpr int TQ_BITS = GA_TQENS_BITS;  // 3 bits/coord → 97% fidelity, 10.7× compression

			if (FA->use_tqens) {
				// ── Multi-dimensional QuantizedEnsemble (full TurboQuantProd) ──
				// Energy descriptor: (com, wal, sas, elec) → 4 dimensions
				// Each chromosome's cfstr provides these component values.
				constexpr int TQ_EDIM = GA_TQENS_ENERGY_DIM;
				turboquant::QuantizedEnsemble qens(TQ_EDIM, TQ_BITS);
				qens.reserve(n_chrom_snapshot);

				// Build energy descriptor vectors from cfstr components
				std::vector<std::array<float, TQ_EDIM>> descriptors(n_chrom_snapshot);
				for (int s = 0; s < n_chrom_snapshot; ++s) {
					const cfstr& cf = (*chrom_snapshot)[s].cf;
					descriptors[s][0] = static_cast<float>(cf.com);
					descriptors[s][1] = static_cast<float>(cf.wal);
					descriptors[s][2] = static_cast<float>(cf.sas);
					descriptors[s][3] = static_cast<float>(cf.elec);
					qens.add_state(std::span<const float>(descriptors[s].data(), TQ_EDIM));
				}

				// Construct beta_E vector: β times a unit energy-weighting direction
				// For the Boltzmann partition function Z = Σ exp(-β·E_total),
				// E_total = com + wal (standard CF).  We project via beta_E = β·(1,1,0,0)
				// so that ⟨beta_E, descriptor⟩ = β·(com + wal) = β·E_total.
				float beta_val = static_cast<float>(1.0 / (statmech::kB_kcal * T_K));
				std::array<float, TQ_EDIM> beta_E = {beta_val, beta_val, 0.0f, 0.0f};

				// Compute approximate partition function via QuantizedEnsemble
				std::vector<float> approx_weights(n_chrom_snapshot);
				float log_Z_approx = qens.compute_partition_function(
					std::span<const float>(beta_E.data(), TQ_EDIM),
					std::span<float>(approx_weights));

				// Compute exact Boltzmann weights from StatMechEngine for comparison
				std::vector<double> exact_bw = sme.boltzmann_weights();

				// Compare approximate vs exact weights
				double sum_bias = 0.0, max_weight_err = 0.0;
				for (int s = 0; s < n_chrom_snapshot; ++s) {
					double err = static_cast<double>(approx_weights[s]) - exact_bw[s];
					sum_bias += err;
					max_weight_err = std::max(max_weight_err, std::abs(err));
				}
				double mean_bias = sum_bias / n_chrom_snapshot;

				// Compute exact log(Z) for comparison
				statmech::Thermodynamics td_exact = sme.compute();
				double log_Z_exact = td_exact.log_Z;
				double pf_err = std::abs(static_cast<double>(log_Z_approx) - log_Z_exact);

				printf("--- TurboQuant QuantizedEnsemble (b=%d, d=%d) ---\n", TQ_BITS, TQ_EDIM);
				printf("  Conformers             N   = %d\n", n_chrom_snapshot);
				printf("  Energy descriptor dim  d   = %d (com, wal, sas, elec)\n", TQ_EDIM);
				printf("  Memory (quantized)         = %zu bytes\n", qens.memory_bytes());
				printf("  Memory (raw float)         = %zu bytes\n",
				       static_cast<size_t>(n_chrom_snapshot) * TQ_EDIM * sizeof(float));
				printf("  Mean Boltzmann weight bias = %+.6e\n", mean_bias);
				printf("  Max  Boltzmann weight err  = %.6e\n", max_weight_err);
				printf("  log(Z) exact               = %.6f\n", log_Z_exact);
				printf("  log(Z) approx              = %.6f\n", static_cast<double>(log_Z_approx));
				printf("  |Δlog(Z)|                  = %.6e\n", pf_err);
			} else {
				// ── Legacy scalar-only diagnostic (d=1, skip TurboQuant which requires d>=2) ──
				constexpr int TQ_DIM = 1;
				size_t raw_bytes = n_chrom_snapshot * sizeof(float);
				size_t quant_bytes = n_chrom_snapshot * ((TQ_DIM * TQ_BITS + 7) / 8 + sizeof(float));
				printf("--- TurboQuant ensemble compression (b=%d, d=%d) ---\n", TQ_BITS, TQ_DIM);
				printf("  Conformers             N   = %d\n", n_chrom_snapshot);
				printf("  Raw size                   = %zu bytes\n", raw_bytes);
				printf("  Quantized size             = %zu bytes\n", quant_bytes);
				printf("  Compression ratio          = %.1f×\n",
				       static_cast<double>(raw_bytes) / quant_bytes);
			}
		}

		// ── Phase 3: TorsionalENM vibrational entropy ────────────────
		tencm::TorsionalENM tencm_model;
		if (FA->is_protein && FA->res_cnt > GA_TENCM_MIN_RESIDUES) {
			tencm_model.build(atoms, residue, FA->res_cnt);
			if (tencm_model.is_built()) {
				// Store mode count on FA for BindingMode vibrational correction
				FA->normal_modes = static_cast<int>(tencm_model.modes().size());

				// Run full ShannonThermoStack: Shannon conf entropy + torsional vib entropy
				shannon_thermo::FullThermoResult ftr =
					shannon_thermo::run_shannon_thermo_stack(
						sme, tencm_model, td.free_energy, T_K);

				printf("--- ShannonThermoStack (vibrational entropy integration) ---\n");
				printf("  Shannon conf entropy    = %10.4f nats\n", ftr.shannonEntropy);
				printf("  Torsional vib entropy   = %10.6f kcal/(mol·K)\n", ftr.torsionalVibEntropy);
				printf("  Entropy contribution    = %10.4f kcal/mol (-TΔS)\n", ftr.entropyContribution);
				printf("  Total ΔG (F + vib corr) = %10.4f kcal/mol\n", ftr.deltaG);
			}
		}
	}

	// NATURaL co-translational / co-transcriptional DualAssembly analysis
	// Skipped when --folded flag or advanced.assume_folded=true (receptor is fully folded)
	if (!FA->assume_folded && FA->resligand && FA->resligand->fatm && FA->resligand->latm) {
		int lig_start   = FA->resligand->fatm[0];
		int lig_end     = FA->resligand->latm[0];
		int n_lig_atoms = lig_end - lig_start + 1;
		if (n_lig_atoms > 0 && FA->MIN_NUM_RESIDUE > 0) {
			natural::NATURaLConfig ncfg = natural::auto_configure(
				&atoms[lig_start], n_lig_atoms,
				residue, FA->MIN_NUM_RESIDUE);
			if (ncfg.enabled) {
				ncfg.temperature_K = (FA->temperature > 0)
				                     ? static_cast<double>(FA->temperature)
				                     : GA_NATURAL_DEFAULT_TEMP;
				natural::DualAssemblyEngine engine(
					ncfg, FA, VC, atoms, residue, FA->MIN_NUM_RESIDUE);
				auto trajectory = engine.run();
				printf("--- NATURaL Co-translational DualAssembly (%zu growth steps) ---\n",
				       trajectory.size());
				if (!trajectory.empty()) {
					printf("  Final ΔG (co-translational) = %10.4f kcal/mol\n",
					       engine.final_deltaG());
					FA->natural_deltaG = engine.final_deltaG();

					int  n_pause        = 0;
					int  n_tm           = 0;
					int  n_burst_events = 0;
					int  max_burst_size = 0;
					int  n_nuc_seeds    = 0;
					std::unordered_set<int> seen_bursts, seen_seeds;

					for (const auto& step : trajectory) {
						if (step.is_pause_site) ++n_pause;
						if (step.tm_inserted)   ++n_tm;
						if (step.burst_unit_id >= 0 &&
						    seen_bursts.insert(step.burst_unit_id).second) {
							++n_burst_events;
							if (step.burst_size > max_burst_size)
								max_burst_size = step.burst_size;
						}
						if (step.nucleation_seed_id >= 0 &&
						    seen_seeds.insert(step.nucleation_seed_id).second)
							++n_nuc_seeds;
					}
					printf("  Pause sites detected        = %d\n",    n_pause);
					printf("  TM insertions               = %d\n",    n_tm);
					printf("  Burst elongation events     = %d (max %d residues/burst)\n",
					       n_burst_events, max_burst_size);
					printf("  Nucleation seeds detected   = %d\n",    n_nuc_seeds);
				}
			}
		}
	}

	return n_chrom_snapshot;
}

void copy_chrom(chromosome* dest, const chromosome* src, int num_genes){

	dest->cf = src->cf;
	dest->evalue = src->evalue;
	dest->app_evalue = src->app_evalue;
	dest->fitnes = src->fitnes;
	dest->boltzmann_weight = src->boltzmann_weight;
	dest->free_energy = src->free_energy;
	dest->status = src->status;

	for(int j=0; j<num_genes; j++){
	        dest->genes[j].to_ic = src->genes[j].to_ic;
		dest->genes[j].to_int32 = src->genes[j].to_int32;
	}
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void save_snapshot(chromosome* chrom_snapshot, const chromosome* chrom, int num_chrom, int num_genes){

	for(int i=0; i<num_chrom; i++)
		copy_chrom(&chrom_snapshot[i],&chrom[i],num_genes);

}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int check_state(char* pausefile, char* abortfile, char* stopfile, int interval){
	FILE* STATE;

	STATE = NULL;

	// try and open pause/stop file
	// (works with PyMOL interface)

	STATE = fopen(pausefile,"r");
	if(STATE != NULL) {
		do {
			fclose(STATE);

# ifdef _WIN32
			Sleep(SLEEP);
# else
			usleep(SLEEP*1000);
# endif
			STATE = fopen(pausefile,"r");

		}while(STATE != NULL);
	}

	STATE = fopen(abortfile,"r");
	if(STATE != NULL) {
		fclose(STATE);
		printf("manual aborting\n");
		return -1;
	}

	STATE = fopen(stopfile,"r");
	if(STATE != NULL) {
		fclose(STATE);
		printf("simulation stopped prematurely\n");
		return 1;
	}

	return 0;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void fitness_stats(GB_Global* GB, const chromosome* chrom,int pop_size){
	int i;
	int flag;

	//calculate fitness max and and average of the whole pop
	GB->fit_max=0.0;
	GB->fit_avg=0.0;

	flag=1;
	for(i=0;i<pop_size;i++){
		if (flag){
			GB->fit_max=chrom[i].fitnes;
			flag=0;
		}

		if(chrom[i].fitnes > GB->fit_max)
			GB->fit_max=chrom[i].fitnes;

		GB->fit_avg+=chrom[i].fitnes;
	}

	GB->fit_avg /= (double)pop_size;

	return;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void adapt_prob(GB_Global* GB,double fit1, double fit2, double* mutp, double* crossp){
	//printf("crossing fit1[%8.3f] with fit2[%8.3f]\n",fit1,fit2);

	//find which crossed individual has higher fitness
	if(fit1 > fit2){
		GB->fit_high=fit1;
		GB->fit_low=fit2;
	}else{
		GB->fit_high=fit2;
		GB->fit_low=fit1;
	}

	//crossp=k1 when high=avg
	//mutp=k2 when high=avg
	//crossp/mutp=0 when high=max

	//calculate new probabilities (pc/pm)
	double denom = GB->fit_max - GB->fit_avg;
	if (denom < GA_FITNESS_DENOM_FLOOR) denom = GA_FITNESS_DENOM_FLOOR;  // prevent division by zero when converged

	if (GB->fit_high > GB->fit_avg) *crossp = GB->k1*(GB->fit_max-GB->fit_high)/denom;
	else *crossp = GB->k3;

	if (GB->fit_low > GB->fit_avg) *mutp = GB->k2*(GB->fit_max-GB->fit_low)/denom;
	else *mutp = GB->k4;

	/*
	  printf("f'=%.1f\tf=%.1f\tfmax=%.1f\tfavg=%.1f\t\tPc=%5.3f\tPm=%5.3f\n",
	  GB->fit_high,GB->fit_low,GB->fit_max,GB->fit_avg,*crossp,*mutp);
	*/

	return;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int reproduce(FA_Global* FA,GB_Global* GB,VC_Global* VC, chromosome* chrom, const genlim* gene_lim,
               atom* atoms,resid* residue,gridpoint* cleftgrid,char* repmodel,
               double mutprob, double crossprob, int print,
	       std::function<int32_t()> & dice,
	       std::unordered_map<size_t, int> & duplicates,
               cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
               GAContext& ctx){

	int& nrejected = ctx.nrejected;

	int i,j,k;
	int nnew,p1,p2;

	gene chrop1_gen[MAX_NUM_GENES];
	gene chrop2_gen[MAX_NUM_GENES];

	int num_genes_wo_sc=0;

	/*
	  std::mt19937 rng;
	  std::uniform_int_distribution<int32_t> one_to_max_int32( 1, MAX_RANDOM_VALUE );
	  std::function<int32_t()> dice = [&](){ return one_to_max_int32(rng); };
	*/

	if(strcmp(repmodel,"STEADY")==0){
		nnew = GB->ssnum;
	}else if(strcmp(repmodel,"BOOM")==0){
		nnew = (int)(GB->pbfrac*(double)GB->num_chrom);
	}else{
		nnew = 0;
	}

	i=0;
	while(i<nnew){

		/************************************/
		/****** SELECTION OF PARENTS ********/
		/************************************/
		p1=roullete_wheel(chrom,GB->num_chrom);
		p2=roullete_wheel(chrom,GB->num_chrom);
		if (GB->adaptive_ga) adapt_prob(GB,chrom[p1].fitnes,chrom[p2].fitnes,&mutprob,&crossprob);

		/************************************/
		/****** CROSSOVER OPERATOR  ********/
		/************************************/
		// create temporary genes
		memcpy(chrop1_gen,chrom[p1].genes,GB->num_genes*sizeof(gene));
		memcpy(chrop2_gen,chrom[p2].genes,GB->num_genes*sizeof(gene));

		if(RandomDouble() < crossprob){
			crossover(chrop1_gen,chrop2_gen,GB->num_genes,GB->intragenes);
		}

		/************************************/
		/****** MUTATION OPERATOR  ********/
		/************************************/
		num_genes_wo_sc = GB->num_genes-FA->nflxsc_real;

		mutate(chrop1_gen,GB->num_genes-FA->nflxsc_real,mutprob);
		k=0;
		for(j=0;j<FA->nflxsc;j++){
			if(residue[FA->flex_res[j].inum].trot != 0){
				if(RandomDouble() < FA->flex_res[j].prob){
					mutate(&chrop1_gen[num_genes_wo_sc+k],1,mutprob);
				}
				k++;
			}
		}

		mutate(chrop2_gen,GB->num_genes-FA->nflxsc_real,mutprob);
		k=0;
		for(j=0;j<FA->nflxsc;j++){
			if(residue[FA->flex_res[j].inum].trot != 0){
				if(RandomDouble() < FA->flex_res[j].prob){
					mutate(&chrop2_gen[num_genes_wo_sc+k],1,mutprob);
				}
				k++;
			}
		}

		for(j=0; j<GB->num_genes; j++){
			chrop1_gen[j].to_ic = genetoic(&gene_lim[j],chrop1_gen[j].to_int32);
			chrop2_gen[j].to_ic = genetoic(&gene_lim[j],chrop2_gen[j].to_int32);
		}

		size_t sig1 = hash_genes(chrop1_gen,GB->num_genes);
		size_t sig2 = hash_genes(chrop2_gen,GB->num_genes);

		/************************************/
		/****** CHECK DUPLICATION  ********/
		/************************************/
		if(GB->duplicates || duplicates.find(sig1) == duplicates.end()){

			/*
			  if(!FA->useflexdee ||
			  cmp_chrom2rotlist(FA->psFlexDEENode,chrom,gene_lim,num_genes_wo_sc,
			  FA->nflxsc_real,GB->num_chrom,FA->FlexDEE_Nodes)==0){
			*/

			//nrejected += filter_deelig(FA,GB,chrom,chrop1_gen,GB->num_chrom+i,atoms,gene_lim,dice);
			memcpy(chrom[GB->num_chrom+i].genes,chrop1_gen,GB->num_genes*sizeof(gene));

			chrom[GB->num_chrom+i].cf=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,
								  chrom[GB->num_chrom+i].genes,target);
			chrom[GB->num_chrom+i].evalue=get_cf_evalue(&chrom[GB->num_chrom+i].cf);
			chrom[GB->num_chrom+i].app_evalue=get_apparent_cf_evalue(&chrom[GB->num_chrom+i].cf);
			ccbm_inject_strain(FA, chrom[GB->num_chrom+i], gene_lim);  // CCBM strain
			chrom[GB->num_chrom+i].status='n';

			duplicates[sig1] = 1;
			i++;
		}

		if(i==nnew) break;

		if(GB->duplicates || duplicates.find(sig2) == duplicates.end()){

			/*
			  if(!FA->useflexdee ||
			  cmp_chrom2rotlist(FA->psFlexDEENode,chrom,gene_lim,num_genes_wo_sc,
			  FA->nflxsc_real,GB->num_chrom,FA->FlexDEE_Nodes)==0){
			*/
			//nrejected += filter_deelig(FA,GB,chrom,chrop2_gen,GB->num_chrom+i,atoms,gene_lim,dice);
			memcpy(chrom[GB->num_chrom+i].genes,chrop2_gen,GB->num_genes*sizeof(gene));

			chrom[GB->num_chrom+i].cf=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,
								  chrom[GB->num_chrom+i].genes,target);
			chrom[GB->num_chrom+i].evalue=get_cf_evalue(&chrom[GB->num_chrom+i].cf);
			chrom[GB->num_chrom+i].app_evalue=get_apparent_cf_evalue(&chrom[GB->num_chrom+i].cf);
			ccbm_inject_strain(FA, chrom[GB->num_chrom+i], gene_lim);  // CCBM strain
			chrom[GB->num_chrom+i].status='n';

			duplicates[sig2] = 1;
			i++;
		}
	}

	if(strcmp(repmodel,"STEADY")==0){
		// replace the n individuals from the old population with the new one (elitism)
		QuickSort(chrom,0,GB->num_chrom-1,true);
		for(i=0;i<nnew;i++) chrom[GB->num_chrom-1-i]=chrom[GB->num_chrom+i];
		calculate_fitness(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,
				  GB->fitness_model,GB->num_chrom,print,target,ctx);
	}else if(strcmp(repmodel,"BOOM")==0){
		// merge and sort both merged populations
		calculate_fitness(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,
				  GB->fitness_model,GB->num_chrom+nnew,print,target,ctx);
	}

	//printf("number of conformers rejected: %d\n", nrejected);

	return nrejected;
}

size_t hash_genes(const gene* g, int n){
	size_t h = 0;
	for(int i = 0; i < n; ++i)
		h ^= std::hash<int32_t>{}(static_cast<int32_t>(g[i].to_ic + 0.5)) + 0x9e3779b9 + (h << 6) + (h >> 2);
	return h;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int filter_deelig(FA_Global* FA, GB_Global* GB, chromosome* chrom, gene* genes, int ci, atom* atoms, const genlim* gene_lim,
		   std::function<int32_t()> & dice)
{
	int nrejected = 0;

	if(FA->deelig_flex && FA->nflexbonds){

		int j,deelig_list[GA_MAX_DEELIG_DIHEDRALS];

		for(j=1; j<=FA->resligand->fdih; j++)
			deelig_list[j] = GA_DEELIG_SENTINEL;

		for(j=0; j<GB->num_genes; j++)
			if(FA->map_par[j].typ == 2 && FA->map_par[j].bnd != -1)
				deelig_list[FA->map_par[j].bnd] = (int)(genes[j].to_ic+0.5);

		/*
		printf("searched deelig list = [");
		for(int k=1; k<=FA->resligand->fdih; k++){
			printf("%d,", deelig_list[k]);
		}
		printf("]\n");
		*/

		if(deelig_search(FA->deelig_root_node, deelig_list, FA->resligand->fdih)){
			/*
			printf("conformer rejected:");
			for(j=1; j<=FA->resligand->fdih; j++)
				printf("%d ", deelig_list[j]);
			printf("\n");
			getchar();
			*/

			// generate a new conformer until the conformer
			// has not already been assigned as 'clashing conformer'
			// and is also not a duplicate
			do{
				nrejected++;

				// only generate a new conformer, do not modify other variables
				generate_random_individual(FA,GB,atoms,genes,gene_lim,dice,
							   FA->map_par_flexbond_first_index,
							   FA->map_par_flexbond_first_index+FA->nflexbonds);

				for(j=1; j<=FA->resligand->fdih; j++)
					deelig_list[j] = GA_DEELIG_SENTINEL;

				for(j=0; j<GB->num_genes; j++)
					if(FA->map_par[j].typ == 2 && FA->map_par[j].bnd != -1)
						deelig_list[FA->map_par[j].bnd] = (int)(genes[j].to_ic+0.5);
				/*
				printf("searched do-while deelig list = [");
				for(int k=1; k<=FA->resligand->fdih; k++){
					printf("%d,", deelig_list[k]);
				}
				printf("]\n");
				*/

				/*
				if(deelig_search(FA->deelig_root_node, deelig_list, FA->resligand->fdih)){
					printf("do-while conformer rejected:");
					for(j=1; j<=FA->resligand->fdih; j++)
						printf("%d ", deelig_list[j]);
					printf("\n");
					getchar();
				}
				*/
			}while((!GB->duplicates && cmp_chrom2pop(chrom,genes,GB->num_genes,0,ci)) ||
			       (deelig_search(FA->deelig_root_node, deelig_list, FA->resligand->fdih)));
		}

	}

	return nrejected;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int deelig_search(struct deelig_node_struct* root_node, int* deelig_list, int fdih)
{
	std::map<int, struct deelig_node_struct*>::iterator it;
	struct deelig_node_struct* node = root_node;

	for(int i=1; i<=fdih; i++){
		//printf("[%d]: searching %d\n", i, deelig_list[i]);
		if((it=node->childs.find(deelig_list[i])) != node->childs.end() ||
		   (deelig_list[i] != GA_DEELIG_SENTINEL && (it=node->childs.find(GA_DEELIG_SENTINEL)) != node->childs.end())){
			//printf("found %d\n", it->first);
			node = it->second;
		}else{
			//printf("not found %d\n", deelig_list[i]);
			return(0);
		}
	}

	return(1);
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int roullete_wheel(const chromosome* chrom,int n){
	double r;
	double tot=0.0;
	int i;

	if (n <= 0) return 0;

	for(i=0;i<n;i++){tot += chrom[i].fitnes;}

	// Guard: if total fitness is zero or negative, return random index
	if (tot <= 0.0) return static_cast<int>(RandomDouble() * n) % n;

	r=RandomDouble()*tot;

	i=0;
	tot=0.0;
	while(tot <= r && i < n){
		tot += chrom[i].fitnes;
		i++;
	}
	i--;
	if (i < 0) i = 0;

	return i;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void calculate_fitness(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,
                       atom* atoms,resid* residue,gridpoint* cleftgrid,char method[], int pop_size, int print,
                       cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
                       GAContext& ctx){

	int& gen_id = ctx.gen_id;
	int i;

	// ── TurboQuant QuantizedContactMatrix (TQCM) ─────────────────────────
	// When TQCM is enabled, build a compressed representation of the
	// energy interaction matrix FA->energy_matrix at first call.  The QCM
	// stores 256×256 rows quantized at 2 bits/coordinate, giving 16×
	// compression.  Subsequent calls can use
	// qcm.approximate_score(type_i, type_j) for fast scoring.
	// The compressed matrix is rebuilt whenever ntypes changes.
	auto& s_tqcm = ctx.tqcm;
	auto& s_tqcm_ntypes = ctx.tqcm_ntypes;

	if (FA->use_tqcm && (gen_id == 0 || s_tqcm_ntypes != FA->ntypes)) {
		// Sample the energy_matrix spline functions at midpoint (area=0.5)
		// to build a flat ntypes×ntypes matrix suitable for QCM compression.
		// QuantizedContactMatrix expects exactly 256×256 floats, so we pad
		// to 256 if ntypes < 256 (values beyond ntypes are zero-filled).
		const int nt = FA->ntypes;
		constexpr int QCM_DIM = turboquant::QuantizedContactMatrix::kNumAtomTypes;
		std::vector<float> flat_matrix(QCM_DIM * QCM_DIM, 0.0f);

		for (int t1 = 0; t1 < nt && t1 < QCM_DIM; ++t1) {
			for (int t2 = 0; t2 < nt && t2 < QCM_DIM; ++t2) {
				struct energy_matrix* em = &FA->energy_matrix[t1 * nt + t2];
				if (em->energy_values != NULL) {
					// Sample the density-of-contact curve at area = 0.5
					// This gives the representative interaction strength
					flat_matrix[t1 * QCM_DIM + t2] = static_cast<float>(get_yval(em, GA_TQCM_SAMPLE_AREA));
				}
			}
		}

		delete s_tqcm;
		s_tqcm = new turboquant::QuantizedContactMatrix(/*bit_width=*/GA_TQCM_BIT_WIDTH);
		s_tqcm->build(flat_matrix.data());
		s_tqcm_ntypes = nt;

		printf("--- TurboQuant QuantizedContactMatrix (TQCM) built ---\n");
		printf("  Atom types             = %d (padded to %d)\n", nt, QCM_DIM);
		printf("  Bit width              = %d\n", s_tqcm->bit_width());
		printf("  Compression ratio      = %.1f×\n", s_tqcm->compression_ratio());
		printf("  Memory (compressed)    = %zu bytes\n", s_tqcm->memory_bytes());
		printf("  Memory (original)      = %zu bytes\n",
		       static_cast<size_t>(QCM_DIM * QCM_DIM) * sizeof(float));
		printf("  MSE bound per element  = %.6f\n",
		       s_tqcm->quantizer().theoretical_mse());

		// Validate: spot-check a few type pairs against exact values
		if (nt >= 2) {
			double max_err = 0.0;
			int n_checks = std::min(nt * nt, GA_TQCM_MAX_SPOT_CHECKS);
			for (int c = 0; c < n_checks; ++c) {
				int ti = c / nt, tj = c % nt;
				float exact_val = flat_matrix[ti * QCM_DIM + tj];
				float approx_val = s_tqcm->approximate_score(ti, tj);
				double ae = std::abs(static_cast<double>(exact_val - approx_val));
				max_err = std::max(max_err, ae);
			}
			printf("  Spot-check max |error| = %.6f (over %d pairs)\n", max_err, n_checks);
		}
	}


	// ── Chromosome evaluation ────────────────────────────────────────────────
	// Runtime dispatch: CUDA GPU → Metal GPU → OpenMP CPU (thread-safe).
	// All compiled-in backends are available simultaneously; select_backend()
	// picks the best one at runtime based on detected hardware.

#if defined(FLEXAIDS_USE_CUDA) || defined(FLEXAIDS_USE_METAL)
	// Helper lambda: sample each energy-matrix density function at n_samples
	// evenly-spaced x values in [0, 1] and pack into a flat float array
	// [n_types × n_types × n_samples] for GPU upload.
	// When Eigen is available, the x-value linspace is built via Eigen::ArrayXd
	// for vectorised construction; the get_yval evaluation loop is then
	// auto-vectorisable because it operates on a contiguous double buffer.
	auto build_emat_sampled = [&](int n_types, int n_samples) -> std::vector<float> {
		const size_t total = static_cast<size_t>(n_types) * n_types * n_samples;
		std::vector<float> out(total, 0.0f);

		// Build the x-sample linspace via Eigen (vectorised).
		Eigen::ArrayXd xs = Eigen::ArrayXd::LinSpaced(n_samples, 0.0, 1.0);
		for (int t1 = 0; t1 < n_types; ++t1) {
			for (int t2 = 0; t2 < n_types; ++t2) {
				struct energy_matrix* em = &FA->energy_matrix[t1 * n_types + t2];
				if (em->energy_values == NULL) continue;
				float* dst = &out[(t1 * n_types + t2) * n_samples];
				for (int k = 0; k < n_samples; ++k)
					dst[k] = static_cast<float>(get_yval(em, xs[k]));
			}
		}
		return out;
	};

	// Helper lambda: pack gene internal coordinates into a flat array for GPU.
	auto pack_genes_batch = [&](int n_genes) -> std::vector<double> {
		std::vector<double> h_genes(pop_size * n_genes, 0.0);
		for (int c = 0; c < pop_size; ++c)
			for (int g = 0; g < n_genes; ++g)
				h_genes[c * n_genes + g] = chrom[c].genes[g].to_ic;
		return h_genes;
	};

	// Helper lambda: unpack GPU batch results into chromosome CF structures.
	auto unpack_gpu_results = [&](const std::vector<double>& h_com,
	                              const std::vector<double>& h_wal,
	                              const std::vector<double>& h_sas) {
		for (int c = 0; c < pop_size; ++c) {
			if (chrom[c].status != 'n') {
				chrom[c].cf.com    = h_com[c];
				chrom[c].cf.wal    = h_wal[c];
				chrom[c].cf.sas    = h_sas[c];
				chrom[c].cf.con    = 0.0;
				chrom[c].cf.gist   = 0.0;
				chrom[c].cf.hbond  = 0.0;
				chrom[c].cf.totsas = 0.0;
				chrom[c].cf.rclash = (h_wal[c] > CLASH_THRESHOLD) ? 1 : 0;
				chrom[c].evalue     = get_cf_evalue(&chrom[c].cf);
				chrom[c].app_evalue = get_apparent_cf_evalue(&chrom[c].cf);
				ccbm_inject_strain(FA, chrom[c], gene_lim);  // CCBM strain
				chrom[c].status    = 'n';
			}
		}
	};

	// Helper lambda: prepare GPU atom data arrays from the atoms array.
	struct GPUAtomData {
		std::vector<float> xyz;
		std::vector<int>   type;
		std::vector<float> radius;
		int lig_first;
		int lig_last;
	};
	auto prepare_gpu_atoms = [&]() -> GPUAtomData {
		const int n_atoms = FA->atm_cnt_real;
		GPUAtomData d;
		d.xyz.resize(n_atoms * 3);
		d.type.resize(n_atoms);
		d.radius.resize(n_atoms);
		for (int a = 0; a < n_atoms; ++a) {
			d.xyz[a*3+0] = atoms[a].coor[0];
			d.xyz[a*3+1] = atoms[a].coor[1];
			d.xyz[a*3+2] = atoms[a].coor[2];
			d.type[a]    = atoms[a].type - 1;  // 1-based → 0-based
			d.radius[a]  = atoms[a].radius;
		}
		d.lig_first = (FA->resligand && FA->resligand->fatm)
		            ? FA->resligand->fatm[0] : 0;
		d.lig_last  = (FA->resligand && FA->resligand->latm)
		            ? FA->resligand->latm[0] : 0;
		return d;
	};
#endif  // FLEXAIDS_USE_CUDA || FLEXAIDS_USE_METAL

	// Log dispatch decision on first call.
	[[maybe_unused]] const auto backend = flexaids::select_backend();
	if (!ctx.dispatch_logged) {
		auto report = flexaids::get_dispatch_report();
		fprintf(stderr, "[FlexAIDdS] Hardware dispatch: %s (%s)\n",
		        flexaids::backend_name(static_cast<flexaids::HardwareBackend>(
		            static_cast<uint8_t>(report.selected))), report.reason.c_str());
		ctx.dispatch_logged = true;
	}

	[[maybe_unused]] bool gpu_handled = false;

#ifdef FLEXAIDS_USE_CUDA
	if (backend == flexaids::HardwareBackend::CUDA) {
		const int n_atoms = FA->atm_cnt_real;
		const int n_types = FA->ntypes;
		const int n_genes = GB->num_genes;
		const int ns      = CUDA_EMAT_SAMPLES;

		// Thread-safe GPU context pool — shared across concurrent GA instances
		auto& pool = GPUContextPool::instance();
		auto handle = pool.acquire_cuda(n_atoms, n_types, [&]() {
			auto ad = prepare_gpu_atoms();
			std::vector<float> h_emat = build_emat_sampled(n_types, ns);
			return cuda_eval_init(n_atoms, n_types, MAX_NUM_CHROM,
			                     n_genes, ad.lig_first, ad.lig_last,
			                     FA->permeability,
			                     ad.xyz.data(), ad.type.data(),
			                     ad.radius.data(), h_emat.data());
		});

		std::vector<double> h_genes = pack_genes_batch(n_genes);
		std::vector<double> h_com(pop_size), h_wal(pop_size), h_sas(pop_size);
		cuda_eval_batch(handle.ctx, pop_size, n_genes, h_genes.data(),
		                h_com.data(), h_wal.data(), h_sas.data());
		unpack_gpu_results(h_com, h_wal, h_sas);
		pool.release_cuda(handle);
		gpu_handled = true;
	}
#endif

#ifdef FLEXAIDS_USE_METAL
	if (!gpu_handled && backend == flexaids::HardwareBackend::METAL) {
		const int n_atoms = FA->atm_cnt_real;
		const int n_types = FA->ntypes;
		const int n_genes = GB->num_genes;
		const int ns      = METAL_EMAT_SAMPLES;

		auto& pool = GPUContextPool::instance();
		auto handle = pool.acquire_metal(n_atoms, n_types, [&]() {
			auto ad = prepare_gpu_atoms();
			std::vector<float> h_emat = build_emat_sampled(n_types, ns);
			return metal_eval_init(n_atoms, n_types, MAX_NUM_CHROM,
			                      ad.lig_first, ad.lig_last,
			                      FA->permeability,
			                      ad.xyz.data(), ad.type.data(),
			                      ad.radius.data(), h_emat.data(), ns);
		});

		std::vector<double> h_genes = pack_genes_batch(n_genes);
		std::vector<double> h_com(pop_size), h_wal(pop_size), h_sas(pop_size);
		metal_eval_batch(handle.ctx, pop_size, n_genes, h_genes.data(),
		                 h_com.data(), h_wal.data(), h_sas.data());
		unpack_gpu_results(h_com, h_wal, h_sas);
		pool.release_metal(handle);
		gpu_handled = true;
	}
#endif

	if (!gpu_handled) {
	// ── Thread-safe CPU path (AVX-512 / AVX2 / OpenMP / scalar) ─────────
	// Each OpenMP thread receives its own private copies of every data
	// structure that Vcontacts/vcfunction/ic2cf writes to:
	//   • atoms[]        – internal coords (dis/ang/dih) and Cartesian (coor)
	//   • residue[]      – rotamer index (.rot)
	//   • FA scratch     – contacts[], contributions[], optres[].cf
	//   • VC workspace   – Calc[], Calclist[], ca_index[], ca_rec[],
	//                      seed[], contlist[], ptorder[], centerpt[],
	//                      poly[], cont[], vedge[]
	// Read-only fields (energy_matrix, map_par, …) are shared.
	// The DEE linked-list update in ic2cf is skipped in parallel mode
	// (guarded by omp_in_parallel() in ic2cf.cpp) to avoid concurrent
	// linked-list corruption; DEE pruning still operates in serial calls.
	{
#ifdef _OPENMP
		const int n_thr = omp_get_max_threads();
#else
		const int n_thr = 1;
#endif
		const int natm  = FA->atm_cnt;
		const int natmr = FA->atm_cnt_real;
		const int nres  = FA->res_cnt;
		const int nopt  = FA->num_optres;
		const int nctb  = FA->ntypes * FA->ntypes;

		// ── Dirty-tracking optimisation ─────────────────────────────────
		// ic2cf only modifies atoms belonging to optimizable residues
		// (ligand + flex sidechains) and buildcc rebuilds their Cartesian
		// coords. vcfunction writes .acs for these same atoms.
		// When normal modes are disabled, we restore only these "dirty"
		// atoms per chromosome instead of copying the entire atom array.
		// This reduces per-eval memory bandwidth by 90%+ for typical systems.
		bool has_normal_modes = false;
		for (int p = 0; p < FA->npar; ++p) {
			if (FA->map_par[p].typ == 3) { has_normal_modes = true; break; }
		}

		// Build sorted unique list of atom indices modified by ic2cf.
		// Sources: mov[] lists (buildcc targets) + map_par[].atm (IC targets).
		std::vector<int> dirty_atm;
		std::vector<int> dirty_res_idx;
		if (!has_normal_modes) {
			// Atoms in mov[] rebuild lists (ligand + flex sidechain Cartesian)
			for (int r = 0; r < FA->nors; ++r)
				for (int m = 0; m < FA->nmov[r]; ++m)
					dirty_atm.push_back(FA->mov[r][m]);
			// Atoms directly referenced by map_par (IC fields: dis/ang/dih)
			for (int p = 0; p < FA->npar; ++p)
				dirty_atm.push_back(FA->map_par[p].atm);
			// Cascade dihedral atoms (atoms whose .dih depends on a flex bond)
			for (int p = 0; p < FA->npar; ++p) {
				if (FA->map_par[p].typ == 2) {
					int j = FA->map_par[p].atm;
					int cat = atoms[j].rec[3];
					while (cat != 0 && cat != FA->map_par[p].atm) {
						dirty_atm.push_back(cat);
						j = cat;
						cat = atoms[j].rec[3];
					}
				}
			}
			// Sort and deduplicate
			std::sort(dirty_atm.begin(), dirty_atm.end());
			dirty_atm.erase(std::unique(dirty_atm.begin(), dirty_atm.end()),
			                dirty_atm.end());

			// Residue indices with rotamer genes (typ==4 modifies .rot)
			for (int p = 0; p < FA->npar; ++p) {
				if (FA->map_par[p].typ == 4)
					dirty_res_idx.push_back(atoms[FA->map_par[p].atm].ofres);
			}
			std::sort(dirty_res_idx.begin(), dirty_res_idx.end());
			dirty_res_idx.erase(
				std::unique(dirty_res_idx.begin(), dirty_res_idx.end()),
				dirty_res_idx.end());
		}
		const bool use_selective = !has_normal_modes &&
		    static_cast<int>(dirty_atm.size()) < natm / 2;
		const int n_dirty_atm = static_cast<int>(dirty_atm.size());
		const int n_dirty_res = static_cast<int>(dirty_res_idx.size());

		// Per-thread mutable atom arrays.
		std::vector<std::vector<atom>>  tl_atoms(n_thr,
		    std::vector<atom>(atoms, atoms + natm + 1));
		// Per-thread residue arrays (pointer fields shared read-only; .rot private).
		std::vector<std::vector<resid>> tl_res(n_thr,
		    std::vector<resid>(residue, residue + nres + 1));
		// Per-thread FA copies with redirected mutable scratch buffers.
		std::vector<FA_Global>           tl_fa(n_thr, *FA);
		std::vector<std::vector<int>>    tl_contacts(n_thr, std::vector<int>(MAX_ATOM_NUMBER, 0));
		std::vector<std::vector<float>>  tl_contrib(n_thr, std::vector<float>(nctb, 0.0f));
		std::vector<std::vector<OptRes>> tl_optres(n_thr,
		    std::vector<OptRes>(FA->optres, FA->optres + nopt));
		// Per-thread VC workspace (Vcontacts writes all these each call).
		std::vector<VC_Global>               tl_vc(n_thr, *VC);
		std::vector<std::vector<atomsas>>    tl_calc(n_thr, std::vector<atomsas>(natmr));
		std::vector<std::vector<int>>        tl_calclist(n_thr, std::vector<int>(natmr));
		std::vector<std::vector<int>>        tl_caidx(n_thr, std::vector<int>(natmr, -1));
		std::vector<std::vector<ca_struct>>  tl_carec(n_thr,
		    std::vector<ca_struct>(VC->ca_recsize));
		std::vector<std::vector<int>>        tl_seed(n_thr,
		    std::vector<int>(3 * natmr));
		std::vector<std::vector<contactlist>> tl_contlist(n_thr,
		    std::vector<contactlist>(GA_CONTLIST_SIZE));
		std::vector<std::vector<ptindex>>    tl_ptorder(n_thr,
		    std::vector<ptindex>(MAX_PT));
		std::vector<std::vector<vertex>>     tl_centerpt(n_thr,
		    std::vector<vertex>(MAX_PT));
		std::vector<std::vector<vertex>>     tl_poly(n_thr,
		    std::vector<vertex>(MAX_POLY));
		std::vector<std::vector<plane>>      tl_cont(n_thr,
		    std::vector<plane>(MAX_PT));
		std::vector<std::vector<edgevector>> tl_vedge(n_thr,
		    std::vector<edgevector>(MAX_POLY));

		for (int t = 0; t < n_thr; ++t) {
			// Redirect FA mutable scratch to per-thread buffers.
			tl_fa[t].contacts      = tl_contacts[t].data();
			tl_fa[t].contributions = tl_contrib[t].data();
			tl_fa[t].optres        = tl_optres[t].data();
			// Redirect VC mutable workspace to per-thread buffers.
			tl_vc[t].Calc      = tl_calc[t].data();
			tl_vc[t].Calclist  = tl_calclist[t].data();
			tl_vc[t].ca_index  = tl_caidx[t].data();
			tl_vc[t].ca_rec    = tl_carec[t].data();
			tl_vc[t].seed      = tl_seed[t].data();
			tl_vc[t].contlist  = tl_contlist[t].data();
			tl_vc[t].ptorder   = tl_ptorder[t].data();
			tl_vc[t].centerpt  = tl_centerpt[t].data();
			tl_vc[t].poly      = tl_poly[t].data();
			tl_vc[t].cont      = tl_cont[t].data();
			tl_vc[t].vedge     = tl_vedge[t].data();
			// box is shared: if vindex==1 it's pre-built read-only;
			// if vindex==0 Vcontacts will malloc/vcfunction will free per call.
		}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(chrom, pop_size, GB, gene_lim, cleftgrid, target, \
	       atoms, residue, FA, VC, \
	       tl_atoms, tl_res, tl_fa, tl_optres, tl_vc, \
	       natm, nres, nopt, \
	       use_selective, dirty_atm, dirty_res_idx, n_dirty_atm, n_dirty_res)
#endif
		for (int ii = 0; ii < pop_size; ++ii) {
			if (chrom[ii].status == 'n') continue;
#ifdef _OPENMP
			const int tid = omp_get_thread_num();
#else
			const int tid = 0;
#endif
			// Reset per-thread state to the reference protein configuration.
			// When normal modes are off, only restore the atoms/residues that
			// ic2cf + vcfunction actually modify (typically <10% of total).
			if (use_selective) {
				for (int d = 0; d < n_dirty_atm; ++d) {
					const int ai = dirty_atm[d];
					tl_atoms[tid][ai] = atoms[ai];
				}
				for (int d = 0; d < n_dirty_res; ++d) {
					const int ri = dirty_res_idx[d];
					tl_res[tid][ri] = residue[ri];
				}
			} else {
				std::copy(atoms,   atoms + natm + 1,   tl_atoms[tid].begin());
				std::copy(residue, residue + nres + 1, tl_res[tid].begin());
			}
			// Redirect per-thread atom optres pointers to per-thread optres array.
			// atoms[j].optres points to FA->optres (original); redirect to tl_optres[tid]
			// so vcfunction scoring writes to (and ic2cf reads from) the same buffer.
			for (int ai = 1; ai <= natm; ++ai) {
				atom& a = tl_atoms[tid][ai];
				if (a.optres) {
					ptrdiff_t oidx = a.optres - FA->optres;
					a.optres = &tl_optres[tid][oidx];
				}
			}
			// optres cf fields are cleared by vcfunction itself; pre-clear for safety.
			for (int o = 0; o < nopt; ++o) {
				tl_optres[tid][o].cf.com    = 0.0;
				tl_optres[tid][o].cf.wal    = 0.0;
				tl_optres[tid][o].cf.sas    = 0.0;
				tl_optres[tid][o].cf.totsas = 0.0;
				tl_optres[tid][o].cf.con    = 0.0;
				tl_optres[tid][o].cf.gist   = 0.0;
				tl_optres[tid][o].cf.elec   = 0.0;
				tl_optres[tid][o].cf.hbond  = 0.0;
				tl_optres[tid][o].cf.gist_desolv = 0.0;
				tl_optres[tid][o].cf.rclash = 0;
			}
			tl_vc[tid].numcarec = 0;

			chrom[ii].cf = eval_chromosome(
			    &tl_fa[tid], GB, &tl_vc[tid], gene_lim,
			    tl_atoms[tid].data(), tl_res[tid].data(),
			    cleftgrid, chrom[ii].genes, target);
			chrom[ii].evalue     = get_cf_evalue(&chrom[ii].cf);
			chrom[ii].app_evalue = get_apparent_cf_evalue(&chrom[ii].cf);
			ccbm_inject_strain(FA, chrom[ii], gene_lim);  // CCBM strain
			chrom[ii].status     = 'n';
		}
	}
	}  // !gpu_handled

	QuickSort(chrom,0,pop_size-1,true);

	//print_par(chrom,gene_lim,5,GB->num_genes);
	//PAUSE;
	//chrom_hpsort(pop_size,0,chrom);

	if(strcmp(method,"LINEAR")==0){
		/* the fitness value is a number between 0 and num_chrom.
		   each chromosome is assigned an integer value that
		   corresponds to its position in index_map.
		*/
		for(i=0;i<GB->num_chrom;i++){
			chrom[i].fitnes=(double)(GB->num_chrom-i);
		}
	}

	if(strcmp(method,"PSHARE")==0){
		/* the fitness value is a number between 0 and num_chrom.
		   each chromosome is assigned an integer value that
		   corresponds to its position in index_map. Moreover,
		   each chromosome's fitness is lowered by sharing.
		   The niche count (share) must be accumulated over ALL j before
		   dividing — fixed from the previous per-j assignment bug.
		   The outer loop is data-race free (each i writes only chrom[i].fitnes)
		   and is parallelised with OpenMP.
		*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(chrom, GB, FA, cleftgrid)
#endif
		for(int pi=0; pi<GB->num_chrom; pi++){
			double pshare = 0.0;
			for(int pj=0; pj<GB->num_chrom; pj++){
				double prmsp = calc_rmsp(GB->num_genes,
				                         chrom[pi].genes, chrom[pj].genes,
				                         FA->map_par, cleftgrid);
				if(prmsp <= GB->sig_share){
					pshare += (1.0 - pow((prmsp/GB->sig_share), GB->alpha));
				}
			}
			// Assign fitness AFTER accumulating the full niche count.
			chrom[pi].fitnes = (double)(GB->num_chrom - pi) / pshare;
		}
	}

	if(strcmp(method,"SMFREE")==0){
		/* SMFREE — StatMech Free-energy-weighted fitness with niche sharing.
		   Uses the StatMechEngine to compute Boltzmann weights from the
		   current population's energies. Fitness blends rank-based fitness
		   with thermodynamic Boltzmann probability:
		     fitness_i = [(1-w) * rank_component + w * boltzmann_component] / share_i
		   where w = entropy_weight ∈ [0,1].
		   This biases selection toward thermodynamically favorable poses
		   (low free energy) while maintaining diversity via niche sharing.
		*/
		if (FA->temperature > 0) {
			const double T = static_cast<double>(FA->temperature);
			statmech::StatMechEngine engine(T);

			// Feed all chromosome energies into the engine.
			for (int si = 0; si < GB->num_chrom; si++) {
				engine.add_sample(chrom[si].evalue);
			}

			// Compute ensemble thermodynamics and Boltzmann weights.
			auto thermo = engine.compute();
			auto bweights = engine.boltzmann_weights();

			// Store Boltzmann weights and free energy on each chromosome.
			for (int si = 0; si < GB->num_chrom; si++) {
				chrom[si].boltzmann_weight = bweights[static_cast<size_t>(si)];
				chrom[si].free_energy = thermo.free_energy;
			}

			// Find max Boltzmann weight for normalisation of the Boltzmann component.
			double max_bw = 0.0;
			for (int si = 0; si < GB->num_chrom; si++) {
				if (chrom[si].boltzmann_weight > max_bw)
					max_bw = chrom[si].boltzmann_weight;
			}
			if (max_bw <= 0.0) max_bw = 1.0;

			const double w = GB->entropy_weight;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(chrom, GB, FA, cleftgrid, max_bw, w)
#endif
			for (int pi = 0; pi < GB->num_chrom; pi++) {
				// Niche sharing (same as PSHARE).
				double pshare = 0.0;
				for (int pj = 0; pj < GB->num_chrom; pj++) {
					double prmsp = calc_rmsp(GB->num_genes,
					                         chrom[pi].genes, chrom[pj].genes,
					                         FA->map_par, cleftgrid);
					if (prmsp <= GB->sig_share) {
						pshare += (1.0 - pow((prmsp / GB->sig_share), GB->alpha));
					}
				}

				// Rank component: normalised to [0, 1].
				double rank_component = static_cast<double>(GB->num_chrom - pi) /
				                        static_cast<double>(GB->num_chrom);

				// Boltzmann component: normalised to [0, 1] by max weight.
				double boltz_component = chrom[pi].boltzmann_weight / max_bw;

				// Blended fitness divided by niche count.
				double blended = (1.0 - w) * rank_component + w * boltz_component;
				chrom[pi].fitnes = blended * static_cast<double>(GB->num_chrom) / pshare;
			}

			// Log ensemble thermodynamics periodically.
			if (gen_id % GA_SMFREE_LOG_INTERVAL == 0) {
				fprintf(stderr, "[SMFREE] gen=%d  F=%.3f  <E>=%.3f  S=%.6f  Cv=%.4f  σ_E=%.3f\n",
				        gen_id, thermo.free_energy, thermo.mean_energy,
				        thermo.entropy, thermo.heat_capacity, thermo.std_energy);
			}
		} else {
			// Temperature = 0: fall back to rank-only (same as LINEAR).
			for (i = 0; i < GB->num_chrom; i++) {
				chrom[i].fitnes = static_cast<double>(GB->num_chrom - i);
				chrom[i].boltzmann_weight = 0.0;
				chrom[i].free_energy = 0.0;
			}
		}
	}

	if(print){

		FILE* outfile_ptr = get_update_file_ptr(FA);

		if(outfile_ptr == NULL){
			fprintf(stderr,"ERROR: The NRGsuite failed to update within the timeout.\n");
			Terminate(10);
		}

		fprintf(outfile_ptr, "Generation: %5d\n", gen_id);
		fprintf(outfile_ptr, "best by energy\n");

		print_par(chrom,gene_lim,GB->num_print,GB->num_genes, outfile_ptr);

		fflush(outfile_ptr);

		close_update_file_ptr(FA, outfile_ptr);

	}

	gen_id++;

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void close_update_file_ptr(FA_Global* FA, FILE* outfile_ptr)
{

	if(FA->nrg_suite){
		fclose(outfile_ptr);
	}

}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
FILE* get_update_file_ptr(FA_Global* FA)
{

	if(!FA->nrg_suite){
		return stdout;
	}

	FILE* outfile_ptr = NULL;
	char UPDATEFILE[MAX_PATH__];
	long long timeout = 0;

#ifdef _WIN32
	snprintf(UPDATEFILE,MAX_PATH__,"%s\\.update",FA->state_path);
#else
	snprintf(UPDATEFILE,MAX_PATH__,"%s/.update",FA->state_path);
#endif

	outfile_ptr = fopen(UPDATEFILE,"r");
	if(outfile_ptr != NULL) {
		do {
			fclose(outfile_ptr);

# ifdef _WIN32
			Sleep(SLEEP);
# else
			usleep(SLEEP*1000);
# endif

			timeout += SLEEP;
			if(timeout >= FA->nrg_suite_timeout*1000){
				return NULL;
			}

			outfile_ptr = fopen(UPDATEFILE,"r");

		}while(outfile_ptr != NULL);
	}

	outfile_ptr = fopen(UPDATEFILE,"w");
	if(outfile_ptr == NULL){
		fprintf(stderr,"ERROR: Cannot open update file '%s' for reading.\n", UPDATEFILE);
		Terminate(10);
	}

	return outfile_ptr;

}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
cfstr eval_chromosome(FA_Global* FA,GB_Global* GB,VC_Global* VC,const genlim* gene_lim,
		      atom* atoms,resid* residue,gridpoint* cleftgrid,gene* john,
		      cfstr (*function)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*)){

	double icv[MAX_NUM_GENES] = {0};

	for(int i=0;i<GB->num_genes;i++){
		if(john[i].to_ic > gene_lim[i].max) {
			fprintf(stderr, "Exceptional out of bounds error at: max: %.5lf when ic: %.5lf\n", gene_lim[i].max, john[i].to_ic);
			john[i].to_ic = gene_lim[i].max;
		}else if(john[i].to_ic < gene_lim[i].min) {
			fprintf(stderr, "Exceptional out of bounds error at: min: %.5lf when ic: %.5lf\n", gene_lim[i].max, john[i].to_ic);
			john[i].to_ic = gene_lim[i].min;
		}

		icv[i] = john[i].to_ic;
	}

	return (*function)(FA,VC,atoms,residue,cleftgrid,GB->num_genes,icv);
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void generate_random_individual(FA_Global* FA, GB_Global* GB, atom* atoms, gene* genes, const genlim* gene_lim,
				std::function<int32_t()> & dice,
				int from_gene, int to_gene)
{
	for(int j=from_gene;j<to_gene;j++)
	{
		// side-chain optimization
		if(FA->map_par[j].typ == 4)
		{
			int l=0;
			while(FA->flex_res[l].inum != atoms[FA->map_par[j].atm].ofres){
				l++;
			};

			//printf("probability of atom[%d].ofres[%d]\t flex_res[%d](%s).inum[%d]= %.3f\n", FA->map_par[j].atm, atoms[FA->map_par[j].atm].ofres, l, FA->flex_res[l].name, FA->flex_res[l].inum, FA->flex_res[l].prob);

			if(RandomDouble() < FA->flex_res[l].prob)
			{
				genes[j].to_int32 = dice();
			}else{
				genes[j].to_int32 = 0;
			}
		}else{
			genes[j].to_int32 = dice();
		}

		genes[j].to_ic = genetoic(&gene_lim[j],genes[j].to_int32);
	}

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void populate_chromosomes(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,
                          atom* atoms,resid* residue,gridpoint* cleftgrid,char method[],
                          cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
                          char file[], long int at, int popoffset, int print,
                          std::function<int32_t()> & dice,
                          std::unordered_map<size_t, int> & duplicates){

	int i,j;

	/*
	  std::mt19937 rng;
	  std::uniform_int_distribution<int32_t> one_to_max_int32( 1, MAX_RANDOM_VALUE );
	  std::function<int32_t()> dice = [&](){ return one_to_max_int32(rng); };
	*/

	FILE* infile_ptr = NULL;

	// initialise genes to zero
	for(i=popoffset;i<GB->num_chrom;i++){
		for(j=0;j<GB->num_genes;j++){
			chrom[i].genes[j].to_int32=0;
			chrom[i].genes[j].to_ic=0.0;
		}
	}

	//------------------------------------------------------------------------------
	// use method to create new genes
	if(strcmp(method,"RANDOM")==0){
		printf("generating random population...\n");
		//printf("num_chrom=%d num_genes=%d\n",GB->num_chrom,GB->num_genes);

		int gener=0;
		size_t sig = 0;

		i=popoffset;
		while(i<GB->num_chrom){
			while(1){
				generate_random_individual(FA,GB,atoms,chrom[i].genes,gene_lim,dice,0,GB->num_genes);

				// ── MIF-weighted or RefLig seeding override for gene 0 ──
				if (FA->reflig_nearest_count > 0 &&
				    i < popoffset + static_cast<int>(FA->reflig_seed_fraction *
				        static_cast<float>(GB->num_chrom - popoffset))) {
					// RefLig seeding: distribute K nearest grid points across seeded fraction
					int k = std::abs(chrom[i].genes[0].to_int32) % FA->reflig_nearest_count;
					int grid_idx = FA->reflig_nearest_grid[k];
					chrom[i].genes[0].to_ic = static_cast<double>(grid_idx);
					chrom[i].genes[0].to_int32 = ictogene(&gene_lim[0],
					                                       static_cast<double>(grid_idx));
				} else if (FA->mif_enabled && FA->mif_cdf && FA->mif_count > 0) {
					// MIF-weighted Boltzmann sampling
					double u = RandomDouble(dice());
					auto it = std::lower_bound(FA->mif_cdf,
					                           FA->mif_cdf + FA->mif_count, u);
					int idx = static_cast<int>(std::distance(FA->mif_cdf, it));
					idx = std::clamp(idx, 0, FA->mif_count - 1);
					int grid_idx = FA->mif_sorted[idx];
					chrom[i].genes[0].to_ic = static_cast<double>(grid_idx);
					chrom[i].genes[0].to_int32 = ictogene(&gene_lim[0],
					                                       static_cast<double>(grid_idx));
				}

				sig = hash_genes(chrom[i].genes,GB->num_genes);
				if(GB->duplicates || duplicates.find(sig) == duplicates.end()){
					break;
				}
			}

			gener++;
			i++;
			duplicates[sig] = 1;
		}

		printf("generated %d randomized individuals\n", gener);

	}

	//------------------------------------------------------------------------------

	if(strcmp(method,"IPFILE")==0){
		printf("generating population from file...\n");

		if(!OpenFile_B(file,"rb",&infile_ptr)){
			fprintf(stderr,"ERROR: Cannot open file '%s' for reading.\n", file);
			Terminate(8);
		}

		fseek(infile_ptr, at, SEEK_SET);

		i=0;
		j=0;
		while(i<GB->num_chrom && fread(&chrom[i].genes[j].to_int32, 1, sizeof(int32_t), infile_ptr))
		{
			chrom[i].genes[j].to_ic = genetoic(&gene_lim[j],chrom[i].genes[j].to_int32);

			j++;
			if(j==GB->num_genes){
				i++;
				j=0;
			}
		}

		CloseFile_B(&infile_ptr,"r");

		printf("generated %d individuals from file\n", i);

		// reset to RANDOM afterwards
		strcpy(method,"RANDOM");

		// complete remaining population when necessary
		return populate_chromosomes(FA,GB,VC,chrom,gene_lim,atoms,residue,
					    cleftgrid,GB->pop_init_method,target,
					    GB->pop_init_file,at,i,print,dice,duplicates);
	}

	//------------------------------------------------------------------------------

	// calculate evalue for each chromosome — thread-safe OpenMP parallel eval.
	// Uses the same per-thread VC/atoms/FA workspace strategy as calculate_fitness.
	{
#ifdef _OPENMP
		const int n_thr = omp_get_max_threads();
#else
		const int n_thr = 1;
#endif
		const int natm  = FA->atm_cnt;
		const int natmr = FA->atm_cnt_real;
		const int nres  = FA->res_cnt;
		const int nopt  = FA->num_optres;
		const int nctb  = FA->ntypes * FA->ntypes;
		const int range = GB->num_chrom - popoffset;

		std::vector<std::vector<atom>>   p_atoms(n_thr, std::vector<atom>(atoms, atoms + natm + 1));
		std::vector<std::vector<resid>>  p_res(n_thr, std::vector<resid>(residue, residue + nres + 1));
		std::vector<FA_Global>           p_fa(n_thr, *FA);
		std::vector<std::vector<int>>    p_contacts(n_thr, std::vector<int>(MAX_ATOM_NUMBER, 0));
		std::vector<std::vector<float>>  p_contrib(n_thr, std::vector<float>(nctb, 0.0f));
		std::vector<std::vector<OptRes>> p_optres(n_thr,
		    std::vector<OptRes>(FA->optres, FA->optres + nopt));
		std::vector<VC_Global>               p_vc(n_thr, *VC);
		std::vector<std::vector<atomsas>>    p_calc(n_thr, std::vector<atomsas>(natmr));
		std::vector<std::vector<int>>        p_calclist(n_thr, std::vector<int>(natmr));
		std::vector<std::vector<int>>        p_caidx(n_thr, std::vector<int>(natmr, -1));
		std::vector<std::vector<ca_struct>>  p_carec(n_thr,
		    std::vector<ca_struct>(VC->ca_recsize));
		std::vector<std::vector<int>>        p_seed(n_thr, std::vector<int>(3 * natmr));
		std::vector<std::vector<contactlist>> p_contlist(n_thr, std::vector<contactlist>(GA_CONTLIST_SIZE));
		std::vector<std::vector<ptindex>>    p_ptorder(n_thr, std::vector<ptindex>(MAX_PT));
		std::vector<std::vector<vertex>>     p_centerpt(n_thr, std::vector<vertex>(MAX_PT));
		std::vector<std::vector<vertex>>     p_poly(n_thr, std::vector<vertex>(MAX_POLY));
		std::vector<std::vector<plane>>      p_cont(n_thr, std::vector<plane>(MAX_PT));
		std::vector<std::vector<edgevector>> p_vedge(n_thr, std::vector<edgevector>(MAX_POLY));

		for (int t = 0; t < n_thr; ++t) {
			p_fa[t].contacts      = p_contacts[t].data();
			p_fa[t].contributions = p_contrib[t].data();
			p_fa[t].optres        = p_optres[t].data();
			p_vc[t].Calc      = p_calc[t].data();
			p_vc[t].Calclist  = p_calclist[t].data();
			p_vc[t].ca_index  = p_caidx[t].data();
			p_vc[t].ca_rec    = p_carec[t].data();
			p_vc[t].seed      = p_seed[t].data();
			p_vc[t].contlist  = p_contlist[t].data();
			p_vc[t].ptorder   = p_ptorder[t].data();
			p_vc[t].centerpt  = p_centerpt[t].data();
			p_vc[t].poly      = p_poly[t].data();
			p_vc[t].cont      = p_cont[t].data();
			p_vc[t].vedge     = p_vedge[t].data();
		}

		(void)range;  // suppress unused warning when _OPENMP not defined

		// ── Dirty-tracking optimisation (same logic as main eval loop) ───
		bool p_has_normal_modes = false;
		for (int p = 0; p < FA->npar; ++p) {
			if (FA->map_par[p].typ == 3) { p_has_normal_modes = true; break; }
		}
		std::vector<int> p_dirty_atm;
		std::vector<int> p_dirty_res_idx;
		if (!p_has_normal_modes) {
			for (int r = 0; r < FA->nors; ++r)
				for (int m = 0; m < FA->nmov[r]; ++m)
					p_dirty_atm.push_back(FA->mov[r][m]);
			for (int p = 0; p < FA->npar; ++p)
				p_dirty_atm.push_back(FA->map_par[p].atm);
			for (int p = 0; p < FA->npar; ++p) {
				if (FA->map_par[p].typ == 2) {
					int j = FA->map_par[p].atm;
					int cat = atoms[j].rec[3];
					while (cat != 0 && cat != FA->map_par[p].atm) {
						p_dirty_atm.push_back(cat);
						j = cat;
						cat = atoms[j].rec[3];
					}
				}
			}
			std::sort(p_dirty_atm.begin(), p_dirty_atm.end());
			p_dirty_atm.erase(std::unique(p_dirty_atm.begin(), p_dirty_atm.end()),
			                   p_dirty_atm.end());
			for (int p = 0; p < FA->npar; ++p) {
				if (FA->map_par[p].typ == 4)
					p_dirty_res_idx.push_back(atoms[FA->map_par[p].atm].ofres);
			}
			std::sort(p_dirty_res_idx.begin(), p_dirty_res_idx.end());
			p_dirty_res_idx.erase(
				std::unique(p_dirty_res_idx.begin(), p_dirty_res_idx.end()),
				p_dirty_res_idx.end());
		}
		const bool p_use_selective = !p_has_normal_modes &&
		    static_cast<int>(p_dirty_atm.size()) < natm / 2;
		const int p_n_dirty_atm = static_cast<int>(p_dirty_atm.size());
		const int p_n_dirty_res = static_cast<int>(p_dirty_res_idx.size());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(chrom, FA, GB, VC, gene_lim, atoms, residue, cleftgrid, target, \
	       popoffset, p_atoms, p_res, p_fa, p_optres, p_vc, natm, nres, nopt, \
	       p_use_selective, p_dirty_atm, p_dirty_res_idx, p_n_dirty_atm, p_n_dirty_res)
#endif
		for(i=popoffset;i<GB->num_chrom;i++){
#ifdef _OPENMP
			const int tid = omp_get_thread_num();
#else
			const int tid = 0;
#endif
			if (p_use_selective) {
				for (int d = 0; d < p_n_dirty_atm; ++d) {
					const int ai = p_dirty_atm[d];
					p_atoms[tid][ai] = atoms[ai];
				}
				for (int d = 0; d < p_n_dirty_res; ++d) {
					const int ri = p_dirty_res_idx[d];
					p_res[tid][ri] = residue[ri];
				}
			} else {
				std::copy(atoms,   atoms + natm + 1,   p_atoms[tid].begin());
				std::copy(residue, residue + nres + 1, p_res[tid].begin());
			}
			// Redirect per-thread atom optres pointers to per-thread optres array.
			for (int ai = 1; ai <= natm; ++ai) {
				atom& a = p_atoms[tid][ai];
				if (a.optres) {
					ptrdiff_t oidx = a.optres - FA->optres;
					a.optres = &p_optres[tid][oidx];
				}
			}
			for (int o = 0; o < nopt; ++o) {
				p_optres[tid][o].cf.com    = 0.0;
				p_optres[tid][o].cf.wal    = 0.0;
				p_optres[tid][o].cf.sas    = 0.0;
				p_optres[tid][o].cf.totsas = 0.0;
				p_optres[tid][o].cf.con    = 0.0;
				p_optres[tid][o].cf.gist   = 0.0;
				p_optres[tid][o].cf.elec   = 0.0;
				p_optres[tid][o].cf.hbond  = 0.0;
				p_optres[tid][o].cf.gist_desolv = 0.0;
				p_optres[tid][o].cf.rclash = 0;
			}
			p_vc[tid].numcarec = 0;

			chrom[i].cf = eval_chromosome(
			    &p_fa[tid], GB, &p_vc[tid], gene_lim,
			    p_atoms[tid].data(), p_res[tid].data(),
			    cleftgrid, chrom[i].genes, target);
			chrom[i].evalue     = get_cf_evalue(&chrom[i].cf);
			chrom[i].app_evalue = get_apparent_cf_evalue(&chrom[i].cf);
			chrom[i].status     = 'n';
			ccbm_inject_strain(FA, chrom[i], gene_lim);  // CCBM strain
		}
	}

	// sort and calculate fitness (use a local GAContext for initial population)
	GAContext pop_ctx;
	calculate_fitness(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,GB->fitness_model,GB->num_chrom,print,target,pop_ctx);

	return;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int cmp_chrom2rotlist(psFlexDEE_Node psFlexDEE_INI_Node, const chromosome* chrom, const genlim* gene_lim,
                      int gene_offset, int num_genes, int tot, int num_nodes){

	int   par[GA_MAX_FLEXDEE_PARAMS];
	//int* genes = NULL;
	sFlexDEE_Node sFlexDEENode;

	memset(&par,0,sizeof(par));

	if ( psFlexDEE_INI_Node == NULL ) { return 0; }


	sFlexDEENode.rotlist = par;

	for(int i=0;i<tot;i++){
		//genes = &chrom[i].genesic[gene_offset];

		psFlexDEE_INI_Node = psFlexDEE_INI_Node->last;

		if ( dee_pivot(&sFlexDEENode,&psFlexDEE_INI_Node,1,num_nodes,(num_nodes+1)/2,num_nodes,num_genes) == 0 ) { return 1; }

	}

	return 0;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int cmp_chrom2pop(const chromosome* chrom,const gene* genes, int num_genes,int start, int last){
	int i,j,flag;

	for(i=start;i<last;i++){
		flag=0;
		for(j=0;j<num_genes;j++){
			//printf("individuals[%d][%d].gene[%d]=%.3f\t%.3f\n", start-1, i, j,
			//       genes[j].to_ic, chrom[i].genes[j].to_ic);
			flag += abs(genes[j].to_ic - chrom[i].genes[j].to_ic) < GA_GENE_MATCH_TOLERANCE;
		}

		//printf("flag=%d\n",flag);
		if(flag == num_genes){return 1;}
	}

	return 0;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int cmp_chrom2pop_int(const chromosome* chrom,const gene* genes, int num_genes,int start, int last){
	int i,j,flag;

	for(i=start;i<last;i++){
		flag=0;
		for(j=0;j<num_genes;j++){
			//printf("comparing %u to %u\n",c->genes[j],chrom[i].genes[j]);
			flag += ( genes[j].to_int32 == chrom[i].genes[j].to_int32 );
		}

		//printf("flag=%d\n",flag);
		if(flag == num_genes){return 1;}
	}

	return 0;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void validate_dups(GB_Global* GB, genlim* gene_lim, int num_genes){

	double n_poss = calc_poss(gene_lim, num_genes);

	if(n_poss < (double)GB->num_chrom && !GB->duplicates){
		fprintf(stderr,"Too many chromosomes for the number of possibilites (%.1lf) when no duplicates allowed.\n", n_poss);
		fprintf(stderr,"Duplicates are then allowed.\n");
		GB->duplicates = 1;
	}

	return;
}

double calc_poss(genlim* gene_lim, int num_genes){

	double n_poss = 0.0;

	for(int i=0; i<num_genes; i++){
		if(n_poss > 0.0){
			n_poss *= gene_lim[i].nbin;
		}else{
			n_poss = gene_lim[i].nbin;
		}
	}

	return n_poss;
}

void set_bins(genlim* gene_lim, int num_genes){

	for(int i=0; i<num_genes; i++){
		double nbin = (gene_lim[i].max - gene_lim[i].min) / gene_lim[i].del;
		if(nbin - (int)nbin > 0.0){ nbin += 1.0; }
		if(gene_lim[i].map){ nbin += 1.0; }

		gene_lim[i].bin = 1.0/nbin;
		gene_lim[i].nbin = nbin;
	}

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void set_bins(genlim* gene_lim){

	double nbin = (gene_lim->max - gene_lim->min) / gene_lim->del;
	if(nbin - (int)nbin > 0.0){ nbin += 1.0; }
	if(gene_lim->map){ nbin += 1.0; }

	gene_lim->bin = 1.0/nbin;

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void read_gainputs(FA_Global* FA,GB_Global* GB,int* gen_int,int* sz_part,char file[]){

	FILE *infile_ptr;        /* pointer to input file */
	char buffer[MAX_PATH__];         /* a line from the INPUT file */
	char field[9];           /* field names on INPUT file */

	// Direct mode: GA params already set by apply_config — skip file reading
	if(file[0] == '\0'){
		printf("read_gainputs: using pre-configured GA parameters (direct mode)\n");
		return;
	}

	//printf("file here is <%s>\n",file);
	// In direct mode (no .ga.inp file), all GA params are set via
	// apply_config().  Skip file parsing when the path is empty.
	if (file[0] == '\0') {
		printf("read_gainputs: no GA input file — using config defaults\n");
		return;
	}
	infile_ptr=NULL;
	if(!OpenFile_B(file,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Cannot find file '%s'.\n", file);
		Terminate(8);
	}

	while (fgets(buffer, sizeof(buffer),infile_ptr)){
		size_t blen = strlen(buffer);
		if (blen > 0 && buffer[blen-1] == '\n')
			buffer[--blen] = '\0';
		if (blen > 0 && buffer[blen-1] == '\r')
			buffer[--blen] = '\0';


		if(strncmp(buffer,"NUMCHROM",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->num_chrom);
		}else if(strncmp(buffer,"OPTIGRID",8) == 0){
			sscanf(buffer,"%s %d %d %d",field,&FA->opt_grid,gen_int,sz_part);
		}else if(strncmp(buffer,"NUMGENER",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->max_generations);
		}else if(strncmp(buffer,"ADAPTVGA",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->adaptive_ga);
		}else if(strncmp(buffer,"ADAPTKCO",8) == 0){
			//adaptive response parameters
			//k1-k4 are values ranging from 0.0-1.0 inclusively
			sscanf(buffer,"%s %lf %lf %lf %lf",field,&GB->k1,&GB->k2,&GB->k3,&GB->k4);
		}else if(strncmp(buffer,"CROSRATE",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->cross_rate);
		}else if(strncmp(buffer,"MUTARATE",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->mut_rate);
		}else if(strncmp(buffer,"INTRAGEN",8) == 0){
			GB->intragenes = 1;
		}else if(strncmp(buffer,"INIMPROB",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->ini_mut_prob);
		}else if(strncmp(buffer,"ENDMPROB",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->end_mut_prob);
		}else if(strncmp(buffer,"POPINIMT",8) == 0){
			sscanf(buffer,"%s %8s",field,GB->pop_init_method);
			//0         1         2
			//012345678901234567890123456789
			//POPINIMT IPFILE file.dat
			if(strcmp(GB->pop_init_method,"IPFILE") == 0 && blen > 16){
				strncpy(GB->pop_init_file,&buffer[16],MAX_PATH__-1);
				GB->pop_init_file[MAX_PATH__-1]='\0';
			}
		}else if(strncmp(buffer,"FITMODEL",8) == 0){
			sscanf(buffer,"%s %8s",field,GB->fitness_model);
		}else if(strncmp(buffer,"REPMODEL",8) == 0){
			sscanf(buffer,"%s %8s",field,GB->rep_model);
		}else if(strncmp(buffer,"DUPLICAT",8) == 0){
			GB->duplicates = 1;
		}else if(strncmp(buffer,"BOOMFRAC",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->pbfrac);
		}else if(strncmp(buffer,"STEADNUM",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->ssnum);
		}else if(strncmp(buffer,"SHAREALF",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->alpha);
		}else if(strncmp(buffer,"SHAREPEK",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->peaks);
		}else if(strncmp(buffer,"SHARESCL",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->scale);
		}else if(strncmp(buffer,"OUTGENER",8) == 0){
			GB->outgen = 1;
		}else if(strncmp(buffer,"STRTSEED",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->seed);
		}else if(strncmp(buffer,"PRINTCHR",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->num_print);
		}else if(strncmp(buffer,"PRINTINT",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->print_int);
		}else if(strncmp(buffer,"PRINTRRG",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->rrg_skip);
		}else if(strncmp(buffer,"ENTRCNVG",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->entropy_convergence);
		}else if(strncmp(buffer,"ENTRCHKI",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->entropy_check_interval);
		}else if(strncmp(buffer,"ENTRWIND",8) == 0){
			sscanf(buffer,"%s %d",field,&GB->entropy_window);
		}else if(strncmp(buffer,"ENTRTHRS",8) == 0){
			sscanf(buffer,"%s %lf",field,&GB->entropy_rel_threshold);
		}else if(strncmp(buffer,"MIFWEIGH",8) == 0){
			sscanf(buffer,"%*s %d", &FA->mif_enabled);
		}else if(strncmp(buffer,"MIFTEMPR",8) == 0){
			sscanf(buffer,"%*s %f", &FA->mif_temperature);
		}else if(strncmp(buffer,"GRIDPRIO",8) == 0){
			sscanf(buffer,"%*s %f", &FA->grid_prio_percent);
		}else if(strncmp(buffer,"REFLGFIL",8) == 0){
			sscanf(buffer,"%*s %s", FA->reflig_file);
		}else if(strncmp(buffer,"REFLGSED",8) == 0){
			sscanf(buffer,"%*s %f", &FA->reflig_seed_fraction);
		}else if(strncmp(buffer,"REFLGKNN",8) == 0){
			sscanf(buffer,"%*s %d", &FA->reflig_k_nearest);
		}else if(strncmp(buffer,"REFLGHTM",8) == 0){
			sscanf(buffer,"%*s %d", &FA->reflig_hetatm_fallback);
		}else if(strncmp(buffer,"AUTOFLXE",8) == 0){
			sscanf(buffer,"%*s %d", &FA->autoflex_enabled);
		}else if(strncmp(buffer,"AUTOFLXN",8) == 0){
			sscanf(buffer,"%*s %d", &FA->autoflex_max);
		}else{
			// ...
		}

	}

	CloseFile_B(&infile_ptr,"r");

}

long int read_pop_init_file(FA_Global* FA, GB_Global* GB, genlim* gene_lim, char* pop_init_file)
{

	long int at = 0;
	FILE* infile_ptr = NULL;

	if(!OpenFile_B(pop_init_file,"rb",&infile_ptr)){
		fprintf(stderr,"ERROR: Cannot open file '%s' for reading.\n", pop_init_file);
		Terminate(8);
	}

	char genes_tag[6];
	fread(&genes_tag[0], 1, sizeof(genes_tag)-1, infile_ptr);
	genes_tag[5] = '\0';
	//printf("genes_tag=%s\n", genes_tag);

	if(strcmp(genes_tag,"genes") == 0){

		int i=0;
		while(i < GB->num_genes){
			fread(&gene_lim[i], 1, sizeof(genlim), infile_ptr);
			i++;
		}

		char chrom_tag[6];
		fread(&chrom_tag[0], 1, sizeof(chrom_tag)-1, infile_ptr);
		chrom_tag[5] = '\0';
		//printf("chrom_tag=%s\n", chrom_tag);

		if(strcmp(chrom_tag,"chrom") == 0){
			at = ftell(infile_ptr);
		}

	}

	CloseFile_B(&infile_ptr, "r");

	return at;
}

void set_gene_lim(FA_Global* FA, GB_Global* GB, genlim* gene_lim)
{

	for(int ngenes=0;ngenes<GB->num_genes;ngenes++){
		gene_lim[ngenes].min=FA->min_opt_par[ngenes];
		gene_lim[ngenes].max=FA->max_opt_par[ngenes];
		gene_lim[ngenes].del=FA->del_opt_par[ngenes];
		gene_lim[ngenes].map=FA->map_opt_par[ngenes];

		printf("gene %d: min: %10.2f max: %10.2f delta: %10.2f map: %d\n", ngenes,
		       gene_lim[ngenes].min,
		       gene_lim[ngenes].max,
		       gene_lim[ngenes].del,
		       gene_lim[ngenes].map);

	}

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void crossover(gene *john,gene *mary,int num_genes, int intragenes){

	/* john and mary are two chromosomes to be crossover at points a and b
	 */

	int i,j;
	unsigned int optr;
	int temp;
	int gen_a,gen_b,aux_gen;
	int pnt_a,pnt_b,aux_pnt;

	gen_a=(int)(RandomDouble()*(double)num_genes);
	gen_b=(int)(RandomDouble()*(double)num_genes);
	if (gen_a >= num_genes) gen_a = num_genes - 1;
	if (gen_b >= num_genes) gen_b = num_genes - 1;
	//printf("gen_a=%d\tgen_b=%d\n",gen_a,gen_b);



	if(intragenes){
		pnt_a=(int)(RandomDouble()*(double)(MAX_GEN_LENGTH));
		pnt_b=(int)(RandomDouble()*(double)(MAX_GEN_LENGTH));

		if(gen_a > gen_b){
			aux_gen=gen_a;
			aux_pnt=pnt_a;
			gen_a=gen_b;
			pnt_a=pnt_b;
			gen_b=aux_gen;
			pnt_b=aux_pnt;
		}

		if(gen_a == gen_b && pnt_a < pnt_b){
			aux_pnt=pnt_b;
			pnt_b=pnt_a;
			pnt_a=aux_pnt;
		}
	}else{
		if(gen_a > gen_b){
			aux_gen=gen_a;
			gen_a=gen_b;
			gen_b=aux_gen;
		}

		if(gen_a != gen_b){
			// find left of right bound of gene a
			if(RandomDouble() < 0.5){
				pnt_a=MAX_GEN_LENGTH;
			}else{
				pnt_a=0;
			}

			// find left of right bound of gene b
			if(RandomDouble() < 0.5){
				pnt_b=MAX_GEN_LENGTH;
			}else{
				pnt_b=0;
			}

			if((gen_b - gen_a) == 1 && pnt_a == 0 && pnt_b == MAX_GEN_LENGTH){
				pnt_b=pnt_a;
			}
		}else{
			pnt_a=MAX_GEN_LENGTH; pnt_b=0;
		}
	}

	//printf("gen_a=%d\tpnt_a=%d\tgen_b=%d\tpnt_b=%d\n",gen_a,pnt_a,gen_b,pnt_b);

	for(j=gen_a;j<=gen_b;j++){
		optr=1u;
		aux_pnt = (j==gen_a)?pnt_a:(MAX_GEN_LENGTH);
		for(i=0;i<aux_pnt;i++) optr |= (optr << 1);
		unsigned int uj = static_cast<unsigned int>(john[j].to_int32);
		unsigned int um = static_cast<unsigned int>(mary[j].to_int32);
		john[j].to_int32 = static_cast<int32_t>((uj & ~optr) | (um &  optr));
		mary[j].to_int32 = static_cast<int32_t>((uj &  optr) | (um & ~optr));
	}

	if(pnt_b > 0){
		optr=1u;
		for(i=0;i<pnt_b-1;i++) optr |= (optr << 1);
		unsigned int uj = static_cast<unsigned int>(john[gen_b].to_int32);
		unsigned int um = static_cast<unsigned int>(mary[gen_b].to_int32);
		john[gen_b].to_int32 = static_cast<int32_t>((uj & ~optr) | (um &  optr));
		mary[gen_b].to_int32 = static_cast<int32_t>((uj &  optr) | (um & ~optr));
	}

	/*
	  printf("john after:\n");
	  print_chrom(john,num_genes,0);
	  printf("mary after:\n");
	  print_chrom(mary,num_genes,0);
	*/

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void mutate(gene *john,int num_genes,double mut_rate){
	/* creates an operator with 1's with rate= mut_rate
	   uses it to mutate john.
	*/
	int i,j;
	unsigned int optr;
	unsigned int test;

	for(j=0;j<num_genes;j++){
		optr=0u;
		test=1u;
		for(i=0;i<32;i++){
			if(RandomDouble() < mut_rate){
				optr |= test;
			}
			test <<= 1;
		}
		john[j].to_int32 ^= static_cast<int32_t>(optr);
	}

	return;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void bin_print(int dec, int len){
	int i,val;
	int test=0;
	int op=1;
	op <<= len-1;
	//printf("op=%u\n",op);
	//printf("dec=%u len=%d\n",dec,len);
	for(i=len-1;i>=0;i--){
		test = (int)pow(2.0f,i);
		//printf("\n[%u]&[%u]=%u test=%u: ",dec,op,dec&op,test);
		val=0;
		if((dec&op) == test) val=1;
		printf("%1d",val);
		op >>= 1;
	}
	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_chrom(chromosome *x, chromosome *y){
	chromosome t=*x;*x=*y;*y=t;
}

void QuickSort(chromosome* list, int beg, int end, bool energy)
{
    QS_TYPE piv;

    int  l,r,p;

    while (beg<end)    // This while loop will avoid the second recursive call
    {
        l = beg; p = beg + (end-beg)/2; r = end;

		if(energy)
			piv = list[p].evalue;
		else
			piv = list[p].fitnes;

        while (1)
        {
            while ( (l<=r) && ( ( energy && QS_ASC(list[l].evalue,piv) <= 0 ) ||
								( !energy && QS_DSC(list[l].fitnes,piv) <= 0 ) ) ) l++;
            while ( (l<=r) && ( ( energy && QS_ASC(list[r].evalue,piv) > 0 ) ||
								( !energy && QS_DSC(list[r].fitnes,piv) > 0 ) ) ) r--;

            if (l>r) break;

			swap_chrom(&list[l],&list[r]);

            if (p==r) p=l;

            l++; r--;
        }

		swap_chrom(&list[p],&list[r]);
        //list[p]=list[r]; list[r].evalue=piv;
        r--;

        // Recursion on the shorter side & loop (with new indexes) on the longer
        if ((r-beg)<(end-l))
        {
            QuickSort(list, beg, r, energy);
            beg=l;
        }
        else
        {
            QuickSort(list, l, end, energy);
            end=r;
        }
    }
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int remove_dups(chromosome* chrom, int num_chrom, int num_genes){

	int i=0;
	int j;
	if (num_chrom<=1) return num_chrom;

	for (j=1;j<num_chrom;j++)
	{
		int flag = 0;
		for(int l=0;l<num_genes;l++){
			flag += abs(chrom[j].genes[l].to_ic - chrom[i].genes[l].to_ic) < 0.1;
		}
		if(flag != num_genes)
		{
			copy_chrom(&chrom[++i],&chrom[j],num_genes);
		}
	}

	return i+1;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void print_par(const chromosome* chrom,const genlim* gene_lim,int num_chrom,int num_genes, FILE* outfile_ptr){
	for(int i=0;i<num_chrom;i++){
		fprintf(outfile_ptr, "%4d (",i);
		for(int j=0;j<num_genes;j++) fprintf(outfile_ptr, "%10.2f ", chrom[i].genes[j].to_ic);
		fprintf(outfile_ptr, ") ");
		fprintf(outfile_ptr, " cf=%9.3f cf.app=%9.3f fitnes=%9.3f\n",
			chrom[i].evalue, chrom[i].app_evalue, chrom[i].fitnes);
	}

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void write_par(const chromosome* chrom,const genlim* gene_lim,int ger,char* outfile,int num_chrom,int num_genes){
	int i,j;
	FILE *outfile_ptr;

	outfile_ptr=NULL;
	if(!OpenFile_B(outfile,"wb",&outfile_ptr)){
		Terminate(6);
	}else{

		char genes_tag[5] = { 'g' , 'e' , 'n' , 'e' , 's' };

		fwrite(&genes_tag[0], 1, sizeof(genes_tag), outfile_ptr);
		for(j=0;j<num_genes;j++){
			fwrite(&gene_lim[j], 1, sizeof(genlim), outfile_ptr);
		}

		char chrom_tag[5] = { 'c' , 'h' , 'r' , 'o' , 'm' };

		fwrite(&chrom_tag[0], 1, sizeof(chrom_tag), outfile_ptr);
		for(i=0;i<num_chrom;i++)
		{
			for(j=0;j<num_genes;j++)
			{
				fwrite(&chrom[i].genes[j].to_int32, 1, sizeof(int32_t), outfile_ptr);
			}
		}

	}

	CloseFile_B(&outfile_ptr,"w");

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void print_pop(const chromosome* chrom,const genlim* gene_lim,int numc, int numg){
	int i,j;

	for(i=0;i<numc;i++){
		printf("%2d (",i);
		for(j=0;j<numg;j++){printf(" %10d",chrom[i].genes[j].to_int32);}
		printf(") ");
		for(j=0;j<numg;j++){printf(" "),bin_print(chrom[i].genes[j].to_int32,(MAX_GEN_LENGTH));}
		printf("\n");
	}
	return;
}
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void print_chrom(const chromosome* chrom, int num_genes, int real_flag){
	int j;
	//int i;

	printf("(");
	for(j=0;j<num_genes;j++){
		if(real_flag){
			printf(" %10.5f",chrom->genes[j].to_ic);
		}else{
			printf(" %10d",chrom->genes[j].to_int32);
		}
	}
	printf(") ");
	printf("\n");

	return;
}

/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
void print_chrom(const gene* genes, int num_genes, int real_flag){
	int j;
	//int i;

	printf("(");
	for(j=0;j<num_genes;j++){
		if(real_flag){
			printf(" %10.5f",genes[j].to_ic);
		}else{
			printf(" %10d",genes[j].to_int32);
		}
	}
	printf(") ");
	printf("\n");

	return;
}

/********************************************************************************
 * This function calculates the RSMD between atomic coordinates of the atoms in *
 * the register ori_ligatm and those for the atoms of the ligand in             *
 * residue[opt_res[0]] after reconstructing the coordinates using opt_par       *
 ********************************************************************************/

double calc_rmsp(int npar, const gene* g1, const gene* g2, const optmap* map_par, gridpoint* cleftgrid){
	// Vectorised RMSP using Eigen strided Map over the to_ic field.
	// gene_struct lays out {int32_t to_int32; double to_ic}, so stride = sizeof(gene).
	using EMap = Eigen::Map<const Eigen::VectorXd,
	                        Eigen::Unaligned,
	                        Eigen::InnerStride<sizeof(gene)/sizeof(double)>>;
	// to_ic is the second field; offset by one int32_t (4 bytes / 8 bytes per double = 0.5 — not aligned).
	// Fall back to a plain gather instead to avoid alignment issues.
	Eigen::VectorXd diff(npar);
	for (int ii = 0; ii < npar; ++ii) diff[ii] = g1[ii].to_ic - g2[ii].to_ic;
	return std::sqrt(diff.squaredNorm() / (double)npar);
}

double genetoic(const genlim* gene_lim, int32_t gene){

	int i=0;
	double tot=gene_lim->bin;

	while(tot < RandomDouble(gene))
	{
		tot += gene_lim->bin;
		i++;
	}

	double ic = gene_lim->min + gene_lim->del * (double)i;

	/* printf("ic=%.1f gene=%d randdouble=%.8f min=%.3f del=%.3f bin=%.8f\n",
		  ic, gene, RandomDouble(gene),
		  gene_lim->min, gene_lim->del, gene_lim->bin);
	*/

	return(ic);
}

int ictogene(const genlim* gene_lim, double ic){

	int i = (int)((ic - gene_lim->min) / gene_lim->del);

	double tot = 1.0;

	while(i > 0){
		tot -= gene_lim->bin;
		i--;
	}

	int gene = RandomInt(tot);

        /*
	  printf("ic=%.3f gene=%d randdouble=%.5f min=%.3f del=%.3f bin=%.3f\n",
	  ic, gene, RandomDouble(gene),
	  gene_lim->min, gene_lim->del, gene_lim->bin);
	*/

	return(gene);
}


int RandomInt(double frac){
	double raw = frac * ((double)RAND_MAX + 1.0);
	if (raw >= (double)RAND_MAX + 1.0) return RAND_MAX;
	if (raw < 0.0) return 0;
	return (int)raw;
}

double RandomDouble(int32_t gene){
	return gene/((double)MAX_RANDOM_VALUE+1.0);
}

double RandomDouble(){
	// Thread-safe RNG (replaces non-reentrant rand())
	thread_local std::mt19937 tl_rng(std::random_device{}());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(tl_rng);
}