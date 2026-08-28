#include "gaboom.h"
#include "Vcontacts.h"
#include "fileio.h"
#include "coarse_init.h"
#include "flexaid_exception.h"
#include "ga_constants.h"
#include "UnifiedHardwareDispatch.h"
#include "MIFGrid.h"
#include "CavityDetect/SpatialGrid.h"
#include "RngSeed.h"
#include "ensemble_pipeline.h"
#include "ProtocolConfig.h"
#include "niche_distance.h"
#include "niche_hash.h"
#include "new_search_arch.h"

#include <random>
#include <functional>
#include <climits>     // SIZE_MAX
#include <cstdint>
#include <vector>
#include <algorithm>
#include <array>
#include <memory>
#include <span>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <cstdlib>
#include <cstring>

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
#include "ThermalExtrapolation.h"
#include "ga_diversity.h"
#include "TurboQuant.h"
#include "GAContext.h"
#include "GPUContextPool.h"
#include "fast_optics.hpp"
#include "NATURaL/NATURaLDualAssembly.h"
#include "InStreamClustering.h"
#include "VibEntropy.h"
#include "LigandRingFlex/LigandRingFlex.h"   // Phase 2: ring pucker GA genes

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

// ═══ Ring pucker side-channel GA operators (LigandRingFlex Phase 2) ═══════════
// These operate directly on the chromosome's POD ring arrays (ring_phases /
// ring_six / ring_five). Counts come from FA->ring_flex_template (per-complex,
// detected at read time). All are no-ops when ring flex is inactive, so the
// default-OFF path and ring-free ligands are untouched. Conformer-index genes
// (six/five) are carried for completeness; the module's apply() currently
// deforms only furanose sugar rings (6-ring conformer apply is a documented
// no-op in LigandRingFlex.cpp), so only ring_phases changes Cartesian coords.
static inline void ring_randomise_chrom(const FA_Global* FA, chromosome* c) {
    if (!FA->ring_flex_active || !FA->ring_flex_template) return;
    const ligand_ring_flex::RingFlexGenes& t = *FA->ring_flex_template;
    auto& lib = ring_flex::RingConformerLibrary::instance();
    const int n6 = static_cast<int>(t.conformer_indices.size());
    const int n5 = static_cast<int>(t.five_conformer_indices.size());
    const int ns = FA->ring_n_sugars;
    for (int i = 0; i < ns && i < MAX_RING_FLEX; ++i)
        c->ring_phases[i] = static_cast<float>(RandomDouble() * 360.0);
    for (int i = 0; i < n6 && i < MAX_RING_FLEX; ++i)
        c->ring_six[i] = static_cast<uint8_t>(RandomDouble() * lib.n_six());
    for (int i = 0; i < n5 && i < MAX_RING_FLEX; ++i)
        c->ring_five[i] = static_cast<uint8_t>(RandomDouble() * lib.n_five());
}

static inline void ring_mutate_chrom(const FA_Global* FA, chromosome* c,
                                     double ring_mut_prob = 0.05,
                                     double pucker_mut_prob = 0.12) {
    if (!FA->ring_flex_active || !FA->ring_flex_template) return;
    const ligand_ring_flex::RingFlexGenes& t = *FA->ring_flex_template;
    auto& lib = ring_flex::RingConformerLibrary::instance();
    const int n6 = static_cast<int>(t.conformer_indices.size());
    const int n5 = static_cast<int>(t.five_conformer_indices.size());
    const int ns = FA->ring_n_sugars;
    for (int i = 0; i < ns && i < MAX_RING_FLEX; ++i)
        if (RandomDouble() < pucker_mut_prob)
            c->ring_phases[i] = sugar_pucker::mutate_phase(c->ring_phases[i]);
    for (int i = 0; i < n6 && i < MAX_RING_FLEX; ++i)
        if (RandomDouble() < ring_mut_prob)
            c->ring_six[i] = static_cast<uint8_t>(RandomDouble() * lib.n_six());
    for (int i = 0; i < n5 && i < MAX_RING_FLEX; ++i)
        if (RandomDouble() < ring_mut_prob)
            c->ring_five[i] = static_cast<uint8_t>(RandomDouble() * lib.n_five());
}

// Single-point crossover swapping the tail of each ring-gene array between two
// children (mirrors the standard gene crossover that produced them).
static inline void ring_crossover_chrom(const FA_Global* FA,
                                        chromosome* a, chromosome* b) {
    if (!FA->ring_flex_active || !FA->ring_flex_template) return;
    const ligand_ring_flex::RingFlexGenes& t = *FA->ring_flex_template;
    const int n6 = static_cast<int>(t.conformer_indices.size());
    const int ns = FA->ring_n_sugars;
    if (ns > 1) {
        int pt = 1 + static_cast<int>(RandomDouble() * (ns - 1));
        for (int i = pt; i < ns && i < MAX_RING_FLEX; ++i)
            std::swap(a->ring_phases[i], b->ring_phases[i]);
    }
    if (n6 > 1) {
        int pt = 1 + static_cast<int>(RandomDouble() * (n6 - 1));
        for (int i = pt; i < n6 && i < MAX_RING_FLEX; ++i)
            std::swap(a->ring_six[i], b->ring_six[i]);
    }
}

// Copy a chromosome's furanose pucker phases into FA->ring_cur_phases so the
// next ic2cf() reconstructs that chromosome's puckered ring. Must be called on
// the FA instance ic2cf will see (the per-thread FA copy in the GA eval loops).
static inline void ring_load_chrom_to_fa(FA_Global* FA, const chromosome* c) {
    if (!FA->ring_flex_active) return;
    const int ns = FA->ring_n_sugars;
    for (int i = 0; i < ns && i < MAX_RING_FLEX; ++i)
        FA->ring_cur_phases[i] = c->ring_phases[i];
}

// ── receptor_chain_normalizer ───────────────────────────────────────────────
// Experimental diagnostic only. Blindly dividing GA evalue/app_evalue by the
// number of receptor chains changes the effective selection temperature and
// regresses multichain Astex cases. Keep production/default scoring extensive;
// enable this only with FLEXAIDDS_CHAIN_NORM=1 for targeted ablations.
static int count_receptor_chains(FA_Global* FA, const resid* residue,
                                 bool chain_norm_enabled = false) {
    // When chain_norm_enabled is false (default), preserve the historical
    // extensive scoring path (return 1). Callers that want the diagnostic
    // multi-chain normalisation pass ProtocolConfig::chain_norm.
    if (!chain_norm_enabled) return 1;

    std::unordered_set<char> seen;
    for (int r = 0; r < FA->res_cnt; ++r) {
        if (residue[r].type == 0 && residue[r].chn != '\0' && residue[r].chn != ' ')
            seen.insert(residue[r].chn);
    }
    return seen.empty() ? 1 : static_cast<int>(seen.size());
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

	// Typed protocol snapshot (env remains the compatibility adapter).
	const flexaids::ProtocolConfig proto = flexaids::ProtocolConfig::from_env();

	// ── OMP thread default: 2/worker if OMP_NUM_THREADS not set in environment.
	// Leaves headroom for Metal dispatch and OS scheduling on M3 Pro 11 P-cores.
	// OMP_NUM_THREADS is a platform OpenMP contract — left as raw getenv.
#ifdef _OPENMP
	if (!std::getenv("OMP_NUM_THREADS")) {
		omp_set_num_threads(2);
	}
#endif

	int i;
	int print=0;

	// ── Experimental multi-chain VCT normalisation ────────────────────────────
	// Default is 1. FLEXAIDDS_CHAIN_NORM=1 enables a diagnostic ablation that
	// divides every evalue/app_evalue by the number of explicit receptor chain
	// IDs. Do not enable for production benchmarks without a target-specific
	// justification: it changes GA selection pressure and regressed Astex smoke.
	const int n_receptor_chains = count_receptor_chains(FA, residue, proto.chain_norm);
	if (n_receptor_chains > 1)
		fprintf(stderr, "[CHAIN] %d receptor chains detected — "
		        "normalising VCT evalue by chain count (diagnostic)\n", n_receptor_chains);

	// ── Level-3 H(ω) diagnostic env override ──────────────────────────────────
	// FLEXAIDDS_USE_SHANNON=1 enables the ligand vibrational-mode Shannon-entropy
	// monitor (per-generation [HVIB] lines) even on the bare binary, mirroring the
	// DatasetRunner config toggle.  Engine-side env wins, matching FLEXAIDDS_N_ELITE.
	// Purely diagnostic — does NOT enter CF or fitness.  Default stays OFF.
	if (proto.use_shannon) {
		if (!GB->use_shannon)
			fprintf(stderr, "[HVIB] FLEXAIDDS_USE_SHANNON set: enabling ligand "
			        "vibrational-entropy H(ω) monitor (diagnostic only)\n");
		GB->use_shannon = 1;
	}

	//char tmp_rrgfile[MAX_PATH__];
	//int rrg_flag;
	//int rrg_skip=100;

	char outfile[MAX_PATH__];
	int n_chrom_snapshot=0;
	char gridfile[MAX_PATH__];
	char gridfilename[MAX_PATH__];
	(void)gridfilename; // reserved for future grid-naming use

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

	// Entropy convergence — enabled by default; ENTRCNVG=0 in gainp to disable.
	// These MUST be set unconditionally: the no-GA-input-file path (JSON
	// --config benchmark runs) does not read a gainp file, so leaving
	// entropy_check_interval at 0 caused an integer modulo-by-zero (SIGFPE)
	// at generation 0 in the H-plateau / entropy-convergence checks below.
	// The gainp path may still override these values via read_gainputs().
	GB->entropy_convergence    = 1;
	GB->entropy_check_interval = GA_DEFAULT_ENTROPY_CHECK_INTERVAL;
	GB->entropy_window         = GA_DEFAULT_ENTROPY_WINDOW;
	GB->entropy_rel_threshold  = GA_DEFAULT_ENTROPY_REL_THRESHOLD;

	if (gainpfile[0] != '\0') {
		//GB->rrg_skip=0;
		GB->adaptive_ga=0;
		GB->num_print=GA_DEFAULT_NUM_PRINT;
		GB->print_int=GA_DEFAULT_PRINT_INT;
		GB->seed = GA_DEFAULT_SEED;

		printf("file in GA is <%s>\n",gainpfile);

		read_gainputs(FA,GB,&geninterval,&popszpartition,gainpfile);
	} else {
		printf("No GA input file — using pre-configured parameters\n");
	}

	// The search budget determines the result and, until this line, left no
	// witness on a SUCCESSFUL run: num_chrom appeared only in the
	// chrom_snapshot-overflow fprintf below, which fires solely when the value
	// is <= 0.  So "what population did this run use?" was answerable from the
	// source but never from the run's own output — and the three suppliers
	// disagree: the gate takes config_defaults.h's 1000, the campaign multiplies
	// by ceil(n_genes/4), and the legacy path reads NUMCHROM/NUMGENER from gainp.
	//
	// Placed AFTER the gainp block deliberately.  Printed before it, this line
	// would report the pre-file value on the legacy `./FlexAID cfg.inp ga.inp`
	// path — right for the two harnesses we benchmark and silently wrong for the
	// third, which is precisely the stale-witness failure this PR exists to
	// remove.  Here it reports what the GA actually runs with, on every path.
	printf("num_chrom=%d max_generations=%d\n",
	       GB->num_chrom, GB->max_generations);

	// Defensive clamp: a zero (or negative) check interval reaches three
	// integer division / modulo sites in the generation loop below and would
	// raise SIGFPE. Guarantee a sane value no matter how it was configured.
	if (GB->entropy_check_interval <= 0) {
		GB->entropy_check_interval = GA_DEFAULT_ENTROPY_CHECK_INTERVAL;
	}
	unsigned int tt;
	if (GB->seed == 0) {
		if (flexaids_rng::has_master_seed()) {
			tt = static_cast<unsigned int>(flexaids_rng::master_seed());
		} else {
			std::uint64_t env_seed = 0;
			if (flexaids_rng::env_seed(env_seed)) {
				tt = static_cast<unsigned int>(env_seed);
			} else {
				fprintf(stderr,
				        "ERROR: GA seed is 0 and FLEXAID_SEED is unset. "
				        "Refusing time(0) fallback (non-reproducible). "
				        "Set GB->seed or FLEXAID_SEED.\n");
				Terminate(1);
			}
		}
	} else {
		tt = GB->seed;
	}
	printf("srand=%u\n", tt);
	flexaids_rng::set_master_seed(static_cast<std::uint64_t>(tt));
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
		// Zero ring pucker side-channel (POD; malloc leaves it indeterminate).
		memset((*chrom)[i].ring_phases, 0, sizeof((*chrom)[i].ring_phases));
		memset((*chrom)[i].ring_six,    0, sizeof((*chrom)[i].ring_six));
		memset((*chrom)[i].ring_five,   0, sizeof((*chrom)[i].ring_five));
	}

	// *** chrom_snapshot
	// Use std::size_t for the product to prevent 32-bit signed overflow
	// (e.g., num_chrom=10000 × max_generations=100000 overflows int).
	const std::size_t snap_count = static_cast<std::size_t>(GB->num_chrom)
	                              * static_cast<std::size_t>(GB->max_generations);
	if (GB->num_chrom <= 0 || GB->max_generations <= 0 ||
	    snap_count > (SIZE_MAX / sizeof(chromosome)))
	{
		fprintf(stderr,"ERROR: chrom_snapshot size overflow (num_chrom=%d, max_generations=%d).\n",
		        GB->num_chrom, GB->max_generations);
		Terminate(2);
	}
	(*chrom_snapshot) = (chromosome*)malloc(snap_count * sizeof(chromosome));
	if(!(*chrom_snapshot))
	{
		fprintf(stderr,"ERROR: memory allocation error for chrom_snapshot (requested %zu bytes).\n",
		        snap_count * sizeof(chromosome));
		Terminate(2);
	}

	for(std::size_t snap_i = 0; snap_i < snap_count; ++snap_i)
	{
		(*chrom_snapshot)[snap_i].genes = (gene*)malloc(GB->num_genes*sizeof(gene));

		if(!(*chrom_snapshot)[snap_i].genes){
			fprintf(stderr,"ERROR: memory allocation error for chrom_snapshot[%zu].genes.\n",snap_i);
			Terminate(2);
		}

		(*chrom_snapshot)[snap_i].app_evalue = 0.0;
		(*chrom_snapshot)[snap_i].evalue = 0.0;
		(*chrom_snapshot)[snap_i].fitnes = 0.0;
		(*chrom_snapshot)[snap_i].boltzmann_weight = 0.0;
		(*chrom_snapshot)[snap_i].free_energy = 0.0;
		(*chrom_snapshot)[snap_i].status = ' ';
		memset((*chrom_snapshot)[snap_i].ring_phases, 0, sizeof((*chrom_snapshot)[snap_i].ring_phases));
		memset((*chrom_snapshot)[snap_i].ring_six,    0, sizeof((*chrom_snapshot)[snap_i].ring_six));
		memset((*chrom_snapshot)[snap_i].ring_five,   0, sizeof((*chrom_snapshot)[snap_i].ring_five));
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
	// G4.2: Cartesian ligand heavy-atom niche (env-OFF default). Gate lives in
	// niche_distance.h (unit-tested). Gene-space calc_rmsp mixes cleft-grid
	// ordinal (gene 0) with angular genes — PHASE4_GATES_ACTUALIZED defect.
	if (flexaids::niche_cartesian_env_enabled()) {
		GB->sig_share = flexaids::niche_cartesian_sigma_ang(2.0);
		fprintf(stderr,
		        "[NICHE-CART] enabled: sigma_share=%.4f A (ligand heavy-atom RMSD); "
		        "gene-space RMSP niche OFF\n",
		        GB->sig_share);
		printf("[NICHE-CART] SIGMA_SHARE=%.4f A (Cartesian ligand RMSD)\n", GB->sig_share);
	}
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
	// prev_best_fitness removed: SMFREE stagnation now tracks CF (evalue), not fit_max
	// For PSHARE fit_max is always num_chrom (rank 0); track best evalue instead.
	double prev_best_evalue = 1e30;    // PSHARE: best CF seen (lower=better)
	int    stagnation_count  = 0;
	bool   ga_stagnant = false;

	// ── [P5-ADAPTIVE-GEN] Adaptive-generation early convergence state (BEGIN) ──
	// Off by default. Opt in with FLEXAIDDS_ADAPTIVE_GENERATIONS=<K> (patience in
	// generations); the GA stops once best-CF improves by < eps for K consecutive
	// generations instead of always running max_generations. Independent of the
	// SMFREE/entropy termination paths above (no sec_may_terminate gate) so it
	// works for the classic CF-only fast-docking mode. Purely a wall-clock lever;
	// ranking/poses are unchanged when the flag is unset.
	int    ag_patience   = 0;       // 0 = disabled
	double ag_eps        = 1.0;     // kcal/mol improvement threshold (best-CF)
	double ag_prev_best  = 1e30;    // best CF observed at previous check
	int    ag_plateau    = 0;       // consecutive plateau generations
	{
		const char* ag_e = std::getenv("FLEXAIDDS_ADAPTIVE_GENERATIONS");
		if (ag_e && *ag_e) {
			ag_patience = std::atoi(ag_e);
			if (ag_patience < 0) ag_patience = 0;
		}
		const char* ag_eps_e = std::getenv("FLEXAIDDS_ADAPTIVE_EPS");
		if (ag_eps_e && *ag_eps_e) {
			const double v = std::atof(ag_eps_e);
			if (v > 0.0) ag_eps = v;
		}
		if (ag_patience > 0) {
			printf("[P5-ADAPTIVE-GEN] adaptive generations ON: patience=%d gens, "
			       "eps=%.4f kcal/mol (best-CF plateau early stop)\n",
			       ag_patience, ag_eps);
		}
	}
	// ── [P5-ADAPTIVE-GEN] state (END) ──

	// ── Always-on H plateau early exit: ring buffer over last 20 checks ─────
	// Fires independently of GB->entropy_convergence (that flag controls the
	// soft/hard/plateau checks above).  ε = 0.001 nats matches the thermal noise
	// floor of the Shannon estimate at convergence.
	constexpr int    kHPlateauWindow = 20;
	constexpr double kHPlateauEps    = 0.001;  // nats; ~0.00144 bits
	std::array<double, kHPlateauWindow> h_plateau_ring{};
	h_plateau_ring.fill(0.0);
	int  h_plateau_head   = 0;
	int  h_plateau_filled = 0;

	// ── P8: dual SEC termination gate (energy + gene-space joint convergence) ──
	// The Shannon-Entropy-Collapse (SEC) checks below detect collapse of the
	// *energy* histogram.  But an energy plateau can coexist with a population
	// that is still spread across multiple distinct binding-mode hypotheses in
	// *gene* space — terminating then throws away diversity that better cluster
	// selection (or further search) could still exploit.  When diversity
	// monitoring is enabled we therefore require JOINT collapse: an energy SEC
	// trigger only terminates if gene-space allele entropy has *also* collapsed
	// below diversity_collapse_threshold.  With monitoring off, behaviour is
	// unchanged (energy SEC alone terminates).
	auto sec_may_terminate = [&](int gen) -> bool {
		// SEC min-gen guard reverted (v29 ablation: early SEC collapse in oracle mode
		// is correct convergence behaviour, not premature death — gate was net harmful).
		if (!GB->diversity_monitoring) return true;  // gate disabled → legacy behaviour
		auto dm = ga_diversity::compute_diversity(
			*chrom, GB->num_chrom, GB->num_genes, *gene_lim,
			GB->diversity_collapse_threshold);
		if (!dm.collapse_detected) {
			printf("SEC fired at gen %d but gene-space still diverse "
			       "(allele_H=%.3f >= %.3f, min_gene_H=%.3f) — deferring termination\n",
			       gen, dm.allele_entropy, GB->diversity_collapse_threshold,
			       dm.min_gene_entropy);
		}
		return dm.collapse_detected;
	};

	// ── InStreamClustering: online medoid clustering during GA ──
	flexaids::InStreamCluster instream_cluster(
	    GA_INSTREAM_RMSD_THRESHOLD,
	    GA_INSTREAM_MAX_MEDOIDS,
	    GB->num_genes);

	// Merge/H(ω) cadence: defaults to GA_INSTREAM_INTERVAL.  FLEXAIDDS_INSTREAM_INTERVAL
	// overrides it (clamped >= 1) so short diagnostic runs can populate the medoid
	// set — and emit the Level-3 [HVIB] monitor — before the GA's entropy-collapse
	// early exit (which can trigger well before generation 100).  Default behaviour
	// for production benchmarks is unchanged.
	int instream_interval = GA_INSTREAM_INTERVAL;
	if (proto.instream_interval >= 1) {
		instream_interval = proto.instream_interval;
		fprintf(stderr, "[INSTREAM] merge/H(ω) cadence overridden to every %d "
		        "generation(s) via FLEXAIDDS_INSTREAM_INTERVAL\n", instream_interval);
	}

	////// Genetic Algorithm ///////
	////////////////////////////////

	// FLEXAIDDS_NO_SEC (presence) disables early-exit paths (stagnation +
	// entropy/SEC) so the *full* generation budget is always used.
	//
	// FLEXAIDDS_BENCHMARK is NOT equivalent. The two env vars set two distinct
	// fields (ProtocolConfig.cpp: cfg.no_sec / cfg.benchmark_mode):
	//
	//   FLEXAIDDS_NO_SEC     -> disables BOTH entropy/SEC exits AND stagnation
	//   FLEXAIDDS_BENCHMARK  -> disables the STAGNATION plateau exit only
	//                           (`benchmark_full`, below); the entropy/SEC
	//                           guards test `no_sec` alone and ignore it
	//
	// This comment previously read "(or FLEXAIDDS_BENCHMARK=1)", generalising
	// the stagnation behaviour to both exits. That is the half that matters:
	// all 80 restarts of the Jul-28 arm-A pilot terminated by ENTROPY
	// convergence (stagnation fired 0/80) at a median of 165 of 2000
	// generations -- ~8% of the intended budget -- and the results were
	// recorded as full-budget.
	//
	// If you need the full budget, set FLEXAIDDS_NO_SEC and verify the
	// "[SEC] ... DISABLED" line below appears on STDERR (not stdout).
	// During benchmarking this ensures equal search effort vs other methods.
	// "Spare" generations after a plateau are used with boosted mutation/exploration
	// (see stagnation handling) to search conformational space more effectively.
	const bool no_sec = proto.no_sec;
	if (no_sec)
		fprintf(stderr, "[SEC] All entropy-convergence early exits DISABLED "
		        "(FLEXAIDDS_NO_SEC=1) — GA will run to max_generations.\n");
	const unsigned int target_temperature_K = FA->temperature;

	// ── True GA elitism (v27) snapshot buffers ──
	// The n_elite lowest-CF individuals are deep-copied each generation BEFORE
	// boom injection / reproduce()/sharing, then restored over the worst of the
	// freshly reproduced population so the running best can never be ejected by
	// diversity pressure (boom injection or niche-sharing fitness reduction).
	const int n_elite = (GB->n_elite > 0)
	                    ? std::min(GB->n_elite, GB->num_chrom)
	                    : 0;
	std::vector<gene>   elite_genes_buf(static_cast<size_t>(n_elite) * GB->num_genes);
	std::vector<cfstr>  elite_cf_buf(n_elite);
	std::vector<double> elite_eval_buf(n_elite), elite_app_buf(n_elite);
	if (n_elite > 0)
		fprintf(stderr, "[ELITE] GA-internal elitism active: protecting %d "
		        "lowest-CF individual(s) per generation\n", n_elite);

	// ── Temperature annealing (FLEXAIDDS_T_HOT via ProtocolConfig) ────────────
	// Exponential decay T_hot -> target K:
	//   T(alpha) = T_hot*exp(-5alpha) + T_target*(1-exp(-5alpha)).
	// Affects SMFREE Boltzmann-weight selection only; post-GA thermodynamics
	// use the configured target temperature even when the GA stops early.
	// arm3b ablation (5000K constant, Fable 5) was net-neutral in oracle mode
	// → native basin gravitationally dominant.  Annealing targets near-miss
	// false-minimum escape early in the run while native seeds lock in.
	// Useful calibration range: 500–2000 K.
	const double t_hot_anneal = proto.t_hot;
	const double target_temperature_d = static_cast<double>(target_temperature_K);
	const bool do_anneal = (target_temperature_K > 0) &&
	                        (t_hot_anneal > target_temperature_d);
	if (do_anneal) {
		fprintf(stderr, "[ANNEAL] Temperature annealing enabled: "
		        "T_hot=%.0f K -> %.0f K over %d generations (exp-5 schedule)\n",
		        t_hot_anneal, target_temperature_d, GB->max_generations);
		// Prime initial temperature so gen-0 SMFREE sees T_hot
		FA->temperature = static_cast<unsigned int>(std::round(t_hot_anneal));
		FA->beta        = 1.0 / t_hot_anneal;
	}

	// ── ThermodynamicEngine: unbound reference H(ω) (before GA, atoms[] = initial state) ──
	if (FA->thermo_engine_enabled && FA->thermo_engine != nullptr) {
		// (a) Receptor backbone torsional H(ω) — CA positions are fixed throughout docking
		if (FA->is_protein && FA->res_cnt > GA_TENCM_MIN_RESIDUES) {
			tencm::TorsionalENM rec_enm;
			rec_enm.build(atoms, residue, FA->res_cnt);
			if (rec_enm.is_built()) {
				std::vector<double> rec_eigs;
				rec_eigs.reserve(rec_enm.modes().size());
				for (const auto& nm : rec_enm.modes())
					if (nm.eigenvalue > 0.0) rec_eigs.push_back(nm.eigenvalue);
				if (!rec_eigs.empty()) {
					const std::vector<std::vector<double>> single = { rec_eigs };
					FA->H_rep_receptor_ref = static_cast<float>(
						vibentropy::compute_vib_entropy_collapse(single).H_pop);
				}
			}
		}
		// (b) Free-ligand H(ω) — initial (unbound) conformation, before GA modifies atoms[]
		{
			const int lig_start     = (FA->resligand && FA->resligand->fatm)
			                          ? FA->resligand->fatm[0] : -1;
			const int lig_end_incl  = (FA->resligand && FA->resligand->latm)
			                          ? FA->resligand->latm[0] : -1;
			if (lig_start >= 0 && lig_end_incl >= lig_start) {
				tencm::TorsionalENM lig_enm_free;
				lig_enm_free.build_from_ligand(atoms, lig_start, lig_end_incl + 1);
				if (lig_enm_free.is_built()) {
					std::vector<double> eigs;
					eigs.reserve(lig_enm_free.modes().size());
					for (const auto& nm : lig_enm_free.modes())
						if (nm.eigenvalue > 0.0) eigs.push_back(nm.eigenvalue);
					if (!eigs.empty()) {
						const std::vector<std::vector<double>> single = { eigs };
						FA->H_rep_ligand_ref = static_cast<float>(
							vibentropy::compute_vib_entropy_collapse(single).H_pop);
					}
				}
			}
		}
		FA->thermo_engine->set_unbound_reference(FA->H_rep_receptor_ref, FA->H_rep_ligand_ref);
	}

	// ── Per-generation timing (bench) ──
	double _sum_gen_ms = 0.0;
	int    _n_gen_timed = 0;
	for(i=0;i<GB->max_generations;i++)
	{
		///////////////////////////////////////////////////

		state=check_state(PAUSEFILE,ABORTFILE,STOPFILE,INTERVAL);

		if(state == -1){
			return(state);
		}else if(state == 1){
			break;
		}

		auto _t0_gen = std::chrono::steady_clock::now();

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

			// Recompute MIF for adapted grid.
			// Allocate new buffers BEFORE freeing old ones so OOM leaves the
			// previous MIF intact rather than crashing on null deref.
			if (FA->mif_enabled || FA->grid_prio_percent < 100.0f) {
				std::vector<atom> protein_atoms(atoms, atoms + FA->atm_cnt_real);
				cavity_detect::SpatialGrid sg;
				sg.build(protein_atoms);
				auto mif = mif::compute_mif(*cleftgrid, FA->num_grd,
				                             atoms, FA->atm_cnt_real, sg);
				mif::build_sampling_cdf(mif, FA->mif_temperature);

				const std::size_t n_energies = mif.energies.size();
				const std::size_t n_sorted   = mif.sorted_indices.size();
				const std::size_t n_cdf      = mif.cdf.size();
				float*  new_energies = static_cast<float*>(malloc(n_energies * sizeof(float)));
				int*    new_sorted   = static_cast<int*>(malloc(n_sorted * sizeof(int)));
				double* new_cdf      = static_cast<double*>(malloc(n_cdf * sizeof(double)));

				if (!new_energies || !new_sorted || !new_cdf) {
					fprintf(stderr,
					        "ERROR: MIF allocation failed at generation %d "
					        "(energies=%zu sorted=%zu cdf=%zu) — keeping old MIF.\n",
					        i+1, n_energies, n_sorted, n_cdf);
					free(new_energies); free(new_sorted); free(new_cdf);
				} else {
					std::copy_n(mif.energies.data(), n_energies, new_energies);
					std::copy_n(mif.sorted_indices.data(), n_sorted, new_sorted);
					std::copy_n(mif.cdf.data(), n_cdf, new_cdf);

					free(FA->mif_energies); free(FA->mif_sorted); free(FA->mif_cdf);
					FA->mif_count    = static_cast<int>(n_sorted);
					FA->mif_energies = new_energies;
					FA->mif_sorted   = new_sorted;
					FA->mif_cdf      = new_cdf;
				}
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

		// Stagnation detection: track best CF (evalue), not SMFREE fitness.
		// SMFREE fitness overflows to cap=1000 for all chromosomes immediately
		// (exp(-CF/kT) underflows for typical CF -50 to -200 kcal/mol at kT=0.596),
		// so tracking fit_max was a bug: stagnation fired at gen ~100-300 while
		// gene-space allele_H was still 43-86% — search murdered by saturated proxy.
		// Fable-5 analysis: preserve _ps discriminator for _tol (PSHARE converges
		// in sub-kcal increments during late descent; 1e-3 sensitivity is correct there).
		// Joint termination gate: stagnation only terminates if BOTH CF is stuck AND
		// gene space has collapsed (sec_may_terminate checks allele_H < 0.300).
		// For SMFREE+N=1000: allele_H drifts to 0.300 in ~2100 gens (mutation-drift
		// equilibrium), so 1GPK-class targets (H=0.863) run to gen 2000 — correct.
		if ((i + 1) % STAGNATION_WINDOW == 0 && i > 0) {
			const bool _ps = (strcmp(GB->fitness_model,"PSHARE")==0);
			const double _cur  = (*chrom)[0].evalue;  // always CF; elitism guarantees monotonic
			const double _prev = prev_best_evalue;
			const double _tol  = _ps ? 1e-3 : 1.0;   // SMFREE: 1 kcal/mol (~1.68 kT) threshold
			const bool stagnant = (_prev - _cur) < _tol;  // improvement = _prev - _cur
			if (stagnant) {
				stagnation_count += STAGNATION_WINDOW;
				if (stagnation_count >= STAGNATION_LIMIT) {
					const bool benchmark_full = proto.no_sec || proto.benchmark_mode;
					if (benchmark_full) {
						printf("GA plateau at gen %d (stagnant %d gens, best_CF=%.4f); "
						       "BENCHMARK mode: continuing with exploration boost for remaining gens.\n",
						       i+1, stagnation_count, _cur);
						GB->mut_rate = std::min(0.25, GB->mut_rate * 3.0);
						stagnation_count = 0;
					} else if (sec_may_terminate(i + 1)) {
						// Joint condition: CF stagnant AND gene space collapsed.
						// If allele_H is still high (e.g. 0.863 for 1GPK-class targets
						// under SMFREE+drift), sec_may_terminate defers -> keep running.
						printf("GA terminated: CF stagnant for %d gens (best_CF=%.4f) "
						       "with gene-space collapsed\n", stagnation_count, _cur);
						ga_stagnant = true;
						break;
					}
					// CF stagnant but gene space diverse -> continue search
				}
			} else {
				stagnation_count = 0;
			}
			prev_best_evalue = _cur;  // always update (was conditional on _ps)
		}

		// ── [P5-ADAPTIVE-GEN] Best-CF plateau early stop (BEGIN) ────────────
		// Localized generation-loop convergence check. Guarded by ag_patience>0
		// (FLEXAIDDS_ADAPTIVE_GENERATIONS); when disabled this is a single
		// predictable branch and behavior is identical to legacy. Uses the
		// elitism-guaranteed best chromosome (*chrom)[0].evalue (CF, lower=better).
		if (ag_patience > 0 && i > 0) {
			const double ag_cur = (*chrom)[0].evalue;
			if ((ag_prev_best - ag_cur) < ag_eps) {
				if (++ag_plateau >= ag_patience) {
					printf("[P5-ADAPTIVE-GEN] GA converged: best-CF plateau for %d "
					       "gens (best_CF=%.4f) at gen %d — early stop "
					       "(max_generations=%d)\n",
					       ag_plateau, ag_cur, i + 1, GB->max_generations);
					break;
				}
			} else {
				ag_plateau = 0;
			}
			ag_prev_best = ag_cur;
		}
		// ── [P5-ADAPTIVE-GEN] Best-CF plateau early stop (END) ──────────────

		// ── Always-on H plateau early exit ─────────────────────────────────
		// Every entropy_check_interval generations, sample H of the current
		// population and push into a 20-slot ring buffer.  If the absolute
		// difference between the newest and oldest slot < kHPlateauEps nats,
		// the distribution has stopped collapsing → write best pose and stop.
		if (!no_sec && !entropy_converged && !ga_stagnant &&
		    ((i + 1) % GB->entropy_check_interval == 0)) {
			std::vector<double> _hp_energies(GB->num_chrom);
			for (int _c = 0; _c < GB->num_chrom; ++_c)
				_hp_energies[_c] = (*chrom)[_c].evalue;
			const double H_now = shannon_thermo::compute_shannon_entropy(
				_hp_energies, shannon_thermo::DEFAULT_HIST_BINS);
			h_plateau_ring[h_plateau_head] = H_now;
			h_plateau_head = (h_plateau_head + 1) % kHPlateauWindow;
			if (h_plateau_filled < kHPlateauWindow) ++h_plateau_filled;
			if (h_plateau_filled == kHPlateauWindow) {
				// oldest entry is now at h_plateau_head (ring has wrapped)
				const double delta = std::abs(H_now - h_plateau_ring[h_plateau_head]);
				if (delta < kHPlateauEps && sec_may_terminate(i + 1)) {
					printf("Early exit at gen %d: H plateau < %.4f "
					       "(H_now=%.6f nats, delta=%.6f nats)\n",
					       i + 1, kHPlateauEps, H_now, delta);
					entropy_converged = true;
					break;
				}
			}
		}

		// Entropy convergence check (opt-in via ENTRCNVG config keyword)
		if (!no_sec && GB->entropy_convergence &&
		    ((i + 1) % GB->entropy_check_interval == 0)) {
			std::vector<double> pop_energies(GB->num_chrom);
			for (int c = 0; c < GB->num_chrom; ++c) {
				pop_energies[c] = (*chrom)[c].evalue;
			}
			double H = shannon_thermo::compute_shannon_entropy(
				pop_energies, shannon_thermo::DEFAULT_HIST_BINS);
			entropy_history.push_back(H);

			if (H <= shannon_thermo::kHSC_soft_nats && sec_may_terminate(i + 1)) {
				printf("Entropy collapse convergence at generation %d "
				       "(H=%.4f nats <= %.4f nats / %.1f bits)\n",
				       i + 1, H,
				       shannon_thermo::kHSC_soft_nats,
				       shannon_thermo::kHSC_soft_bits);
				entropy_converged = true;
				break;
			}

			if (shannon_thermo::detect_entropy_plateau(
			        entropy_history, GB->entropy_window,
			        GB->entropy_rel_threshold) && sec_may_terminate(i + 1)) {
				printf("Entropy convergence at generation %d "
				       "(H=%.4f nats, stable for %d checks)\n",
				       i + 1, H, GB->entropy_window);
				entropy_converged = true;
				break;
			}

			// Hard-zone variance plateau: H < ln2 AND σ² < (0.005 nats)²
			// Once in the collapsed zone, tiny variance confirms convergence.
			if (H < shannon_thermo::kHSC_hard_nats &&
			    static_cast<int>(entropy_history.size()) >= 20) {
				constexpr int    kHardWindow = 20;
				constexpr double kPlateauVarThresh = 0.005 * 0.005; // nats²
				const size_t h_start = entropy_history.size() - kHardWindow;
				double h_sum = 0.0, h_sum2 = 0.0;
				for (size_t k = h_start; k < entropy_history.size(); ++k) {
					h_sum  += entropy_history[k];
					h_sum2 += entropy_history[k] * entropy_history[k];
				}
				const double h_mean = h_sum / kHardWindow;
				const double h_var  = h_sum2 / kHardWindow - h_mean * h_mean;
				if (h_var < kPlateauVarThresh && sec_may_terminate(i + 1)) {
					printf("Shannon entropy converged at generation %d "
					       "(H=%.4f < %.4f nats, var=%.6f nats²) — early stop\n",
					       i + 1, H, shannon_thermo::kHSC_hard_nats, h_var);
					entropy_converged = true;
					break;
				}
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
				if (flexaids::new_search::basin_reinject_enabled()) {
					// B: basin-aware reinject — re-randomize worst fraction until
					// Cartesian ligand RMSD vs current best exceeds sigma (default 2 Å).
					static bool s_basin_logged = false;
					if (!s_basin_logged) {
						s_basin_logged = true;
						std::fprintf(stderr,
						             "[NEW-SEARCH-ARCH] basin_reinject=1: catastrophic "
						             "reinject prefers Cartesian RMSD > sigma vs best "
						             "(FLEXAIDDS_BASIN_SIGMA_ANG, default 2.0)\n");
					}
					const double sigma = flexaids::new_search::basin_sigma_ang(2.0);
					const int n_chrom = GB->num_chrom;
					const int n_genes = GB->num_genes;
					int n_mutate = static_cast<int>(
					    std::ceil(GB->catastrophic_mutation_fraction * n_chrom));
					if (n_mutate > n_chrom) n_mutate = n_chrom;
					// Sort indices by evalue ascending (best first).
					std::vector<int> order(static_cast<size_t>(n_chrom));
					for (int q = 0; q < n_chrom; ++q) order[static_cast<size_t>(q)] = q;
					std::sort(order.begin(), order.end(), [&](int a, int b) {
						return (*chrom)[a].evalue < (*chrom)[b].evalue;
					});
					const int best_i = order[0];
					// Coords for current best.
					constexpr int kCoordStride = MAX_ATM_HET * 3;
					std::vector<float> best_xyz(static_cast<size_t>(kCoordStride), 0.0f);
					std::vector<float> trial_xyz(static_cast<size_t>(kCoordStride), 0.0f);
					calc_rmsd_chrom(FA, GB, *chrom, *gene_lim, atoms, residue, *cleftgrid,
					                n_genes, best_i, best_i, best_xyz.data(), nullptr, false);
					const int lres = atoms[FA->map_par[0].atm].ofres;
					const int rot = residue[lres].rot;
					int n_lig = residue[lres].latm[rot] - residue[lres].fatm[rot] + 1;
					if (n_lig < 1) n_lig = 1;
					int n_accepted = 0;
					for (int k = 0; k < n_mutate; ++k) {
						const int ci = order[static_cast<size_t>(n_chrom - 1 - k)];
						bool ok = false;
						for (int attempt = 0; attempt < 24; ++attempt) {
							generate_random_individual(FA, GB, atoms, (*chrom)[ci].genes,
							                           *gene_lim, dice, 0, n_genes);
							for (int g = 0; g < n_genes; ++g)
								(*chrom)[ci].genes[g].to_ic =
								    genetoic(&(*gene_lim)[g], (*chrom)[ci].genes[g].to_int32);
							calc_rmsd_chrom(FA, GB, *chrom, *gene_lim, atoms, residue,
							                *cleftgrid, n_genes, ci, ci, trial_xyz.data(),
							                nullptr, false);
							const double rmsd = flexaids::niche_cartesian_rmsd(
							    best_xyz.data(), trial_xyz.data(), n_lig);
							if (flexaids::new_search::outside_basin(rmsd, sigma) ||
							    attempt == 23) {
								ok = flexaids::new_search::outside_basin(rmsd, sigma);
								break;
							}
						}
						if (ok) ++n_accepted;
						ring_randomise_chrom(FA, &(*chrom)[ci]);
						ring_load_chrom_to_fa(FA, &(*chrom)[ci]);
						(*chrom)[ci].cf = eval_chromosome(FA, GB, VC, *gene_lim, atoms,
						                                  residue, *cleftgrid,
						                                  (*chrom)[ci].genes, target);
						(*chrom)[ci].evalue =
						    get_cf_evalue(&(*chrom)[ci].cf, FA) / n_receptor_chains;
						(*chrom)[ci].app_evalue =
						    get_apparent_cf_evalue(&(*chrom)[ci].cf) / n_receptor_chains;
						ccbm_inject_strain(FA, (*chrom)[ci], *gene_lim);
						(*chrom)[ci].status = 'n';
					}
					GB->catastrophic_mutation_count++;
					fprintf(stderr,
					        "[BASIN-REINJECT] #%d gen %d n=%d outside_basin=%d "
					        "sigma=%.2fA (H_allele=%.3f min_gene=%.3f)\n",
					        GB->catastrophic_mutation_count, i + 1, n_mutate, n_accepted,
					        sigma, dm.allele_entropy, dm.min_gene_entropy);
				} else {
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
		}

		// ── True GA elitism (v27): snapshot the n_elite lowest-CF individuals ──
		// Taken every generation BEFORE boom injection and reproduce()/sharing so
		// the global best survives both.  Restored over the worst of the new
		// population right after reproduce() below.  evalue is the (non-apparent)
		// CF; lower = better pose.
		if (n_elite > 0) {
			// Simple CF-minimum elite selection (v31+).
			// ε-tiebreaker reverted: IC-space distance to opt_par is unreliable for
			// already-converged poses (IC≠Cartesian proximity; nonlinear FK map +
			// angular degeneracy). v30 ablation: net-neutral with 6 regressions on
			// previously-perfect targets (0.00→2-6Å). Revert to plain CF sort.
			std::vector<int> eidx(GB->num_chrom);
			for (int q = 0; q < GB->num_chrom; ++q) eidx[q] = q;
			std::partial_sort(eidx.begin(), eidx.begin() + n_elite, eidx.end(),
				[&](int a, int b){
					return (*chrom)[a].evalue < (*chrom)[b].evalue;
				});
			for (int e = 0; e < n_elite; ++e) {
				const chromosome& src = (*chrom)[eidx[e]];
				elite_cf_buf[e]   = src.cf;
				elite_eval_buf[e] = src.evalue;
				elite_app_buf[e]  = src.app_evalue;
				for (int g = 0; g < GB->num_genes; ++g)
					elite_genes_buf[static_cast<size_t>(e) * GB->num_genes + g] = src.genes[g];
			}
		}

		// ── P5: periodic BOOM random injection (diversity insurance) ──
		// Every boom_inject_interval generations, replace the worst
		// (boom_inject_fraction × num_chrom/2) chromosomes with FRESH random
		// individuals (not seeds), always preserving the better half.  Unlike the
		// collapse-triggered catastrophic mutation above (which mutates existing
		// genes only in the first half of the run), this fires unconditionally on
		// a fixed cadence and resets genes to pure random — forcing the population
		// to keep exploring alternative binding-mode basins instead of all
		// converging onto the single deepest (often false) VCT minimum.  The base
		// population [0,num_chrom) already carries valid CF/evalue from the prior
		// generation, so "worst" is well-defined; fresh randoms are scored here so
		// the upcoming reproduce()/fitness pass sees correct energies.
		if (GB->boom_inject_interval > 0 && GB->boom_inject_fraction > 0.0 &&
		    ((i + 1) % GB->boom_inject_interval == 0) && (i + 1) < GB->max_generations) {
			const int half = GB->num_chrom / 2;
			int n_inject = (int)(GB->boom_inject_fraction * (double)half);
			if (n_inject > half) n_inject = half;
			if (n_inject > 0) {
				// Worst n_inject chromosomes by evalue (higher evalue = worse pose).
				std::vector<int> bidx(GB->num_chrom);
				for (int q = 0; q < GB->num_chrom; ++q) bidx[q] = q;
				std::partial_sort(bidx.begin(), bidx.begin() + n_inject, bidx.end(),
					[&](int a, int b){ return (*chrom)[a].evalue > (*chrom)[b].evalue; });
				for (int q = 0; q < n_inject; ++q) {
					int ci = bidx[q];
					generate_random_individual(FA, GB, atoms, (*chrom)[ci].genes,
					                           *gene_lim, dice, 0, GB->num_genes);
					ring_randomise_chrom(FA, &(*chrom)[ci]);
					ring_load_chrom_to_fa(FA, &(*chrom)[ci]);
					(*chrom)[ci].cf = eval_chromosome(FA, GB, VC, *gene_lim, atoms,
					                                  residue, *cleftgrid,
					                                  (*chrom)[ci].genes, target);
					(*chrom)[ci].evalue     = get_cf_evalue(&(*chrom)[ci].cf, FA) / n_receptor_chains;
					(*chrom)[ci].app_evalue = get_apparent_cf_evalue(&(*chrom)[ci].cf) / n_receptor_chains;
					ccbm_inject_strain(FA, (*chrom)[ci], *gene_lim);
					(*chrom)[ci].status = 'n';
				}
				GB->boom_inject_count++;
				fprintf(stderr, "[BOOM] injection #%d at gen %d: re-randomized worst "
				        "%d/%d chromosomes (fresh random, better half preserved)\n",
				        GB->boom_inject_count, i + 1, n_inject, GB->num_chrom);
			}
		}

		// ── Per-generation temperature update for SMFREE annealing ───────────
		// Must fire before reproduce() so the SMFREE Boltzmann-weight selection
		// inside reproduce() sees the correct annealed temperature.
		if (do_anneal && GB->max_generations > 1) {
			const double alpha = static_cast<double>(i) /
			                     static_cast<double>(GB->max_generations - 1);
			const double T_now = t_hot_anneal * std::exp(-5.0 * alpha)
			                   + target_temperature_d * (1.0 - std::exp(-5.0 * alpha));
			FA->temperature = static_cast<unsigned int>(std::round(T_now));
			FA->beta        = 1.0 / T_now;
			if (i % 200 == 0) {
				fprintf(stderr, "[ANNEAL] gen=%4d  T=%7.1f K  α=%.4f\n",
				        i + 1, T_now, alpha);
			}
		}

		// G4.3: FLEXAIDDS_MUTATION_GRANULAR uses phenotype-live ±1-bin steps.
		// hash_genes() keys the lifetime `duplicates` map on to_ic (phenotype).
		// Local steps exhaust the neighborhood; reproduce() then stalls in
		// rejection sampling (~0.08 gen/s observed on 1L7F after gen ~4100).
		// Clear once per generation under granular mode so uniqueness is
		// within-gen only. Classic (env off) keeps lifetime uniqueness.
		{
			const char* gran_e = std::getenv("FLEXAIDDS_MUTATION_GRANULAR");
			if (gran_e && (gran_e[0] == '1' || gran_e[0] == 'y' || gran_e[0] == 'Y' ||
			               gran_e[0] == 't' || gran_e[0] == 'T')) {
				duplicates.clear();
				static bool s_dup_clear_logged = false;
				if (!s_dup_clear_logged) {
					s_dup_clear_logged = true;
					std::fprintf(stderr,
					             "[MUT-GRAN] per-generation phenotype-duplicate clear "
					             "(lifetime map incompatible with ±1-bin local search)\n");
				}
			}
		}

		nrejected = reproduce(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),
				      GB->rep_model,GB->mut_rate,GB->cross_rate,print,dice,duplicates,target,*ctx);

		// ── True GA elitism (v27): restore snapshotted elites ──
		// reproduce() has rebuilt the population (selection/crossover/mutation
		// over sharing-reduced fitness).  Overwrite the n_elite WORST individuals
		// (highest evalue) of the new generation with the elites captured before
		// boom/sharing, guaranteeing the running best is carried forward intact.
		if (n_elite > 0) {
			std::vector<int> widx(GB->num_chrom);
			for (int q = 0; q < GB->num_chrom; ++q) widx[q] = q;
			std::partial_sort(widx.begin(), widx.begin() + n_elite, widx.end(),
				[&](int a, int b){ return (*chrom)[a].evalue > (*chrom)[b].evalue; });
			for (int e = 0; e < n_elite; ++e) {
				chromosome& dst = (*chrom)[widx[e]];
				dst.cf              = elite_cf_buf[e];
				dst.evalue          = elite_eval_buf[e];
				dst.app_evalue      = elite_app_buf[e];
				dst.fitnes          = 0.0;   // recomputed next reproduce()
				dst.boltzmann_weight = 0.0;
				dst.free_energy     = 0.0;
				dst.status          = 'n';   // CF is valid (deep-copied) — no re-eval
				for (int g = 0; g < GB->num_genes; ++g)
					dst.genes[g] = elite_genes_buf[static_cast<size_t>(e) * GB->num_genes + g];
			}
		}

		// Fix 6: write snapshots COMPACTLY (stride = save_num_chrom), so the
		// writer's layout matches what the post-GA thermo reader consumes.
		// Previously the writer strided by GB->num_chrom while n_chrom_snapshot
		// (and every downstream reader: QuickSort, StatMechEngine, FastOPTICS,
		// TQENS) advanced by save_num_chrom and read the compact prefix
		// [0, n_chrom_snapshot). With save_num_chrom = 5% of num_chrom, only
		// generation 0's records landed in that prefix; every later generation
		// was written far past it, leaving the prefix filled with zero-init
		// chromosomes → predicted_dH = predicted_TdS = 0 in every run. Compact
		// stride packs gen i at [i*save_num_chrom, (i+1)*save_num_chrom), exactly
		// the range the reader scans. The allocation (num_chrom*max_generations)
		// already dwarfs the compact extent (save_num_chrom*max_generations).
		save_snapshot(&(*chrom_snapshot)[i*save_num_chrom],(*chrom),save_num_chrom,GB->num_genes);
		n_chrom_snapshot += save_num_chrom;


		// ── InStreamClustering: merge top-K elites every N generations ──
		if (((i + 1) % instream_interval == 0) && save_num_chrom > 0) {
			const int top_k = std::min(GA_INSTREAM_TOP_K, save_num_chrom);
			std::vector<float> elite_genes(static_cast<size_t>(top_k) * GB->num_genes);
			std::vector<double> elite_scores(top_k);
			for (int ek = 0; ek < top_k; ++ek) {
				for (int g = 0; g < GB->num_genes; ++g) {
					elite_genes[static_cast<size_t>(ek) * GB->num_genes + g] =
						static_cast<float>((*chrom)[ek].genes[g].to_ic);
				}
				elite_scores[ek] = (*chrom)[ek].app_evalue;
			}
			instream_cluster.merge_elites(
				elite_genes.data(), elite_scores.data(),
				top_k, i + 1, GB->num_genes);

			// ── Level-3 H(ω): ligand vibrational-mode Shannon entropy ──────────
			// For each current cluster representative, materialise its ligand
			// conformation (eval_chromosome → ic2cf writes Cartesian coords into
			// the shared atoms[] buffer), build a Cartesian ANM over the ligand
			// heavy atoms, and harvest the mode eigenvalues.  Each rep has a
			// different ligand conformation → different stiffness matrix →
			// different spectrum, so the pooled distribution measures vibrational
			// (ω-space) collapse across the GA population.  Purely diagnostic —
			// eigenvalues never enter CF or fitness.  Gated on use_shannon, which
			// defaults OFF so existing benchmarks are unaffected.
			if (GB->use_shannon) {
				const int lig_start = (FA->resligand && FA->resligand->fatm)
				                      ? FA->resligand->fatm[0] : -1;
				const int lig_end_incl = (FA->resligand && FA->resligand->latm)
				                      ? FA->resligand->latm[0] : -1;
				if (lig_start >= 0 && lig_end_incl >= lig_start) {
					const auto& reps = instream_cluster.snapshot();
					std::vector<std::vector<double>> rep_eigs;
					rep_eigs.reserve(reps.size());
					std::vector<gene>     tmp_genes(GB->num_genes);
					tencm::TorsionalENM   lig_enm;
					for (const auto& med : reps) {
						if (static_cast<int>(med.genes_ic.size()) < GB->num_genes)
							continue;
						for (int g = 0; g < GB->num_genes; ++g)
							tmp_genes[g].to_ic = med.genes_ic[g];
						// Materialise coords into shared atoms[] (intentionally
						// overwrites it — the GA carries genes, not coords, so the
						// next generation re-materialises from genes anyway).
						eval_chromosome(FA, GB, VC, *gene_lim, atoms, residue,
						                *cleftgrid, tmp_genes.data(), target);
						// Half-open [lig_start, lig_end); latm[0] is inclusive.
						lig_enm.build_from_ligand(atoms, lig_start, lig_end_incl + 1);
						if (!lig_enm.is_built()) continue;
						std::vector<double> eigs;
						eigs.reserve(lig_enm.modes().size());
						for (const auto& nm : lig_enm.modes())
							eigs.push_back(nm.eigenvalue);
						rep_eigs.push_back(std::move(eigs));
					}
					if (!rep_eigs.empty()) {
						vibentropy::VibEntropyResult vr =
							vibentropy::compute_vib_entropy_collapse(rep_eigs);
						fprintf(stderr,
						        "[HVIB] gen=%d H_pop=%.6f H_rep_mean=%.6f "
						        "D_vib=%.6f n_reps=%d\n",
						        i + 1, vr.H_pop, vr.H_rep_mean, vr.D_vib,
						        vr.n_reps);
					}
				}
			}
		}

		if(strcmp(GB->fitness_model,"PSHARE")==0){
			QuickSort((*chrom),0,GB->num_chrom-1,false);

			if(print){
				printf("best by fitnes\n");
				print_par((*chrom),(*gene_lim),GB->num_print,GB->num_genes, stdout);
			}
		}

		// ── Mid-run H_shannon snapshots for entropy-collapse causality test ──
		// Fires at gen 500 and gen 1000 (if reached) when thermo is enabled.
		// Uses the same 256-bin gene-space histogram as the post-GA
		// ThermodynamicEngine::shannon_entropy() call: normalize each gene via
		// to_int32 / MAX_RANDOM_VALUE → [0,1], bin to [0,255], compute Shannon H,
		// then scale by T_eff to match the TdS_shannon units in [THERMO].
		if (FA->thermo_engine_enabled && FA->thermo_engine != nullptr) {
			const int cur_gen = i + 1;
			if (cur_gen == 500 || cur_gen == 1000) {
				float H_snap = 0.0f;
				std::array<int, 256> hist{};
				for (int g = 0; g < GB->num_genes; ++g) {
					hist.fill(0);
					for (int c = 0; c < GB->num_chrom; ++c) {
						const int bin = std::clamp(
							static_cast<int>(
								static_cast<float>((*chrom)[c].genes[g].to_int32)
								/ static_cast<float>(MAX_RANDOM_VALUE) * 255.0f),
							0, 255);
						++hist[bin];
					}
					const float n_inv = 1.0f / static_cast<float>(GB->num_chrom);
					for (int b = 0; b < 256; ++b) {
						if (hist[b] > 0) {
							const float p = static_cast<float>(hist[b]) * n_inv;
							H_snap -= p * std::log2(p);
						}
					}
				}
				const float TdS_snap = FA->thermo_T_eff * H_snap;
				printf("[THERMO_SNAP gen=%d] TdS_shannon=%.6f\n", cur_gen, TdS_snap);
			}
		}

		// ── IP-6: per-generation any-pose RMSD trace (.gentrace.tsv) ────────────
		// WHY THIS EXISTS
		//   The terminal .pop.tsv / .rrd dumps are a snapshot of the FINAL
		//   population only. A restart that samples a near-native pose at
		//   generation 12 and loses it to selection pressure by generation 50 is
		//   indistinguishable, in those files, from a restart that never sampled
		//   one at all. On the 2026-08-06 Astex-85 dock-once campaign that
		//   ambiguity covers 72 of 84 targets: they are reported "search-limited"
		//   when the true failure may be convergence/retention inside the GA.
		//   1R1H is the existence proof -- its 5 restarts have terminal population
		//   minima 11.23 / 9.86 / 1.78 / 10.24 / 10.79 A, so the basin is reachable
		//   but usually lost.
		//
		// WHAT IT RECORDS
		//   Every generation, the MINIMUM RMSD over the ENTIRE live population,
		//   plus the CF of that pose and the CF of the population's best-scoring
		//   member. Deliberately unfiltered: selecting by CF first would hide
		//   exactly the case under test (a near-native pose the objective ranks
		//   badly). best_rmsd_so_far makes "was it ever visited" a single column.
		//
		// COST
		//   2 calc_rmsd calls per chromosome per generation (raw + symmetry).
		//   Gated OFF by default; opt in with FLEXAIDDS_GENTRACE=1, and thin with
		//   FLEXAIDDS_GENTRACE_EVERY=N to trace every Nth generation. Requires
		//   refstructure==1 (a native pose to measure against) -- audit only,
		//   never the benchmark hot path.
		//
		// ACTIVATION (all three are required -- GENTRACE alone is silently inert)
		//   FLEXAIDDS_GENTRACE=1
		//   FLEXAIDDS_DUMP_POP=1        <-- sets refstructure; without it the gate
		//                                   below is false and nothing is written
		//   FLEXAIDDS_RMSDST=<crystal.sdf>
		//   refstructure is only ever set by load_dump_pop_refstructure()
		//   (top.cpp), which no-ops unless DUMP_POP is truthy. Setting GENTRACE on
		//   its own produces no file and no warning.
		if (FA->refstructure == 1) {
			static const int gentrace_every = [](){
				const char* on = std::getenv("FLEXAIDDS_GENTRACE");
				if (!on || on[0] == '\0' || std::atoi(on) == 0) return 0;
				const char* ev = std::getenv("FLEXAIDDS_GENTRACE_EVERY");
				const int n = (ev && std::atoi(ev) > 0) ? std::atoi(ev) : 1;
				return n;
			}();
			if (gentrace_every > 0 && ((i + 1) % gentrace_every == 0 || i == 0)) {
				// thread_local, not static: GA() runs inside an OpenMP parallel-for
				// under --parallel-dock (ParallelDock.cpp:115), where plain statics
				// would be a data race and two threads would clobber one file.
				// The thread id in the filename keeps per-subgrid traces separate.
				// RAII so the handle is closed at thread exit rather than leaked
				// for the process lifetime (one FD per GA thread under
				// --parallel-dock).
				struct GtFile {
					FILE* fp = nullptr;
					~GtFile() { if (fp) fclose(fp); }
				};
				thread_local GtFile gt_file;
				FILE*& gt_fp = gt_file.fp;
				thread_local double gt_best_raw  = 1e9;
				thread_local double gt_best_sym  = 1e9;
				thread_local int    gt_best_gen  = -1;
				if (!gt_fp) {
					char gt_path[MAX_PATH__];
#ifdef _OPENMP
					const int gt_tid = omp_in_parallel() ? omp_get_thread_num() : -1;
#else
					const int gt_tid = -1;
#endif
					if (gt_tid >= 0)
						snprintf(gt_path, MAX_PATH__, "%s.t%d.gentrace.tsv",
						         FA->rrgfile, gt_tid);
					else
						snprintf(gt_path, MAX_PATH__, "%s.gentrace.tsv", FA->rrgfile);
					gt_fp = fopen(gt_path, "w");
					if (gt_fp) {
						fprintf(gt_fp, "generation\tn_chrom\tmin_rmsd_raw\tmin_rmsd_sym\t"
						               "cf_of_min_rmsd\tbest_cf\trmsd_of_best_cf\t"
						               "best_rmsd_raw_so_far\tbest_rmsd_sym_so_far\t"
						               "gen_of_best\n");
						fprintf(stdout, "[GENTRACE] writing %s (every %d gen)\n",
						        gt_path, gentrace_every);
					} else {
						fprintf(stderr, "WARNING: [GENTRACE] cannot open %s\n", gt_path);
					}
				}
				if (gt_fp) {
					// STATE HYGIENE. calc_rmsd() is not a pure function: it writes
					// atoms[].dis/ang/dih, calls buildcc() to rebuild atoms[].coor,
					// may call alter_mode() on atoms[]/residue[], and writes
					// residue[].rot for rotamer genes. cluster.cpp's DUMP_POP can
					// ignore that ("runs last: only frees follow"); here the GA
					// continues for another generation, so every mutated field must
					// be put back. GB->num_genes == FA->npar (gaboom.cpp:281), so
					// one length covers opt_par.
					std::vector<double> saved_par(FA->opt_par,
					                              FA->opt_par + GB->num_genes);
					std::vector<int> saved_rot(static_cast<size_t>(FA->res_cnt) + 1);
					for (int r = 1; r <= FA->res_cnt; ++r) saved_rot[r] = residue[r].rot;
					double min_raw = 1e9, min_sym = 1e9, cf_at_min = 0.0;
					double best_cf = 1e9, rmsd_at_best_cf = 0.0;
					for (int c = 0; c < GB->num_chrom; ++c) {
						for (int g = 0; g < GB->num_genes; ++g)
							FA->opt_par[g] = (*chrom)[c].genes[g].to_ic;
						bool Hung = false;
						const double rr = calc_rmsd(FA, atoms, residue, *cleftgrid,
						                            FA->npar, FA->opt_par, Hung);
						Hung = true;
						const double rs = calc_rmsd(FA, atoms, residue, *cleftgrid,
						                            FA->npar, FA->opt_par, Hung);
						const double cf = (*chrom)[c].evalue;
						if (rr < min_raw) { min_raw = rr; min_sym = rs; cf_at_min = cf; }
						if (cf < best_cf) { best_cf = cf; rmsd_at_best_cf = rr; }
					}
					if (min_raw < gt_best_raw) {
						gt_best_raw = min_raw; gt_best_sym = min_sym; gt_best_gen = i + 1;
					}
					// best_rmsd_sym_so_far is the symmetry-corrected counterpart of
					// best_rmsd_raw_so_far, recorded at the same generation. It is
					// the value to compare against the 2 A success criterion, which
					// is symmetry-corrected.
					fprintf(gt_fp,
					        "%d\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\n",
					        i + 1, GB->num_chrom, min_raw, min_sym, cf_at_min,
					        best_cf, rmsd_at_best_cf, gt_best_raw, gt_best_sym,
					        gt_best_gen);
					fflush(gt_fp);   // partial trace survives a timeout kill

					// Restore in reverse order of mutation: opt_par, then rotamers,
					// then re-derive atoms[] geometry from the restored opt_par via
					// the same code path that perturbed it. The final calc_rmsd call
					// is discarded -- it is invoked for its buildcc() side effect.
					std::copy(saved_par.begin(), saved_par.end(), FA->opt_par);
					for (int r = 1; r <= FA->res_cnt; ++r) residue[r].rot = saved_rot[r];
					(void)calc_rmsd(FA, atoms, residue, *cleftgrid,
					                FA->npar, FA->opt_par, false);
				}
			}
		}

		// ── Record generation wall-clock ──
		{
			auto _t1_gen = std::chrono::steady_clock::now();
			double _ms = std::chrono::duration<double,std::milli>(_t1_gen - _t0_gen).count();
			_sum_gen_ms += _ms;
			++_n_gen_timed;
			if (i < 20)
				fprintf(stderr, "TIMING GEN %4d: %.2f ms  (~%.3f us/eval, %d chrom)\n",
				        i+1, _ms, _ms*1000.0/(2.0*GB->num_chrom), GB->num_chrom);
		}
	}

	if (do_anneal) {
		FA->temperature = target_temperature_K;
		FA->beta        = 1.0 / target_temperature_d;
		fprintf(stderr,
		        "[ANNEAL] restored target temperature %.0f K for post-GA thermodynamics\n",
		        target_temperature_d);
	}

	// ── Timing summary ──
	if (_n_gen_timed > 0) {
		double _avg = _sum_gen_ms / _n_gen_timed;
		fprintf(stderr,
		        "TIMING SUMMARY: %d gens timed, avg %.2f ms/gen, "
		        "~%.2f us/eval (2x-pop est), "
		        "est %.1f s for %d-gen x %d-chrom run\n",
		        _n_gen_timed, _avg,
		        _avg * 1000.0 / (2.0 * GB->num_chrom),
		        _avg * GB->max_generations / 1000.0,
		        GB->max_generations, GB->num_chrom);
	}

	printf("%d ligand conformers rejected\n", nrejected);
	if (entropy_converged)
		printf("GA terminated early by entropy convergence\n");
	if (ga_stagnant)
		printf("GA terminated early by fitness stagnation\n");

	// Wave 3.4: consume FA->use_memetic (set only when MEMETIC+WALL_PILOT_PASS).
	// Full burial-safe local-refine kernel is E5 and stays off until wall PASS.
	// When armed, record a durable post-GA marker on FA (observable by receipts /
	// diagnostics) so the gate is not warn-only — refine kernel hooks here later.
	if (FA->use_memetic) {
		// Sticky arm marker: non-zero means the GA path actually read use_memetic.
		// Value records max_generations budget (diagnostics only; not a refine step).
		FA->memetic_armed_at_gen =
		    GB->max_generations > 0 ? GB->max_generations : 1;
		fprintf(stderr,
		        "[MEMETIC] use_memetic=1 armed_at_gen=%d: post-GA local refine "
		        "ARMED (kernel deferred until wall PASS + E5 design)\n",
		        FA->memetic_armed_at_gen);
	} else {
		FA->memetic_armed_at_gen = 0;
	}

	// Print H_final for two-pass benchmark script parsing
	{
		std::vector<double> hfinal_energies(GB->num_chrom);
		for (int c = 0; c < GB->num_chrom; ++c)
			hfinal_energies[c] = (*chrom)[c].evalue;
		const double H_final_val = shannon_thermo::compute_shannon_entropy(
			hfinal_energies, shannon_thermo::DEFAULT_HIST_BINS);
		printf("H_final = %.6f\n", H_final_val);
	}

	QuickSort((*chrom),0,GB->num_chrom-1,true);

	// ── Legacy mixed-domain thermodynamic diagnostics (proxy-only) ──
	if (FA->thermo_engine_enabled && FA->thermo_engine != nullptr) {
		std::vector<std::vector<float>> gene_pop(GB->num_chrom,
			std::vector<float>(GB->num_genes));
		std::vector<float> cf_pop(GB->num_chrom);

		for (int c = 0; c < GB->num_chrom; ++c) {
			cf_pop[c] = static_cast<float>((*chrom)[c].evalue);
			for (int g = 0; g < GB->num_genes; ++g)
				gene_pop[c][g] = static_cast<float>((*chrom)[c].genes[g].to_int32)
				                 / static_cast<float>(MAX_RANDOM_VALUE);
		}

		// (c) Bound complex H(ω) — materialise rank-0 chromosome, compute ligand ANM
		int thermo_n_heavy = 0;
		{
			const int lig_start    = (FA->resligand && FA->resligand->fatm)
			                         ? FA->resligand->fatm[0] : -1;
			const int lig_end_incl = (FA->resligand && FA->resligand->latm)
			                         ? FA->resligand->latm[0] : -1;

			// Materialise rank-0 Cartesian coordinates into atoms[].coor[].
			// eval_chromosome → ic2cf → buildcc writes the docked pose.
			std::vector<gene> rank0_genes(GB->num_genes);
			for (int g = 0; g < GB->num_genes; ++g)
				rank0_genes[g] = (*chrom)[0].genes[g];
			eval_chromosome(FA, GB, VC, *gene_lim, atoms, residue,
			                *cleftgrid, rank0_genes.data(), target);

			if (lig_start >= 0 && lig_end_incl >= lig_start) {
				// Count heavy atoms (not H) for per-heavy-atom H_vct normalisation.
				auto is_h = [](const atom& a) noexcept {
					const char* e = a.element;
					while (*e == ' ') ++e;
					return e[0] == 'H' && (e[1] == '\0' || e[1] == ' ');
				};
				for (int ai = lig_start; ai <= lig_end_incl; ++ai)
					if (!is_h(atoms[ai])) ++thermo_n_heavy;

				fprintf(stderr,
				        "[THERMO-DBG] lig_start=%d lig_end=%d n_heavy=%d "
				        "coor0=(%.3f,%.3f,%.3f)\n",
				        lig_start, lig_end_incl, thermo_n_heavy,
				        atoms[lig_start].coor[0],
				        atoms[lig_start].coor[1],
				        atoms[lig_start].coor[2]);

				tencm::TorsionalENM lig_enm_bound;
				lig_enm_bound.build_from_ligand(atoms, lig_start, lig_end_incl + 1);
				if (lig_enm_bound.is_built()) {
					std::vector<double> eigs;
					eigs.reserve(lig_enm_bound.modes().size());
					for (const auto& nm : lig_enm_bound.modes())
						if (nm.eigenvalue > 0.0) eigs.push_back(nm.eigenvalue);
					if (!eigs.empty()) {
						const std::vector<std::vector<double>> single = { eigs };
						FA->H_rep_bound_complex = static_cast<float>(
							vibentropy::compute_vib_entropy_collapse(single).H_pop);
					}
				}
				fprintf(stderr,
				        "[THERMO-DBG] is_built=%d H_rep_bound=%.6f H_rep_ref=%.6f\n",
				        (int)lig_enm_bound.is_built(),
				        FA->H_rep_bound_complex,
				        FA->H_rep_receptor_ref + FA->H_rep_ligand_ref);
			}
		}

		FA->thermo_result = FA->thermo_engine->compute(
			gene_pop, cf_pop, FA->H_rep_bound_complex, thermo_n_heavy,
			FA->thermo_report_T);

		printf("[THERMO] claim_validity=proxy_only energy_domain=cf_arbitrary_units "
		       "ensemble_measure=optimizer_samples G_bind=%.6f H_vct=%.6f H_vct_raw=%.6f n_heavy=%d "
		       "TdS_shannon=%.6f TdS_vib=%.6f D_vib=%.6f compensation=%.6f\n",
		       FA->thermo_result.G_bind,
		       FA->thermo_result.H_vct,
		       FA->thermo_result.H_vct_raw,
		       FA->thermo_result.n_heavy_atoms,
		       FA->thermo_result.TdS_shannon,
		       FA->thermo_result.TdS_vib,
		       FA->H_rep_bound_complex,
		       FA->thermo_result.compensation);

		// Reporting-only whiteboard diagnostics — computed at thermo_report_T
		// (default 21.0 = kT_ISMB), independent of thermo_T_eff/G_bind above.
		// Left-hand labelling per whiteboard: G_bind_T21/I_ES/CF_r2s/regime are
		// all "(T=21)"-defined quantities, not RHS substitutions. Additive
		// output only; nothing here is read back into scoring or GA state.
		printf("[THERMO2] report_T=%.6f I_ES=%.6f CF_r2s=%.6f regime=%s\n",
		       FA->thermo_result.report_T,
		       FA->thermo_result.I_ES,
		       FA->thermo_result.CF_r2s,
		       FA->thermo_result.binding_regime.c_str());

		// dG_eff = <CF> - T*H over optimizer records at two score-scale
		// parameters. Even when FLEXAIDDS_THERMO_SCORE=1 computes the legacy
		// sentinel, the later exact-CF rescore/clustering path does not consume it.
		printf("[THERMO3] dG_eff=%.6f mean_CF=%.6f H=%.6f T_eff=%.6f | "
		       "dG_eff_T21=%.6f mean_CF_T21=%.6f H_T21=%.6f report_T=%.6f | "
		       "n_poses=%d thermo_score=%d\n",
		       FA->thermo_result.dG_eff,
		       FA->thermo_result.mean_CF,
		       FA->thermo_result.H_pose,
		       FA->thermo_result.T_eff_used,
		       FA->thermo_result.dG_eff_T21,
		       FA->thermo_result.mean_CF_T21,
		       FA->thermo_result.H_pose_T21,
		       FA->thermo_result.report_T,
		       GB->num_chrom,
		       flexaids::thermo_score_enabled() ? 1 : 0);

		if (flexaids::thermo_score_enabled()) {
			printf("[THERMO_GATE_SUMMARY] enforced_in_final_election=0 impossible=%d n_impossible_poses=%d "
			       "gate_dS=%.4f dG_eff=%.6f\n",
			       FA->thermo_result.thermo_impossible ? 1 : 0,
			       FA->thermo_result.n_impossible_poses,
			       FA->thermo_result.gate_dS_used,
			       FA->thermo_result.dG_eff);
		}
	}

	snprintf(outfile,MAX_PATH__,"%s_par.res",FA->rrgfile);
	if (FA->htpmode == false) {write_par((*chrom),(*gene_lim),i+1,outfile,GB->num_chrom,GB->num_genes);}

	printf("sorting chrom_snapshot\n");

	// ── InStreamClustering: finalize and report ──
	{
		auto medoids = instream_cluster.finalize();
		printf("--- InStreamClustering: %d clusters from %lld total merges ---\n",
		       instream_cluster.cluster_count(),
		       (long long)instream_cluster.total_merged());
		if (!medoids.empty()) {
			printf("  Best cluster score: %.4f (members: %d, first seen: gen %d)\n",
			       medoids[0].best_score,
			       medoids[0].member_count,
			       medoids[0].first_seen_gen);
		}
	}
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
		statmech::StatMechEngine sme(
			T_K, statmech::make_contact_function_optimizer_provenance());
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
				statmech::StatMechEngine sme_filtered(
					T_K, statmech::make_contact_function_optimizer_provenance());
				for (size_t idx : sc_indices)
					sme_filtered.add_sample((*chrom_snapshot)[idx].evalue);

				printf("--- SuperCluster pre-filter: %zu / %d poses selected ---\n",
				       sc_indices.size(), n_chrom_snapshot);
				sme = sme_filtered;
			}
		}

		statmech::Thermodynamics td = sme.compute();
		printf("--- CF-proxy ensemble diagnostics (T parameter = %.1f K, N = %d records) ---\n",
		       td.temperature, n_chrom_snapshot);
		printf("  claim_validity            = proxy_only\n");
		printf("  F-like proxy          F~  = %10.4f [legacy transform]\n", td.free_energy);
		printf("  Mean CF proxy       <CF>  = %10.4f [CF units]\n", td.mean_energy);
		printf("  CF std dev          sigma = %10.4f [CF units]\n", td.std_energy);
		printf("  C_v-like diagnostic       = %10.4f [proxy scale]\n", td.heat_capacity);
		printf("  S-like diagnostic         = %10.6f [proxy scale/K]\n", td.entropy);

		// ── Enthalpy-Entropy Index (Williams et al. 2017, Drug Discov. Today) ──
		// I_EE = (ΔH + T·ΔS) / ΔG   — diagnostic only, never for ranking
		{
			const statmech::ThermodynamicBreakdown bd = sme.compute_breakdown();
			if (bd.has_I_EE) {
				printf("  Legacy H/S ratio      I_EE= %10.4f  [proxy_only; diagnostic]\n",
				       bd.I_EE);
			}
		}

		// ── Kirchhoff ΔG(T) extrapolation (Robertson & Murphy 1997) ────────
		// Activated only when DSF/TSA Tm has been supplied via dsf_Tm_K.
		// ΔG(T) = ΔHm(1 − T/Tm) − ΔCp[(Tm − T) + T·ln(T/Tm)]
		if (FA->dsf_Tm_K > 0.0 && td.allows_canonical_physical_claim()) {
			thermal_extrap::KirchhoffInput kin;
			kin.Tm_K     = FA->dsf_Tm_K;
			kin.delta_Hm = (FA->dsf_delta_Hm != 0.0)
			                 ? FA->dsf_delta_Hm
			                 : td.mean_energy; // fallback: use computed ⟨E⟩ at Tm
			// ΔCp: requires two-temperature runs; use Cv as provisional stand-in
			// with a warning. Will be replaced once compute_delta_Cp() is wired
			// to a paired apo run.
			kin.delta_Cp = td.heat_capacity; // provisional — Cv, not ΔCp
			const double T_target = (FA->temperature > 0)
			                         ? static_cast<double>(FA->temperature)
			                         : 298.15;
			const auto kr = thermal_extrap::kirchhoff_deltaG(kin, T_target);
			printf("--- Kirchhoff ΔG(T) [Robertson & Murphy 1997] ---\n");
			printf("  Tm (DSF/TSA)          = %10.2f K  (%.1f °C)\n",
			       kin.Tm_K, kin.Tm_K - 273.15);
			printf("  ΔHm (at Tm)           = %10.4f kcal/mol\n", kin.delta_Hm);
			printf("  ΔCp (provisional Cv)  = %10.6f kcal/(mol·K)  [⚠ use compute_delta_Cp]\n",
			       kin.delta_Cp);
			printf("  ΔG(%.1f K)             = %10.4f kcal/mol\n", T_target, kr.delta_G);
			printf("  ΔH(%.1f K)             = %10.4f kcal/mol\n", T_target, kr.delta_H);
			printf("  T·ΔS(%.1f K)           = %10.4f kcal/mol\n", T_target, kr.T_delta_S);
		} else if (FA->dsf_Tm_K > 0.0) {
			printf("--- Kirchhoff extrapolation unavailable: CF proxy lacks calibrated energy/measure provenance ---\n");
		}

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

				printf("--- ShannonThermoStack legacy mixed-domain diagnostic (proxy_only) ---\n");
				printf("  Shannon conf entropy    = %10.4f nats\n", ftr.shannonEntropy);
				printf("  Torsional vib diagnostic= %10.6f [model scale/K]\n", ftr.torsionalVibEntropy);
				printf("  Entropy contribution    = %10.4f [mixed proxy scale]\n", ftr.entropyContribution);
				printf("  Legacy deltaG field     = %10.4f [proxy; configurational entropy is double-counted]\n", ftr.deltaG);
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
					printf("  Final co-translational diagnostic = %10.4f [model scale; proxy_only]\n",
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

	// Ring pucker side-channel travels with the chromosome (POD copy) so
	// snapshots and dedup-compaction keep ring genes synced with standard genes.
	memcpy(dest->ring_phases, src->ring_phases, sizeof(dest->ring_phases));
	memcpy(dest->ring_six,    src->ring_six,    sizeof(dest->ring_six));
	memcpy(dest->ring_five,   src->ring_five,   sizeof(dest->ring_five));
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

	// ── GA-PARALLEL-EVAL (P3, finding #1) ──────────────────────────────────
	// Default ON (drift allowed — see docs/FLEXAID_FAST_DOCKING_PLAN.md §P3 and
	// OPTIMIZATION_KNOWN_ISSUES.md): offspring CF evaluation is deferred so that
	// BOTH parents and offspring are scored in the single OpenMP-parallel loop
	// inside calculate_fitness() (called unconditionally at the end of this
	// function for both STEADY and BOOM). Each new offspring's status is set to
	// ' ' ("needs eval") — the STALE-STATUS FIX — because chrom[num_chrom+i] is
	// REUSED memory whose status is typically 'n' from the previous generation,
	// and calculate_fitness()'s eval loop SKIPS status=='n'. Without the explicit
	// reset the deferred offspring would keep the prior occupant's stale CF
	// (the CF -5.23 vs serial -51.93 defect documented in OPTIMIZATION_KNOWN_ISSUES.md).
	// No new parallel region and no change to eval_chromosome()/vcfunction()
	// math: this flag only decides WHERE (serial reproduce() vs parallel
	// calculate_fitness()) the exact same eval is called from.
	//
	// Override: FLEXAIDDS_PARALLEL_REPRODUCE=0 forces the legacy serial inline
	// path (bit-reproducible reference); any other explicit value keeps it ON.
	// The CI A/B reproducibility mode is FLEXAID_DETERMINISTIC (see
	// calculate_fitness()), which pins a serial-equivalent reduction order for
	// the parallel eval loop rather than disabling the deferred path here.
	// Default OFF per METHODOLOGY §1: this flag was gated off for reproducibility
	// (the ~0.2% multi-thread chromosome divergence in OPTIMIZATION_KNOWN_ISSUES.md),
	// and §1 requires an intended behaviour change to be opt-in behind a flag that
	// defaults OFF with parity holding when it is OFF. The drift allowance that
	// would unblock it is a maintainer decision, not one this change can grant
	// itself. Set FLEXAIDDS_PARALLEL_REPRODUCE=1 to opt in and benchmark it.
	static const bool parallel_reproduce_eval = flexaids::parallel_reproduce_enabled();

	// Multi-chain VCT normalisation (see GA() comment for rationale)
	const int n_receptor_chains = count_receptor_chains(FA, residue);
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

		const bool did_cross = RandomDouble() < crossprob;
		if(did_cross){
			crossover(chrop1_gen,chrop2_gen,GB->num_genes,GB->intragenes);
		}

		/************************************/
		/****** MUTATION OPERATOR  ********/
		/************************************/
		num_genes_wo_sc = GB->num_genes-FA->nflxsc_real;

		mutate(chrop1_gen,GB->num_genes-FA->nflxsc_real,mutprob,gene_lim);
		k=0;
		for(j=0;j<FA->nflxsc;j++){
			if(residue[FA->flex_res[j].inum].trot != 0){
				if(RandomDouble() < FA->flex_res[j].prob){
					mutate(&chrop1_gen[num_genes_wo_sc+k],1,mutprob,
					       gene_lim ? &gene_lim[num_genes_wo_sc+k] : nullptr);
				}
				k++;
			}
		}

		mutate(chrop2_gen,GB->num_genes-FA->nflxsc_real,mutprob,gene_lim);
		k=0;
		for(j=0;j<FA->nflxsc;j++){
			if(residue[FA->flex_res[j].inum].trot != 0){
				if(RandomDouble() < FA->flex_res[j].prob){
					mutate(&chrop2_gen[num_genes_wo_sc+k],1,mutprob,
					       gene_lim ? &gene_lim[num_genes_wo_sc+k] : nullptr);
				}
				k++;
			}
		}

		for(j=0; j<GB->num_genes; j++){
			chrop1_gen[j].to_ic = genetoic(&gene_lim[j],chrop1_gen[j].to_int32);
			chrop2_gen[j].to_ic = genetoic(&gene_lim[j],chrop2_gen[j].to_int32);
		}

		// ── Ring pucker crossover + mutation (mirrors the standard-gene flow) ──
		// Derive each child's ring genes from its parent, cross the two children
		// when the standard genes were crossed, then mutate. rc1/rc2 are scratch
		// chromosomes (only ring_* fields used; genes pointer stays null/unused).
		chromosome rc1{}, rc2{};
		if (FA->ring_flex_active) {
			memcpy(rc1.ring_phases, chrom[p1].ring_phases, sizeof(rc1.ring_phases));
			memcpy(rc1.ring_six,    chrom[p1].ring_six,    sizeof(rc1.ring_six));
			memcpy(rc1.ring_five,   chrom[p1].ring_five,   sizeof(rc1.ring_five));
			memcpy(rc2.ring_phases, chrom[p2].ring_phases, sizeof(rc2.ring_phases));
			memcpy(rc2.ring_six,    chrom[p2].ring_six,    sizeof(rc2.ring_six));
			memcpy(rc2.ring_five,   chrom[p2].ring_five,   sizeof(rc2.ring_five));
			if (did_cross) ring_crossover_chrom(FA, &rc1, &rc2);
			ring_mutate_chrom(FA, &rc1);
			ring_mutate_chrom(FA, &rc2);
		}

		// A: phenotype-bin uniqueness when enabled (classic path); default keeps
		// historical hash_genes on rounded to_ic.
		const bool pheno_unique = flexaids::new_search::phenotype_unique_enabled();
		size_t sig1 = pheno_unique && gene_lim
		                  ? flexaids::new_search::hash_phenotype_bins(
		                        chrop1_gen, GB->num_genes, gene_lim)
		                  : hash_genes(chrop1_gen, GB->num_genes);
		size_t sig2 = pheno_unique && gene_lim
		                  ? flexaids::new_search::hash_phenotype_bins(
		                        chrop2_gen, GB->num_genes, gene_lim)
		                  : hash_genes(chrop2_gen, GB->num_genes);

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

			if (FA->ring_flex_active) {
				memcpy(chrom[GB->num_chrom+i].ring_phases, rc1.ring_phases, sizeof(rc1.ring_phases));
				memcpy(chrom[GB->num_chrom+i].ring_six,    rc1.ring_six,    sizeof(rc1.ring_six));
				memcpy(chrom[GB->num_chrom+i].ring_five,   rc1.ring_five,   sizeof(rc1.ring_five));
				// Deferred to calculate_fitness()'s per-thread FA copy when the
				// eval itself is deferred (see parallel_reproduce_eval below):
				// loading into the shared FA here would race with the other
				// offspring evals still running serially in this while-loop.
				if (!parallel_reproduce_eval)
					ring_load_chrom_to_fa(FA, &chrom[GB->num_chrom+i]);
			}

			if (parallel_reproduce_eval) {
				// GA-PARALLEL-EVAL: defer CF eval to calculate_fitness()'s
				// OpenMP loop. STALE-STATUS FIX — explicitly mark this reused
				// slot as ' ' ("needs eval"); it typically holds 'n' from the
				// previous generation, which the eval loop skips.
				chrom[GB->num_chrom+i].status=' ';
			} else {
				chrom[GB->num_chrom+i].cf=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,
									  chrom[GB->num_chrom+i].genes,target);
				chrom[GB->num_chrom+i].evalue=get_cf_evalue(&chrom[GB->num_chrom+i].cf, FA) / n_receptor_chains;
				chrom[GB->num_chrom+i].app_evalue=get_apparent_cf_evalue(&chrom[GB->num_chrom+i].cf) / n_receptor_chains;
				ccbm_inject_strain(FA, chrom[GB->num_chrom+i], gene_lim);  // CCBM strain
				chrom[GB->num_chrom+i].status='n';
			}

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

			if (FA->ring_flex_active) {
				memcpy(chrom[GB->num_chrom+i].ring_phases, rc2.ring_phases, sizeof(rc2.ring_phases));
				memcpy(chrom[GB->num_chrom+i].ring_six,    rc2.ring_six,    sizeof(rc2.ring_six));
				memcpy(chrom[GB->num_chrom+i].ring_five,   rc2.ring_five,   sizeof(rc2.ring_five));
				// See offspring #1 above: deferred under parallel_reproduce_eval.
				if (!parallel_reproduce_eval)
					ring_load_chrom_to_fa(FA, &chrom[GB->num_chrom+i]);
			}

			if (parallel_reproduce_eval) {
				// GA-PARALLEL-EVAL: see offspring #1 above. STALE-STATUS FIX.
				chrom[GB->num_chrom+i].status=' ';
			} else {
				chrom[GB->num_chrom+i].cf=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,
									  chrom[GB->num_chrom+i].genes,target);
				chrom[GB->num_chrom+i].evalue=get_cf_evalue(&chrom[GB->num_chrom+i].cf, FA) / n_receptor_chains;
				chrom[GB->num_chrom+i].app_evalue=get_apparent_cf_evalue(&chrom[GB->num_chrom+i].cf) / n_receptor_chains;
				ccbm_inject_strain(FA, chrom[GB->num_chrom+i], gene_lim);  // CCBM strain
				chrom[GB->num_chrom+i].status='n';
			}

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

	// Multi-chain VCT normalisation (see GA() comment for rationale)
	const int n_receptor_chains = count_receptor_chains(FA, residue);

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
	// Runtime dispatch keeps production GA fitness on the CPU/OpenMP path.
	// The legacy GPU evaluators below remain compiled for explicit experiments,
	// but default backend selection does not use them until their raw gene
	// decoding and CF terms match the full ic2cf/vcfunction CPU path.

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

	// Helper lambda: unpack full-fidelity GPU batch results (P4) into chromosome
	// CF structures. The GPU now fills every scoring channel that feeds
	// get_cf_evalue(): com/wal/sas/con/elec/hbond/gist_desolv/pb_clash. On the
	// Metal MULTI (screening) path only com/wal/sas are populated — the other
	// vectors are left zero by the caller, which is the documented reduced
	// fidelity for that path. Channels with NO GPU implementation (cf.gist is the
	// non-scoring diagnostic channel; metal_coord/entropy/h_rep) are zeroed here
	// so a stale host-buffer value can never leak into get_cf_evalue().
	auto unpack_gpu_results = [&](const std::vector<double>& h_com,
	                              const std::vector<double>& h_wal,
	                              const std::vector<double>& h_sas,
	                              const std::vector<double>& h_con,
	                              const std::vector<double>& h_elec,
	                              const std::vector<double>& h_hbond,
	                              const std::vector<double>& h_gist,
	                              const std::vector<double>& h_pb) {
		for (int c = 0; c < pop_size; ++c) {
			if (chrom[c].status != 'n') {
				chrom[c].cf.com         = h_com[c];
				chrom[c].cf.wal         = h_wal[c];
				chrom[c].cf.sas         = h_sas[c];
				chrom[c].cf.con         = h_con[c];
				chrom[c].cf.elec        = h_elec[c];
				chrom[c].cf.hbond       = h_hbond[c];
				chrom[c].cf.gist_desolv = h_gist[c];   // scoring GIST channel (get_cf_evalue sums gist_desolv, not cf.gist)
				chrom[c].cf.pb_clash    = h_pb[c];
				// Terms with no GPU implementation — zeroed (see divergence guard).
				chrom[c].cf.gist        = 0.0;
				chrom[c].cf.metal_coord = 0.0;
				chrom[c].cf.entropy     = 0.0;
				chrom[c].cf.h_rep       = 0.0;
				chrom[c].cf.totsas = 0.0;
				chrom[c].cf.rclash = (h_wal[c] > CLASH_THRESHOLD) ? 1 : 0;
				chrom[c].evalue     = get_cf_evalue(&chrom[c].cf, FA) / n_receptor_chains;
				chrom[c].app_evalue = get_apparent_cf_evalue(&chrom[c].cf) / n_receptor_chains;
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

	// Build the per-batch full-fidelity CF scalar params from FA (+ the CPU env
	// gates it mirrors), so the GPU CF assembly matches get_cf_evalue().
	auto build_gpu_params = [&]() -> GpuCfParams {
		GpuCfParams P; std::memset(&P, 0, sizeof(P));
		// Distance weighting: mirror vcfunction's resolution exactly.
		double dw = 0.0;
		if (std::getenv("FLEXAIDDS_DIST_WEIGHT_CON") != nullptr) {
			float r0 = 3.5f;
			if (const char* s = std::getenv("FLEXAIDDS_CON_R0")) {
				float pr = strtof(s, nullptr); if (pr > 0.0f) r0 = pr;
			}
			dw = r0;
		} else if (FA->vct_dist_weight_r0 > 0.0) {
			dw = FA->vct_dist_weight_r0;
		}
		P.dw_r0       = (float)dw;
		P.elec_on     = FA->use_elec ? 1 : 0;
		P.dielectric  = FA->dielectric;
		// H-bond enters GA fitness only under the same gate as vcfunction.cpp.
		const bool hb_on = FA->use_hbond && FA->use_hbond_search && !FA->hbond_rank_rescore;
		P.hbond_weight      = hb_on ? (float)FA->hbond_weight : 0.0f;
		P.hbond_salt_weight = hb_on ? (float)FA->hbond_salt_bridge_weight : 0.0f;
		P.hbond_opt_dist    = (float)FA->hbond_optimal_dist;
		P.hbond_sigma_dist  = (float)FA->hbond_sigma_dist;
		// Representative angle term (drift model — the GPU omits the true
		// virtual-H directionality; see cuda_eval.cuh error note).
		P.hbond_angle_repr  = 0.7f;
		P.pb_clash_weight   = (float)FA->pb_clash_weight;
		P.pb_clash_ratio    = (float)FA->pb_clash_ratio;
		P.pb_clash_exponent = (float)FA->pb_clash_exponent;
		P.pb_pocket_weight  = (float)FA->pb_pocket_weight;
		P.pb_pocket_radius  = (float)FA->pb_pocket_radius;
		P.kdist       = (float)KDIST;
		P.sas_weight  = (float)FA->sas_weight;
		P.solvent_flat= (FA->solventterm != 0.0f) ? 1 : 0;
		P.solventterm = (float)FA->solventterm;
		return P;
	};

	// Build the rigid full-fidelity static device inputs (pb-vdw radii, charges,
	// H-bond donor/acceptor+heavy flags, covalent constraints). GIST grid upload
	// is intentionally left unwired (GISTGrid exposes no grid accessors) so GPU
	// gist_desolv stays 0 — covered by the divergence guard. Same atom indexing
	// convention as prepare_gpu_atoms (atoms[a], a in [0, atm_cnt_real)).
	struct GPUExtraData {
		std::vector<float> pbvdw, charge;
		std::vector<int>   hflags;
		std::vector<int>   cons_i, cons_j;
		std::vector<float> cons_bl, cons_md;
		GpuCfExtraStatic   extra;
	};
	auto prepare_gpu_extra = [&]() -> GPUExtraData {
		const int n_atoms = FA->atm_cnt_real;
		GPUExtraData d;
		d.pbvdw.resize(n_atoms);
		d.charge.resize(n_atoms);
		d.hflags.resize(n_atoms);
		for (int a = 0; a < n_atoms; ++a) {
			d.pbvdw[a]  = (float)atoms[a].pb_vdw_radius;
			d.charge[a] = atoms[a].has_resp ? atoms[a].resp_charge : atoms[a].charge;
			// type256 bit layout (atom_typing_256.h): bit7 = H-bond donor,
			// bit6 = H-bond acceptor. bit2 (local) = heavy (non-hydrogen).
			int f = 0;
			const uint8_t t256 = atoms[a].type256;
			if ((t256 >> 7) & 0x1) f |= 0x1;   // donor
			if ((t256 >> 6) & 0x1) f |= 0x2;   // acceptor
			const bool is_h = (atoms[a].element[0] == 'H') ||
			                  (atoms[a].element[0] == ' ' && atoms[a].element[1] == 'H');
			if (!is_h) f |= 0x4;               // heavy
			d.hflags[a] = f;
		}
		// Covalent (type-1) constraints only; use internal atom indices (inum).
		for (int c = 0; c < FA->num_constraints; ++c) {
			if (FA->constraints[c].type != 1) continue;
			const int i = FA->constraints[c].inum1;
			const int j = FA->constraints[c].inum2;
			if (i < 0 || i >= n_atoms || j < 0 || j >= n_atoms) continue;
			d.cons_i.push_back(i);
			d.cons_j.push_back(j);
			d.cons_bl.push_back(FA->constraints[c].bond_len);
			d.cons_md.push_back(FA->constraints[c].max_dist);
		}
		std::memset(&d.extra, 0, sizeof(d.extra));
		d.extra.atom_pbvdw  = d.pbvdw.data();
		d.extra.atom_charge = d.charge.data();
		d.extra.atom_hflags = d.hflags.data();
		d.extra.n_cons = (int)d.cons_i.size();
		if (d.extra.n_cons > 0) {
			d.extra.cons_i       = d.cons_i.data();
			d.extra.cons_j       = d.cons_j.data();
			d.extra.cons_bondlen = d.cons_bl.data();
			d.extra.cons_maxdist = d.cons_md.data();
		}
		d.extra.gist_nx = 0;   // GPU GIST upload not wired (see comment above)
		return d;
	};
#endif  // FLEXAIDS_USE_CUDA || FLEXAIDS_USE_METAL

	// Log dispatch decision on first call.
	// FLEXAIDDS_FORCE_CPU (truthy, default unset): pin fitness evaluation to a CPU
	// backend. Exists so a run that needs the PoseBust terms (pb_clash / pb_pocket),
	// which are CPU-only, has a way to get them — see the divergence guard below.
	// Unset → select_backend() exactly as before, bit-identical.
	[[maybe_unused]] const auto backend = []() {
		const char* fc = std::getenv("FLEXAIDDS_FORCE_CPU");
		if (fc && fc[0] != '\0' && std::strcmp(fc, "0") != 0)
			return flexaids::select_cpu_backend();
		return flexaids::select_backend();
	}();
	if (!ctx.dispatch_logged) {
		auto report = flexaids::get_dispatch_report();
		fprintf(stderr, "[FlexAIDdS] Hardware dispatch: %s (%s)\n",
		        flexaids::backend_name(static_cast<flexaids::HardwareBackend>(
		            static_cast<uint8_t>(report.selected))), report.reason.c_str());
		ctx.dispatch_logged = true;
	}

	// ── Accelerated-path divergence guard (P4) ────────────────────────────────
	// The GPU CF now reproduces com/wal/sas/con/elec/hbond/gist_desolv/pb_clash
	// (PoseBust clash + pocket now INCLUDED) within the ranking-preserving drift
	// tolerance. A few CPU terms still have NO GPU implementation and are zeroed
	// in unpack_gpu_results() so no stale value leaks into get_cf_evalue():
	//   • cf.metal_coord   (Gaussian metal-coordination Morse term)
	//   • cf.entropy       (Shannon contact-type vct-entropy penalty)
	//   • cf.h_rep         (tENCoM vibrational Shannon entropy, tencom_weight)
	//   • cf.gist_desolv   when use_gist is set: the GPU GIST grid upload is not
	//                       yet wired (GISTGrid has no grid accessors) → 0.
	// Additionally the Metal MULTI (screening) path computes only com/wal/sas.
	// Warn once if any active weight would make the accelerated score diverge.
	if (backend == flexaids::HardwareBackend::CUDA ||
	    backend == flexaids::HardwareBackend::METAL) {
		static bool gpu_div_warned = false;
		const bool metal_c = FA->use_metal_coord != 0;
		const bool entropy = FA->vct_entropy_weight > 0.0;
		const bool tencom  = FA->tencom_weight > 0.0f;
		const bool gist    = FA->use_gist != 0;   // GPU GIST upload not wired yet
		if (!gpu_div_warned && (metal_c || entropy || tencom || gist)) {
			gpu_div_warned = true;
			fprintf(stderr,
			        "[GPU_CF] WARNING: %s accelerated path does not compute "
			        "metal_coord=%d vct_entropy(w=%.4g) tENCoM(w=%.4g) gist=%d; these "
			        "are zeroed and scoring will DIVERGE from the CPU backend. Set "
			        "FLEXAIDDS_FORCE_CPU=1 to score them on CPU.\n",
			        flexaids::backend_name(backend),
			        (int)metal_c, FA->vct_entropy_weight, FA->tencom_weight, (int)gist);
		}
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
			CudaEvalCtx* c = cuda_eval_init(n_atoms, n_types, MAX_NUM_CHROM,
			                     n_genes, ad.lig_first, ad.lig_last,
			                     FA->permeability,
			                     ad.xyz.data(), ad.type.data(),
			                     ad.radius.data(), h_emat.data());
			// Upload the rigid full-fidelity static arrays once per context.
			auto ex = prepare_gpu_extra();
			cuda_eval_set_extra(c, &ex.extra);
			return c;
		});

		std::vector<double> h_genes = pack_genes_batch(n_genes);
		std::vector<double> h_com(pop_size), h_wal(pop_size), h_sas(pop_size),
		                    h_con(pop_size), h_elec(pop_size), h_hbond(pop_size),
		                    h_gist(pop_size), h_pb(pop_size);
		GpuCfParams  params = build_gpu_params();
		GpuCfResults res;
		res.com = h_com.data(); res.wal = h_wal.data(); res.sas = h_sas.data();
		res.con = h_con.data(); res.elec = h_elec.data(); res.hbond = h_hbond.data();
		res.gist_desolv = h_gist.data(); res.pb_clash = h_pb.data();
		cuda_eval_batch(handle.ctx, pop_size, n_genes, h_genes.data(), &params, &res);
		unpack_gpu_results(h_com, h_wal, h_sas, h_con, h_elec, h_hbond, h_gist, h_pb);
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

		// Batch size: how many complexes to collect before one GPU dispatch.
		// GB->metal_batch_n == 1 → instant dispatch (same as single-batch).
		// GB->metal_batch_n == N → N concurrent GA workers share one kernel.
		const int batch_n = (GB->metal_batch_n > 0) ? GB->metal_batch_n : 1;

		// max_pop for the context must accommodate all batched chromosomes.
		const int ctx_max_pop = pop_size * batch_n;

		auto& pool = GPUContextPool::instance();
		auto handle = pool.acquire_metal(n_atoms, n_types, ctx_max_pop, [&]() {
			auto ad = prepare_gpu_atoms();
			std::vector<float> h_emat = build_emat_sampled(n_types, ns);
			MetalEvalCtx* c = metal_eval_init(n_atoms, n_types, ctx_max_pop,
			                      ad.lig_first, ad.lig_last,
			                      FA->permeability,
			                      ad.xyz.data(), ad.type.data(),
			                      ad.radius.data(), h_emat.data(), ns);
			// Upload the rigid full-fidelity static arrays once per context.
			if (c) { auto ex = prepare_gpu_extra(); metal_eval_set_extra(c, &ex.extra); }
			return c;
		});

		if (handle.ctx) {
			std::vector<double> h_genes = pack_genes_batch(n_genes);
			// All eight channels are zero-initialised; the multi (screening) path
			// fills only com/wal/sas, leaving the rest 0 (its documented fidelity).
			std::vector<double> h_com(pop_size), h_wal(pop_size), h_sas(pop_size),
			                    h_con(pop_size, 0.0), h_elec(pop_size, 0.0),
			                    h_hbond(pop_size, 0.0), h_gist(pop_size, 0.0),
			                    h_pb(pop_size, 0.0);
			GpuCfParams params = build_gpu_params();
			// Stash params so the multi-complex kernel (fixed GPUContextPool API)
			// can read dw_r0/sas_weight/solvent for its com/wal/sas evaluation.
			metal_eval_set_params(handle.ctx, &params);

			if (batch_n <= 1) {
				// Single-complex fast path — full-fidelity CF.
				GpuCfResults res;
				res.com = h_com.data(); res.wal = h_wal.data(); res.sas = h_sas.data();
				res.con = h_con.data(); res.elec = h_elec.data(); res.hbond = h_hbond.data();
				res.gist_desolv = h_gist.data(); res.pb_clash = h_pb.data();
				metal_eval_batch(handle.ctx, pop_size, n_genes, h_genes.data(),
				                 &params, &res);
			} else {
				// Multi-complex path: pack per-complex atom data and queue for
				// batched dispatch.  When N concurrent workers all reach this
				// point in the same generation, they share one GPU kernel launch
				// (N × pop_size chromosomes per dispatch). Screening pre-filter:
				// only com/wal/sas are returned (see metal_eval.h).
				auto ad = prepare_gpu_atoms();
				MetalMultiBatchEntry entry;
				entry.h_genes      = h_genes.data();
				entry.h_atom_xyz    = ad.xyz.data();
				entry.h_atom_type   = ad.type.data();
				entry.h_atom_radius = ad.radius.data();
				entry.n_atoms       = n_atoms;
				entry.lig_first     = ad.lig_first;
				entry.lig_last      = ad.lig_last;
				entry.n_types       = n_types;
				entry.perm          = FA->permeability;
				entry.h_com_out     = h_com.data();
				entry.h_wal_out     = h_wal.data();
				entry.h_sas_out     = h_sas.data();
				pool.submit_metal_batch_multi(handle.ctx, pop_size, n_genes,
				                              entry, batch_n);
			}

			unpack_gpu_results(h_com, h_wal, h_sas, h_con, h_elec, h_hbond, h_gist, h_pb);
			gpu_handled = true;
		} else {
			fprintf(stderr,
			        "[FlexAIDdS] Metal backend selected but no runtime device was available; falling back to CPU\n");
		}
		pool.release_metal(handle);
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
				for (int m = 0; m < FA->nmov[r]; ++m) {
					const int ai = FA->mov[r][m];
					if (ai >= 1 && ai <= natm) dirty_atm.push_back(ai);
				}
			// Atoms directly referenced by map_par (IC fields: dis/ang/dih)
			// Direct-mode rigid-body parameters use pseudo indices >= 90000;
			// those are not entries in atoms[] and must never enter copy lists.
			for (int p = 0; p < FA->npar; ++p) {
				const int ai = FA->map_par[p].atm;
				if (ai >= 1 && ai <= natm) dirty_atm.push_back(ai);
			}
			// Cascade dihedral atoms (atoms whose .dih depends on a flex bond)
			for (int p = 0; p < FA->npar; ++p) {
				if (FA->map_par[p].typ == 2 &&
				    FA->map_par[p].atm >= 1 && FA->map_par[p].atm <= natm) {
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
				if (FA->map_par[p].typ == 4 &&
				    FA->map_par[p].atm >= 1 && FA->map_par[p].atm <= natm)
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

		// ── P3: resident per-thread workspace (lean receptor copy) ───────────
		// The receptor is READ-ONLY during scoring, so rather than re-cloning
		// the full atom/residue arrays and every Voronoi scratch buffer on each
		// generation (the dominant memory-bandwidth cost of the OpenMP path —
		// O(generations × threads × natm)), we allocate them ONCE and keep them
		// resident across calculate_fitness() calls. The resident receptor atoms
		// stay valid because a rigid receptor never moves and flexible side-chain
		// atoms are always in the per-chromosome dirty set restored below; only
		// the small mutable ligand/pose/flex state is refreshed each generation.
		// The cheap live snapshots (FA/VC scalars, optres records) are re-synced
		// each call so nothing goes stale. The cache is fully rebuilt whenever
		// the problem shape or any base pointer changes, so multiple sequential
		// docks in one process are safe. Vcontacts/vcfunction already reset their
		// scratch per eval (contacts is epoch-stamped, ca_index re-seeded, etc.),
		// so carrying those buffers across calls matches the existing across-eval
		// reuse within a single call — no behavioural change to the CF math.
		struct ParEvalWS {
			int n_thr=0, natm=-1, nres=-1, natmr=-1, nopt=-1, nctb=-1, ca_recsize=-1;
			const void *fa=nullptr,*vc=nullptr,*atoms=nullptr,*res=nullptr,*optres=nullptr;
			std::vector<std::vector<atom>>        tl_atoms;
			std::vector<std::vector<resid>>       tl_res;
			std::vector<FA_Global>                tl_fa;
			std::vector<std::vector<int>>         tl_contacts;
			std::vector<std::vector<float>>       tl_contrib;
			std::vector<std::vector<OptRes>>      tl_optres;
			std::vector<VC_Global>                tl_vc;
			std::vector<std::vector<atomsas>>     tl_calc;
			std::vector<std::vector<int>>         tl_calclist;
			std::vector<std::vector<int>>         tl_caidx;
			std::vector<std::vector<ca_struct>>   tl_carec;
			std::vector<std::vector<int>>         tl_seed;
			std::vector<std::vector<contactlist>> tl_contlist;
			std::vector<std::vector<ptindex>>     tl_ptorder;
			std::vector<std::vector<vertex>>      tl_centerpt;
			std::vector<std::vector<vertex>>      tl_poly;
			std::vector<std::vector<plane>>       tl_cont;
			std::vector<std::vector<edgevector>>  tl_vedge;
			// Per-thread Calc[] indices with score==true (sibling VC_Global
			// fields: scorable_list, n_scorable, scorable_cap, fastpath_used).
			std::vector<std::vector<int>>         tl_scorable;
		};
		thread_local ParEvalWS ws;  // resident across generations on THIS thread.
		                       // Must not be process-wide static: --parallel-dock
		                       // runs GA() under an outer OpenMP parallel-for
		                       // (ParallelDock.cpp). A shared static races on
		                       // rebuild and tl_* buffers. Inner eval still
		                       // shares this thread's tl_* via references
		                       // captured before the inner omp for, so the
		                       // serial claim path is unchanged.

		const bool ws_valid =
			ws.n_thr == n_thr && ws.natm == natm && ws.nres == nres &&
			ws.natmr == natmr && ws.nopt == nopt && ws.nctb == nctb &&
			ws.ca_recsize == VC->ca_recsize &&
			ws.fa == (const void*)FA && ws.vc == (const void*)VC &&
			ws.atoms == (const void*)atoms && ws.res == (const void*)residue &&
			ws.optres == (const void*)FA->optres;

		if (!ws_valid) {
			// Full (re)allocation — the receptor is cloned exactly once per shape.
			ws = ParEvalWS{};
			ws.n_thr = n_thr; ws.natm = natm; ws.nres = nres; ws.natmr = natmr;
			ws.nopt = nopt;   ws.nctb = nctb; ws.ca_recsize = VC->ca_recsize;
			ws.fa = FA; ws.vc = VC; ws.atoms = atoms; ws.res = residue;
			ws.optres = FA->optres;
			ws.tl_atoms.assign(n_thr, std::vector<atom>(atoms, atoms + natm + 1));
			ws.tl_res.assign(n_thr, std::vector<resid>(residue, residue + nres + 1));
			ws.tl_fa.assign(n_thr, *FA);
			// CONTACTS_BUFFER_SIZE: stamps + the trailing epoch slot (flexaid.h).
			// This buffer is RESIDENT across generations while tl_fa[t] below is
			// re-snapshotted from *FA every generation; the epoch must therefore
			// live here, in the buffer, and not in the FA_Global copy.
			ws.tl_contacts.assign(n_thr, std::vector<int>(CONTACTS_BUFFER_SIZE, 0));
			ws.tl_contrib.assign(n_thr, std::vector<float>(nctb, 0.0f));
			ws.tl_optres.assign(n_thr, std::vector<OptRes>(FA->optres, FA->optres + nopt));
			ws.tl_vc.assign(n_thr, *VC);
			ws.tl_calc.assign(n_thr, std::vector<atomsas>(natmr));
			ws.tl_calclist.assign(n_thr, std::vector<int>(natmr));
			ws.tl_caidx.assign(n_thr, std::vector<int>(natmr, -1));
			ws.tl_carec.assign(n_thr, std::vector<ca_struct>(VC->ca_recsize));
			ws.tl_seed.assign(n_thr, std::vector<int>(3 * natmr));
			ws.tl_contlist.assign(n_thr, std::vector<contactlist>(GA_CONTLIST_SIZE));
			ws.tl_ptorder.assign(n_thr, std::vector<ptindex>(MAX_PT));
			ws.tl_centerpt.assign(n_thr, std::vector<vertex>(MAX_PT));
			ws.tl_poly.assign(n_thr, std::vector<vertex>(MAX_POLY));
			ws.tl_cont.assign(n_thr, std::vector<plane>(MAX_PT));
			ws.tl_vedge.assign(n_thr, std::vector<edgevector>(MAX_POLY));
			ws.tl_scorable.assign(n_thr, std::vector<int>(natmr, 0));
		}

		// Aliases keep the eval loop below identical to the previous per-call form.
		auto& tl_atoms   = ws.tl_atoms;    auto& tl_res      = ws.tl_res;
		auto& tl_fa      = ws.tl_fa;       auto& tl_contacts = ws.tl_contacts;
		auto& tl_contrib = ws.tl_contrib;  auto& tl_optres   = ws.tl_optres;
		auto& tl_vc      = ws.tl_vc;       auto& tl_calc     = ws.tl_calc;
		auto& tl_calclist= ws.tl_calclist; auto& tl_caidx    = ws.tl_caidx;
		auto& tl_carec   = ws.tl_carec;    auto& tl_seed     = ws.tl_seed;
		auto& tl_contlist= ws.tl_contlist; auto& tl_ptorder  = ws.tl_ptorder;
		auto& tl_centerpt= ws.tl_centerpt; auto& tl_poly     = ws.tl_poly;
		auto& tl_cont    = ws.tl_cont;     auto& tl_vedge    = ws.tl_vedge;
		auto& tl_scorable= ws.tl_scorable;

		// Wire per-thread scorable-list scratch. Sibling adds these four
		// fields at the end of VC_Global; the generic lambda keeps this
		// block compiling if the header has not landed yet.
		auto wire_scorable = [](auto& vc, int* buf, int cap) {
			if constexpr (requires {
				vc.scorable_list; vc.n_scorable;
				vc.scorable_cap; vc.fastpath_used;
			}) {
				vc.scorable_list = buf;
				vc.n_scorable = 0;
				vc.scorable_cap = cap;
				vc.fastpath_used = 0;
			}
		};

		for (int t = 0; t < n_thr; ++t) {
			// Refresh the cheap live FA snapshot each generation (scalar state may
			// change between generations), then redirect FA scratch to the
			// resident per-thread buffers.
			tl_fa[t] = *FA;
			// NOTE: this snapshot rewinds every scalar in tl_fa[t] to the master
			// FA's value once per generation. Nothing that must stay monotonic
			// against a RESIDENT buffer may live in FA_Global — the contacts
			// epoch lives inside tl_contacts[t] itself for exactly this reason
			// (flexaid.h, tests/test_contacts_epoch.cpp).
			tl_fa[t].contacts      = tl_contacts[t].data();
			tl_fa[t].contributions = tl_contrib[t].data();
			tl_fa[t].optres        = tl_optres[t].data();
			// Keep optres non-cf fields in sync with the reference (cf fields are
			// cleared per-chromosome in the eval loop below).
			std::copy(FA->optres, FA->optres + nopt, tl_optres[t].begin());
			// Refresh the cheap live VC snapshot, then redirect VC scratch.
			tl_vc[t] = *VC;
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
			wire_scorable(tl_vc[t], tl_scorable[t].data(), natmr);
			// Keep the reference-calculation retry path enabled in GA workers.
			// The direct native probe uses recalc=1; forcing 0 here caused the
			// same pose to fall into the non-convergence penalty path.
			tl_vc[t].recalc    = 1;
			// box is shared: if vindex==1 it's pre-built read-only;
			// if vindex==0 Vcontacts will malloc/vcfunction will free per call.
		}

		// A4a: precompute list of atoms that have optres pointers so the
		// per-chromosome redirect loop skips the >90% of atoms with optres==NULL.
		struct OptresAtomEntry { int ai; ptrdiff_t oidx; };
		std::vector<OptresAtomEntry> optres_atom_list;
		optres_atom_list.reserve(nopt * 4);
		for(int ai = 1; ai <= natm; ++ai) {
			if(atoms[ai].optres)
				optres_atom_list.push_back({ai, atoms[ai].optres - FA->optres});
		}
		const int n_optres_atoms = (int)optres_atom_list.size();

		// ── FLEXAID_DETERMINISTIC (P3 · CI A/B reproducibility) ──────────────
		// When set, the parallel CF-eval loop runs on a single thread so its
		// reduction order is serial-equivalent and bit-reproducible run-to-run
		// (each chromosome writes only its own chrom[ii].cf, so 1-thread ==
		// deterministic order). Unset (default) uses the fast multi-thread path,
		// where the ~0.2% chromosome-level numeric drift documented in
		// OPTIMIZATION_KNOWN_ISSUES.md is accepted (see FLEXAID_FAST_DOCKING_PLAN
		// §P3/§4.3). Enabled by the compile-time macro -DFLEXAID_DETERMINISTIC
		// OR by the env var FLEXAID_DETERMINISTIC=<non-empty, not "0">.
		static const bool deterministic_eval = [](){
#ifdef FLEXAID_DETERMINISTIC
			return true;
#else
			const char* env = std::getenv("FLEXAID_DETERMINISTIC");
			return env && env[0] != '\0' && env[0] != '0';
#endif
		}();
#ifdef _OPENMP
		const int eval_threads = deterministic_eval ? 1 : n_thr;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(eval_threads) default(none) \
	shared(chrom, pop_size, GB, gene_lim, cleftgrid, target, \
	       atoms, residue, FA, VC, eval_threads, \
	       tl_atoms, tl_res, tl_fa, tl_optres, tl_vc, \
	       natm, nres, nopt, n_receptor_chains, \
	       use_selective, dirty_atm, dirty_res_idx, n_dirty_atm, n_dirty_res, \
	       optres_atom_list, n_optres_atoms)
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
			// A4a: redirect optres pointers using precomputed index list
			// (skips atoms with optres==NULL — typically >90% of atoms).
			for (int oa = 0; oa < n_optres_atoms; ++oa) {
				const auto& e = optres_atom_list[oa];
				tl_atoms[tid][e.ai].optres = &tl_optres[tid][e.oidx];
			}
			// optres cf fields are cleared by vcfunction itself; pre-clear for safety.
			for (int o = 0; o < nopt; ++o) {
				tl_optres[tid][o].cf.com         = 0.0;
				tl_optres[tid][o].cf.wal         = 0.0;
				tl_optres[tid][o].cf.sas         = 0.0;
				tl_optres[tid][o].cf.totsas      = 0.0;
				tl_optres[tid][o].cf.con         = 0.0;
				tl_optres[tid][o].cf.gist        = 0.0;
				tl_optres[tid][o].cf.elec        = 0.0;
				tl_optres[tid][o].cf.hbond       = 0.0;
				tl_optres[tid][o].cf.metal_coord = 0.0;
				tl_optres[tid][o].cf.gist_desolv = 0.0;
				tl_optres[tid][o].cf.pb_clash    = 0.0;
				tl_optres[tid][o].cf.rclash      = 0;
			}
			tl_vc[tid].numcarec = 0;

			// Load this chromosome's ring pucker phases into the per-thread FA
			// so ic2cf reconstructs its puckered ring (no-op when inactive).
			ring_load_chrom_to_fa(&tl_fa[tid], &chrom[ii]);

			chrom[ii].cf = eval_chromosome(
			    &tl_fa[tid], GB, &tl_vc[tid], gene_lim,
			    tl_atoms[tid].data(), tl_res[tid].data(),
			    cleftgrid, chrom[ii].genes, target);
			chrom[ii].evalue     = get_cf_evalue(&chrom[ii].cf, FA) / n_receptor_chains;
			chrom[ii].app_evalue = get_apparent_cf_evalue(&chrom[ii].cf) / n_receptor_chains;
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
		// G4.2: optional Cartesian ligand RMSD niche (precompute coords once).
		// Distance math: flexaids::niche_* in niche_distance.h (unit-tested).
		const bool niche_cart = flexaids::niche_cartesian_env_enabled();
		constexpr int kCoordStride = MAX_ATM_HET * 3;
		std::vector<float> lig_xyz;
		int n_lig_atoms = 0;
		if (niche_cart) {
			lig_xyz.assign(static_cast<size_t>(GB->num_chrom) * kCoordStride, 0.0f);
			for (int c = 0; c < GB->num_chrom; ++c) {
				calc_rmsd_chrom(FA, GB, chrom, gene_lim, atoms, residue, cleftgrid,
				                GB->num_genes, c, c,
				                &lig_xyz[static_cast<size_t>(c) * kCoordStride],
				                nullptr, false);
			}
			const int lres = atoms[FA->map_par[0].atm].ofres;
			const int rot = residue[lres].rot;
			n_lig_atoms = residue[lres].latm[rot] - residue[lres].fatm[rot] + 1;
			if (n_lig_atoms < 1) n_lig_atoms = 1;
		}
		std::vector<double> pshare_out(static_cast<size_t>(GB->num_chrom), 1.0);
		const bool niche_hash = niche_cart && flexaids::niche_hash_enabled();
		std::vector<float> niche_cents;
		std::unordered_map<flexaids::NicheCell, std::vector<int>, flexaids::NicheCellHash> niche_map;
		if (niche_hash) {
			niche_cents.assign(static_cast<size_t>(GB->num_chrom) * 3, 0.f);
			for (int c = 0; c < GB->num_chrom; ++c) {
				const float* xyz = &lig_xyz[static_cast<size_t>(c) * kCoordStride];
				float sx = 0.f, sy = 0.f, sz = 0.f;
				for (int a = 0; a < n_lig_atoms; ++a) {
					sx += xyz[a * 3 + 0];
					sy += xyz[a * 3 + 1];
					sz += xyz[a * 3 + 2];
				}
				const float inv = n_lig_atoms > 0 ? 1.f / static_cast<float>(n_lig_atoms) : 0.f;
				niche_cents[static_cast<size_t>(c) * 3 + 0] = sx * inv;
				niche_cents[static_cast<size_t>(c) * 3 + 1] = sy * inv;
				niche_cents[static_cast<size_t>(c) * 3 + 2] = sz * inv;
			}
			niche_map = flexaids::niche_hash_build(
				niche_cents.data(), GB->num_chrom, static_cast<float>(GB->sig_share));
		}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(chrom, GB, FA, cleftgrid, niche_cart, lig_xyz, n_lig_atoms, pshare_out, \
	       niche_hash, niche_cents, niche_map)
#endif
		for(int pi=0; pi<GB->num_chrom; pi++){
			double pshare = 0.0;
			std::vector<int> neigh;
			if (niche_hash) {
				const float* c = &niche_cents[static_cast<size_t>(pi) * 3];
				flexaids::niche_hash_neighbors(
					niche_map,
					flexaids::niche_cell_of(c[0], c[1], c[2],
					                       static_cast<float>(GB->sig_share)),
					neigh);
			}
			const int n_j = niche_hash ? static_cast<int>(neigh.size()) : GB->num_chrom;
			for(int jk=0; jk<n_j; jk++){
				const int pj = niche_hash ? neigh[static_cast<size_t>(jk)] : jk;
				double prmsp = 0.0;
				if (niche_cart) {
					const float* a = &lig_xyz[static_cast<size_t>(pi) * kCoordStride];
					const float* b = &lig_xyz[static_cast<size_t>(pj) * kCoordStride];
					prmsp = flexaids::niche_cartesian_rmsd(a, b, n_lig_atoms);
				} else {
					prmsp = calc_rmsp(GB->num_genes,
					                 chrom[pi].genes, chrom[pj].genes,
					                 FA->map_par, cleftgrid,
					                 GB->sig_share * GB->sig_share);  // A4b early exit
				}
				if(prmsp <= GB->sig_share){
					pshare += (1.0 - pow((prmsp/GB->sig_share), GB->alpha));
				}
			}
			pshare_out[static_cast<size_t>(pi)] = pshare;
			// Assign fitness AFTER accumulating the full niche count.
			// v27 elitism: the top n_elite (lowest evalue → smallest pi after the
			// ascending QuickSort above) are exempt from the sharing reduction so
			// niching can never demote the running best out of the selection pool.
			if (pi < GB->n_elite)
				chrom[pi].fitnes = (double)(GB->num_chrom - pi);
			else
				chrom[pi].fitnes = (double)(GB->num_chrom - pi) / pshare;
		}
		if (niche_cart && (gen_id % 50 == 0)) {
			double sum_ps = 0.0;
			int n_lonely = 0;
			for (int pi = 0; pi < GB->num_chrom; ++pi) {
				const double ps = pshare_out[static_cast<size_t>(pi)];
				sum_ps += ps;
				if (ps <= 1.0 + 1e-9) ++n_lonely;
			}
			const double mean_ps = sum_ps / std::max(1, GB->num_chrom);
			// n_niches proxy: chromosomes that are alone in their niche (pshare≈1)
			// plus a soft count of crowded niches via mean share.
			fprintf(stderr,
			        "[NICHE-CART] gen=%d n_lonely=%d/%d mean_pshare=%.3f sigma=%.3fA n_lig=%d\n",
			        gen_id, n_lonely, GB->num_chrom, mean_ps, GB->sig_share, n_lig_atoms);
		}
	}

	if(strcmp(method,"SMFREE")==0){
		/* SMFREE — soft-β CF sampling with niche sharing (ensemble layer 3).
		   Selection uses β_sel = 1/T (same as ACF clustering), NOT physical
		   1/(kB·T). Niche share is gene-space calc_rmsp unless FLEXAIDDS_NICHE_CARTESIAN=1
		   (G4.2 Cartesian ligand heavy-atom RMSD). The kB-based CF transform is
		   proxy-only and diagnostic. Reproducibility: same β_sel as election.
		*/
		if (FA->temperature > 0) {
			const double T = static_cast<double>(FA->temperature);
			double beta_sel = 0.0;
			(void)ensemble::soft_selection_beta(T, &beta_sel);
			statmech::StatMechEngine engine(
				T, statmech::make_contact_function_optimizer_provenance());

			// Feed all chromosome energies into the engine.
			for (int si = 0; si < GB->num_chrom; si++) {
				engine.add_sample(chrom[si].evalue);
			}

			// Compute the legacy CF ensemble transform (β_num = 1/kBT) and
			// SELECTION weights (β_sel = 1/T, matching the clustering
			// convention FA->beta). Using β_num here would collapse
			// selection to a zero-temperature argmax (e^{βΔCF} with β≈1.68),
			// killing the thermal diversity SMFREE is meant to inject. See P1.
			auto thermo = engine.compute();
			auto bweights = engine.selection_weights();

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

			const bool niche_cart = flexaids::niche_cartesian_env_enabled();
			constexpr int kCoordStride = MAX_ATM_HET * 3;
			std::vector<float> lig_xyz;
			int n_lig_atoms = 0;
			if (niche_cart) {
				lig_xyz.assign(static_cast<size_t>(GB->num_chrom) * kCoordStride, 0.0f);
				for (int c = 0; c < GB->num_chrom; ++c) {
					calc_rmsd_chrom(FA, GB, chrom, gene_lim, atoms, residue, cleftgrid,
					                GB->num_genes, c, c,
					                &lig_xyz[static_cast<size_t>(c) * kCoordStride],
					                nullptr, false);
				}
				const int lres = atoms[FA->map_par[0].atm].ofres;
				const int rot = residue[lres].rot;
				n_lig_atoms = residue[lres].latm[rot] - residue[lres].fatm[rot] + 1;
				if (n_lig_atoms < 1) n_lig_atoms = 1;
			}
			std::vector<double> pshare_out(static_cast<size_t>(GB->num_chrom), 1.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) default(none) \
	shared(chrom, GB, FA, cleftgrid, max_bw, w, niche_cart, lig_xyz, n_lig_atoms, pshare_out)
#endif
			for (int pi = 0; pi < GB->num_chrom; pi++) {
				// Niche sharing via flexaids::niche_* (gene RMSP or Cartesian RMSD).
				double pshare = 0.0;
				for (int pj = 0; pj < GB->num_chrom; pj++) {
					double prmsp = 0.0;
					if (niche_cart) {
						const float* a = &lig_xyz[static_cast<size_t>(pi) * kCoordStride];
						const float* b = &lig_xyz[static_cast<size_t>(pj) * kCoordStride];
						prmsp = flexaids::niche_cartesian_rmsd(a, b, n_lig_atoms);
					} else {
						prmsp = calc_rmsp(GB->num_genes,
						                 chrom[pi].genes, chrom[pj].genes,
						                 FA->map_par, cleftgrid,
						                 GB->sig_share * GB->sig_share);  // A4b early exit
					}
					if (prmsp <= GB->sig_share) {
						pshare += (1.0 - pow((prmsp / GB->sig_share), GB->alpha));
					}
				}
				pshare_out[static_cast<size_t>(pi)] = pshare;

				// Rank component: normalised to [0, 1].
				double rank_component = static_cast<double>(GB->num_chrom - pi) /
				                        static_cast<double>(GB->num_chrom);

				// Boltzmann component: normalised to [0, 1] by max weight.
				double boltz_component = chrom[pi].boltzmann_weight / max_bw;

				// Blended fitness divided by niche count.
				// v27 elitism: top n_elite (lowest evalue → smallest pi after the
				// ascending QuickSort) are exempt from the sharing reduction.
				double blended = (1.0 - w) * rank_component + w * boltz_component;
				if (pi < GB->n_elite)
					chrom[pi].fitnes = blended * static_cast<double>(GB->num_chrom);
				else
					chrom[pi].fitnes = blended * static_cast<double>(GB->num_chrom) / pshare;
			}

			if (niche_cart && (gen_id % 50 == 0)) {
				double sum_ps = 0.0;
				int n_lonely = 0;
				for (int pi = 0; pi < GB->num_chrom; ++pi) {
					const double ps = pshare_out[static_cast<size_t>(pi)];
					sum_ps += ps;
					if (ps <= 1.0 + 1e-9) ++n_lonely;
				}
				fprintf(stderr,
				        "[NICHE-CART] gen=%d n_lonely=%d/%d mean_pshare=%.3f sigma=%.3fA n_lig=%d (SMFREE)\n",
				        gen_id, n_lonely, GB->num_chrom,
				        sum_ps / std::max(1, GB->num_chrom), GB->sig_share, n_lig_atoms);
			}

			// Log selection β + thermo periodically (greppable reproducibility audit).
			if (gen_id % GA_SMFREE_LOG_INTERVAL == 0) {
				fprintf(stderr,
					"[SMFREE] gen=%d  beta_sel=%.6f  T=%.1f  F=%.3f  <E>=%.3f  "
					"S=%.6f  Cv=%.4f  σ_E=%.3f\n",
					gen_id, beta_sel, T, thermo.free_energy, thermo.mean_energy,
					thermo.entropy, thermo.heat_capacity, thermo.std_energy);
			}
		} else {
			// Temperature = 0: rank-only. Loud warn — product claims entropy but
			// sampling is not soft-β (non-reproducible vs T>0 runs).
			static int smfree_t0_warned = 0;
			if (!smfree_t0_warned) {
				smfree_t0_warned = 1;
				fprintf(stderr,
					"[SMFREE] WARN: temperature=0 → rank-only fitness "
					"(no soft-β sampling). Set thermodynamics.temperature>0 "
					"for reproducible ensemble layer 3 (β=1/T).\n");
				if (flexaids::ProtocolConfig::from_env().smfree_require_t) {
					fprintf(stderr,
						"[SMFREE] FATAL: FLEXAIDDS_SMFREE_REQUIRE_T set and T=0\n");
					Terminate(2);
				}
			}
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

	// Safety net: catch any Terminate() thrown deep inside the scoring pipeline
	// (e.g., residual Vcontacts bounding-box checks) and return maximum penalty
	// so the chromosome survives non-competitively rather than aborting the run.
	try {
		return (*function)(FA,VC,atoms,residue,cleftgrid,GB->num_genes,icv);
	} catch (const FlexAIDException&) {
		cfstr cf_penalty{};
		cf_penalty.com = 99999.0;
		return cf_penalty;
	}
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

	// Multi-chain VCT normalisation (see GA() comment for rationale)
	const int n_receptor_chains = count_receptor_chains(FA, residue);

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

		// ── Coarse-init pocket scan (autonomous / blinded mode) ───────────
		// Run once per populate_chromosomes() call (coarse_seeds_count==0 guard
		// prevents re-runs on boom/SEC repopulations after the first gen).
		if (FA->coarse_init_enabled && FA->coarse_seeds_count == 0 &&
		    FA->num_grd > 0 && cleftgrid != nullptr) {
			run_coarse_pocket_scan(FA, VC, GB, atoms, residue, cleftgrid,
			                       gene_lim, dice);
		}
		// Inject pre-screened coarse seeds at the front of the population.
		int coarse_offset = 0;
		if (FA->coarse_init_enabled && FA->coarse_seeds_count > 0) {
			coarse_offset = std::min(FA->coarse_seeds_count,
			                         GB->num_chrom - popoffset);
			for (int si = 0; si < coarse_offset; si++) {
				const int ci = popoffset + si;
				const double grid_ic = static_cast<double>(FA->coarse_seeds_grid[si]);
				chrom[ci].genes[0].to_ic    = grid_ic;
				chrom[ci].genes[0].to_int32 = ictogene(&gene_lim[0], grid_ic);
				for (int g = 1; g < GB->num_genes; g++) {
					double ic = static_cast<double>(
					    FA->coarse_seeds_genes[si * (GB->num_genes - 1) + (g - 1)]);
					chrom[ci].genes[g].to_ic    = ic;
					chrom[ci].genes[g].to_int32 = ictogene(&gene_lim[g], ic);
				}
				// Preserve the exact genes that were scored by coarse_init. Ring
				// randomisation here changed a screened pose after its score was
				// accepted, breaking the seed-to-gen-0 correspondence.
				const size_t csig = hash_genes(chrom[ci].genes, GB->num_genes);
				duplicates[csig] = 1;
			}
			printf("[COARSE-INIT] Injected %d pre-screened seeds into gen-0\n",
			       coarse_offset);
		}
		// ── End coarse-init injection ─────────────────────────────────────

		i = popoffset + coarse_offset;
		while(i<GB->num_chrom){
			while(1){
				generate_random_individual(FA,GB,atoms,chrom[i].genes,gene_lim,dice,0,GB->num_genes);

				// ── MIF-weighted or RefLig seeding override for gene 0 ──
				// Oracle-ceiling / re-dock: when pose_seed_enabled and seed_fraction>0,
				// inject the crystal IC (gene0=0 + opt_par orientations/torsions) into
				// the seeded fraction.  Historically this worked when reflig_file was
				// empty (native_direct_seed).  Modern DatasetRunner always sets
				// reference_ligand.file to the crystal SDF for RMSD/PB; that must NOT
				// disable orientation seeding or the seeded fraction stays random
				// (Astex 1GPK: gen-0 CF≈+32 with no native, hist CF≈−74 with native flood).
				const int seed_budget = static_cast<int>(
				    FA->reflig_seed_fraction *
				    static_cast<float>(GB->num_chrom - popoffset));
				const bool want_pose_seed =
				    FA->reflig_pose_seed_enabled &&
				    FA->opt_par != nullptr &&
				    FA->map_par != nullptr &&
				    seed_budget > 0;
				const bool reflig_seeded =
				    (FA->reflig_nearest_count > 0 || want_pose_seed) &&
				    i < popoffset + seed_budget;
				if (reflig_seeded) {
					// Prefer native grid anchor (index 0) whenever the IC frame
					// allows it and pose seeding is on.  Nearest-grid cycling is a
					// fallback only when gene0 cannot be 0 (or pose seed is off).
					const bool native_direct_seed =
					    FA->resligand != NULL && gene_lim[0].min <= 0.0 &&
					    (FA->reflig_file[0] == '\0' || FA->reflig_pose_seed_enabled);
					int grid_idx = 0;
					if (!native_direct_seed && FA->reflig_nearest_count > 0 &&
					    FA->reflig_nearest_grid) {
						// Explicit RefLig grid bias without native pose seed.
						int k = (i - popoffset) % FA->reflig_nearest_count;
						grid_idx = FA->reflig_nearest_grid[k];
					}
					chrom[i].genes[0].to_ic = static_cast<double>(grid_idx);
					chrom[i].genes[0].to_int32 = ictogene(&gene_lim[0],
					                                       static_cast<double>(grid_idx));
					if (want_pose_seed) {
						for (int g = 1; g < GB->num_genes; g++) {
							// typ: -1 translation (gene0 only), 1 angle, 2 dihedral,
							// 3 normal mode (skip), 4 other special.
							const int typ = FA->map_par[g].typ;
							if (typ == 3) continue;
							double ref_ic = FA->opt_par[g];
							// Tiny per-chromosome jitter keeps near-native diversity
							// when GB->duplicates is off, without leaving the basin.
							if (i > popoffset && (typ == 1 || typ == 2)) {
								const double jitter =
								    (static_cast<int>(dice() % 5) - 2) * 0.25; // −0.5..0.5°
								ref_ic += jitter;
							}
							chrom[i].genes[g].to_ic = ref_ic;
							chrom[i].genes[g].to_int32 =
							    ictogene(&gene_lim[g], ref_ic);
						}
					}
					if (i == popoffset) {
						printf("[REFLIG-SEED] native_anchor=%d pose_seed=%d frac=%.2f "
						       "budget=%d nearest=%d ngenes=%d first=(",
						       native_direct_seed ? 1 : 0,
						       FA->reflig_pose_seed_enabled,
						       FA->reflig_seed_fraction, seed_budget,
						       FA->reflig_nearest_count, GB->num_genes);
						for (int g = 0; g < GB->num_genes && g < 8; g++)
							printf("%s%.3f", g ? "," : "",
							       chrom[i].genes[g].to_ic);
						printf(")\n");
					}
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
				if(reflig_seeded || GB->duplicates || duplicates.find(sig) == duplicates.end()){
					break;
				}
			}

			// Gen-0 ring pucker randomisation (no-op unless ring flex is active).
			ring_randomise_chrom(FA, &chrom[i]);

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
		// CONTACTS_BUFFER_SIZE: stamps + the trailing epoch slot (flexaid.h).
		std::vector<std::vector<int>>    p_contacts(n_thr, std::vector<int>(CONTACTS_BUFFER_SIZE, 0));
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
		std::vector<std::vector<int>>        p_scorable(n_thr, std::vector<int>(natmr, 0));

		for (int t = 0; t < n_thr; ++t) {
			p_fa[t].contacts      = p_contacts[t].data();
			p_fa[t].contributions = p_contrib[t].data();
			p_fa[t].optres        = p_optres[t].data();
			p_vc[t].Calc      = p_calc[t].data();
			p_vc[t].Calclist  = p_calclist[t].data();
			p_vc[t].ca_index  = p_caidx[t].data();
			p_vc[t].scorable_list = p_scorable[t].data();
			p_vc[t].n_scorable = 0;
			p_vc[t].scorable_cap = natmr;
			p_vc[t].fastpath_used = 0;
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
				for (int m = 0; m < FA->nmov[r]; ++m) {
					const int ai = FA->mov[r][m];
					if (ai >= 1 && ai <= natm) p_dirty_atm.push_back(ai);
				}
			for (int p = 0; p < FA->npar; ++p) {
				const int ai = FA->map_par[p].atm;
				if (ai >= 1 && ai <= natm) p_dirty_atm.push_back(ai);
			}
			for (int p = 0; p < FA->npar; ++p) {
				if (FA->map_par[p].typ == 2 &&
				    FA->map_par[p].atm >= 1 && FA->map_par[p].atm <= natm) {
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
				if (FA->map_par[p].typ == 4 &&
				    FA->map_par[p].atm >= 1 && FA->map_par[p].atm <= natm)
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
	       n_receptor_chains, p_use_selective, p_dirty_atm, p_dirty_res_idx, \
	       p_n_dirty_atm, p_n_dirty_res)
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
				p_optres[tid][o].cf.pb_clash = 0.0;
				p_optres[tid][o].cf.rclash = 0;
			}
			p_vc[tid].numcarec = 0;

			ring_load_chrom_to_fa(&p_fa[tid], &chrom[i]);

			chrom[i].cf = eval_chromosome(
			    &p_fa[tid], GB, &p_vc[tid], gene_lim,
			    p_atoms[tid].data(), p_res[tid].data(),
			    cleftgrid, chrom[i].genes, target);
			chrom[i].evalue     = get_cf_evalue(&chrom[i].cf, FA) / n_receptor_chains;
			chrom[i].app_evalue = get_apparent_cf_evalue(&chrom[i].cf) / n_receptor_chains;
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
		(void)temp; // reserved for swap in crossover extensions
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
	mutate(john, num_genes, mut_rate, /*gene_lim=*/nullptr);
}

void mutate(gene *john,int num_genes,double mut_rate,const genlim* gene_lim){
	/* creates an operator with 1's with rate= mut_rate
	   uses it to mutate john.

	   Default (gene_lim null or FLEXAIDDS_MUTATION_GRANULAR unset): classic
	   per-bit flip across all 32 bits.  Low-order bits are often "dead" for
	   decoding: genetoic bins by gene/2^31, so flips below ~2^31/nbin do not
	   change the IC phenotype (PHASE4 G4.3).

	   G4.3 mode (FLEXAIDDS_MUTATION_GRANULAR=1 and gene_lim provided): when a
	   gene is selected for mutation, apply a ±1-bin step in gene integer
	   space (phenotype-changing small move). L4: one-shot [MUT-GRAN] on stderr.

	   S4-A FLEXAIDDS_PHENOTYPE_UNIQUE=1 (classic path only): if classic bit-flips
	   leave phenotype bins unchanged, force one ±1-bin phenotype step.
	*/
	// Env is process-lifetime for dock runs. Tests setenv before first mutate.
	// Skip getenv when gene_lim is null (classic overload path).
	bool env_on = false;
	if (gene_lim != nullptr) {
		const char* e = std::getenv("FLEXAIDDS_MUTATION_GRANULAR");
		env_on = e && (e[0] == '1' || e[0] == 'y' || e[0] == 'Y' || e[0] == 't' ||
		               e[0] == 'T');
	}
	const bool use_granular = env_on && (gene_lim != nullptr);
	if (use_granular) {
		static bool s_logged = false;
		if (!s_logged) {
			s_logged = true;
			std::fprintf(stderr,
			             "[MUT-GRAN] FLEXAIDDS_MUTATION_GRANULAR=1: bin-aware ±1-bin gene steps "
			             "(phenotype-live mutations; classic dead low bits avoided)\n");
		}
	}

	const bool force_pheno =
	    (gene_lim != nullptr) && !use_granular &&
	    flexaids::new_search::phenotype_unique_enabled();
	std::size_t h_before = 0;
	if (force_pheno) {
		static bool s_pheno_logged = false;
		if (!s_pheno_logged) {
			s_pheno_logged = true;
			std::fprintf(stderr,
			             "[NEW-SEARCH-ARCH] phenotype_unique=1: classic mutate forces "
			             "±1-bin phenotype step when bit-flips leave bins unchanged\n");
		}
		for (int j = 0; j < num_genes; ++j)
			john[j].to_ic = genetoic(&gene_lim[j], john[j].to_int32);
		h_before = flexaids::new_search::hash_phenotype_bins(john, num_genes, gene_lim);
	}

	for (int j = 0; j < num_genes; j++) {
		if (use_granular) {
			// Gate once per gene at mut_rate (not per bit).
			if (RandomDouble() >= mut_rate) continue;
			const double nbin = gene_lim[j].nbin > 1.0 ? gene_lim[j].nbin : 2.0;
			const int32_t step = static_cast<int32_t>(std::max(
			    1.0, std::floor((static_cast<double>(MAX_RANDOM_VALUE) + 1.0) / nbin)));
			// ±1 bin; occasional ±2 for slightly larger local moves (still small).
			int k = 1;
			if (RandomDouble() < 0.25) k = 2;
			const int sign = (RandomDouble() < 0.5) ? 1 : -1;
			int64_t ng = static_cast<int64_t>(john[j].to_int32) +
			             static_cast<int64_t>(sign) * static_cast<int64_t>(step) * k;
			if (ng < 0) ng = 0;
			if (ng > static_cast<int64_t>(MAX_RANDOM_VALUE))
				ng = static_cast<int64_t>(MAX_RANDOM_VALUE);
			john[j].to_int32 = static_cast<int32_t>(ng);
			continue;
		}

		unsigned int optr = 0u;
		unsigned int test = 1u;
		for (int i = 0; i < 32; i++) {
			if (RandomDouble() < mut_rate) {
				optr |= test;
			}
			test <<= 1;
		}
		john[j].to_int32 ^= static_cast<int32_t>(optr);
	}

	if (force_pheno) {
		for (int j = 0; j < num_genes; ++j)
			john[j].to_ic = genetoic(&gene_lim[j], john[j].to_int32);
		const std::size_t h_after =
		    flexaids::new_search::hash_phenotype_bins(john, num_genes, gene_lim);
		if (h_after == h_before && num_genes > 0) {
			const int j = static_cast<int>(RandomDouble() * static_cast<double>(num_genes)) %
			              num_genes;
			const double nbin = gene_lim[j].nbin > 1.0 ? gene_lim[j].nbin : 2.0;
			const int sign = (RandomDouble() < 0.5) ? 1 : -1;
			const int k = (RandomDouble() < 0.25) ? 2 : 1;
			flexaids::new_search::apply_phenotype_bin_step(
			    &john[j], nbin, sign, k, MAX_RANDOM_VALUE);
			john[j].to_ic = genetoic(&gene_lim[j], john[j].to_int32);
		}
	}
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

/* A4b: early-exit RMSP — pass early_exit_sq > 0 to short-circuit once the
   partial sum of squared differences exceeds early_exit_sq * npar
   (i.e., partial RMSP already exceeds the threshold).
   Default 0.0 = no early exit (full computation, same as before). */
double calc_rmsp(int npar, const gene* g1, const gene* g2, const optmap* map_par,
                 gridpoint* cleftgrid, double early_exit_sq)
{
	double sum_sq = 0.0;
	const double threshold = early_exit_sq * (double)npar;  // 0 → no early exit
	for (int ii = 0; ii < npar; ++ii) {
		const double d = g1[ii].to_ic - g2[ii].to_ic;
		sum_sq += d * d;
		if (threshold > 0.0 && sum_sq > threshold)
			return std::sqrt(sum_sq / (double)npar);  // guaranteed > early_exit_sq
	}
	return std::sqrt(sum_sq / (double)npar);
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

	// genetoic() decodes a gene by counting bins UP from `bin` and returns
	// index i_dec = ceil(frac/bin) - 1.  The previous encoding here used a
	// DECREASING index (tot = 1 - i*bin), which made genetoic(ictogene(ic))
	// reflect to i_dec = nbin - i instead of round-tripping.  Encode the
	// fractional position at the CENTER of bin i so the decoder recovers the
	// same index: frac = (i + 0.5)*bin  =>  i_dec == i.  Same gene range and
	// bin count — only the encode/decode correspondence is fixed.
	double tot = ((double)i + 0.5) * gene_lim->bin;

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
	// Thread-safe RNG tied to ga.seed / FLEXAID_SEED via lazy_thread_rng.
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
}
