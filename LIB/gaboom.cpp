#include "gaboom.h"
#include "Vcontacts.h"
#include "fileio.h"
#include "flexaid_exception.h"
#include "ga_constants.h"
#include "UnifiedHardwareDispatch.h"
#include "MIFGrid.h"
#include "CavityDetect/SpatialGrid.h"
#include "RngSeed.h"

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

	// ── OMP thread default: 2/worker if OMP_NUM_THREADS not set in environment.
	// Leaves headroom for Metal dispatch and OS scheduling on M3 Pro 11 P-cores.
#ifdef _OPENMP
	if (!std::getenv("OMP_NUM_THREADS")) {
		omp_set_num_threads(2);
	}
#endif

	int i;
	int print=0;

	// ── Level-3 H(ω) diagnostic env override ──────────────────────────────────
	// FLEXAIDDS_USE_SHANNON=1 enables the ligand vibrational-mode Shannon-entropy
	// monitor (per-generation [HVIB] lines) even on the bare binary, mirroring the
	// DatasetRunner config toggle.  Engine-side env wins, matching FLEXAIDDS_N_ELITE.
	// Purely diagnostic — does NOT enter CF or fitness.  Default stays OFF.
	if (const char* _hv = std::getenv("FLEXAIDDS_USE_SHANNON")) {
		if (_hv[0] != '\0' && _hv[0] != '0') {
			if (!GB->use_shannon)
				fprintf(stderr, "[HVIB] FLEXAIDDS_USE_SHANNON set: enabling ligand "
				        "vibrational-entropy H(ω) monitor (diagnostic only)\n");
			GB->use_shannon = 1;
		}
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

	// Defensive clamp: a zero (or negative) check interval reaches three
	// integer division / modulo sites in the generation loop below and would
	// raise SIGFPE. Guarantee a sane value no matter how it was configured.
	if (GB->entropy_check_interval <= 0) {
		GB->entropy_check_interval = GA_DEFAULT_ENTROPY_CHECK_INTERVAL;
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
	if (const char* _ii = std::getenv("FLEXAIDDS_INSTREAM_INTERVAL")) {
		int v = std::atoi(_ii);
		if (v >= 1) {
			instream_interval = v;
			fprintf(stderr, "[INSTREAM] merge/H(ω) cadence overridden to every %d "
			        "generation(s) via FLEXAIDDS_INSTREAM_INTERVAL\n", instream_interval);
		}
	}

	////// Genetic Algorithm ///////
	////////////////////////////////
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

	// ── Temperature annealing (FLEXAIDDS_T_HOT) ──────────────────────────────
	// Exponential decay T_hot → 298 K: T(α) = T_hot·exp(−5α) + 298·(1−exp(−5α))
	// where α = gen/(max_gen−1) ∈ [0,1].  Affects SMFREE Boltzmann-weight
	// selection only; post-GA thermodynamics use the final temperature (≈298 K).
	// arm3b ablation (5000K constant, Fable 5) was net-neutral in oracle mode
	// → native basin gravitationally dominant.  Annealing targets near-miss
	// false-minimum escape early in the run while native seeds lock in.
	// Useful calibration range: 500–2000 K.
	const double t_hot_anneal = []() -> double {
		const char* env = std::getenv("FLEXAIDDS_T_HOT");
		return (env && env[0] != '\0') ? std::atof(env) : 0.0;
	}();
	const bool do_anneal = (t_hot_anneal > 298.0) && (FA->temperature > 0);
	if (do_anneal) {
		fprintf(stderr, "[ANNEAL] Temperature annealing enabled: "
		        "T_hot=%.0f K → 298 K over %d generations (exp-5 schedule)\n",
		        t_hot_anneal, GB->max_generations);
		// Prime initial temperature so gen-0 SMFREE sees T_hot
		FA->temperature = static_cast<unsigned int>(std::round(t_hot_anneal));
		FA->beta        = 1.0 / t_hot_anneal;
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

		// ── Always-on H plateau early exit ─────────────────────────────────
		// Every entropy_check_interval generations, sample H of the current
		// population and push into a 20-slot ring buffer.  If the absolute
		// difference between the newest and oldest slot < kHPlateauEps nats,
		// the distribution has stopped collapsing → write best pose and stop.
		if (!entropy_converged && !ga_stagnant &&
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
		if (GB->entropy_convergence &&
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
					(*chrom)[ci].evalue     = get_cf_evalue(&(*chrom)[ci].cf);
					(*chrom)[ci].app_evalue = get_apparent_cf_evalue(&(*chrom)[ci].cf);
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
			                   + 298.0 * (1.0 - std::exp(-5.0 * alpha));
			FA->temperature = static_cast<unsigned int>(std::round(T_now));
			FA->beta        = 1.0 / T_now;
			if (i % 200 == 0) {
				fprintf(stderr, "[ANNEAL] gen=%4d  T=%7.1f K  α=%.4f\n",
				        i + 1, T_now, alpha);
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

		save_snapshot(&(*chrom_snapshot)[i*GB->num_chrom],(*chrom),save_num_chrom,GB->num_genes);
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

		// ── Enthalpy-Entropy Index (Williams et al. 2017, Drug Discov. Today) ──
		// I_EE = (ΔH + T·ΔS) / ΔG   — diagnostic only, never for ranking
		{
			const statmech::ThermodynamicBreakdown bd = sme.compute_breakdown();
			if (bd.has_I_EE) {
				printf("  Enthalpy-Entropy Idx  I_EE= %10.4f  [Williams 2017]\n", bd.I_EE);
				const char* regime =
					(bd.I_EE > 1.05)  ? "entropy-assisted" :
					(bd.I_EE < 0.95 && bd.I_EE >= 0.0) ? "entropy-opposed" :
					(bd.I_EE < 0.0)   ? "entropy-driven (rare)" :
					                    "pure enthalpy";
				printf("                        (%s)\n", regime);
			}
		}

		// ── Kirchhoff ΔG(T) extrapolation (Robertson & Murphy 1997) ────────
		// Activated only when DSF/TSA Tm has been supplied via dsf_Tm_K.
		// ΔG(T) = ΔHm(1 − T/Tm) − ΔCp[(Tm − T) + T·ln(T/Tm)]
		if (FA->dsf_Tm_K > 0.0) {
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

			if (FA->ring_flex_active) {
				memcpy(chrom[GB->num_chrom+i].ring_phases, rc1.ring_phases, sizeof(rc1.ring_phases));
				memcpy(chrom[GB->num_chrom+i].ring_six,    rc1.ring_six,    sizeof(rc1.ring_six));
				memcpy(chrom[GB->num_chrom+i].ring_five,   rc1.ring_five,   sizeof(rc1.ring_five));
				ring_load_chrom_to_fa(FA, &chrom[GB->num_chrom+i]);
			}

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

			if (FA->ring_flex_active) {
				memcpy(chrom[GB->num_chrom+i].ring_phases, rc2.ring_phases, sizeof(rc2.ring_phases));
				memcpy(chrom[GB->num_chrom+i].ring_six,    rc2.ring_six,    sizeof(rc2.ring_six));
				memcpy(chrom[GB->num_chrom+i].ring_five,   rc2.ring_five,   sizeof(rc2.ring_five));
				ring_load_chrom_to_fa(FA, &chrom[GB->num_chrom+i]);
			}

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
			return metal_eval_init(n_atoms, n_types, ctx_max_pop,
			                      ad.lig_first, ad.lig_last,
			                      FA->permeability,
			                      ad.xyz.data(), ad.type.data(),
			                      ad.radius.data(), h_emat.data(), ns);
		});

		if (handle.ctx) {
			std::vector<double> h_genes = pack_genes_batch(n_genes);
			std::vector<double> h_com(pop_size), h_wal(pop_size), h_sas(pop_size);

			if (batch_n <= 1) {
				// Single-complex fast path — no batching overhead.
				metal_eval_batch(handle.ctx, pop_size, n_genes, h_genes.data(),
				                 h_com.data(), h_wal.data(), h_sas.data());
			} else {
				// Multi-complex path: pack per-complex atom data and queue for
				// batched dispatch.  When N concurrent workers all reach this
				// point in the same generation, they share one GPU kernel launch
				// (N × pop_size chromosomes per dispatch).
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

			unpack_gpu_results(h_com, h_wal, h_sas);
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
			// Keep the reference-calculation retry path enabled in GA workers.
			// The direct native probe uses recalc=1; forcing 0 here caused the
			// same pose to fall into the non-convergence penalty path.
			tl_vc[t].recalc    = 1;
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

			// Load this chromosome's ring pucker phases into the per-thread FA
			// so ic2cf reconstructs its puckered ring (no-op when inactive).
			ring_load_chrom_to_fa(&tl_fa[tid], &chrom[ii]);

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
			// v27 elitism: the top n_elite (lowest evalue → smallest pi after the
			// ascending QuickSort above) are exempt from the sharing reduction so
			// niching can never demote the running best out of the selection pool.
			if (pi < GB->n_elite)
				chrom[pi].fitnes = (double)(GB->num_chrom - pi);
			else
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

			// Compute ensemble thermodynamics (physical β = 1/kBT) and
			// SELECTION weights (β_sel = 1/T, matching the clustering
			// convention FA->beta). Using the physical β here would collapse
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
				// v27 elitism: top n_elite (lowest evalue → smallest pi after the
				// ascending QuickSort) are exempt from the sharing reduction.
				double blended = (1.0 - w) * rank_component + w * boltz_component;
				if (pi < GB->n_elite)
					chrom[pi].fitnes = blended * static_cast<double>(GB->num_chrom);
				else
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
				const bool reflig_seeded =
				    FA->reflig_nearest_count > 0 &&
				    i < popoffset + static_cast<int>(FA->reflig_seed_fraction *
				        static_cast<float>(GB->num_chrom - popoffset));
				if (reflig_seeded) {
					// Direct-mode native fallback: grid index 0 is the input pose
					// anchor. Keep it exactly so redocking starts from a physically
					// valid native-like chromosome instead of a nearby clash point.
					bool native_direct_seed =
					    FA->reflig_file[0] == '\0' && FA->reflig_hetatm_fallback &&
					    FA->resligand != NULL && gene_lim[0].min <= 0.0;
					int grid_idx = 0;
					if (!native_direct_seed) {
						// Explicit RefLig seeding: distribute K nearest grid points.
						int k = (i - popoffset) % FA->reflig_nearest_count;
						grid_idx = FA->reflig_nearest_grid[k];
					}
					chrom[i].genes[0].to_ic = static_cast<double>(grid_idx);
					chrom[i].genes[0].to_int32 = ictogene(&gene_lim[0],
					                                       static_cast<double>(grid_idx));
					for (int g = 1; g < GB->num_genes; g++) {
						if (!FA->map_par || !FA->opt_par) break;
						if (FA->map_par[g].typ == 1 || FA->map_par[g].typ == 2) {
							double ref_ic = FA->opt_par[g];
							chrom[i].genes[g].to_ic = ref_ic;
							chrom[i].genes[g].to_int32 = ictogene(&gene_lim[g], ref_ic);
						}
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
				} else if (FA->num_grd > 0) {
					// ── Cleft-biased GPA0 seeding ──────────────────────────────────
					// Seed gene[0] (rigid-body grid index) near cleftgrid[0], the
					// highest Voronoi contact density point.  Box-Muller Gaussian
					// with σ = max(3, num_grd/10) indices ≈ 2 Å spread across the
					// pocket.  Dramatically reduces the 479 k clashes/run caused by
					// blind uniform placement.  Chromosomes that still clash receive
					// the OOB penalty as normal and die through selection pressure.
					const int sigma_idx = std::max(3, FA->num_grd / 10);
					const double u1 = std::max(1e-10, RandomDouble(dice()));
					const double u2 = RandomDouble(dice());
					// Box-Muller N(0,1) → N(0, sigma_idx)
					const double z = std::sqrt(-2.0 * std::log(u1))
					                 * std::cos(2.0 * M_PI * u2);
					const int grid_idx = std::clamp(
					    static_cast<int>(std::round(z * sigma_idx)),
					    0, FA->num_grd - 1);
					chrom[i].genes[0].to_ic    = static_cast<double>(grid_idx);
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

			ring_load_chrom_to_fa(&p_fa[tid], &chrom[i]);

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
		// EMap typedef removed (unused) — plain gather loop used below
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
	// Thread-safe RNG (replaces non-reentrant rand())
	thread_local std::mt19937 tl_rng = flexaids_rng::make_thread_rng(0x9A800DULL);
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(tl_rng);
}
