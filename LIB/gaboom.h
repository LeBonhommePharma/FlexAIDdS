#ifndef GABOOM_H
#define GABOOM_H

#ifdef _WIN32
# ifndef _CRT_SECURE_NO_WARNINGS
#  define _CRT_SECURE_NO_WARNINGS
# endif
#endif

#include <random>
#include <functional>
#include <cstdint>
#include <cmath>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>

#include "flexaid.h"
#include "Vcontacts.h"

// Forward declaration for re-entrant GA context
struct GAContext;

#define MAX_NUM_GENES 100
#define MAX_NUM_CHROM 1000
#define MAX_GEN_LENGTH 32                         // in number of bits
#define MAX_RANDOM_VALUE 2147483647              // upper bound of 32-bit integer
#define SAVE_CHROM_FRACTION 0.05

#define QS_TYPE double
#define QS_ASC(a,b) ((a)-(b))
#define QS_DSC(a,b) ((b)-(a))
#define K(i,j,n) (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
//#define K(i,j,n) ( (i < j) ? (i*n+j) : (j*n+i) )


typedef std::mt19937 RNGType;

struct genelimits_struct{
	double max;
	double min;
	double del;
	double bin;
	double nbin;
	int map;              // mapping gene (maps into an array)
};
typedef struct genelimits_struct genlim;

struct gene_struct{
	int32_t  to_int32;
	double          to_ic;
};
typedef struct gene_struct gene;

struct chromosome_struct{
	gene* genes;    // pointer to genes array
	cfstr  cf;       // function evaluation of chromosome
	double evalue;   // NOT the apparent cf
	double app_evalue;   // THE apparent cf
	double fitnes;   // fitness score of chromosome
	double boltzmann_weight;  // Boltzmann probability p_i = exp(-β E_i) / Z
	double free_energy;       // per-pose free energy contribution (kcal/mol)
	char   status;   /* status, n -> eval is correct
			    o -> need to recalculate eval
			 */

	// ── Ring pucker side-channel (LigandRingFlex Phase 2) ──
	// Plain-old-data arrays so the chromosome stays trivially copyable: it is
	// allocated with malloc(), bulk-swapped by swap_chrom (struct copy), and
	// shallow-assigned in the STEADY replacement path — none of which would be
	// safe for a heap-owning std::vector member. Ring topology (atom indices,
	// sugar types, counts) is shared per-complex on FA->ring_flex_template; only
	// the variable conformer choice travels with each chromosome here.
	float   ring_phases[MAX_RING_FLEX]; // furanose sugar pucker phases (deg)
	uint8_t ring_six[MAX_RING_FLEX];    // 6-membered ring conformer indices
	uint8_t ring_five[MAX_RING_FLEX];   // 5-membered (non-sugar) conformer indices
};
typedef struct chromosome_struct chromosome;

struct GB_Global_struct{
	//long long    seed;

	int          num_chrom;
	int          num_genes;
	int          max_generations;

	double        alpha;
	double        peaks;
	double        scale;
	double        sig_share;

	double        mut_rate;
	double        cross_rate;
	double        ini_mut_prob;
	double        end_mut_prob;

	double        pbfrac;
	int           ssnum;
	int           intragenes;

	// adapt Genetic Algorithm operators (mutation,crossover)
	int          adaptive_ga;
	double       fit_avg;
	double       fit_max;
	double       fit_low;
	double       fit_high;
	double       k1,k2,k3,k4;

	int          outgen;
	int          rrg_skip;

	int          seed;
	int          num_print;
	int	     	print_int;

	// Entropy convergence (opt-in early termination)
	int          entropy_convergence;    // 0=off (default), 1=on
	int          entropy_check_interval; // check every N generations (default: 10)
	int          entropy_window;         // plateau window size (default: 5)
	double       entropy_rel_threshold;  // relative change threshold (default: 0.01)

	char         pop_init_method[9];
	char         pop_init_file[MAX_PATH__];
	char         fitness_model[9];
	char         rep_model[9];
	int          duplicates;

	// ── Entropy-driven fitness parameters ──
	double       entropy_weight;      // blending: 0.0 = pure rank, 1.0 = full Boltzmann (default 0.5)
	int          entropy_interval;    // compute ensemble thermo every N generations (0 = off)
	int          use_shannon;         // Level-3 H(ω) ligand vibrational-entropy monitor (diagnostic only; FLEXAIDDS_USE_SHANNON). 0=off (default)

	// ── Diversity monitoring (entropy collapse mitigation) ──
	int          diversity_monitoring;            // 0=off (default), 1=on
	int          diversity_check_interval;        // check every N generations (default: 10)
	double       diversity_collapse_threshold;    // normalized allele entropy threshold (default: 0.3)
	double       catastrophic_mutation_fraction;  // fraction of population to re-randomize (default: 0.2)
	int          catastrophic_mutation_count;     // counter: how many times triggered

	// ── P5: periodic BOOM random injection (diversity insurance) ──
	// Every boom_inject_interval generations, replace the worst
	// (boom_inject_fraction × num_chrom/2) chromosomes with fresh random
	// individuals (not seeds), preserving the better half. Distinct from the
	// BOOM *reproduction* model and from collapse-triggered catastrophic mutation:
	// this fires unconditionally on a fixed cadence to keep multiple binding-mode
	// hypotheses alive against the deepest (often false) VCT minimum. 0 = off.
	int          boom_inject_interval;            // generations between injections (default: 100)
	double       boom_inject_fraction;            // fraction of the worst half to replace (default: 1.0)
	int          boom_inject_count;               // counter: how many times triggered

	// ── True GA elitism (v27) ──
	// Number of lowest-CF individuals snapshotted before boom injection and
	// niche-sharing, then restored over the worst of the reproduced population
	// so the global best is never lost to diversity pressure.  Elites are also
	// exempt from the PSHARE/SMFREE sharing fitness reduction.  0 = disabled,
	// default 1 (single-elite). Configurable via FLEXAIDDS_N_ELITE.
	int          n_elite;

	// ── Metal GPU multi-complex batching ──
	// When > 1 concurrent GA workers share a process, they queue up N entries
	// before issuing a single N×pop_size GPU kernel dispatch.
	// 1 = single-complex fast path (default, zero overhead).
	int          metal_batch_n;

};
typedef struct GB_Global_struct GB_Global;

// Density Peaks clustering algorithm data structure definitions
struct ClusterChrom
{
	uint index;						// original index in chrommose* chrom (at the input of density_cluster function)
	bool isHalo;					// is part of cluster core?
	bool isCenter;					// is cluster center?
	bool isBorder;					// is part of the border region
	bool isClustered;
	chromosome* Chromosome;			// Chromosomes list
	int Cluster;					// Assigned Cluster
	int Density;					// Density of points in distance cut-off
	double CF;						// Complementarity Function value
	float PiDi;						// Density x Distance
	float Distance;					// Nearest highest density peak distance
	float Coord[3*MAX_ATM_HET];		// Cartesian Coordinates
	struct ClusterChrom* DP;		// Nearest Density Peak (point of higher density)
}; typedef struct ClusterChrom ClusterChrom;

struct Cluster_struct
{
       int ID;						 // assigned cluster number (ID)
       int Frequency;				 // observation frequency of this cluster (number of representatives in cluster)
       double totCF;
       double lowestCF;
       // Pointer to best CF value in cluster
       ClusterChrom* Representative; // Pointer to the ClusterChrom individual with the lowest CF in cluster
       // Pointer to Cluster Center
       ClusterChrom* Center;		 // Queue of ClusterChrom (first element is the cluster center)
};
typedef struct Cluster_struct DPcluster;
/***********************************************************************/
/* 1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/* 1         2         3         4         5         6         7*/
/***********************************************************************/
int   GA(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome** chrom,chromosome** chrom_snapshot,genlim** gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,char gainpfile[], int* memchrom, cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int, double*), GAContext* ctx = nullptr);
int   check_state(char* pausefile, char* abortfile, char* stopfile, int interval);
void  QuickSort(chromosome*, int, int, bool);
void  QuickSort_Clusters(int*, int*, double*, double*, int*, int, int);
void  swap_clusters(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*);
void  crossover(gene *john,gene *mary,int num_genes, int intragenes);
void  mutate(gene *john,int num_genes,double mut_rate);
void  bin_print(int dec,int len);
void  read_gainputs(FA_Global* FA,GB_Global* GB,int*,int*,char file[]);
int   deelig_search(struct deelig_node_struct* root_node, int* deelig_list, int fdih);
int   filter_deelig(FA_Global* FA, GB_Global* GB, chromosome* chrom, gene* genes, int ci, atom* atoms, const genlim* gene_lim,
		   std::function<int32_t()> & dice);

void   set_gene_lim(FA_Global* FA, GB_Global* GB, genlim* gene_lim);
long int read_pop_init_file(FA_Global* FA, GB_Global* GB, genlim* gene_lim, char* pop_init_file);
void  set_bins(genlim* gene_lim);
void  set_bins(genlim* gene_lim, int num_genes);
double calc_poss(genlim* gene_lim, int num_genes);
void validate_dups(GB_Global* GB, genlim* gene_lim, int num_genes);
double genetoic(const genlim* gene_lim, int32_t gene);
int ictogene(const genlim* gene_lim, double ic);

int 	RandomInt(double frac);
double 	RandomDouble();
double 	RandomDouble(int32_t dice);

void  	swap_chrom(chromosome * x, chromosome * y);
void  	copy_chrom(chromosome* dest, const chromosome* src, int num_genes);
int   	remove_dups(chromosome* list, int num_chrom, int num_genes);

FILE* get_update_file_ptr(FA_Global* FA);
void 	close_update_file_ptr(FA_Global* FA, FILE* outfile_ptr);

void  	generate_random_individual(FA_Global* FA, GB_Global* GB, atom* atoms, gene* genes, const genlim* gene_lim,
				 std::function<int32_t()> &, int from_gene, int to_gene);
void populate_chromosomes(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom,
                        const genlim* gene_lim, atom* atoms,resid* residue,
                        gridpoint* cleftgrid,char method[],
                        cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*),
                        char file[], long int at, int offset, int print,
                        std::function<int32_t()> &,
                        std::unordered_map<size_t, int> &);
cfstr 	eval_chromosome(FA_Global* FA,GB_Global* GB,VC_Global* VC,const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,gene* john, cfstr (*function)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));
void  	calculate_fitness(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,char method[],int pop_size, int print, cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int, double*));
int reproduce(FA_Global* FA,GB_Global* GB,VC_Global* VC, chromosome* chrom,
             const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,
             char rmodel[], double mutprob, double crossprob, int print,
             std::function<int32_t()> &,
             std::unordered_map<size_t, int> &,
             cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));
void  	print_pop(const chromosome* chrom,const genlim* gene_lim,int numc, int numg);
void  	print_chrom(const chromosome* chrom, int num_genes, int real_flag);
void  	print_chrom(const gene* genes, int num_genes, int real_flag);
void  	print_par(const chromosome* chrom,const genlim* gene_lim,int num_chrom,int num_genes, FILE* outfile_ptr);
int   	roullete_wheel(const chromosome* chrom,int n);
int   	cmp_chrom2pop(chromosome* chrom,chromosome* c, int num_genes,int start, int last);
int   	cmp_chrom2pop_int(const chromosome* chrom,const gene* genes, int num_genes,int start, int last);
int   	cmp_chrom2rotlist(psFlexDEE_Node psFlexDEE_INI_Node, const chromosome* chrom, const genlim* gene_lim,int gene_offset, int num_genes, int tot, int num_nodes);
int   	cmp_chrom2pop(const chromosome* chrom,const gene* genes, int num_genes,int start, int last);
void  	save_snapshot(chromosome* chrom_snapshot, const chromosome* chrom, int num_chrom, int num_genes);
void  	cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC,chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue,gridpoint* cleftgrid, int memchrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
void  	DensityPeak_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int memchrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
void 	FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
//long long time_seed();
/* A4b: early_exit_sq = sig_share² to prune distant pairs; 0.0 = full compute */
double 	calc_rmsp(int npar, const gene* g1, const gene* g2, const optmap* map_par, gridpoint* cleftgrid, double early_exit_sq = 0.0);
void 	write_par(const chromosome* chrom,const genlim* gene_lim,int ger, char* outfile,int num_chrom,int num_genes);
void 	adapt_prob(GB_Global* GB,double fitnes1,double fitnes2, double* mut_prob, double* cross_prob);
void 	fitness_stats(GB_Global* GB, const chromosome* chrom,int nchrom);
float  	calc_rmsd_chrom(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, int chrom_a, int chrom_b, float*, float*, bool calc_rmsd, int* out_n_atoms = nullptr); // calculates RMSD between chromossomes
int    	write_rrg(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue, gridpoint* cleftgrid, char* outfile);        // writes GA output during simulation
int    	write_rrd(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int* Clus_GAPOP,float* Clus_RMSDT,char outfile[]);
int 	write_DensityPeak_rrd(FA_Global* FA, GB_Global* GB, const chromosome* chrom, const genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, ClusterChrom* Chrom, DPcluster* Clust, float* RMSD, char outfile[]);
void   	partition_grid(FA_Global* FA,chromosome* chrom,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,int pop_size, int expansion);        // partition grid size where favorable conformations are found
void   	slice_grid(FA_Global* FA,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid);                      // slice grid symmetrically in half

// Density Peaks Clustering algorithm function declarations
void 	QuickSort_Cluster_by_CF(DPcluster* Clust, bool Entropic, int beg, int end);
void 	swap_clusters(DPcluster* xClust, DPcluster* yClust);
float 	getDistanceCutoff(float* RMSD, int num_chrom);
void 	QuickSort_ChromCluster_by_CF(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void 	QuickSort_ChromCluster_by_higher_Density(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void 	QuickSort_ChromCluster_by_lower_Density(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void 	swap_elements(ClusterChrom* Chrom, ClusterChrom* ChromX, ClusterChrom* ChromY, int num_chrom);
float 	calculate_stddev(ClusterChrom* Chrom, int num_chrom);
float 	calculate_mean(ClusterChrom* Chrom, int num_chrom);
int 	DistanceComparator(const void*, const void*);

size_t	hash_genes(const gene* g, int n);

#endif // include guard