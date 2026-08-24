#ifndef BINDINGMODE_H
#define BINDINGMODE_H

#include "gaboom.h"
#include "fileio.h"
#include "statmech.h"  // statistical mechanics engine
#include "encom.h"     // ENCoM vibrational entropy (Phase 3)
#include "ShannonThermoStack/ShannonThermoStack.h"  // Shannon configurational entropy

#include <cstddef>

//#define UNDEFINED_DIST FLT_MAX // Defined in FOPTICS as > than +INF
#define UNDEFINED_DIST -0.1f // Defined in FOPTICS as > than +INF
#define isUndefinedDist(a) ((a - UNDEFINED_DIST) <= FLT_EPSILON)

class BindingPopulation; // forward-declaration in order to access BindingPopulation* Population pointer

/// Write DatasetRunner-compatible member-CF sidecar next to a FO/cluster PDB.
/// Format matches cluster.cpp: one app_evalue per line (head first); finite only.
/// Replaces trailing `.pdb` with `.mcf`. Returns true if the file was written.
bool write_mcf_sidecar(const char* pdb_path,
                       const std::vector<double>& app_evalues_head_first);

/// Always emit a vibrational REMARK, including fail-closed 0.0.
/// Surface honesty only: this does not invent an eigenvalue channel or elect poses.
void format_vibrational_diagnostic_remark(char* buf, size_t buflen, double vib_corr);

/*****************************************\
			  Pose
\*****************************************/
struct Pose
{
	// friend class BindingPopulation;
	
	// public constructor :
	Pose(chromosome* chrom, int chrom_index, int order, float dist, uint temperature, std::vector<float>);
	~Pose();
	// public (default behavior when struct is used instead of class)
	int chrom_index;
	int order;
	float reachDist;
	chromosome* chrom;
	double CF;
	statmech::EnergyComponents energy_components;
	double boltzmann_weight;  // ← DEPRECATED: now computed by StatMechEngine
	std::vector<float> vPose;

	// ═══ CONFORMER-COUPLED BINDING MODE (CCBM) FIELDS ═══
	int    model_index;       // receptor conformer index (0 = single model)
	float* model_coords;      // pointer to this conformer's atom coordinates (NOT owned)
	double receptor_strain;   // E_strain(r) = E_conformer(r) - E_conformer(r_ref)
	//
	// The TOTAL energy of this microstate is:
	//   E_total = CF + receptor_strain
	// where CF is the ligand-receptor interaction + ligand internal energy
	// and receptor_strain is the conformational penalty for adopting conformer r.
	//
	// This means a high-affinity pose on a strained receptor conformer may
	// be disfavored compared to a medium-affinity pose on the relaxed conformer.
	double total_energy() const { return CF + receptor_strain; }

	inline bool const operator< (const Pose& rhs);
};

struct PoseClassifier
{
   inline bool operator() ( const Pose& Pose1, const Pose& Pose2 )
   {
       	if(Pose1.order < Pose2.order) return true;
       	else if(Pose1.order > Pose2.order) return false;
		if(Pose1.reachDist < Pose2.reachDist) return true;
		else if(Pose1.reachDist > Pose2.reachDist) return false;
		if(Pose1.chrom_index < Pose2.chrom_index) return true;
		else if(Pose1.chrom_index > Pose2.chrom_index) return false;
		
		return false;

		
   }
};
/*****************************************\
			  BindingMode
\*****************************************/
class BindingMode // aggregation of poses (Cluster)
{
	friend class BindingPopulation;
	
	public:
		explicit 	 					BindingMode(BindingPopulation*); // public constructor (explicitely needs a pointer to a BindingPopulation of type BindingPopulation*)

			void	 					add_Pose(Pose&);
			void	 					clear_Poses();
			int	  						get_BindingMode_size() const;
			const Pose&					get_pose(int index) const;
			
			// ═══ LEGACY INTERFACE (backward compatibility) ═══
			double	 					compute_energy() const;      // returns Helmholtz free energy F = H - TS
			double	 					compute_entropy() const;     // returns configurational entropy S
			double	 					compute_enthalpy() const;    // returns Boltzmann-weighted ⟨E⟩
			
			// ═══ NEW STATMECH API ═══
			statmech::Thermodynamics	get_thermodynamics() const;  // full thermo struct (F, S, H, Cv, σ_E)
			statmech::ThermodynamicBreakdown get_thermodynamic_breakdown() const; // explicit config/correction ledger
			statmech::ComponentAverages get_component_averages() const; // Boltzmann-weighted component diagnostics
			double	 					get_free_energy() const;     // alias for compute_energy()
			double	 					get_heat_capacity() const;   // heat capacity C_v
			std::vector<double>	 		get_boltzmann_weights() const; // weights for all poses
			double	 					delta_G_relative_to(const BindingMode& reference) const;  // ΔG between modes
			
			// ═══ ADVANCED: WHAM FREE ENERGY PROFILES ═══
			std::vector<statmech::WHAMBin> free_energy_profile(
				const std::vector<double>& coordinates,
				int nbins = 20
			) const;  // 1D FE profile along arbitrary coordinate
			
		std::vector<Pose>::const_iterator elect_Representative(bool useOPTICSordering) const;
			inline bool const 			operator<(const BindingMode&);

			// ═══ PUBLIC ACCESSORS (for bindings) ═══
			const std::vector<Pose>&	get_poses() const { return Poses; }

			// Cached free energy populated by set_energy(). Used by sort
			// comparators to avoid recomputing the StatMechEngine on each
			// comparison (was the source of O(N² log N) sort cost).
			double	get_cached_energy() const noexcept { return energy; }

			// ═══ CONFORMER-COUPLED ENSEMBLE THERMODYNAMICS (CCBM) ═══

			/// Receptor conformer population weights:
			/// p(r) = Σ_i w(r,i) / Σ_{r',i} w(r',i)
			/// where w(r,i) = exp(-β·E_total(r,i)) / Z
			/// This tells you which receptor conformer is most populated
			/// in the BOUND state — conformational selection vs induced fit!
			std::vector<double> conformer_populations() const;

			/// Shannon entropy of the receptor conformer distribution:
			/// S_receptor = -kB Σ_r p(r) ln p(r)
			/// Quantifies how much receptor flexibility contributes to binding entropy.
			/// If S_receptor ≈ 0: one dominant conformer (conformational selection)
			/// If S_receptor is high: many conformers contribute (induced fit / population shift)
			double receptor_conformational_entropy() const;

			/// Mutual information between ligand pose and receptor conformer:
			/// I(ligand; receptor) = S_ligand + S_receptor - S_joint
			/// This quantifies the coupling between ligand pose choice and
			/// receptor conformer selection. High MI = strongly coupled binding.
			double ligand_receptor_mutual_information() const;

			/// Decompose total entropy into:
			///   S_total = S_ligand_pose + S_receptor_conformer - I(ligand; receptor)
			/// This is the entropy decomposition that NO other docking tool can do.
			struct EntropyDecomposition {
				double S_total;          // total configurational entropy
				double S_ligand;         // marginal ligand pose entropy
				double S_receptor;       // marginal receptor conformer entropy
				double I_mutual;         // mutual information (coupling)
				double S_vibrational;    // fail-closed 0.0 on ranking path (no eigenvalue channel)
			};
			EntropyDecomposition decompose_entropy() const;

 	protected:
		std::vector<Pose> Poses;
		BindingPopulation* Population; // used to access the BindingPopulation
		
		// ═══ STATMECH ENGINE (replaces manual Boltzmann summation) ═══
		mutable statmech::StatMechEngine engine_;  // mutable: allows lazy evaluation in const methods
		mutable bool thermo_cache_valid_;          // track if engine_ matches current Poses
		mutable double vib_correction_cache_;      // cached vibrational correction value
		mutable bool   vib_cache_valid_;           // track if vib cache matches current state

		void	set_energy();                         // updates cached energy value
		void	rebuild_engine() const;               // populates engine_ from Poses (called on-demand)
		double	compute_vibrational_correction() const; // fail-closed 0.0; atom::eigen is eigenvectors
		// Classic FlexAID soft-β ranking when T>0 && !FA->force_cf_rank_emission.
		// Live objective is configurational (soft-β) only. Vib is called but
		// fail-closes to 0.0 and does not elect. Rollback: force_cf_rank_emission.
		bool	use_classic_entropy_ranking() const noexcept;

	private:
		void 	output_BindingMode(int num_result, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints);
		void	output_dynamic_BindingMode(int nBindingMode, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints);
		double energy;  // cached free energy value
};

/*****************************************\
			BindingPopulation  
\*****************************************/
class BindingPopulation
{
	friend class BindingMode;
	public:
		// Temperature is used for energy calculations of BindingModes
		unsigned int Temperature;
		
		// explicit 	BindingPopulation(unsigned int);// public constructor (explicitely needs an int representative of Temperature)
		explicit 	BindingPopulation(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom);
		 	// add new binding mode to population
		 	void	add_BindingMode(BindingMode&);
		 	// return the number of BindinMonde (size getter)
		 	int	 	get_Population_size();
		 	// access individual BindingMode by index (bounds-checked)
		 	const BindingMode&	get_binding_mode(int index) const;
		 	BindingMode&		get_binding_mode(int index);
		 	// output BindingMode up to nResults results
		 	void	output_Population(int nResults, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints);
		 	
		 	// ═══ NEW STATMECH API ═══
		 	/// Compute ΔG between two binding modes (relative binding free energy)
		 	double	compute_delta_G(const BindingMode& mode1, const BindingMode& mode2) const;
		 	/// Get global ensemble StatMechEngine aggregating all binding modes
		 	statmech::StatMechEngine get_global_ensemble() const;
		 	/// Convenience: log_Z from the full population ensemble (P1 grand canonical hook)
		 	/// Returns the canonical log partition function for use with TargetServer / GrandPartitionFunction.
		 	double get_log_Z() const;
		 	/// Get canonical engine (alias for get_global_ensemble for clarity in GPF paths)
		 	statmech::StatMechEngine get_canonical_engine() const { return get_global_ensemble(); }
		 	/// Get super-cluster filtered ensemble (dominant energy basin only)
		 	statmech::StatMechEngine get_super_cluster_ensemble() const;

		 	// ═══ POPULATION-LEVEL SHANNON ENTROPY ═══
		 	/// Shannon configurational entropy across all binding modes: S = -kB * sum(p_i * ln(p_i))
		 	double	get_shannon_entropy() const;
		 	/// ΔG matrix between all pairs of binding modes (upper triangle, row-major)
		 	std::vector<std::vector<double>> get_deltaG_matrix() const;
		 	// ═══ PUBLIC ACCESSORS (for bindings) ═══
		 	// Returns insertion order; call output_Population() or use indexed
		 	// get_binding_mode() for energy-sorted access.
		 	const std::vector<BindingMode>& get_binding_modes() const { return BindingModes; }

		 	// ═══ POPULATION-LEVEL SHANNON ENTROPY ═══
		 	/// Shannon configurational entropy across all binding modes: S = -kB * sum(p_i * ln(p_i))
		 	/// ΔG matrix between all pairs of binding modes (upper triangle, row-major)

	protected:
		double PartitionFunction;	// sum of all Boltzmann_weight (DEPRECATED: use StatMechEngine)
		int nChroms;	 			// n_chrom_snapshot input to clustergin function

		// FlexAID pointer
		FA_Global* 	FA;	 			// pointer to FA_Global struct
		GB_Global* 	GB;	 			// pointer to GB_Global struct
		VC_Global* 	VC;	 			// pointer to VC_Global struct
		chromosome* chroms;	 		// pointer to chromosomes' array
		genlim* gene_lim;	 		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		atom* atoms;	 			// pointer to atoms' array
		resid* residue;	 			// pointer to residues' array
		gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)
	
	private:

		std::vector< BindingMode > 	BindingModes;	// BindingMode container
		mutable double				shannonS_population_;	// cached Shannon entropy
		mutable bool				shannon_cache_valid_;	// cache validity flag
		bool						sorted_ = false;        // true iff BindingModes is currently sorted by energy

		void 	 					Entropize(); 	// Sort BindinModes by cached free energy (no-op if already sorted)

		struct EnergyComparator
		{
			inline bool operator() ( const BindingMode& BindingMode1, const BindingMode& BindingMode2 )
			{
				// Use the cached energy populated by set_energy() before
				// sorting. Calling compute_energy() here re-runs the
				// StatMechEngine on each comparison, which made
				// std::sort take O(N log N × N_poses) per Entropize().
				return (BindingMode1.get_cached_energy() < BindingMode2.get_cached_energy());
			}
		};
};
#endif
