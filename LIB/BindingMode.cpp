#include "BindingMode.h"
#include "fast_optics.hpp"
#include "SoftBetaFreeEnergy.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

// DatasetRunner reads sibling `.mcf` files (one app_evalue per line, head first)
// to build multi-member Shannon G̃. Format matches LIB/cluster.cpp exactly.
bool write_mcf_sidecar(const char* pdb_path,
                       const std::vector<double>& app_evalues_head_first)
{
	if (pdb_path == nullptr || pdb_path[0] == '\0') return false;
	std::string mcf_path(pdb_path);
	if (mcf_path.size() > 4 &&
	    mcf_path.substr(mcf_path.size() - 4) == ".pdb")
		mcf_path = mcf_path.substr(0, mcf_path.size() - 4) + ".mcf";
	else
		mcf_path += ".mcf";
	FILE* mf = std::fopen(mcf_path.c_str(), "w");
	if (!mf) return false;
	for (double v : app_evalues_head_first) {
		if (std::isfinite(v))
			std::fprintf(mf, "%.6f\n", v);
	}
	std::fclose(mf);
	return true;
}

/*****************************************\\
			BindingPopulation  
\\*****************************************/
// public (explicitely requires unsigned int temperature) *non-overloadable*
// BindingPopulation::BindingPopulation(unsigned int temp) : Temperature(temp)
BindingPopulation::BindingPopulation(FA_Global* pFA, GB_Global* pGB, VC_Global* pVC, chromosome* pchrom, genlim* pgene_lim, atom* patoms, resid* presidue, gridpoint* pcleftgrid, int num_chrom)
	: Temperature(pFA->temperature),
	  PartitionFunction(0.0),
	  nChroms(num_chrom),
	  FA(pFA),
	  GB(pGB),
	  VC(pVC),
	  chroms(pchrom),
	  gene_lim(pgene_lim),
	  atoms(patoms),
	  residue(presidue),
	  cleftgrid(pcleftgrid),
	  shannonS_population_(0.0),
	  shannon_cache_valid_(false)
{
}


void BindingPopulation::add_BindingMode(BindingMode& mode)
{
	for (std::vector<Pose>::iterator pose = mode.Poses.begin(); pose != mode.Poses.end(); ++pose)
	{
		this->PartitionFunction += pose->boltzmann_weight;
	}
	mode.set_energy();
	this->BindingModes.push_back(mode);
	this->shannon_cache_valid_ = false;  // Invalidate Shannon cache
	this->sorted_ = false;               // mark for lazy re-sort

	// Previously this called Entropize() on every insertion, re-sorting all
	// existing modes (each comparison invokes compute_energy() and rebuilds
	// the StatMechEngine). Inserting M modes was O(M² log M × N_poses).
	// Sorting is now lazy — done in output_Population() and accessor paths.
}


void BindingPopulation::Entropize()
{
	if (this->sorted_) return;

	for (std::vector<BindingMode>::iterator it = this->BindingModes.begin(); it != this->BindingModes.end(); ++it)
	{
		it->set_energy();
	}
	std::sort(this->BindingModes.begin(), this->BindingModes.end(), BindingPopulation::EnergyComparator());
	this->sorted_ = true;
}


int BindingPopulation::get_Population_size() { return this->BindingModes.size(); }



const BindingMode& BindingPopulation::get_binding_mode(int index) const
{
	if (index < 0 || index >= static_cast<int>(this->BindingModes.size()))
	{
		throw std::out_of_range(
			"BindingPopulation::get_binding_mode: index " +
			std::to_string(index) + " out of range [0, " +
			std::to_string(this->BindingModes.size()) + ")"
		);
	}
	return this->BindingModes[index];
}


BindingMode& BindingPopulation::get_binding_mode(int index)
{
	this->Entropize();
	if (index < 0 || index >= static_cast<int>(this->BindingModes.size()))
	{
		throw std::out_of_range(
			"BindingPopulation::get_binding_mode: index " +
			std::to_string(index) + " out of range [0, " +
			std::to_string(this->BindingModes.size()) + ")"
		);
	}
	return this->BindingModes[index];
}


// output BindingMode up to nResults results
void BindingPopulation::output_Population(int nResults, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints)
{
	// Output Population information ~= output clusters informations (*.cad)
	// Ensure the population is sorted by free energy before emitting results
	// (lazy sort: previously every add_BindingMode triggered a full re-sort).
	this->Entropize();

	// P1.5 additive diagnostic: explicit log when 0 binding modes (common friction for "best BindingMode")
	if (this->BindingModes.empty()) {
		std::cerr << "[P1 DIAGNOSTIC] BindingPopulation::output_Population: 0 binding modes after clustering (minPoints=" << minPoints << "). No best mode will be emitted. Check GA params, data quality, FastOPTICS/OPTICS settings, or RMSD clustering thresholds. This run may produce only placeholder 999 RMSD results.\n";
	}

	// Looping through BindingModes
	int num_result = 0;
	if (!nResults) nResults = this->get_Population_size() - 1; // if 0 is sent to this function, output all
	for (std::vector<BindingMode>::iterator mode = this->BindingModes.begin(); mode != this->BindingModes.end() && nResults > 0; ++mode, --nResults, ++num_result)
	{
		mode->output_BindingMode(num_result, end_strfile, tmp_end_strfile, dockinp, gainp, minPoints);
		// mode->output_dynamic_BindingMode(num_result,end_strfile, tmp_end_strfile, dockinp, gainp, minPoints);
	}
}


/// === NEW: Binding Population thermodynamic APIs ===
double BindingPopulation::compute_delta_G(const BindingMode& mode1, const BindingMode& mode2) const
{
	return mode1.get_free_energy() - mode2.get_free_energy();
}


statmech::StatMechEngine BindingPopulation::get_global_ensemble() const
{
	statmech::StatMechEngine global_engine(static_cast<double>(this->Temperature));

	// Build the population-level ensemble from raw microstates exactly once.
	// Do not feed per-mode Boltzmann probabilities back as multiplicities here:
	// that would reweight each pose by exp(-beta E) twice and distort the global
	// partition function, entropy, and any downstream thermodynamic diagnostic.
	for (const auto& mode : this->BindingModes)
	{
		for (const auto& pose : mode.Poses)
		{
			global_engine.add_sample(pose.total_energy(), 1.0);
		}
	}

	return global_engine;
}


statmech::StatMechEngine BindingPopulation::get_super_cluster_ensemble() const
{
	// Collect all pose energies across all binding modes
	std::vector<double> all_energies;
	for (const auto& mode : this->BindingModes)
		for (const auto& pose : mode.Poses)
			all_energies.push_back(pose.CF);

	if (all_energies.size() <= 4)
		return get_global_ensemble();  // too few poses for meaningful filtering

	// Build 1D energy points for lightweight FastOPTICS
	std::vector<fast_optics::Point> energy_pts(all_energies.size());
	for (size_t i = 0; i < all_energies.size(); ++i)
		energy_pts[i].coords = { all_energies[i] };

	fast_optics::FastOPTICS sc_optics(energy_pts,
		std::max(4, static_cast<int>(all_energies.size()) / 20));
	auto sc_indices = sc_optics.extractSuperCluster(
		fast_optics::ClusterMode::SUPER_CLUSTER_ONLY);

	if (sc_indices.empty())
		return get_global_ensemble();  // fallback if extraction yields nothing

	statmech::StatMechEngine sc_engine(static_cast<double>(this->Temperature));
	for (size_t idx : sc_indices)
		sc_engine.add_sample(all_energies[idx], 1.0);

	return sc_engine;
}


double BindingPopulation::get_shannon_entropy() const
{
	if (shannon_cache_valid_)
	{
		return shannonS_population_;
	}

	auto global_engine = get_global_ensemble();
	if (global_engine.size() == 0)
	{
		shannonS_population_ = 0.0;
		shannon_cache_valid_ = true;
		return 0.0;
	}

	const std::vector<double> weights = global_engine.boltzmann_weights();
	double shannon_nats = 0.0;
	for (double w : weights)
	{
		if (w > 0.0) shannon_nats -= w * std::log(w);
	}

	// Convert from dimensionless nats to thermodynamic units: S = kB * H
	shannonS_population_ = statmech::kB_kcal * shannon_nats;
	shannon_cache_valid_ = true;
	return shannonS_population_;
}


std::vector<std::vector<double>> BindingPopulation::get_deltaG_matrix() const
{
	int n = static_cast<int>(this->BindingModes.size());
	std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));

	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
		{
			double dg = compute_delta_G(this->BindingModes[i], this->BindingModes[j]);
			matrix[i][j] = dg;
			matrix[j][i] = -dg;
		}
	}

	return matrix;
}


/*****************************************\\
			  BindingMode
\\*****************************************/

// public constructor *non-overloadable*
BindingMode::BindingMode(BindingPopulation* pop)
	: Population(pop),
	  engine_(pop ? static_cast<double>(pop->Temperature) : 298.15),
	  thermo_cache_valid_(false),
	  vib_correction_cache_(0.0),
	  vib_cache_valid_(false),
	  energy(0.0)
{
}


// public method for pose addition
void BindingMode::add_Pose(Pose& pose)
{
	this->Poses.push_back(pose);
	this->thermo_cache_valid_ = false;
	this->vib_cache_valid_ = false;
}


/// === Cache rebuild infrastructure (Phase 1 + CCBM) ===
/// Uses pose.total_energy() (= CF + receptor_strain) for the true
/// multi-conformer free energy: F = -kT ln Σ_{r,i} exp(-β(E_CF(r,i) + E_strain(r)))
/// Physical-kB path used for diagnostic ledger + force_cf_rank_emission ranking.
void BindingMode::rebuild_engine() const
{
	if (thermo_cache_valid_)
	{
		return;
	}

	engine_.clear();
	for (const auto& pose : Poses)
	{
		engine_.add_sample(pose.total_energy(), 1.0);
	}
	thermo_cache_valid_ = true;
}


bool BindingMode::use_classic_entropy_ranking() const noexcept
{
	// Classic FlexAID: soft-β free energy over mode members (local Z) elects modes
	// when T>0 (SoftBetaFreeEnergy.h ≡ ACF). force_cf_rank_emission restores
	// physical per-mode StatMech ranking (rollback).
	if (!Population || Population->Temperature == 0)
		return false;
	if (Population->FA && Population->FA->force_cf_rank_emission)
		return false;
	return true;
}


// Collect mode-member CF values for SoftBeta (3Dsig / DatasetRunner identity).
// Election pb_clash re-weight (selection/search decoupling).
// FLEXAIDDS_PB_CLASH_ELECT_WEIGHT, when set, re-ranks modes as if the pb_clash
// penalty had been scored at the ELECTION weight instead of the SEARCH weight,
// WITHOUT changing the GA search fitness (which stays at pb_clash_weight so the
// near-native pose is still SAMPLED). Each member's stored CF already includes
// pb_clash at the search weight (chrom->cf.pb_clash = search_w * raw). To lift it
// to elect_w we add (elect_w/search_w - 1) * chrom->cf.pb_clash to the member
// energy. Validated offline on 1SJ0: search w=100 samples 2.47A but elects 6.43A;
// re-ranking at elect_w>=250 elects 2.67A (near-native). Default: no adjustment.
static double election_pb_scale(const FA_Global* FA)
{
	const char* e = std::getenv("FLEXAIDDS_PB_CLASH_ELECT_WEIGHT");
	if (!e || e[0] == '\0') return 0.0;               // 0.0 => disabled (no adjustment)
	const double elect_w = atof(e);
	const double search_w = (FA != nullptr) ? FA->pb_clash_weight : 0.0;
	if (elect_w <= 0.0 || search_w <= 0.0) return 0.0; // need both weights active
	return (elect_w / search_w) - 1.0;                 // additive factor on stored cf.pb_clash
}

static std::vector<double> soft_beta_mode_energies(const std::vector<Pose>& poses,
                                                   double pb_elect_factor)
{
	std::vector<double> energies;
	energies.reserve(poses.size());
	for (const auto& pose : poses) {
		double e = static_cast<double>(pose.CF);
		// Election re-weight: lift pb_clash from search weight to election weight.
		if (pb_elect_factor != 0.0 && pose.chrom != nullptr)
			e += pb_elect_factor * static_cast<double>(pose.chrom->cf.pb_clash);
		if (std::isfinite(e))
			energies.push_back(e);
	}
	return energies;
}


double BindingMode::compute_enthalpy() const
{
	if (!use_classic_entropy_ranking())
	{
		rebuild_engine();
		return engine_.compute().mean_energy;
	}
	// Soft-β H̃ over mode members only (≡ SoftBetaFreeEnergy / cluster ACF / DatasetRunner).
	// Local re-normalization — not global PartitionFunction weights.
	if (!Population || Poses.empty())
		return 0.0;
	const double T = static_cast<double>(Population->Temperature);
	return flexaids::soft_beta::free_energy(soft_beta_mode_energies(Poses, election_pb_scale(Population->FA)), T).H;
}


double BindingMode::compute_entropy() const
{
	if (!use_classic_entropy_ranking())
	{
		rebuild_engine();
		return engine_.compute().entropy;
	}
	// Soft-β S̃ = −Σ p_i ln p_i over mode members (nats). Same T as ACF / 3Dsig.
	if (!Population || Poses.empty())
		return 0.0;
	const double T = static_cast<double>(Population->Temperature);
	return flexaids::soft_beta::free_energy(soft_beta_mode_energies(Poses, election_pb_scale(Population->FA)), T).S;
}


double BindingMode::compute_energy() const
{
	if (!use_classic_entropy_ranking())
	{
		rebuild_engine();
		double nat_dg = (Population && Population->FA) ? Population->FA->natural_deltaG : 0.0;
		return engine_.compute().free_energy + compute_vibrational_correction() + nat_dg;
	}
	// Ranking objective: G̃ = H̃ − T·S̃ over mode members (SoftBetaFreeEnergy.h).
	// Identity: G̃ ≡ ACF ≡ E_min − T ln Z_local. Same formula as cluster.cpp emission
	// order and DatasetRunner S1 election. FlexAIDdS adds vib + optional NATURaL:
	//   F_rank = G̃_conf + (−T·S_vib) + ΔG_NATURaL
	// Vib is additive; it does not replace soft-β configurational ranking.
	if (!Population)
		return 0.0;
	const double T = static_cast<double>(Population->Temperature);
	const auto fe = flexaids::soft_beta::free_energy(soft_beta_mode_energies(Poses, election_pb_scale(Population->FA)), T);
	const double nat_dg =
		(Population->FA) ? Population->FA->natural_deltaG : 0.0;
	return fe.G + compute_vibrational_correction() + nat_dg;
}


/// === Thermodynamic APIs ===
statmech::Thermodynamics BindingMode::get_thermodynamics() const
{
	rebuild_engine();
	statmech::Thermodynamics td = engine_.compute();
	// Phase 3: include vibrational free energy correction in reported free energy
	td.free_energy += compute_vibrational_correction();
	// NATURaL: include co-translational ΔG (0.0 if assume_folded or not computed)
	td.free_energy += (Population && Population->FA) ? Population->FA->natural_deltaG : 0.0;
	return td;
}

statmech::ThermodynamicBreakdown BindingMode::get_thermodynamic_breakdown() const
{
	rebuild_engine();
	const double vib = compute_vibrational_correction();
	const double natural = (Population && Population->FA) ? Population->FA->natural_deltaG : 0.0;
	statmech::ThermodynamicBreakdown breakdown = engine_.compute_breakdown(
		vib,
		natural,
		0.0,
		vib != 0.0,
		natural != 0.0,
		false);
	if (!Poses.empty()) {
		breakdown.components = get_component_averages();
		breakdown.has_components = true;
	}
	return breakdown;
}

statmech::ComponentAverages BindingMode::get_component_averages() const
{
	rebuild_engine();
	std::vector<statmech::EnergyComponents> components;
	components.reserve(Poses.size());
	for (const auto& pose : Poses) {
		statmech::EnergyComponents c = pose.energy_components;
		c.total = pose.total_energy();
		c.cf = pose.CF;
		c.receptor_strain = pose.receptor_strain;
		components.push_back(c);
	}
	return engine_.component_averages(components);
}


double BindingMode::get_free_energy() const
{
	return get_thermodynamics().free_energy;
}


double BindingMode::get_heat_capacity() const
{
	return get_thermodynamics().heat_capacity;
}


std::vector<double> BindingMode::get_boltzmann_weights() const
{
	rebuild_engine();
	return engine_.boltzmann_weights();
}


double BindingMode::delta_G_relative_to(const BindingMode& reference) const
{
	return this->get_free_energy() - reference.get_free_energy();
}


std::vector<statmech::WHAMBin> BindingMode::free_energy_profile(
	const std::vector<double>& coordinates,
	int nbins
) const
{
	if (coordinates.size() != Poses.size())
	{
		throw std::invalid_argument(
			"free_energy_profile: coordinate vector size (" +
			std::to_string(coordinates.size()) +
			") must match number of poses (" +
			std::to_string(Poses.size()) + ")"
		);
	}

	if (nbins < 2 || nbins > 1000)
	{
		throw std::invalid_argument(
			"free_energy_profile: nbins must be in [2, 1000]"
		);
	}

	// Use total_energy() = CF + receptor_strain to match rebuild_engine().
	// Previously this used pose.CF only, so the FE profile and the main
	// thermodynamics were computed from different energies in CCBM mode.
	std::vector<double> energies;
	energies.reserve(Poses.size());
	for (const auto& pose : Poses)
	{
		energies.push_back(pose.total_energy());
	}

	return statmech::StatMechEngine::boltzmann_pmf(
		energies,
		coordinates,
		static_cast<double>(Population->Temperature),
		nbins
	);
}


int BindingMode::get_BindingMode_size() const { return this->Poses.size(); }



const Pose& BindingMode::get_pose(int index) const
{
	if (index < 0 || index >= static_cast<int>(this->Poses.size()))
	{
		throw std::out_of_range(
			"BindingMode::get_pose: index " +
			std::to_string(index) + " out of range [0, " +
			std::to_string(this->Poses.size()) + ")"
		);
	}
	return this->Poses[index];
}


void BindingMode::clear_Poses()
{
	this->Poses.clear();
	this->engine_.clear();
	this->thermo_cache_valid_ = false;
	this->vib_cache_valid_ = false;
}


void BindingMode::set_energy()
{
	this->energy = this->compute_energy();
}


void BindingMode::output_BindingMode(int num_result, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints)
{
	cfstr CF;
	resid* pRes = NULL;
	cfstr* pCF = NULL;

	char sufix[25];
	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];

	std::vector<Pose>::const_iterator Rep_lowCF = this->elect_Representative(false);
	std::vector<Pose>::const_iterator Rep_lowOPTICS = this->elect_Representative(true);

	// Elect the OPTICS density-center pose (lowest reachability distance) as the
	// Binding Mode representative, not the lowest-CF member. The center is the
	// most representative of the cluster's density; emitting the lowest-CF pose
	// instead made Fast OPTICS output indistinguishable from the CF algorithm.
	// Fall back to the lowest-CF pose only if no valid OPTICS center exists
	// (e.g. a singleton mode where every reachDist is undefined).
	std::vector<Pose>::const_iterator Rep = Rep_lowOPTICS;
	if (Rep == this->Poses.end() || isUndefinedDist(Rep->reachDist)) Rep = Rep_lowCF;

	// P2 (feature-flagged): PB-aware member promotion within this BindingMode.
	// When FLEXAIDDS_PB_AWARE_PROMOTION=1, elect the member with lowest stored
	// intermolecular wall energy (CF.wal) among poses near the density center.
	// Mode ranking / BindingPopulation sort keys are untouched — emission only.
	// Metric is crystal-blind (no RMSD-to-ref). Default OFF.
	bool pb_promotion = false;
	if (const char* env_pb = std::getenv("FLEXAIDDS_PB_AWARE_PROMOTION")) {
		pb_promotion = (env_pb[0] != '\0' && env_pb[0] != '0');
	}
	if (pb_promotion && this->Poses.size() > 1 && Rep != this->Poses.end()) {
		const float center_reach = Rep->reachDist;
		const double center_cf = Rep->CF;
		// Cap: keep candidates in the same basin (reachDist or CF window).
		const float reach_cap = 2.0f;  // same scale as typical cluster_rmsd
		const double cf_window = 5.0;  // kcal proxy; wide enough for near-center members
		std::vector<Pose>::const_iterator best = Rep;
		double best_wal = (Rep->chrom != nullptr) ? Rep->chrom->cf.wal : 1.0e300;
		double best_cf = Rep->CF;
		for (std::vector<Pose>::const_iterator it = this->Poses.begin();
		     it != this->Poses.end(); ++it) {
			if (it->chrom == nullptr) continue;
			// Basin filter: OPTICS reachDist when defined, else CF proximity.
			bool in_basin = true;
			if (!isUndefinedDist(center_reach) && !isUndefinedDist(it->reachDist)) {
				if (it->reachDist > center_reach + reach_cap) in_basin = false;
			} else if ((it->CF - center_cf) > cf_window) {
				in_basin = false;
			}
			if (!in_basin) continue;
			const double wal = it->chrom->cf.wal;
			const double cf = it->CF;
			// Lexicographic: min wall, then min CF (same mode, lower clash preferred).
			if (wal < best_wal - 1e-9 ||
			    (std::fabs(wal - best_wal) <= 1e-9 && cf < best_cf - 1e-9)) {
				best = it;
				best_wal = wal;
				best_cf = cf;
			}
		}
		Rep = best;
	}

	for (int k = 0; k < this->Population->GB->num_genes; ++k) this->Population->FA->opt_par[k] = Rep->chrom->genes[k].to_ic;

	// Ring pucker (LigandRingFlex Phase 2): emit the representative pose's
	// furanose conformation alongside its standard genes. No-op when inactive.
	if (this->Population->FA->ring_flex_active) {
		for (int s = 0; s < this->Population->FA->ring_n_sugars && s < MAX_RING_FLEX; ++s)
			this->Population->FA->ring_cur_phases[s] = Rep->chrom->ring_phases[s];
	}

	// Rebuild coordinates for output + populate per-optres CF breakdown.
	// Returned re-score is discarded for the REMARK CF line: at emission
	// Vcontacts can early-return a raw uncapped clash penalty that bypasses
	// the per-contact wall cap and blows up the reported CF. Report the stored
	// chromosome evalue — what the GA actually optimized.
	CF = ic2cf(this->Population->FA, this->Population->VC, this->Population->atoms, this->Population->residue, this->Population->cleftgrid, this->Population->GB->num_genes, this->Population->FA->opt_par);

	size_t remark_len = 0;
	remark[0] = '\0';
	safe_remark_cat(remark, "REMARK optimized structure\n", &remark_len);

	if (pb_promotion) {
		snprintf(tmpremark, MAX_REMARK,
			"REMARK emission_rep_policy=pb_aware_promotion (min CF.wal within mode; FLEXAIDDS_PB_AWARE_PROMOTION=1)\n");
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK,
			"REMARK emission_promoted=1 optics_center_cf=%8.5f elected_cf=%8.5f elected_wal=%8.5f\n",
			(Rep_lowOPTICS != this->Poses.end() ? Rep_lowOPTICS->CF : Rep->CF),
			Rep->CF,
			(Rep->chrom != nullptr ? Rep->chrom->cf.wal : 0.0));
		safe_remark_cat(remark, tmpremark, &remark_len);
	} else {
		snprintf(tmpremark, MAX_REMARK, "REMARK Fast OPTICS clustering algorithm used to output the OPTICS density center as Binding Mode representative\n");
		safe_remark_cat(remark, tmpremark, &remark_len);
	}

	snprintf(tmpremark, MAX_REMARK, "REMARK CF=%8.5f\n", Rep->chrom->evalue);
	safe_remark_cat(remark, tmpremark, &remark_len);
	snprintf(tmpremark, MAX_REMARK, "REMARK CF.app=%8.5f\n", Rep->chrom->app_evalue);
	safe_remark_cat(remark, tmpremark, &remark_len);

	for (int j = 0; j < this->Population->FA->num_optres; ++j)
	{
		pRes = &this->Population->residue[this->Population->FA->optres[j].rnum];
		pCF = &this->Population->FA->optres[j].cf;

		snprintf(tmpremark, MAX_REMARK, "REMARK optimizable residue %s %c %d\n", pRes->name, pRes->chn, pRes->number);
		safe_remark_cat(remark, tmpremark, &remark_len);

		snprintf(tmpremark, MAX_REMARK, "REMARK CF.com=%8.5f\n", pCF->com);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.sas=%8.5f\n", pCF->sas);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.wal=%8.5f\n", pCF->wal);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.pb_clash=%8.5f\n", pCF->pb_clash);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.con=%8.5f\n", pCF->con);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.gist=%8.5f\n", pCF->gist);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.hbond=%8.5f\n", pCF->hbond);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK Residue has an overall SAS of %.3f\n", pCF->totsas);
		safe_remark_cat(remark, tmpremark, &remark_len);
	}

	snprintf(tmpremark, MAX_REMARK, "REMARK Binding Mode:%d Best CF in Binding Mode:%8.5f OPTICS Center (CF):%8.5f Binding Mode Total CF:%8.5f Binding Mode Frequency:%d\n",
		num_result, Rep_lowCF->CF, Rep_lowOPTICS->CF, this->compute_energy(), this->get_BindingMode_size());
	safe_remark_cat(remark, tmpremark, &remark_len);
	{
		double vib_corr = this->compute_vibrational_correction();
		if (std::abs(vib_corr) > 1e-12) {
			snprintf(tmpremark, MAX_REMARK, "REMARK Vibrational correction (ENCoM heuristic): %10.4f kcal/mol\n", vib_corr);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}
	}
	// Canonical thermodynamic ledger for plugin / load_results (audited fields only)
	{
		const statmech::Thermodynamics td = this->get_thermodynamics();
		double T = 300.0;
		if (this->Population && this->Population->FA && this->Population->FA->temperature > 0) {
			T = static_cast<double>(this->Population->FA->temperature);
		} else if (this->Population && this->Population->Temperature > 0) {
			T = static_cast<double>(this->Population->Temperature);
		}
		snprintf(tmpremark, MAX_REMARK, "REMARK free_energy = %.6f\n", td.free_energy);
		safe_remark_cat(remark, tmpremark, &remark_len);
		// Ranking objective (soft-β G̃); free_energy above is physical StatMech ledger.
		snprintf(tmpremark, MAX_REMARK, "REMARK soft_beta_G = %.6f\n", this->compute_energy());
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK enthalpy = %.6f\n", td.mean_energy);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK entropy = %.8f\n", td.entropy);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK heat_capacity = %.8f\n", td.heat_capacity);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK temperature = %.2f\n", T);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK binding_mode = %d\n", num_result);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK pose_rank = 1\n");
		safe_remark_cat(remark, tmpremark, &remark_len);
	}
	for (int j = 0; j < this->Population->FA->npar; ++j)
	{
		snprintf(tmpremark, MAX_REMARK, "REMARK [%8.3f]\n", this->Population->FA->opt_par[j]);
		safe_remark_cat(remark, tmpremark, &remark_len);
	}

	if (this->Population->FA->refstructure == 1)
	{
		bool Hungarian = false;
		const double rmsd_raw = calc_rmsd(this->Population->FA, this->Population->atoms, this->Population->residue, this->Population->cleftgrid, this->Population->FA->npar, this->Population->FA->opt_par, Hungarian);
		snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n", rmsd_raw);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK rmsd_raw = %.5f\n", rmsd_raw);
		safe_remark_cat(remark, tmpremark, &remark_len);
		Hungarian = true;
		const double rmsd_sym = calc_rmsd(this->Population->FA, this->Population->atoms, this->Population->residue, this->Population->cleftgrid, this->Population->FA->npar, this->Population->FA->opt_par, Hungarian);
		snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n", rmsd_sym);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK rmsd_sym = %.5f\n", rmsd_sym);
		safe_remark_cat(remark, tmpremark, &remark_len);
	}
	snprintf(tmpremark, MAX_REMARK, "REMARK inputs: %s & %s\n", dockinp, gainp);
	safe_remark_cat(remark, tmpremark, &remark_len);
	snprintf(sufix, sizeof(sufix), "_%d_%d.pdb", minPoints, num_result);
	snprintf(tmp_end_strfile, MAX_PATH__, "%s%s", end_strfile, sufix);
	write_pdb(this->Population->FA, this->Population->atoms, this->Population->residue, tmp_end_strfile, remark);

	// ── Write member-CF sidecar (.mcf) for DatasetRunner Shannon G̃ ──
	// Same format as cluster.cpp: one app_evalue per line; head first, then
	// remaining finite members. Without this, FO path collapses to S̃=0, G̃=CF.
	{
		std::vector<double> mcf_vals;
		mcf_vals.reserve(this->Poses.size());
		if (Rep != this->Poses.end() && Rep->chrom != nullptr &&
		    std::isfinite(Rep->chrom->app_evalue)) {
			mcf_vals.push_back(Rep->chrom->app_evalue);
		}
		for (std::vector<Pose>::const_iterator it = this->Poses.begin();
		     it != this->Poses.end(); ++it) {
			if (it == Rep) continue;
			if (it->chrom == nullptr) continue;
			if (!std::isfinite(it->chrom->app_evalue)) continue;
			mcf_vals.push_back(it->chrom->app_evalue);
		}
		if (!mcf_vals.empty())
			(void)write_mcf_sidecar(tmp_end_strfile, mcf_vals);
	}
}


void BindingMode::output_dynamic_BindingMode(int num_result, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints)
{
	cfstr CF;
	resid* pRes = NULL;
	cfstr* pCF = NULL;

	char sufix[25];
	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];
	int nModel = 1;
	for (std::vector<Pose>::iterator Pose = this->Poses.begin(); Pose != this->Poses.end(); ++Pose, ++nModel)
	{
		for (int k = 0; k < this->Population->GB->num_genes; ++k) this->Population->FA->opt_par[k] = Pose->chrom->genes[k].to_ic;

		// Ring pucker (LigandRingFlex Phase 2): emit each pose's furanose
		// conformation alongside its standard genes. No-op when inactive.
		if (this->Population->FA->ring_flex_active) {
			for (int s = 0; s < this->Population->FA->ring_n_sugars && s < MAX_RING_FLEX; ++s)
				this->Population->FA->ring_cur_phases[s] = Pose->chrom->ring_phases[s];
		}

		// Re-score is discarded for the REMARK CF line (emission clash
		// early-return blows up the raw penalty); report the stored chromosome
		// evalue — what the GA optimized.
		CF = ic2cf(this->Population->FA, this->Population->VC, this->Population->atoms, this->Population->residue, this->Population->cleftgrid, this->Population->GB->num_genes, this->Population->FA->opt_par);

		size_t remark_len = 0;
		remark[0] = '\0';
		safe_remark_cat(remark, "REMARK optimized structure\n", &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK Fast OPTICS clustering algorithm used to output the lowest OPTICS ordering as Binding Mode representative\n");
		safe_remark_cat(remark, tmpremark, &remark_len);

		snprintf(tmpremark, MAX_REMARK, "REMARK CF=%8.5f\n", Pose->chrom->evalue);
		safe_remark_cat(remark, tmpremark, &remark_len);
		snprintf(tmpremark, MAX_REMARK, "REMARK CF.app=%8.5f\n", Pose->chrom->app_evalue);
		safe_remark_cat(remark, tmpremark, &remark_len);

		for (int j = 0; j < this->Population->FA->num_optres; ++j)
		{
			pRes = &this->Population->residue[this->Population->FA->optres[j].rnum];
			pCF = &this->Population->FA->optres[j].cf;

			snprintf(tmpremark, MAX_REMARK, "REMARK optimizable residue %s %c %d\n", pRes->name, pRes->chn, pRes->number);
			safe_remark_cat(remark, tmpremark, &remark_len);

			snprintf(tmpremark, MAX_REMARK, "REMARK CF.com=%8.5f\n", pCF->com);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.sas=%8.5f\n", pCF->sas);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.wal=%8.5f\n", pCF->wal);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.pb_clash=%8.5f\n", pCF->pb_clash);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.con=%8.5f\n", pCF->con);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.gist=%8.5f\n", pCF->gist);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK CF.hbond=%8.5f\n", pCF->hbond);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK Residue has an overall SAS of %.3f\n", pCF->totsas);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}

		// Mode-level ledger + per-pose identity (pose_rank = MODEL index)
		{
			const statmech::Thermodynamics td = this->get_thermodynamics();
			double T = 300.0;
			if (this->Population && this->Population->FA && this->Population->FA->temperature > 0) {
				T = static_cast<double>(this->Population->FA->temperature);
			} else if (this->Population && this->Population->Temperature > 0) {
				T = static_cast<double>(this->Population->Temperature);
			}
			snprintf(tmpremark, MAX_REMARK, "REMARK Binding Mode:%d Best CF in Binding Mode:%8.5f Binding Mode Frequency:%d\n",
				num_result, this->Poses.empty() ? 0.0 : this->Poses.front().CF, this->get_BindingMode_size());
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK free_energy = %.6f\n", td.free_energy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK soft_beta_G = %.6f\n", this->compute_energy());
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK enthalpy = %.6f\n", td.mean_energy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK entropy = %.8f\n", td.entropy);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK heat_capacity = %.8f\n", td.heat_capacity);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK temperature = %.2f\n", T);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK binding_mode = %d\n", num_result);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK pose_rank = %d\n", nModel);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}

		for (int j = 0; j < this->Population->FA->npar; ++j)
		{
			snprintf(tmpremark, MAX_REMARK, "REMARK [%8.3f]\n", this->Population->FA->opt_par[j]);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}

		if (this->Population->FA->refstructure == 1)
		{
			bool Hungarian = false;
			const double rmsd_raw = calc_rmsd(this->Population->FA, this->Population->atoms, this->Population->residue, this->Population->cleftgrid, this->Population->FA->npar, this->Population->FA->opt_par, Hungarian);
			snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n", rmsd_raw);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK rmsd_raw = %.5f\n", rmsd_raw);
			safe_remark_cat(remark, tmpremark, &remark_len);
			Hungarian = true;
			const double rmsd_sym = calc_rmsd(this->Population->FA, this->Population->atoms, this->Population->residue, this->Population->cleftgrid, this->Population->FA->npar, this->Population->FA->opt_par, Hungarian);
			snprintf(tmpremark, MAX_REMARK, "REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n", rmsd_sym);
			safe_remark_cat(remark, tmpremark, &remark_len);
			snprintf(tmpremark, MAX_REMARK, "REMARK rmsd_sym = %.5f\n", rmsd_sym);
			safe_remark_cat(remark, tmpremark, &remark_len);
		}
		snprintf(tmpremark, MAX_REMARK, "REMARK inputs: %s & %s\n", dockinp, gainp);
		safe_remark_cat(remark, tmpremark, &remark_len);

		snprintf(sufix, sizeof(sufix), "_%d_MODEL_%d.pdb", minPoints, num_result);
		snprintf(tmp_end_strfile, MAX_PATH__, "%s%s", end_strfile, sufix);
		if (Pose == this->Poses.begin() && Pose + 1 == this->Poses.end())
		{
			write_MODEL_pdb(true, true, nModel, this->Population->FA, this->Population->atoms, this->Population->residue, tmp_end_strfile, remark);
		}
		else if (Pose == this->Poses.begin())
		{
			write_MODEL_pdb(true, false, nModel, this->Population->FA, this->Population->atoms, this->Population->residue, tmp_end_strfile, remark);
		}
		else if (Pose + 1 == this->Poses.end())
		{
			write_MODEL_pdb(false, true, nModel, this->Population->FA, this->Population->atoms, this->Population->residue, tmp_end_strfile, remark);
		}
		else
		{
			write_MODEL_pdb(false, false, nModel, this->Population->FA, this->Population->atoms, this->Population->residue, tmp_end_strfile, remark);
		}
	}
}


std::vector<Pose>::const_iterator BindingMode::elect_Representative(bool useOPTICSorder) const
{
	std::vector<Pose>::const_iterator Rep = this->Poses.begin();
	for (std::vector<Pose>::const_iterator it = this->Poses.begin(); it != this->Poses.end(); ++it)
	{
		if (!useOPTICSorder && (Rep->CF - it->CF) > DBL_EPSILON) Rep = it;
		if (useOPTICSorder && it->reachDist < Rep->reachDist && !isUndefinedDist(it->reachDist)) Rep = it;
	}
	return Rep;
}


inline bool const BindingMode::operator<(const BindingMode& rhs) { return (this->get_cached_energy() < rhs.get_cached_energy()); }


/*****************************************\\
				  Pose
\\*****************************************/
Pose::Pose(chromosome* chrom, int index, int iorder, float dist, uint temperature, std::vector<float> vec)
	: chrom_index(index),
	  order(iorder),
	  reachDist(dist),
	  chrom(chrom),
	  CF(chrom->app_evalue),
	  energy_components(),
	  boltzmann_weight(0.0),
	  vPose(vec),
	  model_index(0),
	  model_coords(nullptr),
	  receptor_strain(0.0)
{
	energy_components.total = total_energy();
	energy_components.cf = CF;
	energy_components.hbond = chrom->cf.hbond;
	energy_components.gist = chrom->cf.gist_desolv;
	energy_components.metal = chrom->cf.metal_coord;
	energy_components.water = chrom->cf.gist;
	energy_components.complete = false;
	// Classic FlexAID soft-β: w = exp(−E/T). Physical 1/(kB·T) collapsed ranking
	// entropy to zero on CF magnitudes. PartitionFunction + classic BindingMode F
	// use these weights. Physical ledger uses StatMechEngine (unchanged).
	// Rollback of ranking product is FA->force_cf_rank_emission, not this weight.
	if (temperature > 0) {
		this->boltzmann_weight = std::exp(
			-static_cast<double>(chrom->app_evalue) / static_cast<double>(temperature));
	} else {
		this->boltzmann_weight = 0.0;
	}
}


Pose::~Pose() {}


inline bool const Pose::operator<(const Pose& rhs)
{
	if (this->order < rhs.order) return true;
	else if (this->order > rhs.order) return false;

	if (this->reachDist < rhs.reachDist) return true;
	else if (this->reachDist > rhs.reachDist) return false;

	if (this->chrom_index < rhs.chrom_index) return true;
	else if (this->chrom_index > rhs.chrom_index) return false;

	return false;
}


/*****************************************\\
     ENCoM vibrational correction (Phase 3)
\\*****************************************/

/*****************************************\\
     CCBM: Conformer-Coupled Binding Modes
\\*****************************************/

/// Helper: compute log-sum-exp of a vector of log-weights for numerical stability.
static double ccbm_log_sum_exp(const std::vector<double>& x) {
	if (x.empty()) return -std::numeric_limits<double>::infinity();
	double mx = *std::max_element(x.begin(), x.end());
	if (mx == -std::numeric_limits<double>::infinity()) return mx;
	double s = 0.0;
	for (double v : x)
		s += std::exp(v - mx);
	return mx + std::log(s);
}

/// Compute the maximum model_index across all poses in this binding mode.
static int ccbm_max_model_index(const std::vector<Pose>& poses) {
	int mx = 0;
	for (const auto& p : poses)
		if (p.model_index > mx) mx = p.model_index;
	return mx;
}

std::vector<double> BindingMode::conformer_populations() const
{
	if (Poses.empty()) return {};

	const double beta = 1.0 / (statmech::kB_kcal * static_cast<double>(Population->Temperature));
	const int n_models = ccbm_max_model_index(Poses) + 1;

	// Compute log-weights per model: ln w(r) = ln Σ_i exp(-β E_total(r,i))
	std::vector<std::vector<double>> model_log_weights(n_models);
	for (const auto& pose : Poses) {
		int r = pose.model_index;
		if (r >= 0 && r < n_models)
			model_log_weights[r].push_back(-beta * pose.total_energy());
	}

	// log Z_r = log-sum-exp of all pose weights for model r
	std::vector<double> log_Z_r(n_models, -std::numeric_limits<double>::infinity());
	for (int r = 0; r < n_models; ++r) {
		if (!model_log_weights[r].empty())
			log_Z_r[r] = ccbm_log_sum_exp(model_log_weights[r]);
	}

	// log Z_total = log-sum-exp over all models
	double log_Z_total = ccbm_log_sum_exp(log_Z_r);

	// p(r) = exp(log_Z_r - log_Z_total)
	std::vector<double> populations(n_models, 0.0);
	for (int r = 0; r < n_models; ++r) {
		if (log_Z_r[r] > -std::numeric_limits<double>::infinity())
			populations[r] = std::exp(log_Z_r[r] - log_Z_total);
	}

	return populations;
}


double BindingMode::receptor_conformational_entropy() const
{
	auto pops = conformer_populations();
	if (pops.empty()) return 0.0;

	// S_receptor = -kB Σ_r p(r) ln p(r)
	double S = 0.0;
	for (double p : pops) {
		if (p > 0.0)
			S -= p * std::log(p);
	}
	return statmech::kB_kcal * S;
}


double BindingMode::ligand_receptor_mutual_information() const
{
	if (Poses.empty()) return 0.0;

	const double beta = 1.0 / (statmech::kB_kcal * static_cast<double>(Population->Temperature));
	[[maybe_unused]] const int n_models = ccbm_max_model_index(Poses) + 1;	const int n_poses = static_cast<int>(Poses.size());

	// Compute joint Boltzmann weights: w_{r,i} = exp(-β E_total(r,i))
	std::vector<double> log_w(n_poses);
	for (int i = 0; i < n_poses; ++i)
		log_w[i] = -beta * Poses[i].total_energy();

	double log_Z = ccbm_log_sum_exp(log_w);

	// Joint probabilities: p(r,i) = exp(log_w_i - log_Z)
	std::vector<double> p_joint(n_poses);
	for (int i = 0; i < n_poses; ++i)
		p_joint[i] = std::exp(log_w[i] - log_Z);

	// Joint entropy: S_joint = -kB Σ_{r,i} p(r,i) ln p(r,i)
	double S_joint = 0.0;
	for (int i = 0; i < n_poses; ++i) {
		if (p_joint[i] > 0.0)
			S_joint -= p_joint[i] * std::log(p_joint[i]);
	}
	S_joint *= statmech::kB_kcal;

	// Marginal receptor entropy (already computed)
	double S_receptor = receptor_conformational_entropy();

	// Marginal ligand entropy: S_ligand = -kB Σ_i p(i) ln p(i)
	// where p(i) = Σ_r p(r,i_in_cluster) — for individual poses, each is unique
	// so p_ligand(i) = p_joint(i) summed over models that share the same pose index.
	// In practice, each (r,i) pair is a unique microstate, so ligand marginals
	// are: p_L(i) = sum over all models r of p(r, pose_i)
	// For CCBM, each pose IS a unique (r,i) tuple. The "ligand identity" is
	// the vPose coordinates. To compute marginal ligand entropy properly we
	// need to group poses by their ligand coordinates. However since each
	// GA chromosome produces a unique pose, the simplest correct approach is:
	//
	// For the mutual information decomposition, we use:
	//   I(L;R) = S_L + S_R - S_joint
	//
	// where S_L = -kB Σ_l p_L(l) ln p_L(l), with p_L(l) = Σ_r p(r,l)
	// Since poses are indexed per-conformer, the "ligand pose" marginal
	// groups all (r,*) that have the same ligand coordinates. But in the
	// GA each chromosome generates a unique ligand pose, so we identify
	// poses by their chrom_index (ligand identity).

	// Build marginal ligand distribution by chrom_index
	std::map<int, double> ligand_marginals;
	for (int i = 0; i < n_poses; ++i) {
		ligand_marginals[Poses[i].chrom_index] += p_joint[i];
	}

	double S_ligand = 0.0;
	for (const auto& [_, p_l] : ligand_marginals) {
		if (p_l > 0.0)
			S_ligand -= p_l * std::log(p_l);
	}
	S_ligand *= statmech::kB_kcal;

	// I(L;R) = S_L + S_R - S_joint
	double MI = S_ligand + S_receptor - S_joint;
	// Mutual information should be non-negative; clamp numerical errors
	if (MI < 0.0) MI = 0.0;

	return MI;
}


BindingMode::EntropyDecomposition BindingMode::decompose_entropy() const
{
	EntropyDecomposition ed;
	ed.S_vibrational = 0.0;

	if (Poses.empty()) {
		ed.S_total = 0.0;
		ed.S_ligand = 0.0;
		ed.S_receptor = 0.0;
		ed.I_mutual = 0.0;
		return ed;
	}

	const double beta = 1.0 / (statmech::kB_kcal * static_cast<double>(Population->Temperature));
	const int n_poses = static_cast<int>(Poses.size());

	// Joint Boltzmann distribution
	std::vector<double> log_w(n_poses);
	for (int i = 0; i < n_poses; ++i)
		log_w[i] = -beta * Poses[i].total_energy();

	double log_Z = ccbm_log_sum_exp(log_w);

	std::vector<double> p_joint(n_poses);
	for (int i = 0; i < n_poses; ++i)
		p_joint[i] = std::exp(log_w[i] - log_Z);

	// S_total (joint entropy) = -kB Σ p ln p
	double H_joint = 0.0;
	for (int i = 0; i < n_poses; ++i) {
		if (p_joint[i] > 0.0)
			H_joint -= p_joint[i] * std::log(p_joint[i]);
	}
	ed.S_total = statmech::kB_kcal * H_joint;

	// S_receptor
	ed.S_receptor = receptor_conformational_entropy();

	// S_ligand marginal
	std::map<int, double> ligand_marginals;
	for (int i = 0; i < n_poses; ++i)
		ligand_marginals[Poses[i].chrom_index] += p_joint[i];

	double H_ligand = 0.0;
	for (const auto& [_, p_l] : ligand_marginals) {
		if (p_l > 0.0)
			H_ligand -= p_l * std::log(p_l);
	}
	ed.S_ligand = statmech::kB_kcal * H_ligand;

	// I(L;R) = S_L + S_R - S_joint
	ed.I_mutual = ed.S_ligand + ed.S_receptor - ed.S_total;
	if (ed.I_mutual < 0.0) ed.I_mutual = 0.0;

	// Vibrational correction from ENCoM modes (if available)
	ed.S_vibrational = 0.0;
	if (Population && Population->FA && Population->FA->normal_modes) {
		double T = static_cast<double>(Population->Temperature);
		double vib_corr = compute_vibrational_correction();
		if (T > 0.0 && std::abs(vib_corr) > 1e-12)
			ed.S_vibrational = -vib_corr / T;  // correction = -T*S_vib, so S_vib = -correction/T
	}

	return ed;
}


double BindingMode::compute_vibrational_correction() const
{
	// ── Entropy-ablation hook: FLEXAIDDS_NO_TENCOM ───────────────────────────
	// Zeroes the tENCoM/ENCoM vibrational free-energy correction (-T·S_vib).
	// NOTE: this correction is derived from the *receptor* elastic-network
	// normal modes (atoms[0].eigen) and is cached pose-independently, so it is
	// an identical additive constant across every binding mode of a given
	// receptor — it cannot change intra-receptor pose ranking and is therefore a
	// structural NULL on RMSD success by construction. This arm exists to
	// confirm that null empirically (it can only shift the reported predicted_dG
	// column, not which pose is emitted).
	static const bool no_tencom = (std::getenv("FLEXAIDDS_NO_TENCOM") != nullptr);
	if (no_tencom) return 0.0;

	if (!this->Population->FA->normal_modes) return 0.0;

	// Return cached value if still valid
	if (this->vib_cache_valid_) return this->vib_correction_cache_;

	std::vector<encom::NormalMode> modes;
	int mode_count = this->Population->FA->normal_modes;
	const atom* atoms = this->Population->atoms;

	if (atoms && atoms[0].eigen)
	{
		for (int m = 0; m < mode_count; ++m)
		{
			if (!atoms[0].eigen[m]) continue;
			encom::NormalMode mode;
			mode.index = m + 1;
			mode.eigenvalue = static_cast<double>(atoms[0].eigen[m][0]);
			mode.frequency = std::sqrt(std::abs(mode.eigenvalue));
			modes.push_back(mode);
		}
	}

	if (modes.empty()) {
		this->vib_correction_cache_ = 0.0;
		this->vib_cache_valid_ = true;
		return 0.0;
	}

	double T = static_cast<double>(this->Population->Temperature);
	encom::VibrationalEntropy vs = encom::ENCoMEngine::compute_vibrational_entropy(modes, T);

	this->vib_correction_cache_ = -T * vs.S_vib_kcal_mol_K;
	this->vib_cache_valid_ = true;
	return this->vib_correction_cache_;
}
