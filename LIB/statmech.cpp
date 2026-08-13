// statmech.cpp — Statistical Mechanics Engine implementation
//
// Notation:
//   β  = 1/(kB T)
//   Z  = Σ_i  n_i exp(−β E_i)        (canonical partition function)
//   F  = −kT ln Z                      (Helmholtz free energy)
//   ⟨E⟩ = (1/Z) Σ_i  n_i E_i exp(−β E_i)
//   C_v = (⟨E²⟩ − ⟨E⟩²) / (kT²)      (heat capacity)
//   S  = (⟨E⟩ − F) / T                 (entropy)
//
// All sums use log-sum-exp for numerical stability when energies span
// hundreds of kcal/mol (common in docking).
//
// Hardware dispatch (runtime via UnifiedHardwareDispatch layer):
//   1. AVX-512 16-wide SIMD (+ OpenMP)
//   2. Eigen3 vectorised array ops (auto-vectorises to AVX2/AVX-512)
//   3. OpenMP parallel reductions for large ensembles
//   4. Scalar fallback (always available)

#include "statmech.h"
#include "UnifiedHardwareDispatch.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <map>
#include <vector>

#include <Eigen/Dense>

#ifdef _OPENMP
#  include <omp.h>
#endif

#if defined(__AVX512F__) && defined(__AVX512DQ__)
#  include <immintrin.h>
#  define STATMECH_HAS_AVX512 1
#else
#  define STATMECH_HAS_AVX512 0
#endif

namespace statmech {

// Threshold above which OpenMP parallelisation pays off for reductions.
[[maybe_unused]] static constexpr std::size_t OMP_THRESHOLD = 4096;

// ─── scientific provenance / claim validity ─────────────────────────────────

namespace {

bool has_artifact_sha256(const std::string& value) noexcept {
    constexpr char prefix[] = "sha256:";
    constexpr char historical_filler[] =
        "3f7a9c2b1e4d5f6a7b8c9d0e1f2a3b4c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a";
    constexpr char empty_artifact_sha256[] =
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855";
    constexpr std::size_t prefix_length = sizeof(prefix) - 1;
    if (value.size() != prefix_length + 64)
        return false;
    for (std::size_t index = 0; index < prefix_length; ++index) {
        if (value[index] != prefix[index])
            return false;
    }

    bool seen[16] = {};
    std::size_t distinct = 0;
    bool is_historical_filler = true;
    bool is_empty_artifact = true;
    for (std::size_t index = 0; index < 64; ++index) {
        const char character = value[prefix_length + index];
        unsigned int nibble = 0;
        if (character >= '0' && character <= '9')
            nibble = static_cast<unsigned int>(character - '0');
        else if (character >= 'a' && character <= 'f')
            nibble = static_cast<unsigned int>(character - 'a' + 10);
        else if (character >= 'A' && character <= 'F')
            nibble = static_cast<unsigned int>(character - 'A' + 10);
        else
            return false;

        if (!seen[nibble]) {
            seen[nibble] = true;
            ++distinct;
        }
        const char filler_character = historical_filler[index];
        const unsigned int filler_nibble =
            filler_character <= '9'
                ? static_cast<unsigned int>(filler_character - '0')
                : static_cast<unsigned int>(filler_character - 'a' + 10);
        is_historical_filler = is_historical_filler && nibble == filler_nibble;
        const char empty_character = empty_artifact_sha256[index];
        const unsigned int empty_nibble =
            empty_character <= '9'
                ? static_cast<unsigned int>(empty_character - '0')
                : static_cast<unsigned int>(empty_character - 'a' + 10);
        is_empty_artifact = is_empty_artifact && nibble == empty_nibble;
    }
    return distinct >= 3 && !is_historical_filler && !is_empty_artifact;
}

ScientificProvenance provenance_for_breakdown(
    const ScientificProvenance& source,
    double G_vib_kcal_mol,
    double G_natural_kcal_mol,
    double G_other_kcal_mol,
    bool has_vib,
    bool has_natural,
    bool has_other)
{
    // Correction terms do not yet carry independent artifact receipts. Any
    // supplied correction therefore makes the aggregate ledger proxy-only,
    // even when the configurational ensemble itself has physical provenance.
    // The numerical fields remain unchanged.
    const bool has_unreceipted_correction =
        has_vib || has_natural || has_other ||
        G_vib_kcal_mol != 0.0 ||
        G_natural_kcal_mol != 0.0 ||
        G_other_kcal_mol != 0.0;
    return has_unreceipted_correction ? ScientificProvenance{} : source;
}

} // namespace

bool ScientificProvenance::allows_canonical_physical_claim() const noexcept {
    if (schema_version != kScientificProvenanceSchemaVersion)
        return false;

    const bool calibrated_energy =
        energy_domain == EnergyDomain::CalibratedKcalPerMol &&
        has_artifact_sha256(energy_provenance);
    const bool physical_measure =
        (ensemble_measure == EnsembleMeasure::EnumeratedMicrostates ||
         ensemble_measure == EnsembleMeasure::WeightedQuadrature) &&
        has_artifact_sha256(measure_provenance);

    return calibrated_energy && physical_measure;
}

bool ScientificProvenance::allows_binding_physical_claim() const noexcept {
    return allows_canonical_physical_claim() &&
           reference_state == ReferenceState::MatchedAssociationCycle &&
           has_artifact_sha256(reference_provenance);
}

bool ScientificProvenance::is_proxy_only() const noexcept {
    return claim_validity() == ClaimValidity::ProxyOnly;
}

ClaimValidity ScientificProvenance::claim_validity() const noexcept {
    if (allows_binding_physical_claim())
        return ClaimValidity::BindingPhysical;
    if (allows_canonical_physical_claim())
        return ClaimValidity::CanonicalPhysical;
    return ClaimValidity::ProxyOnly;
}

ScientificProvenance make_contact_function_optimizer_provenance(
    ReferenceState reference_state)
{
    ScientificProvenance provenance;
    provenance.energy_domain = EnergyDomain::ContactFunctionArbitraryUnits;
    provenance.ensemble_measure = EnsembleMeasure::OptimizerSamples;
    provenance.reference_state = reference_state;
    provenance.energy_provenance =
        "FlexAID Voronoi/contact-function score; no kcal/mol calibration";
    provenance.measure_provenance =
        "optimizer-selected, deduplicated and/or clustered GA pose records";
    provenance.reference_provenance =
        reference_state == ReferenceState::BoundOnly
            ? "bound pose ensemble only"
            : "no matched association-cycle artifact";
    return provenance;
}

// ─── construction ────────────────────────────────────────────────────────────

StatMechEngine::StatMechEngine(double temperature_K)
    : StatMechEngine(temperature_K, ScientificProvenance{})
{
}

StatMechEngine::StatMechEngine(double temperature_K,
                               ScientificProvenance provenance)
    : T_(temperature_K)
    , beta_(1.0 / (kB_kcal * temperature_K))
    , beta_selection_(1.0 / temperature_K)   // 1/T convention for GA/cluster selection (P1)
    , provenance_(provenance)
{
    if (temperature_K <= 0.0)
        throw std::invalid_argument("StatMechEngine: temperature must be > 0");
}

// ─── add_sample ──────────────────────────────────────────────────────────────

void StatMechEngine::add_sample(double energy, double multiplicity) {
    ensemble_.push_back({energy, multiplicity});
}

// ─── log_sum_exp ─────────────────────────────────────────────────────────────

double StatMechEngine::log_sum_exp(std::span<const double> x) {
    // Delegate to the unified hardware dispatch layer which handles
    // AVX-512, Eigen, OpenMP, and scalar paths with runtime selection.
    return flexaids::log_sum_exp_dispatch(x);
}

// ─── compute ─────────────────────────────────────────────────────────────────

Thermodynamics StatMechEngine::compute() const {
    if (ensemble_.empty())
        throw std::runtime_error("StatMechEngine::compute: empty ensemble");

    const std::size_t N = ensemble_.size();
    bool any_positive_multiplicity = false;
    for (std::size_t i = 0; i < N; ++i) {
        if (!std::isfinite(ensemble_[i].energy))
            throw std::runtime_error("StatMechEngine::compute: non-finite sample energy");
        if (std::isfinite(ensemble_[i].count) && ensemble_[i].count > 0.0) {
            any_positive_multiplicity = true;
        }
    }
    if (!any_positive_multiplicity)
        throw std::runtime_error("StatMechEngine::compute: no positive-multiplicity samples");

    // Build array of log-weights:  w_i = ln(n_i) − β E_i
    std::vector<double> log_w(N);

    // Eigen-vectorised log-weight construction
    {
        Eigen::ArrayXd counts(static_cast<Eigen::Index>(N));
        Eigen::ArrayXd energies(static_cast<Eigen::Index>(N));
        for (std::size_t i = 0; i < N; ++i) {
            double n = ensemble_[i].count;
            if (!std::isfinite(n) || n < 0.0) n = 0.0;
            counts(static_cast<Eigen::Index>(i))   = n;
            energies(static_cast<Eigen::Index>(i)) = ensemble_[i].energy;
        }
        Eigen::ArrayXd lw = counts.log() - beta_ * energies;
        Eigen::Map<Eigen::ArrayXd>(log_w.data(), static_cast<Eigen::Index>(N)) = lw;
    }

    double lnZ = log_sum_exp(log_w);

    // ⟨E⟩  = (1/Z) Σ n_i E_i exp(−β E_i)
    //       = exp(−lnZ) Σ E_i exp(log_w_i)
    // To keep stability: ⟨E⟩ = Σ E_i exp(log_w_i − lnZ)
    double E_avg  = 0.0;
    double E2_avg = 0.0;

    // Build contiguous energy array for vectorised paths.
    std::vector<double> energies_vec(N);
    for (std::size_t i = 0; i < N; ++i)
        energies_vec[i] = ensemble_[i].energy;

#if STATMECH_HAS_AVX512
    // AVX-512 path: 8-wide fused probability × energy moment accumulation.
    if (N >= 16) {
        __m512d v_Eavg  = _mm512_setzero_pd();
        __m512d v_E2avg = _mm512_setzero_pd();
        __m512d v_lnZ   = _mm512_set1_pd(lnZ);

        std::size_t i = 0;
        for (; i + 7 < N; i += 8) {
            __m512d v_lw = _mm512_loadu_pd(log_w.data() + i);
            __m512d v_E  = _mm512_loadu_pd(energies_vec.data() + i);

            // p_i = exp(log_w_i - lnZ)
            __m512d v_arg = _mm512_sub_pd(v_lw, v_lnZ);
            alignas(64) double tmp_exp[8];
            _mm512_storeu_pd(tmp_exp, v_arg);
            for (int k = 0; k < 8; ++k) tmp_exp[k] = std::exp(tmp_exp[k]);
            __m512d v_p = _mm512_loadu_pd(tmp_exp);

            // Accumulate p_i * E_i and p_i * E_i^2
            v_Eavg  = _mm512_fmadd_pd(v_p, v_E, v_Eavg);
            v_E2avg = _mm512_fmadd_pd(v_p, _mm512_mul_pd(v_E, v_E), v_E2avg);
        }
        E_avg  = _mm512_reduce_add_pd(v_Eavg);
        E2_avg = _mm512_reduce_add_pd(v_E2avg);

        // Scalar tail
        for (; i < N; ++i) {
            double p_i = std::exp(log_w[i] - lnZ);
            double Ei  = energies_vec[i];
            E_avg  += p_i * Ei;
            E2_avg += p_i * Ei * Ei;
        }
    } else
#endif

    // Eigen vectorised path: auto-vectorises to AVX2/AVX-512 via Eigen's backend.
    if (N >= 16) {
        Eigen::Map<const Eigen::ArrayXd> lw(log_w.data(), static_cast<Eigen::Index>(N));
        Eigen::Map<const Eigen::ArrayXd> E(energies_vec.data(), static_cast<Eigen::Index>(N));

        Eigen::ArrayXd probs = (lw - lnZ).exp();
        E_avg  = (probs * E).sum();
        E2_avg = (probs * E * E).sum();
    } else
    {
        for (std::size_t i = 0; i < N; ++i) {
            double p_i = std::exp(log_w[i] - lnZ);
            double Ei  = energies_vec[i];
            E_avg  += p_i * Ei;
            E2_avg += p_i * Ei * Ei;
        }
    }

    double kT  = kB_kcal * T_;
    double var = E2_avg - E_avg * E_avg;

    Thermodynamics th;
    th.temperature    = T_;
    th.log_Z          = lnZ;
    th.free_energy    = -kT * lnZ;
    th.mean_energy    = E_avg;
    th.mean_energy_sq = E2_avg;
    th.heat_capacity  = std::max(0.0, var) / (kB_kcal * T_ * T_);
    th.entropy        = (E_avg - th.free_energy) / T_;
    th.std_energy     = std::sqrt(std::max(0.0, var));
    th.provenance     = provenance_;
    return th;
}

// ─── compute_at_temperature ───────────────────────────────────────────────────
// Re-evaluate the current ensemble at T_K without touching the stored T_ / beta_.
// Mirrors the Eigen branch of compute() exactly; only beta differs.

Thermodynamics StatMechEngine::compute_at_temperature(double T_K) const
{
    if (ensemble_.empty())
        throw std::runtime_error("StatMechEngine::compute_at_temperature: empty ensemble");
    if (T_K <= 0.0)
        throw std::invalid_argument("StatMechEngine::compute_at_temperature: T_K must be > 0");

    const double beta_T = 1.0 / (kB_kcal * T_K);
    const std::size_t N = ensemble_.size();
    bool any_positive_multiplicity = false;
    for (std::size_t i = 0; i < N; ++i) {
        if (!std::isfinite(ensemble_[i].energy))
            throw std::runtime_error("StatMechEngine::compute_at_temperature: non-finite sample energy");
        if (std::isfinite(ensemble_[i].count) && ensemble_[i].count > 0.0) {
            any_positive_multiplicity = true;
        }
    }
    if (!any_positive_multiplicity)
        throw std::runtime_error("StatMechEngine::compute_at_temperature: no positive-multiplicity samples");

    // log-weights at the requested temperature
    std::vector<double> log_w(N);
    {
        Eigen::ArrayXd counts(static_cast<Eigen::Index>(N));
        Eigen::ArrayXd energies(static_cast<Eigen::Index>(N));
        for (std::size_t i = 0; i < N; ++i) {
            double n = ensemble_[i].count;
            if (!std::isfinite(n) || n < 0.0) n = 0.0;
            counts(static_cast<Eigen::Index>(i))   = n;
            energies(static_cast<Eigen::Index>(i)) = ensemble_[i].energy;
        }
        Eigen::Map<Eigen::ArrayXd>(log_w.data(), static_cast<Eigen::Index>(N)) =
            counts.log() - beta_T * energies;
    }

    const double lnZ = log_sum_exp(log_w);

    // ⟨E⟩ and ⟨E²⟩ via Eigen vectorised path
    double E_avg = 0.0, E2_avg = 0.0;
    {
        Eigen::Map<const Eigen::ArrayXd> lw(log_w.data(), static_cast<Eigen::Index>(N));
        std::vector<double> ev(N);
        for (std::size_t i = 0; i < N; ++i) ev[i] = ensemble_[i].energy;
        Eigen::Map<const Eigen::ArrayXd> E(ev.data(), static_cast<Eigen::Index>(N));

        const Eigen::ArrayXd probs = (lw - lnZ).exp();
        E_avg  = (probs * E).sum();
        E2_avg = (probs * E * E).sum();
    }

    const double kT  = kB_kcal * T_K;
    const double var = E2_avg - E_avg * E_avg;

    Thermodynamics th;
    th.temperature    = T_K;
    th.log_Z          = lnZ;
    th.free_energy    = -kT * lnZ;
    th.mean_energy    = E_avg;
    th.mean_energy_sq = E2_avg;
    th.heat_capacity  = std::max(0.0, var) / (kB_kcal * T_K * T_K);
    th.entropy        = (E_avg - th.free_energy) / T_K;
    th.std_energy     = std::sqrt(std::max(0.0, var));
    th.provenance     = provenance_;
    return th;
}

ThermodynamicBreakdown StatMechEngine::compute_breakdown(
    double G_vib_kcal_mol,
    double G_natural_kcal_mol,
    double G_other_kcal_mol,
    bool has_vib,
    bool has_natural,
    bool has_other) const
{
    const Thermodynamics th = compute();

    ThermodynamicBreakdown b;
    b.temperature_K = th.temperature;
    b.logZ_config = th.log_Z;
    b.G_config_kcal_mol = th.free_energy;
    b.H_eff_kcal_mol = th.mean_energy;
    b.S_config_kcal_mol_K = th.entropy;
    b.minus_T_S_config_kcal_mol = th.free_energy - th.mean_energy;
    b.Cv_kcal_mol_K = th.heat_capacity;
    b.sigma_E_kcal_mol = th.std_energy;
    b.provenance = provenance_for_breakdown(
        th.provenance,
        G_vib_kcal_mol,
        G_natural_kcal_mol,
        G_other_kcal_mol,
        has_vib,
        has_natural,
        has_other);
    b.G_vib_kcal_mol = G_vib_kcal_mol;
    b.G_natural_kcal_mol = G_natural_kcal_mol;
    b.G_other_kcal_mol = G_other_kcal_mol;
    b.G_total_kcal_mol = b.G_config_kcal_mol
                       + b.G_vib_kcal_mol
                       + b.G_natural_kcal_mol
                       + b.G_other_kcal_mol;
    b.has_vib = has_vib;
    b.has_natural = has_natural;
    b.has_other = has_other;

    // I_E-E index (Williams et al. 2017) — diagnostic, never for ranking
    fill_IEE(b);

    return b;
}

// ─── compute_delta_Cp ────────────────────────────────────────────────────────
// Central finite-difference ΔCp of binding:
//   ΔCp ≈ [ΔH(T+dT) − ΔH(T−dT)] / (2·dT)
//
// Consistency check via entropy path:
//   ΔCp ≈ T_ref × [ΔS(T+dT) − ΔS(T−dT)] / (2·dT)
//
// Both routes use the SAME GA ensemble re-evaluated at T±dT.
// No additional sampling is required.

DeltaCpResult compute_delta_Cp(
    const StatMechEngine& bound,
    const StatMechEngine& unbound,
    double T_ref_K,
    double dT_K)
{
    if (dT_K <= 0.0)
        throw std::invalid_argument("compute_delta_Cp: dT_K must be > 0");
    if (T_ref_K - dT_K <= 0.0)
        throw std::invalid_argument("compute_delta_Cp: T_ref_K - dT_K must be > 0");

    const Thermodynamics bnd_lo  = bound.compute_at_temperature(T_ref_K - dT_K);
    const Thermodynamics bnd_hi  = bound.compute_at_temperature(T_ref_K + dT_K);
    const Thermodynamics ref_lo  = unbound.compute_at_temperature(T_ref_K - dT_K);
    const Thermodynamics ref_hi  = unbound.compute_at_temperature(T_ref_K + dT_K);

    const double dH_lo = bnd_lo.mean_energy - ref_lo.mean_energy;
    const double dH_hi = bnd_hi.mean_energy - ref_hi.mean_energy;
    const double dS_lo = bnd_lo.entropy     - ref_lo.entropy;
    const double dS_hi = bnd_hi.entropy     - ref_hi.entropy;

    DeltaCpResult r;
    r.T_ref_K              = T_ref_K;
    r.dT_K                 = dT_K;
    r.delta_H_lo           = dH_lo;
    r.delta_H_hi           = dH_hi;
    r.delta_S_lo           = dS_lo;
    r.delta_S_hi           = dS_hi;
    r.delta_Cp             = (dH_hi - dH_lo) / (2.0 * dT_K);
    r.delta_Cp_from_entropy = T_ref_K * (dS_hi - dS_lo) / (2.0 * dT_K);

    // Fractional consistency between the two ΔCp estimates
    const double abs_mean = 0.5 * (std::abs(r.delta_Cp) + std::abs(r.delta_Cp_from_entropy));
    r.consistency_check = std::abs(r.delta_Cp - r.delta_Cp_from_entropy)
                          / (abs_mean + 1e-9);
    r.consistent = (r.consistency_check < 0.05);

    return r;
}

ComponentAverages StatMechEngine::component_averages(
    std::span<const EnergyComponents> components) const
{
    if (components.empty())
        throw std::invalid_argument("StatMechEngine::component_averages: empty component list");
    if (components.size() != ensemble_.size())
        throw std::invalid_argument("StatMechEngine::component_averages: component count must match ensemble size");

    const std::vector<double> weights = boltzmann_weights();
    ComponentAverages avg;
    avg.component_completeness_flag = true;

    for (std::size_t i = 0; i < components.size(); ++i) {
        const double p = weights[i];
        const EnergyComponents& c = components[i];
        avg.mean_CF_kcal_mol += p * c.cf;
        avg.mean_receptor_strain_kcal_mol += p * c.receptor_strain;
        avg.mean_ligand_internal_kcal_mol += p * c.ligand_internal;
        avg.mean_hbond_kcal_mol += p * c.hbond;
        avg.mean_gist_kcal_mol += p * c.gist;
        avg.mean_metal_kcal_mol += p * c.metal;
        avg.mean_water_kcal_mol += p * c.water;
        avg.mean_other_kcal_mol += p * c.other;
        avg.component_completeness_flag = avg.component_completeness_flag && c.complete;
    }

    avg.component_sum_kcal_mol = avg.mean_CF_kcal_mol
        + avg.mean_receptor_strain_kcal_mol
        + avg.mean_ligand_internal_kcal_mol
        + avg.mean_hbond_kcal_mol
        + avg.mean_gist_kcal_mol
        + avg.mean_metal_kcal_mol
        + avg.mean_water_kcal_mol
        + avg.mean_other_kcal_mol;
    avg.component_status = avg.component_completeness_flag
        ? ComponentStatus::Available
        : ComponentStatus::IncludedInOther;
    return avg;
}

// ─── boltzmann_weights ───────────────────────────────────────────────────────

std::vector<double> StatMechEngine::boltzmann_weights() const {
    if (ensemble_.empty()) return {};

    const std::size_t N = ensemble_.size();

    // Build raw energy array and use the unified dispatch layer.
    std::vector<double> energies(N);
    for (std::size_t i = 0; i < N; ++i)
        energies[i] = ensemble_[i].energy;

    auto result = flexaids::compute_boltzmann_batch(energies, beta_);

    // Normalise weights accounting for multiplicities.
    std::vector<double> w(N);
    double Z_with_mult = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        Z_with_mult += ensemble_[i].count * result.weights[i];

    if (Z_with_mult > 0.0) {
        for (std::size_t i = 0; i < N; ++i)
            w[i] = ensemble_[i].count * result.weights[i] / Z_with_mult;
    }
    return w;
}

// ─── selection_weights ───────────────────────────────────────────────────────
// Identical to boltzmann_weights() but uses β_sel = 1/T (NOT the kB-folded
// physical β). For GA/cluster selection only — see header for the rationale.

std::vector<double> StatMechEngine::selection_weights() const {
    if (ensemble_.empty()) return {};

    const std::size_t N = ensemble_.size();

    std::vector<double> energies(N);
    for (std::size_t i = 0; i < N; ++i)
        energies[i] = ensemble_[i].energy;

    auto result = flexaids::compute_boltzmann_batch(energies, beta_selection_);

    std::vector<double> w(N);
    double Z_with_mult = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        Z_with_mult += ensemble_[i].count * result.weights[i];

    if (Z_with_mult > 0.0) {
        for (std::size_t i = 0; i < N; ++i)
            w[i] = ensemble_[i].count * result.weights[i] / Z_with_mult;
    }
    return w;
}

// ─── delta_G ─────────────────────────────────────────────────────────────────
// ΔG = F_this − F_ref = −kT (ln Z_this − ln Z_ref)

double StatMechEngine::delta_G(const StatMechEngine& reference) const {
    auto this_th = this->compute();
    auto ref_th  = reference.compute();
    double kT = kB_kcal * T_;
    return -kT * (this_th.log_Z - ref_th.log_Z);
}

// ─── Helmholtz convenience ───────────────────────────────────────────────────

double StatMechEngine::helmholtz(std::span<const double> energies, double T) {
    if (energies.empty())
        throw std::invalid_argument("helmholtz: empty energy list");
    double beta = 1.0 / (kB_kcal * T);

    // Use unified dispatch for the Boltzmann batch computation.
    auto result = flexaids::compute_boltzmann_batch(energies, beta);
    // F = -kT * ln(Z) where log_Z already accounts for E_min shift.
    return -(kB_kcal * T) * result.log_Z;
}

// ─── replica exchange ────────────────────────────────────────────────────────

std::vector<Replica>
StatMechEngine::init_replicas(std::span<const double> temperatures) {
    std::vector<Replica> reps;
    reps.reserve(temperatures.size());
    int id = 0;
    for (double T : temperatures) {
        Replica r;
        r.id             = id++;
        r.temperature    = T;
        r.beta           = 1.0 / (kB_kcal * T);
        r.current_energy = 0.0;
        reps.push_back(r);
    }
    return reps;
}

bool StatMechEngine::attempt_swap(Replica& a, Replica& b, std::mt19937& rng) {
    // Metropolis criterion:
    //   Δ = (β_a − β_b)(E_a − E_b)
    //   P_accept = min(1, exp(Δ))
    double delta = (a.beta - b.beta) * (a.current_energy - b.current_energy);
    if (delta >= 0.0) {
        std::swap(a.current_energy, b.current_energy);
        return true;
    }
    std::uniform_real_distribution<double> U(0.0, 1.0);
    if (U(rng) < std::exp(delta)) {
        std::swap(a.current_energy, b.current_energy);
        return true;
    }
    return false;
}

// ─── Boltzmann-reweighted PMF ────────────────────────────────────────────────
// Single-window post-hoc reweighting of an ensemble onto a 1D collective
// coordinate. NOT multi-window WHAM (Kumar et al. 1992) — that requires
// biased simulations with per-window offsets, neither of which are
// available here. The historical name `wham()` survives as a deprecated
// alias in the header (see statmech.h).

std::vector<WHAMBin> StatMechEngine::boltzmann_pmf(
    std::span<const double> energies,
    std::span<const double> coordinates,
    double temperature,
    int    n_bins,
    int    max_iter,
    double tolerance)
{
    if (energies.size() != coordinates.size())
        throw std::invalid_argument("boltzmann_pmf: energies and coordinates size mismatch");
    if (energies.empty() || n_bins <= 0)
        throw std::invalid_argument("boltzmann_pmf: invalid input");

    const std::size_t N = energies.size();
    double beta = 1.0 / (kB_kcal * temperature);

    // Find coordinate range (single pass)
    auto [cmin_it, cmax_it] = std::minmax_element(coordinates.begin(), coordinates.end());
    double cmin = *cmin_it;
    double cmax = *cmax_it;
    double bin_w = (cmax - cmin) / n_bins;
    if (bin_w <= 0.0) bin_w = 1.0;

    // Histogram + Boltzmann-weighted histogram.
    // Use the dispatch layer for the Boltzmann weight computation,
    // then bin the pre-computed weights for O(N) histogramming.
    auto boltz_result = flexaids::compute_boltzmann_batch(energies, beta);

    std::vector<double> raw_count(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> boltz_sum(static_cast<std::size_t>(n_bins), 0.0);
    [[maybe_unused]] double inv_bw = 1.0 / bin_w;

#ifdef _OPENMP
    // OpenMP parallel histogram with per-thread private bins
    if (N >= OMP_THRESHOLD) {
        int n_threads = omp_get_max_threads();
        std::vector<std::vector<double>> t_raw(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));
        std::vector<std::vector<double>> t_boltz(n_threads,
            std::vector<double>(static_cast<std::size_t>(n_bins), 0.0));

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(N); ++i) {
            int tid = omp_get_thread_num();
            int b = static_cast<int>((coordinates[i] - cmin) * inv_bw);
            b = std::min(std::max(b, 0), n_bins - 1);
            t_raw[tid][static_cast<std::size_t>(b)]  += 1.0;
            t_boltz[tid][static_cast<std::size_t>(b)] += boltz_result.weights[i];
        }
        // Reduce thread-private histograms
        for (auto& tr : t_raw)
            for (int b = 0; b < n_bins; ++b)
                raw_count[static_cast<std::size_t>(b)] += tr[static_cast<std::size_t>(b)];
        for (auto& tb : t_boltz)
            for (int b = 0; b < n_bins; ++b)
                boltz_sum[static_cast<std::size_t>(b)] += tb[static_cast<std::size_t>(b)];
    } else
#endif
    {
        for (std::size_t i = 0; i < N; ++i) {
            int b = static_cast<int>((coordinates[i] - cmin) / bin_w);
            if (b < 0) b = 0;
            if (b >= n_bins) b = n_bins - 1;
            raw_count[static_cast<std::size_t>(b)] += 1.0;
            boltz_sum[static_cast<std::size_t>(b)] += boltz_result.weights[i];
        }
    }

    // Free energy per bin: F_b = −kT ln( weighted_count_b / raw_count_b )
    // Iterative self-consistency (single-window simplification)
    std::vector<double> f_old(static_cast<std::size_t>(n_bins), 0.0);
    std::vector<double> f_new(static_cast<std::size_t>(n_bins), 0.0);

    for (int iter = 0; iter < max_iter; ++iter) {
        // Eigen vectorised WHAM self-consistency update
        Eigen::Map<const Eigen::ArrayXd> rc(raw_count.data(), n_bins);
        Eigen::Map<const Eigen::ArrayXd> bs(boltz_sum.data(), n_bins);
        Eigen::Map<Eigen::ArrayXd> fn(f_new.data(), n_bins);
        Eigen::Map<Eigen::ArrayXd> fo(f_old.data(), n_bins);

        // F_b = -kT * ln(boltz_sum_b / raw_count_b) where raw_count_b > 0
        auto occupied = (rc > 0.0);
        Eigen::ArrayXd safe_rc = occupied.select(rc, Eigen::ArrayXd::Ones(n_bins));
        fn = occupied.select(
            -(kB_kcal * temperature) * (bs / safe_rc).log(),
            Eigen::ArrayXd::Zero(n_bins));

        // Shift so minimum = 0
        fn -= fn.minCoeff();

        // Check convergence
        {
            double maxdiff = (fn - fo).abs().maxCoeff();
            fo = fn;
            if (maxdiff < tolerance) break;
        }
    }

    // Build output
    std::vector<WHAMBin> result(static_cast<std::size_t>(n_bins));
    for (int b = 0; b < n_bins; ++b) {
        result[static_cast<std::size_t>(b)].coord_center = cmin + (b + 0.5) * bin_w;
        result[static_cast<std::size_t>(b)].count        = raw_count[static_cast<std::size_t>(b)];
        result[static_cast<std::size_t>(b)].free_energy  = f_new[static_cast<std::size_t>(b)];
    }
    return result;
}

// ─── thermodynamic integration ───────────────────────────────────────────────
// ΔG = ∫₀¹ ⟨∂V/∂λ⟩_λ dλ   (trapezoidal rule)

double StatMechEngine::thermodynamic_integration(std::span<const TIPoint> points) {
    if (points.size() < 2)
        throw std::invalid_argument("TI requires at least 2 points");

    double integral = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        double dl = points[i].lambda - points[i-1].lambda;
        integral += 0.5 * dl * (points[i].dV_dlambda + points[i-1].dV_dlambda);
    }
    return integral;
}

// ─── BoltzmannLUT ────────────────────────────────────────────────────────────

BoltzmannLUT::BoltzmannLUT(double beta, double e_min, double e_max, int n_bins)
    : beta_(beta)
    , e_min_(e_min)
    , n_bins_(n_bins)
    , table_(static_cast<std::size_t>(n_bins))
{
    double range = e_max - e_min;
    if (range <= 0.0) range = 1.0;
    inv_bin_width_ = n_bins / range;

    // Eigen-vectorised LUT initialisation
    {
        Eigen::ArrayXd idx = Eigen::ArrayXd::LinSpaced(n_bins, 0.5, n_bins - 0.5);
        Eigen::ArrayXd E = e_min + idx * (range / n_bins);
        Eigen::Map<Eigen::ArrayXd>(table_.data(), n_bins) = (-beta * E).exp();
    }
}

double BoltzmannLUT::operator()(double energy) const noexcept {
    int idx = static_cast<int>((energy - e_min_) * inv_bin_width_);
    if (idx < 0) idx = 0;
    if (idx >= n_bins_) idx = n_bins_ - 1;
    return table_[static_cast<std::size_t>(idx)];
}

// ─── ensemble merging (parallel grid-decomposed docking) ─────────────────────

void StatMechEngine::merge(const StatMechEngine& other) {
    if (std::fabs(other.T_ - T_) > 1e-6)
        throw std::invalid_argument("Cannot merge engines at different temperatures");
    ensemble_.insert(ensemble_.end(),
                     other.ensemble_.begin(), other.ensemble_.end());

    // A merged physical claim is safe only when both source ensembles carry
    // exactly the same provenance witness. Numeric merging is unchanged; a
    // mismatch only downgrades interpretation to the fail-closed default.
    if (provenance_ != other.provenance_)
        provenance_ = ScientificProvenance{};
}

void StatMechEngine::merge_samples(std::span<const double> energies,
                                    std::span<const double> multiplicities) {
    if (energies.size() != multiplicities.size())
        throw std::invalid_argument("energies and multiplicities must have same size");
    for (size_t i = 0; i < energies.size(); ++i)
        ensemble_.push_back({energies[i], multiplicities[i]});

    // Raw arrays carry no provenance witness. merge(const StatMechEngine&)
    // can compare two witnesses and downgrade on mismatch; this transport
    // (MPI / deserialization) has nothing to compare against, so a receiving
    // engine cannot attest that the injected samples came from the calibrated
    // source its own witness describes. Fail closed rather than let unattested
    // energies inherit an authorizing witness.
    //
    // Deliberately scoped to engines that could actually authorize a claim: an
    // already-proxy-only witness keeps its descriptive domain/measure strings,
    // so every current CF/contact-function path — including ParallelDock's
    // aggregation — is byte-for-byte unaffected in both numerics and emitted
    // metadata. Callers that can independently attest the transported samples
    // must re-declare the witness explicitly via set_provenance().
    if (!energies.empty() && !provenance_.is_proxy_only())
        provenance_ = ScientificProvenance{};
}

std::vector<double> StatMechEngine::serialize_energies() const {
    std::vector<double> out(ensemble_.size());
    for (size_t i = 0; i < ensemble_.size(); ++i)
        out[i] = ensemble_[i].energy;
    return out;
}

std::vector<double> StatMechEngine::serialize_multiplicities() const {
    std::vector<double> out(ensemble_.size());
    for (size_t i = 0; i < ensemble_.size(); ++i)
        out[i] = ensemble_[i].count;
    return out;
}

// ─── make_breakdown (Task 1 ledger) ──────────────────────────────────────────
// Derives the full audited ThermodynamicBreakdown from a live engine.
// All identities from thermo_invariants.md are enforced by construction here
// (G_config = -kT logZ, S = (H-G)/T, minus_TS = G-H, G_total = sum of parts).
// Corrections are passed in by the caller (BindingMode for vib/natural, etc.).
// No ranking side-effects. Safe for use in tests and future JSON paths.
ThermodynamicBreakdown StatMechEngine::make_breakdown(
    const StatMechEngine& engine,
    double G_vib_kcal_mol,     bool has_vib,
    double G_natural_kcal_mol, bool has_natural,
    double G_other_kcal_mol,   bool has_other)
{
    ThermodynamicBreakdown b;
    // Reuse compute()'s fail-closed empty-ensemble contract. A zero-filled
    // ledger is not a valid substitute for missing thermodynamic moments.
    const auto th = engine.compute();

    b.temperature_K = th.temperature;

    b.logZ_config = th.log_Z;
    b.G_config_kcal_mol = th.free_energy;
    b.H_eff_kcal_mol = th.mean_energy;
    b.S_config_kcal_mol_K = th.entropy;
    b.minus_T_S_config_kcal_mol = th.free_energy - th.mean_energy;
    b.Cv_kcal_mol_K = th.heat_capacity;
    b.sigma_E_kcal_mol = th.std_energy;
    b.provenance = provenance_for_breakdown(
        th.provenance,
        G_vib_kcal_mol,
        G_natural_kcal_mol,
        G_other_kcal_mol,
        has_vib,
        has_natural,
        has_other);

    b.G_vib_kcal_mol = G_vib_kcal_mol;
    b.G_natural_kcal_mol = G_natural_kcal_mol;
    b.G_other_kcal_mol = G_other_kcal_mol;

    b.G_total_kcal_mol = th.free_energy + G_vib_kcal_mol + G_natural_kcal_mol + G_other_kcal_mol;

    b.has_vib = has_vib;
    b.has_natural = has_natural;
    b.has_other = has_other;

    return b;
}

// ─── Task 3: Component-wise weighted averages ────────────────────────────────

EnergyComponents StatMechEngine::compute_weighted_components(
    std::span<const double> weights,
    std::span<const EnergyComponents> components)
{
    const size_t n = weights.size();
    if (n == 0 || n != components.size())
        return {};

    double sum_w = 0.0;
    EnergyComponents result{};

    for (size_t i = 0; i < n; ++i) {
        const double w = weights[i];
        sum_w += w;

        const auto& c = components[i];
        result.cf               += w * c.cf;
        result.receptor_strain  += w * c.receptor_strain;
        result.ligand_internal  += w * c.ligand_internal;
        result.hbond            += w * c.hbond;
        result.gist             += w * c.gist;
        result.metal            += w * c.metal;
        result.water            += w * c.water;
        result.other            += w * c.other;
        result.total            += w * c.total;
    }

    if (sum_w > 1e-300) {
        const double inv = 1.0 / sum_w;
        result.cf              *= inv;
        result.receptor_strain *= inv;
        result.ligand_internal *= inv;
        result.hbond           *= inv;
        result.gist            *= inv;
        result.metal           *= inv;
        result.water           *= inv;
        result.other           *= inv;
        result.total           *= inv;
    }

    // component_sum is the sum of the averaged pieces (diagnostic)
    // Note: this may legitimately differ from H_eff if not all energy was decomposed.
    return result;
}

ThermodynamicBreakdown StatMechEngine::make_breakdown_with_components(
    const StatMechEngine& engine,
    std::span<const EnergyComponents> components,
    double G_vib_kcal_mol,     bool has_vib,
    double G_natural_kcal_mol, bool has_natural,
    double G_other_kcal_mol,   bool has_other)
{
    ThermodynamicBreakdown b = make_breakdown(
        engine, G_vib_kcal_mol, has_vib,
        G_natural_kcal_mol, has_natural,
        G_other_kcal_mol, has_other);

    if (!components.empty() && components.size() == engine.size()) {
        auto weights = engine.boltzmann_weights();
        b.component_means = compute_weighted_components(weights, components);

        // Compute component_sum for convenience / diagnostics
        const auto& m = b.component_means;
        b.component_sum_kcal_mol =
            m.cf + m.receptor_strain + m.ligand_internal + m.hbond +
            m.gist + m.metal + m.water + m.other;

        // Heuristic completeness: if the two biggest terms (CF + strain) are
        // marked Available, we consider the decomposition "reasonably complete".
        b.components_complete =
            (m.cf_status == ComponentStatus::Available) &&
            (m.receptor_strain_status == ComponentStatus::Available ||
             m.receptor_strain_status == ComponentStatus::NotComputed); // allow single-conformer case
    }

    return b;
}

// ─── Task 4: Diagnostic metric implementations (on the ledger) ───────────────

double ThermodynamicBreakdown::entropy_fraction() const {
    return statmech::entropy_fraction(H_eff_kcal_mol, minus_T_S_config_kcal_mol);
}

double ThermodynamicBreakdown::enthalpy_fraction() const {
    return statmech::enthalpy_fraction(H_eff_kcal_mol, minus_T_S_config_kcal_mol);
}

double ThermodynamicBreakdown::compensation_score() const {
    return statmech::compensation_score(G_config_kcal_mol, H_eff_kcal_mol, minus_T_S_config_kcal_mol);
}

// ─── Task 5: Joint Receptor–Ligand Ensemble (EXPERIMENTAL) ───────────────────

namespace {

double safe_entropy(double p) {
    if (p <= 0.0) return 0.0;
    return -p * std::log(p);
}

} // anonymous

JointEnsembleResult StatMechEngine::compute_joint_ensemble(
    std::span<const JointMicrostate> microstates,
    double temperature_K)
{
    JointEnsembleResult result;
    result.temperature_K = temperature_K;
    result.experimental = true;

    if (microstates.empty()) {
        return result;
    }

    const double beta = 1.0 / (kB_kcal * temperature_K);

    // Step 1: Compute log-weights for all microstates
    const size_t N = microstates.size();
    std::vector<double> log_w(N);
    for (size_t i = 0; i < N; ++i) {
        const auto& m = microstates[i];
        log_w[i] = m.log_multiplicity - beta * m.energy.total;
    }

    double logZ = log_sum_exp(log_w);
    result.logZ = logZ;
    result.G_kcal_mol = -kB_kcal * temperature_K * logZ;

    // Step 2: Compute probabilities p(r,i) and accumulate for marginals + moments
    std::map<int, double> p_receptor;
    std::map<int, double> p_ligand;
    double H = 0.0;

    std::vector<double> p(N);

    for (size_t i = 0; i < N; ++i) {
        p[i] = std::exp(log_w[i] - logZ);
        const auto& m = microstates[i];

        H += p[i] * m.energy.total;

        p_receptor[m.receptor_conformer_id] += p[i];
        p_ligand[m.ligand_pose_id] += p[i];
    }

    result.H_kcal_mol = H;

    // Step 3: Convert maps to vectors (sorted by id for determinism)
    std::vector<int> receptor_ids;
    for (auto& kv : p_receptor) receptor_ids.push_back(kv.first);
    std::sort(receptor_ids.begin(), receptor_ids.end());

    std::vector<int> ligand_ids;
    for (auto& kv : p_ligand) ligand_ids.push_back(kv.first);
    std::sort(ligand_ids.begin(), ligand_ids.end());

    result.receptor_population.resize(receptor_ids.size());
    result.ligand_population.resize(ligand_ids.size());

    for (size_t r = 0; r < receptor_ids.size(); ++r) {
        result.receptor_population[r] = p_receptor[receptor_ids[r]];
    }
    for (size_t l = 0; l < ligand_ids.size(); ++l) {
        result.ligand_population[l] = p_ligand[ligand_ids[l]];
    }

    // Step 4: Entropies (in kcal/mol/K, using kB_kcal)
    double S_joint = 0.0;
    for (double prob : p) {
        S_joint += safe_entropy(prob);
    }
    result.S_joint_kcal_mol_K = kB_kcal * S_joint;

    double S_receptor = 0.0;
    for (double pr : result.receptor_population) {
        S_receptor += safe_entropy(pr);
    }
    result.S_receptor_kcal_mol_K = kB_kcal * S_receptor;

    double S_ligand = 0.0;
    for (double pi : result.ligand_population) {
        S_ligand += safe_entropy(pi);
    }
    result.S_ligand_kcal_mol_K = kB_kcal * S_ligand;

    // Step 5: Mutual information I(R;L) = S_R + S_L − S_joint (nats).
    // The previous S_joint − S_R − S_L sign was the negation of I.
    result.mutual_information_dimensionless = S_receptor + S_ligand - S_joint;
    if (result.mutual_information_dimensionless < 0.0)
        result.mutual_information_dimensionless = 0.0;

    // Fallback detection
    bool has_real_receptor_ids = false;
    for (const auto& m : microstates) {
        if (m.receptor_conformer_id >= 0) { has_real_receptor_ids = true; break; }
    }
    if (!has_real_receptor_ids) {
        result.fallback_single_receptor = true;
        result.S_receptor_kcal_mol_K = 0.0;
        result.mutual_information_dimensionless = 0.0;
    }

    return result;
}

// ─── Task 6: Standard-State Affinity Calibration (safe utilities) ────────────

double deltaG_standard_to_Kd_M(double deltaG_kcal_mol, double T_K, double c0_M) {
    if (T_K <= 0.0) {
        throw std::invalid_argument("Temperature must be > 0 K for affinity conversion");
    }
    if (c0_M <= 0.0) {
        throw std::invalid_argument("Standard state concentration c0_M must be > 0");
    }

    // ΔG° (kcal/mol) = RT ln(Kd / c0)   →   Kd (M) = c0 * exp(ΔG° / RT_kcal)
    const double RT_kcal = kB_kcal * T_K;
    return c0_M * std::exp(deltaG_kcal_mol / RT_kcal);
}

double Kd_M_to_deltaG_standard(double Kd_M, double T_K, double c0_M) {
    if (T_K <= 0.0) {
        throw std::invalid_argument("Temperature must be > 0 K for affinity conversion");
    }
    if (Kd_M <= 0.0) {
        throw std::invalid_argument("Kd must be > 0 M");
    }
    if (c0_M <= 0.0) {
        throw std::invalid_argument("Standard state concentration c0_M must be > 0");
    }

    const double RT_kcal = kB_kcal * T_K;
    // ΔG° = RT ln(Kd / c0)
    return RT_kcal * std::log(Kd_M / c0_M);
}

}  // namespace statmech
