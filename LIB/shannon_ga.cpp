/**
 * @file shannon_ga.cpp
 * @brief Implementation of the Shannon Thermodynamic GA
 * @author Le Bonhomme Pharma / FlexAIDdS
 */
#include "shannon_ga.h"
#include "RngSeed.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace FlexAID {
namespace StatMech {

// ─────────────────────────────────────────────────────────────────────────────
ShannonThermodynamicGA::ShannonThermodynamicGA(size_t pop_size,
                                               double mutation_rate,
                                               double temperature_K)
    : pop_size_(pop_size)
    , mutation_rate_(mutation_rate)
    , temperature_K_(temperature_K)
    , current_entropy_(0.0)
    , current_log_Z_(0.0)
    , engine_(temperature_K)
    , last_thermo_{}
    , rng_(flexaids_rng::seed_from_env_or_random(0x534841ULL))
{
    (void)pop_size_;  // reserved for future population-level heuristics
}

// ─────────────────────────────────────────────────────────────────────────────
void ShannonThermodynamicGA::evaluate_enthalpy_batch(chromosome* chrom,
                                                      int num_chrom)
{
    engine_.clear();

    #ifdef _OPENMP
    // Collect energies in parallel, add to engine sequentially
    // (StatMechEngine::add_sample is not thread-safe)
    std::vector<double> energies(num_chrom);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < num_chrom; ++i)
        energies[i] = chrom[i].evalue;

    for (int i = 0; i < num_chrom; ++i)
        engine_.add_sample(energies[i]);
    #else
    for (int i = 0; i < num_chrom; ++i)
        engine_.add_sample(chrom[i].evalue);
    #endif
}

// ─────────────────────────────────────────────────────────────────────────────
GenerationThermo ShannonThermodynamicGA::collapse_thermodynamics(
    chromosome* chrom, int num_chrom, int generation)
{
    // Ensure engine is populated
    if (engine_.size() == 0)
        evaluate_enthalpy_batch(chrom, num_chrom);

    // Compute full thermodynamics via LSE partition function
    last_thermo_ = engine_.compute();
    current_log_Z_ = last_thermo_.log_Z;

    // Get Boltzmann weights and write back onto chromosomes
    auto bweights = engine_.boltzmann_weights();
    for (int i = 0; i < num_chrom; ++i) {
        chrom[i].boltzmann_weight = bweights[i];
        chrom[i].free_energy = last_thermo_.free_energy;
    }

    // Compute Shannon entropy of the energy distribution
    std::vector<double> pop_energies(num_chrom);
    for (int i = 0; i < num_chrom; ++i)
        pop_energies[i] = chrom[i].evalue;

    double H = shannon_thermo::compute_shannon_entropy(
        pop_energies, shannon_thermo::DEFAULT_HIST_BINS);
    current_entropy_ = H;
    entropy_history_.push_back(H);

    // Check for plateau
    bool plateau = false;
    if (entropy_history_.size() >= 2) {
        // Use the existing detect_entropy_plateau with a minimal window
        int window = std::min(static_cast<int>(entropy_history_.size()), 5);
        plateau = shannon_thermo::detect_entropy_plateau(
            entropy_history_, window, 0.01);
    }

    GenerationThermo snap;
    snap.generation       = generation;
    snap.free_energy      = last_thermo_.free_energy;
    snap.mean_energy      = last_thermo_.mean_energy;
    snap.entropy          = last_thermo_.entropy;
    snap.heat_capacity    = last_thermo_.heat_capacity;
    snap.shannon_H        = H;
    snap.plateau_detected = plateau;

    return snap;
}

// ─────────────────────────────────────────────────────────────────────────────
int ShannonThermodynamicGA::select_parent(const chromosome* chrom,
                                           int num_chrom)
{
    // Boltzmann-weighted roulette wheel selection
    double total = 0.0;
    for (int i = 0; i < num_chrom; ++i)
        total += chrom[i].boltzmann_weight;

    if (total <= 0.0) {
        // Fallback: uniform random
        std::uniform_int_distribution<int> dist(0, num_chrom - 1);
        return dist(rng_);
    }

    std::uniform_real_distribution<double> dist(0.0, total);
    double r = dist(rng_);
    double cumulative = 0.0;
    for (int i = 0; i < num_chrom; ++i) {
        cumulative += chrom[i].boltzmann_weight;
        if (r <= cumulative)
            return i;
    }
    return num_chrom - 1;
}

// ─────────────────────────────────────────────────────────────────────────────
void ShannonThermodynamicGA::mutate_and_crossover(chromosome* chrom,
                                                    int num_chrom,
                                                    int num_genes,
                                                    const genlim* gene_lim)
{
    // Adaptive mutation: boost when Shannon entropy is low (population collapsed)
    // Normalized entropy: H / ln(num_bins) ∈ [0, 1]
    double H_max = std::log(static_cast<double>(shannon_thermo::DEFAULT_HIST_BINS));
    double H_norm = (H_max > 0.0) ? (current_entropy_ / H_max) : 1.0;

    // When H_norm → 0 (collapsed), adaptive_rate → 3× base rate (capped at 0.5)
    // When H_norm → 1 (diverse),  adaptive_rate → base rate
    double entropy_boost = 1.0 + 2.0 * (1.0 - H_norm);
    double adaptive_rate = std::min(mutation_rate_ * entropy_boost, 0.5);

    // Apply to bottom half of population (sorted by energy, worst individuals)
    int half = num_chrom / 2;
    for (int i = half; i < num_chrom; ++i) {
        // Select two parents via Boltzmann-weighted roulette
        int p1 = select_parent(chrom, half);  // select from top half
        int p2 = select_parent(chrom, half);
        while (p2 == p1 && half > 1)
            p2 = select_parent(chrom, half);

        // Crossover: copy genes from parents
        crossover(chrom[p1].genes, chrom[p2].genes, num_genes, 0);

        // Copy result into offspring slot
        for (int g = 0; g < num_genes; ++g)
            chrom[i].genes[g] = chrom[p1].genes[g];

        // Mutate with entropy-adaptive rate
        mutate(chrom[i].genes, num_genes, adaptive_rate);

        // Decode genes
        for (int g = 0; g < num_genes; ++g)
            chrom[i].genes[g].to_ic = genetoic(&gene_lim[g],
                                                 chrom[i].genes[g].to_int32);

        chrom[i].status = 'o';  // mark for re-evaluation
    }
}

// ─────────────────────────────────────────────────────────────────────────────
CollapseResult ShannonThermodynamicGA::run(chromosome* chrom, int num_chrom,
                                            int num_genes,
                                            const genlim* gene_lim,
                                            int max_generations,
                                            double entropy_tolerance,
                                            int check_interval,
                                            int plateau_window)
{
    CollapseResult result;
    result.converged = false;
    result.final_generation = max_generations;

    entropy_history_.clear();

    for (int gen = 0; gen < max_generations; ++gen) {

        // Sort population by energy (ascending = best first)
        QuickSort(chrom, 0, num_chrom - 1, true);

        // Periodically compute thermodynamics and check convergence
        if ((gen + 1) % check_interval == 0) {

            // Feed energies into engine and compute thermo
            evaluate_enthalpy_batch(chrom, num_chrom);
            GenerationThermo snap = collapse_thermodynamics(
                chrom, num_chrom, gen + 1);
            result.history.push_back(snap);

            std::fprintf(stderr,
                "[ShannonGA] gen=%d  F=%.3f  <E>=%.3f  S=%.6f  H=%.4f\n",
                gen + 1, snap.free_energy, snap.mean_energy,
                snap.entropy, snap.shannon_H);

            if (snap.shannon_H <= shannon_thermo::kHSC_soft_nats) {
                std::fprintf(stderr,
                    "[ShannonGA] Entropy collapse at generation %d "
                    "(H=%.4f nats <= %.4f nats / %.1f bits)\n",
                    gen + 1, snap.shannon_H,
                    shannon_thermo::kHSC_soft_nats,
                    shannon_thermo::kHSC_soft_bits);
                result.converged = true;
                result.final_generation = gen + 1;
                break;
            }

            // Check for entropy plateau
            if (static_cast<int>(entropy_history_.size()) >= plateau_window) {
                bool plateau = shannon_thermo::detect_entropy_plateau(
                    entropy_history_, plateau_window, entropy_tolerance);
                if (plateau) {
                    std::fprintf(stderr,
                        "[ShannonGA] Entropy plateau at generation %d "
                        "(H=%.4f, tol=%.4f)\n",
                        gen + 1, snap.shannon_H, entropy_tolerance);
                    result.converged = true;
                    result.final_generation = gen + 1;
                    break;
                }
            }
        }

        // Apply entropy-aware GA operators
        mutate_and_crossover(chrom, num_chrom, num_genes, gene_lim);
    }

    // Final thermodynamic snapshot
    if (!result.history.empty()) {
        const auto& last = result.history.back();
        result.final_free_energy = last.free_energy;
        result.final_entropy     = last.entropy;
        result.final_shannon_H   = last.shannon_H;
    } else {
        // Compute once if we never hit a check interval
        evaluate_enthalpy_batch(chrom, num_chrom);
        auto snap = collapse_thermodynamics(chrom, num_chrom, max_generations);
        result.history.push_back(snap);
        result.final_free_energy = snap.free_energy;
        result.final_entropy     = snap.entropy;
        result.final_shannon_H   = snap.shannon_H;
    }

    return result;
}

} // namespace StatMech
} // namespace FlexAID
