/**
 * @file shannon_ga.h
 * @brief Shannon Thermodynamic Genetic Algorithm Core
 * @author Le Bonhomme Pharma / FlexAIDdS
 *
 * Unified wrapper that combines the StatMechEngine (LSE partition function),
 * ShannonThermoStack (configurational entropy), and GA diversity monitoring
 * into a single entropy-collapse GA driver.
 *
 * This does NOT replace gaboom.cpp; it provides a reusable component that
 * gaboom's main loop can delegate to for thermodynamic-aware generations.
 */
#pragma once

#include "statmech.h"
#include "ShannonThermoStack/ShannonThermoStack.h"
#include "gaboom.h"

#include <vector>
#include <random>
#include <functional>

namespace FlexAID {
namespace StatMech {

// ─── per-generation thermodynamic snapshot ───────────────────────────────────
struct GenerationThermo {
    int    generation;
    double free_energy;       // F = -kT ln Z  (kcal/mol)
    double mean_energy;       // <E>
    double entropy;           // S = (<E> - F) / T
    double heat_capacity;     // C_v
    double shannon_H;         // Shannon entropy of energy distribution (nats)
    bool   plateau_detected;  // true if entropy plateau triggered
};

// ─── convergence result ──────────────────────────────────────────────────────
struct CollapseResult {
    bool   converged;             // true if entropy collapse threshold or plateau reached
    int    final_generation;      // generation at convergence (or max)
    double final_free_energy;     // F at termination
    double final_entropy;         // S at termination
    double final_shannon_H;       // Shannon H at termination
    std::vector<GenerationThermo> history;  // per-check-interval snapshots
};

// ─── main class ──────────────────────────────────────────────────────────────

class ShannonThermodynamicGA {
public:
    /**
     * @param pop_size       GA population size
     * @param mutation_rate  base mutation probability
     * @param temperature_K  thermodynamic temperature (default 298.15 K)
     */
    ShannonThermodynamicGA(size_t pop_size,
                           double mutation_rate,
                           double temperature_K = 298.15);

    // ── batch enthalpy evaluation ────────────────────────────────────────────
    // Reads evalue from each chromosome and feeds the StatMechEngine.
    void evaluate_enthalpy_batch(chromosome* chrom, int num_chrom);

    // ── thermodynamic collapse step ──────────────────────────────────────────
    // Computes partition function (LSE), Boltzmann weights, Shannon entropy,
    // and writes boltzmann_weight / free_energy back onto chromosomes.
    // Returns the GenerationThermo snapshot for this generation.
    GenerationThermo collapse_thermodynamics(chromosome* chrom, int num_chrom,
                                             int generation);

    // ── entropy-weighted selection ────────────────────────────────────────────
    // Selects a parent index via Boltzmann-weighted roulette wheel.
    int select_parent(const chromosome* chrom, int num_chrom);

    // ── GA operators with entropy feedback ───────────────────────────────────
    // Applies crossover + mutation with adaptive rates based on current entropy.
    // When entropy is low (population collapsed), mutation rate is boosted.
    void mutate_and_crossover(chromosome* chrom, int num_chrom,
                              int num_genes, const genlim* gene_lim);

    // ── full run (standalone driver) ─────────────────────────────────────────
    // Runs the entropy-collapse GA loop for up to max_generations.
    // Stops early if Shannon entropy crosses the HSC nats threshold or plateaus
    // within entropy_tolerance.
    // The caller must have already populated and evaluated the initial population.
    CollapseResult run(chromosome* chrom, int num_chrom,
                       int num_genes, const genlim* gene_lim,
                       int max_generations,
                       double entropy_tolerance = 0.01,
                       int check_interval = 10,
                       int plateau_window = 5);

    // ── accessors ────────────────────────────────────────────────────────────
    double current_entropy()     const noexcept { return current_entropy_; }
    double current_partition_Z() const noexcept { return current_log_Z_; }
    double temperature()         const noexcept { return temperature_K_; }
    const statmech::Thermodynamics& last_thermo() const noexcept { return last_thermo_; }

private:
    [[maybe_unused]] size_t pop_size_;
    double mutation_rate_;
    double temperature_K_;
    double current_entropy_;
    double current_log_Z_;
    statmech::StatMechEngine engine_;
    statmech::Thermodynamics last_thermo_;
    std::vector<double> entropy_history_;
    std::mt19937 rng_;
};

} // namespace StatMech
} // namespace FlexAID
