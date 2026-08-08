// ReferenceEntropy.cpp — Reference state entropy corrections
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.

#include "ReferenceEntropy.h"
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>
#include <cstdio>
#include <span>

#include <Eigen/Dense>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace reference_entropy {

// ─── Receptor conformational entropy from model ensemble ─────────────────────
ReceptorEntropyResult compute_receptor_entropy(
    const std::vector<double>& model_energies,
    double temperature_K)
{
    ReceptorEntropyResult r;
    r.n_models = static_cast<int>(model_energies.size());

    if (r.n_models <= 1) {
        r.S_conf = 0.0;
        r.T_S_conf = 0.0;
        r.energy_spread = 0.0;
        return r;
    }

    // Compute Boltzmann weights from model energies
    // Eigen-accelerated path: vectorized exp/log operations
    const double beta = 1.0 / (kB_kcal * temperature_K);

    // Eigen vectorized: all exp/log/sum operations auto-vectorize to AVX/SSE
    Eigen::Map<const Eigen::ArrayXd> E(model_energies.data(), r.n_models);
    const double e_min = E.minCoeff();
    const double e_max = E.maxCoeff();
    r.energy_spread = e_max - e_min;

    // Shifted exponents for log-sum-exp stability
    Eigen::ArrayXd shifted = -beta * (E - e_min);
    const double max_s = shifted.maxCoeff();
    const double log_Z = max_s + std::log((shifted - max_s).exp().sum());

    // Boltzmann probabilities and Shannon entropy (vectorized)
    Eigen::ArrayXd log_pi = shifted - log_Z;
    Eigen::ArrayXd pi = log_pi.exp();
    // H = -sum(p_i * log(p_i)), mask zeros
    Eigen::ArrayXd safe_log = (pi > 1e-15).select(log_pi, Eigen::ArrayXd::Zero(r.n_models));
    Eigen::ArrayXd safe_pi  = (pi > 1e-15).select(pi, Eigen::ArrayXd::Zero(r.n_models));
    double H = -(safe_pi * safe_log).sum();

    // Convert from nats to kcal/mol·K via: S = kB * H
    r.S_conf = kB_kcal * H;
    r.T_S_conf = -temperature_K * r.S_conf;

    return r;
}

// ─── Ligand solution entropy ─────────────────────────────────────────────────
LigandEntropyResult compute_ligand_solution_entropy(
    int n_rotatable_bonds,
    double temperature_K)
{
    LigandEntropyResult r;
    r.n_rotatable = n_rotatable_bonds;

    // Each rotatable bond contributes ~R·ln(3) to conformational entropy
    // (3-fold rotational barrier approximation)
    r.S_rot = n_rotatable_bonds * R_kcal * std::log(3.0);
    r.T_S_rot = -temperature_K * r.S_rot;

    return r;
}

// ─── Cratic entropy ─────────────────────────────────────────────────────────
CraticResult compute_cratic_entropy(double temperature_K)
{
    CraticResult r;

    // Standard concentration: 1M = 1 mol/L = 6.022e23 / L
    // Accessible volume per molecule at 1M: V = 1/N_A = 1.66e-27 m³ = 1660 Å³
    // ΔS_cratic = -R·ln(V_site/V_standard) ≈ -R·ln(1/1660) = R·ln(1660)
    r.S_cratic = -R_kcal * std::log(1.0 / 1660.0);  // positive (entropy loss)
    r.T_S_cratic = -temperature_K * (-r.S_cratic);   // positive penalty in ΔG

    return r;
}

// ─── Combined reference correction ──────────────────────────────────────────
ReferenceEntropyCorrection compute_reference_correction(
    const std::vector<double>& receptor_model_energies,
    int n_rotatable_bonds,
    double temperature_K)
{
    ReferenceEntropyCorrection c;
    c.temperature = temperature_K;

    auto rec = compute_receptor_entropy(receptor_model_energies, temperature_K);
    auto lig = compute_ligand_solution_entropy(n_rotatable_bonds, temperature_K);
    auto cra = compute_cratic_entropy(temperature_K);

    c.T_dS_receptor = rec.T_S_conf;
    c.T_dS_ligand   = lig.T_S_rot;
    c.T_dS_cratic   = cra.T_S_cratic;
    c.T_dS_total    = c.T_dS_receptor + c.T_dS_ligand + c.T_dS_cratic;

    // Claim firewall (see LIB/statmech.h ScientificProvenance).
    //
    // `receptor_model_energies` are uncalibrated model scores — in the only
    // wired caller (ParallelCampaign) they are literally surrogate placeholder
    // values. Feeding them through a kB/T expression produces a MODEL-SCALE
    // diagnostic, not a calibrated kcal/mol reference-state correction, so no
    // physical unit may be printed here. Numbers and line order are unchanged;
    // only the labels and the unit wording are proxy-scoped.
    printf("--- Reference-state entropy proxy correction (T parameter = %.0f K) ---\n",
           temperature_K);
    printf("  claim_validity: proxy_only (uncalibrated model scores in, model scale out)\n");
    printf("  Receptor conf term     -T*dS_rec-like = %8.3f [model scale] (%d models)\n",
           c.T_dS_receptor, static_cast<int>(receptor_model_energies.size()));
    printf("  Ligand solution term   -T*dS_lig-like = %8.3f [model scale] (%d rot bonds)\n",
           c.T_dS_ligand, n_rotatable_bonds);
    printf("  Cratic term            -T*dS_cra-like = %8.3f [model scale]\n",
           c.T_dS_cratic);
    printf("  Total correction       -T*dS_ref-like = %8.3f [model scale]\n",
           c.T_dS_total);

    return c;
}

// ─── Ensemble consensus scoring ──────────────────────────────────────────────
EnsembleConsensusResult compute_ensemble_consensus(
    const std::vector<double>& per_model_dG,
    double temperature_K)
{
    EnsembleConsensusResult r;
    r.n_models = static_cast<int>(per_model_dG.size());
    r.per_model_dG = per_model_dG;

    if (r.n_models == 0) {
        r.dG_consensus = 0.0;
        r.dG_best = r.dG_worst = r.dG_mean = r.dG_stddev = 0.0;
        r.best_model_idx = -1;
        return r;
    }

    if (r.n_models == 1) {
        r.dG_consensus = per_model_dG[0];
        r.dG_best = r.dG_worst = r.dG_mean = per_model_dG[0];
        r.dG_stddev = 0.0;
        r.best_model_idx = 0;
        return r;
    }

    double beta = 1.0 / (kB_kcal * temperature_K);

    // Find best (most negative) and worst
    r.dG_best  =  std::numeric_limits<double>::max();
    r.dG_worst = -std::numeric_limits<double>::max();
    r.best_model_idx = 0;

    for (int i = 0; i < r.n_models; i++) {
        if (per_model_dG[i] < r.dG_best) {
            r.dG_best = per_model_dG[i];
            r.best_model_idx = i;
        }
        if (per_model_dG[i] > r.dG_worst)
            r.dG_worst = per_model_dG[i];
    }

    // Mean and stddev
    double sum = 0.0;
    for (double g : per_model_dG) sum += g;
    r.dG_mean = sum / r.n_models;

    double var = 0.0;
    for (double g : per_model_dG) var += (g - r.dG_mean) * (g - r.dG_mean);
    r.dG_stddev = std::sqrt(var / r.n_models);

    // Boltzmann-weighted consensus:
    // ΔG_consensus = -kT · ln[ (1/N) · Σ exp(-βΔG_i) ]
    // Using log-sum-exp for stability
    std::vector<double> exponents(r.n_models);
    for (int i = 0; i < r.n_models; i++)
        exponents[i] = -beta * per_model_dG[i];

    double max_exp = *std::max_element(exponents.begin(), exponents.end());
    double lse = 0.0;
    for (double e : exponents)
        lse += std::exp(e - max_exp);
    double log_avg = max_exp + std::log(lse / r.n_models);

    r.dG_consensus = -log_avg / beta;

    // Claim firewall: `per_model_dG` is an uncalibrated per-model score, so the
    // log-sum-exp below is a soft-min over model scale — it is not a binding
    // free energy and must never be printed as ΔG / kcal/mol. Arithmetic and
    // line order are unchanged; only the labels and units are proxy-scoped.
    printf("--- Ensemble consensus over per-model proxy scores (T parameter = %.0f K, %d models) ---\n",
           temperature_K, r.n_models);
    printf("  Best model [%d]      score-like = %8.3f [model scale]\n",
           r.best_model_idx + 1, r.dG_best);
    printf("  Worst model         score-like = %8.3f [model scale]\n", r.dG_worst);
    printf("  Arithmetic mean     score-like = %8.3f ± %.3f [model scale]\n",
           r.dG_mean, r.dG_stddev);
    printf("  Boltzmann-weighted  score-like = %8.3f [model scale; proxy_only]\n",
           r.dG_consensus);

    return r;
}

} // namespace reference_entropy
