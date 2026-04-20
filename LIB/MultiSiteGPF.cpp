// MultiSiteGPF.cpp — Multi-site grand canonical partition function
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "MultiSiteGPF.h"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace target {

MultiSiteGPF::MultiSiteGPF(double temperature_K)
    : T_(temperature_K)
{
}

// ── Site management ────────────────────────────────────────────────────

int MultiSiteGPF::add_site(const std::string& site_name)
{
    if (site_name_to_idx_.count(site_name))
        throw std::invalid_argument("Site '" + site_name + "' already registered");

    int idx = static_cast<int>(sites_.size());
    sites_.push_back(std::make_unique<GrandPartitionFunction>(T_));
    site_name_to_idx_[site_name] = idx;

    // Extend coupling matrix (N×N → (N+1)×(N+1))
    for (auto& row : coupling_) {
        row.push_back(1.0);  // default: independent
    }
    coupling_.push_back(std::vector<double>(idx + 1, 1.0));
    coupling_[idx][idx] = 1.0;  // self-coupling = 1

    return idx;
}

int MultiSiteGPF::num_sites() const
{
    return static_cast<int>(sites_.size());
}

int MultiSiteGPF::site_index(const std::string& name) const
{
    auto it = site_name_to_idx_.find(name);
    return (it != site_name_to_idx_.end()) ? it->second : -1;
}

GrandPartitionFunction& MultiSiteGPF::site(int idx)
{
    validate_site_idx(idx);
    return *sites_[idx];
}

const GrandPartitionFunction& MultiSiteGPF::site(int idx) const
{
    validate_site_idx(idx);
    return *sites_[idx];
}

// ── Ligand registration ────────────────────────────────────────────────

void MultiSiteGPF::add_ligand(int site_idx, const std::string& ligand_name,
                               double log_Z, double concentration_M)
{
    validate_site_idx(site_idx);
    sites_[site_idx]->add_ligand(ligand_name, log_Z, concentration_M);
}

void MultiSiteGPF::add_ligand(const std::string& site_name,
                               const std::string& ligand_name,
                               double log_Z, double concentration_M)
{
    int idx = site_index(site_name);
    if (idx < 0)
        throw std::invalid_argument("Site '" + site_name + "' not found");
    add_ligand(idx, ligand_name, log_Z, concentration_M);
}

// ── Cooperativity ──────────────────────────────────────────────────────

void MultiSiteGPF::set_cooperativity(int site_a, int site_b, double omega)
{
    if (omega <= 0.0)
        throw std::invalid_argument("Cooperativity factor must be positive");
    validate_site_idx(site_a);
    validate_site_idx(site_b);

    std::lock_guard<std::mutex> lock(coupling_mtx_);
    coupling_[site_a][site_b] = omega;
    coupling_[site_b][site_a] = omega;  // symmetric
}

double MultiSiteGPF::cooperativity(int site_a, int site_b) const
{
    validate_site_idx(site_a);
    validate_site_idx(site_b);
    std::lock_guard<std::mutex> lock(coupling_mtx_);
    return coupling_[site_a][site_b];
}

// ── Thermodynamic queries ──────────────────────────────────────────────

double MultiSiteGPF::log_Xi() const
{
    if (sites_.empty()) return 0.0;

    // For independent sites: Ξ_total = Π_s Ξ_s  →  ln Ξ = Σ_s ln Ξ_s
    double log_xi_sum = 0.0;
    for (const auto& s : sites_) {
        log_xi_sum += s->log_Xi();
    }

    // Check if any coupling deviates from 1.0
    bool has_coupling = false;
    {
        std::lock_guard<std::mutex> lock(coupling_mtx_);
        int n = static_cast<int>(coupling_.size());
        for (int a = 0; a < n && !has_coupling; ++a) {
            for (int b = a + 1; b < n; ++b) {
                if (std::abs(coupling_[a][b] - 1.0) > 1e-12) {
                    has_coupling = true;
                    break;
                }
            }
        }
    }

    if (!has_coupling) return log_xi_sum;

    // With cooperativity:
    // Ξ_total = Π_s Ξ_s + Σ_{pairs} (ω_ab - 1) · (Ξ_a - 1)(Ξ_b - 1) + ...
    // First-order pairwise approximation:
    //   ΔΞ = Σ_{a<b} (ω_ab - 1) · Π_{s≠a,b} Ξ_s · (Ξ_a - 1)(Ξ_b - 1)
    //
    // In log space: compute each pair's contribution and add via log_sum_exp

    int n = static_cast<int>(sites_.size());
    std::vector<double> log_xi_per_site(n);
    for (int s = 0; s < n; ++s) {
        log_xi_per_site[s] = sites_[s]->log_Xi();
    }

    // Collect all terms in log space for log_sum_exp
    // Term 0: the independent product (base term)
    std::vector<double> log_terms = {log_xi_sum};  // Π_s Ξ_s

    {
        std::lock_guard<std::mutex> lock(coupling_mtx_);
        for (int a = 0; a < n; ++a) {
            for (int b = a + 1; b < n; ++b) {
                double omega = coupling_[a][b];
                if (std::abs(omega - 1.0) < 1e-12) continue;

                // Coupling contribution:
                // (ω_ab - 1) · Π_{s≠a,b} Ξ_s · (Ξ_a - 1)(Ξ_b - 1)
                // In log: ln|ω - 1| + Σ_{s≠a,b} ln Ξ_s + ln|Ξ_a - 1| + ln|Ξ_b - 1|

                double log_product_other = 0.0;
                for (int s = 0; s < n; ++s) {
                    if (s == a || s == b) continue;
                    log_product_other += log_xi_per_site[s];
                }

                // Ξ_a - 1 = Σ_i z_i·Z_i (all ligands at site a, no empty term)
                double xi_a = std::exp(log_xi_per_site[a]);
                double xi_b = std::exp(log_xi_per_site[b]);

                if (xi_a <= 1.0 || xi_b <= 1.0) continue;  // no ligands at one site

                double log_xi_a_minus_1 = std::log(xi_a - 1.0);
                double log_xi_b_minus_1 = std::log(xi_b - 1.0);

                double log_coupling = std::log(std::abs(omega - 1.0));
                double log_term = log_coupling + log_product_other
                                + log_xi_a_minus_1 + log_xi_b_minus_1;

                // Only add positive contributions to first-order approx
                if (omega > 1.0) {
                    log_terms.push_back(log_term);
                }
            }
        }
    }

    // log_sum_exp over all terms
    double max_val = *std::max_element(log_terms.begin(), log_terms.end());
    double sum = 0.0;
    for (double lt : log_terms) {
        sum += std::exp(lt - max_val);
    }
    return max_val + std::log(sum);
}

double MultiSiteGPF::empty_probability() const
{
    return std::exp(-log_Xi());
}

double MultiSiteGPF::binding_probability(int site_idx,
                                          const std::string& ligand_name) const
{
    validate_site_idx(site_idx);
    // Marginal probability: p(ligand at site s)
    // For independent sites: just the single-site probability
    // For coupled sites: first-order approximation using site occupancy
    double p_site = sites_[site_idx]->binding_probability(ligand_name);

    // Check for coupling adjustments
    int n = num_sites();
    bool has_coupling = false;
    {
        std::lock_guard<std::mutex> lock(coupling_mtx_);
        for (int b = 0; b < n; ++b) {
            if (b == site_idx) continue;
            if (std::abs(coupling_[site_idx][b] - 1.0) > 1e-12) {
                has_coupling = true;
                break;
            }
        }
    }

    if (!has_coupling) return p_site;

    // Coupled approximation:
    // p(i,s) ≈ p_independent(i,s) · [1 + Σ_{s'≠s} (ω_{ss'} - 1) · occupancy(s')]
    double occupancy_factor = 1.0;
    {
        std::lock_guard<std::mutex> lock(coupling_mtx_);
        for (int b = 0; b < n; ++b) {
            if (b == site_idx) continue;
            double omega = coupling_[site_idx][b];
            if (std::abs(omega - 1.0) > 1e-12) {
                double occ_b = sites_[b]->empty_probability();
                double p_bound_b = 1.0 - occ_b;  // total occupancy at site b
                occupancy_factor += (omega - 1.0) * p_bound_b;
            }
        }
    }

    return std::min(1.0, p_site * occupancy_factor);
}

double MultiSiteGPF::site_occupancy(int site_idx) const
{
    validate_site_idx(site_idx);
    return 1.0 - sites_[site_idx]->empty_probability();
}

double MultiSiteGPF::selectivity(int site_idx,
                                  const std::string& a,
                                  const std::string& b) const
{
    validate_site_idx(site_idx);
    return sites_[site_idx]->selectivity(a, b);
}

// ── Cross-site analysis ────────────────────────────────────────────────

std::vector<MultiSiteGPF::CrossSiteSelectivity>
MultiSiteGPF::cross_site_analysis() const
{
    // Collect all ligand names across all sites
    std::unordered_map<std::string, std::vector<int>> ligand_sites;
    for (int s = 0; s < num_sites(); ++s) {
        for (const auto& [name, log_Z] : sites_[s]->all_log_Z()) {
            ligand_sites[name].push_back(s);
        }
    }

    std::vector<CrossSiteSelectivity> results;
    for (const auto& [ligand_name, site_list] : ligand_sites) {
        CrossSiteSelectivity css;
        css.ligand_name = ligand_name;
        css.site_probabilities.resize(num_sites(), 0.0);

        double best_prob = -1.0;
        for (int s : site_list) {
            double p = binding_probability(s, ligand_name);
            css.site_probabilities[s] = p;
            if (p > best_prob) {
                best_prob = p;
                css.best_site_idx = s;
            }
        }

        // Compute selectivity vs next-best site
        css.best_site_selectivity = 0.0;
        for (int s : site_list) {
            if (s == css.best_site_idx) continue;
            [[maybe_unused]] double sel = selectivity(css.best_site_idx, ligand_name,
                                     ligand_name);            // selectivity of same ligand at different sites = Z_best / Z_other
            // Cross-site selectivity: ratio of binding probabilities
            if (css.site_probabilities[s] > 0.0) {
                double ratio = best_prob / css.site_probabilities[s];
                if (ratio > css.best_site_selectivity)
                    css.best_site_selectivity = ratio;
            }
        }

        results.push_back(std::move(css));
    }

    return results;
}

// ── Validation ─────────────────────────────────────────────────────────

void MultiSiteGPF::validate_site_idx(int idx) const
{
    if (idx < 0 || idx >= static_cast<int>(sites_.size()))
        throw std::invalid_argument("Site index " + std::to_string(idx)
                                    + " out of range [0, "
                                    + std::to_string(sites_.size()) + ")");
}

} // namespace target
