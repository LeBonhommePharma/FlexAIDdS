// FibrilGrowthOracle.cpp — see header for design.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include "FibrilGrowthOracle.h"

#include <cmath>
#include <stdexcept>

namespace natural {

FibrilGrowthOracle::FibrilGrowthOracle(double temperature_K,
                                       double acceptance_threshold)
    : T_K_(temperature_K),
      acceptance_threshold_(acceptance_threshold),
      Xi_(std::make_unique<target::GrandPartitionFunction>(temperature_K))
{
    if (!(T_K_ > 0.0))
        throw std::invalid_argument("FibrilGrowthOracle: temperature must be > 0");
    if (acceptance_threshold_ < 0.0 || acceptance_threshold_ > 1.0)
        throw std::invalid_argument("FibrilGrowthOracle: acceptance_threshold must be in [0,1]");
}

ElongationDecision FibrilGrowthOracle::gate(
    const statmech::StatMechEngine& monomer_engine,
    double                          c_monomer_M)
{
    if (!(c_monomer_M > 0.0))
        throw std::invalid_argument("FibrilGrowthOracle: c_monomer_M must be > 0");

    // ln Z = -F / (kT), with kT in kcal/mol.
    const auto th = monomer_engine.compute();
    const double kT_kcal = statmech::kB_kcal * T_K_;
    const double log_Z   = -th.free_energy / kT_kcal;

    // Replace any prior "monomer" entry with the fresh ensemble at the current
    // concentration. add_or_overwrite is thread-safe and avoids the TOCTOU race.
    Xi_->add_or_overwrite("monomer", log_Z, c_monomer_M);

    const double p_elong  = Xi_->binding_probability("monomer");
    // Apparent concentration-corrected elongation free energy:
    //   ΔG_app = -kT ln Z + kT ln(c°/c)
    // GrandPartitionFunction::delta_G_bind() is concentration-independent by
    // design, so apply the cratic term explicitly here.
    const double F_bound = Xi_->delta_G_bind("monomer", /*F_ref=*/0.0);
    const double dG_elong = F_bound + kT_kcal * std::log(target::c_standard / c_monomer_M);

    return ElongationDecision{
        .p_elong  = p_elong,
        .dG_elong = dG_elong,
        .gated_in = (p_elong >= acceptance_threshold_),
    };
}

} // namespace natural
