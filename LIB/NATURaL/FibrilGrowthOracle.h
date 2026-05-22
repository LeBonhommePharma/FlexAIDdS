// FibrilGrowthOracle.h — Grand canonical fibril elongation gate
//
// For protofibril_n + monomer ⇌ protofibril_{n+1}, the grand partition function
//   Ξ = 1 + z · Z,  z = c_monomer / c°
// gives p(elongation) = z·Z / Ξ and ΔG_elong = −kT ln Z + kT ln(c°/c_monomer).
//
// Reuses target::GrandPartitionFunction (see LIB/GrandPartitionFunction.h:46), which
// already provides log-space arithmetic and concentration-aware fugacities.
//
// References: docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md §3
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#pragma once

#include "../GrandPartitionFunction.h"
#include "../statmech.h"

#include <memory>

namespace natural {

struct ElongationDecision {
    double p_elong;   // ∈ [0, 1)
    double dG_elong;  // kcal/mol — apparent ΔG at configured monomer concentration
    bool   gated_in;  // p_elong ≥ acceptance_threshold
};

class FibrilGrowthOracle {
public:
    explicit FibrilGrowthOracle(double temperature_K       = 310.15,
                                double acceptance_threshold = 0.5);

    // Returns the elongation decision for the current monomer ensemble.
    //
    // monomer_engine: StatMechEngine populated by a Sim C GA run, with the protofibril
    //                 as receptor and a single free monomer as ligand. Internally we
    //                 derive ln Z from the engine's free energy: ln Z = -F / (kT).
    // c_monomer_M:    free monomer concentration in Molar. Default = 1 µM, the upper
    //                 bound for physiological cytosolic free Aβ42 (Bjorkdahl 2008).
    ElongationDecision gate(const statmech::StatMechEngine& monomer_engine,
                            double                          c_monomer_M = 1.0e-6);

    double temperature_K()        const noexcept { return T_K_; }
    double acceptance_threshold() const noexcept { return acceptance_threshold_; }

private:
    double T_K_;
    double acceptance_threshold_;
    // Single-ligand GrandPartitionFunction. We re-register the monomer entry on each
    // gate() call so a fresh Sim C ensemble overwrites the previous Z.
    std::unique_ptr<target::GrandPartitionFunction> Xi_;
};

} // namespace natural
