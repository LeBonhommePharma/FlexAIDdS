// DiFTGAAdapter.h — bridge between DiFT torsional potentials and the GA fitness
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// ─────────────────────────────────────────────────────────────────────────────
// This header is the single, clean integration surface between the DiFT engine
// and the FlexAID∆S genetic algorithm. It is header-only and self-contained so
// it can be unit-tested in isolation and called from the fitness evaluator
// (LIB/spfunction.cpp / LIB/cffunction.cpp) behind an opt-in flag without
// disturbing any existing scoring path.
//
// ONE FFT, TWO PAYOFFS — applied per ligand rotatable bond:
//   energy   = Σ_b  V_tors,b(φ_b)            → adds to the complementarity score
//   minus_TS = Σ_b  −T·S_tors,b              → the torsional ΔS free-energy
//                                              penalty (confinement cost)
// The GA minimizes  CF + score_torsional().total(), so a pose that parks a
// rotatable bond in a deep DiFT well pays the well's energy AND its entropy
// penalty, both derived from the same Fourier spectrum.
// ─────────────────────────────────────────────────────────────────────────────
#pragma once

#include "DiFT.h"

#include <span>
#include <vector>

namespace dift {

// A ligand rotatable bond with its DiFT-parametrized torsional potential.
struct RotatableBondTorsion {
    TorsionalPotential potential;     // DiFT V_tors for this bond
    int                gene_index = -1;  // GA dihedral gene that drives φ for this bond
};

// Decomposed torsional contribution to the binding free energy (kcal/mol).
struct TorsionalScore {
    double energy   = 0.0;   // Σ V_tors,b relative to each well minimum
    double minus_TS = 0.0;   // Σ −T·S_tors,b — torsional confinement penalty
    int    n_bonds  = 0;     // rotatable bonds accounted for

    // Total free-energy contribution the GA should add to the CF score.
    double total() const noexcept { return energy + minus_TS; }
};

// Score a pose's torsional state.
//   bonds                — per-rotatable-bond DiFT potentials
//   dihedral_angles_rad  — current dihedral value (radians) of each bond,
//                          in the SAME order as `bonds`
//   temperature_K        — ensemble temperature for the entropy term
//
// `energy` uses the potential RELATIVE to its minimum (only relative torsional
// energies are physical); `minus_TS` is the temperature-weighted excess
// torsional entropy from DiFTEngine::thermodynamics — a free rotor contributes
// exactly zero, a confined bond contributes a positive penalty.
inline TorsionalScore score_torsional(std::span<const RotatableBondTorsion> bonds,
                                      std::span<const double> dihedral_angles_rad,
                                      double temperature_K = 300.0) {
    TorsionalScore s;
    if (bonds.size() != dihedral_angles_rad.size())
        return s;   // caller contract violated — return a zero contribution

    DiFTEngine engine(temperature_K);
    for (std::size_t b = 0; b < bonds.size(); ++b) {
        const TorsionalPotential& pot = bonds[b].potential;
        s.energy   += pot.relative(dihedral_angles_rad[b]);
        s.minus_TS += engine.thermodynamics(pot).minus_TS;
        ++s.n_bonds;
    }
    return s;
}

// Convenience: build a RotatableBondTorsion by parametrizing a raw torsional
// profile (QM scan or Boltzmann-inverted CG histogram) on the spot.
inline RotatableBondTorsion make_bond_torsion(std::span<const double> profile,
                                              int    gene_index,
                                              double temperature_K   = 300.0,
                                              int    max_multiplicity = 6) {
    DiFTEngine engine(temperature_K);
    RotatableBondTorsion rbt;
    rbt.potential   = engine.parametrize(profile, max_multiplicity);
    rbt.gene_index  = gene_index;
    return rbt;
}

} // namespace dift
