// LigandRingFlex.h — unified ring flexibility interface for the GA
//
// Combines non-aromatic ring conformers (RingConformerLibrary) and
// furanose sugar pucker (SugarPucker) under a single API that the GA
// calls during initialisation, mutation, and pose evaluation.
#pragma once

#include "RingConformerLibrary.h"
#include "SugarPucker.h"
#include <vector>
#include <cstdint>

// flexaid.h is already included via RingConformerLibrary.h

namespace ligand_ring_flex {

// ─── chromosome extensions ───────────────────────────────────────────────────
// These are appended to the GA chromosome conceptually; in the FlexAID
// codebase they live as supplementary side-channel arrays alongside the
// existing `gene` array (not binary-encoded — stored as raw indices/floats).

struct RingFlexGenes {
    // One entry per detected non-aromatic ring in the ligand.
    // Value: index into RingConformerLibrary (six_conformers or five_conformers).
    std::vector<uint8_t> conformer_indices; // 6-membered rings
    std::vector<uint8_t> five_conformer_indices; // 5-membered non-sugar rings

    // Sugar pucker phases (degrees) — one per detected furanose ring
    std::vector<float>                       sugar_phases;
    std::vector<sugar_pucker::SugarType>     sugar_types;
    std::vector<std::vector<int>>            sugar_ring_indices; // atom idx lists

    // True when at least one flexible ring (6-ring, non-sugar 5-ring, or
    // furanose sugar) was detected — i.e. ring-flex genes are present.
    bool has_rings() const {
        return !conformer_indices.empty() ||
               !five_conformer_indices.empty() ||
               !sugar_phases.empty();
    }
};

// ─── ring detection at startup ────────────────────────────────────────────────
// Analyses the ligand topology stored in atoms[] and populates a RingFlexGenes
// template (0-initialised) ready for population initialisation.
RingFlexGenes detect_rings(const atom* atoms, int n_lig_atoms,
                            const int*  lig_atom_offset);

// ─── population initialisation ───────────────────────────────────────────────
// Randomise all ring conformer indices for a new individual.
void randomise(RingFlexGenes& genes);

// ─── mutation ────────────────────────────────────────────────────────────────
// With probability `ring_mut_prob` flip a random ring conformer.
// With probability `pucker_mut_prob` mutate a random sugar pucker phase.
void mutate(RingFlexGenes& genes,
            double ring_mut_prob   = 0.05,
            double pucker_mut_prob = 0.12);

// ─── crossover ───────────────────────────────────────────────────────────────
// Single-point crossover of ring indices and sugar phases between two individuals.
void crossover(RingFlexGenes& parent_a, RingFlexGenes& parent_b);

// ─── pose application ────────────────────────────────────────────────────────
// Apply the ring-flex genes to the atom array (updates internal coordinates).
// Called inside calculate_fitness() before the CF evaluation.
void apply(atom* atoms, const RingFlexGenes& genes);

// ─── Shannon entropy contribution ────────────────────────────────────────────
// Compute bits of configurational entropy from ring conformer diversity
// in a population of RingFlexGenes (one entry per GA individual).
double compute_ring_entropy(const std::vector<RingFlexGenes>& population);

} // namespace ligand_ring_flex
