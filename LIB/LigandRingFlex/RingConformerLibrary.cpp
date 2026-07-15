// RingConformerLibrary.cpp — ring conformer tables and detection
#include "RingConformerLibrary.h"
#include "../RngSeed.h"
#include "../flexaid.h"   // atom_struct (bond graph), MBNDS
#include <cstdlib>
#include <cstring>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

namespace ring_flex {

// ─── RingConformerLibrary ────────────────────────────────────────────────────

RingConformerLibrary& RingConformerLibrary::instance() {
    static RingConformerLibrary lib;
    return lib;
}

RingConformerLibrary::RingConformerLibrary() {
    build_six_conformers();
    build_five_conformers();
}

void RingConformerLibrary::build_six_conformers() {
    // Reference dihedral tables from Cremer & Pople (JACS 1975) for pyranose rings.
    // Angles in degrees; ν0 = C1-C2-C3-C4, ν1 = C2-C3-C4-C5, …
    six_ = {
        { SixConformerType::Chair1,   "4C1",
          { 55.9f, -55.9f,  55.9f, -55.9f,  55.9f, -55.9f } },
        { SixConformerType::Chair2,   "1C4",
          {-55.9f,  55.9f, -55.9f,  55.9f, -55.9f,  55.9f } },
        { SixConformerType::BoatA,    "2,5B",
          {  0.0f,  60.0f,   0.0f, -60.0f,   0.0f,  60.0f } },
        { SixConformerType::BoatB,    "B2,5",
          {  0.0f, -60.0f,   0.0f,  60.0f,   0.0f, -60.0f } },
        { SixConformerType::TwistBoatA, "2SO",
          { 30.0f,  30.0f, -60.0f,  30.0f,  30.0f, -60.0f } },
        { SixConformerType::TwistBoatB, "OS2",
          {-30.0f, -30.0f,  60.0f, -30.0f, -30.0f,  60.0f } },
        { SixConformerType::HalfChairA, "3H4",
          { 45.0f, -20.0f, -25.0f,  45.0f, -45.0f,  20.0f } },
        { SixConformerType::HalfChairB, "4H3",
          {-45.0f,  20.0f,  25.0f, -45.0f,  45.0f, -20.0f } },
    };
}

void RingConformerLibrary::build_five_conformers() {
    // Five-membered ring envelope conformers (Altona & Sundaralingam, JACS 1972).
    // Phase angles: P = 0° (C3'-endo), 36° (C4'-exo), 72° (O4'-endo),
    //               108° (C1'-exo), 144° (C2'-endo) for nucleoside furanoses.
    five_ = {
        { FiveConformerType::E0, "C3'-endo",
          {  0.0f, -36.0f,  36.0f, -36.0f,  36.0f }, 0.0f },
        { FiveConformerType::E1, "C4'-exo",
          { -36.0f,  36.0f, -36.0f,  36.0f,   0.0f }, 36.0f },
        { FiveConformerType::E2, "O4'-endo",
          {  36.0f, -36.0f,  36.0f,   0.0f, -36.0f }, 72.0f },
        { FiveConformerType::E3, "C1'-exo",
          { -36.0f,  36.0f,   0.0f, -36.0f,  36.0f }, 108.0f },
        { FiveConformerType::E4, "C2'-endo",
          {  36.0f,   0.0f, -36.0f,  36.0f, -36.0f }, 144.0f },
    };
}

int RingConformerLibrary::random_six_index() const {
    auto& rng = flexaids_rng::lazy_thread_rng(0x516ULL);
    return std::uniform_int_distribution<int>(0, n_six() - 1)(rng);
}

int RingConformerLibrary::random_five_index() const {
    auto& rng = flexaids_rng::lazy_thread_rng(0x515ULL);
    return std::uniform_int_distribution<int>(0, n_five() - 1)(rng);
}

// ─── detect_non_aromatic_rings ───────────────────────────────────────────────
// Finds 5- and 6-membered rings in the ligand by DFS on the bond graph,
// then classifies each as aromatic or non-aromatic using an sp2 heuristic.
// Only non-aromatic rings are returned (aromatic rings are rigid).

// Check if a ring is likely aromatic: all ring atoms are sp2-hybridised
// (≤ 3 bonds, element is C/N/O/S) and the ring size is 5 or 6.
static bool is_likely_aromatic(const std::vector<int>& ring,
                               const ::atom* atoms) {
    for (int idx : ring) {
        int nbonds = atoms[idx].bond[0];
        if (nbonds > 3) return false;
        char elem = atoms[idx].element[0];
        // sp2 atoms in aromatic rings: C, N, O, S with 2-3 bonds
        if (elem != 'C' && elem != 'N' && elem != 'O' && elem != 'S')
            return false;
        if (nbonds < 2) return false;
    }
    return true;
}

std::vector<RingDescriptor> detect_non_aromatic_rings(
    const int* atom_indices, int n_atoms,
    const ::atom* atoms)
{
    std::vector<RingDescriptor> result;
    if (!atom_indices || n_atoms <= 0 || !atoms) return result;

    // Build a set of ligand atom indices for fast membership lookup
    std::unordered_set<int> lig_set(atom_indices, atom_indices + n_atoms);

    // Build local adjacency restricted to ligand atoms
    std::unordered_map<int, std::vector<int>> adj;
    for (int i = 0; i < n_atoms; ++i) {
        int ai = atom_indices[i];
        int nb = atoms[ai].bond[0];
        for (int b = 1; b <= nb && b <= 6; ++b) {
            int neighbor = atoms[ai].bond[b];
            if (lig_set.count(neighbor))
                adj[ai].push_back(neighbor);
        }
    }

    // Find small rings (size 5 and 6) via bounded DFS from each atom.
    // For each edge (u,v), do a BFS/DFS from v (without going through u)
    // looking for a path back to u of length 4 or 5 (giving rings of 5 or 6).
    std::vector<std::vector<int>> found_rings;
    std::unordered_set<int> visited_global;

    for (int i = 0; i < n_atoms; ++i) {
        int start = atom_indices[i];
        for (int nb : adj[start]) {
            // Search for paths from nb back to start, max depth 5
            // (yielding rings of size 3..6; we keep only 5 and 6)
            // Use iterative DFS with path tracking
            struct Frame { int node; int parent; std::vector<int> path; };
            std::vector<Frame> stack;
            stack.push_back({nb, start, {start, nb}});

            while (!stack.empty()) {
                Frame f = std::move(stack.back());
                stack.pop_back();

                if ((int)f.path.size() > 6) continue; // max ring size 6

                for (int next : adj[f.node]) {
                    if (next == f.parent && (int)f.path.size() < 4) continue;
                    if (next == start && (int)f.path.size() >= 5) {
                        // Found a ring of size path.size()
                        std::vector<int> ring = f.path;
                        // Canonicalize: rotate so smallest index is first
                        auto it = std::min_element(ring.begin(), ring.end());
                        std::rotate(ring.begin(), it, ring.end());
                        // Normalize direction
                        if (ring.size() > 2 && ring[1] > ring.back())
                            std::reverse(ring.begin() + 1, ring.end());

                        // Check for duplicate
                        bool dup = false;
                        for (const auto& existing : found_rings) {
                            if (existing == ring) { dup = true; break; }
                        }
                        if (!dup) found_rings.push_back(ring);
                        continue;
                    }
                    // Don't revisit atoms already in path
                    bool in_path = false;
                    for (int p : f.path) {
                        if (p == next) { in_path = true; break; }
                    }
                    if (in_path) continue;
                    if ((int)f.path.size() < 6) {
                        std::vector<int> new_path = f.path;
                        new_path.push_back(next);
                        stack.push_back({next, f.node, std::move(new_path)});
                    }
                }
            }
        }
    }

    // Convert to RingDescriptors, filtering out aromatic rings
    for (auto& ring : found_rings) {
        int sz = static_cast<int>(ring.size());
        if (sz != 5 && sz != 6) continue;

        bool aromatic = is_likely_aromatic(ring, atoms);
        if (aromatic) continue; // skip aromatic rings (rigid)

        RingDescriptor desc;
        desc.atom_indices = std::move(ring);
        desc.size = (sz == 5) ? RingSize::Five : RingSize::Six;
        desc.is_aromatic = false;
        result.push_back(std::move(desc));
    }

    return result;
}

} // namespace ring_flex
