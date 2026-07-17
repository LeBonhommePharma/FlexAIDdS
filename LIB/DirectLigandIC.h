#pragma once

#include "flexaid.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <utility>
#include <vector>

namespace direct_ligand_ic {

using BondGraph = std::vector<std::vector<std::pair<int, int>>>;

struct ReconstructionTree {
    std::array<int, 3> gpa{{0, 0, 0}};  // local atom indices
    std::vector<int> parent;             // local parent, -1 at root
    std::vector<bool> visited;
};

inline bool is_bridge(const BondGraph& graph, int u, int v) {
    std::vector<bool> seen(graph.size(), false);
    std::queue<int> pending;
    pending.push(u);
    seen[u] = true;
    while (!pending.empty()) {
        const int current = pending.front();
        pending.pop();
        for (const auto& [next, order] : graph[current]) {
            (void)order;
            if ((current == u && next == v) || (current == v && next == u)) continue;
            if (!seen[next]) {
                seen[next] = true;
                pending.push(next);
            }
        }
    }
    return !seen[v];
}

inline int bond_order(const BondGraph& graph, int u, int v) {
    for (const auto& [next, order] : graph[u]) {
        if (next == v) return order;
    }
    return 0;
}

inline int heavy_degree(const BondGraph& graph,
                        const std::vector<bool>& is_heavy,
                        int atom_index) {
    int degree = 0;
    for (const auto& [next, order] : graph[atom_index]) {
        (void)order;
        if (is_heavy[next]) ++degree;
    }
    return degree;
}

// Element token from atom.element / atom name (H, C, N, O, …).
inline char element_letter(const atom& a) {
    if (a.element[0] != '\0')
        return static_cast<char>(std::toupper(static_cast<unsigned char>(a.element[0])));
    // Fallback: first non-space of PDB-style name.
    for (int i = 0; i < 5 && a.name[i]; ++i) {
        if (a.name[i] != ' ')
            return static_cast<char>(std::toupper(static_cast<unsigned char>(a.name[i])));
    }
    return '?';
}

// Resonance-locked C–N (amide / urea / carbamate / guanidine / amidine).
//
// MUST NOT use VCT type==11 as amide proof: SDF generic N and MOL2 N.1/N.2/N.3
// all map to type 11 for matrix scoring, so type 11 freezes ordinary amine C–N.
// Only:
//   1) explicit MOL2 amide bond order 10 (am → 10 in Mol2Reader), or
//   2) structural conjugation: C of the C–N has a double bond (order ≥ 2) to O or N.
inline bool is_resonance_locked_cn(const BondGraph& graph,
                                   const atom* atoms,
                                   int first_atom,
                                   int u,
                                   int v) {
    // Explicit MOL2 amide bond order encoding (am → 10).
    if (bond_order(graph, u, v) == 10) return true;

    const atom& au = atoms[first_atom + u];
    const atom& av = atoms[first_atom + v];
    const char eu = element_letter(au);
    const char ev = element_letter(av);
    int c_local = -1;
    if (eu == 'C' && ev == 'N') { c_local = u; }
    else if (eu == 'N' && ev == 'C') { c_local = v; }
    else return false;

    // Structural C(=O/N)–N conjugation only (not atom.type).
    for (const auto& [nb, order] : graph[c_local]) {
        if (order < 2) continue;
        const char en = element_letter(atoms[first_atom + nb]);
        if (en == 'O' || en == 'N') return true;
    }
    return false;
}

inline bool is_rotatable(const BondGraph& graph,
                         const std::vector<bool>& is_heavy,
                         int u,
                         int v) {
    const int order = bond_order(graph, u, v);
    // order 1 = ordinary single; order 10 = MOL2 amide (never a rotor).
    if (!(order == 1)) return false;
    if (heavy_degree(graph, is_heavy, u) < 2 ||
        heavy_degree(graph, is_heavy, v) < 2) return false;
    if (!is_bridge(graph, u, v)) return false;
    return true;
}

// Atom-aware overload used by configure_rotatable_bonds / choose_gpa scoring.
inline bool is_rotatable(const BondGraph& graph,
                         const std::vector<bool>& is_heavy,
                         const atom* atoms,
                         int first_atom,
                         int u,
                         int v) {
    if (!is_rotatable(graph, is_heavy, u, v)) return false;
    if (atoms && is_resonance_locked_cn(graph, atoms, first_atom, u, v))
        return false;
    return true;
}

inline double frame_quality(const atom* atoms, int first_atom,
                            int a, int b, int c) {
    const atom& aa = atoms[first_atom + a];
    const atom& ab = atoms[first_atom + b];
    const atom& ac = atoms[first_atom + c];
    const double ux = aa.coor[0] - ab.coor[0];
    const double uy = aa.coor[1] - ab.coor[1];
    const double uz = aa.coor[2] - ab.coor[2];
    const double vx = ac.coor[0] - ab.coor[0];
    const double vy = ac.coor[1] - ab.coor[1];
    const double vz = ac.coor[2] - ab.coor[2];
    const double un = std::sqrt(ux * ux + uy * uy + uz * uz);
    const double vn = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (un <= 1e-8 || vn <= 1e-8) return 0.0;
    const double cx = uy * vz - uz * vy;
    const double cy = uz * vx - ux * vz;
    const double cz = ux * vy - uy * vx;
    return std::sqrt(cx * cx + cy * cy + cz * cz) / (un * vn);
}

// Choose GPA frame for IC reconstruction.
//
// Preferred: ≥3 non-collinear **heavy** atoms (full flexible docking).
// Safe rigid / virtual-frame path (parse + rigid docking preserved):
//   - water / methane / 1–2 heavy: use all atoms including H for a non-collinear
//     triple when possible (H₂O: O–H–H; CH₄: C–H–H);
//   - single atom or collinear-only: rigid virtual frame with FA->ori grandparents
//     (duplicate GPA slots allowed; build_tree already handles g1==g0 / g2==g*).
// Empty molecule only is a hard fail. Never rejects parse of valid small ligands.
inline bool choose_gpa(const atom* atoms,
                       int first_atom,
                       const BondGraph& graph,
                       const std::vector<bool>& is_heavy,
                       std::array<int, 3>& best) {
    const int n = static_cast<int>(graph.size());
    if (n <= 0) {
        std::fprintf(stderr, "ERROR [DIRECT-IC]: empty ligand — no GPA\n");
        best = {{-1, -1, -1}};
        return false;
    }

    std::vector<int> heavy;
    std::vector<int> all;
    heavy.reserve(static_cast<size_t>(n));
    all.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        all.push_back(i);
        if (is_heavy[i]) heavy.push_back(i);
    }

    auto try_triple = [&](int a, int b, int c) -> bool {
        if (a < 0 || b < 0 || c < 0) return false;
        if (a == b || b == c || a == c) return false;
        return frame_quality(atoms, first_atom, a, b, c) >= 1e-3;
    };

    // ── Phase 1: preferred bonded heavy GPA (docking-quality frame) ──────────
    if (heavy.size() >= 3) {
        best = {{heavy[0], heavy[1], heavy[2]}};
        double best_score = -1.0;

        for (int center = 0; center < n; ++center) {
            if (!is_heavy[center]) continue;
            for (const auto& [left_raw, left_order] : graph[center]) {
                (void)left_order;
                if (!is_heavy[left_raw]) continue;
                for (const auto& [right_raw, right_order] : graph[center]) {
                    (void)right_order;
                    if (left_raw >= right_raw || !is_heavy[right_raw]) continue;
                    const int left = left_raw;
                    const int right = right_raw;
                    const double quality = frame_quality(atoms, first_atom,
                                                         left, center, right);
                    if (quality < 1e-3) continue;
                    const int rigid_edges =
                        (!is_rotatable(graph, is_heavy, atoms, first_atom, left, center) ? 1 : 0) +
                        (!is_rotatable(graph, is_heavy, atoms, first_atom, center, right) ? 1 : 0);
                    const int local_degree = heavy_degree(graph, is_heavy, center) +
                                             heavy_degree(graph, is_heavy, left) +
                                             heavy_degree(graph, is_heavy, right);
                    const double score = 10000.0 * rigid_edges +
                                         100.0 * quality + local_degree;
                    if (score > best_score) {
                        best_score = score;
                        best = {{left, center, right}};
                    }
                }
            }
        }
        if (best_score >= 0.0 && try_triple(best[0], best[1], best[2])) {
            return true;
        }

        // Unbonded / collinear bonded frames: any non-collinear heavy triple.
        for (size_t i = 0; i < heavy.size(); ++i) {
            for (size_t j = i + 1; j < heavy.size(); ++j) {
                for (size_t k = j + 1; k < heavy.size(); ++k) {
                    if (try_triple(heavy[i], heavy[j], heavy[k])) {
                        best = {{heavy[i], heavy[j], heavy[k]}};
                        return true;
                    }
                }
            }
        }
    }

    // ── Phase 2: rigid frame including H (water, methane, 1–2 heavy) ─────────
    // Reader support is preserved; ligands with <3 non-collinear heavies dock as
    // rigid bodies (fdih typically 0 — no bridge heavy–heavy rotors).
    for (size_t i = 0; i < all.size(); ++i) {
        for (size_t j = i + 1; j < all.size(); ++j) {
            for (size_t k = j + 1; k < all.size(); ++k) {
                if (try_triple(all[i], all[j], all[k])) {
                    best = {{all[i], all[j], all[k]}};
                    std::fprintf(stderr,
                        "WARN [DIRECT-IC]: rigid/virtual-frame GPA (heavy=%zu total=%d) "
                        "local=%d,%d,%d — parse OK; dock as rigid if fdih=0\n",
                        heavy.size(), n, best[0] + 1, best[1] + 1, best[2] + 1);
                    return true;
                }
            }
        }
    }

    // ── Phase 3: virtual FA->ori frame for 1–2 atoms / fully collinear ────────
    // build_tree uses rec=0 → FA->ori for missing grandparents. Allow duplicate
    // GPA slots so single-atom (Xe) and diatomic (CO) parse and type cleanly.
    if (n == 1) {
        best = {{0, 0, 0}};
        std::fprintf(stderr,
            "WARN [DIRECT-IC]: single-atom rigid virtual-frame GPA (FA->ori)\n");
        return true;
    }
    if (n == 2) {
        best = {{0, 1, 0}};
        std::fprintf(stderr,
            "WARN [DIRECT-IC]: diatomic rigid virtual-frame GPA (FA->ori)\n");
        return true;
    }
    // ≥3 collinear atoms: still assign first three (CoordBuilder/buildcc
    // tolerate degenerate frames via perpendicular fallback).
    best = {{0, 1, std::min(2, n - 1)}};
    if (best[2] == best[1]) best[2] = 0;
    std::fprintf(stderr,
        "WARN [DIRECT-IC]: collinear/degenerate rigid virtual-frame GPA "
        "local=%d,%d,%d\n",
        best[0] + 1, best[1] + 1, best[2] + 1);
    return true;
}

// Returns false if GPA selection fails (ligand load must abort).
inline bool build_tree(atom* atoms,
                       resid& ligand,
                       int first_atom,
                       const BondGraph& graph,
                       const std::vector<bool>& is_heavy,
                       ReconstructionTree& tree) {
    const int n = static_cast<int>(graph.size());
    tree.parent.assign(n, -1);
    tree.visited.assign(n, false);
    if (!choose_gpa(atoms, first_atom, graph, is_heavy, tree.gpa)) {
        return false;
    }

    if (!ligand.gpa) {
        ligand.gpa = static_cast<int*>(std::malloc(3 * sizeof(int)));
        if (!ligand.gpa) {
            std::fprintf(stderr, "ERROR: direct-ligand GPA allocation failed\n");
            std::abort();
        }
    }
    for (int i = 0; i < 3; ++i) ligand.gpa[i] = first_atom + tree.gpa[i];

    const int g0 = tree.gpa[0];
    const int g1 = tree.gpa[1];
    const int g2 = tree.gpa[2];
    std::queue<int> pending;
    tree.visited[g0] = true;
    pending.push(g0);
    if (g1 != g0) {
        tree.visited[g1] = true;
        tree.parent[g1] = g0;
        pending.push(g1);
    }
    if (g2 != g0 && g2 != g1) {
        tree.visited[g2] = true;
        tree.parent[g2] = g1;
        pending.push(g2);
    }

    while (!pending.empty()) {
        const int current = pending.front();
        pending.pop();
        for (const auto& [next, order] : graph[current]) {
            (void)order;
            if (!tree.visited[next]) {
                tree.visited[next] = true;
                tree.parent[next] = current;
                pending.push(next);
            }
        }
    }

    auto absolute = [first_atom](int local) {
        return local >= 0 ? first_atom + local : 0;
    };
    auto moving_reference = [&](int current, int r0, int r1) {
        for (int candidate : tree.gpa) {
            if (candidate != current && candidate != r0 && candidate != r1)
                return candidate;
        }
        return -1;
    };

    for (int local = 0; local < n; ++local) {
        atom& current = atoms[first_atom + local];
        current.recs = 'm';
        if (local == g0) {
            current.rec[0] = current.rec[1] = current.rec[2] = 0;
            continue;
        }
        if (local == g1) {
            current.rec[0] = absolute(g0);
            current.rec[1] = current.rec[2] = 0;
            continue;
        }
        if (local == g2) {
            current.rec[0] = absolute(g1);
            current.rec[1] = absolute(g0);
            current.rec[2] = 0;
            continue;
        }

        int r0 = tree.parent[local];
        int r1 = r0 >= 0 ? tree.parent[r0] : -1;
        int r2 = r1 >= 0 ? tree.parent[r1] : -1;
        if (r0 < 0) r0 = g0;
        if (r1 < 0 || r1 == local || r1 == r0)
            r1 = moving_reference(local, r0, -1);
        if (r2 < 0 || r2 == local || r2 == r0 || r2 == r1)
            r2 = moving_reference(local, r0, r1);
        current.rec[0] = absolute(r0);
        current.rec[1] = absolute(r1);
        current.rec[2] = absolute(r2);
    }

    std::fprintf(stderr,
                 "[DIRECT-IC] GPA topology atoms=%d,%d,%d local=%d,%d,%d\n",
                 atoms[ligand.gpa[0]].number,
                 atoms[ligand.gpa[1]].number,
                 atoms[ligand.gpa[2]].number,
                 g0 + 1, g1 + 1, g2 + 1);
    return true;
}

inline int configure_rotatable_bonds(atom* atoms,
                                     resid& ligand,
                                     int first_atom,
                                     const BondGraph& graph,
                                     const std::vector<bool>& is_heavy,
                                     const ReconstructionTree& tree) {
    const int n = static_cast<int>(graph.size());
    for (int local = 0; local < n; ++local) {
        atoms[first_atom + local].rec[3] = 0;
        atoms[first_atom + local].shift = 0.0f;
    }

    int fdih = 0;
    for (int child = 0; child < n; ++child) {
        const int parent = tree.parent[child];
        if (parent < 0) continue;
        if (child == tree.gpa[1] || child == tree.gpa[2]) continue;
        if (!is_rotatable(graph, is_heavy, atoms, first_atom, parent, child)) continue;

        std::vector<int> controls;
        for (int candidate = 0; candidate < n; ++candidate) {
            if (tree.parent[candidate] == child) controls.push_back(candidate);
        }
        if (controls.empty()) {
            std::fprintf(stderr,
                         "WARNING: no IC control atom for rotatable bond %d-%d\n",
                         atoms[first_atom + parent].number,
                         atoms[first_atom + child].number);
            continue;
        }

        std::sort(controls.begin(), controls.end());
        ligand.bond[++fdih] = first_atom + controls.front();
        if (controls.size() > 1) {
            for (std::size_t i = 0; i < controls.size(); ++i) {
                const int current = first_atom + controls[i];
                const int next = first_atom + controls[(i + 1) % controls.size()];
                atoms[current].rec[3] = next;
                atoms[next].shift = atoms[next].dih - atoms[current].dih;
            }
        }
    }
    ligand.fdih = fdih;
    return fdih;
}

}  // namespace direct_ligand_ic
