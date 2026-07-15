#pragma once

#include "flexaid.h"

#include <algorithm>
#include <array>
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

inline bool is_rotatable(const BondGraph& graph,
                         const std::vector<bool>& is_heavy,
                         int u,
                         int v) {
    return bond_order(graph, u, v) == 1 &&
           heavy_degree(graph, is_heavy, u) >= 2 &&
           heavy_degree(graph, is_heavy, v) >= 2 &&
           is_bridge(graph, u, v);
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

inline std::array<int, 3> choose_gpa(const atom* atoms,
                                     int first_atom,
                                     const BondGraph& graph,
                                     const std::vector<bool>& is_heavy) {
    const int last = std::max(0, static_cast<int>(graph.size()) - 1);
    std::array<int, 3> best{{0, std::min(1, last), std::min(2, last)}};
    double best_score = -1.0;

    for (int center = 0; center < static_cast<int>(graph.size()); ++center) {
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
                    (!is_rotatable(graph, is_heavy, left, center) ? 1 : 0) +
                    (!is_rotatable(graph, is_heavy, center, right) ? 1 : 0);
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
    return best;
}

inline ReconstructionTree build_tree(atom* atoms,
                                     resid& ligand,
                                     int first_atom,
                                     const BondGraph& graph,
                                     const std::vector<bool>& is_heavy) {
    ReconstructionTree tree;
    const int n = static_cast<int>(graph.size());
    tree.parent.assign(n, -1);
    tree.visited.assign(n, false);
    tree.gpa = choose_gpa(atoms, first_atom, graph, is_heavy);

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
    return tree;
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
        if (!is_rotatable(graph, is_heavy, parent, child)) continue;

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
