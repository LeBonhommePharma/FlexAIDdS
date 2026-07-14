// ChecksGeometry.cpp — clean-room geometry plausibility for FlexAIDdS PoseBust
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "ChecksGeometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <queue>
#include <sstream>
#include <utility>

namespace flexaids::posebust {
namespace {

// ---------------------------------------------------------------------------
// Covalent radii (Å) — Cordero et al. single-bond style table (common organics)
// ---------------------------------------------------------------------------

float covalent_radius(int Z) noexcept {
    switch (Z) {
        case 1:  return 0.31f; // H
        case 5:  return 0.84f; // B
        case 6:  return 0.76f; // C
        case 7:  return 0.71f; // N
        case 8:  return 0.66f; // O
        case 9:  return 0.57f; // F
        case 14: return 1.11f; // Si
        case 15: return 1.07f; // P
        case 16: return 1.05f; // S
        case 17: return 1.02f; // Cl
        case 33: return 1.19f; // As
        case 34: return 1.20f; // Se
        case 35: return 1.20f; // Br
        case 53: return 1.39f; // I
        default: return 1.00f;
    }
}

// ---------------------------------------------------------------------------
// Graph helpers
// ---------------------------------------------------------------------------

std::vector<std::vector<int>> build_adjacency(const Molecule& mol) {
    const int n = static_cast<int>(mol.atoms.size());
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(n));
    for (const Bond& b : mol.bonds) {
        if (b.a < 0 || b.b < 0 || b.a >= n || b.b >= n || b.a == b.b) {
            continue;
        }
        adj[static_cast<std::size_t>(b.a)].push_back(b.b);
        adj[static_cast<std::size_t>(b.b)].push_back(b.a);
    }
    for (auto& nbrs : adj) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
    return adj;
}

/// All-pairs graph distance; disconnected → large sentinel.
std::vector<std::vector<int>> all_pairs_graph_distance(
    const std::vector<std::vector<int>>& adj) {
    const int n = static_cast<int>(adj.size());
    constexpr int kInf = 1'000'000;
    std::vector<std::vector<int>> dist(
        static_cast<std::size_t>(n),
        std::vector<int>(static_cast<std::size_t>(n), kInf));

    for (int s = 0; s < n; ++s) {
        dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(s)] = 0;
        std::queue<int> q;
        q.push(s);
        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            const int du = dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(u)];
            for (int v : adj[static_cast<std::size_t>(u)]) {
                auto& dv = dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(v)];
                if (dv > du + 1) {
                    dv = du + 1;
                    q.push(v);
                }
            }
        }
    }
    return dist;
}

// ---------------------------------------------------------------------------
// Angle / torsion
// ---------------------------------------------------------------------------

float bond_angle_degrees(const Vec3& i, const Vec3& j, const Vec3& k) noexcept {
    Vec3 v1 = i - j;
    Vec3 v2 = k - j;
    const float n1 = norm(v1);
    const float n2 = norm(v2);
    if (n1 < 1e-12f || n2 < 1e-12f) {
        return 0.f;
    }
    v1 = v1 * (1.f / n1);
    v2 = v2 * (1.f / n2);
    float c = dot(v1, v2);
    c = std::clamp(c, -1.f, 1.f);
    return std::acos(c) * (180.f / std::numbers::pi_v<float>);
}

/// Ideal valence angle from hybridization heuristic on graph degree of j.
float ideal_angle_for_degree(int degree) noexcept {
    if (degree >= 4) {
        return 109.5f; // tetrahedral / sp3
    }
    if (degree == 3) {
        return 120.0f; // trigonal / sp2
    }
    if (degree == 2) {
        return 180.0f; // linear / sp
    }
    return 0.f; // no angle at degree < 2
}

/// Signed torsion angle (degrees) for a-b-c-d (IUPAC convention via atan2).
float torsion_degrees(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) noexcept {
    const Vec3 b1 = b - a;
    const Vec3 b2 = c - b;
    const Vec3 b3 = d - c;
    const Vec3 n1 = cross(b1, b2);
    const Vec3 n2 = cross(b2, b3);
    const float n1n = norm(n1);
    const float n2n = norm(n2);
    const float b2n = norm(b2);
    if (n1n < 1e-12f || n2n < 1e-12f || b2n < 1e-12f) {
        return 0.f;
    }
    const Vec3 n1u = n1 * (1.f / n1n);
    const Vec3 n2u = n2 * (1.f / n2n);
    const Vec3 b2u = b2 * (1.f / b2n);
    const Vec3 m1  = cross(n1u, b2u);
    const float x  = dot(n1u, n2u);
    const float y  = dot(m1, n2u);
    return std::atan2(y, x) * (180.f / std::numbers::pi_v<float>);
}

// ---------------------------------------------------------------------------
// Best-fit plane distance (SVD-free: covariance eigen via analytic 3×3)
// ---------------------------------------------------------------------------

/// Max absolute distance of points to the best-fit plane (centroid + PCA normal).
float max_out_of_plane_distance(const std::vector<Vec3>& pts) {
    if (pts.size() < 3) {
        return 0.f;
    }
    Vec3 centroid{0.f, 0.f, 0.f};
    for (const Vec3& p : pts) {
        centroid = centroid + p;
    }
    centroid = centroid * (1.f / static_cast<float>(pts.size()));

    // Covariance matrix (symmetric)
    double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
    for (const Vec3& p : pts) {
        const double dx = p.x - centroid.x;
        const double dy = p.y - centroid.y;
        const double dz = p.z - centroid.z;
        xx += dx * dx;
        xy += dx * dy;
        xz += dx * dz;
        yy += dy * dy;
        yz += dy * dz;
        zz += dz * dz;
    }

    // Power iteration for smallest eigenvector of cov (plane normal).
    // Start with cross product of two long edges for a stable seed.
    Vec3 n = cross(pts[1] - pts[0], pts[2] - pts[0]);
    const float n0 = norm(n);
    if (n0 * n0 < 1e-12f) {
        n = {0.f, 0.f, 1.f};
    } else {
        n = n * (1.f / n0);
    }

    // Inverse-power style: repeatedly apply (trace*I − cov) to amplify small eig.
    const double tr = xx + yy + zz;
    for (int iter = 0; iter < 32; ++iter) {
        const double nx = n.x, ny = n.y, nz = n.z;
        // v' = (tr I − C) v
        double vx = (tr - xx) * nx - xy * ny - xz * nz;
        double vy = -xy * nx + (tr - yy) * ny - yz * nz;
        double vz = -xz * nx - yz * ny + (tr - zz) * nz;
        const double nn = std::sqrt(vx * vx + vy * vy + vz * vz);
        if (nn < 1e-18) {
            break;
        }
        n = {static_cast<float>(vx / nn), static_cast<float>(vy / nn),
             static_cast<float>(vz / nn)};
    }

    float max_d = 0.f;
    for (const Vec3& p : pts) {
        const float d = std::fabs(dot(p - centroid, n));
        max_d = std::max(max_d, d);
    }
    return max_d;
}

// ---------------------------------------------------------------------------
// Simple cycle enumeration (size 5–6), de-duplicated
// ---------------------------------------------------------------------------

using Cycle = std::vector<int>;

Cycle normalize_cycle(Cycle c) {
    if (c.empty()) {
        return c;
    }
    // Rotate so minimum atom index is first.
    const auto min_it = std::min_element(c.begin(), c.end());
    std::rotate(c.begin(), min_it, c.end());
    // Choose orientation with smaller second element.
    Cycle rev = c;
    std::reverse(rev.begin() + 1, rev.end());
    if (rev < c) {
        return rev;
    }
    return c;
}

void dfs_cycles(int start,
                int u,
                int depth,
                int max_depth,
                const std::vector<std::vector<int>>& adj,
                std::vector<int>& path,
                std::vector<char>& on_path,
                std::vector<Cycle>& out) {
    if (depth > max_depth) {
        return;
    }
    for (int v : adj[static_cast<std::size_t>(u)]) {
        if (depth >= 4 && depth <= 5 && v == start) {
            // Closed cycle of size depth+1 ∈ {5,6}
            out.push_back(normalize_cycle(path));
            continue;
        }
        if (on_path[static_cast<std::size_t>(v)]) {
            continue;
        }
        // Only explore v ≥ start to cut rotations early (still normalize).
        if (v < start && depth > 0) {
            continue;
        }
        on_path[static_cast<std::size_t>(v)] = 1;
        path.push_back(v);
        dfs_cycles(start, v, depth + 1, max_depth, adj, path, on_path, out);
        path.pop_back();
        on_path[static_cast<std::size_t>(v)] = 0;
    }
}

std::vector<Cycle> find_simple_cycles_5_6(const std::vector<std::vector<int>>& adj) {
    const int n = static_cast<int>(adj.size());
    std::vector<Cycle> raw;
    std::vector<int> path;
    std::vector<char> on_path(static_cast<std::size_t>(n), 0);

    for (int s = 0; s < n; ++s) {
        path.clear();
        path.push_back(s);
        on_path[static_cast<std::size_t>(s)] = 1;
        dfs_cycles(s, s, 0, 5, adj, path, on_path, raw);
        on_path[static_cast<std::size_t>(s)] = 0;
    }

    std::sort(raw.begin(), raw.end());
    raw.erase(std::unique(raw.begin(), raw.end()), raw.end());
    // Keep only size 5–6 (normalize preserves size).
    raw.erase(std::remove_if(raw.begin(), raw.end(),
                             [](const Cycle& c) {
                                 return c.size() < 5 || c.size() > 6;
                             }),
              raw.end());
    return raw;
}

bool majority_aromatic_or_sp2(const Molecule& /*mol*/,
                              const Cycle& cyc,
                              const std::vector<std::vector<int>>& adj) {
    int votes = 0;
    for (int idx : cyc) {
        const int deg = static_cast<int>(adj[static_cast<std::size_t>(idx)].size());
        // Degree-3 trigonal / sp2 heuristic (no aromatic flag on Atom in Types.h).
        if (deg == 3) {
            ++votes;
        }
    }
    return votes * 2 > static_cast<int>(cyc.size()); // strict majority
}

// ---------------------------------------------------------------------------
// Emit helper (Types.h: key / label / passed / detail / metric / …)
// ---------------------------------------------------------------------------

CheckItem make_item(std::string key,
                    std::string label,
                    bool passed,
                    int n_checked,
                    int n_failed,
                    double metric,
                    double threshold,
                    const std::string& detail) {
    CheckItem item;
    item.key       = std::move(key);
    item.label     = std::move(label);
    item.passed    = passed;
    item.n_checked = n_checked;
    item.n_failed  = n_failed;
    item.metric    = static_cast<float>(metric);
    item.threshold = static_cast<float>(threshold);
    item.detail    = detail;
    return item;
}

} // namespace

// ---------------------------------------------------------------------------
// Public: Bondi-ish vdW radii
// ---------------------------------------------------------------------------

float vdw_radius(int Z) noexcept {
    // Bondi (1964) van der Waals radii (Å); common organics + defaults.
    switch (Z) {
        case 1:  return 1.20f; // H
        case 5:  return 1.92f; // B (approx.)
        case 6:  return 1.70f; // C
        case 7:  return 1.55f; // N
        case 8:  return 1.52f; // O
        case 9:  return 1.47f; // F
        case 14: return 2.10f; // Si
        case 15: return 1.80f; // P
        case 16: return 1.80f; // S
        case 17: return 1.75f; // Cl
        case 33: return 1.85f; // As
        case 34: return 1.90f; // Se
        case 35: return 1.85f; // Br
        case 53: return 1.98f; // I
        default: return 2.00f;
    }
}

// ---------------------------------------------------------------------------
// Distance geometry suite
// ---------------------------------------------------------------------------

void check_distance_geometry(const Molecule& pred, std::vector<CheckItem>& out) {
    const auto adj = build_adjacency(pred);
    const int n    = static_cast<int>(pred.atoms.size());

    // --- bond_lengths_within_bounds -----------------------------------------
    {
        int n_checked = 0;
        int n_failed  = 0;
        double worst_rel = 0.0; // |d-ideal|/ideal
        for (const Bond& b : pred.bonds) {
            if (b.a < 0 || b.b < 0 || b.a >= n || b.b >= n) {
                continue;
            }
            const Atom& a1 = pred.atoms[static_cast<std::size_t>(b.a)];
            const Atom& a2 = pred.atoms[static_cast<std::size_t>(b.b)];
            const float r1 = covalent_radius(a1.atomic_num);
            const float r2 = covalent_radius(a2.atomic_num);
            const float ideal = r1 + r2;
            if (ideal < 1e-6f) {
                continue;
            }
            const float d = dist(a1.pos(), a2.pos());
            const float half = kThresholdBadBondLength * ideal;
            const float rel  = std::fabs(d - ideal) / ideal;
            worst_rel = std::max(worst_rel, static_cast<double>(rel));
            ++n_checked;
            if (std::fabs(d - ideal) > half + 1e-6f) {
                ++n_failed;
            }
        }
        const bool passed = (n_failed == 0);
        std::ostringstream oss;
        oss << n_checked << " bonds; " << n_failed << " outside covalent sum ± "
            << kThresholdBadBondLength << "*(r1+r2); worst |Δ|/ideal="
            << worst_rel;
        out.push_back(make_item("bond_lengths_within_bounds",
                                "Bond lengths within bounds",
                                passed,
                                n_checked,
                                n_failed,
                                worst_rel,
                                static_cast<double>(kThresholdBadBondLength),
                                oss.str()));
    }

    // --- bond_angles_within_bounds ------------------------------------------
    {
        int n_checked = 0;
        int n_failed  = 0;
        double worst_rel = 0.0; // |Δθ|/θ_ideal
        for (int j = 0; j < n; ++j) {
            const Atom& center = pred.atoms[static_cast<std::size_t>(j)];
            if (center.is_h) {
                continue; // angles at H ignored
            }
            const auto& nbrs = adj[static_cast<std::size_t>(j)];
            const int degree = static_cast<int>(nbrs.size());
            const float theta_ideal = ideal_angle_for_degree(degree);
            if (theta_ideal <= 0.f || degree < 2) {
                continue;
            }
            for (std::size_t a = 0; a < nbrs.size(); ++a) {
                for (std::size_t b = a + 1; b < nbrs.size(); ++b) {
                    const int i = nbrs[a];
                    const int k = nbrs[b];
                    const Atom& ai = pred.atoms[static_cast<std::size_t>(i)];
                    const Atom& ak = pred.atoms[static_cast<std::size_t>(k)];
                    // Optional: skip pure-H wings (H–X–H and angles with only H).
                    if (ai.is_h && ak.is_h) {
                        continue;
                    }
                    const float theta = bond_angle_degrees(ai.pos(), center.pos(), ak.pos());
                    const float rel   = std::fabs(theta - theta_ideal) / theta_ideal;
                    worst_rel = std::max(worst_rel, static_cast<double>(rel));
                    ++n_checked;
                    if (rel > kThresholdBadAngle + 1e-6f) {
                        ++n_failed;
                    }
                }
            }
        }
        const bool passed = (n_failed == 0);
        std::ostringstream oss;
        oss << n_checked << " angles; " << n_failed
            << " with |Δθ|/θ_ideal > " << kThresholdBadAngle
            << " (deg4→109.5, deg3→120, deg2→180); worst rel=" << worst_rel;
        out.push_back(make_item("bond_angles_within_bounds",
                                "Bond angles within bounds",
                                passed,
                                n_checked,
                                n_failed,
                                worst_rel,
                                static_cast<double>(kThresholdBadAngle),
                                oss.str()));
    }

    // --- no_internal_clash --------------------------------------------------
    {
        const auto gdist = all_pairs_graph_distance(adj);
        int n_checked = 0;
        int n_failed  = 0;
        double worst_ratio = 0.0; // min d / (scale*(ri+rj))  — lower is worse; we track min
        double min_ratio   = std::numeric_limits<double>::infinity();

        for (int i = 0; i < n; ++i) {
            const Atom& ai = pred.atoms[static_cast<std::size_t>(i)];
            if (ai.is_h) {
                continue;
            }
            for (int j = i + 1; j < n; ++j) {
                const Atom& aj = pred.atoms[static_cast<std::size_t>(j)];
                if (aj.is_h) {
                    continue;
                }
                const int gd = gdist[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
                if (gd < kClashMinGraphDistance) {
                    continue; // 1–2, 1–3, 1–4 (and self) excluded
                }
                const float ri = vdw_radius(ai.atomic_num);
                const float rj = vdw_radius(aj.atomic_num);
                const float cutoff = kClashVdwScale * (ri + rj);
                const float d = dist(ai.pos(), aj.pos());
                const double ratio = (cutoff > 1e-12f)
                                         ? static_cast<double>(d / cutoff)
                                         : 1.0;
                min_ratio = std::min(min_ratio, ratio);
                ++n_checked;
                if (d + 1e-6f < cutoff) {
                    ++n_failed;
                }
            }
        }
        if (std::isfinite(min_ratio)) {
            worst_ratio = min_ratio; // metric: smallest d/cutoff (fail if < 1)
        } else {
            worst_ratio = 1.0;
        }
        const bool passed = (n_failed == 0);
        std::ostringstream oss;
        oss << n_checked << " heavy pairs (graph dist≥" << kClashMinGraphDistance
            << "); " << n_failed << " with d < " << kClashVdwScale
            << "*(vdW_i+vdW_j); min d/cutoff=" << worst_ratio;
        out.push_back(make_item("no_internal_clash",
                                "No internal steric clash",
                                passed,
                                n_checked,
                                n_failed,
                                worst_ratio,
                                1.0, // pass when metric ≥ 1
                                oss.str()));
    }
}

// ---------------------------------------------------------------------------
// Flatness suite
// ---------------------------------------------------------------------------

void check_flatness(const Molecule& pred, std::vector<CheckItem>& out) {
    const auto adj = build_adjacency(pred);
    const int n    = static_cast<int>(pred.atoms.size());

    // --- flatness_passes (aromatic / sp2 rings) -----------------------------
    {
        const auto cycles = find_simple_cycles_5_6(adj);
        int n_checked = 0;
        int n_failed  = 0;
        double worst_oop = 0.0;

        for (const Cycle& cyc : cycles) {
            if (!majority_aromatic_or_sp2(pred, cyc, adj)) {
                continue;
            }
            std::vector<Vec3> pts;
            pts.reserve(cyc.size());
            for (int idx : cyc) {
                if (idx < 0 || idx >= n) {
                    continue;
                }
                pts.push_back(pred.atoms[static_cast<std::size_t>(idx)].pos());
            }
            if (pts.size() < 5) {
                continue;
            }
            const float oop = max_out_of_plane_distance(pts);
            worst_oop = std::max(worst_oop, static_cast<double>(oop));
            ++n_checked;
            if (oop > kThresholdFlatnessAngstrom + 1e-6f) {
                ++n_failed;
            }
        }

        const bool passed = (n_failed == 0); // vacuously true if none checked
        std::ostringstream oss;
        oss << n_checked << " aromatic/sp2 rings (size 5–6); " << n_failed
            << " with max OOP > " << kThresholdFlatnessAngstrom
            << " Å; worst OOP=" << worst_oop << " Å";
        out.push_back(make_item("flatness_passes",
                                "Aromatic ring flatness",
                                passed,
                                n_checked,
                                n_failed,
                                worst_oop,
                                static_cast<double>(kThresholdFlatnessAngstrom),
                                oss.str()));
    }

    // --- flat_double_bonds (separate flat key) ------------------------------
    {
        int n_checked = 0;
        int n_failed  = 0;
        double worst_sin = 0.0;

        for (const Bond& b : pred.bonds) {
            if (b.a < 0 || b.b < 0 || b.a >= n || b.b >= n) {
                continue;
            }
            // order ≥ 2 (double/triple); aromatic rings handled above
            if (b.order < 2) {
                continue;
            }

            const int a = b.a;
            const int c = b.b;
            const auto& na = adj[static_cast<std::size_t>(a)];
            const auto& nc = adj[static_cast<std::size_t>(c)];

            // First substituent on each end (prefer heavy atoms).
            auto pick_sub = [](const std::vector<int>& nbrs, int partner,
                               const Molecule& mol) -> int {
                int heavy = -1;
                int any   = -1;
                for (int v : nbrs) {
                    if (v == partner) {
                        continue;
                    }
                    if (any < 0) {
                        any = v;
                    }
                    if (!mol.atoms[static_cast<std::size_t>(v)].is_h) {
                        heavy = v;
                        break;
                    }
                }
                return heavy >= 0 ? heavy : any;
            };

            const int sub_a = pick_sub(na, c, pred);
            const int sub_c = pick_sub(nc, a, pred);
            if (sub_a < 0 || sub_c < 0) {
                // Terminal double bond with insufficient substituents — skip.
                continue;
            }

            const float phi = torsion_degrees(
                pred.atoms[static_cast<std::size_t>(sub_a)].pos(),
                pred.atoms[static_cast<std::size_t>(a)].pos(),
                pred.atoms[static_cast<std::size_t>(c)].pos(),
                pred.atoms[static_cast<std::size_t>(sub_c)].pos());
            const float s = std::fabs(std::sin(phi * std::numbers::pi_v<float> / 180.f));
            worst_sin = std::max(worst_sin, static_cast<double>(s));
            ++n_checked;
            if (s > kThresholdDoubleBondSinPhi + 1e-6f) {
                ++n_failed;
            }
        }

        const bool passed = (n_failed == 0);
        std::ostringstream oss;
        oss << n_checked << " bonds with order≥2; " << n_failed
            << " with |sin(φ)| > " << kThresholdDoubleBondSinPhi
            << "; worst |sin(φ)|=" << worst_sin;
        out.push_back(make_item("flat_double_bonds",
                                "Double-bond planarity",
                                passed,
                                n_checked,
                                n_failed,
                                worst_sin,
                                static_cast<double>(kThresholdDoubleBondSinPhi),
                                oss.str()));
    }
}

} // namespace flexaids::posebust
