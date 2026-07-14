// ChecksProtein.cpp — Protein–ligand geometry checks for FlexAIDdS PoseBust
//
// Clean-room C++26 (Apache-2.0). No third-party PoseBusters source was used.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "ChecksProtein.h"
#include "ChecksGeometry.h"  // shared vdw_radius() definition

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace flexaids::posebust {
namespace {

// ── thresholds (dock-suite style) ──────────────────────────────────────────
constexpr float kMinDistanceA            = 1.5f;   // Å hard clash floor
constexpr float kVdwClashScale           = 0.75f;  // soft vdW clash factor
constexpr float kCellSizeA               = 5.0f;   // Å spatial-hash cell
constexpr float kPocketPresenceA         = 5.0f;   // Å pocket presence cutoff
constexpr float kGridSpacingA            = 0.5f;   // Å voxel pitch
constexpr float kBBoxExpandA             = 2.0f;   // Å ligand bbox pad
constexpr double kMaxVolumeOverlap       = 0.075;  // 7.5 %
constexpr std::size_t kBruteForcePairCap = 5'000'000;

// ── helpers ────────────────────────────────────────────────────────────────

[[nodiscard]] inline bool is_heavy_atom(const Atom& a) noexcept {
    if (a.is_h) return false;
    return a.atomic_num > 1;
}

[[nodiscard]] inline Vec3 atom_pos(const Atom& a) noexcept {
    return Vec3{a.x, a.y, a.z};
}

struct HeavyRef {
    Vec3  pos{};
    int   Z    = 0;
    float rvdw = 0.f;
};

[[nodiscard]] std::vector<HeavyRef> collect_heavy(const Molecule& mol) {
    std::vector<HeavyRef> out;
    out.reserve(mol.atoms.size());
    for (const Atom& a : mol.atoms) {
        if (!is_heavy_atom(a)) continue;
        HeavyRef h;
        h.pos  = atom_pos(a);
        h.Z    = a.atomic_num;
        h.rvdw = vdw_radius(a.atomic_num);
        out.push_back(h);
    }
    return out;
}

// Sparse 3-D cell list over heavy atoms (cell edge = kCellSizeA).
class CellList {
public:
    void build(const std::vector<HeavyRef>& atoms) {
        atoms_ = &atoms;
        if (atoms.empty()) {
            dim_[0] = dim_[1] = dim_[2] = 0;
            return;
        }

        float bmin[3] = {atoms[0].pos.x, atoms[0].pos.y, atoms[0].pos.z};
        float bmax[3] = {bmin[0], bmin[1], bmin[2]};
        for (const auto& a : atoms) {
            bmin[0] = std::min(bmin[0], a.pos.x);
            bmin[1] = std::min(bmin[1], a.pos.y);
            bmin[2] = std::min(bmin[2], a.pos.z);
            bmax[0] = std::max(bmax[0], a.pos.x);
            bmax[1] = std::max(bmax[1], a.pos.y);
            bmax[2] = std::max(bmax[2], a.pos.z);
        }
        for (int d = 0; d < 3; ++d) {
            origin_[d] = bmin[d] - kCellSizeA;
            dim_[d]    = static_cast<int>(std::floor(
                          (bmax[d] + kCellSizeA - origin_[d]) / kCellSizeA)) +
                      1;
            dim_[d] = std::max(dim_[d], 1);
        }

        const std::size_t ncells = static_cast<std::size_t>(dim_[0]) *
                                   static_cast<std::size_t>(dim_[1]) *
                                   static_cast<std::size_t>(dim_[2]);
        cells_.assign(ncells, {});
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            const int c = cell_index(atoms[i].pos);
            if (c >= 0) cells_[static_cast<std::size_t>(c)].push_back(i);
        }
    }

    template <class Fn>
    void for_neighbors(const Vec3& p, Fn&& fn) const {
        if (!atoms_ || atoms_->empty() || dim_[0] == 0) return;

        const int cx =
            static_cast<int>(std::floor((p.x - origin_[0]) / kCellSizeA));
        const int cy =
            static_cast<int>(std::floor((p.y - origin_[1]) / kCellSizeA));
        const int cz =
            static_cast<int>(std::floor((p.z - origin_[2]) / kCellSizeA));

        for (int dx = -1; dx <= 1; ++dx) {
            const int ix = cx + dx;
            if (ix < 0 || ix >= dim_[0]) continue;
            for (int dy = -1; dy <= 1; ++dy) {
                const int iy = cy + dy;
                if (iy < 0 || iy >= dim_[1]) continue;
                for (int dz = -1; dz <= 1; ++dz) {
                    const int iz = cz + dz;
                    if (iz < 0 || iz >= dim_[2]) continue;
                    const std::size_t ci = flat(ix, iy, iz);
                    for (std::size_t idx : cells_[ci]) {
                        fn((*atoms_)[idx]);
                    }
                }
            }
        }
    }

private:
    [[nodiscard]] std::size_t flat(int ix, int iy, int iz) const noexcept {
        return static_cast<std::size_t>(ix) +
               static_cast<std::size_t>(dim_[0]) *
                   (static_cast<std::size_t>(iy) +
                    static_cast<std::size_t>(dim_[1]) *
                        static_cast<std::size_t>(iz));
    }

    [[nodiscard]] int cell_index(const Vec3& p) const noexcept {
        const int ix =
            static_cast<int>(std::floor((p.x - origin_[0]) / kCellSizeA));
        const int iy =
            static_cast<int>(std::floor((p.y - origin_[1]) / kCellSizeA));
        const int iz =
            static_cast<int>(std::floor((p.z - origin_[2]) / kCellSizeA));
        if (ix < 0 || iy < 0 || iz < 0 || ix >= dim_[0] || iy >= dim_[1] ||
            iz >= dim_[2]) {
            return -1;
        }
        return static_cast<int>(flat(ix, iy, iz));
    }

    const std::vector<HeavyRef>*          atoms_ = nullptr;
    float                                 origin_[3] = {0.f, 0.f, 0.f};
    int                                   dim_[3]    = {0, 0, 0};
    std::vector<std::vector<std::size_t>> cells_;
};

struct DistClashResult {
    float min_dist     = std::numeric_limits<float>::infinity();
    bool  has_pair     = false;
    bool  no_vdw_clash = true;
    int   n_pairs      = 0;
    int   n_clashes    = 0;
};

[[nodiscard]] DistClashResult compute_distance_clash(
    const std::vector<HeavyRef>& lig,
    const std::vector<HeavyRef>& prot) {
    DistClashResult r;
    if (lig.empty() || prot.empty()) {
        r.min_dist     = std::numeric_limits<float>::infinity();
        r.has_pair     = false;
        r.no_vdw_clash = true;
        return r;
    }

    const std::size_t n = lig.size();
    const std::size_t m = prot.size();
    const bool use_brute =
        (m > 0 && n <= (std::numeric_limits<std::size_t>::max() / m) &&
         n * m < kBruteForcePairCap);

    auto consider = [&](const HeavyRef& L, const HeavyRef& P) {
        const float d = dist(L.pos, P.pos);
        r.has_pair    = true;
        ++r.n_pairs;
        if (d < r.min_dist) r.min_dist = d;
        const float thresh = kVdwClashScale * (L.rvdw + P.rvdw);
        if (d < thresh) {
            r.no_vdw_clash = false;
            ++r.n_clashes;
        }
    };

    if (use_brute) {
        for (const auto& L : lig) {
            for (const auto& P : prot) {
                consider(L, P);
            }
        }
    } else {
        // Cell list on the larger set. With 5 Å cells any pair closer than
        // ~5 Å (hard 1.5 Å clash / soft vdW clash) lies in the 3×3×3 hood.
        const bool protein_is_base = m >= n;
        CellList   cells;
        if (protein_is_base) {
            cells.build(prot);
            for (const auto& L : lig) {
                cells.for_neighbors(L.pos, [&](const HeavyRef& P) {
                    consider(L, P);
                });
            }
        } else {
            cells.build(lig);
            for (const auto& P : prot) {
                cells.for_neighbors(P.pos, [&](const HeavyRef& L) {
                    consider(L, P);
                });
            }
        }
        if (!r.has_pair) {
            // Nothing inside any neighbourhood ⇒ separation ≳ cell size.
            r.min_dist     = kCellSizeA;
            r.no_vdw_clash = true;
        }
    }

    return r;
}

[[nodiscard]] double compute_volume_overlap_fraction(
    const std::vector<HeavyRef>& lig,
    const std::vector<HeavyRef>& prot) {
    if (lig.empty()) return 0.0;

    float bmin[3] = {lig[0].pos.x, lig[0].pos.y, lig[0].pos.z};
    float bmax[3] = {bmin[0], bmin[1], bmin[2]};
    for (const auto& a : lig) {
        bmin[0] = std::min(bmin[0], a.pos.x);
        bmin[1] = std::min(bmin[1], a.pos.y);
        bmin[2] = std::min(bmin[2], a.pos.z);
        bmax[0] = std::max(bmax[0], a.pos.x);
        bmax[1] = std::max(bmax[1], a.pos.y);
        bmax[2] = std::max(bmax[2], a.pos.z);
    }
    for (int d = 0; d < 3; ++d) {
        bmin[d] -= kBBoxExpandA;
        bmax[d] += kBBoxExpandA;
    }

    const int nx = std::max(
        1, static_cast<int>(std::ceil((bmax[0] - bmin[0]) / kGridSpacingA)) + 1);
    const int ny = std::max(
        1, static_cast<int>(std::ceil((bmax[1] - bmin[1]) / kGridSpacingA)) + 1);
    const int nz = std::max(
        1, static_cast<int>(std::ceil((bmax[2] - bmin[2]) / kGridSpacingA)) + 1);

    CellList prot_cells;
    if (!prot.empty()) prot_cells.build(prot);

    std::uint64_t lig_voxels  = 0;
    std::uint64_t both_voxels = 0;

    for (int ix = 0; ix < nx; ++ix) {
        const float x = bmin[0] + static_cast<float>(ix) * kGridSpacingA;
        for (int iy = 0; iy < ny; ++iy) {
            const float y = bmin[1] + static_cast<float>(iy) * kGridSpacingA;
            for (int iz = 0; iz < nz; ++iz) {
                const float z = bmin[2] + static_cast<float>(iz) * kGridSpacingA;
                const Vec3  v{x, y, z};

                bool in_lig = false;
                for (const auto& L : lig) {
                    if (dist2(v, L.pos) <= L.rvdw * L.rvdw) {
                        in_lig = true;
                        break;
                    }
                }
                if (!in_lig) continue;
                ++lig_voxels;

                bool in_prot = false;
                if (!prot.empty()) {
                    prot_cells.for_neighbors(v, [&](const HeavyRef& P) {
                        if (in_prot) return;
                        if (dist2(v, P.pos) <= P.rvdw * P.rvdw) in_prot = true;
                    });
                    // Full scan when protein is modest (correctness over speed).
                    if (!in_prot && prot.size() < 20'000ull) {
                        for (const auto& P : prot) {
                            if (dist2(v, P.pos) <= P.rvdw * P.rvdw) {
                                in_prot = true;
                                break;
                            }
                        }
                    }
                }
                if (in_prot) ++both_voxels;
            }
        }
    }

    if (lig_voxels == 0) return 0.0;
    return static_cast<double>(both_voxels) / static_cast<double>(lig_voxels);
}

[[nodiscard]] bool has_pocket_presence(const std::vector<HeavyRef>& lig,
                                       const std::vector<HeavyRef>& prot) {
    if (lig.empty() || prot.empty()) return false;
    const float cut2 = kPocketPresenceA * kPocketPresenceA;

    const std::size_t n = lig.size();
    const std::size_t m = prot.size();
    const bool use_brute =
        m > 0 && n <= (std::numeric_limits<std::size_t>::max() / m) &&
        n * m < kBruteForcePairCap;

    if (use_brute) {
        for (const auto& L : lig) {
            for (const auto& P : prot) {
                if (dist2(L.pos, P.pos) <= cut2) return true;
            }
        }
        return false;
    }

    CellList cells;
    cells.build(prot);
    for (const auto& L : lig) {
        bool found = false;
        cells.for_neighbors(L.pos, [&](const HeavyRef& P) {
            if (!found && dist2(L.pos, P.pos) <= cut2) found = true;
        });
        if (found) return true;
    }
    return false;
}

}  // namespace

// ── public API ─────────────────────────────────────────────────────────────
// vdw_radius() is defined once in ChecksGeometry.cpp

void check_intermolecular_distance(const Molecule& ligand,
                                   const Molecule& protein,
                                   std::vector<CheckItem>& out) {
    const auto lig  = collect_heavy(ligand);
    const auto prot = collect_heavy(protein);
    const DistClashResult r = compute_distance_clash(lig, prot);

    // Dock-suite key: no_clashes — pass when min lig–prot heavy distance is OK.
    const bool min_ok =
        !r.has_pair ||
        (std::isfinite(r.min_dist) && r.min_dist >= kMinDistanceA);

    {
        std::ostringstream oss;
        if (!r.has_pair) {
            oss << "no heavy-heavy pair (lig_heavy=" << lig.size()
                << ", prot_heavy=" << prot.size() << "); vacuous pass";
        } else {
            oss << "min_dist=" << r.min_dist << " A (threshold "
                << kMinDistanceA << " A); pairs_examined=" << r.n_pairs
                << "; soft_vdw_clashes=" << r.n_clashes
                << " (scale=" << kVdwClashScale << ")";
        }
        CheckItem item;
        item.key       = "no_clashes";
        item.label     = "No ligand–protein clashes";
        item.passed    = min_ok;
        item.detail    = oss.str();
        item.metric    = r.has_pair ? r.min_dist : std::numeric_limits<float>::quiet_NaN();
        item.threshold = kMinDistanceA;
        item.n_checked = r.n_pairs;
        item.n_failed  = min_ok ? 0 : 1;
        out.push_back(std::move(item));
    }
}

void check_volume_overlap(const Molecule& ligand,
                          const Molecule& protein,
                          std::vector<CheckItem>& out) {
    const auto lig  = collect_heavy(ligand);
    const auto prot = collect_heavy(protein);

    const double frac   = compute_volume_overlap_fraction(lig, prot);
    const bool   vol_ok = frac <= kMaxVolumeOverlap;
    {
        std::ostringstream oss;
        oss << "overlap_fraction=" << frac << " (threshold " << kMaxVolumeOverlap
            << "); grid=" << kGridSpacingA << " A; bbox_pad=" << kBBoxExpandA
            << " A; lig_heavy=" << lig.size() << "; prot_heavy=" << prot.size();
        CheckItem item;
        item.key       = "no_volume_clash";
        item.label     = "No ligand–protein volume clash";
        item.passed    = vol_ok;
        item.detail    = oss.str();
        item.metric    = static_cast<float>(frac);
        item.threshold = static_cast<float>(kMaxVolumeOverlap);
        out.push_back(std::move(item));
    }

    // not_too_far_away: pocket presence — some protein heavy within 5 Å.
    const bool pocket = has_pocket_presence(lig, prot);
    {
        std::ostringstream oss;
        oss << (pocket ? "protein heavy atom within " : "no protein heavy atom within ")
            << kPocketPresenceA << " A of ligand; lig_heavy=" << lig.size()
            << "; prot_heavy=" << prot.size();
        CheckItem item;
        item.key       = "not_too_far_away";
        item.label     = "Ligand not too far from protein";
        item.passed    = pocket;
        item.detail    = oss.str();
        item.metric    = pocket ? 0.f : kPocketPresenceA;
        item.threshold = kPocketPresenceA;
        out.push_back(std::move(item));
    }
}

}  // namespace flexaids::posebust
