// niche_hash.h — spatial hash for Cartesian niche sharing (default OFF)
//
// FLEXAIDDS_NICHE_HASH=1 plus Cartesian niche: only compare chromosomes whose
// ligand centroids fall in neighboring cells of size sig_share. Cell size ==
// sigma so any pair with RMSD <= sigma has centroids in the 27-neighborhood.
// Unset → callers keep the O(pop²) loop.
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "EnvFlags.h"

#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace flexaids {

inline bool niche_hash_enabled() noexcept
{
    return env_bool("FLEXAIDDS_NICHE_HASH", false);
}

struct NicheCell {
    int x = 0, y = 0, z = 0;
    bool operator==(const NicheCell& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};

struct NicheCellHash {
    std::size_t operator()(const NicheCell& c) const noexcept {
        const std::uint64_t h =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.x)) * 73856093ull)
            ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.y)) * 19349663ull)
            ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.z)) * 83492791ull);
        return static_cast<std::size_t>(h);
    }
};

inline NicheCell niche_cell_of(float cx, float cy, float cz, float cell)
{
    if (!(cell > 0.f)) cell = 1.f;
    return {static_cast<int>(std::floor(cx / cell)),
            static_cast<int>(std::floor(cy / cell)),
            static_cast<int>(std::floor(cz / cell))};
}

/// Build cell → chromosome indices from per-chrom centroids (xyz stride 3).
inline std::unordered_map<NicheCell, std::vector<int>, NicheCellHash>
niche_hash_build(const float* centroids_xyz, int n_chrom, float cell)
{
    std::unordered_map<NicheCell, std::vector<int>, NicheCellHash> map;
    map.reserve(static_cast<std::size_t>(n_chrom));
    for (int i = 0; i < n_chrom; ++i) {
        const float* p = centroids_xyz + static_cast<std::size_t>(i) * 3;
        map[niche_cell_of(p[0], p[1], p[2], cell)].push_back(i);
    }
    return map;
}

/// Append 27-neighborhood occupants of cell c into `out` (may include `self`).
inline void niche_hash_neighbors(
    const std::unordered_map<NicheCell, std::vector<int>, NicheCellHash>& map,
    NicheCell c, std::vector<int>& out)
{
    out.clear();
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
                NicheCell n{c.x + dx, c.y + dy, c.z + dz};
                auto it = map.find(n);
                if (it == map.end()) continue;
                out.insert(out.end(), it->second.begin(), it->second.end());
            }
}

}  // namespace flexaids
