#pragma once

#include "RngSeed.h"
#include <cstdint>
#include <mutex>

namespace flexaids {

// Region GAs still reach nested shared state and process-global seed epochs.
// This is the production manager's dispatcher, not a claim that arbitrary
// callers of GA() elsewhere in the process are safe to run concurrently.
inline std::mutex& region_execution_mutex() {
    static std::mutex mutex;
    return mutex;
}

template<class RunRegion>
void for_each_serial_region(int begin, int end, int stride, RunRegion&& run_region) {
    std::lock_guard<std::mutex> lock(region_execution_mutex());
    for (int region = begin; region < end; region += stride)
        run_region(region);
}

// A prior region GA advances the process master seed. An explicit parent
// protocol must take precedence over that mutable fallback on every run.
inline bool resolve_region_parent_seed(int configured_seed, std::uint64_t& seed) {
    if (configured_seed != 0) {
        seed = static_cast<unsigned int>(configured_seed);
        return true;
    }
    if (flexaids_rng::env_seed(seed)) return true;
    if (!flexaids_rng::has_master_seed()) return false;
    seed = flexaids_rng::master_seed();
    return true;
}

// A pure function of parent seed and global region index: invariant to worker
// scheduling and MPI rank assignment. GB_Global::seed needs a nonzero int.
inline int region_seed(std::uint64_t parent_seed, int region_index) {
    const auto mixed = flexaids_rng::splitmix64(
        parent_seed + static_cast<std::uint64_t>(region_index));
    return static_cast<int>(mixed % 2147483647ULL) + 1;
}

}  // namespace flexaids
