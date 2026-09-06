#pragma once

#include "RngSeed.h"

#include <cstdint>
#include <string>

// CMake supplies these values. A build without that information must not
// claim a known revision or a clean tree (0=clean, 1=dirty, 2=unknown).
#ifndef FLEXAIDS_GIT_COMMIT
#define FLEXAIDS_GIT_COMMIT "unknown"
#endif
#ifndef FLEXAIDS_GIT_DIRTY
#define FLEXAIDS_GIT_DIRTY 2
#endif

namespace flexaids::pose_provenance {

inline std::string format_remark(const char* commit, int dirty,
                                 std::uint64_t seed)
{
    return std::string("REMARK FLEXAID.commit=") + commit +
           " FLEXAID.dirty=" + std::to_string(dirty) +
           " FLEXAID.seed=" + std::to_string(seed) + "\n";
}

// Read the effective GA seed; do not initialize or draw from any RNG here.
// Keep the existing BindingMode fallback of zero when no master seed is set.
inline std::string remark()
{
    return format_remark(FLEXAIDS_GIT_COMMIT, FLEXAIDS_GIT_DIRTY,
                         flexaids_rng::has_master_seed()
                             ? flexaids_rng::master_seed() : 0);
}

// Classic and density-peak writers already fill a bounded scientific REMARK
// buffer. Insert metadata only after that buffer is complete: adding it to
// the bounded buffer would displace/truncate scientific lines at its tail.
// The expanded string preserves every input byte, including an incomplete
// final line, and keeps provenance near the head of the output.
inline std::string add_to_remarks(std::string scientific_remarks)
{
    const auto newline = scientific_remarks.find('\n');
    scientific_remarks.insert(newline == std::string::npos ? 0 : newline + 1,
                              remark());
    return scientific_remarks;
}

} // namespace flexaids::pose_provenance
