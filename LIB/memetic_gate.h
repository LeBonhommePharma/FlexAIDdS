// Wave 3.4 memetic enable gate (header-only, unit-testable).
// Product path stays OFF unless BOTH env vars are truthy integers.
// WALL_PILOT_PASS must only be set after score-only wall oracle PASS.
#pragma once

#include <cstdlib>

namespace flexaids {

/// Returns 1 iff FLEXAIDDS_MEMETIC is set to a non-zero integer AND
/// FLEXAIDDS_WALL_PILOT_PASS is set to a non-zero integer; else 0.
inline int resolve_use_memetic_from_env() {
    const char* e = std::getenv("FLEXAIDDS_MEMETIC");
    const bool want = e != nullptr && e[0] != '\0' && std::atoi(e) != 0;
    const char* wall_ok = std::getenv("FLEXAIDDS_WALL_PILOT_PASS");
    const bool wall_pass =
        wall_ok != nullptr && wall_ok[0] != '\0' && std::atoi(wall_ok) != 0;
    return (want && wall_pass) ? 1 : 0;
}

}  // namespace flexaids
