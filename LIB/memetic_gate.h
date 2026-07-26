// Wave 3.4 memetic enable gate (header-only, unit-testable).
// Product path stays OFF unless FLEXAIDDS_MEMETIC is truthy AND a burial/scoring
// pilot unlock is set.
//
// Unlock keys (either is sufficient; both may be set):
//   FLEXAIDDS_PB_CLASH_PHASE2_PASS=1  — revised Phase 2 (ROADMAP_v2): SCORING-LOCKED
//                                      pb_clash oracle with magnitude floor PASS
//   FLEXAIDDS_WALL_PILOT_PASS=1      — legacy wall oracle (structurally unpassable;
//                                      kept so old docs/tests still parse; do not
//                                      set from WAL evidence)
//
// Never enable memetic from micro-ΔdCF SEARCH-MISS oracles or sign-only noise.
#pragma once

#include <cstdlib>

namespace flexaids {

/// Returns 1 iff FLEXAIDDS_MEMETIC is a non-zero integer AND at least one pilot
/// unlock env is a non-zero integer; else 0.
inline int resolve_use_memetic_from_env() {
    const char* e = std::getenv("FLEXAIDDS_MEMETIC");
    const bool want = e != nullptr && e[0] != '\0' && std::atoi(e) != 0;
    if (!want) {
        return 0;
    }
    const char* phase2 = std::getenv("FLEXAIDDS_PB_CLASH_PHASE2_PASS");
    const bool phase2_pass =
        phase2 != nullptr && phase2[0] != '\0' && std::atoi(phase2) != 0;
    const char* wall_ok = std::getenv("FLEXAIDDS_WALL_PILOT_PASS");
    const bool wall_pass =
        wall_ok != nullptr && wall_ok[0] != '\0' && std::atoi(wall_ok) != 0;
    return (phase2_pass || wall_pass) ? 1 : 0;
}

}  // namespace flexaids
