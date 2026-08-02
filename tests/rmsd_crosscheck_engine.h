// Cross-check fixture: expose the ENGINE's symmetry-corrected RMSD
// (LIB/calc_rmsd.cpp calc_Hungarian_RMSD) behind a clean, element-labelled
// interface so it can be compared against dataset::hungarian_rmsd in a gtest
// translation unit WITHOUT pulling in flexaid.h (which #defines the macro E).
#pragma once

#include <array>
#include <string>
#include <vector>

namespace crosscheck {

struct AtomSpec {
    int type;                        // FlexAID/SYBYL atom type (>= 1; DUMMY row is one type)
    std::string element;             // element symbol ("C", "N", "O", ...)
    std::array<float, 3> coor;       // docked pose coordinate
    std::array<float, 3> coor_ref;   // crystal / reference coordinate
};

// Symmetry-corrected RMSD via the engine's calc_Hungarian_RMSD, grouping atoms
// by AtomSpec::type (the force-field/SYBYL type). Defined in
// rmsd_crosscheck_engine.cpp, the only TU that includes flexaid.h.
float engine_hungarian_rmsd(const std::vector<AtomSpec>& atoms);

}  // namespace crosscheck
