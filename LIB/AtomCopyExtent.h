#pragma once

// FlexAID atom[] and residue[] are 1-based: live entries occupy [1, count].
// A half-open C++ copy of `count` elements is [0, count) and DROPS the last
// live entry. Gaboom's ParEvalWS copies `atoms + natm + 1` (LIB/gaboom.cpp).
// ParallelDock create_workspace must use the same extent.
inline int flexaid_one_based_copy_n(int count) noexcept {
    return count + 1;
}
