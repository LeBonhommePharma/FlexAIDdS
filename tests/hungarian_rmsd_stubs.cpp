// Slim stubs for test_hungarian_rmsd_bounds.
// Links the real LIB/calc_rmsd.cpp; only stub symbols that calc_rmsd.cpp
// references outside the Hungarian path under test.

#include "flexaid.h"

#include <cstdio>
#include <cstdlib>

void alter_mode(atom*, resid*, float*, int, int) {}
void buildcc(FA_Global*, atom*, int, int[]) {}

void Terminate(int status) {
    std::fprintf(stderr, "Terminate(%d) called from hungarian_rmsd test\n", status);
    std::exit(status == 0 ? 1 : status);
}
