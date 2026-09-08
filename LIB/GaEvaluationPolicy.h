#pragma once

#include "EnvFlags.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace flexaids {

// Shared by generation zero and every later CPU fitness evaluation. Limit
// allocation as well as scheduling, so deterministic runs need one workspace.
inline int ga_evaluation_threads() {
#if defined(FLEXAID_DETERMINISTIC) || !defined(_OPENMP)
    return 1;
#else
    return env_bool("FLEXAID_DETERMINISTIC", false) ? 1 : omp_get_max_threads();
#endif
}

}  // namespace flexaids
