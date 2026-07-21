// cmaes_mock_seams_stub.cpp — five engine-seam stubs for mock-only CMA-ES tests
// Apache-2.0 © 2026 Le Bonhomme Pharma
//
// cmaes_search.cpp always defines cmaes_run_dock, which references the real
// engine seams. Mock unit tests never call cmaes_run_dock, but the linker still
// needs definitions. Prefer these lightweight stubs over linking flexaid_core.
//
// If LIB/cmaes_search.cpp is compiled with -DFLEXAIDS_CMAES_MOCK_ONLY (future
// chunk1 enhancement that #if-outs the dock path), this TU may be omitted.

#include "gaboom.h"
#include "flexaid.h"

void set_gene_lim(FA_Global* /*FA*/, GB_Global* /*GB*/, genlim* /*gene_lim*/) {}

void set_bins(genlim* /*gene_lim*/, int /*num_genes*/) {}

// One-arg overload also declared in gaboom.h (legacy).
void set_bins(genlim* /*gene_lim*/) {}

cfstr eval_chromosome(FA_Global* /*FA*/, GB_Global* /*GB*/, VC_Global* /*VC*/,
                      const genlim* /*gene_lim*/, atom* /*atoms*/, resid* /*residue*/,
                      gridpoint* /*cleftgrid*/, gene* /*john*/,
                      cfstr (* /*function*/)(FA_Global*, VC_Global*, atom*, resid*,
                                             gridpoint*, int, double*)) {
    return cfstr{};
}

double get_cf_evalue(cfstr* cf, FA_Global* /*FA*/) {
    if (!cf) return 0.0;
    return cf->com + cf->wal + cf->sas + cf->con + cf->elec;
}

double get_apparent_cf_evalue(cfstr* cf) {
    if (!cf) return 0.0;
    return cf->com + cf->wal + cf->sas + cf->elec;
}
