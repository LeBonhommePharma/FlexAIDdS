// coarse_init.h — Coarse pocket-scan seeding for gen-0 (no-seed docking)
//
// Scores a coarse grid over the binding cleft with random orientations before
// the GA initialisation loop, then injects the lowest apparent-CF
// placements as guaranteed gen-0 seeds so the population starts with real
// pocket contacts rather than floating or sterically invalid chromosomes.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma / NRGlab, Université de Montréal

#pragma once

#include <functional>
#include <cstdint>

#include "flexaid.h"
#include "gaboom.h"

// Run the coarse pocket scan and populate FA->coarse_seeds_grid / coarse_seeds_genes.
//
// Must be called after reflig_nearest_grid has been populated (top.cpp direct-mode
// seeding block) and after gene_lim has been initialised (set_gene_lim in gaboom.cpp).
// The caller is populate_chromosomes() in gaboom.cpp — all required data is available.
//
// On return FA->coarse_seeds_count contains the number of seeds ready for injection
// (≤ FA->coarse_init_n_seeds). A zero count is safe — the caller falls through to
// the normal RANDOM gen-0 loop unchanged.
void run_coarse_pocket_scan(
    FA_Global*                  FA,
    VC_Global*                  VC,
    GB_Global*                  GB,
    atom*                       atoms,
    resid*                      residue,
    gridpoint*                  cleftgrid,
    const genlim*               gene_lim,
    std::function<int32_t()>&   dice);
