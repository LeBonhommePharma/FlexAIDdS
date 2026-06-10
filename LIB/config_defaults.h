// config_defaults.h — Single source of truth for all FlexAIDdS parameters
//
// Every key has a sensible default — the user only needs to override
// what they want to change via a JSON config file (-c flag).
//
// Design principle: FULL FLEXIBILITY ON by default.
// Use --rigid to disable all flexibility for fast screening.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "json_value.h"

inline json::Value flexaid_default_config() {
    using V = json::Value;
    using O = json::Object;
    using A = json::Array;

    return V(O{
        // ── Scoring ──────────────────────────────────────────────
        {"scoring", V(O{
            {"function",         V("VCT")},     // Voronoi contact scoring (best accuracy)
            {"self_consistency", V("MAX")},      // A→B and B→A contact self-consistency
            {"plane_definition", V("X")},        // plane definition mode
            {"normalize_area",   V(false)},
            {"vct_dist_weight_r0", V(4.0)},      // P9: per-contact exp(-r/r0) decay length (Å); <=0 disables
            {"accessible_surface", V(false)},    // ACS weighting
            {"acs_weight",       V(1.0)},
            {"solvent_penalty",  V(0.0)},
            // Angular-dependent hydrogen bond potential
            {"hbond_enabled",           V(false)},
            {"hbond_optimal_distance",  V(2.8)},
            {"hbond_optimal_angle",     V(180.0)},
            {"hbond_sigma_distance",    V(0.4)},
            {"hbond_sigma_angle",       V(30.0)},
            {"hbond_weight",            V(-2.5)},
            {"hbond_salt_bridge_weight",V(-5.0)},
            // Metal ion coordination potential
            {"metal_coord_enabled",     V(false)},
            {"metal_coord_weight",      V(1.0)},
            {"metal_coord_sigma",       V(0.45)},
            {"metal_coord_cn_weight",   V(0.5)},
            // GIST desolvation grid
            {"gist_enabled",  V(false)},
            {"gist_dx_file",  V("")},
            {"gist_weight",   V(1.0)},
        })},

        // ── Optimization step sizes ──────────────────────────────
        {"optimization", V(O{
            {"translation_step", V(0.25)},       // Angstroms
            {"angle_step",       V(5.0)},        // degrees
            {"dihedral_step",    V(5.0)},        // degrees
            {"flexible_step",    V(10.0)},       // degrees for ligand flex bonds
            {"grid_spacing",     V(0.375)},      // grid spacer length
        })},

        // ── Flexibility (FULL by default) ────────────────────────
        {"flexibility", V(O{
            {"ligand_torsions",             V(true)},   // DEEFLX: dead-end elimination for ligand flex
            {"intramolecular",              V(true)},   // consider intramolecular forces
            {"intramolecular_fraction",     V(1.0)},
            {"permeability",                V(1.0)},    // atom permeability
            {"rotamer_permeability",        V(0.8)},    // rotamer acceptance VDW permeability
            {"ring_conformers",             V(true)},   // LigandRingFlex sampling
            {"chirality",                   V(true)},   // ChiralCenter R/S discrimination
            {"binding_site_conformations",  V(1)},      // pbloops
            {"bonded_loops",                V(2)},      // bloops: exclude interactions n bonds away
            {"use_flexdee",                 V(false)},  // dead-end elimination for sidechains
            {"dee_clash",                   V(0.5)},
            {"multi_model",                 V(false)},  // CCBM: multi-conformer receptor docking
        })},

        // ── Thermodynamics ───────────────────────────────────────
        {"thermodynamics", V(O{
            {"temperature",           V(300)},   // Kelvin (0 = entropy off)
            {"clustering_algorithm",  V("CF")},  // CF, DP, or FO
            {"cluster_rmsd",          V(2.0)},   // RMSD threshold for clustering
        })},

        // ── Genetic Algorithm ────────────────────────────────────
        {"ga", V(O{
            {"num_chromosomes",      V(1000)},
            {"num_generations",      V(2000)},  // P6: 500→2000 base budget (scaled by ceil(n_genes/4) in benchmark)
            {"crossover_rate",       V(0.8)},
            {"mutation_rate",        V(0.03)},
            {"fitness_model",        V("SMFREE")},  // SMFREE = entropy-aware (StatMech Free energy + sharing)
            {"reproduction_model",   V("BOOM")},
            {"boom_fraction",        V(1.0)},
            {"population_init",      V("RANDOM")},
            {"seed",                 V(0)},          // 0 = time-based
            {"adaptive",             V(true)},   // P6: ADAPTVGA on (Srinivas-Patnaik adaptive Pc/Pm)
            // k1=max crossover (0.95 at avg fitness, →0 at max), k2=max mutation,
            // k3=below-avg crossover (full disruption), k4=below-avg mutation.
            {"adaptive_k",           V(A{V(0.95), V(0.10), V(1.0), V(0.05)})},
            {"sharing_alpha",        V(4.0)},   // P5: niche-sharing exponent (was 1.0; 4 = original FlexAID)
            {"sharing_peaks",        V(5.0)},
            {"sharing_scale",        V(10.0)},
            {"intragenes",           V(false)},
            {"duplicates",           V(false)},
            {"initial_mutation_prob", V(0.0)},
            {"end_mutation_prob",    V(0.0)},
            {"steady_state_num",     V(0)},
            {"entropy_weight",       V(0.5)},   // SMFREE blending: 0=rank-only, 1=pure Boltzmann
            {"entropy_interval",     V(0)},     // log ensemble thermo every N gens (0=auto)
            {"use_shannon",          V(false)},  // include Shannon configurational entropy
            // Diversity monitoring (entropy collapse mitigation)
            // P8: ON by default — drives both catastrophic-mutation rescue and the
            // dual SEC termination gate (energy SEC only terminates when gene-space
            // allele entropy has also collapsed; see gaboom.cpp sec_may_terminate).
            {"diversity_monitoring",             V(true)},
            {"diversity_check_interval",         V(10)},
            {"diversity_collapse_threshold",     V(0.3)},
            {"catastrophic_mutation_fraction",   V(0.2)},
            // P5: periodic BOOM random injection (diversity insurance)
            {"boom_inject_interval",             V(100)},
            {"boom_inject_fraction",             V(1.0)},
        })},

        // ── Distributed Computing ───────────────────────────────
        {"distributed", V(O{
            {"backend", V("thread")},  // "thread" (default), "mpi"
        })},

        // ── Output ───────────────────────────────────────────────
        {"output", V(O{
            {"max_results",        V(10)},
            {"scored_only",        V(false)},
            {"score_ligand_only",  V(false)},
            {"htp_mode",           V(false)},
            {"print_chromosomes",  V(10)},
            {"print_interval",     V(1)},
            {"rrg_skip",           V(0)},
            {"output_generations", V(false)},
            {"output_range",       V(false)},
            {"rotamer_output",     V(false)},
        })},

        // ── Protein ──────────────────────────────────────────────
        {"protein", V(O{
            {"is_protein",                  V(true)},
            {"exclude_het",                 V(false)},
            {"remove_water",                V(true)},
            {"keep_ions",                   V(true)},
            {"keep_structural_waters",      V(true)},
            {"structural_water_bfactor_max",V(20.0f)},
            {"omit_buried",                 V(false)},
        })},

        // ── Reference Ligand Seeding ─────────────────────────────
        {"reference_ligand", V(O{
            {"file",             V("")},
            {"seed_fraction",    V(0.25)},
            {"k_nearest",        V(10)},
            {"hetatm_fallback",  V(true)},
        })},

        // ── Advanced ─────────────────────────────────────────────
        {"advanced", V(O{
            {"vcontacts_index",    V(false)},
            {"supernode",          V(false)},
            {"force_interaction",  V(false)},
            {"interaction_factor", V(5.0)},
            {"assume_folded",      V(false)},  // skip NATURaL co-translational chain growth
        })},
    });
}

inline json::Value flexaid_rigid_overrides() {
    using V = json::Value;
    using O = json::Object;
    return V(O{
        {"flexibility", V(O{
            {"ligand_torsions",      V(false)},
            {"intramolecular",       V(false)},
            {"ring_conformers",      V(false)},
            {"chirality",            V(false)},
            {"use_flexdee",          V(false)},
            {"permeability",         V(1.0)},
            {"rotamer_permeability", V(1.0)},
        })},
        {"thermodynamics", V(O{
            {"temperature", V(0)},
        })},
    });
}
