// config_parser.cpp — JSON config loader & applier for FlexAIDdS
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#include "config_parser.h"
#include "config_defaults.h"
#include "flexaid.h"
#include "gaboom.h"
#include "memetic_gate.h"
#include "statmech.h"

#include <cstring>
#include <stdexcept>

// ─── Helper: safely get a value from section.key with fallback ───────────

static bool jbool(const json::Value& cfg, const char* section, const char* key, bool fb) {
    return cfg[section][key].is_null() ? fb : cfg[section][key].as_bool(fb);
}
static int jint(const json::Value& cfg, const char* section, const char* key, int fb) {
    return cfg[section][key].is_null() ? fb : cfg[section][key].as_int(fb);
}
static double jdbl(const json::Value& cfg, const char* section, const char* key, double fb) {
    return cfg[section][key].is_null() ? fb : cfg[section][key].as_double(fb);
}
static float jflt(const json::Value& cfg, const char* section, const char* key, float fb) {
    return cfg[section][key].is_null() ? fb : cfg[section][key].as_float(fb);
}
static std::string jstr(const json::Value& cfg, const char* section, const char* key, const char* fb) {
    return cfg[section][key].is_null() ? std::string(fb) : cfg[section][key].as_string(fb);
}

// ─── Public API ──────────────────────────────────────────────────────────

json::Value load_config(const std::string& config_path) {
    json::Value defaults = flexaid_default_config();

    if (config_path.empty())
        return defaults;

    json::Value user_config = json::parse_file(config_path);
    return json::merge(defaults, user_config);
}

void apply_config(const json::Value& config, FA_Global* FA, GB_Global* GB,
                  const flexaids::ProtocolConfig* protocol) {
    // Single snapshot at entry: avoids dual-protocol if env mutates mid-apply.
    const flexaids::ProtocolConfig local_proto =
        protocol ? flexaids::ProtocolConfig{} : flexaids::ProtocolConfig::from_env();
    const flexaids::ProtocolConfig& proto = protocol ? *protocol : local_proto;

    // ── Scoring ──
    {
        auto complf = jstr(config, "scoring", "function", "VCT");
        std::strncpy(FA->complf, complf.c_str(), sizeof(FA->complf) - 1);

        auto sc = jstr(config, "scoring", "self_consistency", "MAX");
        std::strncpy(FA->vcontacts_self_consistency, sc.c_str(), sizeof(FA->vcontacts_self_consistency) - 1);

        auto pd = jstr(config, "scoring", "plane_definition", "X");
        FA->vcontacts_planedef = pd.empty() ? 'X' : pd[0];

        FA->normalize_area = jbool(config, "scoring", "normalize_area", false) ? 1 : 0;
        // P9: VCT distance-weighted contacts — exp(-r/r0) per-contact decay.
        // Default 7.0 Å (ON); set <= 0 to restore legacy equal-weight contacts.
        FA->vct_dist_weight_r0 = jdbl(config, "scoring", "vct_dist_weight_r0", 7.0);
        // Lever 2: intensive CF.com (divide by contact count). Default OFF
        // keeps the extensive score so existing arms stay byte-for-byte stable.
        FA->vct_normalize_contacts = jbool(config, "scoring", "vct_normalize_contacts", false) ? 1 : 0;
        FA->vct_entropy_weight = jdbl(config, "scoring", "vct_entropy_weight", 0.0);
        if (proto.vct_entropy_weight_set)
            FA->vct_entropy_weight = proto.vct_entropy_weight;
        FA->useacs          = jbool(config, "scoring", "accessible_surface", false) ? 1 : 0;
        FA->acsweight       = jflt(config, "scoring", "acs_weight", 1.0f);
        FA->solventterm     = jflt(config, "scoring", "solvent_penalty", 0.0f);
        FA->sas_weight      = jdbl(config, "scoring", "sas_weight", 1.0);
        FA->pb_clash_weight   = jdbl(config, "scoring", "pb_clash_weight", 0.0);
        FA->pb_clash_exponent = jdbl(config, "scoring", "pb_clash_exponent", 3.0);
        FA->pb_clash_ratio    = jdbl(config, "scoring", "pb_clash_ratio", 0.75);
        FA->pb_pocket_weight  = jdbl(config, "scoring", "pb_pocket_weight", 0.0);
        FA->pb_pocket_radius  = jdbl(config, "scoring", "pb_pocket_radius", 6.0);
        // Env overrides, applied after the JSON read and in the same order and with
        // the same validity rules as the legacy .inp path (top.cpp). Without these
        // the JSON path silently clobbered the env values with its own defaults, so
        // exporting FLEXAIDDS_PB_CLASH_WEIGHT into a DatasetRunner campaign — the
        // ONLY documented way to enable pb_clash, since the generated config omits
        // the key — was a no-op. Unset env leaves the JSON/default value untouched.
        if (const char* e = std::getenv("FLEXAIDDS_PB_CLASH_WEIGHT"))  { FA->pb_clash_weight = atof(e); }
        if (const char* e = std::getenv("FLEXAIDDS_PB_CLASH_EXP"))     { double v = atof(e); if (v > 0.0) FA->pb_clash_exponent = v; }
        if (const char* e = std::getenv("FLEXAIDDS_PB_CLASH_RATIO"))   { double v = atof(e); if (v > 0.0) FA->pb_clash_ratio = v; }
        if (const char* e = std::getenv("FLEXAIDDS_PB_POCKET_WEIGHT")) { FA->pb_pocket_weight = atof(e); }
        if (const char* e = std::getenv("FLEXAIDDS_PB_POCKET_RADIUS")) { double v = atof(e); if (v > 0.0) FA->pb_pocket_radius = v; }

        // Angular-dependent hydrogen bond potential
        const bool hbond_on = jbool(config, "scoring", "hbond_enabled", false);
        FA->use_hbond              = hbond_on ? 1 : 0;
        // v58 split: search (GA fitness) vs rank (post-GA re-score).  When the
        // split keys are absent, both follow hbond_enabled for legacy configs.
        FA->use_hbond_search       = jbool(config, "scoring", "hbond_search_enabled", hbond_on) ? 1 : 0;
        FA->use_hbond_rank         = jbool(config, "scoring", "hbond_rank_enabled", hbond_on) ? 1 : 0;
        FA->hbond_optimal_dist     = jdbl(config, "scoring", "hbond_optimal_distance", 2.8);
        FA->hbond_optimal_angle    = jdbl(config, "scoring", "hbond_optimal_angle", 180.0);
        FA->hbond_sigma_dist       = jdbl(config, "scoring", "hbond_sigma_distance", 0.4);
        FA->hbond_sigma_angle      = jdbl(config, "scoring", "hbond_sigma_angle", 30.0);
        FA->hbond_weight           = jdbl(config, "scoring", "hbond_weight", -2.5);
        FA->hbond_salt_bridge_weight = jdbl(config, "scoring", "hbond_salt_bridge_weight", -5.0);

        // Metal ion coordination potential
        FA->use_metal_coord      = jbool(config, "scoring", "metal_coord_enabled", false) ? 1 : 0;
        FA->metal_coord_weight   = jdbl(config, "scoring", "metal_coord_weight", 1.0);
        FA->metal_coord_sigma    = jdbl(config, "scoring", "metal_coord_sigma", 0.45);
        FA->metal_coord_cn_weight = jdbl(config, "scoring", "metal_coord_cn_weight", 0.5);

        // Coulomb electrostatics (E7 / Wave 1.2). Default OFF — USEELC legacy
        // keyword still enables; modern JSON was unreachable before this key.
        // Env FLEXAIDDS_USE_ELEC=1 forces ON; =0 forces OFF after JSON.
        FA->use_elec = jbool(config, "scoring", "electrostatics_enabled", false) ? 1 : 0;
        if (jbool(config, "scoring", "use_elec", false))
            FA->use_elec = 1;
        if (const char* e = std::getenv("FLEXAIDDS_USE_ELEC")) {
            FA->use_elec = (e[0] != '\0' && std::atoi(e) != 0) ? 1 : 0;
        }

        {
            double tw = jdbl(config, "scoring", "tencom_weight", 0.0);
            if (tw < 0.0) tw = 0.0;
            if (tw > 2.0) tw = 2.0;
            FA->tencom_weight = static_cast<float>(tw);
        }

        // GIST desolvation grid — HARD-DISABLED for all runs until the
        // evaluator / grid-type confusion is repaired (audit 2026-07-17).
        // Strict claims must not consume GIST scores; ignore JSON enable.
        const bool gist_requested = jbool(config, "scoring", "gist_enabled", false);
        FA->use_gist = 0;
        FA->gist_weight = 0.0;
        FA->gist_evaluator = nullptr;
        if (gist_requested) {
            fprintf(stderr,
                "WARN [GIST]: gist_enabled requested but HARD-DISABLED until "
                "evaluator/grid type confusion is repaired (strict claims)\n");
        }
        // Note: re-enable only behind a new validated gate + tests; never via
        // gist_enabled alone until that repair lands.
    }

    // ── Optimization ──
    {
        FA->delta_angstron = jdbl(config, "optimization", "translation_step", 0.25);
        FA->delta_angle    = jdbl(config, "optimization", "angle_step", 5.0);
        FA->delta_dihedral = jdbl(config, "optimization", "dihedral_step", 5.0);
        FA->delta_flexible = jdbl(config, "optimization", "flexible_step", 10.0);
        FA->spacer_length  = jflt(config, "optimization", "grid_spacing", 0.375f);
    }

    // ── Flexibility ──
    {
        FA->deelig_flex          = jbool(config, "flexibility", "ligand_torsions", true) ? 1 : 0;
        FA->intramolecular       = jbool(config, "flexibility", "intramolecular", true) ? 1 : 0;
        FA->intrafraction        = jflt(config, "flexibility", "intramolecular_fraction", 1.0f);
        FA->permeability         = jflt(config, "flexibility", "permeability", 1.0f);
        FA->rotamer_permeability = jflt(config, "flexibility", "rotamer_permeability", 0.8f);
        FA->soft_wall_cutoff     = jflt(config, "flexibility", "soft_wall_cutoff", 0.40f);
        FA->intermolecular_clash_ratio =
            jflt(config, "flexibility", "intermolecular_clash_ratio", 0.0f);
        FA->pbloops              = jint(config, "flexibility", "binding_site_conformations", 1);
        FA->bloops               = jint(config, "flexibility", "bonded_loops", 2);
        FA->useflexdee           = jbool(config, "flexibility", "use_flexdee", false) ? 1 : 0;
        FA->dee_clash            = jflt(config, "flexibility", "dee_clash", 0.5f);
    }

    // ── Thermodynamics ──
    {
        FA->temperature = static_cast<unsigned int>(jint(config, "thermodynamics", "temperature", 300));
        if (FA->temperature > 0) {
            // β = 1 / T — in FlexAID's partition function Z = Σ exp(−CF_i/T), T is an
            // effective softmax temperature over the (unitless) CF landscape, NOT a
            // physical kT in kcal/mol (see Morency, 3Dsig 2017). Folding in kB_kcal
            // made β ~503× larger, collapsing the Boltzmann weights to a single-pose
            // delta → entropy → 0 → ΔG reverts to raw CF ranking. Keep β = 1/T.
            FA->beta = 1.0 / static_cast<double>(FA->temperature);
        } else {
            FA->beta = 0.0;
        }
        FA->cluster_rmsd = jflt(config, "thermodynamics", "cluster_rmsd", 2.0f);
        FA->use_super_cluster = jbool(config, "thermodynamics", "use_super_cluster", false);
        // Classic FlexAID entropy ranking is the default product when T>0.
        // force_cf_rank_emission restores P3b lowest-CF emission (easy rollback).
        // classic_entropy_ranking=false is an alias that also forces CF emission.
        FA->force_cf_rank_emission = jbool(config, "thermodynamics", "force_cf_rank_emission", false);
        if (!jbool(config, "thermodynamics", "classic_entropy_ranking", true)) {
            FA->force_cf_rank_emission = true;
        }
        // ProtocolConfig overrides (snapshot); nullopt = leave JSON-derived value.
        if (proto.force_cf_rank_emission.has_value()) {
            FA->force_cf_rank_emission = *proto.force_cf_rank_emission;
        }
        if (proto.classic_entropy_ranking.has_value()) {
            // classic=false → CF emission; classic=true → entropy ranking (default).
            FA->force_cf_rank_emission = !*proto.classic_entropy_ranking;
        }
        FA->use_tqcm  = jbool(config, "turboquant", "compressed_contact_matrix", false);
        FA->use_tqens = jbool(config, "turboquant", "ensemble_compression", false);
        FA->use_tqnn  = jbool(config, "turboquant", "compressed_nn", false);

        auto ca = jstr(config, "thermodynamics", "clustering_algorithm", "CF");
        std::strncpy(FA->clustering_algorithm, ca.c_str(), sizeof(FA->clustering_algorithm) - 1);
    }

    // ── GA ──
    {
        GB->num_chrom        = jint(config, "ga", "num_chromosomes", 1000);
        GB->max_generations  = jint(config, "ga", "num_generations", 2000);  // P6: 500→2000 base budget
        GB->cross_rate       = jdbl(config, "ga", "crossover_rate", 0.8);
        GB->mut_rate         = jdbl(config, "ga", "mutation_rate", 0.03);

        auto fm = jstr(config, "ga", "fitness_model", "PSHARE");
        std::strncpy(GB->fitness_model, fm.c_str(), sizeof(GB->fitness_model) - 1);

        auto rm = jstr(config, "ga", "reproduction_model", "BOOM");
        std::strncpy(GB->rep_model, rm.c_str(), sizeof(GB->rep_model) - 1);

        GB->pbfrac = jdbl(config, "ga", "boom_fraction", 1.0);

        auto pi = jstr(config, "ga", "population_init", "RANDOM");
        std::strncpy(GB->pop_init_method, pi.c_str(), sizeof(GB->pop_init_method) - 1);

        GB->seed        = jint(config, "ga", "seed", 0);
        GB->adaptive_ga = jbool(config, "ga", "adaptive", true) ? 1 : 0;  // P6: ADAPTVGA default ON

        // adaptive_k array (P6 ADAPTVGA defaults: crossover 0.95 max, mutation 0.10
        // max, full crossover for below-avg individuals, 0.05 below-avg mutation).
        // Fall through to these when the config omits the array so adaptive Pc/Pm
        // never run with zero-initialised k1..k4 (which would disable crossover).
        const auto& ak = config["ga"]["adaptive_k"];
        if (ak.is_array() && ak.size() >= 4) {
            GB->k1 = ak[static_cast<size_t>(0)].as_double(0.95);
            GB->k2 = ak[static_cast<size_t>(1)].as_double(0.10);
            GB->k3 = ak[static_cast<size_t>(2)].as_double(1.0);
            GB->k4 = ak[static_cast<size_t>(3)].as_double(0.05);
        } else {
            GB->k1 = 0.95;
            GB->k2 = 0.10;
            GB->k3 = 1.0;
            GB->k4 = 0.05;
        }

        GB->alpha       = jdbl(config, "ga", "sharing_alpha", 4.0);  // P5: niche-sharing exponent (was 1.0)
        GB->peaks       = jdbl(config, "ga", "sharing_peaks", 5.0);
        GB->scale       = jdbl(config, "ga", "sharing_scale", 10.0);
        GB->intragenes  = jbool(config, "ga", "intragenes", false) ? 1 : 0;
        GB->duplicates  = jbool(config, "ga", "duplicates", false) ? 1 : 0;
        GB->ini_mut_prob = jdbl(config, "ga", "initial_mutation_prob", 0.0);
        GB->end_mut_prob = jdbl(config, "ga", "end_mutation_prob", 0.0);
        GB->ssnum            = jint(config, "ga", "steady_state_num", 0);
        GB->entropy_weight   = jdbl(config, "ga", "entropy_weight", 0.5);
        GB->entropy_interval = jint(config, "ga", "entropy_interval", 0);
        GB->use_shannon      = jbool(config, "ga", "use_shannon", false) ? 1 : 0;

        // Diversity monitoring (entropy collapse mitigation)
        GB->diversity_monitoring          = jbool(config, "ga", "diversity_monitoring", true) ? 1 : 0;  // P8: default ON (dual SEC gate)
        GB->diversity_check_interval      = jint(config, "ga", "diversity_check_interval", 10);
        GB->diversity_collapse_threshold  = jdbl(config, "ga", "diversity_collapse_threshold", 0.3);
        GB->catastrophic_mutation_fraction = jdbl(config, "ga", "catastrophic_mutation_fraction", 0.2);
        GB->catastrophic_mutation_count   = 0;

        // P5: periodic BOOM random injection (diversity insurance)
        GB->boom_inject_interval          = jint(config, "ga", "boom_inject_interval", 100);
        GB->boom_inject_fraction          = jdbl(config, "ga", "boom_inject_fraction", 1.0);
        GB->boom_inject_count             = 0;
        // Wave 3 / P1 anti-collapse: env overrides (default OFF = leave JSON).
        // FLEXAIDDS_BOOM_INTERVAL=<gens>  (0 disables periodic injection)
        // FLEXAIDDS_BOOM_FRAC=<0..1>
        // FLEXAIDDS_SIGMA_SCALE=<k> multiplies GB->scale (niche radius control)
        if (const char* e = std::getenv("FLEXAIDDS_BOOM_INTERVAL")) {
            int v = std::atoi(e);
            if (v >= 0) GB->boom_inject_interval = v;
        }
        if (const char* e = std::getenv("FLEXAIDDS_BOOM_FRAC")) {
            double v = std::atof(e);
            if (v >= 0.0 && v <= 1.0) GB->boom_inject_fraction = v;
        }
        if (const char* e = std::getenv("FLEXAIDDS_SIGMA_SCALE")) {
            double v = std::atof(e);
            if (v > 0.0) GB->scale *= v;
        }

        // ── True GA elitism (v27) ──
        GB->n_elite                       = jint(config, "ga", "n_elite", 1);

        // ProtocolConfig overrides for v27 GA-internal elitism knobs. These win
        // over JSON so a single binary can sweep them without re-emitting
        // dock_config.json (env adapter via ProtocolConfig::from_env).
        if (proto.n_elite_set)
            GB->n_elite = proto.n_elite;
        if (proto.sharing_alpha.has_value())
            GB->alpha = *proto.sharing_alpha;
        if (proto.boom_frac.has_value())
            GB->boom_inject_fraction = *proto.boom_frac;

        // Entropy-ablation hooks (unset → byte-identical to JSON defaults).
        if (proto.entropy_weight.has_value())
            GB->entropy_weight = *proto.entropy_weight;
        if (proto.diversity_monitoring.has_value())
            GB->diversity_monitoring = *proto.diversity_monitoring ? 1 : 0;
    }

    // ── Output ──
    {
        FA->max_results        = jint(config, "output", "max_results", 10);
        FA->output_scored_only = jbool(config, "output", "scored_only", false) ? 1 : 0;
        FA->score_ligand_only  = jbool(config, "output", "score_ligand_only", false) ? 1 : 0;
        FA->htpmode            = jbool(config, "output", "htp_mode", false);
        GB->num_print          = jint(config, "output", "print_chromosomes", 10);
        GB->print_int          = jint(config, "output", "print_interval", 1);
        GB->rrg_skip           = jint(config, "output", "rrg_skip", 0);
        GB->outgen             = jbool(config, "output", "output_generations", false) ? 1 : 0;
        FA->rotout             = jbool(config, "output", "rotamer_output", false) ? 1 : 0;
    }

    // ── Protein ──
    {
        FA->is_protein   = jbool(config, "protein", "is_protein", true) ? 1 : 0;
        FA->exclude_het  = jbool(config, "protein", "exclude_het", false) ? 1 : 0;
        FA->remove_water = jbool(config, "protein", "remove_water", true) ? 1 : 0;
        FA->keep_ions    = jbool(config, "protein", "keep_ions", true) ? 1 : 0;
        FA->keep_structural_waters =
            jbool(config, "protein", "keep_structural_waters", true) ? 1 : 0;
        FA->structural_water_bfactor_max =
            jflt(config, "protein", "structural_water_bfactor_max", 20.0f);
        FA->omit_buried  = jbool(config, "protein", "omit_buried", false) ? 1 : 0;
    }

    // ── Reference Ligand Seeding ──
    {
        auto rf = jstr(config, "reference_ligand", "file", "");
        std::strncpy(FA->reflig_file, rf.c_str(), sizeof(FA->reflig_file) - 1);
        FA->reflig_seed_fraction =
            jflt(config, "reference_ligand", "seed_fraction", 0.25f);
        FA->reflig_pose_seed_enabled =
            jbool(config, "reference_ligand", "pose_seed_enabled", true) ? 1 : 0;
        FA->reflig_k_nearest =
            jint(config, "reference_ligand", "k_nearest", 10);
        FA->reflig_hetatm_fallback =
            jbool(config, "reference_ligand", "hetatm_fallback", true) ? 1 : 0;
    }

    // ── Cavity-only MIF seeding ──
    {
        FA->mif_enabled =
            jbool(config, "seeding", "mif_enabled", false) ? 1 : 0;
        FA->mif_temperature =
            jflt(config, "seeding", "mif_temperature", 300.0f);
        FA->grid_prio_percent =
            jflt(config, "seeding", "grid_prio_percent", 100.0f);
    }

    // ── Advanced ──
    {
        FA->vindex             = jbool(config, "advanced", "vcontacts_index", false) ? 1 : 0;
        FA->supernode          = jbool(config, "advanced", "supernode", false) ? 1 : 0;
        FA->force_interaction  = jbool(config, "advanced", "force_interaction", false) ? 1 : 0;
        FA->interaction_factor = jflt(config, "advanced", "interaction_factor", 5.0f);
        FA->assume_folded      = jbool(config, "advanced", "assume_folded", false) ? 1 : 0;
    }

    // ── Coarse-init pocket scan ──
    {
        FA->coarse_init_enabled   = jbool(config, "coarse_init", "enabled",       false);
        FA->coarse_init_grid_step = jflt (config, "coarse_init", "grid_step",     3.0f);
        FA->coarse_init_n_seeds   = jint (config, "coarse_init", "n_seeds",       25);
        FA->coarse_init_n_orient  = jint (config, "coarse_init", "n_orientations",64);
        // Wave 3: FLEXAIDDS_COARSE_ORIENTATIONS overrides n_orientations (pilot A/B).
        if (const char* e = std::getenv("FLEXAIDDS_COARSE_ORIENTATIONS")) {
            int v = std::atoi(e);
            if (v > 0 && v <= 4096) FA->coarse_init_n_orient = v;
        }
        // coarse_seeds_* arrays are populated at runtime by run_coarse_pocket_scan()
        FA->coarse_seeds_grid     = nullptr;
        FA->coarse_seeds_genes    = nullptr;
        FA->coarse_seeds_count    = 0;
    }

    // Wave 3.4 memetic: real enable flag FA->use_memetic (default 0).
    // Requires FLEXAIDDS_MEMETIC=1 AND a pilot unlock:
    //   preferred: FLEXAIDDS_PB_CLASH_PHASE2_PASS=1 (ROADMAP_v2 SCORING-LOCKED
    //              magnitude-floor PASS), or
    //   legacy:    FLEXAIDDS_WALL_PILOT_PASS=1 (structurally unpassable; do not set).
    // MEMETIC alone cannot arm the feature. Logic in memetic_gate.h (unit-tested).
    {
        const char* e = std::getenv("FLEXAIDDS_MEMETIC");
        const bool want = e != nullptr && e[0] != '\0' && std::atoi(e) != 0;
        FA->use_memetic = flexaids::resolve_use_memetic_from_env();
        if (FA->use_memetic) {
            fprintf(stderr,
                "[MEMETIC] enabled (FLEXAIDDS_MEMETIC=1 and "
                "PB_CLASH_PHASE2_PASS or WALL_PILOT_PASS) use_memetic=%d\n",
                FA->use_memetic);
        } else if (want) {
            fprintf(stderr,
                "WARN [MEMETIC]: FLEXAIDDS_MEMETIC=1 ignored — set "
                "FLEXAIDDS_PB_CLASH_PHASE2_PASS=1 only after revised Phase 2 "
                "SCORING-LOCKED pb_clash oracle PASS (magnitude floor ≥1.0 kcal "
                "+ sign flip ≥2/3; ROADMAP_v2); use_memetic=0\n");
        }
    }

    // ── ThermodynamicEngine ──
    {
        FA->thermo_engine_enabled = jbool(config, "thermo_engine", "enabled",      false);
        FA->thermo_T_eff          = jflt (config, "thermo_engine", "T_eff",        0.596f);
        FA->thermo_tencom_scale   = jflt (config, "thermo_engine", "tencom_scale", 1.0f);
        // Reporting-only T for whiteboard diagnostics (I_ES/regime/CF_r2s); does
        // NOT feed thermo_T_eff/G_bind. ISMB 2017 calibrated default = 21.0
        // (kT_ISMB) — the LEFT-hand-defining constant behind ΔG₂₁/P_i(T=21).
        FA->thermo_report_T       = jflt (config, "thermo_engine", "report_T",     21.0f);
        FA->H_rep_bound_complex   = 0.0f;
        FA->H_rep_receptor_ref    = 0.0f;
        FA->H_rep_ligand_ref      = 0.0f;
        FA->thermo_result         = {};
        FA->thermo_engine         = nullptr;

        if (FA->thermo_engine_enabled) {
            FA->thermo_engine = new ThermodynamicEngine(FA->thermo_T_eff, FA->thermo_tencom_scale);
            FA->thermo_engine->set_unbound_reference(FA->H_rep_receptor_ref, FA->H_rep_ligand_ref);
        }
    }

    // Always GA
    std::strcpy(FA->metopt, "GA");
}
