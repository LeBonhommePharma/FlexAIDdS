// NATURaLDualAssembly.h — Native Assembly of co-Transcriptionally/co-Translationally
//                          Unified Receptor–Ligand (NATURaL) module
//
// Auto-detects when the ligand is nucleotide/sugar-containing or the receptor
// is a nucleic acid chain, then activates co-translational DualAssembly mode.
//
// In DualAssembly, the receptor chain grows residue-by-residue at ribosome speed
// (Zhao 2011 master equation; see RibosomeElongation.h) while the ligand is
// present from the start, computing incremental CF and Shannon entropy at each
// growth step to capture co-translational stereochemical selection.
//
// Also supports transmembrane (TM) domain insertion via the Sec translocon
// (see TransloconInsertion model below).
//
// Fully integrated with:
//   – RibosomeElongation (Zhao 2011 master equation, codon-dependent rates)
//   – ShannonThermoStack (entropy per growth step, time-weighted)
//   – SugarPuckerGene (activated for nucleotide ligands)
//   – AlphaShape Contact Function (incremental ΔSASA)
//   – TransloconInsertion (TM helix lateral gating, von Heijne 2007)
#pragma once

#include "../flexaid.h"
#include "../Vcontacts.h"
#include "../statmech.h"
#include "RibosomeElongation.h"
#include "TransloconInsertion.h"
#include "NucleationDetector.h"
#include <memory>
#include <string>
#include <vector>
#include <functional>

namespace natural {

// ─── detection ───────────────────────────────────────────────────────────────

// Returns true if any ligand atom has an O2' or ribose-pattern atom name,
// indicating a nucleotide or nucleoside ligand.
bool is_nucleotide_ligand(const atom* atoms, int n_lig_atoms);

// Returns true if the receptor contains nucleic acid residues
// (deoxyribose/ribose backbone: DA, DG, DC, DT, A, G, C, U, etc.)
bool is_nucleic_acid_receptor(const resid* residues, int n_residues);

// ─── NATURaL configuration ───────────────────────────────────────────────────
struct NATURaLConfig {
    bool             enabled                = false;
    bool             co_translational_growth = false;
    bool             sugar_pucker_auto       = true;
    double           temperature_K           = 298.15;
    int              max_growth_steps        = -1; // -1 = full sequence length
    ribosome::Organism organism              = ribosome::Organism::EcoliK12;
    bool             use_ribosome_speed      = true;  // use Zhao 2011 rates
    bool             model_tm_insertion      = true;  // model TM translocon

    // ── Ion-dependent RNA folding ──────────────────────────────────────────
    // For RNA receptors (ribozymes, riboswitches), Mg²⁺ drives tertiary folding.
    // Secondary structure (hairpins/stems) folds fast (~k_fold_rna_secondary).
    // Active-site/tertiary conformation is Mg²⁺-dependent (Hill equation).
    // References: Penedo 2004 RNA; Martick & Scott 2006 Cell 126:309.
    double           mg_concentration_mM    = 2.0;   // physiological [Mg²⁺] (mM)
    bool             ion_dependent_folding   = true;  // use Mg²⁺ Hill eq for RNA tertiary
    double           k_fold_rna_secondary    = ribosome::K_FOLD_RNA_SECONDARY; // s⁻¹ hairpin
    double           k_fold_rna_tertiary     = ribosome::K_FOLD_RNA_TERTIARY;  // s⁻¹ at K_d

    // ── Burst elongation detection ─────────────────────────────────────────
    // Consecutive fast codons produce bursts where no co-translational
    // folding can compete (dwell_i << 1/k_fold). k_el > burst_threshold×hmean
    // marks a residue as part of a burst unit.
    // Reference: Pechmann 2013 Nat Struct Mol Biol 20:237; Thanaraj 1996.
    double           burst_threshold        = 1.5;   // k_el/hmean ratio → burst member

    // ── Nucleation seed event detection ───────────────────────────────────
    // Sequence contexts that strongly promote co-translational structure
    // initiation. Detected seeds apply a folding_rate_boost multiplier on
    // k_fold at the matching positions in the growth loop.
    // RNA: Watson-Crick hairpins (Turner 2004 ΔG) + G-quadruplexes.
    // Protein: Chou-Fasman α-helix propensity + hydrophobic clusters.
    bool             detect_nucleation_seeds = true;
    int              min_hairpin_stem_bp     = 4;     // RNA: min stem length (bp)
    double           helix_propensity_thresh = 1.03;  // protein: Chou-Fasman P_α mean
    int              hydrophobic_run_min     = 4;     // protein: min ILVFMW run length
};

// Auto-configure from receptor and ligand properties.
NATURaLConfig auto_configure(const atom*  atoms,
                              int          n_lig_atoms,
                              const resid* residues,
                              int          n_residues);

// ─── Human in-vivo protofibril protocol planning ────────────────────────────
//
// A protofibril is best treated as the fixed docking scaffold when the nascent
// chain is the growing, high-entropy object. Reciprocal role swaps are still
// useful controls because FlexAID target/ligand machinery is not perfectly
// symmetric.
enum class GrowthProcess {
    Transcription,  // RNAP II, nucleotide-by-nucleotide
    Translation     // 80S ribosome, amino-acid-by-amino-acid
};

enum class DockingRolePolicy {
    ProtofibrilAsTarget,
    ReciprocalControl
};

struct InVivoAssemblyTrack {
    std::string       name;
    GrowthProcess     process = GrowthProcess::Translation;
    DockingRolePolicy role_policy = DockingRolePolicy::ProtofibrilAsTarget;
    int               chain_length = 0;       // nt for transcription, aa for translation
    double            mean_elongation_rate = 0.0; // nt/s or aa/s
    double            initiation_rate = 0.0;  // s^-1
    double            tunnel_length = 0.0;    // RNAP hybrid nt or ribosome tunnel aa
    bool              protofibril_is_target = true;
    bool              direct_encounter_allowed = true;
    std::string       compartment;

    double dwell_time_s() const noexcept;
    double completion_time_s() const noexcept;
    double first_exposed_time_s() const noexcept;
};

struct ParallelGrowthEvent {
    int           track_index = -1;
    int           unit_index = 0;       // 1-based nt/aa index
    double        t_arrival = 0.0;      // seconds since initiation gate
    int           exposed_units = 0;    // units outside polymerase/ribosome tunnel
    GrowthProcess process = GrowthProcess::Translation;
    DockingRolePolicy role_policy = DockingRolePolicy::ProtofibrilAsTarget;
};

std::vector<InVivoAssemblyTrack> make_human_protofibril_tracks(
    int  transcript_nt,
    int  peptide_aa,
    bool include_reciprocal_controls = true);

std::vector<ParallelGrowthEvent> build_parallel_growth_schedule(
    const std::vector<InVivoAssemblyTrack>& tracks);

// Shannon entropy of the Contact Function growth trajectory, returned in nats.
// This is intentionally testable because it feeds kT*H free-energy arithmetic.
double growth_entropy_nats(const std::vector<double>& cf_trajectory);

// ─── DualAssembly engine ─────────────────────────────────────────────────────
class DualAssemblyEngine {
public:
    explicit DualAssemblyEngine(const NATURaLConfig& cfg,
                                 FA_Global* FA, VC_Global* VC,
                                 atom* atoms, resid* residues,
                                 int n_residues);

    // Run the incremental co-translational growth simulation.
    // At each step: grow one residue, compute CF + Shannon entropy.
    // Returns the thermodynamic score trajectory.
    struct GrowthStep {
        int    residue_idx;
        double t_arrival;         // time since initiation (s) from master equation
        double k_el;              // elongation rate at this codon (s⁻¹)
        double dwell_time;        // 1/k_el (s) — time available for folding
        bool   is_pause_site;     // rate < RIBOSOME_PAUSE_THRESHOLD of mean → folding window
        bool   in_tunnel;         // residue still inside ribosomal exit tunnel
        double cf_score;          // Contact Function (kcal/mol)
        double shannon_entropy;   // dimensionless nats
        double p_cotrans_folded;  // k_fold/(k_fold+k_el) co-translational folding prob
        double cumulative_deltaG; // kcal/mol (time-weighted mean)
        bool   tm_inserted;       // true if TM helix has been laterally gated
        double tm_insertion_dG;   // ΔG of translocon-mediated insertion (kcal/mol)

        // ── Burst elongation annotation ──────────────────────────────────
        int    burst_unit_id;     // ≥0 index into burst table; -1 if not in a burst
        int    burst_size;        // n_residues in the burst unit (1 if not in burst)
        bool   follows_pause;     // true if this burst immediately follows a pause site

        // ── Nucleation seed annotation ───────────────────────────────────
        int    nucleation_seed_id;    // ≥0 index into seed list; -1 if no seed at this pos
        double nucleation_seed_boost; // k_fold multiplier applied (1.0 if no seed)
    };
    std::vector<GrowthStep> run();

    // Final ΔG after full co-translational folding
    double final_deltaG() const noexcept { return final_deltaG_; }

    // Whether NATURaL mode was active
    bool is_active() const noexcept { return config_.enabled; }

private:
    NATURaLConfig config_;
    FA_Global*    FA_;
    VC_Global*    VC_;
    atom*         atoms_;
    resid*        residues_;
    int           n_residues_;
    double        final_deltaG_ = 0.0;
    bool          is_rna_receptor_ = false;  // set in constructor via is_nucleic_acid_receptor()

    // Compute a lightweight CF score for the current partial complex
    double compute_partial_cf(int n_grown_residues) const;

    // Compute Shannon entropy over accumulated growth ensemble, returned in nats.
    double compute_growth_entropy(const std::vector<double>& cf_trajectory) const;

    // Hill equation: k_eff = k_max × [Mg]ⁿ / (K_d ⁿ + [Mg]ⁿ)
    // Used to scale RNA tertiary folding rate by Mg²⁺ availability.
    double mg_hill_factor() const noexcept;
};

} // namespace natural
