// RibosomeElongation.h — Ribosome-speed co-translational elongation model
//
// Implements the master equation framework of Zhao et al. (2011)
// J. Phys. Chem. B 115, 3987–3997 (bibliographic citation only;
// no DOI was invented — doi-check: documented).
// 10.1021/jp109255g is Werner et al. J. Phys. Chem. C 2010, not Zhao.
//
//   dP_n(t)/dt = k_{n-1} · P_{n-1}(t) − k_n · P_n(t)
//
// where P_n(t) = probability that the ribosome has just incorporated residue n
// at time t, and k_n = codon-specific elongation rate (s⁻¹).
//
// Key physical parameters:
//   – Ribosomal exit tunnel: ~30–35 aa (Goldman et al. 2010, Cell 143:92)
//   – Mean elongation rate: 10–20 aa/s in vivo (Wohlgemuth et al. 2008)
//   – Codon-dependent rates: proportional to cognate tRNA abundance
//     (Dong et al. 1996; Ikemura 1985 calibration for E. coli / Human)
//   – Co-translational folding competition:
//       P_fold(n) = k_fold(n) / (k_fold(n) + k_el(n))
//
// Integration with NATURaL DualAssembly:
//   – Growth steps are time-stamped (t_n = Σ_{i≤n} 1/k_i)
//   – Slower codons (pausing) are flagged as co-translational folding windows
//   – Shannon entropy and CF contributions are time-weighted
//
// Reference implementation of the Zhao 2011 master equation is validated
// against the analytic mean elongation time T = 1/k_ini + Σ_n 1/k_n.
#pragma once

#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <functional>

namespace ribosome {

// ─── physical constants ───────────────────────────────────────────────────────
inline constexpr double TUNNEL_LENGTH_AA   = 34.0;  // residues occluded in tunnel
inline constexpr double MEAN_EL_RATE_ECOLI = 16.5;  // aa/s (Wohlgemuth 2008)
inline constexpr double MEAN_EL_RATE_HUMAN = 5.6;   // aa/s (Ingolia 2011)
inline constexpr double K_TERM_DEFAULT     = 4.0;   // termination rate (s⁻¹)
inline constexpr double K_INI_DEFAULT      = 0.1;   // initiation rate (s⁻¹)
inline constexpr double K_FOLD_DEFAULT     = 1.0;   // baseline folding rate s⁻¹ (protein / Mg-dep. tertiary)

// ─── RNA-specific folding rates ───────────────────────────────────────────────
// RNA hairpin/stem (secondary structure) folds in microseconds (k ~ 10³–10⁶ s⁻¹).
// References: Porschke 1974; Ninio 1987; Liphardt 2001; Woodside 2006 PNAS.
inline constexpr double K_FOLD_RNA_SECONDARY = 1.0e4; // s⁻¹ RNA hairpin/stem formation
// Mg²⁺-dependent tertiary folding (active-site conformation) — slow, at K_d = 1 mM.
// Reference: Penedo 2004 RNA; Martick & Scott 2006 Cell 126:309.
inline constexpr double K_FOLD_RNA_TERTIARY  = 1.0;   // s⁻¹ at Mg²⁺ = K_d (1 mM)
inline constexpr double KD_MG_RNA_MM         = 1.0;   // mM  Mg²⁺ K_d for tertiary folding
inline constexpr double N_HILL_MG            = 2.0;   // Hill coefficient (cooperative Mg²⁺)

// ─── RNA polymerase constants (co-transcriptional folding) ────────────────────
// Same Zhao 2011 master equation applied to nucleotide-by-nucleotide synthesis.
// Reference: Uptain 1997 Annu. Rev. Biochem. 66:117; Jonkers 2014 Mol. Cell 54:591.
inline constexpr double RNAP_TUNNEL_NT       = 8.0;   // nt in RNAP RNA:DNA hybrid
// Corrected from 50 nt/s: genome-wide mRNA in vivo rate is 21–25 nt/s at 37°C.
// Reference: Sci. Reports 2017 (doi:10.1038/s41598-017-17408-9); rRNA is 80–90 nt/s.
inline constexpr double MEAN_NT_RATE_ECOLI   = 25.0;  // nt/s mRNA in vivo (Sci. Reports 2017)
inline constexpr double MEAN_NT_RATE_HUMAN   = 25.0;  // nt/s (Jonkers 2014)
inline constexpr double K_RNAP_INI_DEFAULT   = 0.05;  // initiation rate (s⁻¹)
inline constexpr double K_RNAP_TERM_DEFAULT  = 5.0;   // termination rate (s⁻¹)
// Pause threshold for ribosome: rate < mean × 0.30 → co-translational folding window
// (Pechmann & Frydman 2013 Nat. Struct. Mol. Biol.)
inline constexpr double RIBOSOME_PAUSE_THRESHOLD = 0.30;
// Pause threshold for RNAP: rate < mean × 0.15 → co-transcriptional folding window
// (RNAP pauses are more extreme than ribosome pauses; Neuman 2003 Science 298:1152)
inline constexpr double RNAP_PAUSE_THRESHOLD = 0.15;

// ─── organism context ─────────────────────────────────────────────────────────
enum class Organism { EcoliK12, HumanHEK293 };

// ─── burst elongation unit ────────────────────────────────────────────────────
// A burst unit is a maximal run of consecutive fast codons (k_el > threshold×mean)
// after tunnel exit. Multiple residues in a burst elongate quasi-simultaneously
// before co-translational folding can compete (dwell_i << 1/k_fold).
struct BurstUnit {
    int    start_residue;    // 0-based first residue of burst
    int    end_residue;      // 0-based last residue (inclusive)
    int    n_residues;       // monomers effectively added "at once"
    double total_dwell;      // Σ 1/k_el[i] for i in burst (s)
    double mean_k;           // harmonic mean elongation rate in burst (s⁻¹)
    bool   follows_pause;    // true if immediately preceded by a pause site
};

// ─── codon rate table ─────────────────────────────────────────────────────────
// Maps 3-letter codon string → elongation rate (s⁻¹).
// Rates are calibrated from tRNA abundance × ribosome A-site occupancy times
// measured by ribosome profiling (Ingolia 2012; Li 2012).
// Default: E. coli K-12 codon usage (Dong 1996 tRNA gene copy numbers).
class CodonRateTable {
public:
    static CodonRateTable build_ecoli();
    static CodonRateTable build_human();

    // ── RNA polymerase rate tables (same master equation, nt-by-nt synthesis) ──
    // "Codons" are single nucleotides: "A", "U"/"T", "G", "C".
    // Rates reflect NTP incorporation kinetics and RNAP pause propensity.
    // Reference: Uptain 1997; Neuman 2003; Jonkers 2014.
    static CodonRateTable build_rnap_ecoli();   // mean = 50 nt/s
    static CodonRateTable build_rnap_human();   // mean = 25 nt/s

    // Elongation rate for a given codon (s⁻¹); returns mean rate if unknown
    double rate(const std::string& codon) const;

    // Mean elongation rate over a codon sequence
    double mean_rate(const std::vector<std::string>& codons) const;

    // Identify "pause sites": codons with rate < mean_rate × threshold
    std::vector<int> pause_sites(const std::vector<std::string>& codons,
                                  double threshold = RIBOSOME_PAUSE_THRESHOLD) const;

    // Most-Frequent Codon (MFC) mapping: amino-acid 1-letter code → dominant codon.
    // Derived from the highest-rate codon per AA in each organism's rate table.
    // Enables per-residue elongation rates for protein chains without explicit
    // DNA sequence (replaces uniform mean-rate fallback).
    static std::unordered_map<char, std::string> aa_to_mfc(Organism org);

    Organism organism;
    double   mean_rate_aa_per_s;

private:
    std::unordered_map<std::string, double> rates_; // codon → rate (s⁻¹)
    explicit CodonRateTable(Organism org, double mean_rate);
    void populate_ecoli();
    void populate_human();
    void populate_rnap_ecoli();  // single-nucleotide RNAP rates, E. coli
    void populate_rnap_human();  // single-nucleotide RNAP rates, Human
};

// ─── master equation state ────────────────────────────────────────────────────
// P_n(t) = probability density at ribosome position n (0 = initiation complex)
// Solved by explicit Euler integration with adaptive step size.
struct MasterEqState {
    int                 n_residues;
    double              t_current;           // s
    std::vector<double> P;                   // P[n] = probability at position n
    std::vector<double> t_arrival;           // leading-edge time (P[n] first > 1e-3) (s)
    std::vector<double> k_el;               // elongation rates (s⁻¹)
    double              k_ini;
    double              k_ter;
    double              t_median_terminal = -1.0; // time when P[N] ≥ 0.5 (median)

    // Fraction of ribosomes that have reached position n at time t_current
    double fraction_reached(int n) const;
};

// ─── RibosomeElongation ───────────────────────────────────────────────────────
class RibosomeElongation {
public:
    // Construct from amino-acid sequence (1-letter) + codon sequence.
    // If codons is empty, mean rates are used for every residue.
    RibosomeElongation(const std::string&              aa_sequence,
                       const std::vector<std::string>& codons,
                       const CodonRateTable&            rate_table,
                       double                           k_ini  = K_INI_DEFAULT,
                       double                           k_ter  = K_TERM_DEFAULT);

    // ── Analytical solution (Zhao 2011 Eq. 7) ────────────────────────────────
    // Mean first-passage time to incorporate residue n:
    //   <T_n> = 1/k_ini + Σ_{i=1}^{n} 1/k_i
    double mean_arrival_time(int n) const;

    // Total mean translation time for the full chain
    double mean_total_time() const;

    // ── Master equation integration ───────────────────────────────────────────
    // Integrate dP_n/dt = k_{n-1}*P_{n-1} − k_n*P_n from t=0 to t_max.
    // Returns time-resolved occupancy trajectories.
    MasterEqState integrate(double t_max,
                             double dt    = 1e-3,  // s (Euler step)
                             bool   adaptive = true) const;

    // ── Co-translational folding windows ─────────────────────────────────────
    struct FoldingWindow {
        int    residue_idx;    // first residue outside tunnel at this pause
        double t_available;    // time available for folding (s) = 1/k_el(n)
        double k_fold;         // estimated folding rate at this window (s⁻¹)
        double p_folded_cotrans; // k_fold / (k_fold + k_el)
        bool   is_pause_site;   // rate < mean × RIBOSOME_PAUSE_THRESHOLD
    };
    std::vector<FoldingWindow> folding_windows(double k_fold_base = K_FOLD_DEFAULT) const;

    // ── Time-weighted thermodynamic score ────────────────────────────────────
    // Integrates a per-residue score function S(n) weighted by dwell time 1/k_n.
    // Models the time each intermediate is sampled at the ribosome exit.
    double time_weighted_score(const std::function<double(int)>& score_fn) const;

    // ── Accessors ─────────────────────────────────────────────────────────────
    int    n_residues()   const noexcept { return static_cast<int>(k_el_.size()); }
    double k_ini()        const noexcept { return k_ini_; }
    double k_ter()        const noexcept { return k_ter_; }
    const  std::vector<double>& elongation_rates() const noexcept { return k_el_; }
    const  std::vector<int>&    pause_sites()      const noexcept { return pause_sites_; }

private:
    std::vector<double> k_el_;         // k_el_[n] = elongation rate for residue n
    double              k_ini_;
    double              k_ter_;
    std::vector<int>    pause_sites_;  // indices where k_el < RIBOSOME_PAUSE_THRESHOLD * mean
    double              mean_rate_;
};

// ─── validation ───────────────────────────────────────────────────────────────
// Validate the model against:
//   1. Analytic mean elongation time (should match integrate() steady-state)
//   2. Known E. coli rate: ~16.5 aa/s → T ≈ 60 aa / 16.5 ≈ 3.6 s for a 60-aa protein
// Returns true if validation passes within tolerance.
struct ValidationResult {
    bool   passed;
    double analytic_T;    // s (Σ 1/k_i + 1/k_ini)
    double ode_T;         // s (from integrate() half-max arrival)
    double relative_error;
    std::string message;
};
ValidationResult validate_master_equation(int n_residues = 60,
                                           Organism org = Organism::EcoliK12);

} // namespace ribosome
