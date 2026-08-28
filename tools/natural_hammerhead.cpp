// tools/natural_hammerhead.cpp
// NATURaL co-transcriptional folding simulation of the minimal hammerhead ribozyme
//
// Uses the Zhao 2011 master equation framework (NATURaL module) in RNAP mode to model
// nucleotide-by-nucleotide synthesis of the Schistosoma mansoni minimal hammerhead
// ribozyme and compute co-transcriptional folding kinetics.
//
// Scientific references:
//   [1] Zhao et al. (2011) J. Phys. Chem. B 115:3987 — master equation framework
//   [2] Martick & Scott (2006) Science 313:1514    — HHRz crystal structure (PDB: 2OEU)
//   [3] Uptain et al. (1997) Annu. Rev. Biochem.  — E. coli RNAP rates (~50 nt/s)
//   [4] Neuman et al. (2003) Science 298:1152      — RNAP pausing
//   [5] Pechmann & Frydman (2013) Nat. Struct. Biol. — pause-site folding boost (3×)
//
// Build (via CMake):
//   cmake -DENABLE_NATURAL_HAMMERHEAD=ON -S . -B build && cmake --build build -t natural_hammerhead
//
// The hammerhead ribozyme (HHRz) is a self-cleaving RNA. Upon transcription by RNA
// polymerase, it folds co-transcriptionally into three stems (I, II, III) plus an
// invariant catalytic core. This simulation tracks the kinetics of each nucleotide
// incorporation event, identifies RNAP pause sites (co-transcriptional folding
// windows), and computes the probability that each stem element folds as it emerges
// from the RNAP RNA:DNA hybrid channel.

#include "NATURaL/RibosomeElongation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "TuiColor.h"

// ─── Terminal colours ─────────────────────────────────────────────────────────
// These were raw 8-colour ANSI behind a COMPILE-TIME `NO_COLOUR` guard and no
// isatty() anywhere, so `natural_hammerhead > run.log` wrote literal \033[3xm
// bytes into the log. They now resolve through TuiColor.h, which returns "" the
// moment stdout is not a terminal — redirects and pipes come out clean.
//
// The hues are the design system's, and each keeps ONE meaning:
//   tunnel / Stem III  -> aqua        (was blue)
//   pause site / Loop  -> strawberry  (was yellow — a retired colour)
//   folded / PASS      -> mint        (was green)
//   FAIL               -> fail text   (was red)
//   run complete       -> tangerine   (was cyan — also retired)
// The legend below names the colour it shows, so those words move too — a
// legend that says "Yellow" beside a strawberry swatch is just wrong.
#define COL_RESET  tui::reset()
#define COL_BOLD   tui::bold()
#define COL_RED    tui::failtext()
#define COL_GREEN  tui::mint()
#define COL_YELLOW tui::strawberry()
#define COL_BLUE   tui::aqua()
#define COL_CYAN   tui::tangerine()
#define COL_WHITE  tui::fg()

// ─── Hammerhead Ribozyme ────────────────────────────────────────────────────────
//
// Schistosoma mansoni HHRz (minimal, 40 nt; Martick & Scott 2006 Science 313:1514)
// This is the catalytic strand of the full-length minimal hammerhead ribozyme.
// The molecule folds into:
//
//   Stem I  (nt  1– 5 / 36–40):  GpGpCpUpA / UpApGpCpC  — 5 bp
//   Stem III (nt  6–12 / 23–29):  GpGpCpCpUpA / UpApGpGpCpC — 7 bp
//   Catalytic core (nt 13–22):    AUGAGUU UGAA  — 10 nt (invariant)
//   Stem II (nt 30–35):           CGAAAq  — GNRA-type tetraloop closing 6-bp stem
//   Cleavage site: C17 (5'-..CpGpA..-3', cleavage 3' of C17)
//
// Dot-bracket secondary structure (approximate):
//   (((((((..((((....)))).....((((....)))).....))))))):
//
static constexpr const char* HHRZ_SEQUENCE =
    // 1234567890123456789012345678901234567890
      "GGCUAUGGCCUAAUGAGUUUGAACCGAAACGUUCACUGCC";

static constexpr int HHRZ_LENGTH = 40;

// Structural annotation for each position (0-indexed)
// Region codes: S1=Stem-I, S2=Stem-II, S3=Stem-III, CR=Core, L2=Loop-II, L3=Loop-III
static constexpr std::array<const char*, 40> HHRZ_REGIONS = {
    "S1","S1","S1","S1","S1",          //  1- 5  Stem I (closing pair)
    "S3","S3","S3","S3","S3","S3","S3",//  6-12  Stem III
    "CR","CR","CR","CR","CR",          // 13-17  Catalytic core (cleavage at 17)
    "CR","CR","CR",                    // 18-20  Catalytic core continued
    "S3","S3","S3","S3","S3","S3","S3",// 21-27  Stem III (closing strand)
    "S2","S2","S2",                    // 28-30  Stem II
    "L2","L2","L2","L2",               // 31-34  Loop II (GAAA tetraloop)
    "S2","S2","S2",                    // 35-37  Stem II (closing strand)
    "S1","S1","S1"                     // 38-40  Stem I (outer closing strand)
};

// Brief structural role per region code
static const char* region_name(const char* code) {
    if (std::strcmp(code, "S1") == 0) return "Stem I";
    if (std::strcmp(code, "S2") == 0) return "Stem II";
    if (std::strcmp(code, "S3") == 0) return "Stem III";
    if (std::strcmp(code, "CR") == 0) return "Catalytic core";
    if (std::strcmp(code, "L2") == 0) return "Loop II";
    if (std::strcmp(code, "L3") == 0) return "Loop III";
    return "Unknown";
}

// ─── Shannon entropy over a real-valued trajectory ───────────────────────────
// Mirrors compute_growth_entropy() from NATURaLDualAssembly.cpp.
// Returns Shannon entropy (bits) over a 32-bin histogram of the trajectory.
static double shannon_entropy_bits(const std::vector<double>& traj) {
    if (traj.empty()) return 0.0;
    double mn = *std::min_element(traj.begin(), traj.end());
    double mx = *std::max_element(traj.begin(), traj.end());
    if (mx - mn < 1e-8) return 0.0;

    constexpr int BINS = 32;
    double bw = (mx - mn) / BINS + 1e-10;
    std::array<int, BINS> cnt{};
    for (double v : traj) {
        int b = static_cast<int>((v - mn) / bw);
        cnt[std::clamp(b, 0, BINS - 1)]++;
    }
    double H = 0.0, N = static_cast<double>(traj.size());
    const double log2_inv = 1.0 / std::log(2.0);
    for (int c : cnt) {
        if (c > 0) {
            double p = c / N;
            H -= p * std::log(p) * log2_inv;
        }
    }
    return H;
}

// ─── GrowthStep — mirrors NATURaLDualAssembly::DualAssemblyEngine::GrowthStep ─
struct GrowthStep {
    int    residue_idx;
    char   nucleotide;        // A / G / C / U
    const char* region;       // structural region code
    double t_arrival;         // s — mean first-passage time (Zhao 2011 Eq. 7)
    double k_el;              // nt/s  — RNAP rate at this position
    double dwell_time;        // 1/k_el (s) — time available for folding
    bool   is_pause_site;     // rate < 15% of mean (RNAP threshold)
    bool   in_tunnel;         // inside RNAP RNA:DNA hybrid (first ~8 nt)
    double p_cotrans_folded;  // k_fold / (k_fold + k_el)
    double shannon_entropy;   // bits over dwell-time trajectory so far
    double dwell_cumulative;  // cumulative dwell time up to this nt (s)
};

// ─── Run NATURaL RNAP simulation on a nucleotide sequence ─────────────────────
// Implements the DualAssemblyEngine::run() loop for RNAP / co-transcriptional mode,
// using RibosomeElongation directly. Equivalent to NATURaLDualAssembly with
//   cfg.use_ribosome_speed = false  (RNAP rates)
//   cfg.model_tm_insertion  = false (N/A for RNA)
//   FA_ = VC_ = nullptr            (CF scoring disabled; pure kinetics)
static std::vector<GrowthStep> run_rnap_natural(
    const std::string&              nt_sequence,
    const std::array<const char*, 40>& regions,
    ribosome::Organism              organism = ribosome::Organism::EcoliK12)
{
    using namespace ribosome;

    int N = static_cast<int>(nt_sequence.size());

    // ── Build RNAP rate table (same Zhao 2011 master equation, nt-by-nt) ──────
    CodonRateTable rnap_table =
        (organism == Organism::EcoliK12)
            ? CodonRateTable::build_rnap_ecoli()
            : CodonRateTable::build_rnap_human();

    // "Codons" are single nucleotides ("A", "G", "C", "U")
    std::vector<std::string> codons;
    codons.reserve(N);
    for (char nt : nt_sequence)
        codons.emplace_back(1, nt);

    // ── RibosomeElongation (RNAP kinetics) ───────────────────────────────────
    RibosomeElongation rnap(nt_sequence, codons, rnap_table,
                            K_RNAP_INI_DEFAULT, K_RNAP_TERM_DEFAULT);

    const auto& k_el = rnap.elongation_rates();

    // Harmonic mean rate (for pause detection at 15% threshold per Neuman 2003)
    // Filter NaN/Inf/zero values to prevent division-by-zero propagation
    double inv_sum = 0.0;
    int valid_count = 0;
    for (double k : k_el) {
        if (std::isfinite(k) && k > 1e-9) {
            inv_sum += 1.0 / k;
            ++valid_count;
        }
    }
    double hmean = (valid_count > 0 && inv_sum > 1e-15)
                   ? static_cast<double>(valid_count) / inv_sum
                   : MEAN_NT_RATE_ECOLI;

    // ── Pre-compute analytical arrival times (Zhao 2011 Eq. 7) ─────────────
    // <T_n> = 1/k_ini + Σ_{i=0}^{n-1} 1/k_i
    std::vector<double> t_arrival(N + 1);
    t_arrival[0] = 1.0 / rnap.k_ini();
    for (int n = 1; n <= N; ++n) {
        double k = (n - 1 < (int)k_el.size()) ? k_el[n - 1] : hmean;
        double inv_k = (std::isfinite(k) && k > 1e-9) ? 1.0 / k : 1.0 / hmean;
        t_arrival[n] = t_arrival[n - 1] + inv_k;
    }

    // ── Main growth loop (mirrors DualAssemblyEngine::run()) ─────────────────
    std::vector<GrowthStep> trajectory;
    trajectory.reserve(N);

    std::vector<double> dwell_traj;        // trajectory of dwell times for entropy
    dwell_traj.reserve(N);
    double dwell_cumulative = 0.0;

    for (int step = 0; step < N; ++step) {
        double k_n       = (step < (int)k_el.size()) ? k_el[step] : hmean;
        double dwell     = (std::isfinite(k_n) && k_n > 1e-9)
                           ? 1.0 / k_n
                           : 1.0 / hmean;
        double t_arr     = (step < (int)t_arrival.size()) ? t_arrival[step] : 0.0;

        // RNAP tunnel: first RNAP_TUNNEL_NT nt are inside the RNA:DNA hybrid
        bool   in_tunnel = (step < static_cast<int>(RNAP_TUNNEL_NT));

        // Pause detection (RNAP threshold = 15% of harmonic mean, per Neuman 2003)
        bool   is_pause  = !in_tunnel && (k_n < RNAP_PAUSE_THRESHOLD * hmean);

        // Co-translational folding probability: P_fold = k_fold/(k_fold + k_el)
        // Pechmann & Frydman (2013): pause sites get 3× k_fold boost
        double k_fold = K_FOLD_DEFAULT;
        if (is_pause) k_fold *= 3.0;
        double p_cotrans = k_fold / (k_fold + k_n);

        // Shannon entropy over dwell-time trajectory so far
        dwell_traj.push_back(dwell);
        dwell_cumulative += dwell;
        double S_dwell = shannon_entropy_bits(dwell_traj);

        const char* reg = (step < (int)regions.size()) ? regions[step] : "?";

        trajectory.push_back({
            step,
            nt_sequence[step],
            reg,
            t_arr,
            k_n,
            dwell,
            is_pause,
            in_tunnel,
            p_cotrans,
            S_dwell,
            dwell_cumulative
        });
    }

    return trajectory;
}

// ─── ASCII bar chart of elongation rates ─────────────────────────────────────
static std::string rate_bar(double k, double k_max, int width = 30) {
    int filled = static_cast<int>(std::round(k / k_max * width));
    filled = std::clamp(filled, 0, width);
    std::string bar(filled, '#');
    bar += std::string(width - filled, '.');
    return bar;
}

// ─── Print full trajectory table ─────────────────────────────────────────────
static void print_trajectory(const std::vector<GrowthStep>& traj,
                              double mean_rate, double hmean_rate)
{
    const int W = 120;
    std::string sep(W, '-');

    std::cout << COL_BOLD << "\n  NATURaL Co-Transcriptional Growth Trajectory\n"
              << COL_RESET << "  " << sep << "\n";

    std::cout << COL_BOLD
              << std::setw(4)  << "nt"
              << std::setw(4)  << "nt"
              << std::setw(10) << "Region"
              << std::setw(10) << "t_arr(s)"
              << std::setw(10) << "k_el(1/s)"
              << std::setw(11) << "dwell(ms)"
              << std::setw(8)  << "Tunnel?"
              << std::setw(8)  << "Pause?"
              << std::setw(10) << "P_fold"
              << std::setw(10) << "Hd(bits)"
              << "  Rate bar\n"
              << COL_RESET
              << "  " << sep << "\n";

    double k_max = 0.0;
    for (const auto& s : traj) k_max = std::max(k_max, s.k_el);

    for (const auto& s : traj) {
        bool pause = s.is_pause_site;
        bool tunnel = s.in_tunnel;

        std::string col = tunnel  ? COL_BLUE   :
                          pause   ? COL_YELLOW  :
                          (s.p_cotrans_folded > 0.05) ? COL_GREEN : COL_WHITE;

        std::string nt_display(1, s.nucleotide);

        std::cout << col
                  << std::setw(4)  << (s.residue_idx + 1)
                  << std::setw(4)  << nt_display
                  << std::setw(10) << region_name(s.region)
                  << std::setw(10) << std::fixed << std::setprecision(4) << s.t_arrival
                  << std::setw(10) << std::fixed << std::setprecision(2) << s.k_el
                  << std::setw(11) << std::fixed << std::setprecision(3) << (s.dwell_time * 1000.0)
                  << std::setw(8)  << (tunnel ? "yes" : " ")
                  << std::setw(8)  << (pause  ? "PAUSE" : " ")
                  << std::setw(10) << std::fixed << std::setprecision(4) << s.p_cotrans_folded
                  << std::setw(10) << std::fixed << std::setprecision(4) << s.shannon_entropy
                  << "  [" << rate_bar(s.k_el, k_max, 30) << "]\n"
                  << COL_RESET;
    }

    std::cout << "  " << sep << "\n";

    // Legend
    std::cout << "\n  Legend:\n"
              << "    " << COL_BLUE   << "Aqua"   << COL_RESET << "   = inside RNAP RNA:DNA hybrid tunnel (first "
                        << static_cast<int>(ribosome::RNAP_TUNNEL_NT) << " nt)\n"
              << "    " << COL_YELLOW << "Strawberry" << COL_RESET << " = RNAP pause site (k_el < "
                        << static_cast<int>(ribosome::RNAP_PAUSE_THRESHOLD * 100) << "% of harmonic mean "
                        << std::fixed << std::setprecision(1) << hmean_rate << " nt/s) → co-transcriptional folding window\n"
              << "    " << COL_GREEN  << "Mint"  << COL_RESET << "  = significant co-transcriptional folding probability (P_fold > 5%)\n"
              << "\n";
}

// ─── Summary statistics ───────────────────────────────────────────────────────
static void print_summary(const std::vector<GrowthStep>& traj,
                           const ribosome::RibosomeElongation& rnap)
{
    using namespace ribosome;

    std::cout << COL_BOLD << "\n  Simulation Summary\n" << COL_RESET;
    std::cout << "  " << std::string(60, '=') << "\n";

    int n_pause      = 0;
    int n_tunnel     = 0;
    int n_fold_win   = 0;  // folding windows (pause + outside tunnel)
    double total_dwell = 0.0;
    double max_p_fold  = 0.0;
    int    max_p_fold_idx = -1;

    for (const auto& s : traj) {
        if (s.in_tunnel)    ++n_tunnel;
        if (s.is_pause_site) { ++n_pause; ++n_fold_win; }
        total_dwell += s.dwell_time;
        if (s.p_cotrans_folded > max_p_fold) {
            max_p_fold     = s.p_cotrans_folded;
            max_p_fold_idx = s.residue_idx;
        }
    }

    double T_total = rnap.mean_total_time();
    auto valid     = validate_master_equation(
        std::min(static_cast<int>(traj.size()), 40), Organism::EcoliK12);

    std::cout << "  Sequence length      : " << traj.size() << " nt\n"
              << "  Total transcription  : " << std::fixed << std::setprecision(3)
              << T_total << " s (Zhao 2011 analytic, Eq. 7)\n"
              << "  Mean RNAP rate       : " << std::fixed << std::setprecision(1)
              << MEAN_NT_RATE_ECOLI << " nt/s (E. coli, Uptain 1997)\n"
              << "  Pause sites (RNAP)   : " << n_pause << " / " << traj.size()
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * n_pause / traj.size()) << "%)\n"
              << "  Tunnel residues      : " << n_tunnel
              << " (RNA:DNA hybrid, " << static_cast<int>(RNAP_TUNNEL_NT) << " nt)\n"
              << "  Co-trans fold windows: " << n_fold_win << "\n";

    if (max_p_fold_idx >= 0) {
        const auto& s = traj[max_p_fold_idx];
        std::cout << "  Best folding window  : nt " << (max_p_fold_idx + 1)
                  << " (" << s.nucleotide << ", " << region_name(s.region)
                  << ", P_fold=" << std::fixed << std::setprecision(4) << max_p_fold
                  << ", k_el=" << std::fixed << std::setprecision(2) << s.k_el << " nt/s)\n";
    }

    // Master equation validation
    std::cout << "\n  Master Equation Validation (Zhao 2011 Eq. 7):\n"
              << "    Analytic T  : " << std::fixed << std::setprecision(4) << valid.analytic_T << " s\n"
              << "    ODE T       : " << std::fixed << std::setprecision(4) << valid.ode_T      << " s\n"
              << "    Relative err: " << std::fixed << std::setprecision(2)
              << (valid.relative_error * 100.0) << "%  → "
              << (valid.passed ? COL_GREEN : COL_RED) << (valid.passed ? "PASS" : "FAIL")
              << COL_RESET << "\n";

    std::cout << "  " << std::string(60, '=') << "\n";
}

// ─── Folding window report ─────────────────────────────────────────────────────
static void print_folding_windows(const std::vector<GrowthStep>& traj,
                                   const ribosome::RibosomeElongation& rnap)
{
    auto fw = rnap.folding_windows(ribosome::K_FOLD_DEFAULT);
    if (fw.empty()) {
        std::cout << "\n  No co-transcriptional folding windows detected.\n";
        return;
    }

    std::cout << COL_BOLD << "\n  Co-Transcriptional Folding Windows (Zhao 2011 + Pechmann 2013)\n"
              << COL_RESET;
    std::cout << "  " << std::string(80, '-') << "\n";
    std::cout << COL_BOLD
              << std::setw(6)  << "nt"
              << std::setw(4)  << "nt"
              << std::setw(14) << "Region"
              << std::setw(12) << "t_avail(s)"
              << std::setw(10) << "k_fold"
              << std::setw(12) << "P_fold"
              << std::setw(8)  << "Pause?"
              << "\n" << COL_RESET
              << "  " << std::string(80, '-') << "\n";

    for (const auto& w : fw) {
        int idx = w.residue_idx;
        char nt  = (idx < static_cast<int>(std::strlen(HHRZ_SEQUENCE)))
                   ? HHRZ_SEQUENCE[idx] : '?';
        const char* reg = (idx < static_cast<int>(HHRZ_REGIONS.size()))
                          ? HHRZ_REGIONS[idx] : "?";

        std::string col = w.is_pause_site ? COL_YELLOW : COL_GREEN;

        std::cout << col
                  << std::setw(6)  << (idx + 1)
                  << std::setw(4)  << nt
                  << std::setw(14) << region_name(reg)
                  << std::setw(12) << std::fixed << std::setprecision(5) << w.t_available
                  << std::setw(10) << std::fixed << std::setprecision(3) << w.k_fold
                  << std::setw(12) << std::fixed << std::setprecision(4) << w.p_folded_cotrans
                  << std::setw(8)  << (w.is_pause_site ? "YES" : "no")
                  << "\n" << COL_RESET;
    }
    std::cout << "  " << std::string(80, '-') << "\n";
}

// ─── ASCII secondary structure map ───────────────────────────────────────────
static void print_structure_map(const std::vector<GrowthStep>& traj) {
    std::cout << COL_BOLD << "\n  HHRz Sequence & Structural Annotation\n" << COL_RESET;
    std::cout << "  " << std::string(80, '-') << "\n";

    // Print nucleotide sequence in groups of 10
    std::cout << "  Sequence:  ";
    for (int i = 0; i < HHRZ_LENGTH; ++i) {
        if (i > 0 && i % 10 == 0) std::cout << " ";
        std::cout << HHRZ_SEQUENCE[i];
    }
    std::cout << "\n";

    // Print region codes aligned
    std::cout << "  Regions:   ";
    for (int i = 0; i < HHRZ_LENGTH; ++i) {
        if (i > 0 && i % 10 == 0) std::cout << " ";
        const char* r = HHRZ_REGIONS[i];
        char ch = r[0];   // S or C or L
        // Colour per region
        if      (ch == 'S' && r[1] == '1') std::cout << COL_CYAN   << "1" << COL_RESET;
        else if (ch == 'S' && r[1] == '2') std::cout << COL_GREEN  << "2" << COL_RESET;
        else if (ch == 'S' && r[1] == '3') std::cout << COL_BLUE   << "3" << COL_RESET;
        else if (ch == 'C')                std::cout << COL_RED    << "C" << COL_RESET;
        else if (ch == 'L')                std::cout << COL_YELLOW << "L" << COL_RESET;
        else                               std::cout << "?";
    }
    std::cout << "\n";

    // Print ruler
    std::cout << "  Position:  ";
    for (int i = 0; i < HHRZ_LENGTH; ++i) {
        if (i > 0 && i % 10 == 0) std::cout << " ";
        if ((i + 1) % 10 == 0)      std::cout << (i + 1) / 10;
        else if ((i + 1) % 5 == 0)  std::cout << '|';
        else                         std::cout << '.';
    }
    std::cout << "\n\n";

    // Colour legend
    std::cout << "  "
              << COL_CYAN   << "1" << COL_RESET << "=Stem I  "
              << COL_GREEN  << "2" << COL_RESET << "=Stem II  "
              << COL_BLUE   << "3" << COL_RESET << "=Stem III  "
              << COL_RED    << "C" << COL_RESET << "=Catalytic core  "
              << COL_YELLOW << "L" << COL_RESET << "=Loop\n";

    // Mark RNAP tunnel and pause sites
    std::cout << "\n  RNAP markers:\n  ";
    for (int i = 0; i < HHRZ_LENGTH; ++i) {
        if (i > 0 && i % 10 == 0) std::cout << " ";
        bool in_tunnel = (i < static_cast<int>(ribosome::RNAP_TUNNEL_NT));
        bool is_pause  = i < static_cast<int>(traj.size()) && traj[i].is_pause_site;
        if      (in_tunnel) std::cout << COL_BLUE << "~" << COL_RESET;
        else if (is_pause)  std::cout << COL_YELLOW << "^" << COL_RESET;
        else                std::cout << " ";
    }
    std::cout << "\n";
    std::cout << "  (" << COL_BLUE << "~" << COL_RESET << "=tunnel  "
              << COL_YELLOW << "^" << COL_RESET << "=pause site)\n";

    std::cout << "  " << std::string(80, '-') << "\n";
}

// ─── Folding probability heat-map ─────────────────────────────────────────────
static void print_pfold_heatmap(const std::vector<GrowthStep>& traj) {
    std::cout << COL_BOLD << "\n  Co-Transcriptional Folding Probability (P_fold per nucleotide)\n"
              << COL_RESET;
    std::cout << "  " << std::string(80, '-') << "\n";

    // Show bar chart: each nt = one row with horizontal bar
    double p_max = 0.0;
    for (const auto& s : traj) p_max = std::max(p_max, s.p_cotrans_folded);

    for (int i = 0; i < static_cast<int>(traj.size()); ++i) {
        const auto& s = traj[i];

        // Bin probability into 0..20 bar width
        int bar_len = static_cast<int>(std::round(s.p_cotrans_folded / p_max * 20.0));
        bar_len = std::clamp(bar_len, 0, 20);

        std::string col = s.in_tunnel  ? COL_BLUE   :
                          s.is_pause_site ? COL_YELLOW  : COL_GREEN;
        std::cout << col
                  << std::setw(4)  << (i + 1)
                  << std::setw(2)  << s.nucleotide
                  << std::setw(14) << region_name(s.region)
                  << "  |";
        std::cout << std::string(bar_len, '|') << std::string(20 - bar_len, ' ')
                  << "| " << std::fixed << std::setprecision(4) << s.p_cotrans_folded
                  << COL_RESET << "\n";
    }
    std::cout << "  " << std::string(80, '-') << "\n";
}

// ─── Time-weighted thermodynamic score ────────────────────────────────────────
static void print_time_weighted_score(const ribosome::RibosomeElongation& rnap,
                                       const std::vector<GrowthStep>& traj)
{
    // Time-weighted average of p_cotrans_folded (proxy for co-translational efficiency)
    double tw_p = rnap.time_weighted_score([&](int n) {
        return (n < static_cast<int>(traj.size())) ? traj[n].p_cotrans_folded : 0.0;
    });

    std::cout << COL_BOLD << "\n  Time-Weighted Analysis (Zhao 2011 Eq. 12)\n" << COL_RESET;
    std::cout << "  " << std::string(60, '-') << "\n";
    std::cout << "  Time-weighted mean P_fold : " << std::fixed << std::setprecision(4) << tw_p << "\n";
    std::cout << "  Mean total transcription  : " << std::fixed << std::setprecision(3)
              << rnap.mean_total_time() << " s\n";
    std::cout << "  1/k_ini (initiation time) : " << std::fixed << std::setprecision(3)
              << (1.0 / rnap.k_ini()) << " s\n";
    std::cout << "  1/k_ter (termination time): " << std::fixed << std::setprecision(3)
              << (1.0 / rnap.k_ter()) << " s\n";
    std::cout << "  " << std::string(60, '-') << "\n";
}

// ─── main ─────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    using namespace ribosome;

    (void)argc; (void)argv;

    // ── Header ───────────────────────────────────────────────────────────────
    std::cout << "\n"
              << COL_BOLD << COL_CYAN
              << "  ╔══════════════════════════════════════════════════════════╗\n"
              << "  ║   NATURaL — Co-Transcriptional Folding Simulation        ║\n"
              << "  ║   Hammerhead Ribozyme (Schistosoma mansoni, HHRz-40)     ║\n"
              << "  ║   Zhao 2011 master equation  ·  E. coli RNAP rates       ║\n"
              << "  ╚══════════════════════════════════════════════════════════╝\n"
              << COL_RESET << "\n";

    std::cout << "  Sequence  : " << COL_BOLD << HHRZ_SEQUENCE << COL_RESET
              << "  (" << HHRZ_LENGTH << " nt)\n"
              << "  Reference : Martick & Scott (2006) Science 313:1514  [PDB: 2OEU]\n"
              << "  Mode      : co-transcriptional (RNAP, E. coli K-12, ~"
              << static_cast<int>(MEAN_NT_RATE_ECOLI) << " nt/s)\n"
              << "  Kinetics  : Zhao 2011 master equation, RNAP rate table\n"
              << "              (Uptain 1997; Neuman 2003; Jonkers 2014)\n"
              << "  Tunnel    : first " << static_cast<int>(RNAP_TUNNEL_NT)
              << " nt occluded in RNAP RNA:DNA hybrid (Nudler 2012)\n"
              << "  Pausing   : threshold = " << static_cast<int>(RNAP_PAUSE_THRESHOLD * 100)
              << "% of mean rate (RNAP; Neuman 2003)\n"
              << "  k_fold    : " << K_FOLD_DEFAULT << " s⁻¹ base (×3 at pause sites; Pechmann 2013)\n\n";

    // ── Build the elongation model ──────────────────────────────────────────
    std::string nt_seq(HHRZ_SEQUENCE);
    std::vector<std::string> codons;
    codons.reserve(HHRZ_LENGTH);
    for (char nt : nt_seq) codons.emplace_back(1, nt);

    CodonRateTable rnap_table = CodonRateTable::build_rnap_ecoli();
    RibosomeElongation rnap(nt_seq, codons, rnap_table,
                            K_RNAP_INI_DEFAULT, K_RNAP_TERM_DEFAULT);

    // Harmonic mean for threshold computation
    const auto& k_el_vec = rnap.elongation_rates();
    double inv_sum = 0.0;
    for (double k : k_el_vec) if (k > 1e-9) inv_sum += 1.0 / k;
    double hmean = HHRZ_LENGTH / inv_sum;

    // ── Run NATURaL simulation ───────────────────────────────────────────────
    auto traj = run_rnap_natural(nt_seq, HHRZ_REGIONS, Organism::EcoliK12);

    // ── Print outputs ────────────────────────────────────────────────────────
    print_structure_map(traj);
    print_trajectory(traj, rnap_table.mean_rate_aa_per_s, hmean);
    print_pfold_heatmap(traj);
    print_folding_windows(traj, rnap);
    print_time_weighted_score(rnap, traj);
    print_summary(traj, rnap);

    // ── Master equation ODE integration (full trajectory) ───────────────────
    {
        std::cout << COL_BOLD << "\n  Master Equation ODE Integration (full)\n" << COL_RESET;
        std::cout << "  " << std::string(60, '-') << "\n";
        double t_max = rnap.mean_total_time() * 3.0;
        auto state = rnap.integrate(t_max, 1e-4, /*adaptive=*/true);
        std::cout << "  Integration time: t_max=" << std::fixed << std::setprecision(3)
                  << t_max << " s  (3 × T_mean)\n";
        std::cout << "  Ribosome occupancy at completion (P[N]): "
                  << std::fixed << std::setprecision(4) << state.P.back() << "\n";
        std::cout << "  Fraction reached last nt: "
                  << std::fixed << std::setprecision(4)
                  << state.fraction_reached(HHRZ_LENGTH) << "\n";

        // Print the first-passage times for key structural elements
        std::cout << "\n  Analytical mean first-passage times (Zhao 2011 Eq. 7):\n";

        // Stem I (nt 5 is the last nt of Stem I)
        // Stem III end (nt 12)
        // Catalytic core (nt 17, cleavage site)
        // Stem II (nt 35)
        // Full molecule (nt 40)
        const struct { int nt; const char* label; } landmarks[] = {
            {  5, "Stem I    nt  5 (last Stem-I opener)"},
            { 12, "Stem III  nt 12 (Stem-III complete)"},
            { 17, "Core      nt 17 (C17 cleavage site)"},
            { 22, "Core end  nt 22 (core complete)"},
            { 27, "Stem III  nt 27 (Stem-III closure)"},
            { 35, "Stem II   nt 35 (Stem-II closure)"},
            { 40, "Full HHRz nt 40 (complete ribozyme)"},
        };
        for (const auto& lm : landmarks) {
            int idx = lm.nt - 1;
            if (idx >= 0 && idx < HHRZ_LENGTH) {
                double t = rnap.mean_arrival_time(idx);
                std::cout << "    " << lm.label << "  →  t=" << std::fixed
                          << std::setprecision(4) << t << " s\n";
            }
        }
        std::cout << "  " << std::string(60, '-') << "\n";
    }

    // ── Final message ────────────────────────────────────────────────────────
    std::cout << "\n"
              << COL_BOLD << COL_CYAN
              << "  NATURaL simulation complete.\n" << COL_RESET
              << "  The hammerhead ribozyme folds co-transcriptionally across "
              << HHRZ_LENGTH << " RNAP steps\n"
              << "  using the Zhao 2011 master equation in E. coli RNAP mode.\n"
              << "  RNAP pause sites mark windows where RNA secondary structure\n"
              << "  (Stems I, II, III and the catalytic core) can form before\n"
              << "  downstream sequence is transcribed.\n\n";

    return 0;
}
