// RibosomeElongation.cpp — Zhao 2011 master equation implementation
//
// Reference:
//   Zhao et al. (2011) "Mechanism of Ribosomal Translocation"
//   J. Phys. Chem. B 115, 3987–3997.  DOI: 10.1021/jp109255g
//
// Codon rates calibrated from:
//   Dong et al. (1996) J Mol Biol 260:649 (E. coli tRNA gene copies)
//   Ingolia et al. (2011) Cell 147:789 (human ribosome profiling)
//   Wohlgemuth et al. (2008) EMBO J 27:1458 (mean elongation rate)
//   Goldman et al. (2010) Cell 143:92 (tunnel length 34–36 aa)

#include "RibosomeElongation.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace ribosome {

// ─── CodonRateTable ───────────────────────────────────────────────────────────

CodonRateTable CodonRateTable::build_ecoli() {
    CodonRateTable t(Organism::EcoliK12, MEAN_EL_RATE_ECOLI);
    t.populate_ecoli();
    return t;
}

CodonRateTable CodonRateTable::build_human() {
    CodonRateTable t(Organism::HumanHEK293, MEAN_EL_RATE_HUMAN);
    t.populate_human();
    return t;
}

CodonRateTable::CodonRateTable(Organism org, double mean_rate)
    : organism(org), mean_rate_aa_per_s(mean_rate) {}

double CodonRateTable::rate(const std::string& codon) const {
    auto it = rates_.find(codon);
    return (it != rates_.end()) ? it->second : mean_rate_aa_per_s;
}

double CodonRateTable::mean_rate(const std::vector<std::string>& codons) const {
    if (codons.empty()) return mean_rate_aa_per_s;
    double sum = 0.0;
    for (const auto& c : codons) sum += rate(c);
    return sum / static_cast<double>(codons.size());
}

std::vector<int> CodonRateTable::pause_sites(
    const std::vector<std::string>& codons, double threshold) const
{
    double mr = mean_rate(codons);
    std::vector<int> pauses;
    for (int i = 0; i < (int)codons.size(); ++i)
        if (rate(codons[i]) < mr * threshold)
            pauses.push_back(i);
    return pauses;
}

// ─── E. coli K-12 codon rate table ────────────────────────────────────────────
// Rates (s⁻¹) derived from tRNA gene copy numbers (Dong 1996 Table 1) scaled
// so that the harmonic mean equals MEAN_EL_RATE_ECOLI = 16.5 aa/s.
// Rare codons (AGG, AGA, AUA, CUA, CGA, CCC) have rates ≤ 3 s⁻¹.
// Common codons (AAA, GAA, CTG) have rates ≈ 40–60 s⁻¹.
void CodonRateTable::populate_ecoli() {
    // Format: codon → rate (s⁻¹)
    // Data from Dong (1996) tRNA gene copies, calibrated to Wohlgemuth (2008)
    rates_ = {
        // Phe
        {"TTT", 12.1}, {"TTC", 24.2},
        // Leu
        {"TTA",  8.5}, {"TTG", 14.0}, {"CTT", 10.5}, {"CTC", 12.0},
        {"CTA",  2.8}, {"CTG", 52.0},
        // Ile
        {"ATT", 28.0}, {"ATC", 35.0}, {"ATA",  2.5},
        // Met
        {"ATG", 28.0},
        // Val
        {"GTT", 28.0}, {"GTC", 18.0}, {"GTA", 15.0}, {"GTG", 30.0},
        // Ser
        {"TCT", 18.0}, {"TCC", 20.0}, {"TCA", 12.0}, {"TCG",  8.0},
        {"AGT", 10.0}, {"AGC", 22.0},
        // Pro
        {"CCT", 14.0}, {"CCC",  3.5}, {"CCA", 16.0}, {"CCG", 30.0},
        // Thr
        {"ACT", 18.0}, {"ACC", 32.0}, {"ACA", 10.0}, {"ACG", 20.0},
        // Ala
        {"GCT", 22.0}, {"GCC", 24.0}, {"GCA", 18.0}, {"GCG", 30.0},
        // Tyr
        {"TAT", 16.0}, {"TAC", 28.0},
        // Stop (not used in elongation)
        {"TAA",  0.0}, {"TAG",  0.0}, {"TGA",  0.0},
        // His
        {"CAT", 12.0}, {"CAC", 16.0},
        // Gln
        {"CAA", 20.0}, {"CAG", 42.0},
        // Asn
        {"AAT", 18.0}, {"AAC", 35.0},
        // Lys
        {"AAA", 50.0}, {"AAG", 28.0},
        // Asp
        {"GAT", 30.0}, {"GAC", 35.0},
        // Glu
        {"GAA", 58.0}, {"GAG", 28.0},
        // Cys
        {"TGT",  8.0}, {"TGC", 12.0},
        // Trp
        {"TGG", 14.0},
        // Arg
        {"CGT", 38.0}, {"CGC", 30.0}, {"CGA",  2.5}, {"CGG",  5.0},
        {"AGA",  2.0}, {"AGG",  1.5},
        // Gly
        {"GGT", 30.0}, {"GGC", 38.0}, {"GGA",  8.0}, {"GGG", 10.0},
    };
    // Rescale so harmonic mean = MEAN_EL_RATE_ECOLI
    // Compute current harmonic mean over elongating codons (exclude stops)
    double inv_sum = 0.0; int cnt = 0;
    for (auto& kv : rates_) {
        if (kv.second > 0.1) { inv_sum += 1.0 / kv.second; ++cnt; }
    }
    if (cnt > 0) {
        double current_hmean = cnt / inv_sum;
        double scale = MEAN_EL_RATE_ECOLI / current_hmean;
        for (auto& kv : rates_)
            if (kv.second > 0.1) kv.second *= scale;
    }
    mean_rate_aa_per_s = MEAN_EL_RATE_ECOLI;
}

// ─── Human codon rate table ───────────────────────────────────────────────────
// Calibrated from Ingolia (2011) ribosome profiling A-site density.
// Mean rate: 5.6 aa/s (Ingolia 2011 Fig. 3).
void CodonRateTable::populate_human() {
    // Relative rates from Ingolia (2011) ribosome density (inverse occupancy)
    // scaled so harmonic mean = MEAN_EL_RATE_HUMAN = 5.6 aa/s
    rates_ = {
        // Phe
        {"TTT",  4.2}, {"TTC",  7.8},
        // Leu
        {"TTA",  2.1}, {"TTG",  4.5}, {"CTT",  4.0}, {"CTC",  5.5},
        {"CTA",  1.8}, {"CTG", 14.0},
        // Ile
        {"ATT",  7.0}, {"ATC", 10.5}, {"ATA",  2.8},
        // Met
        {"ATG",  8.0},
        // Val
        {"GTT",  5.5}, {"GTC",  7.0}, {"GTA",  3.5}, {"GTG",  9.0},
        // Ser
        {"TCT",  5.0}, {"TCC",  7.5}, {"TCA",  4.2}, {"TCG",  2.5},
        {"AGT",  4.0}, {"AGC",  8.0},
        // Pro
        {"CCT",  5.5}, {"CCC",  6.5}, {"CCA",  5.8}, {"CCG",  2.0},
        // Thr
        {"ACT",  5.2}, {"ACC", 10.0}, {"ACA",  4.5}, {"ACG",  2.8},
        // Ala
        {"GCT",  6.5}, {"GCC", 10.5}, {"GCA",  5.0}, {"GCG",  2.2},
        // Tyr
        {"TAT",  5.0}, {"TAC",  9.0},
        // Stop
        {"TAA",  0.0}, {"TAG",  0.0}, {"TGA",  0.0},
        // His
        {"CAT",  4.5}, {"CAC",  7.5},
        // Gln
        {"CAA",  5.5}, {"CAG", 13.0},
        // Asn
        {"AAT",  5.5}, {"AAC", 10.5},
        // Lys
        {"AAA", 12.0}, {"AAG", 14.5},
        // Asp
        {"GAT",  8.5}, {"GAC", 12.0},
        // Glu
        {"GAA", 14.5}, {"GAG", 13.0},
        // Cys
        {"TGT",  3.5}, {"TGC",  5.5},
        // Trp
        {"TGG",  5.0},
        // Arg
        {"CGT",  3.5}, {"CGC",  6.5}, {"CGA",  2.5}, {"CGG",  5.0},
        {"AGA",  6.5}, {"AGG",  5.5},
        // Gly
        {"GGT",  5.0}, {"GGC",  9.0}, {"GGA",  6.5}, {"GGG",  5.5},
    };
    // Rescale to human mean
    double inv_sum = 0.0; int cnt = 0;
    for (auto& kv : rates_)
        if (kv.second > 0.1) { inv_sum += 1.0 / kv.second; ++cnt; }
    if (cnt > 0) {
        double current_hmean = cnt / inv_sum;
        double scale = MEAN_EL_RATE_HUMAN / current_hmean;
        for (auto& kv : rates_)
            if (kv.second > 0.1) kv.second *= scale;
    }
    mean_rate_aa_per_s = MEAN_EL_RATE_HUMAN;
}

// ─── RNAP rate tables (co-transcriptional, same Zhao 2011 master equation) ────
// The Zhao 2011 framework is identical for RNA polymerase: each nucleotide
// incorporation event n has rate k_n (nt/s), giving dP_n/dt = k_{n-1}*P_{n-1}
// - k_n*P_n.  "Codons" are single nucleotides: "A", "U"/"T", "G", "C".
//
// E. coli RNAP NTP incorporation rates (Uptain 1997 Annu. Rev. Biochem. 66:117;
// Neuman 2003 Science 298:1152; Bai 2004 Cell 119:785):
//   – Mean elongation: ~50 nt/s at saturating NTP
//   – Purine (A, G) slightly faster than pyrimidine (U, C)
//   – Pause signals at U-rich regions (~15% threshold vs 30% for ribosomes)
CodonRateTable CodonRateTable::build_rnap_ecoli() {
    CodonRateTable t(Organism::EcoliK12, MEAN_NT_RATE_ECOLI);
    t.populate_rnap_ecoli();
    return t;
}

CodonRateTable CodonRateTable::build_rnap_human() {
    CodonRateTable t(Organism::HumanHEK293, MEAN_NT_RATE_HUMAN);
    t.populate_rnap_human();
    return t;
}

void CodonRateTable::populate_rnap_ecoli() {
    // E. coli RNAP: single-nucleotide incorporation rates (s⁻¹)
    // Purines are ~10% faster; U is slowest (Uptain 1997; Bai 2004).
    // Raw rates (before rescaling to harmonic mean = 50 nt/s):
    rates_ = {
        {"A",  56.0},  // fast — abundant ATP
        {"G",  52.0},  // fast — GTP
        {"C",  46.0},  // slightly slower — CTP
        {"U",  44.0},  // slowest — UTP; also T for DNA template
        {"T",  44.0},  // DNA-mode alias for U
        // Di-nucleotide pause signals (RNAP can pause at specific sequence motifs)
        // These are treated as separate "codons" for pause detection:
        // GGU, TGT → near-pause sites observed by Neuman 2003
        {"GGU", 7.0},  // pause sequence (RNA hairpin precursor)
        {"TGT", 6.5},  // pause sequence (DNA template)
    };
    // Rescale single-nucleotide rates to harmonic mean = MEAN_NT_RATE_ECOLI
    double inv_sum = 0.0; int cnt = 0;
    for (auto& kv : rates_) {
        if (kv.first.size() == 1 && kv.second > 0.1) {
            inv_sum += 1.0 / kv.second;
            ++cnt;
        }
    }
    if (cnt > 0) {
        double current_hmean = cnt / inv_sum;
        double scale = MEAN_NT_RATE_ECOLI / current_hmean;
        for (auto& kv : rates_)
            if (kv.first.size() == 1 && kv.second > 0.1)
                kv.second *= scale;
    }
    mean_rate_aa_per_s = MEAN_NT_RATE_ECOLI;
}

void CodonRateTable::populate_rnap_human() {
    // Human RNAP II: slower than prokaryotic RNAP (Jonkers 2014 Mol. Cell 54:591)
    // Mean ~25 nt/s; pause propensity higher at CpG dinucleotides and
    // splice sites (Carrillo Oesterreich 2010 Science 327:1914).
    rates_ = {
        {"A",  28.0},
        {"G",  26.0},
        {"C",  22.0},  // CpG pausing modelled via separate di-nt pause signals
        {"U",  24.0},
        {"T",  24.0},
        // Splice-site pause proxies (intronic pyrimidine tracts):
        {"CG", 3.5},   // CpG: RNAP II pausing at methylated CpG (Carrillo 2010)
        {"GT", 4.0},   // GT splice donor signal
    };
    double inv_sum = 0.0; int cnt = 0;
    for (auto& kv : rates_) {
        if (kv.first.size() == 1 && kv.second > 0.1) {
            inv_sum += 1.0 / kv.second;
            ++cnt;
        }
    }
    if (cnt > 0) {
        double current_hmean = cnt / inv_sum;
        double scale = MEAN_NT_RATE_HUMAN / current_hmean;
        for (auto& kv : rates_)
            if (kv.first.size() == 1 && kv.second > 0.1)
                kv.second *= scale;
    }
    mean_rate_aa_per_s = MEAN_NT_RATE_HUMAN;
}

// ─── CodonRateTable::aa_to_mfc ───────────────────────────────────────────────
// Most-Frequent Codon (MFC) per amino acid, derived from highest-rate codons in
// each organism's rate table.
// E. coli K-12: based on Ikemura 1985 tRNA gene dosage; fast codons correlate with
//   highly expressed proteins (CAI-optimal codons).
// Human HEK293: based on Grantham 1980 / Nakamura 2000 codon usage database.
std::unordered_map<char, std::string> CodonRateTable::aa_to_mfc(Organism org)
{
    if (org == Organism::EcoliK12) {
        return {
            {'A', "GCG"}, // Ala — GCG most abundant in E. coli (Ikemura 1985)
            {'R', "CGT"}, // Arg — CGT/CGC dominant (CGG/AGA/AGG are rare)
            {'N', "AAC"}, // Asn — AAC > AAT in E. coli
            {'D', "GAT"}, // Asp — GAT > GAC
            {'C', "TGC"}, // Cys — TGC > TGT
            {'Q', "CAG"}, // Gln — CAG dominant
            {'E', "GAA"}, // Glu — GAA fast; GAG rare
            {'G', "GGC"}, // Gly — GGC/GGT dominate
            {'H', "CAC"}, // His — CAC > CAT
            {'I', "ATC"}, // Ile — ATC fast; ATA is rare codon (2.5 s⁻¹)
            {'L', "CTG"}, // Leu — CTG is the dominant fast Leu codon
            {'K', "AAA"}, // Lys — AAA fast (50 s⁻¹); AAG slightly slower
            {'M', "ATG"}, // Met — only one codon
            {'F', "TTT"}, // Phe — TTT > TTC in E. coli
            {'P', "CCG"}, // Pro — CCG dominant in E. coli
            {'S', "AGC"}, // Ser — AGC common; TCG/TCA rare
            {'T', "ACC"}, // Thr — ACC > ACG > ACA; ACT rare
            {'W', "TGG"}, // Trp — only one codon
            {'Y', "TAT"}, // Tyr — TAT > TAC
            {'V', "GTG"}, // Val — GTG/GTC fast; GTA rare (2.5 s⁻¹)
            {'X', ""   }, // unknown → use mean
        };
    } else {
        // Human HEK293 — C-ending codons dominate (Grantham 1980 genome hypothesis)
        return {
            {'A', "GCC"}, // Ala — GCC dominant in human
            {'R', "AGA"}, // Arg — AGA/AGG more common in human than CGN
            {'N', "AAC"}, // Asn — AAC > AAT
            {'D', "GAC"}, // Asp — GAC > GAT in human
            {'C', "TGC"}, // Cys — TGC > TGT
            {'Q', "CAG"}, // Gln — CAG dominant
            {'E', "GAG"}, // Glu — GAG > GAA in human
            {'G', "GGC"}, // Gly — GGC common
            {'H', "CAC"}, // His — CAC > CAT
            {'I', "ATC"}, // Ile — ATC dominant
            {'L', "CTG"}, // Leu — CTG dominant across eukaryotes
            {'K', "AAG"}, // Lys — AAG > AAA in human (Nakamura 2000)
            {'M', "ATG"}, // Met — only one codon
            {'F', "TTC"}, // Phe — TTC > TTT in human
            {'P', "CCC"}, // Pro — CCC common in human
            {'S', "AGC"}, // Ser — AGC most common in human
            {'T', "ACC"}, // Thr — ACC dominant
            {'W', "TGG"}, // Trp — only one codon
            {'Y', "TAC"}, // Tyr — TAC > TAT in human
            {'V', "GTG"}, // Val — GTG dominant
            {'X', ""   }, // unknown → use mean
        };
    }
}

// ─── RibosomeElongation ───────────────────────────────────────────────────────

RibosomeElongation::RibosomeElongation(
    const std::string&              aa_sequence,
    const std::vector<std::string>& codons,
    const CodonRateTable&           rate_table,
    double k_ini, double k_ter)
    : k_ini_(k_ini), k_ter_(k_ter)
{
    int N = static_cast<int>(aa_sequence.size());
    k_el_.resize(N);

    for (int n = 0; n < N; ++n) {
        if (n < (int)codons.size())
            k_el_[n] = rate_table.rate(codons[n]);
        else
            k_el_[n] = rate_table.mean_rate_aa_per_s;
    }

    mean_rate_ = 0.0;
    double inv_sum = 0.0;
    int valid_count = 0;
    for (double k : k_el_) {
        if (std::isfinite(k) && k > 1e-9) {
            inv_sum += 1.0 / k;
            ++valid_count;
        }
    }
    mean_rate_ = (valid_count > 0 && inv_sum > 1e-15)
                 ? static_cast<double>(valid_count) / inv_sum
                 : rate_table.mean_rate_aa_per_s;

    // Flag pause sites: rate < RIBOSOME_PAUSE_THRESHOLD of mean
    double threshold = RIBOSOME_PAUSE_THRESHOLD * mean_rate_;
    for (int n = 0; n < N; ++n)
        if (k_el_[n] < threshold)
            pause_sites_.push_back(n);
}

// ─── Analytical mean first-passage time (Zhao 2011 Eq. 7) ────────────────────
// <T_n> = 1/k_ini + Σ_{i=0}^{n-1} 1/k_el[i]
double RibosomeElongation::mean_arrival_time(int n) const {
    double t = 1.0 / k_ini_;
    for (int i = 0; i < n && i < (int)k_el_.size(); ++i) {
        double ki = k_el_[i];
        t += (std::isfinite(ki) && ki > 1e-9) ? 1.0 / ki : 1.0 / mean_rate_;
    }
    return t;
}

double RibosomeElongation::mean_total_time() const {
    return mean_arrival_time(static_cast<int>(k_el_.size())) + 1.0 / k_ter_;
}

// ─── Master equation integration (explicit Euler, adaptive dt) ───────────────
// Zhao 2011 Eqs. 3–5:
//   dP_0/dt = k_ini · δ(t) − k_el[0] · P_0
//   dP_n/dt = k_el[n-1] · P_{n-1} − k_el[n] · P_n   (1 ≤ n < N)
//   dP_N/dt = k_el[N-1] · P_{N-1} − k_ter · P_N
//
// Uses an absorbing post-termination state so probability is conserved across
// initiation, elongation, termination, and completion.
MasterEqState RibosomeElongation::integrate(double t_max, double dt, bool adaptive) const {
    int N = static_cast<int>(k_el_.size());
    MasterEqState state;
    state.n_residues = N;
    state.t_current  = 0.0;
    state.k_el       = k_el_;
    state.k_ini      = k_ini_;
    state.k_ter      = k_ter_;
    // P[0] = pre-initiation, P[1..N] = elongation stages,
    // P[N+1] = termination stage, P[N+2] = absorbing completed chain.
    state.P.assign(N + 3, 0.0);
    state.t_arrival.resize(N + 1, -1.0);

    // Inject probability at position 0 at t=0 (initiation)
    // Model: ribosome starts outside; rate k_ini brings it to position 1
    // Here we use the discrete-chain formulation with P[0] = pre-initiation
    std::vector<double> P_init(N + 3, 0.0);
    P_init[0] = 1.0; // all ribosomes at initiation complex
    state.P = P_init;

    // Adaptive step size: dt_max = 0.5 / max(k_el)
    double k_max = *std::max_element(k_el_.begin(), k_el_.end());
    k_max = std::max(k_max, std::max(k_ini_, k_ter_));
    if (adaptive) dt = std::min(dt, 0.5 / k_max);

    double t = 0.0;
    std::vector<double> dP(N + 3, 0.0);

    // Use Eigen for the tridiagonal ODE vector update
    Eigen::VectorXd P_vec(N + 3), dP_vec(N + 3);
    for (int i = 0; i <= N + 2; ++i) P_vec(i) = P_init[i];

    while (t < t_max) {
        // dP[0] = -k_ini * P[0]
        dP_vec(0) = -k_ini_ * P_vec(0);
        // dP[1] = k_ini*P[0] - k_el[0]*P[1]
        dP_vec(1) = k_ini_ * P_vec(0) - k_el_[0] * P_vec(1);
        // dP[n] = k_el[n-2]*P[n-1] - k_el[n-1]*P[n]  (n = 2..N)
        for (int n = 2; n <= N; ++n)
            dP_vec(n) = k_el_[n - 2] * P_vec(n - 1) - k_el_[n - 1] * P_vec(n);
        // dP[N+1] = k_el[N-1]*P[N] - k_ter*P[N+1]
        dP_vec(N + 1) = k_el_[N - 1] * P_vec(N) - k_ter_ * P_vec(N + 1);
        // dP[N+2] = k_ter*P[N+1]  (absorbing completed chain)
        dP_vec(N + 2) = k_ter_ * P_vec(N + 1);

        P_vec += dt * dP_vec;
        // Clamp non-negative
        P_vec = P_vec.cwiseMax(0.0);
        t += dt;

        // Leading-edge: record when intermediate P[n] first exceeds 1e-3
        for (int n = 0; n < N; ++n)
            if (state.t_arrival[n] < 0 && P_vec(n) > 1e-3)
                state.t_arrival[n] = t;
        // Terminal absorbing state: record the median completion time.
        if (state.t_median_terminal < 0 && P_vec(N + 2) >= 0.5)
            state.t_median_terminal = t;
    }
    for (int i = 0; i <= N + 2; ++i) state.P[i] = P_vec(i);

    state.t_current = t;

    // Fill any unrecorded arrival times with analytic estimate
    for (int n = 0; n <= N; ++n)
        if (state.t_arrival[n] < 0)
            state.t_arrival[n] = mean_arrival_time(n);

    return state;
}

double MasterEqState::fraction_reached(int n) const {
    if (n < 0 || n > n_residues) return 0.0;
    // Sum of P[n..N] = probability that ribosome has passed position n
    double sum = 0.0;
    for (int i = n; i < static_cast<int>(P.size()); ++i) sum += P[i];
    return std::min(sum, 1.0);
}

// ─── Co-translational folding windows ────────────────────────────────────────
// A folding window exists at residue n if:
//   1. n > TUNNEL_LENGTH_AA (residue has emerged from ribosome tunnel)
//   2. The dwell time at position n (1/k_el[n]) is long enough for folding
// Competitive folding fraction (Hartl 2011; co-translational = k_fold/(k_fold+k_el)):
std::vector<RibosomeElongation::FoldingWindow>
RibosomeElongation::folding_windows(double k_fold_base) const
{
    std::vector<FoldingWindow> windows;
    int N = static_cast<int>(k_el_.size());

    for (int n = 0; n < N; ++n) {
        if (n < static_cast<int>(TUNNEL_LENGTH_AA)) continue; // still in tunnel

        double k_el_n = k_el_[n];
        double dwell  = (std::isfinite(k_el_n) && k_el_n > 1e-9)
                        ? 1.0 / k_el_n
                        : 1.0 / mean_rate_;  // fallback to mean dwell time

        // Scale folding rate by secondary structure propensity (simple heuristic:
        // folding is faster at pause sites where the chain has time to equilibrate)
        double k_fold = k_fold_base;
        bool is_pause = (std::isfinite(k_el_n) && k_el_n < RIBOSOME_PAUSE_THRESHOLD * mean_rate_);
        if (is_pause) k_fold *= 3.0; // pause sites have 3× folding rate (Pechmann 2013)

        double p_cotrans = k_fold / (k_fold + k_el_n);

        windows.push_back({
            n,
            dwell,
            k_fold,
            p_cotrans,
            is_pause
        });
    }
    return windows;
}

// ─── Time-weighted thermodynamic score ───────────────────────────────────────
// Integrates S(n) over elongating intermediates, weighted by dwell time 1/k_n.
// Initiation/termination are excluded because score_fn is defined only on chain
// intermediates; a constant score must remain constant after normalization.
double RibosomeElongation::time_weighted_score(
    const std::function<double(int)>& score_fn) const
{
    double weighted_sum = 0.0;
    double dwell_sum = 0.0;
    int N = static_cast<int>(k_el_.size());
    for (int n = 0; n < N; ++n) {
        double ki = k_el_[n];
        double dwell = (std::isfinite(ki) && ki > 1e-9) ? 1.0 / ki : 1.0 / mean_rate_;
        weighted_sum += score_fn(n) * dwell;
        dwell_sum += dwell;
    }
    if (dwell_sum < 1e-12) return 0.0;
    return weighted_sum / dwell_sum; // time-weighted average over scored intermediates
}

// ─── Validation (Zhao 2011 internal consistency check) ───────────────────────
// Analytic T vs ODE T should agree within 5% for well-resolved rate differences.
ValidationResult validate_master_equation(int n_residues, Organism org)
{
    // Build a synthetic uniform-rate protein
    CodonRateTable table = (org == Organism::EcoliK12)
                           ? CodonRateTable::build_ecoli()
                           : CodonRateTable::build_human();

    double k_mean = table.mean_rate_aa_per_s;
    std::string aa_seq(n_residues, 'A'); // poly-Ala as neutral test

    RibosomeElongation model(aa_seq, {}, table, K_INI_DEFAULT, K_TERM_DEFAULT);

    double T_analytic = model.mean_total_time();

    // Numerical mean first-passage time from the same master equation:
    // E[T_abs] = integral_0^inf P(not completed at t) dt.
    // This validates the kinetic ODE against the analytic sum without confusing
    // the median terminal occupancy with the mean of an initiation-dominated
    // sequential process.
    const double t_max = T_analytic * 12.0;
    double dt = 1.0 / (k_mean * 100.0);
    double k_max = std::max({k_mean, K_INI_DEFAULT, K_TERM_DEFAULT});
    dt = std::min(dt, 0.25 / k_max);

    const int N = n_residues;
    Eigen::VectorXd P = Eigen::VectorXd::Zero(N + 3);
    Eigen::VectorXd dP = Eigen::VectorXd::Zero(N + 3);
    P(0) = 1.0;

    double T_ode = 0.0;
    for (double t = 0.0; t < t_max; t += dt) {
        double survival = std::max(0.0, 1.0 - P(N + 2));
        T_ode += survival * dt;

        dP.setZero();
        dP(0) = -K_INI_DEFAULT * P(0);
        dP(1) = K_INI_DEFAULT * P(0) - k_mean * P(1);
        for (int n = 2; n <= N; ++n)
            dP(n) = k_mean * P(n - 1) - k_mean * P(n);
        dP(N + 1) = k_mean * P(N) - K_TERM_DEFAULT * P(N + 1);
        dP(N + 2) = K_TERM_DEFAULT * P(N + 1);

        P += dt * dP;
        P = P.cwiseMax(0.0);
    }

    double rel_err = std::abs(T_ode - T_analytic) / T_analytic;
    bool passed   = rel_err < 0.05; // 5% tolerance

    std::ostringstream msg;
    msg << "RibosomeElongation validation ["
        << (org == Organism::EcoliK12 ? "E.coli" : "Human") << ", N="
        << n_residues << "]: "
        << "T_analytic=" << T_analytic << "s, T_ode=" << T_ode << "s, "
        << "rel_err=" << rel_err * 100.0 << "% — "
        << (passed ? "PASS" : "FAIL");

    return { passed, T_analytic, T_ode, rel_err, msg.str() };
}

} // namespace ribosome
