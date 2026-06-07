// NATURaLDualAssembly.cpp — co-translational / co-transcriptional folding
//
// The same Zhao 2011 master-equation rationale applies to both polypeptide
// and nucleotide chains:
//
//   Polypeptide (ribosome):
//     dP_n/dt = k_{n-1} · P_{n-1} − k_n · P_n
//     k_n  = codon-specific aa elongation rate (s⁻¹)  [Dong 1996 / Ingolia 2011]
//     Tunnel: 34 aa occluded (Goldman 2010 Cell 143:92)
//
//   Nucleotide chain (RNA polymerase):
//     Same master equation, k_n = nucleotide-specific NTP incorporation rate
//     k_n  from Uptain 1997 (E. coli RNAP) / Jonkers 2014 (Human RNAP II)
//     Tunnel: 8 nt (RNAP RNA:DNA hybrid; Nudler 2012 Cell 149:1438)
//
// Both modes are handled by the same RibosomeElongation class (same interface,
// different CodonRateTable) — "same rational for nucleotides than polypeptides".
//
// TransloconInsertion (TM helix lateral gating via Sec61 translocon) is also
// evaluated at each growth step for polypeptide chains, using the Hessa 2007
// empirical ΔG scale.
//
// Hardware acceleration in this file:
//   – AVX-512 / AVX2 / Eigen / OpenMP propagated via TransloconInsertion
//   – RibosomeElongation uses Eigen VectorXd for the tridiagonal ODE
//   – ShannonThermoStack / compute_growth_entropy uses AVX-512 histogramming
#include "NATURaLDualAssembly.h"
#include "RibosomeElongation.h"
#include "TransloconInsertion.h"
#include "NucleationDetector.h"
#include "../gaboom.h"   // vcfunction / cfstr / get_apparent_cf_evalue

#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

namespace natural {

// ─── nucleic acid residue names ───────────────────────────────────────────────
static const char* NUCLEIC_ACID_NAMES[] = {
    // RNA — single-letter codes "A","G","C","U" REMOVED (Bug #4 fix).
    // strncmp prefix-match fires on ALA(A), GLY(G), CYS(C), etc. for every
    // standard amino acid → false NA detection on all protein receptors.
    // 3-letter codes cover RNA unambiguously; DA/DG/DC/DT cover DNA.
    "ADE", "GUA", "CYT", "URA",
    // DNA
    "DA", "DG", "DC", "DT", "DA3", "DG3", "DC3", "DT3",
    "DA5", "DG5", "DC5", "DT5",
    nullptr
};

// ─── is_nucleotide_ligand ────────────────────────────────────────────────────
bool is_nucleotide_ligand(const atom* atoms, int n_lig_atoms) {
    if (!atoms || n_lig_atoms <= 0) return false;
    for (int i = 0; i < n_lig_atoms; ++i) {
        const char* name = atoms[i].name;
        if (!name) continue;
        // O2' or O2* indicates ribose (RNA/nucleoside)
        if (strstr(name, "O2'") || strstr(name, "O2*") || strstr(name, "O2P"))
            return true;
        // Common nucleotide atom names
        if (strstr(name, "N9") || strstr(name, "N1") || strstr(name, "C5'"))
            return true;
    }
    return false;
}

// ─── is_nucleic_acid_receptor ────────────────────────────────────────────────
bool is_nucleic_acid_receptor(const resid* residues, int n_residues) {
    if (!residues || n_residues <= 0) return false;
    for (int i = 0; i < n_residues; ++i) {
        const char* name = residues[i].name;
        if (!name) continue;
        for (int k = 0; NUCLEIC_ACID_NAMES[k]; ++k) {
            if (strcmp(name, NUCLEIC_ACID_NAMES[k]) == 0)  // Bug #4 fix: exact match, not prefix
                return true;
        }
    }
    return false;
}

// ─── auto_configure ──────────────────────────────────────────────────────────
NATURaLConfig auto_configure(const atom*  atoms,
                              int          n_lig_atoms,
                              const resid* residues,
                              int          n_residues)
{
    NATURaLConfig cfg;

    bool nucl_lig = is_nucleotide_ligand(atoms, n_lig_atoms);
    bool nucl_rec = is_nucleic_acid_receptor(residues, n_residues);

    // Enable for any receptor with residues.
    // Protein receptor → ribosomal elongation rates (Zhao 2011 / Dong 1996).
    // Nucleic acid receptor → RNAP rates (Uptain 1997 / Jonkers 2014).
    // Callers opt out by setting FA->assume_folded (advanced.assume_folded in config).
    if (n_residues > 0) {
        cfg.enabled                 = true;
        cfg.co_translational_growth = true;
        cfg.sugar_pucker_auto       = nucl_lig;
        // Pure nucleic acid receptor → RNAP rates; all other cases → ribosome rates
        cfg.use_ribosome_speed    = !(nucl_rec && !nucl_lig);
        // Mg²⁺ Hill equation only applies to RNA tertiary folding
        cfg.ion_dependent_folding = nucl_rec;

        std::cout << "[NATURaL] "
                  << (nucl_rec ? "Nucleic acid" : "Protein")
                  << " receptor detected → co-"
                  << (cfg.use_ribosome_speed ? "translational (ribosome)"
                                             : "transcriptional (RNAP)")
                  << " DualAssembly ["
                  << (cfg.organism == ribosome::Organism::EcoliK12 ? "E.coli" : "Human")
                  << "]\n";
    }
    return cfg;
}

// ─── human protofibril parallel protocol helpers ────────────────────────────

double InVivoAssemblyTrack::dwell_time_s() const noexcept {
    return (std::isfinite(mean_elongation_rate) && mean_elongation_rate > 1e-12)
        ? 1.0 / mean_elongation_rate
        : 0.0;
}

double InVivoAssemblyTrack::completion_time_s() const noexcept {
    if (chain_length <= 0) return 0.0;
    double init_wait = (std::isfinite(initiation_rate) && initiation_rate > 1e-12)
        ? 1.0 / initiation_rate
        : 0.0;
    return init_wait + chain_length * dwell_time_s();
}

double InVivoAssemblyTrack::first_exposed_time_s() const noexcept {
    if (chain_length <= 0) return 0.0;
    double init_wait = (std::isfinite(initiation_rate) && initiation_rate > 1e-12)
        ? 1.0 / initiation_rate
        : 0.0;
    int first_exposed_unit = static_cast<int>(std::floor(tunnel_length)) + 1;
    first_exposed_unit = std::clamp(first_exposed_unit, 1, chain_length);
    return init_wait + first_exposed_unit * dwell_time_s();
}

std::vector<InVivoAssemblyTrack> make_human_protofibril_tracks(
    int transcript_nt,
    int peptide_aa,
    bool include_reciprocal_controls)
{
    if (transcript_nt < 0 || peptide_aa < 0)
        throw std::invalid_argument("NATURaL protofibril tracks require non-negative chain lengths");

    std::vector<InVivoAssemblyTrack> tracks;
    tracks.reserve(include_reciprocal_controls ? 4 : 2);

    auto add_track = [&](std::string name,
                         GrowthProcess process,
                         DockingRolePolicy role,
                         int chain_length) {
        if (chain_length <= 0) return;

        const bool is_tx = (process == GrowthProcess::Transcription);
        InVivoAssemblyTrack tr;
        tr.name = std::move(name);
        tr.process = process;
        tr.role_policy = role;
        tr.chain_length = chain_length;
        tr.mean_elongation_rate = is_tx
            ? ribosome::MEAN_NT_RATE_HUMAN
            : ribosome::MEAN_EL_RATE_HUMAN;
        tr.initiation_rate = is_tx
            ? ribosome::K_RNAP_INI_DEFAULT
            : ribosome::K_INI_DEFAULT;
        tr.tunnel_length = is_tx
            ? ribosome::RNAP_TUNNEL_NT
            : ribosome::TUNNEL_LENGTH_AA;
        tr.protofibril_is_target = (role == DockingRolePolicy::ProtofibrilAsTarget);
        tr.compartment = is_tx ? "nucleus" : "cytosol/ER";

        // Eukaryotic transcription and translation are not bacterial-style
        // coupled. Keep the RNAP clock in the parallel schedule, but do not
        // claim direct physical protofibril encounter during nuclear transcript
        // synthesis unless a caller explicitly overrides this field.
        tr.direct_encounter_allowed = !is_tx;
        tracks.push_back(std::move(tr));
    };

    add_track("protofibril_target__nascent_rna_ligand",
              GrowthProcess::Transcription,
              DockingRolePolicy::ProtofibrilAsTarget,
              transcript_nt);
    add_track("protofibril_target__nascent_peptide_ligand",
              GrowthProcess::Translation,
              DockingRolePolicy::ProtofibrilAsTarget,
              peptide_aa);

    if (include_reciprocal_controls) {
        add_track("reciprocal_control__rna_target__protofibril_ligand",
                  GrowthProcess::Transcription,
                  DockingRolePolicy::ReciprocalControl,
                  transcript_nt);
        add_track("reciprocal_control__peptide_target__protofibril_ligand",
                  GrowthProcess::Translation,
                  DockingRolePolicy::ReciprocalControl,
                  peptide_aa);
    }

    return tracks;
}

std::vector<ParallelGrowthEvent> build_parallel_growth_schedule(
    const std::vector<InVivoAssemblyTrack>& tracks)
{
    std::vector<ParallelGrowthEvent> events;
    int total_units = 0;
    for (const auto& tr : tracks) {
        if (tr.chain_length < 0)
            throw std::invalid_argument("NATURaL parallel schedule received a negative chain length");
        total_units += tr.chain_length;
    }
    events.reserve(total_units);

    for (int ti = 0; ti < static_cast<int>(tracks.size()); ++ti) {
        const auto& tr = tracks[ti];
        if (tr.chain_length <= 0) continue;
        if (!std::isfinite(tr.mean_elongation_rate) || tr.mean_elongation_rate <= 0.0)
            throw std::invalid_argument("NATURaL parallel schedule requires positive elongation rates");

        double init_wait = (std::isfinite(tr.initiation_rate) && tr.initiation_rate > 1e-12)
            ? 1.0 / tr.initiation_rate
            : 0.0;
        double dwell = 1.0 / tr.mean_elongation_rate;
        int hidden_units = static_cast<int>(std::floor(std::max(0.0, tr.tunnel_length)));

        for (int unit = 1; unit <= tr.chain_length; ++unit) {
            int exposed = std::max(0, unit - hidden_units);
            events.push_back({
                ti,
                unit,
                init_wait + unit * dwell,
                exposed,
                tr.process,
                tr.role_policy
            });
        }
    }

    std::stable_sort(events.begin(), events.end(),
        [](const ParallelGrowthEvent& a, const ParallelGrowthEvent& b) {
            if (a.t_arrival != b.t_arrival) return a.t_arrival < b.t_arrival;
            if (a.track_index != b.track_index) return a.track_index < b.track_index;
            return a.unit_index < b.unit_index;
        });

    return events;
}

// ─── growth_entropy_nats ─────────────────────────────────────────────────────
// Shannon entropy of the CF distribution accumulated along the growth trajectory.
// Returns dimensionless nats. Convert to bits only at reporting boundaries.
// Dispatches to AVX-512 → AVX2 → Eigen → scalar for the histogram fill.
double growth_entropy_nats(const std::vector<double>& cf_trajectory)
{
    if (cf_trajectory.empty()) return 0.0;

    double min_cf = *std::min_element(cf_trajectory.begin(), cf_trajectory.end());
    double max_cf = *std::max_element(cf_trajectory.begin(), cf_trajectory.end());
    if (max_cf - min_cf < 1e-8) return 0.0;

    constexpr int BINS = 32;
    double bw  = (max_cf - min_cf) / BINS + 1e-10;
    double inv_bw = 1.0 / bw;

    std::array<int, BINS> counts{};
    int total = static_cast<int>(cf_trajectory.size());

#if defined(FLEXAIDS_USE_AVX512) && defined(__AVX512F__)
    {
        __m512d v_min  = _mm512_set1_pd(min_cf);
        __m512d v_ibw  = _mm512_set1_pd(inv_bw);
        __m512d v_bmax = _mm512_set1_pd(static_cast<double>(BINS - 1));
        int i = 0, n = total;
        for (; i + 8 <= n; i += 8) {
            __m512d v    = _mm512_loadu_pd(cf_trajectory.data() + i);
            __m512d off  = _mm512_sub_pd(v, v_min);
            __m512d bins = _mm512_mul_pd(off, v_ibw);
            bins = _mm512_max_pd(_mm512_setzero_pd(),
                                 _mm512_min_pd(bins, v_bmax));
            __m256i bi = _mm512_cvttpd_epi32(bins);
            alignas(32) int bvals[8];
            _mm256_store_si256(reinterpret_cast<__m256i*>(bvals), bi);
            for (int j = 0; j < 8; ++j) counts[bvals[j]]++;
        }
        for (; i < n; ++i) {
            int b = static_cast<int>((cf_trajectory[i] - min_cf) * inv_bw);
            b = std::clamp(b, 0, BINS - 1);
            counts[b]++;
        }
    }
#elif defined(__AVX2__)
    {
        __m256d v_min  = _mm256_set1_pd(min_cf);
        __m256d v_ibw  = _mm256_set1_pd(inv_bw);
        int i = 0, n = total;
        for (; i + 4 <= n; i += 4) {
            __m256d v   = _mm256_loadu_pd(cf_trajectory.data() + i);
            __m256d off = _mm256_sub_pd(v, v_min);
            __m256d bv  = _mm256_mul_pd(off, v_ibw);
            alignas(32) double bvals[4];
            _mm256_store_pd(bvals, bv);
            for (int j = 0; j < 4; ++j) {
                int b = static_cast<int>(bvals[j]);
                counts[std::clamp(b, 0, BINS - 1)]++;
            }
        }
        for (; i < n; ++i) {
            int b = static_cast<int>((cf_trajectory[i] - min_cf) * inv_bw);
            counts[std::clamp(b, 0, BINS - 1)]++;
        }
    }
#else
    {
        Eigen::Map<const Eigen::ArrayXd> vals(cf_trajectory.data(), total);
        Eigen::ArrayXd bins = ((vals - min_cf) * inv_bw).floor().cwiseMax(0).cwiseMin(BINS - 1);
        for (int i = 0; i < total; ++i)
            counts[static_cast<int>(bins(i))]++;
    }
#endif

    Eigen::ArrayXd prob(BINS);
    for (int i = 0; i < BINS; ++i) prob(i) = (double)counts[i] / total;
    Eigen::ArrayXd safe_p   = (prob > 1e-15).select(prob, Eigen::ArrayXd::Ones(BINS));
    Eigen::ArrayXd log_p    = (prob > 1e-15).select(safe_p.log(), Eigen::ArrayXd::Zero(BINS));
    return -(prob * log_p).sum();
}

// ─── build elongation model from config ──────────────────────────────────────
// Returns a RibosomeElongation built from the receptor sequence using either
// ribosome (polypeptide) or RNAP (nucleotide) rates depending on cfg.
// The same Zhao 2011 master equation is used in both cases.
static ribosome::RibosomeElongation build_elongation_model(
    const NATURaLConfig& cfg,
    const resid*         residues,
    int                  n_residues)
{
    using namespace ribosome;

    // Build appropriate rate table
    CodonRateTable table =
        cfg.use_ribosome_speed
            ? (cfg.organism == Organism::EcoliK12
                   ? CodonRateTable::build_ecoli()
                   : CodonRateTable::build_human())
            : (cfg.organism == Organism::EcoliK12
                   ? CodonRateTable::build_rnap_ecoli()
                   : CodonRateTable::build_rnap_human());

    // Extract 1-letter sequence and codon identifiers from residue array.
    // For nucleotide chains, residue names map to single-nucleotide "codons".
    std::string seq;
    std::vector<std::string> codons;
    seq.reserve(n_residues);
    codons.reserve(n_residues);

    for (int i = 0; i < n_residues; ++i) {
        const char* rname = residues[i].name;
        if (!rname) { seq += 'X'; codons.emplace_back(); continue; }

        if (!cfg.use_ribosome_speed) {
            // Nucleotide mode: 1-letter nucleotide code for RNAP rate lookup
            char nt = 'N';
            if      (rname[0] == 'A' || strncmp(rname, "ADE", 3) == 0) nt = 'A';
            else if (rname[0] == 'G' || strncmp(rname, "GUA", 3) == 0) nt = 'G';
            else if (rname[0] == 'C' || strncmp(rname, "CYT", 3) == 0) nt = 'C';
            else if (rname[0] == 'U' || strncmp(rname, "URA", 3) == 0) nt = 'U';
            else if (rname[0] == 'T' || rname[1] == 'T')               nt = 'T';
            seq  += nt;
            codons.emplace_back(1, nt);  // single-character codon
        } else {
            // Polypeptide mode: convert 3-letter residue name to 1-letter AA code
            static const std::pair<const char*, char> aa_table[] = {
                {"ALA",'A'}, {"ARG",'R'}, {"ASN",'N'}, {"ASP",'D'},
                {"CYS",'C'}, {"GLN",'Q'}, {"GLU",'E'}, {"GLY",'G'},
                {"HIS",'H'}, {"ILE",'I'}, {"LEU",'L'}, {"LYS",'K'},
                {"MET",'M'}, {"PHE",'F'}, {"PRO",'P'}, {"SER",'S'},
                {"THR",'T'}, {"TRP",'W'}, {"TYR",'Y'}, {"VAL",'V'},
                {"SEC",'U'}, {"PYL",'O'}, {"HSD",'H'}, {"HSE",'H'},
                {"HSP",'H'}, {"HIE",'H'}, {"HID",'H'}, {"HIP",'H'},
                {"CYX",'C'}, {"ASH",'D'}, {"GLH",'E'},
            };
            char aa_code = 'X';
            for (const auto& [code3, code1] : aa_table) {
                if (strncmp(rname, code3, 3) == 0) { aa_code = code1; break; }
            }
            seq  += aa_code;
            // Use Most-Frequent Codon (MFC) for per-residue elongation rate.
            // MFC gives the highest-abundance codon per AA for the organism,
            // matching dominant in-vivo tRNA pool (Dong 1996; Ikemura 1985).
            // Falls back to empty string (mean rate) for unknown AAs.
            if (aa_code != 'X') {
                static const auto mfc_ecoli =
                    CodonRateTable::aa_to_mfc(Organism::EcoliK12);
                static const auto mfc_human =
                    CodonRateTable::aa_to_mfc(Organism::HumanHEK293);
                const auto& mfc = (cfg.organism == Organism::EcoliK12)
                                  ? mfc_ecoli : mfc_human;
                auto it = mfc.find(aa_code);
                codons.push_back(it != mfc.end() ? it->second : std::string{});
            } else {
                codons.emplace_back(); // unknown AA → mean rate fallback
            }
        }
    }

    double k_ini = cfg.use_ribosome_speed ? K_INI_DEFAULT  : K_RNAP_INI_DEFAULT;
    double k_ter = cfg.use_ribosome_speed ? K_TERM_DEFAULT  : K_RNAP_TERM_DEFAULT;

    return RibosomeElongation(seq, codons, table, k_ini, k_ter);
}

// ─── extract_sequence ────────────────────────────────────────────────────────
// Builds the 1-letter sequence string from residue array for nucleation-seed
// detection. For RNA receptors, uses A/G/C/U nucleotide codes; for protein,
// uses standard 1-letter amino acid codes.
static std::string extract_sequence(const resid* residues, int n, bool is_rna)
{
    std::string seq;
    seq.reserve(n);
    for (int i = 0; i < n; ++i) {
        const char* rname = residues[i].name;
        if (!rname) { seq += is_rna ? 'N' : 'X'; continue; }
        if (is_rna) {
            char nt = 'N';
            if      (rname[0] == 'A' || strncmp(rname, "ADE", 3) == 0) nt = 'A';
            else if (rname[0] == 'G' || strncmp(rname, "GUA", 3) == 0) nt = 'G';
            else if (rname[0] == 'C' || strncmp(rname, "CYT", 3) == 0) nt = 'C';
            else if (rname[0] == 'U' || strncmp(rname, "URA", 3) == 0) nt = 'U';
            else if (rname[0] == 'T' || rname[1] == 'T')               nt = 'U'; // DNA T → U for RNA fold
            seq += nt;
        } else {
            static const std::pair<const char*, char> aa_table[] = {
                {"ALA",'A'}, {"ARG",'R'}, {"ASN",'N'}, {"ASP",'D'},
                {"CYS",'C'}, {"GLN",'Q'}, {"GLU",'E'}, {"GLY",'G'},
                {"HIS",'H'}, {"ILE",'I'}, {"LEU",'L'}, {"LYS",'K'},
                {"MET",'M'}, {"PHE",'F'}, {"PRO",'P'}, {"SER",'S'},
                {"THR",'T'}, {"TRP",'W'}, {"TYR",'Y'}, {"VAL",'V'},
                {"HSD",'H'}, {"HSE",'H'}, {"HSP",'H'}, {"HIE",'H'},
                {"HID",'H'}, {"HIP",'H'}, {"CYX",'C'}, {"ASH",'D'},
                {"GLH",'E'}, {"SEC",'U'}, {"PYL",'O'},
            };
            char aa = 'X';
            for (const auto& [code3, code1] : aa_table) {
                if (strncmp(rname, code3, 3) == 0) { aa = code1; break; }
            }
            seq += aa;
        }
    }
    return seq;
}

// ─── detect_burst_units ──────────────────────────────────────────────────────
// Identifies runs of consecutive fast codons (k_el > burst_threshold × hmean)
// after the ribosomal/RNAP exit tunnel. These are "burst units" where multiple
// monomers are added faster than co-translational folding can compete.
// Reference: Pechmann & Frydman 2013 Nat Struct Mol Biol 20:237.
static std::vector<ribosome::BurstUnit> detect_burst_units(
    const std::vector<double>& k_el,
    double hmean,
    double tunnel_len,
    double burst_threshold)
{
    std::vector<ribosome::BurstUnit> bursts;
    const int N = static_cast<int>(k_el.size());
    const int tunnel_end = static_cast<int>(tunnel_len);
    const double threshold_k = burst_threshold * hmean;

    int i = tunnel_end; // only detect bursts after tunnel exit
    while (i < N) {
        if (k_el[i] > threshold_k) {
            // Start of a potential burst
            int   start     = i;
            double inv_sum  = 0.0;
            while (i < N && k_el[i] > threshold_k) {
                inv_sum += 1.0 / k_el[i];
                ++i;
            }
            int n_res = i - start;
            if (n_res >= 2) { // burst = ≥2 consecutive fast codons
                double total_dwell  = inv_sum;               // Σ 1/k_el
                double hmean_burst  = (inv_sum > 1e-15)
                                      ? n_res / inv_sum : threshold_k;
                bool   fol_pause    = (start > 0 &&
                                       k_el[start - 1] < 0.3 * hmean);
                bursts.push_back({start, i - 1, n_res,
                                   total_dwell, hmean_burst, fol_pause});
            }
        } else {
            ++i;
        }
    }
    return bursts;
}

// ─── DualAssemblyEngine ──────────────────────────────────────────────────────
//
// The full DualAssemblyEngine implementation depends on the FlexAID Vcontacts/
// vcfunction scoring stack. The standalone `dual_assembly` CLI does NOT use the
// in-process engine (it injects synthetic or external GA callbacks instead) and
// therefore defines FLEXAIDS_OMIT_DUAL_ASSEMBLY_ENGINE to skip the heavy
// dependency graph. Other targets (FlexAID, test_natural) leave the macro
// undefined and link the full engine.
#ifndef FLEXAIDS_OMIT_DUAL_ASSEMBLY_ENGINE

DualAssemblyEngine::DualAssemblyEngine(const NATURaLConfig& cfg,
                                         FA_Global* FA, VC_Global* VC,
                                         atom* atoms, resid* residues,
                                         int n_residues)
    : config_(cfg), FA_(FA), VC_(VC),
      atoms_(atoms), residues_(residues), n_residues_(n_residues),
      is_rna_receptor_(is_nucleic_acid_receptor(residues, n_residues))
{}

// Hill equation: k_eff = k_max × [Mg]ⁿ / (K_d ⁿ + [Mg]ⁿ)
// Returns the fractional saturation of Mg²⁺ binding sites (0–1).
double DualAssemblyEngine::mg_hill_factor() const noexcept {
    double mg  = config_.mg_concentration_mM;
    double kd  = ribosome::KD_MG_RNA_MM;
    double n   = ribosome::N_HILL_MG;
    double mgn = std::pow(mg, n);
    double kdn = std::pow(kd, n);
    return mgn / (kdn + mgn);
}

std::vector<DualAssemblyEngine::GrowthStep> DualAssemblyEngine::run() {
    std::vector<GrowthStep> trajectory;

    if (!config_.enabled || !config_.co_translational_growth) {
        return trajectory;
    }

    int max_steps = (config_.max_growth_steps > 0)
                    ? config_.max_growth_steps
                    : n_residues_;
    trajectory.reserve(max_steps);

    // ── Build elongation model (ribosome or RNAP, same master equation) ──────
    auto elong_model = build_elongation_model(config_, residues_, n_residues_);
    const auto& k_el = elong_model.elongation_rates();

    // Tunnel length depends on mode
    const double tunnel_len = config_.use_ribosome_speed
                              ? ribosome::TUNNEL_LENGTH_AA
                              : ribosome::RNAP_TUNNEL_NT;

    // Pause threshold
    const double pause_threshold = config_.use_ribosome_speed
        ? ribosome::RIBOSOME_PAUSE_THRESHOLD
        : ribosome::RNAP_PAUSE_THRESHOLD;

    // Harmonic mean rate for pause detection (filter NaN/Inf/zero values)
    double hmean_rate = elong_model.elongation_rates().empty()
                        ? 16.5
                        : [&](){
                            double inv_sum = 0.0; int cnt = 0;
                            for (double k : k_el) {
                                if (std::isfinite(k) && k > 1e-9) {
                                    inv_sum += 1.0 / k;
                                    ++cnt;
                                }
                            }
                            return (cnt > 0 && inv_sum > 1e-15)
                                   ? static_cast<double>(cnt) / inv_sum
                                   : 16.5;
                          }();

    // Pre-compute analytical arrival times for each residue
    // <T_n> = 1/k_ini + Σ_{i<n} 1/k_i  (Zhao 2011 Eq. 7)
    std::vector<double> t_arrival(max_steps + 1);
    t_arrival[0] = 1.0 / elong_model.k_ini();
    for (int n = 1; n <= max_steps; ++n) {
        double k = (n - 1 < (int)k_el.size()) ? k_el[n - 1]
                                                : elong_model.elongation_rates().back();
        double inv_k = (std::isfinite(k) && k > 1e-9) ? 1.0 / k : 1.0 / hmean_rate;
        t_arrival[n] = t_arrival[n - 1] + inv_k;
    }

    // ── Burst elongation detection ────────────────────────────────────────────
    // Identifies multi-residue fast-codon runs where folding cannot compete.
    auto burst_units = detect_burst_units(
        k_el, hmean_rate, tunnel_len, config_.burst_threshold);

    // Build per-step burst lookup tables
    std::vector<int> step_burst_id  (max_steps, -1);
    std::vector<int> step_burst_size(max_steps,  1);
    std::vector<bool> step_burst_follows_pause(max_steps, false);
    for (int bi = 0; bi < static_cast<int>(burst_units.size()); ++bi) {
        const auto& bu = burst_units[bi];
        for (int s = bu.start_residue;
             s <= bu.end_residue && s < max_steps; ++s) {
            step_burst_id  [s] = bi;
            step_burst_size[s] = bu.n_residues;
            step_burst_follows_pause[s] = bu.follows_pause;
        }
    }

    // ── Nucleation seed detection ─────────────────────────────────────────────
    // Detect RNA hairpin/G-quad or protein helix/hydrophobic-cluster seeds.
    // Each seed applies a folding_rate_boost multiplier to k_fold at that position.
    std::vector<NucleationSeed> nseeds;
    std::vector<double>         seed_boost_map(max_steps, 1.0);
    std::vector<int>            step_seed_id  (max_steps, -1);

    if (config_.detect_nucleation_seeds) {
        std::string seq_str = extract_sequence(residues_, n_residues_, is_rna_receptor_);
        nseeds = NucleationSeedDetector::detect(
            seq_str, is_rna_receptor_, config_.temperature_K);

        // Fill full-length boost map then truncate to max_steps
        auto full_map = NucleationSeedDetector::position_boost_map(
            nseeds, static_cast<int>(seq_str.size()));
        for (int s = 0; s < max_steps && s < static_cast<int>(full_map.size()); ++s)
            seed_boost_map[s] = full_map[s];

        // Build per-step seed-id lookup (first seed wins if seeds overlap)
        for (int si = 0; si < static_cast<int>(nseeds.size()); ++si) {
            for (int s = nseeds[si].start_pos;
                 s <= nseeds[si].end_pos && s < max_steps; ++s) {
                if (step_seed_id[s] < 0) step_seed_id[s] = si;
            }
        }

        if (!nseeds.empty()) {
            std::cout << "[NATURaL] Nucleation seeds detected: "
                      << nseeds.size() << " (";
            int rna_h=0, rna_g=0, pro_hx=0, pro_hp=0;
            for (const auto& s : nseeds) {
                switch (s.type) {
                    case NucleationSeed::Type::RNA_HAIRPIN:          ++rna_h;  break;
                    case NucleationSeed::Type::RNA_GQUADRUPLEX:      ++rna_g;  break;
                    case NucleationSeed::Type::PROTEIN_HELIX:        ++pro_hx; break;
                    case NucleationSeed::Type::PROTEIN_HYDROPHOBIC:  ++pro_hp; break;
                }
            }
            if (rna_h)  std::cout << rna_h  << " RNA hairpin";
            if (rna_g)  std::cout << (rna_h?", ":"") << rna_g  << " G-quad";
            if (pro_hx) std::cout << (rna_h||rna_g?", ":"") << pro_hx << " helix";
            if (pro_hp) std::cout << (rna_h||rna_g||pro_hx?", ":"") << pro_hp << " hydrophobic";
            std::cout << ")\n";
        }
    }

    if (!burst_units.empty()) {
        int max_bs = 0;
        for (const auto& bu : burst_units)
            max_bs = std::max(max_bs, bu.n_residues);
        std::cout << "[NATURaL] Burst elongation units: "
                  << burst_units.size()
                  << " (max " << max_bs << " residues/burst)\n";
    }

    // ── TransloconInsertion (polypeptide chains only) ─────────────────────────
    std::unique_ptr<translocon::TransloconInsertion> tm_model;
    std::string aa_seq_for_tm;
    if (config_.model_tm_insertion && config_.use_ribosome_speed) {
        tm_model = std::make_unique<translocon::TransloconInsertion>(
            config_.temperature_K, 0.5,
            static_cast<int>(ribosome::TUNNEL_LENGTH_AA));
        // Build 1-letter sequence for TransloconInsertion window scanning
        // Uses same 3-letter→1-letter conversion as elongation model
        static const std::pair<const char*, char> aa_table[] = {
            {"ALA",'A'}, {"ARG",'R'}, {"ASN",'N'}, {"ASP",'D'},
            {"CYS",'C'}, {"GLN",'Q'}, {"GLU",'E'}, {"GLY",'G'},
            {"HIS",'H'}, {"ILE",'I'}, {"LEU",'L'}, {"LYS",'K'},
            {"MET",'M'}, {"PHE",'F'}, {"PRO",'P'}, {"SER",'S'},
            {"THR",'T'}, {"TRP",'W'}, {"TYR",'Y'}, {"VAL",'V'},
        };
        aa_seq_for_tm.reserve(n_residues_);
        for (int r = 0; r < n_residues_; ++r) {
            const char* rname = residues_[r].name;
            char aa_code = 'X';
            if (rname) {
                for (const auto& [code3, code1] : aa_table) {
                    if (strncmp(rname, code3, 3) == 0) { aa_code = code1; break; }
                }
            }
            aa_seq_for_tm += aa_code;
        }
    }

    // ── Folding windows from master equation ─────────────────────────────────
    auto fw = elong_model.folding_windows(ribosome::K_FOLD_DEFAULT);
    // Build a quick lookup: residue_idx → folding window info
    // (fw is already sorted by residue index; use binary search or map)

    // ── Main growth loop ──────────────────────────────────────────────────────
    std::vector<double> cf_trajectory;
    cf_trajectory.reserve(max_steps);
    double cumulative_dG = 0.0;
    double weighted_dG_sum = 0.0;
    double dwell_sum = 0.0;

    for (int step = 0; step < max_steps; ++step) {
        // Elongation kinetics at this residue
        double k_n = (step < (int)k_el.size()) ? k_el[step] : hmean_rate;
        double dwell_time  = (std::isfinite(k_n) && k_n > 1e-9)
                             ? 1.0 / k_n
                             : 1.0 / hmean_rate;
        double t_arr       = (step < (int)t_arrival.size()) ? t_arrival[step] : 0.0;
        bool   in_tunnel   = (step < static_cast<int>(tunnel_len));
        bool   is_pause    = (!in_tunnel) && (k_n < pause_threshold * hmean_rate);

        // Co-translational folding probability: P_fold = k_fold / (k_fold + k_el)
        //
        // RNA receptors use differentiated rates:
        //   • Secondary structure (stems/hairpins) folds in microseconds — k ~ 1e4 s⁻¹
        //     → P_fold ≈ 1 during continuous elongation at 25 nt/s
        //   • Tertiary / active-site conformation is Mg²⁺-dependent (Hill equation)
        //     → k_eff = K_FOLD_RNA_TERTIARY × [Mg]ⁿ/(K_d ⁿ + [Mg]ⁿ), only at pause sites
        // Refs: Woodside 2006 PNAS (hairpin); Penedo 2004 RNA; Martick & Scott 2006 Cell.
        double k_fold_here;
        if (is_rna_receptor_ && config_.ion_dependent_folding) {
            if (is_pause) {
                // Pause site: tertiary folding window, Mg²⁺-gated
                k_fold_here = config_.k_fold_rna_tertiary * mg_hill_factor();
            } else {
                // Continuous elongation: secondary structure (fast)
                k_fold_here = config_.k_fold_rna_secondary;
            }
        } else {
            k_fold_here = ribosome::K_FOLD_DEFAULT;
            if (is_pause) k_fold_here *= 3.0;  // Pechmann 2013 pause-site boost
        }
        // Apply nucleation seed boost: detected hairpin/G-quad/helix/hydrophobic
        // seeds increase local k_fold at their sequence positions.
        {
            double boost = seed_boost_map[step];
            if (boost > 1.0 + 1e-9) k_fold_here *= boost;
        }
        double p_cotrans = k_fold_here / (k_fold_here + k_n);

        // CF score for partial complex
        double cf = compute_partial_cf(step + 1);
        cf_trajectory.push_back(cf);

        // Shannon entropy over the growing CF ensemble, dimensionless nats.
        double S_growth = compute_growth_entropy(cf_trajectory);

        // Instantaneous ΔG estimate: ΔH from CF + (−T·ΔS) entropy term.
        // kT [kcal/mol] × H [nats] stays in kcal/mol.
        const double kT = statmech::kB_kcal * config_.temperature_K;
        double step_dG = cf - kT * S_growth;

        // ── TM translocon insertion check ─────────────────────────────────────
        bool   tm_inserted   = false;
        double tm_insert_dG  = 0.0;
        if (tm_model && !in_tunnel && step >= translocon::TM_WINDOW_LEN - 1) {
            // Check if the most recently emerged window is a TM segment
            int win_start = std::max(0, step - translocon::TM_WINDOW_LEN + 1);
            auto win = tm_model->check_window(aa_seq_for_tm, win_start);
            tm_inserted  = win.is_inserted;
            tm_insert_dG = win.deltaG_insert;
            if (tm_inserted) {
                // Include translocon ΔG in the instantaneous kcal/mol score.
                step_dG += tm_insert_dG;
            }
        }

        // Time-weighted mean over growth intermediates. The accumulator is
        // kcal·s/mol internally, but the reported correction remains kcal/mol.
        weighted_dG_sum += step_dG * dwell_time;
        dwell_sum += dwell_time;
        cumulative_dG = (dwell_sum > 0.0) ? (weighted_dG_sum / dwell_sum) : step_dG;

        trajectory.push_back({
            step,
            t_arr,
            k_n,
            dwell_time,
            is_pause,
            in_tunnel,
            cf,
            S_growth,
            p_cotrans,
            cumulative_dG,
            tm_inserted,
            tm_insert_dG,
            // burst annotation
            step_burst_id[step],
            step_burst_size[step],
            step_burst_follows_pause[step],
            // nucleation seed annotation
            step_seed_id[step],
            (step_seed_id[step] >= 0) ? seed_boost_map[step] : 1.0
        });
    }

    final_deltaG_ = cumulative_dG;

    // Optional: validate master equation against analytic solution
    if (n_residues_ >= 20) {
        auto vr = ribosome::validate_master_equation(
            std::min(n_residues_, 60),
            config_.organism);
        if (!vr.passed) {
            std::cerr << "[NATURaL] WARNING: " << vr.message << "\n";
        }
    }

    return trajectory;
}

// ─── compute_partial_cf ───────────────────────────────────────────────────────
// Evaluates the Contact Function for the partially-grown complex using the
// FlexAID energy_matrix.  For each receptor atom in residues [0, n_grown)
// that is within contact distance of a ligand atom, we accumulate the
// complementarity energy weighted by approximate contact area (sphere-point
// sampling) — the same physics as cffunction() but restricted to the grown
// subset of the receptor.
double DualAssemblyEngine::compute_partial_cf(int n_grown_residues) const {
    if (!FA_ || !VC_ || !atoms_ || !residues_) return 0.0;

    // Temporarily mark grown residues as scorable and the rest as non-scorable
    // by adjusting FA_->num_optres to cover only the grown portion.
    // vcfunction() scores the current Cartesian atom coordinates directly —
    // no IC → CC rebuild needed since DualAssembly keeps atoms_ in Cartesian.
    const int saved_num_optres = FA_->num_optres;
    const int n_score = std::min(n_grown_residues, n_residues_);

    // Limit scoring to the first n_score residues (those grown so far).
    // optres[] is ordered by residue index; trim the active count.
    FA_->num_optres = std::min(n_score, saved_num_optres);

    std::vector<std::pair<int,int>> intraclashes;
    bool error = false;
    double cf_val = vcfunction(FA_, VC_, atoms_, residues_, intraclashes, &error);

    // Restore original optres count
    FA_->num_optres = saved_num_optres;

    if (error) return 0.0;
    return cf_val; // kcal/mol (attractive values are negative)
}

// ─── compute_growth_entropy ───────────────────────────────────────────────────
// Shannon entropy of the CF distribution accumulated along the growth trajectory.
// Returns dimensionless nats. Convert to bits only at reporting boundaries.
// Dispatches to AVX-512 → AVX2 → Eigen → scalar for the histogram fill.
double DualAssemblyEngine::compute_growth_entropy(
    const std::vector<double>& cf_trajectory) const
{
    return growth_entropy_nats(cf_trajectory);
}

#endif // FLEXAIDS_OMIT_DUAL_ASSEMBLY_ENGINE

} // namespace natural
