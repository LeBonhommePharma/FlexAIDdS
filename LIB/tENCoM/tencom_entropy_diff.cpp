// tencom_entropy_diff.cpp — TENCoM Vibrational Entropy Differential Tool
//
// Standalone executable for computing vibrational entropy differentials
// between a reference structure and one or more target structures using
// the Torsional Elastic Network Contact Model (TENCoM).
//
// Usage:
//   tencom_entropy_diff --ref reference.pdb target1.pdb [target2.pdb ...]
//   tencom_entropy_diff --ref reference.pdb --targets list.txt
//
// For each target structure:
//   1. Build TENCoM on reference → eigenvalues/eigenvectors (done once)
//   2. Build TENCoM on target   → eigenvalues/eigenvectors
//   3. Compute eigenvalue differentials:  Δλ_m = λ_m(target) − λ_m(ref)
//   4. Compute eigenvector overlaps:     O_m = |v_m(ref) · v_m(target)|
//   5. Compute vibrational entropy differential: ΔS_vib
//   6. Output FlexibilityMode report (global, ligand-free BindingMode analog)
//
// Full flexibility is always ON by default. No ligand required.
//
// Reference: Delarue & Sanejouand, J. Mol. Biol. (2002)
//            Yang, Song & Cui, Biophys. J. (2009)
//            Frappier et al., Proteins 83(11):2073-82 (2015)

#include "tencm.h"
#include "../statmech.h"
#include "../encom.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

// ─── PDB Cα parser (minimal, standalone) ─────────────────────────────────────

struct CaResidue {
    int    resnum;
    char   chain;
    char   resname[4];
    float  x, y, z;
};

static std::vector<CaResidue> parse_pdb_ca(const std::string& pdb_path)
{
    std::ifstream ifs(pdb_path);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open PDB file: " + pdb_path);

    std::vector<CaResidue> cas;
    std::string line;

    while (std::getline(ifs, line)) {
        if (line.size() < 54) continue;

        // Accept ATOM records only (not HETATM)
        if (line.substr(0, 6) != "ATOM  ") continue;

        // Atom name in columns 13-16 (0-indexed: 12..15)
        std::string atom_name = line.substr(12, 4);
        // Match " CA " (standard PDB Cα naming)
        if (atom_name != " CA ") continue;

        CaResidue ca{};
        ca.chain = line[21];

        // Residue name columns 18-20
        std::strncpy(ca.resname, line.substr(17, 3).c_str(), 3);
        ca.resname[3] = '\0';

        // Residue number columns 23-26
        ca.resnum = std::stoi(line.substr(22, 4));

        // Coordinates: x(31-38), y(39-46), z(47-54)
        ca.x = std::stof(line.substr(30, 8));
        ca.y = std::stof(line.substr(38, 8));
        ca.z = std::stof(line.substr(46, 8));

        cas.push_back(ca);
    }

    return cas;
}

static std::vector<std::array<float,3>> ca_to_coords(const std::vector<CaResidue>& cas)
{
    std::vector<std::array<float,3>> coords;
    coords.reserve(cas.size());
    for (const auto& ca : cas)
        coords.push_back({ca.x, ca.y, ca.z});
    return coords;
}

// ─── FlexibilityMode: ligand-free BindingMode analog ─────────────────────────

struct FlexibilityMode {
    std::string label;                // structure identifier (filename)
    int    n_residues    = 0;
    int    n_modes       = 0;         // total normal modes
    int    n_matched     = 0;         // modes compared (min of ref/target)

    // Per-mode differentials (sorted by ascending eigenvalue)
    std::vector<double> ref_eigenvalues;
    std::vector<double> tgt_eigenvalues;
    std::vector<double> delta_eigenvalues;    // Δλ = λ_tgt − λ_ref
    std::vector<double> eigenvector_overlaps; // |v_ref · v_tgt|

    // Per-residue B-factor differentials
    std::vector<float> ref_bfactors;
    std::vector<float> tgt_bfactors;
    std::vector<float> delta_bfactors;        // ΔB = B_tgt − B_ref

    // Global thermodynamics
    double ref_S_vib    = 0.0;   // kcal mol⁻¹ K⁻¹
    double tgt_S_vib    = 0.0;
    double delta_S_vib  = 0.0;   // ΔS_vib = S_tgt − S_ref

    double ref_F_vib    = 0.0;   // −T·S_vib (kcal/mol)
    double tgt_F_vib    = 0.0;
    double delta_F_vib  = 0.0;   // ΔF_vib

    double temperature  = 300.0; // K
};

// ─── Compute vibrational entropy from TENCoM eigenvalues ─────────────────────

static encom::VibrationalEntropy tencom_vibrational_entropy(
    const std::vector<tencm::NormalMode>& modes,
    double temperature_K,
    int skip_rigid = 6)
{
    // Convert tencm::NormalMode → encom::NormalMode for entropy computation
    std::vector<encom::NormalMode> encom_modes;
    for (int m = skip_rigid; m < static_cast<int>(modes.size()); ++m) {
        if (modes[m].eigenvalue < 1e-8) continue;
        encom::NormalMode em;
        em.index      = m + 1;
        em.eigenvalue  = modes[m].eigenvalue;
        em.frequency   = std::sqrt(std::abs(em.eigenvalue));
        encom_modes.push_back(em);
    }

    if (encom_modes.empty()) {
        encom::VibrationalEntropy zero{};
        zero.temperature = temperature_K;
        return zero;
    }

    return encom::ENCoMEngine::compute_vibrational_entropy(encom_modes, temperature_K);
}

// ─── Compute eigenvector overlap (dot product) ──────────────────────────────

static double eigenvector_overlap(const std::vector<double>& v1,
                                   const std::vector<double>& v2)
{
    if (v1.size() != v2.size() || v1.empty()) return 0.0;
    double dot = 0.0, n1 = 0.0, n2 = 0.0;
    for (std::size_t i = 0; i < v1.size(); ++i) {
        dot += v1[i] * v2[i];
        n1  += v1[i] * v1[i];
        n2  += v2[i] * v2[i];
    }
    double denom = std::sqrt(n1) * std::sqrt(n2);
    return denom > 1e-15 ? std::abs(dot / denom) : 0.0;
}

// ─── Build FlexibilityMode from reference + target ──────────────────────────

static FlexibilityMode compute_flexibility_mode(
    const tencm::TorsionalENM& ref_enm,
    const tencm::TorsionalENM& tgt_enm,
    const std::string& tgt_label,
    double temperature)
{
    FlexibilityMode fm;
    fm.label       = tgt_label;
    fm.temperature = temperature;
    fm.n_residues  = tgt_enm.n_residues();

    const auto& ref_modes = ref_enm.modes();
    const auto& tgt_modes = tgt_enm.modes();
    fm.n_modes   = static_cast<int>(tgt_modes.size());
    fm.n_matched = std::min(static_cast<int>(ref_modes.size()),
                            static_cast<int>(tgt_modes.size()));

    // Eigenvalue differentials and eigenvector overlaps
    const int SKIP = 6; // skip rigid-body modes
    for (int m = SKIP; m < fm.n_matched; ++m) {
        fm.ref_eigenvalues.push_back(ref_modes[m].eigenvalue);
        fm.tgt_eigenvalues.push_back(tgt_modes[m].eigenvalue);
        fm.delta_eigenvalues.push_back(tgt_modes[m].eigenvalue - ref_modes[m].eigenvalue);

        // Eigenvector overlap (only meaningful when dimensions match)
        if (ref_modes[m].eigenvector.size() == tgt_modes[m].eigenvector.size()) {
            fm.eigenvector_overlaps.push_back(
                eigenvector_overlap(ref_modes[m].eigenvector, tgt_modes[m].eigenvector));
        } else {
            fm.eigenvector_overlaps.push_back(0.0);
        }
    }

    // B-factor differentials
    fm.ref_bfactors = ref_enm.bfactors(static_cast<float>(temperature));
    fm.tgt_bfactors = tgt_enm.bfactors(static_cast<float>(temperature));

    int n_bf = std::min(static_cast<int>(fm.ref_bfactors.size()),
                        static_cast<int>(fm.tgt_bfactors.size()));
    fm.delta_bfactors.resize(n_bf);
    for (int i = 0; i < n_bf; ++i)
        fm.delta_bfactors[i] = fm.tgt_bfactors[i] - fm.ref_bfactors[i];

    // Vibrational entropy via ENCoM Schlitter formula
    auto ref_vs = tencom_vibrational_entropy(ref_modes, temperature, SKIP);
    auto tgt_vs = tencom_vibrational_entropy(tgt_modes, temperature, SKIP);

    fm.ref_S_vib   = ref_vs.S_vib_kcal_mol_K;
    fm.tgt_S_vib   = tgt_vs.S_vib_kcal_mol_K;
    fm.delta_S_vib = fm.tgt_S_vib - fm.ref_S_vib;

    fm.ref_F_vib   = -temperature * fm.ref_S_vib;
    fm.tgt_F_vib   = -temperature * fm.tgt_S_vib;
    fm.delta_F_vib = fm.tgt_F_vib - fm.ref_F_vib;

    return fm;
}

// ─── Output: FlexibilityMode report (BindingMode-style) ─────────────────────

static void output_flexibility_mode(const FlexibilityMode& fm,
                                     const std::string& ref_label,
                                     std::ostream& os)
{
    os << "\n"
       << "╔══════════════════════════════════════════════════════════════════╗\n"
       << "║  TENCoM Vibrational Entropy Differential — FlexibilityMode     ║\n"
       << "╚══════════════════════════════════════════════════════════════════╝\n\n";

    os << "  Reference : " << ref_label << "\n"
       << "  Target    : " << fm.label  << "\n"
       << "  Residues  : " << fm.n_residues << "\n"
       << "  Modes cmp : " << fm.n_matched << " (excl. 6 rigid-body)\n"
       << "  Temperature: " << std::fixed << std::setprecision(1) << fm.temperature << " K\n"
       << "\n";

    // ── Global thermodynamics ──
    os << "  ┌─────────────────────────────────────────────────────────────┐\n"
       << "  │ Global Vibrational Thermodynamics                          │\n"
       << "  ├─────────────────────────────────────────────────────────────┤\n"
       << std::fixed << std::setprecision(6)
       << "  │  S_vib (ref)     = " << std::setw(14) << fm.ref_S_vib
       << " kcal mol⁻¹ K⁻¹          │\n"
       << "  │  S_vib (target)  = " << std::setw(14) << fm.tgt_S_vib
       << " kcal mol⁻¹ K⁻¹          │\n"
       << "  │  ΔS_vib          = " << std::setw(14) << fm.delta_S_vib
       << " kcal mol⁻¹ K⁻¹          │\n"
       << "  │                                                             │\n"
       << "  │  F_vib (ref)     = " << std::setw(14) << fm.ref_F_vib
       << " kcal mol⁻¹              │\n"
       << "  │  F_vib (target)  = " << std::setw(14) << fm.tgt_F_vib
       << " kcal mol⁻¹              │\n"
       << "  │  ΔF_vib (−TΔS)  = " << std::setw(14) << fm.delta_F_vib
       << " kcal mol⁻¹              │\n"
       << "  └─────────────────────────────────────────────────────────────┘\n\n";

    // ── Top eigenvalue differentials ──
    int n_show = std::min(static_cast<int>(fm.delta_eigenvalues.size()), 20);
    if (n_show > 0) {
        os << "  Eigenvalue Differentials (lowest " << n_show << " non-rigid modes):\n"
           << "  " << std::left
           << std::setw(8)  << "Mode"
           << std::setw(16) << "λ_ref"
           << std::setw(16) << "λ_tgt"
           << std::setw(16) << "Δλ"
           << std::setw(12) << "|v·v|"
           << "\n"
           << "  " << std::string(68, '-') << "\n";

        for (int i = 0; i < n_show; ++i) {
            os << "  " << std::left << std::setw(8) << (i + 7) // mode 7+ (after 6 rigid)
               << std::right << std::scientific << std::setprecision(4)
               << std::setw(14) << fm.ref_eigenvalues[i] << "  "
               << std::setw(14) << fm.tgt_eigenvalues[i] << "  "
               << std::setw(14) << fm.delta_eigenvalues[i] << "  "
               << std::fixed << std::setprecision(4)
               << std::setw(10) << fm.eigenvector_overlaps[i]
               << "\n";
        }
        os << "\n";
    }

    // ── Per-residue B-factor differentials (top movers) ──
    if (!fm.delta_bfactors.empty()) {
        // Find top 10 residues by |ΔB|
        std::vector<int> idx(fm.delta_bfactors.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::partial_sort(idx.begin(),
                          idx.begin() + std::min(10, static_cast<int>(idx.size())),
                          idx.end(),
                          [&](int a, int b) {
                              return std::abs(fm.delta_bfactors[a]) >
                                     std::abs(fm.delta_bfactors[b]);
                          });

        int n_top = std::min(10, static_cast<int>(idx.size()));
        os << "  Per-Residue B-factor Differential (top " << n_top << " by |ΔB|):\n"
           << "  " << std::left
           << std::setw(10) << "Residue"
           << std::setw(14) << "B_ref (Å²)"
           << std::setw(14) << "B_tgt (Å²)"
           << std::setw(14) << "ΔB (Å²)"
           << "\n"
           << "  " << std::string(52, '-') << "\n";

        for (int i = 0; i < n_top; ++i) {
            int r = idx[i];
            os << "  " << std::left << std::setw(10) << (r + 1)
               << std::right << std::fixed << std::setprecision(3)
               << std::setw(12) << fm.ref_bfactors[r] << "  "
               << std::setw(12) << fm.tgt_bfactors[r] << "  "
               << std::setw(12) << fm.delta_bfactors[r]
               << "\n";
        }
        os << "\n";
    }

    // ── Interpretation ──
    os << "  Interpretation:\n";
    if (fm.delta_S_vib > 0.001)
        os << "    → Target is MORE flexible than reference (ΔS_vib > 0)\n";
    else if (fm.delta_S_vib < -0.001)
        os << "    → Target is LESS flexible than reference (ΔS_vib < 0)\n";
    else
        os << "    → Target has similar flexibility to reference (|ΔS_vib| ≈ 0)\n";

    os << "    → ΔF_vib = " << std::fixed << std::setprecision(4) << fm.delta_F_vib
       << " kcal/mol (vibrational free energy penalty/gain)\n\n";
}

// ─── Output: PDB REMARK block (BindingMode-compatible) ──────────────────────

static void output_pdb_remarks(const FlexibilityMode& fm,
                                const std::string& ref_label,
                                const std::string& output_path)
{
    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << output_path << "\n";
        return;
    }

    ofs << "REMARK   TENCoM Vibrational Entropy Differential\n"
        << "REMARK   FlexibilityMode (ligand-free, global flexibility)\n"
        << "REMARK   Reference: " << ref_label << "\n"
        << "REMARK   Target:    " << fm.label << "\n"
        << "REMARK   Residues:  " << fm.n_residues << "\n"
        << "REMARK   Temperature: " << std::fixed << std::setprecision(1)
        << fm.temperature << " K\n"
        << "REMARK\n"
        << "REMARK   THERMODYNAMICS\n"
        << std::fixed << std::setprecision(6)
        << "REMARK   S_vib_ref    = " << fm.ref_S_vib   << " kcal/mol/K\n"
        << "REMARK   S_vib_tgt    = " << fm.tgt_S_vib   << " kcal/mol/K\n"
        << "REMARK   delta_S_vib  = " << fm.delta_S_vib  << " kcal/mol/K\n"
        << "REMARK   F_vib_ref    = " << fm.ref_F_vib    << " kcal/mol\n"
        << "REMARK   F_vib_tgt    = " << fm.tgt_F_vib    << " kcal/mol\n"
        << "REMARK   delta_F_vib  = " << fm.delta_F_vib  << " kcal/mol\n"
        << "REMARK\n";

    // Per-residue ΔB as B-factor column in pseudo-ATOM records
    ofs << "REMARK   Per-residue B-factor differential (ΔB = B_tgt − B_ref)\n";
    int n = static_cast<int>(fm.delta_bfactors.size());
    for (int i = 0; i < n; ++i) {
        // Pseudo-ATOM: Cα with ΔB in the B-factor column (cols 61-66)
        char buf[82];
        std::snprintf(buf, sizeof(buf),
            "ATOM  %5d  CA  UNK A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
            i + 1, i + 1,
            0.0f, 0.0f, 0.0f,  // placeholder coords
            1.00f,
            static_cast<double>(fm.delta_bfactors[i]));
        ofs << buf << "\n";
    }
    ofs << "END\n";

    std::cout << "  → PDB remarks written: " << output_path << "\n";
}

// ─── TSV output for downstream analysis ─────────────────────────────────────

static void output_tsv(const FlexibilityMode& fm,
                        const std::string& ref_label,
                        const std::string& output_path)
{
    std::ofstream ofs(output_path);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << output_path << "\n";
        return;
    }

    // Header with metadata
    ofs << "# TENCoM Entropy Differential\n"
        << "# Reference: " << ref_label << "\n"
        << "# Target: " << fm.label << "\n"
        << "# Temperature: " << fm.temperature << " K\n"
        << "# delta_S_vib: " << fm.delta_S_vib << " kcal/mol/K\n"
        << "# delta_F_vib: " << fm.delta_F_vib << " kcal/mol\n"
        << "#\n";

    // Per-residue table
    ofs << "residue\tB_ref\tB_tgt\tdelta_B\n";
    int n = static_cast<int>(fm.delta_bfactors.size());
    for (int i = 0; i < n; ++i) {
        ofs << (i + 1) << "\t"
            << std::fixed << std::setprecision(4)
            << fm.ref_bfactors[i] << "\t"
            << fm.tgt_bfactors[i] << "\t"
            << fm.delta_bfactors[i] << "\n";
    }

    std::cout << "  → TSV written: " << output_path << "\n";
}

// ─── Usage ──────────────────────────────────────────────────────────────────

static void print_usage(const char* prog)
{
    std::cerr
        << "\nTENCoM Vibrational Entropy Differential Tool\n"
        << "=============================================\n"
        << "Compute vibrational entropy differentials between a reference\n"
        << "structure and one or more target structures. Full flexibility ON.\n\n"
        << "Usage:\n"
        << "  " << prog << " --ref reference.pdb target1.pdb [target2.pdb ...]\n"
        << "  " << prog << " --ref reference.pdb --targets list.txt\n\n"
        << "Options:\n"
        << "  --ref FILE       Reference PDB (required, computed once)\n"
        << "  --targets FILE   Read target PDB paths from file (one per line)\n"
        << "  --temp FLOAT     Temperature in Kelvin (default: 300.0)\n"
        << "  --cutoff FLOAT   Cα contact cutoff in Å (default: 9.0)\n"
        << "  --k0 FLOAT       Spring constant k₀ (default: 1.0)\n"
        << "  --outdir DIR     Output directory for PDB/TSV files (default: .)\n"
        << "  --quiet          Suppress per-mode/per-residue tables\n"
        << "  --help           Show this help\n\n"
        << "Output (per target):\n"
        << "  Console:  FlexibilityMode report with ΔS_vib, ΔF_vib, Δλ, ΔB\n"
        << "  *.pdb:    REMARK block + pseudo-ATOM with ΔB in B-factor column\n"
        << "  *.tsv:    Per-residue B-factor differentials for analysis\n\n";
}

// ─── main ───────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Parse arguments
    std::string ref_path;
    std::string targets_file;
    std::string outdir = ".";
    double temperature = 300.0;
    float  cutoff      = tencm::DEFAULT_RC;
    float  k0          = tencm::DEFAULT_K0;
    bool   quiet       = false;
    std::vector<std::string> target_paths;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--ref" && i + 1 < argc) {
            ref_path = argv[++i];
        } else if (arg == "--targets" && i + 1 < argc) {
            targets_file = argv[++i];
        } else if (arg == "--temp" && i + 1 < argc) {
            temperature = std::stod(argv[++i]);
        } else if (arg == "--cutoff" && i + 1 < argc) {
            cutoff = std::stof(argv[++i]);
        } else if (arg == "--k0" && i + 1 < argc) {
            k0 = std::stof(argv[++i]);
        } else if (arg == "--outdir" && i + 1 < argc) {
            outdir = argv[++i];
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg[0] != '-') {
            target_paths.push_back(arg);
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    // Read target list from file if provided
    if (!targets_file.empty()) {
        std::ifstream tfs(targets_file);
        if (!tfs.is_open()) {
            std::cerr << "Error: cannot open targets file: " << targets_file << "\n";
            return 1;
        }
        std::string line;
        while (std::getline(tfs, line)) {
            // Trim whitespace
            auto start = line.find_first_not_of(" \t\r\n");
            auto end   = line.find_last_not_of(" \t\r\n");
            if (start == std::string::npos) continue;
            line = line.substr(start, end - start + 1);
            if (line.empty() || line[0] == '#') continue;
            target_paths.push_back(line);
        }
    }

    if (ref_path.empty()) {
        std::cerr << "Error: --ref is required.\n";
        print_usage(argv[0]);
        return 1;
    }

    if (target_paths.empty()) {
        std::cerr << "Error: at least one target PDB is required.\n";
        print_usage(argv[0]);
        return 1;
    }

    // ── Step 1: Build TENCoM on reference (done once) ──────────────────────
    std::cout << "\n=== TENCoM Vibrational Entropy Differential ===\n\n";
    std::cout << "Reference: " << ref_path << "\n";
    std::cout << "Targets:   " << target_paths.size() << " structure(s)\n";
    std::cout << "Temperature: " << temperature << " K\n";
    std::cout << "Cutoff: " << cutoff << " Å, k₀: " << k0 << "\n\n";

    std::vector<CaResidue> ref_cas;
    try {
        ref_cas = parse_pdb_ca(ref_path);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing reference: " << e.what() << "\n";
        return 1;
    }

    if (ref_cas.size() < 3) {
        std::cerr << "Error: reference has fewer than 3 Cα atoms (" << ref_cas.size() << ")\n";
        return 1;
    }

    std::cout << "Reference Cα atoms: " << ref_cas.size() << "\n";

    auto ref_coords = ca_to_coords(ref_cas);
    tencm::TorsionalENM ref_enm;
    ref_enm.build_from_ca(ref_coords, cutoff, k0);

    if (!ref_enm.is_built()) {
        std::cerr << "Error: TENCoM build failed on reference.\n";
        return 1;
    }

    // Compute reference vibrational entropy (report once)
    auto ref_vs = tencom_vibrational_entropy(ref_enm.modes(), temperature);
    std::cout << "Reference: " << ref_enm.n_residues() << " residues, "
              << ref_enm.n_bonds() << " torsional DOFs, "
              << ref_enm.modes().size() << " modes\n"
              << "Reference S_vib = " << std::fixed << std::setprecision(6)
              << ref_vs.S_vib_kcal_mol_K << " kcal mol⁻¹ K⁻¹\n\n";

    // ── Step 2: Process each target ────────────────────────────────────────
    //
    // OpenMP parallelism: each target is fully independent (separate TENCoM
    // build + diagonalisation). Console output is serialised via mutex.
    // File I/O per target is also independent.
    std::atomic<int> n_success{0};
    std::atomic<int> n_fail{0};
    const int n_targets = static_cast<int>(target_paths.size());

#ifdef _OPENMP
    std::mutex io_mutex;
    std::cout << "Processing " << n_targets << " targets with "
              << omp_get_max_threads() << " OpenMP threads\n\n";

    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int ti = 0; ti < n_targets; ++ti) {
        const std::string& tgt_path = target_paths[ti];

        std::vector<CaResidue> tgt_cas;
        try {
            tgt_cas = parse_pdb_ca(tgt_path);
        } catch (const std::exception& e) {
#ifdef _OPENMP
            std::lock_guard<std::mutex> lock(io_mutex);
#endif
            std::cerr << "  Error parsing target: " << e.what() << "\n";
            ++n_fail;
            continue;
        }

        if (tgt_cas.size() < 3) {
#ifdef _OPENMP
            std::lock_guard<std::mutex> lock(io_mutex);
#endif
            std::cerr << "  Skipping " << tgt_path << ": fewer than 3 Cα atoms ("
                      << tgt_cas.size() << ")\n";
            ++n_fail;
            continue;
        }

        // Chain compatibility check (warning only)
        if (tgt_cas.size() != ref_cas.size()) {
#ifdef _OPENMP
            std::lock_guard<std::mutex> lock(io_mutex);
#endif
            std::cerr << "  Warning: " << tgt_path << " Cα count mismatch (ref="
                      << ref_cas.size() << ", tgt=" << tgt_cas.size()
                      << "). Proceeding with mode comparison up to min.\n";
        }

        // TENCoM build + diagonalise (thread-safe, no shared state)
        auto tgt_coords = ca_to_coords(tgt_cas);
        tencm::TorsionalENM tgt_enm;
        tgt_enm.build_from_ca(tgt_coords, cutoff, k0);

        if (!tgt_enm.is_built()) {
#ifdef _OPENMP
            std::lock_guard<std::mutex> lock(io_mutex);
#endif
            std::cerr << "  Error: TENCoM build failed on " << tgt_path << "\n";
            ++n_fail;
            continue;
        }

        // Compute FlexibilityMode (independent per target)
        FlexibilityMode fm = compute_flexibility_mode(ref_enm, tgt_enm, tgt_path, temperature);

        // Console output (serialised)
        {
#ifdef _OPENMP
            std::lock_guard<std::mutex> lock(io_mutex);
#endif
            std::cout << "─── Processing: " << tgt_path << " ───\n";
            if (!quiet) {
                output_flexibility_mode(fm, ref_path, std::cout);
            } else {
                std::cout << "  ΔS_vib = " << std::fixed << std::setprecision(6)
                          << fm.delta_S_vib << " kcal/mol/K"
                          << "  ΔF_vib = " << fm.delta_F_vib << " kcal/mol\n";
            }
        }

        // File outputs (independent per target, no lock needed)
        std::string basename = tgt_path;
        auto slash = basename.find_last_of("/\\");
        if (slash != std::string::npos) basename = basename.substr(slash + 1);
        auto dot = basename.find_last_of('.');
        if (dot != std::string::npos) basename = basename.substr(0, dot);

        output_pdb_remarks(fm, ref_path, outdir + "/" + basename + "_tencom_diff.pdb");
        output_tsv(fm, ref_path, outdir + "/" + basename + "_tencom_diff.tsv");

        ++n_success;
    }

    // ── Summary ────────────────────────────────────────────────────────────
    std::cout << "\n=== Summary ===\n"
              << "  Processed: " << n_success.load() << " / " << target_paths.size() << " targets\n";
    if (n_fail.load() > 0)
        std::cout << "  Failed:    " << n_fail.load() << "\n";
    std::cout << "\n";

    return n_fail.load() > 0 ? 1 : 0;
}
