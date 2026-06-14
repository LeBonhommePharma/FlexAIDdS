// main.cpp — CLI entry point for flexaidds_process_ligand
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Usage:
//   flexaidds_process_ligand [options] <input>
//
// Options:
//   -i, --input <file>       Input file (SDF, MOL2) or SMILES string
//   -f, --format <fmt>       Input format: smiles|sdf|mol2 (default: auto)
//   -o, --output <prefix>    Output prefix for .inp and .ga files
//       --lig-name <name>    3-char residue name for PDB records (default: LIG)
//       --sybyl              Print SYBYL atom-type assignments to stdout
//       --type256            Print 256-type encoding to stdout
//       --rings              Print detected rings to stdout
//       --rotatable          Print rotatable bonds to stdout
//       --validate-only      Validate input without generating output files
//       --strict-valence     Treat valence warnings as errors
//       --allow-macrocycles  Bypass macrocycle size guard
//       --allow-peptides     Bypass peptide backbone guard
//   -v, --verbose            Verbose diagnostic output to stderr
//   -h, --help               Show this help message

#include "ProcessLigand.h"
#include "SybylTyper.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace bonmol;

// ---------------------------------------------------------------------------
// Help text
// ---------------------------------------------------------------------------

static void print_help(const char* prog) {
    std::cout <<
R"(flexaidds_process_ligand — FlexAIDdS Phase 3 ligand preprocessor
Copyright 2026 Le Bonhomme Pharma  |  Apache-2.0

Usage:
  )" << prog << R"( [options] <input>

Input:
  -i, --input <file|smiles>  Input file (SDF, MOL2) or bare SMILES string
  -f, --format <fmt>         Input format: smiles | sdf | mol2
                             (default: auto-detect from file extension)

Output:
  -o, --output <prefix>      Write <prefix>.inp and <prefix>.ga
      --lig-name <name>      3-char PDB residue name (default: LIG)

Reporting:
      --sybyl                Print SYBYL atom-type assignments
      --type256              Print 256-type encoded atom types
      --rings                Print detected SSSR rings
      --rotatable            Print rotatable bond list

Validation:
      --validate-only        Validate and report; skip file output
      --strict-valence       Treat valence warnings as errors
      --allow-macrocycles    Bypass ring-size > 12 rejection
      --allow-peptides       Bypass peptide backbone rejection

Diagnostics:
  -v, --verbose              Verbose pipeline stage output
  -h, --help                 Show this help

Examples:
  # Validate aspirin from SMILES
  flexaidds_process_ligand --validate-only "CC(=O)Oc1ccccc1C(=O)O"

  # Process SDF file, write LIG.inp and LIG.ga
  flexaidds_process_ligand -i ligand.sdf -o LIG --lig-name LIG

  # Process MOL2, print SYBYL types and rotatable bonds
  flexaidds_process_ligand -i ligand.mol2 --sybyl --rotatable -o LIG
)";
}

// ---------------------------------------------------------------------------
// Argument parser
// ---------------------------------------------------------------------------

struct CliArgs {
    std::string   input;
    std::string   format_str  = "auto";
    std::string   output;
    std::string   lig_name    = "LIG";
    bool          show_sybyl  = false;
    bool          show_type256 = false;
    bool          show_rings  = false;
    bool          show_rotatable = false;
    bool          validate_only = false;
    bool          strict_valence = false;
    bool          allow_macrocycles = false;
    bool          allow_peptides = false;
    bool          verbose     = false;
    bool          help        = false;
    bool          ok          = true;
    std::string   error;
};

static CliArgs parse_args(int argc, char** argv) {
    CliArgs args;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                args.ok    = false;
                args.error = "option " + a + " requires an argument";
                return "";
            }
            return argv[++i];
        };

        if (a == "-h" || a == "--help")            { args.help = true; return args; }
        else if (a == "-v" || a == "--verbose")    args.verbose = true;
        else if (a == "-i" || a == "--input")      args.input  = next();
        else if (a == "-f" || a == "--format")     args.format_str = next();
        else if (a == "-o" || a == "--output")     args.output = next();
        else if (a == "--lig-name")                args.lig_name = next();
        else if (a == "--sybyl")                   args.show_sybyl = true;
        else if (a == "--type256")                 args.show_type256 = true;
        else if (a == "--rings")                   args.show_rings  = true;
        else if (a == "--rotatable")               args.show_rotatable = true;
        else if (a == "--validate-only")           args.validate_only = true;
        else if (a == "--strict-valence")          args.strict_valence = true;
        else if (a == "--allow-macrocycles")       args.allow_macrocycles = true;
        else if (a == "--allow-peptides")          args.allow_peptides = true;
        else if (a[0] == '-') {
            args.ok    = false;
            args.error = "unknown option: " + a;
            return args;
        } else {
            positional.push_back(a);
        }

        if (!args.ok) return args;
    }

    // Positional arg: input file or SMILES
    if (args.input.empty() && !positional.empty()) {
        args.input = positional[0];
    }

    if (args.input.empty()) {
        args.ok    = false;
        args.error = "no input specified (use -i or pass as positional argument)";
    }

    return args;
}

// ---------------------------------------------------------------------------
// Format string → InputFormat
// ---------------------------------------------------------------------------

static InputFormat parse_format(const std::string& fmt) {
    std::string lower = fmt;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    if (lower == "smiles" || lower == "smi") return InputFormat::SMILES;
    if (lower == "sdf" || lower == "mol")    return InputFormat::SDF;
    if (lower == "mol2")                     return InputFormat::MOL2;
    return InputFormat::AUTO;
}

// ---------------------------------------------------------------------------
// Print ring information
// ---------------------------------------------------------------------------

static void print_rings(const BonMol& mol) {
    std::cout << "\nDetected rings (" << mol.rings.size() << " total):\n";
    for (size_t i = 0; i < mol.rings.size(); ++i) {
        const Ring& r = mol.rings[i];
        std::cout << "  Ring " << (i + 1) << " (size=" << r.size
                  << ", aromatic=" << (r.is_aromatic ? "yes" : "no") << "): atoms [";
        for (size_t k = 0; k < r.atom_indices.size(); ++k) {
            if (k) std::cout << ",";
            std::cout << (r.atom_indices[k] + 1);
        }
        std::cout << "]\n";
    }
}

// ---------------------------------------------------------------------------
// Print SYBYL types
// ---------------------------------------------------------------------------

static void print_sybyl(const BonMol& mol) {
    std::cout << "\nSYBYL atom-type assignments (" << mol.num_atoms() << " atoms):\n";
    std::cout << std::setw(6) << "Idx" << " "
              << std::setw(5) << "Elem" << " "
              << std::setw(8) << "SYBYL" << " "
              << std::setw(8) << "Type"  << " "
              << std::setw(12) << "Charge" << " "
              << "HBD HBA Arom\n";
    std::cout << std::string(60, '-') << "\n";
    for (int i = 0; i < mol.num_atoms(); ++i) {
        const Atom& a = mol.atoms[i];
        std::cout << std::setw(6) << (i + 1) << " "
                  << std::setw(5) << static_cast<int>(a.element) << " "
                  << std::setw(8) << sybyl::sybyl_type_name(a.sybyl_type) << " "
                  << std::setw(8) << a.sybyl_type << " "
                  << std::setw(12) << std::fixed << std::setprecision(4) << a.partial_charge << " "
                  << std::setw(3) << (a.is_hbond_donor    ? "Y" : "N") << " "
                  << std::setw(3) << (a.is_hbond_acceptor ? "Y" : "N") << " "
                  << std::setw(3) << (a.is_aromatic        ? "Y" : "N") << "\n";
    }
}

// ---------------------------------------------------------------------------
// Print 256-type encoding
// ---------------------------------------------------------------------------

static void print_type256(const BonMol& mol) {
    std::cout << "\n256-type atom encoding (" << mol.num_atoms() << " atoms):\n";
    std::cout << std::setw(6) << "Idx" << " "
              << std::setw(8) << "Type256" << " "
              << "  [Base|Donor|Acceptor]\n";
    std::cout << std::string(40, '-') << "\n";
    for (int i = 0; i < mol.num_atoms(); ++i) {
        const Atom& a = mol.atoms[i];
        uint8_t t    = a.type_256;
        uint8_t base = t & 0x3Fu;
        uint8_t donor = (t >> 7) & 0x01u;
        uint8_t acc   = (t >> 6) & 0x01u;
        std::cout << std::setw(6) << (i + 1) << " "
                  << std::setw(8) << static_cast<int>(t)
                  << "  [" << static_cast<int>(base)
                  << "|" << static_cast<int>(donor)
                  << "|" << static_cast<int>(acc) << "]\n";
    }
}

// ---------------------------------------------------------------------------
// Print rotatable bonds
// ---------------------------------------------------------------------------

static void print_rotatable(const BonMol& mol) {
    std::cout << "\nRotatable bonds:\n";
    int count = 0;
    for (const Bond& b : mol.bonds) {
        if (!b.is_rotatable) continue;
        ++count;
        std::cout << "  Bond " << (b.atom_i + 1) << " -- " << (b.atom_j + 1)
                  << "  (" 
                  << static_cast<int>(mol.atoms[b.atom_i].element)
                  << " -- "
                  << static_cast<int>(mol.atoms[b.atom_j].element)
                  << ")\n";
    }
    if (count == 0) std::cout << "  (none)\n";
}

// ---------------------------------------------------------------------------
// Print pipeline summary
// ---------------------------------------------------------------------------

static void print_summary(const ProcessResult& result) {
    std::cout << "\n=== FlexAIDdS ProcessLigand Summary ===\n";
    std::cout << "  Atoms:            " << result.num_atoms        << "\n";
    std::cout << "  Heavy atoms:      " << result.num_heavy_atoms  << "\n";
    std::cout << "  Mol. weight:      " << std::fixed << std::setprecision(3)
              << result.molecular_weight                           << " Da\n";
    std::cout << "  Rings:            " << result.num_rings        << "\n";
    std::cout << "  Aromatic rings:   " << result.num_arom_rings   << "\n";
    std::cout << "  Rotatable bonds:  " << result.num_rot_bonds    << "\n";
    std::cout << "  Dihedral genes:   " << result.writer_result.num_dihedral_genes << "\n";
    std::cout << "\nPipeline stages:\n";
    for (const auto& s : result.stage_results) {
        std::cout << "  [" << (s.ok ? "OK" : "FAIL") << "] "
                  << std::setw(20) << std::left << s.stage << " " << s.message << "\n";
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    CliArgs args = parse_args(argc, argv);

    if (args.help) {
        print_help(argv[0]);
        return 0;
    }

    if (!args.ok) {
        std::cerr << "Error: " << args.error << "\n";
        std::cerr << "Run with -h for help.\n";
        return 1;
    }

    // Build pipeline options
    ProcessOptions opts;
    opts.input            = args.input;
    opts.format           = parse_format(args.format_str);
    opts.output_prefix    = args.output;
    opts.lig_name         = args.lig_name;
    opts.validate_only    = args.validate_only;
    opts.strict_valence   = args.strict_valence;
    opts.allow_macrocycles = args.allow_macrocycles;
    opts.allow_peptides   = args.allow_peptides;
    opts.verbose          = args.verbose;
    opts.write_inp        = !args.output.empty() && !args.validate_only;
    opts.write_ga         = !args.output.empty() && !args.validate_only;

    // Run pipeline
    ProcessResult result = ProcessLigand::process(opts);

    if (!result.success) {
        std::cerr << "ERROR: " << result.error << "\n";
        // Print any stage-level diagnostics
        for (const auto& s : result.stage_results) {
            if (!s.ok)
                std::cerr << "  Stage [" << s.stage << "]: " << s.message << "\n";
        }
        return 2;
    }

    // Print requested reports
    if (args.show_rings)      print_rings(result.mol);
    if (args.show_sybyl)      print_sybyl(result.mol);
    if (args.show_type256)    print_type256(result.mol);
    if (args.show_rotatable)  print_rotatable(result.mol);

    // Always print summary
    print_summary(result);

    // If output was written, say so
    if (!args.output.empty() && !args.validate_only) {
        std::cout << "\nOutput written:\n";
        std::cout << "  " << args.output << ".inp\n";
        std::cout << "  " << args.output << ".ga\n";
    }

    return 0;
}
