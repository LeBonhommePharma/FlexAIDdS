// tencom_main.cpp — Standalone tENCoM vibrational entropy differential tool
//
// Usage: tENCoM reference.pdb target1.pdb [target2.pdb ...] [options]
//
// Computes TorsionalENM normal modes on a reference structure and one or more
// target structures, outputs eigenvalue/eigenvector differentials and
// vibrational entropy differences (ΔS_vib) in BindingMode-like PDB format.
//
// The reference structure is always required (first positional argument).
// Full flexibility is always on (all torsional modes used).
//
// Options:
//   -T <temperature>   Temperature in Kelvin (default: 300.0)
//   -r <cutoff>        Contact cutoff in Angstroms (default: 9.0)
//   -k <spring_const>  Spring constant k0 (default: 1.0)
//   -o <prefix>        Output file prefix (default: "tencom")
//   -h                 Print help and exit

#include "pdb_calpha.h"
#include "tencom_diff.h"
#include "tencom_output.h"
#include "tencm.h"
#include "../encom.h"

#include "../flexaid_exception.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>

// ─── CLI argument parsing ───────────────────────────────────────────────────

struct Options {
    std::vector<std::string> pdb_files;  // [0] = reference, [1..N] = targets
    double temperature  = 300.0;
    std::vector<double> temperatures;    // filled if -T is a range
    float  cutoff       = 9.0f;
    float  k0           = 1.0f;
    std::string prefix  = "tencom";
    double eigenvalue_cutoff = 1e-6;
    double omega_scale = 1.0;
    bool   has_omega_scale = false;
    std::string calibration_label = "model-scale";
    std::string calibration_provenance =
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied";
    bool   has_calibration_provenance = false;
    encom::FrequencyCalibration frequency_calibration =
        encom::FrequencyCalibration::model_scale();
    bool   output_pdb   = true;
    bool   output_json  = false;
    bool   output_csv   = false;
    int    verbosity    = 1;    // 0=quiet, 1=normal, 2=verbose
};

static void print_usage(const char* progname) {
    std::cout
        << "\nFlexAIDdS tENCoM — Vibrational Entropy Differential Tool\n"
        << "  Torsional Elastic Network Contact Model (TENCoM)\n\n"
        << "Usage: " << progname << " reference.pdb target1.pdb [target2.pdb ...] [options]\n\n"
        << "The reference structure is always required (first positional argument).\n"
        << "Full flexibility is always ON (all torsional normal modes).\n\n"
        << "Options:\n"
        << "  -T <temp>     Temperature in Kelvin, or range start:end:step\n"
        << "                  e.g. -T 300 or -T 200:400:25  (default: 300.0)\n"
        << "  -r <cutoff>   Contact cutoff in Angstroms    (default: 9.0)\n"
        << "  -k <k0>       Spring constant                (default: 1.0)\n"
        << "  -o <prefix>   Output file prefix             (default: tencom)\n"
        << "  -f <format>   Output format: pdb, json, csv, all (default: pdb)\n"
        << "  --omega-scale <s>  Convert sqrt(model eigenvalue) to rad/s\n"
        << "                  by omega = s * sqrt(lambda). Without this, S_vib is\n"
        << "                  reported as model-scale heuristic only.\n"
        << "  --calibration-label <text>  Label for --omega-scale provenance\n"
        << "  --calibration-provenance <text>  Required source note with --omega-scale\n"
        << "  --list <file> Read target PDB paths from file (one per line)\n"
        << "  -v            Verbose output (per-mode eigenvalue lists, timing)\n"
        << "  -q            Quiet output (summary table only)\n"
        << "  -h            Print this help and exit\n\n"
        << "Output:\n"
        << "  <prefix>_mode_N.pdb    PDB files with REMARK thermodynamic metadata\n"
        << "  <prefix>_results.json  Full results in JSON format (-f json|all)\n"
        << "  <prefix>_summary.csv   Summary table as CSV (-f csv|all)\n"
        << "  Summary table          Printed to stdout\n\n";
}

// Safe numeric parser — throws on invalid input
static double parse_double(const char* str, const char* flag) {
    char* endptr = nullptr;
    errno = 0;
    double val = std::strtod(str, &endptr);
    if (errno != 0 || endptr == str || *endptr != '\0') {
        throw FlexAIDException(std::string("Invalid numeric value for ") + flag + ": \"" + str + "\"");
    }
    return val;
}

static Options parse_args(int argc, char* argv[]) {
    Options opts;

    if (argc < 2) {
        print_usage(argv[0]);
        throw FlexAIDException("Insufficient arguments: reference PDB required");
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            throw FlexAIDException("Help requested", 0);
        } else if (arg == "-T") {
            if (i + 1 >= argc) {
                throw FlexAIDException("-T requires a temperature value or range (start:end:step)");
            }
            std::string targ = argv[++i];
            // Check if it's a range: start:end:step
            size_t c1 = targ.find(':');
            if (c1 != std::string::npos) {
                size_t c2 = targ.find(':', c1 + 1);
                if (c2 == std::string::npos) {
                    throw FlexAIDException("Temperature range requires format start:end:step");
                }
                double t_start = parse_double(targ.substr(0, c1).c_str(), "-T start");
                double t_end   = parse_double(targ.substr(c1+1, c2-c1-1).c_str(), "-T end");
                double t_step  = parse_double(targ.substr(c2+1).c_str(), "-T step");
                if (t_start <= 0 || t_end <= 0 || t_step <= 0) {
                    throw FlexAIDException("Temperature range values must be positive");
                }
                if (t_start > t_end) {
                    throw FlexAIDException("Temperature range start must be <= end");
                }
                for (double t = t_start; t <= t_end + t_step * 0.01; t += t_step) {
                    opts.temperatures.push_back(t);
                }
                opts.temperature = t_start;  // default to first
            } else {
                opts.temperature = parse_double(targ.c_str(), "-T");
            }
        } else if (arg == "-r") {
            if (i + 1 >= argc) {
                throw FlexAIDException("-r requires a cutoff value");
            }
            opts.cutoff = static_cast<float>(parse_double(argv[++i], "-r"));
        } else if (arg == "-k") {
            if (i + 1 >= argc) {
                throw FlexAIDException("-k requires a spring constant value");
            }
            opts.k0 = static_cast<float>(parse_double(argv[++i], "-k"));
        } else if (arg == "-o") {
            if (i + 1 >= argc) {
                throw FlexAIDException("-o requires an output prefix");
            }
            opts.prefix = argv[++i];
        } else if (arg == "--omega-scale") {
            if (i + 1 >= argc) {
                throw FlexAIDException("--omega-scale requires a finite positive scale");
            }
            opts.omega_scale = parse_double(argv[++i], "--omega-scale");
            opts.has_omega_scale = true;
        } else if (arg == "--calibration-label") {
            if (i + 1 >= argc) {
                throw FlexAIDException("--calibration-label requires text");
            }
            opts.calibration_label = argv[++i];
        } else if (arg == "--calibration-provenance") {
            if (i + 1 >= argc) {
                throw FlexAIDException("--calibration-provenance requires text");
            }
            opts.calibration_provenance = argv[++i];
            opts.has_calibration_provenance = true;
        } else if (arg == "--list") {
            if (i + 1 >= argc) {
                throw FlexAIDException("--list requires a file path");
            }
            std::string listfile = argv[++i];
            std::ifstream lfs(listfile);
            if (!lfs.is_open()) {
                throw FlexAIDException("Cannot open list file: " + listfile);
            }
            std::string pdb_line;
            while (std::getline(lfs, pdb_line)) {
                // Skip empty lines and comments
                if (pdb_line.empty() || pdb_line[0] == '#') continue;
                // Trim trailing whitespace
                while (!pdb_line.empty() && (pdb_line.back() == ' ' || pdb_line.back() == '\r'))
                    pdb_line.pop_back();
                if (!pdb_line.empty())
                    opts.pdb_files.push_back(pdb_line);
            }
        } else if (arg == "-v" || arg == "--verbose") {
            opts.verbosity = 2;
        } else if (arg == "-q" || arg == "--quiet") {
            opts.verbosity = 0;
        } else if (arg == "-f") {
            if (i + 1 >= argc) {
                throw FlexAIDException("-f requires a format: pdb, json, csv, or all");
            }
            std::string fmt = argv[++i];
            if (fmt == "pdb") {
                opts.output_pdb = true;
            } else if (fmt == "json") {
                opts.output_pdb = false; opts.output_json = true;
            } else if (fmt == "csv") {
                opts.output_pdb = false; opts.output_csv = true;
            } else if (fmt == "all") {
                opts.output_pdb = true; opts.output_json = true; opts.output_csv = true;
            } else {
                throw FlexAIDException("Unknown format '" + fmt + "'. Use pdb, json, csv, or all");
            }
        } else if (arg[0] == '-') {
            print_usage(argv[0]);
            throw FlexAIDException("Unknown option: " + arg);
        } else {
            opts.pdb_files.push_back(arg);
        }
    }

    if (opts.pdb_files.empty()) {
        print_usage(argv[0]);
        throw FlexAIDException("Reference PDB is required");
    }

    // Validate numeric parameters
    if (opts.temperature <= 0.0) {
        throw FlexAIDException("Temperature must be positive (got " + std::to_string(opts.temperature) + " K)");
    }
    if (opts.cutoff <= 0.0f) {
        throw FlexAIDException("Contact cutoff must be positive (got " + std::to_string(opts.cutoff) + " A)");
    }
    if (opts.k0 <= 0.0f) {
        throw FlexAIDException("Spring constant must be positive (got " + std::to_string(opts.k0) + ")");
    }
    if (opts.has_omega_scale) {
        if (!std::isfinite(opts.omega_scale) || opts.omega_scale <= 0.0) {
            throw FlexAIDException("--omega-scale must be finite and positive");
        }
        if (!opts.has_calibration_provenance) {
            throw FlexAIDException(
                "--omega-scale requires --calibration-provenance; a scalar alone is not calibration provenance");
        }
        if (opts.calibration_label == "model-scale") {
            opts.calibration_label = "user-supplied";
        }
        opts.frequency_calibration = encom::FrequencyCalibration::calibrated_scale(
            opts.omega_scale,
            opts.calibration_label,
            opts.calibration_provenance);
    } else {
        opts.frequency_calibration = encom::FrequencyCalibration::model_scale();
    }

    // Verify PDB files exist
    for (const auto& path : opts.pdb_files) {
        std::ifstream test(path);
        if (!test.is_open()) {
            throw FlexAIDException("Cannot open PDB file: " + path);
        }
    }

    return opts;
}

static void apply_calibration_metadata(
    tencom_output::FlexMode& mode,
    const encom::VibrationalEntropy& svib)
{
    mode.eigenvalue_to_omega = svib.eigenvalue_to_omega;
    mode.frequency_calibrated = svib.calibrated;
    mode.calibration_label = svib.calibration_label;
    mode.calibration_provenance = svib.calibration_provenance;
}

// ─── main ───────────────────────────────────────────────────────────────────

// Run analysis at a single temperature, outputting results
static void run_analysis_at_temperature(
    const Options& opts,
    double temperature,
    const tencm::TorsionalENM& ref_enm,
    const std::vector<tencm::TorsionalENM>& tgt_enms,
    std::vector<tencom_pdb::CalphaStructure>& all_structures,
    const std::string& prefix_suffix = "")
{
    std::string prefix = opts.prefix + prefix_suffix;

    // Compute reference vibrational entropy at this temperature
    auto ref_encom_modes = tencom_diff::to_encom_modes(ref_enm.modes());
    auto ref_svib = encom::ENCoMEngine::compute_vibrational_entropy(
        ref_encom_modes, temperature, opts.frequency_calibration,
        opts.eigenvalue_cutoff);

    tencom_output::FlexMode ref_mode;
    ref_mode.mode_id = 0;
    ref_mode.pdb_path = opts.pdb_files[0];
    ref_mode.label = "reference";
    ref_mode.S_vib = ref_svib.S_vib_kcal_mol_K;
    ref_mode.delta_S_vib = 0.0;
    ref_mode.delta_F_vib = 0.0;
    ref_mode.bfactors = ref_enm.bfactors(static_cast<float>(temperature));
    ref_mode.n_modes = ref_svib.n_modes;
    ref_mode.n_residues = ref_enm.n_residues();
    apply_calibration_metadata(ref_mode, ref_svib);

    // Per-residue S_vib decomposition for reference
    {
        double sum_bf = 0.0;
        for (float bf : ref_mode.bfactors) sum_bf += bf;
        ref_mode.per_residue_svib.resize(ref_mode.bfactors.size());
        if (sum_bf > 0.0) {
            for (size_t i = 0; i < ref_mode.bfactors.size(); ++i) {
                ref_mode.per_residue_svib[i] =
                    ref_svib.S_vib_kcal_mol_K * (ref_mode.bfactors[i] / sum_bf);
            }
        }
    }

    tencom_output::FlexPopulation population;
    population.temperature = temperature;
    population.output_prefix = prefix;
    population.modes.push_back(ref_mode);

    for (size_t t = 0; t < tgt_enms.size(); ++t) {
        auto diff = tencom_diff::compute_differential(
            ref_enm, tgt_enms[t],
            opts.pdb_files[0], opts.pdb_files[t + 1],
            temperature, opts.frequency_calibration, opts.eigenvalue_cutoff);

        tencom_output::FlexMode tgt_mode;
        tgt_mode.mode_id = static_cast<int>(t + 1);
        tgt_mode.pdb_path = opts.pdb_files[t + 1];
        tgt_mode.label = opts.pdb_files[t + 1];
        tgt_mode.S_vib = diff.svib_tgt.S_vib_kcal_mol_K;
        tgt_mode.delta_S_vib = diff.delta_S_vib;
        tgt_mode.delta_F_vib = diff.delta_F_vib;
        tgt_mode.mode_data = std::move(diff.mode_comparisons);
        tgt_mode.bfactors = std::move(diff.bfactors_tgt);
        tgt_mode.delta_bfactors = std::move(diff.delta_bfactors);
        tgt_mode.per_residue_svib = std::move(diff.per_residue_svib_tgt);
        tgt_mode.per_residue_delta_svib = std::move(diff.per_residue_delta_svib);
        tgt_mode.n_modes = diff.svib_tgt.n_modes;
        tgt_mode.n_residues = tgt_enms[t].n_residues();
        apply_calibration_metadata(tgt_mode, diff.svib_tgt);

        population.modes.push_back(std::move(tgt_mode));
    }

    population.sort_by_free_energy();
    population.print_summary();

    if (opts.output_pdb) {
        int n = std::min(all_structures.size(), population.modes.size());
        for (int i = 0; i < static_cast<int>(n); ++i) {
            population.write_mode_pdb(population.modes[i], all_structures[i]);
        }
    }
    if (opts.output_json) {
        population.write_json(all_structures);
    }
    if (opts.output_csv) {
        population.write_csv();
    }
}

int main(int argc, char* argv[]) {
  try {
    Options opts = parse_args(argc, argv);

    // If no explicit temperature range, use single temperature
    if (opts.temperatures.empty()) {
        opts.temperatures.push_back(opts.temperature);
    }

    int V = opts.verbosity;  // shorthand

    if (V >= 1) {
        std::cout << "\n=== FlexAIDdS tENCoM — Vibrational Entropy Differential Tool ===\n";
        if (opts.temperatures.size() > 1) {
            std::cout << "Temperature scan: " << opts.temperatures.front()
                      << " - " << opts.temperatures.back() << " K ("
                      << opts.temperatures.size() << " points)\n";
        } else {
            std::cout << "Temperature: " << opts.temperatures[0] << " K\n";
        }
        std::cout << "Contact cutoff: " << opts.cutoff << " A\n";
        std::cout << "Spring constant k0: " << opts.k0 << "\n";
        std::cout << "Frequency calibration: "
                  << (opts.frequency_calibration.calibrated
                      ? "calibrated" : "model-scale heuristic")
                  << " (" << opts.frequency_calibration.label
                  << ", scale=" << opts.frequency_calibration.eigenvalue_to_omega
                  << " rad/s per sqrt(model unit))\n";
        std::cout << "Full flexibility: ON\n\n";
    }

    // ── Step 1: Read reference PDB ──────────────────────────────────────────
    if (V >= 1) std::cout << "Reading reference: " << opts.pdb_files[0] << "\n";
    tencom_pdb::CalphaStructure ref_struct;
    try {
        ref_struct = tencom_pdb::read_pdb_calpha(opts.pdb_files[0]);
    } catch (const std::exception& e) {
        std::cerr << "Error reading reference PDB: " << e.what() << "\n";
        return 1;
    }
    if (V >= 1) std::cout << "  " << ref_struct.res_cnt << " residues (backbone atoms)\n";

    // ── Step 2: Build reference TorsionalENM ────────────────────────────────
    if (V >= 1) std::cout << "Building reference TorsionalENM...\n";
    tencm::TorsionalENM ref_enm;
    ref_enm.build(ref_struct.atoms.data(), ref_struct.residues.data(),
                  ref_struct.res_cnt, opts.cutoff, opts.k0);

    if (!ref_enm.is_built()) {
        std::cerr << "Error: failed to build reference ENM (need >= 3 residues).\n";
        return 1;
    }
    if (V >= 1) std::cout << "  " << ref_enm.modes().size() << " normal modes computed\n";

    if (V >= 2) {
        // Verbose: print eigenvalue list
        std::cout << "  Eigenvalues:";
        const auto& modes = ref_enm.modes();
        int n_show = std::min(static_cast<int>(modes.size()), 20);
        for (int i = 0; i < n_show; ++i) {
            std::cout << " " << std::setprecision(4) << modes[i].eigenvalue;
        }
        if (static_cast<int>(modes.size()) > n_show)
            std::cout << " ... (" << modes.size() - n_show << " more)";
        std::cout << "\n";
    }

    // ── Step 3: Build target TorsionalENMs ──────────────────────────────────
    std::vector<tencm::TorsionalENM> tgt_enms;
    std::vector<tencom_pdb::CalphaStructure> all_structures;
    all_structures.push_back(std::move(ref_struct));

    for (size_t t = 1; t < opts.pdb_files.size(); ++t) {
        if (V >= 1) std::cout << "\nReading target " << t << ": " << opts.pdb_files[t] << "\n";

        tencom_pdb::CalphaStructure tgt_struct;
        try {
            tgt_struct = tencom_pdb::read_pdb_calpha(opts.pdb_files[t]);
        } catch (const std::exception& e) {
            std::cerr << "  Error reading target PDB: " << e.what() << " — skipping.\n";
            continue;
        }
        if (V >= 1) std::cout << "  " << tgt_struct.res_cnt << " residues (backbone atoms)\n";

        tencm::TorsionalENM tgt_enm;
        tgt_enm.build(tgt_struct.atoms.data(), tgt_struct.residues.data(),
                      tgt_struct.res_cnt, opts.cutoff, opts.k0);

        if (!tgt_enm.is_built()) {
            std::cerr << "  Warning: failed to build target ENM — skipping.\n";
            continue;
        }
        if (V >= 1) std::cout << "  " << tgt_enm.modes().size() << " normal modes computed\n";

        tgt_enms.push_back(std::move(tgt_enm));
        all_structures.push_back(std::move(tgt_struct));
    }

    // ── Step 4: Run analysis at each temperature ────────────────────────────
    for (size_t ti = 0; ti < opts.temperatures.size(); ++ti) {
        double T = opts.temperatures[ti];
        std::string suffix = "";
        if (opts.temperatures.size() > 1 && V >= 1) {
            std::cout << "\n──── Temperature: " << T << " K ────\n";
            // Append temperature to output prefix for multi-T runs
            std::ostringstream ss;
            ss << "_T" << static_cast<int>(T);
            suffix = ss.str();
        }

        run_analysis_at_temperature(opts, T, ref_enm, tgt_enms,
                                    all_structures, suffix);
    }

    // ── Temperature scan CSV (multi-T only) ─────────────────────────────────
    if (opts.temperatures.size() > 1) {
        std::string scan_file = opts.prefix + "_tscan.csv";
        std::ofstream ofs(scan_file);
        if (ofs.is_open()) {
            ofs << "temperature";
            for (size_t t = 1; t <= tgt_enms.size(); ++t) {
                ofs << ",delta_S_vib_" << t << ",delta_F_vib_" << t;
            }
            ofs << "\n";

            for (double T : opts.temperatures) {
                ofs << std::fixed << std::setprecision(1) << T;
                auto ref_encom = tencom_diff::to_encom_modes(ref_enm.modes());
                for (size_t t = 0; t < tgt_enms.size(); ++t) {
                    auto diff = tencom_diff::compute_differential(
                        ref_enm, tgt_enms[t],
                        opts.pdb_files[0], opts.pdb_files[t + 1],
                        T, opts.frequency_calibration, opts.eigenvalue_cutoff);
                    ofs << "," << std::setprecision(8) << diff.delta_S_vib
                        << "," << std::setprecision(6) << diff.delta_F_vib;
                }
                ofs << "\n";
            }
            ofs.close();
            if (V >= 1) std::cout << "\n  Wrote temperature scan: " << scan_file << "\n";
        }
    }

    if (V >= 1) std::cout << "\ntENCoM analysis complete.\n";

    return 0;
  } catch (const FlexAIDException& e) {
    if (e.exit_code() == 0) return 0;  // help requested
    std::cerr << "FlexAID Error: " << e.what() << "\n";
    return e.exit_code();
  } catch (const std::exception& e) {
    std::cerr << "Fatal error: " << e.what() << "\n";
    return 1;
  }
}
