// =============================================================================
// benchmark_datasets.cpp — Standalone benchmark dataset runner executable
//
// Usage:
//   benchmark_datasets --benchmark astex [--output results/] [--threads 8]
//   benchmark_datasets --benchmark casf2016
//   benchmark_datasets --benchmark all
//   benchmark_datasets --benchmark doi:10.1021/acs.jcim.3c00817
//   benchmark_datasets --benchmark pdb_list:my_targets.txt
//
// Supported benchmarks:
//   astex, astex_nonnative, hap2, casf2016, posebusters, dude,
//   bindingdb_itc, sampl6, sampl7, pdbbind, all
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

#include "DatasetRunner.h"
#include "BenchmarkRunner.h"
#include "FleetRunner.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;

static void print_usage(const char* progname) {
    printf("FlexAIDdS Benchmark Dataset Runner\n\n");
    printf("Usage:\n");
    printf("  %s --benchmark <dataset> [options]\n\n", progname);
    printf("Datasets:\n");
    printf("  astex            Astex Diverse 85 (Hartshorn et al. 2007)\n");
    printf("  astex_nonnative  Astex Non-Native 1112 (Verdonk et al. 2008)\n");
    printf("  hap2             HAP2 59 targets (Gaudreault & Najmanovich 2015)\n");
    printf("  casf2016         CASF-2016 285 complexes (PDBbind core set)\n");
    printf("  posebusters      PoseBusters 308 (Buttenschoen et al. 2024)\n");
    printf("  dude             DUD-E 102 targets (Mysinger et al. 2012)\n");
    printf("  bindingdb_itc    BindingDB ITC thermodynamic data\n");
    printf("  sampl6           SAMPL6 host-guest 27 (ITC: dG, dH, TdS)\n");
    printf("  sampl7           SAMPL7 host-guest ~30 (ITC: dG, dH, TdS)\n");
    printf("  pdbbind          PDBbind Refined 5316 (v2020)\n");
    printf("  all              Run all standard benchmarks\n");
    printf("  doi:<DOI>              Parse PDB codes from a DOI\n");
    printf("  pdb_list:<file>        Load PDB codes from a text file\n");
    printf("  crossdock_json:<file>  Cross-docking pairs from a JSON spec file\n\n");
    printf("Options:\n");
    printf("  --output <dir>        Output directory (default: benchmark_results/)\n");
    printf("  --threads <N>         Concurrent FlexAIDdS workers (default: 1)\n");
    printf("  --omp-threads <N>     OMP threads per worker (default: auto = hw_cores/workers)\n");
    printf("                        Rule: workers × omp-threads ≤ physical P-cores\n");
    printf("                        M3 Pro optimal: --threads 1 --omp-threads 6\n");
    printf("  --gpu <backend>       Enable GPU (cuda or metal)\n");
    printf("  --cache <dir>         Cache directory (default: ~/.flexaidds/benchmarks/)\n");
    printf("  --force               Re-run even if results already exist\n");
    printf("  --prepare-only        Download and prepare only (no docking)\n");
    printf("  --list-codes          List PDB codes for a dataset and exit\n");
    printf("  --only-codes <list>   Restrict dataset run to comma/space-separated PDB codes, or a file\n");
    printf("  --ga-generations <N>  GA generations (default: 500)\n");
    printf("  --ga-population <N>   GA population size (default: 1000)\n");
    printf("  --grid-spacing <F>    Grid spacing in Å (default: 0.375; use 0.5 for coarse pass)\n");
    printf("  --job-timeout-seconds <N>  Per-complex timeout in s (default: 3600)\n");
    printf("  --fleet               Enable immutable Fleet chunk-result output\n");
    printf("  --campaign-id <ID>    Fleet campaign identifier (required with --fleet)\n");
    printf("  --chunk-id <ID>       Fleet chunk identifier (required with --fleet)\n");
    printf("  --attempt-id <ID>     Fleet attempt identifier (required with --fleet)\n");
    printf("  --worker-id <ID>      Fleet worker identifier (required with --fleet)\n");
    printf("  --manifest-sha256 <H> SHA-256 of immutable campaign manifest\n");
    printf("  --runner-sha256 <H>   SHA-256 of this benchmark runner\n");
    printf("  --engine-sha256 <H>   SHA-256 of the pinned FlexAIDdS engine\n");
    printf("  --output-json <path>  Immutable Fleet chunk JSON destination\n");
    printf("  --mode <mode>         Benchmark protocol (Layer 1):\n");
    printf("                        oracle-ceiling  seed_elitism=ON,  blinding=OFF (ceiling)\n");
    printf("                        defined-cleft-redock  seed_elitism=OFF, blinding=ON, known cleft/site\n");
    printf("                        autonomous      seed_elitism=OFF, blinding=ON  (thesis number)\n");
    printf("                        (default: unset — reads FLEXAIDDS_SEED_ELITISM env var)\n");
    printf("  -h, --help            Show this help\n\n");
    printf("Thread sizing guide (M3 Pro, 6 P-cores):\n");
    printf("  --threads 1 --omp-threads 6   → 6 min/complex, optimal throughput\n");
    printf("  --threads 2 --omp-threads 3   → 10 min/complex, 2× dataset parallelism\n");
    printf("  --threads 2                   → AUTO: 6/2=3 OMP threads/worker (same as above)\n");
    printf("  --threads 2 (no --omp-threads, OMP_NUM_THREADS=6) → BUG: 12 threads on 6 cores\n\n");
    printf("Examples:\n");
    printf("  %s --benchmark astex --threads 1 --omp-threads 6\n", progname);
    printf("  %s --benchmark astex --prepare-only\n", progname);
    printf("  %s --benchmark casf2016 --threads 1 --omp-threads 6\n", progname);
    printf("  %s --benchmark all --threads 1 --omp-threads 6\n", progname);
    printf("  %s --benchmark doi:10.1021/acs.jcim.3c00817\n", progname);
    printf("  %s --benchmark astex --list-codes\n", progname);
}

static void print_publication_table(const dataset::BenchmarkReport& report) {
    // Print a publication-ready summary table matching manuscript format
    printf("\n");
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("  FlexAIDdS Benchmark: %s\n", report.dataset_name.c_str());
    printf("═══════════════════════════════════════════════════════════════════\n");
    printf("\n");
    printf("  ┌─────────────────────────────┬────────────────────┐\n");
    printf("  │ Metric                      │ Value              │\n");
    printf("  ├─────────────────────────────┼────────────────────┤\n");
    printf("  │ Total systems               │ %18d │\n", report.total_systems);
    printf("  │ Successful (RMSD <= 2.0 A)  │ %18d │\n", report.successful);
    printf("  │ Success rate                │ %17.1f%% │\n", report.success_rate * 100.0);
    printf("  │ Valid RMSDs                 │ %18d │\n", report.valid_rmsd_count);
    if (report.valid_rmsd_count > 0 && std::isfinite(report.mean_rmsd) && std::isfinite(report.median_rmsd)) {
        printf("  │ Mean RMSD (Å)               │ %18.2f │\n", report.mean_rmsd);
        printf("  │ Median RMSD (Å)             │ %18.2f │\n", report.median_rmsd);
    } else {
        printf("  │ Mean RMSD (Å)               │ %18s │\n", "NA");
        printf("  │ Median RMSD (Å)             │ %18s │\n", "NA");
    }
    printf("  │ Affinity pairs              │ %18d │\n", report.affinity_pairs);
    if (report.affinity_pairs >= 3 &&
        std::isfinite(report.pearson_r) &&
        std::isfinite(report.spearman_rho) &&
        std::isfinite(report.kendall_tau)) {
        printf("  │ Pearson r                   │ %18.3f │\n", report.pearson_r);
        printf("  │ Spearman ρ                  │ %18.3f │\n", report.spearman_rho);
        printf("  │ Kendall τ                   │ %18.3f │\n", report.kendall_tau);
    } else {
        printf("  │ Pearson r                   │ %18s │\n", "NA");
        printf("  │ Spearman ρ                  │ %18s │\n", "NA");
        printf("  │ Kendall τ                   │ %18s │\n", "NA");
    }
    printf("  └─────────────────────────────┴────────────────────┘\n");
    printf("\n");
}

static void print_itc_table(const dataset::BenchmarkReport& report,
                             const std::vector<dataset::DatasetEntry>& entries) {
    // Print ITC-specific thermodynamic comparison table
    bool has_itc = false;
    for (const auto& e : entries) {
        if (e.has_enthalpy()) { has_itc = true; break; }
    }
    if (!has_itc) return;

    // Claim firewall. The `predicted_*` side of every pair below comes from the
    // uncalibrated contact-function optimizer (EnergyDomain =
    // ContactFunctionArbitraryUnits, EnsembleMeasure = OptimizerSamples, no
    // sha256 receipt ⇒ ProxyOnly per LIB/statmech.h). Correlating a proxy
    // against experiment is a legitimate diagnostic — the correlation
    // coefficients are dimensionless and unaffected by the missing calibration
    // — but the proxy columns themselves must not be presented as physical
    // thermodynamics. Only the labels/units change here; the statistics are
    // computed from exactly the same inputs as before.
    printf("\n");
    printf("  ITC correlation diagnostic — CF-proxy predictions vs experiment\n");
    printf("  claim_validity: proxy_only. Predicted values are uncalibrated\n");
    printf("  contact-function diagnostics in arbitrary units, not kcal/mol; only\n");
    printf("  the dimensionless correlations below are interpretable. No ΔG, ΔH,\n");
    printf("  TΔS, Kd, Ki or affinity value is claimed for the predicted side.\n");
    printf("  ─────────────────────────────────────────────────────────────\n");

    // Collect ITC pairs
    std::vector<double> exp_dG, pred_dG;
    std::vector<double> exp_dH, pred_dH;
    std::vector<double> exp_TdS, pred_TdS;

    for (size_t i = 0; i < entries.size() && i < report.results.size(); ++i) {
        const auto& entry = entries[i];
        const auto& result = report.results[i];

        if (entry.has_affinity() && result.predicted_dG != 0.0f) {
            exp_dG.push_back(-entry.experimental_affinity * 1.3636);
            pred_dG.push_back(result.predicted_dG);
        }
        if (entry.has_enthalpy() && result.predicted_dH != 0.0f) {
            exp_dH.push_back(entry.experimental_dH);
            pred_dH.push_back(result.predicted_dH);
        }
        if (entry.has_entropy() && result.predicted_TdS != 0.0f) {
            exp_TdS.push_back(entry.experimental_TdS);
            pred_TdS.push_back(result.predicted_TdS);
        }
    }

    // Row labels name the PROXY channel that was correlated against the
    // experimental ITC column, e.g. "ΔG-like proxy" = predicted_dG (CF units)
    // vs experimental ΔG (kcal/mol). Cell width is 16 columns; keep the label
    // ≤ 14 display columns so the box stays aligned.
    printf("  ┌────────────────┬──────────┬──────────┬──────────┐\n");
    printf("  │ Proxy channel  │ Pearson  │ Spearman │ Kendall  │\n");
    printf("  ├────────────────┼──────────┼──────────┼──────────┤\n");

    if (exp_dG.size() >= 3) {
        printf("  │ ΔG-like proxy  │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_dG, exp_dG),
               dataset::compute_spearman_rho(pred_dG, exp_dG),
               dataset::compute_kendall_tau(pred_dG, exp_dG));
    }
    if (exp_dH.size() >= 3) {
        printf("  │ ΔH-like proxy  │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_dH, exp_dH),
               dataset::compute_spearman_rho(pred_dH, exp_dH),
               dataset::compute_kendall_tau(pred_dH, exp_dH));
    }
    if (exp_TdS.size() >= 3) {
        printf("  │ TΔS-like proxy │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_TdS, exp_TdS),
               dataset::compute_spearman_rho(pred_TdS, exp_TdS),
               dataset::compute_kendall_tau(pred_TdS, exp_TdS));
    }

    printf("  └────────────────┴──────────┴──────────┴──────────┘\n");
    printf("  Correlation only — proxy columns are in arbitrary CF units.\n");
    printf("\n");
}

static void list_pdb_codes(dataset::BenchmarkSet set) {
    std::vector<std::string> codes;
    switch (set) {
        case dataset::BenchmarkSet::ASTEX_DIVERSE:
            codes = dataset::DatasetRunner::astex_diverse_codes();
            break;
        case dataset::BenchmarkSet::CASF_2016:
            codes = dataset::DatasetRunner::casf2016_codes();
            break;
        case dataset::BenchmarkSet::DUD_E:
            codes = dataset::DatasetRunner::dude_targets();
            break;
        case dataset::BenchmarkSet::HAP2:
            codes = dataset::DatasetRunner::hap2_codes();
            break;
        default:
            printf("No hardcoded PDB list for this dataset. Use --prepare-only to fetch.\n");
            return;
    }

    printf("%s — %zu entries:\n", dataset::benchmark_set_name(set).c_str(), codes.size());
    int col = 0;
    for (const auto& code : codes) {
        printf("%-6s", code.c_str());
        if (++col % 12 == 0) printf("\n");
    }
    if (col % 12 != 0) printf("\n");
}

static std::string uppercase_code(std::string code) {
    std::transform(code.begin(), code.end(), code.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return code;
}

static std::vector<std::string> parse_only_codes(const std::string& spec) {
    std::vector<std::string> codes;
    if (spec.empty()) return codes;

    std::stringstream input;
    if (fs::exists(spec)) {
        std::ifstream ifs(spec);
        if (!ifs) {
            std::cerr << "ERROR: cannot open --only-codes file: " << spec << "\n";
            return codes;
        }
        input << ifs.rdbuf();
    } else {
        input << spec;
    }

    std::string text = input.str();
    for (char& ch : text) {
        if (ch == ',' || ch == ';' || ch == '\n' || ch == '\r' || ch == '\t') ch = ' ';
    }

    std::istringstream iss(text);
    std::string code;
    while (iss >> code) {
        auto hash = code.find('#');
        if (hash == 0) {
            std::string ignored;
            std::getline(iss, ignored);
            continue;
        }
        if (hash != std::string::npos) code.erase(hash);
        code = uppercase_code(code);
        if (code.size() == 4) codes.push_back(code);
    }
    return codes;
}

static void filter_entries_by_code(std::vector<dataset::DatasetEntry>& entries,
                                   const std::vector<std::string>& only_codes) {
    if (only_codes.empty()) return;

    std::set<std::string> wanted;
    for (const auto& code : only_codes) wanted.insert(uppercase_code(code));

    std::vector<dataset::DatasetEntry> filtered;
    filtered.reserve(entries.size());
    for (auto& entry : entries) {
        if (wanted.count(uppercase_code(entry.pdb_id))) {
            filtered.push_back(std::move(entry));
        }
    }

    std::cout << "  --only-codes selected " << filtered.size()
              << " / " << wanted.size() << " requested entries\n";
    entries = std::move(filtered);
}

static dataset::BenchmarkReport run_single_benchmark(const std::string& name,
                                  dataset::DatasetRunner& runner,
                                  const dataset::DockingConfig& config,
                                  bool prepare_only,
                                  bool list_codes_only,
                                  const std::vector<std::string>& only_codes) {
    using BS = dataset::BenchmarkSet;

    // Check for special prefixes: doi: and pdb_list:
    if (name.substr(0, 4) == "doi:") {
        std::string doi = name.substr(4);
        auto entries = runner.prepare_from_doi(doi);
        filter_entries_by_code(entries, only_codes);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
        }
        return {};
    }
    if (name.substr(0, 9) == "pdb_list:") {
        std::string file_path = name.substr(9);
        auto entries = runner.prepare_from_pdb_list(file_path);
        filter_entries_by_code(entries, only_codes);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
        }
        return {};
    }
    if (name.substr(0, 15) == "crossdock_json:") {
        // Cross-docking from a JSON pairs file.
        // JSON format (benchmark_crossdock_85.json):
        //   { "pairs": [ { "receptor_id": "1G9V", "ligand_id": "1GM8",
        //                   "receptor_pdb": "/abs/path/1G9V_apo.pdb",
        //                   "ligand_sdf":   "/abs/path/1GM8_ligand.sdf",
        //                   "oracle_site_pdb": "/abs/path/1G9V_binding_site.pdb",
        //                   "cleft_sphere_file": "/abs/path/1G9V_sph_1.pdb" }, ... ] }
        std::string json_file = name.substr(15);
        // Expand leading ~
        if (!json_file.empty() && json_file[0] == '~') {
            if (const char* home = std::getenv("HOME"))
                json_file = std::string(home) + json_file.substr(1);
        }
        std::ifstream ifs(json_file);
        if (!ifs) {
            fprintf(stderr, "ERROR: cannot open crossdock_json file: %s\n", json_file.c_str());
            return {};
        }
        std::string content((std::istreambuf_iterator<char>(ifs)),
                             std::istreambuf_iterator<char>());
        const fs::path json_base = fs::absolute(fs::path(json_file)).parent_path();

        // Minimal JSON string-field extractor (no external JSON lib needed).
        auto extract_str = [&](const std::string& obj, const std::string& key) -> std::string {
            std::string needle = "\"" + key + "\": \"";
            auto pos = obj.find(needle);
            if (pos == std::string::npos) return "";
            pos += needle.size();
            auto end = obj.find('"', pos);
            return (end != std::string::npos) ? obj.substr(pos, end - pos) : "";
        };
        auto resolve_input_path = [&](const std::string& value) -> std::string {
            if (value.empty()) return value;
            fs::path path(value);
            if (path.is_relative()) path = json_base / path;
            std::error_code ec;
            fs::path normalized = fs::weakly_canonical(path, ec);
            return ec ? path.lexically_normal().string() : normalized.string();
        };

        // Walk JSON locating each pair object by scanning for "receptor_id" keys.
        std::vector<dataset::DatasetEntry> entries;
        std::size_t scan = 0;
        while (true) {
            auto kpos = content.find("\"receptor_id\"", scan);
            if (kpos == std::string::npos) break;
            auto obj_start = content.rfind('{', kpos);
            auto obj_end   = content.find('}', kpos);
            if (obj_start == std::string::npos || obj_end == std::string::npos) break;
            std::string obj = content.substr(obj_start, obj_end - obj_start + 1);

            dataset::DatasetEntry entry;
            entry.pdb_id            = extract_str(obj, "receptor_id");
            entry.receptor_path     = resolve_input_path(extract_str(obj, "receptor_pdb"));
            entry.ligand_path       = resolve_input_path(extract_str(obj, "ligand_sdf"));
            entry.rmsd_reference_path = resolve_input_path(extract_str(obj, "rmsd_ref_sdf"));
            entry.binding_site_path = resolve_input_path(extract_str(obj, "oracle_site_pdb"));
            entry.cleft_sphere_path = resolve_input_path(extract_str(obj, "cleft_sphere_file"));
            entry.source            = "astex_crossdock_85";

            if (!entry.pdb_id.empty() &&
                !entry.receptor_path.empty() &&
                !entry.ligand_path.empty()) {
                entries.push_back(std::move(entry));
            }
            scan = obj_end + 1;
        }

        printf("  crossdock_json: %s — loaded %zu pairs\n",
               json_file.c_str(), entries.size());
        filter_entries_by_code(entries, only_codes);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
            return report;
        }
        return {};
    }

    auto bs = dataset::parse_benchmark_set(name);
    if (!bs.has_value()) {
        fprintf(stderr, "ERROR: Unknown benchmark: '%s'\n", name.c_str());
        fprintf(stderr, "Use --help for available datasets.\n");
        return {};
    }

    if (list_codes_only) {
        list_pdb_codes(*bs);
        return {};
    }

    auto entries = runner.prepare(*bs);
    filter_entries_by_code(entries, only_codes);
    printf("  → %zu entries prepared\n", entries.size());
    // P1 diagnostic (additive, non-behavior): explicit T + progress for early diagnosis of best BindingMode runs
    printf("  [P1] Docking phase starting for best BindingMode search at temperature from config (exact 298/310 K fidelity required). Live progress + run_status sidecar will be emitted during run(entries).\n");

    if (prepare_only) {
        printf("  [prepare-only mode] Skipping docking.\n");
        return {};
    }

    if (!entries.empty()) {
        auto report = runner.run(entries, config);
        print_publication_table(report);
        print_itc_table(report, entries);
        runner.write_report(report, config.output_dir);
        return report;
    }
    return {};
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Keep the benchmark runner and docking engine from the same build tree.
    // An explicit FLEXAIDDS_BINARY still wins, but absent that override the
    // sibling FlexAIDdS executable is authoritative. This prevents a
    // build_lto/benchmark_datasets invocation from silently using build/FlexAIDdS.
    if (std::getenv("FLEXAIDDS_BINARY") == nullptr) {
        std::error_code exe_ec;
        fs::path runner_path = fs::weakly_canonical(fs::absolute(argv[0]), exe_ec);
        if (!exe_ec) {
            const fs::path sibling = runner_path.parent_path() / "FlexAIDdS";
            if (fs::is_regular_file(sibling, exe_ec) && !exe_ec) {
#ifdef _WIN32
                _putenv_s("FLEXAIDDS_BINARY", sibling.string().c_str());
#else
                setenv("FLEXAIDDS_BINARY", sibling.string().c_str(), 0);
#endif
            }
        }
    }

    // Parse arguments
    std::string benchmark_name;
    std::string output_dir = "benchmark_results";
    std::string cache_dir;
    int threads = 1;
    int omp_threads = 0;
    int job_timeout_s = 3600;
    bool use_gpu = false;
    std::string gpu_backend = "cuda";
    bool prepare_only = false;
    bool list_codes_only = false;
    bool force_rerun = false;
    int ga_generations = 0;
    int ga_population = 0;
    double temperature = 0.0;
    double grid_spacing = 0.0;
    std::string clustering;
    std::vector<std::string> only_codes;
    // Fleet mode options
    bool fleet_mode = false;
    std::string campaign_id;
    std::string chunk_id;
    std::string attempt_id;
    std::string worker_id;
    std::string manifest_sha256;
    std::string runner_sha256;
    std::string engine_sha256;
    std::string output_json;
    // Layer 1: benchmark protocol mode
    std::string mode_str;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        if (arg == "--benchmark" && i + 1 < argc) {
            benchmark_name = argv[++i];
            continue;
        }
        if (arg == "--output" && i + 1 < argc) {
            output_dir = argv[++i];
            continue;
        }
        if (arg == "--cache" && i + 1 < argc) {
            cache_dir = argv[++i];
            continue;
        }
        if (arg == "--threads" && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--omp-threads" && i + 1 < argc) {
            omp_threads = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--gpu" && i + 1 < argc) {
            use_gpu = true;
            gpu_backend = argv[++i];
            continue;
        }
        if (arg == "--prepare-only") {
            prepare_only = true;
            continue;
        }
        if (arg == "--list-codes") {
            list_codes_only = true;
            continue;
        }
        if (arg == "--only-codes" && i + 1 < argc) {
            only_codes = parse_only_codes(argv[++i]);
            continue;
        }
        if (arg == "--force") {
            force_rerun = true;
            continue;
        }
        if (arg == "--ga-generations" && i + 1 < argc) {
            ga_generations = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--ga-population" && i + 1 < argc) {
            ga_population = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--grid-spacing" && i + 1 < argc) {
            grid_spacing = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--temperature" && i + 1 < argc) {
            temperature = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--clustering" && i + 1 < argc) {
            clustering = argv[++i];
            continue;
        }
        if (arg == "--job-timeout-seconds" && i + 1 < argc) {
            job_timeout_s = std::atoi(argv[++i]);
            continue;
        }
        if (arg == "--fleet") {
            fleet_mode = true;
            continue;
        }
        if (arg == "--campaign-id" && i + 1 < argc) {
            campaign_id = argv[++i];
            continue;
        }
        if (arg == "--chunk-id" && i + 1 < argc) {
            chunk_id = argv[++i];
            continue;
        }
        if (arg == "--attempt-id" && i + 1 < argc) {
            attempt_id = argv[++i];
            continue;
        }
        if (arg == "--worker-id" && i + 1 < argc) {
            worker_id = argv[++i];
            continue;
        }
        if (arg == "--manifest-sha256" && i + 1 < argc) {
            manifest_sha256 = argv[++i];
            continue;
        }
        if (arg == "--runner-sha256" && i + 1 < argc) {
            runner_sha256 = argv[++i];
            continue;
        }
        if (arg == "--engine-sha256" && i + 1 < argc) {
            engine_sha256 = argv[++i];
            continue;
        }
        if (arg == "--output-json" && i + 1 < argc) {
            output_json = argv[++i];
            continue;
        }
        if (arg == "--mode" && i + 1 < argc) {
            mode_str = argv[++i];
            continue;
        }

        // Fallback: if first positional arg, treat as benchmark name
        if (benchmark_name.empty()) {
            benchmark_name = arg;
        }
    }

    if (benchmark_name.empty()) {
        fprintf(stderr, "ERROR: No benchmark specified. Use --benchmark <name>\n");
        print_usage(argv[0]);
        return 1;
    }

    if (fleet_mode) {
        const bool missing_metadata = campaign_id.empty() || chunk_id.empty() ||
            attempt_id.empty() || worker_id.empty() || manifest_sha256.empty() ||
            runner_sha256.empty() || engine_sha256.empty() || output_json.empty();
        if (missing_metadata) {
            fprintf(stderr, "ERROR: --fleet requires campaign/chunk/attempt/worker IDs, "
                    "manifest/runner/engine SHA-256 values, and --output-json\n");
            return 1;
        }
        if (benchmark_name == "all" || prepare_only || list_codes_only ||
            only_codes.empty() || mode_str.empty() || output_json == "-") {
            fprintf(stderr, "ERROR: Fleet mode requires one explicit benchmark, --only-codes, "
                    "an explicit --mode, and a file output; prepare/list/all modes are unsupported\n");
            return 1;
        }
    }

    // Create runner and config
    dataset::DatasetRunner runner(cache_dir);

    dataset::DockingConfig config;
    config.num_threads            = threads;
    config.omp_threads_per_worker = omp_threads;   // 0 → auto-detect in DatasetRunner
    config.use_gpu                = use_gpu;
    config.gpu_backend            = gpu_backend;
    config.output_dir             = output_dir;
    config.skip_completed         = !force_rerun;
    if (ga_generations > 0)       config.ga_generations    = ga_generations;
    if (ga_population  > 0)       config.ga_population     = ga_population;
    if (grid_spacing   > 0.0)     config.grid_spacing      = static_cast<float>(grid_spacing);
    if (temperature    > 0.0)     config.temperature       = static_cast<float>(temperature);
    if (job_timeout_s  > 0)       config.per_job_timeout_s = job_timeout_s;
    if (!clustering.empty())      config.clustering_algorithm = clustering;

    // Layer 1: explicit benchmark protocol mode
    if (mode_str == "oracle-ceiling") {
        config.mode = dataset::BenchmarkMode::ORACLE_CEILING;
    } else if (mode_str == "defined-cleft-redock" || mode_str == "cognate-redock") {
        config.mode = dataset::BenchmarkMode::DEFINED_CLEFT_REDOCK;
    } else if (mode_str == "autonomous") {
        config.mode = dataset::BenchmarkMode::AUTONOMOUS;
    } else if (!mode_str.empty()) {
        fprintf(stderr, "ERROR: Unknown --mode '%s'. Use 'oracle-ceiling', 'defined-cleft-redock', or 'autonomous'\n",
                mode_str.c_str());
        return 1;
    }

    // Ablation hook: FLEXAIDDS_FORCE_RIGID re-pins legacy rigid-body docking
    // (DatasetRunner writes flexibility.intramolecular=false → engine builds a
    // 4-gene chromosome: translation + rotation only, no ligand torsional DoF).
    // Isolates whether a flexible-docking regression is eval-budget dilution
    // (more genes, same generations) vs the oracle-confinement changes.
    if (const char* fr = std::getenv("FLEXAIDDS_FORCE_RIGID")) {
        if (fr[0] && fr[0] != '0') {
            config.force_rigid = true;
            std::cout << "  FORCE_RIGID:  ON (intramolecular=false, num_genes=4)\n";
        }
    }

    // Override receptor rotamer prep gate (default true since v44).
    // Set FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0 to disable for cross-docking
    // benchmarks where the native ligand is absent and sidechain prep against
    // a ghost occupancy makes no physical sense.
    if (const char* rrp = std::getenv("FLEXAIDDS_RECEPTOR_ROTAMER_PREP")) {
        bool enable = (rrp[0] && rrp[0] != '0' &&
                       std::string(rrp) != "false" && std::string(rrp) != "off");
        config.receptor_rotamer_prep = enable;
        std::cout << "  RECEPTOR_ROTAMER_PREP: " << (enable ? "ON" : "OFF")
                  << " (override via env)\n";
    }

    // Compute effective OMP threads for display (mirrors DatasetRunner logic)
    int effective_omp = config.omp_threads_per_worker;
    if (effective_omp <= 0) {
        const char* env_omp = std::getenv("OMP_NUM_THREADS");
        int base = (env_omp && std::atoi(env_omp) > 0)
            ? std::atoi(env_omp)
            : static_cast<int>(std::thread::hardware_concurrency());
        effective_omp = std::max(1, base / std::max(1, config.num_threads));
    }

    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  FlexAIDdS Benchmark Dataset Runner\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    std::cout << "  Cache:        " << runner.cache_dir() << "\n";
    std::cout << "  Output:       " << output_dir << "\n";
    std::cout << "  Workers:      " << threads << " concurrent FlexAIDdS process(es)\n";
    std::cout << "  OMP/worker:   " << effective_omp << " thread(s)"
              << (config.omp_threads_per_worker > 0 ? " (explicit)" : " (auto)") << "\n";
    std::cout << "  Total threads:" << (threads * effective_omp) << " across "
              << std::thread::hardware_concurrency() << " logical cores\n";
    std::cout << "  Skip done:    " << (config.skip_completed ? "yes (--force to override)" : "no") << "\n";
    if (use_gpu) {
        std::cout << "  GPU:          " << gpu_backend << "\n";
    }
    std::cout << "  GA:           pop=" << config.ga_population
              << "  gen=" << config.ga_generations
              << "  (" << (config.ga_population * config.ga_generations / 1000) << "k evals/complex)\n";
    std::cout << "  Temp:         " << config.temperature << " K\n";
    std::cout << "  Cluster:      " << config.clustering_algorithm << "\n";
    std::cout << "  Timeout/job:  " << config.per_job_timeout_s << " s\n";
    // Layer 1: mode
    {
        const char* mode_label =
            (config.mode == dataset::BenchmarkMode::ORACLE_CEILING)       ? "oracle-ceiling" :
            (config.mode == dataset::BenchmarkMode::DEFINED_CLEFT_REDOCK) ? "defined-cleft-redock" :
            (config.mode == dataset::BenchmarkMode::AUTONOMOUS)           ? "autonomous" :
                                                                            "unset (env-var)";
        std::cout << "  Mode:         " << mode_label << "\n";
    }
    if (fleet_mode) {
        std::cout << "  Fleet:        enabled\n";
        std::cout << "  Campaign:     " << campaign_id << "\n";
        std::cout << "  Chunk/attempt:" << chunk_id << "/" << attempt_id << "\n";
        std::cout << "  JSON:         " << output_json << "\n";
    }
    std::cout << "\n";

    int runtime_exit_code = 0;
    // Handle "all" benchmark
    if (benchmark_name == "all") {
        std::vector<std::string> all_benchmarks = {
            "astex", "astex_nonnative", "hap2", "casf2016",
            "posebusters", "dude", "bindingdb_itc",
            "sampl6", "sampl7"
        };

        for (const auto& name : all_benchmarks) {
            std::cout << "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
            std::cout << "  Running: " << name << "\n";
            std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n";

            auto report = run_single_benchmark(name, runner, config, prepare_only, list_codes_only, only_codes);
            runtime_exit_code = std::max(runtime_exit_code,
                dataset::benchmark_runtime_exit_code(report, prepare_only || list_codes_only));
        }

        // Print combined summary
        std::cout << "\n\n═══════════════════════════════════════════════════════════════\n";
        std::cout << "  All benchmark attempts finished. Results in: " << output_dir << "\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n";
    } else {
        const auto fleet_started = std::chrono::steady_clock::now();
        auto report = run_single_benchmark(benchmark_name, runner, config, prepare_only, list_codes_only, only_codes);

        runtime_exit_code = dataset::benchmark_runtime_exit_code(
            report, prepare_only || list_codes_only);
        if (fleet_mode) {
            const double duration_s = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - fleet_started).count();
            std::ostringstream command;
            for (int i = 0; i < argc; ++i) {
                if (i > 0) command << ' ';
                command << argv[i];
            }
            std::error_code path_error;
            fs::path runner_path = fs::weakly_canonical(fs::absolute(argv[0]), path_error);
            if (path_error) runner_path = fs::absolute(argv[0]);
            const char* engine_env = std::getenv("FLEXAIDDS_BINARY");
            fleet::ChunkMetadata metadata{
                campaign_id,
                chunk_id,
                attempt_id,
                worker_id,
                benchmark_name,
                command.str(),
                runner_path.string(),
                runner_sha256,
                engine_env ? engine_env : "",
                engine_sha256,
                manifest_sha256,
            };
            const std::string payload = fleet::FleetRunner::serialize_chunk_result(
                metadata, report, config, duration_s);
            std::string write_error;
            if (!fleet::FleetRunner::write_chunk_result_atomic(
                    output_json, payload, &write_error)) {
                std::cerr << "ERROR: Fleet result publication failed: " << write_error << "\n";
                return 3;
            }
            std::cout << "  [Fleet] immutable chunk result: " << output_json << "\n";
            if (report.total_systems == 0 || report.results.empty()) {
                std::cerr << "ERROR: Fleet chunk produced no target results\n";
                return 2;
            }
        }
    }

    if (runtime_exit_code != 0)
        std::cerr << "ERROR: Incomplete docking; inspect per-target runtime fields and child logs\n";
    return runtime_exit_code;
}
