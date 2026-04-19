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

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static void print_usage(const char* progname) {
    printf("FlexAIDdS Benchmark Dataset Runner\n\n");
    printf("Usage:\n");
    printf("  %s --benchmark <dataset> [options]\n\n", progname);
    printf("Datasets:\n");
    printf("  astex            Astex Diverse 85 (Hartshorn et al. 2007)\n");
    printf("  astex_nonnative  Astex Non-Native 1112 (Verdonk et al. 2008)\n");
    printf("  hap2             HAP2 59 targets (Gaudreault & Bhatt 2015)\n");
    printf("  casf2016         CASF-2016 285 complexes (PDBbind core set)\n");
    printf("  posebusters      PoseBusters 308 (Buttenschoen et al. 2024)\n");
    printf("  dude             DUD-E 102 targets (Mysinger et al. 2012)\n");
    printf("  bindingdb_itc    BindingDB ITC thermodynamic data\n");
    printf("  sampl6           SAMPL6 host-guest 27 (ITC: dG, dH, TdS)\n");
    printf("  sampl7           SAMPL7 host-guest ~30 (ITC: dG, dH, TdS)\n");
    printf("  pdbbind          PDBbind Refined 5316 (v2020)\n");
    printf("  all              Run all standard benchmarks\n");
    printf("  doi:<DOI>        Parse PDB codes from a DOI\n");
    printf("  pdb_list:<file>  Load PDB codes from a text file\n\n");
    printf("Options:\n");
    printf("  --output <dir>     Output directory (default: benchmark_results/)\n");
    printf("  --threads <N>      Number of threads (default: 1)\n");
    printf("  --gpu <backend>    Enable GPU (cuda or metal)\n");
    printf("  --cache <dir>      Cache directory (default: ~/.flexaidds/benchmarks/)\n");
    printf("  --force            Re-run even if results already exist (ignore cache)\n");
    printf("  --prepare-only     Download and prepare only (no docking)\n");
    printf("  --list-codes       List PDB codes for a dataset and exit\n");
    printf("  -h, --help         Show this help\n\n");
    printf("Examples:\n");
    printf("  %s --benchmark astex --output results/astex/\n", progname);
    printf("  %s --benchmark sampl6 --prepare-only\n", progname);
    printf("  %s --benchmark casf2016 --threads 8 --gpu cuda\n", progname);
    printf("  %s --benchmark all --threads 16\n", progname);
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
    printf("  │ Successful (RMSD < 2.0 Å)   │ %18d │\n", report.successful);
    printf("  │ Success rate                │ %17.1f%% │\n", report.success_rate * 100.0);
    printf("  │ Mean RMSD (Å)               │ %18.2f │\n", report.mean_rmsd);
    printf("  │ Median RMSD (Å)             │ %18.2f │\n", report.median_rmsd);
    printf("  │ Pearson r                   │ %18.3f │\n", report.pearson_r);
    printf("  │ Spearman ρ                  │ %18.3f │\n", report.spearman_rho);
    printf("  │ Kendall τ                   │ %18.3f │\n", report.kendall_tau);
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

    printf("\n");
    printf("  ITC Thermodynamic Validation\n");
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

    printf("  ┌────────────────┬──────────┬──────────┬──────────┐\n");
    printf("  │ Property       │ Pearson  │ Spearman │ Kendall  │\n");
    printf("  ├────────────────┼──────────┼──────────┼──────────┤\n");

    if (exp_dG.size() >= 3) {
        printf("  │ ΔG (kcal/mol)  │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_dG, exp_dG),
               dataset::compute_spearman_rho(pred_dG, exp_dG),
               dataset::compute_kendall_tau(pred_dG, exp_dG));
    }
    if (exp_dH.size() >= 3) {
        printf("  │ ΔH (kcal/mol)  │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_dH, exp_dH),
               dataset::compute_spearman_rho(pred_dH, exp_dH),
               dataset::compute_kendall_tau(pred_dH, exp_dH));
    }
    if (exp_TdS.size() >= 3) {
        printf("  │ TΔS (kcal/mol) │ %8.3f │ %8.3f │ %8.3f │\n",
               dataset::compute_pearson_r(pred_TdS, exp_TdS),
               dataset::compute_spearman_rho(pred_TdS, exp_TdS),
               dataset::compute_kendall_tau(pred_TdS, exp_TdS));
    }

    printf("  └────────────────┴──────────┴──────────┴──────────┘\n");
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

static void run_single_benchmark(const std::string& name,
                                  dataset::DatasetRunner& runner,
                                  const dataset::DockingConfig& config,
                                  bool prepare_only,
                                  bool list_codes_only) {
    using BS = dataset::BenchmarkSet;

    // Check for special prefixes: doi: and pdb_list:
    if (name.substr(0, 4) == "doi:") {
        std::string doi = name.substr(4);
        auto entries = runner.prepare_from_doi(doi);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
        }
        return;
    }
    if (name.substr(0, 9) == "pdb_list:") {
        std::string file_path = name.substr(9);
        auto entries = runner.prepare_from_pdb_list(file_path);
        if (!prepare_only && !entries.empty()) {
            auto report = runner.run(entries, config);
            print_publication_table(report);
            runner.write_report(report, config.output_dir);
        }
        return;
    }

    auto bs = dataset::parse_benchmark_set(name);
    if (!bs.has_value()) {
        fprintf(stderr, "ERROR: Unknown benchmark: '%s'\n", name.c_str());
        fprintf(stderr, "Use --help for available datasets.\n");
        return;
    }

    if (list_codes_only) {
        list_pdb_codes(*bs);
        return;
    }

    auto entries = runner.prepare(*bs);
    printf("  → %zu entries prepared\n", entries.size());

    if (prepare_only) {
        printf("  [prepare-only mode] Skipping docking.\n");
        return;
    }

    if (!entries.empty()) {
        auto report = runner.run(entries, config);
        print_publication_table(report);
        print_itc_table(report, entries);
        runner.write_report(report, config.output_dir);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Parse arguments
    std::string benchmark_name;
    std::string output_dir = "benchmark_results";
    std::string cache_dir;
    int threads = 1;
    bool use_gpu = false;
    std::string gpu_backend = "cuda";
    bool prepare_only = false;
    bool list_codes_only = false;
    bool force_rerun = false;
    int ga_generations = 0;
    int ga_population = 0;
    double temperature = 0.0;
    std::string clustering;

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
        if (arg == "--temperature" && i + 1 < argc) {
            temperature = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--clustering" && i + 1 < argc) {
            clustering = argv[++i];
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

    // Create runner and config
    dataset::DatasetRunner runner(cache_dir);

    dataset::DockingConfig config;
    config.num_threads = threads;
    config.use_gpu = use_gpu;
    config.gpu_backend = gpu_backend;
    config.output_dir = output_dir;
    config.skip_completed = !force_rerun;
    if (ga_generations > 0)    config.ga_generations = ga_generations;
    if (ga_population > 0)     config.ga_population = ga_population;
    if (temperature > 0.0)      config.temperature = static_cast<float>(temperature);
    if (!clustering.empty())    config.clustering_algorithm = clustering;

    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  FlexAIDdS Benchmark Dataset Runner\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    std::cout << "  Cache:   " << runner.cache_dir() << "\n";
    std::cout << "  Output:  " << output_dir << "\n";
    std::cout << "  Threads: " << threads << "\n";
    std::cout << "  Skip existing results: " << (config.skip_completed ? "yes (use --force to override)" : "no") << "\n";
    if (use_gpu) {
        std::cout << "  GPU:     " << gpu_backend << "\n";
    }
    std::cout << "  GA:      pop=" << config.ga_population << " gen=" << config.ga_generations << "\n";
    std::cout << "  Temp:    " << config.temperature << " K\n";
    std::cout << "  Cluster: " << config.clustering_algorithm << "\n";
    std::cout << "\n";

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

            run_single_benchmark(name, runner, config, prepare_only, list_codes_only);
        }

        // Print combined summary
        std::cout << "\n\n═══════════════════════════════════════════════════════════════\n";
        std::cout << "  All benchmarks completed. Results in: " << output_dir << "\n";
        std::cout << "═══════════════════════════════════════════════════════════════\n";
    } else {
        run_single_benchmark(benchmark_name, runner, config, prepare_only, list_codes_only);
    }

    return 0;
}
