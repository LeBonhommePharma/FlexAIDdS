// TwoStageScreen.cpp — Two-stage virtual screening pipeline
//
// Stage 1: NRGRank coarse-grained CF screening
// Stage 2: FlexAIDdS full GA docking on top-N
//
// SPDX-License-Identifier: Apache-2.0

#include "TwoStageScreen.h"
#include "pprop.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <unordered_map>

namespace nrgrank {

TwoStageScreener::TwoStageScreener()  = default;
TwoStageScreener::~TwoStageScreener() = default;

void TwoStageScreener::set_config(const TwoStageConfig& cfg) {
    config_ = cfg;
    screener_.set_config(cfg.coarse);
}

void TwoStageScreener::prepare_target(const std::vector<TargetAtom>& atoms,
                                       const std::vector<BindingSiteSphere>& spheres) {
    screener_.prepare_target(atoms, spheres);
}

bool TwoStageScreener::load_target(const std::string& mol2_path,
                                    const std::string& cleft_pdb) {
    auto atoms   = parse_target_mol2(mol2_path);
    auto spheres = parse_binding_site_pdb(cleft_pdb);

    if (atoms.empty() || spheres.empty()) return false;

    prepare_target(atoms, spheres);
    return screener_.is_prepared();
}

std::vector<TwoStageResult> TwoStageScreener::run(
        const std::vector<ScreenLigand>& ligands) {

    namespace chr = std::chrono;
    auto t0 = chr::steady_clock::now();

    // ── Stage 1: Coarse screening ──
    if (config_.verbose) {
        std::printf("Stage 1: Coarse screening %zu ligands "
                    "(%zu anchors × %d³ rotations)...\n",
                    ligands.size(),
                    screener_.num_anchors(),
                    config_.coarse.rotations_per_axis);
    }

    auto coarse_results = screener_.screen(ligands);

    auto t1 = chr::steady_clock::now();
    double stage1_sec = chr::duration<double>(t1 - t0).count();

    if (config_.verbose) {
        std::printf("Stage 1 complete: %.2f s (%.0f ligands/s)\n",
                    stage1_sec,
                    ligands.size() / std::max(stage1_sec, 1e-9));
    }

    // Write Stage 1 CSV if requested
    if (config_.write_coarse_csv) {
        std::filesystem::create_directories(config_.output_dir);
        std::string csv_path = config_.output_dir + "/coarse_screen.csv";
        write_csv(csv_path, coarse_results);
        if (config_.verbose)
            std::printf("Stage 1 results written to %s\n", csv_path.c_str());
    }

    // ── Build output with rank info ──
    const int n_results = static_cast<int>(coarse_results.size());
    const int n_stage2  = std::min(config_.top_n, n_results);

    std::vector<TwoStageResult> results(n_results);
    for (int i = 0; i < n_results; ++i) {
        results[i].coarse_result = coarse_results[i];
        results[i].coarse_rank   = i + 1;
    }

    // ── Stage 2: Full docking on top-N ──
    if (dock_cb_ && n_stage2 > 0) {
        if (config_.verbose)
            std::printf("Stage 2: Full docking top %d hits...\n", n_stage2);

        auto t2_start = chr::steady_clock::now();

        // Find the corresponding ScreenLigand for each top hit
        // Build name→ligand index map
        std::unordered_map<std::string, size_t> name_to_idx;
        for (size_t i = 0; i < ligands.size(); ++i)
            name_to_idx[ligands[i].name] = i;

        for (int i = 0; i < n_stage2; ++i) {
            if (!flexaids::pprop_keep(results[i].coarse_rank, n_results))
                continue;
            auto it = name_to_idx.find(results[i].coarse_result.name);
            if (it == name_to_idx.end()) continue;

            const auto& lig = ligands[it->second];
            results[i].full_dock_score = dock_cb_(lig, results[i].coarse_result);
            results[i].stage2_kind = config_.stage2_kind_label.empty()
                ? "callback" : config_.stage2_kind_label;
        }

        // Re-sort top-N by full dock score. NaN (pProp-dropped / missing
        // callback) sorts last via stage2_score_less — not operator< on NaN.
        std::sort(results.begin(), results.begin() + n_stage2, stage2_score_less);

        int docked_rank = 0;
        for (int i = 0; i < n_stage2; ++i) {
            if (std::isfinite(results[i].full_dock_score))
                results[i].full_rank = ++docked_rank;
            else
                results[i].full_rank = 0;
        }

        auto t2_end = chr::steady_clock::now();
        if (config_.verbose) {
            double s2_sec = chr::duration<double>(t2_end - t2_start).count();
            std::printf("Stage 2 complete: %.2f s\n", s2_sec);
        }
    }

    return results;
}

void TwoStageScreener::write_csv(const std::string& path,
                                  const std::vector<ScreenResult>& results) {
    std::ofstream fout(path);
    if (!fout.is_open()) return;

    fout << "Rank,Name,Score,BestAnchorX,BestAnchorY,BestAnchorZ\n";
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        fout << (i + 1) << ","
             << "\"" << r.name << "\","
             << r.score << ","
             << r.best_position.x << ","
             << r.best_position.y << ","
             << r.best_position.z << "\n";
    }
}

void TwoStageScreener::write_unified_csv(const std::string& path,
                                         const std::vector<TwoStageResult>& results) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(fs::path(path).parent_path(), ec);
    std::ofstream fout(path);
    if (!fout.is_open()) return;

    fout << "coarse_rank,name,coarse_score,stage2_score,stage2_rmsd,stage2_rank,stage2_kind\n";
    for (const auto& r : results) {
        fout << r.coarse_rank << ","
             << "\"" << r.coarse_result.name << "\","
             << r.coarse_result.score << ",";
        if (std::isfinite(r.full_dock_score)) fout << r.full_dock_score;
        fout << ",";
        if (std::isfinite(r.rmsd)) fout << r.rmsd;
        fout << ","
             << r.full_rank << ","
             << r.stage2_kind << "\n";
    }
}

bool TwoStageScreener::write_screen_receipt(const std::string& output_dir,
                                            int n_ligands,
                                            int top_n,
                                            bool stage2_callback,
                                            const std::string& extra_mode) {
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(output_dir, ec);
    const std::string path = output_dir + "/RUN_RECEIPT.json";
    std::ofstream fout(path);
    if (!fout.is_open()) return false;
    const char* mode = extra_mode.empty()
        ? (stage2_callback ? "two-stage-screen" : "coarse-screen")
        : extra_mode.c_str();
    fout << "{\n"
         << "  \"schema\": \"flexaidds_screen_receipt/1\",\n"
         << "  \"mode\": \"" << mode << "\",\n"
         << "  \"n_ligands\": " << n_ligands << ",\n"
         << "  \"top_n\": " << top_n << ",\n"
         << "  \"stage2_callback\": " << (stage2_callback ? "true" : "false") << ",\n"
         << "  \"score_claim\": \"CF/contact-function scoring proxy (coarse); "
            "stage2 is callback/surrogate unless a real GA dock is wired\",\n"
         << "  \"method_doi\": \"10.1101/2025.02.17.638675\",\n"
         << "  \"zenodo_doi\": \"10.5281/zenodo.16861024\",\n"
         << "  \"matrix\": \"MC_5p_norm_P10_M2_2 (constexpr LIB/nrgrank_matrix.h)\",\n"
         << "  \"license_note\": \"clean-room Apache-2.0; NRGlab/NRGRank GPL source is not a dependency\"\n"
         << "}\n";
    return true;
}

std::vector<std::string> TwoStageScreener::top_names(
        const std::vector<TwoStageResult>& results, int top_n) {
    std::vector<std::string> names;
    const int n = std::min(std::max(top_n, 0), static_cast<int>(results.size()));
    names.reserve(static_cast<size_t>(n));
    const int ntot = static_cast<int>(results.size());
    for (int i = 0; i < ntot && static_cast<int>(names.size()) < n; ++i) {
        if (!flexaids::pprop_keep(results[static_cast<size_t>(i)].coarse_rank, ntot))
            continue;
        names.push_back(results[static_cast<size_t>(i)].coarse_result.name);
    }
    return names;
}

} // namespace nrgrank
