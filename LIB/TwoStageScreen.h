// TwoStageScreen.h — Two-stage virtual screening pipeline
//
// Stage 1: coarse-grained cube screening (CoarseScreen)
//          Fast rigid-body CF scoring with 729 rotations × anchor grid.
//          Filters millions of compounds down to top-N.
//
// Stage 2: optional callback (campaign per-ligand dock, or a test mock).
//          Default: unset → Stage 1 only. Never implied by --screen alone.
//
// Method: DesCôteaux, Mailhot, Najmanovich, bioRxiv 2025.02.17.638675v2
// (doi 10.1101/2025.02.17.638675; Zenodo 10.5281/zenodo.16861024).
// Clean-room per docs/licensing/clean-room-policy.md. Not GPL source.
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "CoarseScreen.h"

#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace nrgrank {

/// Configuration for the two-stage pipeline
struct TwoStageConfig {
    /// Coarse screening parameters
    CoarseScreenConfig coarse;

    /// Number of top hits from Stage 1 to pass to Stage 2
    int top_n = 100;

    /// Normalise scores by atom count before ranking
    bool normalise_by_atoms = false;

    /// Output directory for results
    std::string output_dir = "screen_results";

    /// Write coarse-screen CSV
    bool write_coarse_csv = true;

    /// Verbose progress output
    bool verbose = false;

    /// Written into TwoStageResult::stage2_kind when the callback runs.
    /// Default "callback". Use "surrogate" when the hook is not a real GA CF.
    std::string stage2_kind_label = "callback";
};

/// Result from the full two-stage pipeline
struct TwoStageResult {
    /// Stage 1 coarse screening result
    ScreenResult coarse_result;

    /// Stage 2 full docking score (NaN if not docked)
    float full_dock_score = std::numeric_limits<float>::quiet_NaN();

    /// Stage 2 RMSD to crystal pose (if reference available)
    float rmsd = std::numeric_limits<float>::quiet_NaN();

    /// Rank after Stage 1 (1-based)
    int coarse_rank = 0;

    /// Rank after Stage 2 (1-based, 0 if not docked)
    int full_rank = 0;

    /// How Stage 2 was filled: empty (unset), "callback", or "surrogate".
    /// Never claim this field is a thermodynamic ΔG.
    std::string stage2_kind;
};

/// Total order for Stage-2 scores: finite before NaN; NaNs by coarse_rank.
/// Used by TwoStageScreener::run so pProp-dropped rows cannot break sort.
inline bool stage2_score_less(const TwoStageResult& a, const TwoStageResult& b)
{
    const bool a_nan = !std::isfinite(a.full_dock_score);
    const bool b_nan = !std::isfinite(b.full_dock_score);
    if (a_nan != b_nan) return !a_nan && b_nan;
    if (a_nan) return a.coarse_rank < b.coarse_rank;
    if (a.full_dock_score != b.full_dock_score)
        return a.full_dock_score < b.full_dock_score;
    return a.coarse_rank < b.coarse_rank;
}

/// Callback for Stage 2 docking. Receives the ligand and coarse result,
/// returns the full docking score. This allows plugging in the existing
/// FlexAIDdS GA docking engine without tight coupling.
using FullDockCallback = std::function<float(const ScreenLigand& ligand,
                                              const ScreenResult& coarse)>;

/// Two-stage virtual screening pipeline.
///
/// Usage:
///   TwoStageScreener ts;
///   ts.set_config(cfg);
///   ts.prepare_target(atoms, spheres);
///   ts.set_full_dock_callback(my_dock_fn);
///   auto results = ts.run(ligands);
///
class TwoStageScreener {
public:
    TwoStageScreener();
    ~TwoStageScreener();

    void set_config(const TwoStageConfig& cfg);
    const TwoStageConfig& config() const { return config_; }

    /// Prepare the target (delegates to CoarseScreener)
    void prepare_target(const std::vector<TargetAtom>& atoms,
                        const std::vector<BindingSiteSphere>& spheres);

    /// Load target directly from files
    bool load_target(const std::string& mol2_path, const std::string& cleft_pdb);

    /// Set the Stage 2 docking callback. If not set, only Stage 1 runs.
    void set_full_dock_callback(FullDockCallback cb) { dock_cb_ = std::move(cb); }

    /// Run the full pipeline. Returns results sorted by final score.
    std::vector<TwoStageResult> run(const std::vector<ScreenLigand>& ligands);

    /// Write coarse screening results to CSV
    static void write_csv(const std::string& path,
                          const std::vector<ScreenResult>& results);

    /// Unified Stage-1 + Stage-2 CSV (coarse rank/score + dock CF/RMSD).
    static void write_unified_csv(const std::string& path,
                                  const std::vector<TwoStageResult>& results);

    /// Thin screen receipt (method DOI / Zenodo / counts). Not a claim package.
    static bool write_screen_receipt(const std::string& output_dir,
                                     int n_ligands,
                                     int top_n,
                                     bool stage2_callback,
                                     const std::string& extra_mode = "");

    /// Names of the first min(top_n, size) results (Stage-1 order if undocked).
    static std::vector<std::string> top_names(const std::vector<TwoStageResult>& results,
                                              int top_n);

    /// Access the internal coarse screener (for advanced use)
    const CoarseScreener& coarse_screener() const { return screener_; }

private:
    TwoStageConfig    config_;
    CoarseScreener    screener_;
    FullDockCallback  dock_cb_;
};

} // namespace nrgrank
