// two_stage_screen_main.cpp — CLI driver for NRGRank two-stage screening (P5)
//
// Wires the in-repo NRGRank rigid coarse screen (CoarseScreener /
// TwoStageScreener) as a *pre-filter* in front of the expensive FlexAIDdS GA.
// Stage 1 (this driver) scores an entire ligand library with the fast rigid CF
// and keeps only the top-N candidates; those N are then handed to the full GA
// docking engine (classic `FlexAID`) — an order-of-magnitude throughput win for
// virtual screening with no change to per-dock physics.
//
// This is deliberately a thin driver: it reuses the existing CoarseScreen /
// TwoStageScreen classes rather than reimplementing any scoring, and it does not
// pull in the GA/top.cpp engine (owned elsewhere). Stage 2 is exposed as an
// optional callback hook; by default the driver emits the pruned top-N list
// (names + CSV) that the GA campaign scripts consume.
//
// Usage:
//   flexaid_screen --target target.mol2 --cleft cleft.pdb --ligands lib.mol2 \
//                  [--two-stage] [--coarse-topN N] [--rotations-per-axis R] \
//                  [--no-clash] [--sdf] [--out DIR] [--verbose]
//
// SPDX-License-Identifier: Apache-2.0

#include "CoarseScreen.h"
#include "TwoStageScreen.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace nrgrank;

namespace {

struct CliOptions {
    std::string target;          // target MOL2
    std::string cleft;           // binding-site cleft PDB
    std::string ligands;         // multi-mol2 or SDF library
    std::string out_dir = "screen_results";
    int         top_n   = 100;   // --coarse-topN
    int         rotations_per_axis = 9;
    bool        two_stage = false;   // --two-stage: enable the pipeline mode
    bool        use_clash = true;
    bool        sdf       = false;   // ligand library is SDF (else MOL2)
    bool        verbose   = false;
};

void print_usage(const char* prog) {
    std::printf(
        "FlexAIDdS two-stage screening (NRGRank rigid coarse pre-filter → GA)\n"
        "Usage: %s --target T.mol2 --cleft C.pdb --ligands L.{mol2|sdf} [options]\n\n"
        "Required:\n"
        "  --target FILE          Target receptor MOL2\n"
        "  --cleft FILE           Binding-site cleft/sphere PDB\n"
        "  --ligands FILE         Ligand library (multi-MOL2 by default, or SDF)\n\n"
        "Options:\n"
        "  --two-stage            Enable two-stage MODE (coarse prune + Stage-2 hook)\n"
        "  --coarse-topN N        Keep top-N ligands for the GA (default 100)\n"
        "  --rotations-per-axis R Orientations R^3 in coarse screen (default 9 → 729)\n"
        "  --no-clash             Disable clash filtering in coarse screen\n"
        "  --sdf                  Treat --ligands as SDF (default: MOL2)\n"
        "  --out DIR              Output directory (default screen_results)\n"
        "  --verbose              Verbose progress\n"
        "  -h, --help             This help\n",
        prog);
}

bool parse_args(int argc, char** argv, CliOptions& o) {
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "error: %s requires an argument\n", name);
                return nullptr;
            }
            return argv[++i];
        };
        if (!std::strcmp(a, "--target")) {
            const char* v = need(a); if (!v) return false; o.target = v;
        } else if (!std::strcmp(a, "--cleft")) {
            const char* v = need(a); if (!v) return false; o.cleft = v;
        } else if (!std::strcmp(a, "--ligands")) {
            const char* v = need(a); if (!v) return false; o.ligands = v;
        } else if (!std::strcmp(a, "--coarse-topN")) {
            const char* v = need(a); if (!v) return false; o.top_n = std::atoi(v);
        } else if (!std::strcmp(a, "--rotations-per-axis")) {
            const char* v = need(a); if (!v) return false;
            o.rotations_per_axis = std::atoi(v);
        } else if (!std::strcmp(a, "--out")) {
            const char* v = need(a); if (!v) return false; o.out_dir = v;
        } else if (!std::strcmp(a, "--two-stage")) {
            o.two_stage = true;
        } else if (!std::strcmp(a, "--no-clash")) {
            o.use_clash = false;
        } else if (!std::strcmp(a, "--sdf")) {
            o.sdf = true;
        } else if (!std::strcmp(a, "--verbose")) {
            o.verbose = true;
        } else if (!std::strcmp(a, "-h") || !std::strcmp(a, "--help")) {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            std::fprintf(stderr, "error: unknown argument '%s'\n", a);
            return false;
        }
    }
    if (o.target.empty() || o.cleft.empty() || o.ligands.empty()) {
        std::fprintf(stderr, "error: --target, --cleft and --ligands are required\n");
        return false;
    }
    if (o.top_n < 1) o.top_n = 1;
    if (o.rotations_per_axis < 1) o.rotations_per_axis = 1;
    return true;
}

} // namespace

int main(int argc, char** argv) {
    CliOptions o;
    if (!parse_args(argc, argv, o)) {
        print_usage(argv[0]);
        return 2;
    }

    // ── Load ligand library (reuse existing loaders) ──
    std::vector<ScreenLigand> ligands =
        o.sdf ? CoarseScreener::load_ligands_sdf(o.ligands)
              : CoarseScreener::load_ligands_mol2(o.ligands);
    if (ligands.empty()) {
        std::fprintf(stderr, "error: no ligands parsed from '%s'\n",
                     o.ligands.c_str());
        return 1;
    }
    std::printf("Loaded %zu ligands from %s\n", ligands.size(), o.ligands.c_str());

    // ── Configure and prepare the two-stage screener ──
    TwoStageScreener ts;
    TwoStageConfig cfg;
    cfg.coarse.rotations_per_axis = o.rotations_per_axis;
    cfg.coarse.use_clash          = o.use_clash;
    cfg.coarse.top_n              = o.top_n;
    cfg.top_n                     = o.top_n;
    cfg.output_dir                = o.out_dir;
    cfg.write_coarse_csv          = true;
    cfg.verbose                   = o.verbose;
    ts.set_config(cfg);

    if (!ts.load_target(o.target, o.cleft)) {
        std::fprintf(stderr,
                     "error: failed to prepare target (mol2='%s', cleft='%s')\n",
                     o.target.c_str(), o.cleft.c_str());
        return 1;
    }
    std::printf("Target prepared: %zu anchor points\n",
                ts.coarse_screener().num_anchors());

    // ── Stage 2 hook ──
    // In --two-stage mode we leave the GA docking callback UNSET here: coupling
    // the full FlexAID GA (top.cpp) into a callback belongs to the campaign
    // launcher that owns the engine. This driver's job is the Stage-1 prune; it
    // emits the top-N candidate list the GA then docks. Setting a real callback
    // is a one-liner (ts.set_full_dock_callback(fn)) once wired.
    if (o.two_stage) {
        std::printf("Two-stage MODE: Stage-1 coarse prune → top %d handed to GA\n",
                    o.top_n);
    }

    // ── Run ──
    std::vector<TwoStageResult> results = ts.run(ligands);
    if (results.empty()) {
        std::fprintf(stderr, "error: screening produced no results\n");
        return 1;
    }

    // ── Emit pruned top-N candidate list (fed to the GA) ──
    const int n_keep = std::min<int>(o.top_n, static_cast<int>(results.size()));
    std::filesystem::create_directories(o.out_dir);
    const std::string top_path = o.out_dir + "/stage1_topN.txt";
    std::ofstream top_out(top_path);
    if (top_out.is_open()) {
        top_out << "# rank\tname\tcoarse_score\n";
        for (int i = 0; i < n_keep; ++i) {
            const auto& r = results[i];
            top_out << (i + 1) << '\t' << r.coarse_result.name << '\t'
                    << r.coarse_result.score << '\n';
        }
    }

    TwoStageScreener::write_unified_csv(o.out_dir + "/unified.csv", results);
    TwoStageScreener::write_screen_receipt(
        o.out_dir, static_cast<int>(ligands.size()), o.top_n,
        /*stage2_callback=*/false,
        o.two_stage ? "flexaid_screen-two-stage" : "flexaid_screen");

    std::printf("Screening complete: %zu ligands scored, top %d kept for GA.\n"
                "  coarse CSV : %s/coarse_screen.csv\n"
                "  unified CSV: %s/unified.csv\n"
                "  receipt    : %s/RUN_RECEIPT.json\n"
                "  top-N list : %s\n",
                results.size(), n_keep, o.out_dir.c_str(), o.out_dir.c_str(),
                o.out_dir.c_str(), top_path.c_str());
    std::printf("Best candidate: %s (coarse CF = %.4f)\n",
                results[0].coarse_result.name.c_str(),
                results[0].coarse_result.score);
    return 0;
}
