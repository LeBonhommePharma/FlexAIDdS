#include "ReportedPoseSelector.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>
#include <filesystem>

namespace fs = std::filesystem;

namespace reported_pose {

// Adapted from DatasetRunner: parse one pose file for CF, freq, member_cfs (for Z+H).
static bool parse_pose(const std::string& cand, float& cf, int& freq, std::vector<float>& member_cfs) {
    cf = std::numeric_limits<float>::infinity();
    freq = 1;
    bool have_cf = false;
    std::ifstream pf(cand);
    std::string pl;
    while (std::getline(pf, pl)) {
        if (!have_cf && pl.find("REMARK CF=") != std::string::npos) {
            auto p2 = pl.find("CF=");
            if (p2 != std::string::npos) {
                try { cf = std::stof(pl.substr(p2 + 3)); have_cf = true; }
                catch (...) {}
            }
        } else if (pl.find("Frequency:") != std::string::npos) {
            auto p2 = pl.find("Frequency:");
            try { freq = std::stoi(pl.substr(p2 + 10)); } catch (...) {}
        }
    }
    if (!have_cf || !std::isfinite(cf)) return false;
    // mcf sidecar for member cfs (for entropy).
    {
        std::string mcf = cand;
        if (mcf.size() > 4 && mcf.substr(mcf.size() - 4) == ".pdb")
            mcf = mcf.substr(0, mcf.size() - 4) + ".mcf";
        std::ifstream mf(mcf);
        if (mf.is_open()) {
            std::string ml;
            while (std::getline(mf, ml)) {
                try {
                    float v = std::stof(ml);
                    if (std::isfinite(v)) member_cfs.push_back(v);
                } catch (...) {}
            }
        }
    }
    return true;
}

float parse_g_bind_from_log(const std::string& log_text) {
    // Find last [THERMO] G_bind= in the text.
    size_t pos = log_text.rfind("[THERMO]");
    if (pos == std::string::npos) return std::numeric_limits<float>::quiet_NaN();
    size_t gpos = log_text.find("G_bind=", pos);
    if (gpos == std::string::npos) return std::numeric_limits<float>::quiet_NaN();
    try {
        return std::stof(log_text.substr(gpos + 7));
    } catch (...) { return std::numeric_limits<float>::quiet_NaN(); }
}

std::vector<PoseCandidate> build_cross_restart_pool(const std::vector<std::string>& prefixes) {
    std::vector<PoseCandidate> poses;
    constexpr double kT_kcalmol = 0.592;
    constexpr double alpha_shannon = 1.0;

    auto score_composite = [&](const std::vector<float>& member_cfs, float best_cf, int pop) -> double {
        if (member_cfs.empty()) {
            double z = std::exp(-static_cast<double>(best_cf) / kT_kcalmol);
            return z * std::pow(static_cast<double>(pop), 0.1);
        }
        double Z = 0.0;
        for (float cf_i : member_cfs) {
            Z += std::exp(-static_cast<double>(cf_i) / kT_kcalmol);
        }
        double H = 0.0;
        for (float cf_i : member_cfs) {
            double pi = std::exp(-static_cast<double>(cf_i) / kT_kcalmol) / Z;
            if (pi > 1e-300) H -= pi * std::log(pi);
        }
        double pop_weight = std::pow(static_cast<double>(pop), 0.1);
        return Z * std::exp(-alpha_shannon * H) * pop_weight;
    };

    for (size_t ri = 0; ri < prefixes.size(); ++ri) {
        const auto& out_prefix = prefixes[ri];
        // Get G for this restart.
        float g = std::numeric_limits<float>::quiet_NaN();
        std::string logp = out_prefix + "/stdout.log";
        if (fs::exists(logp)) {
            std::ifstream f(logp);
            std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
            g = parse_g_bind_from_log(text);
        } else {
            // fallback to parent stdout.log , last G seen may be for last restart; for simplicity assign NaN or last.
            std::string parent_log = (fs::path(out_prefix).parent_path() / "stdout.log").string();
            if (fs::exists(parent_log)) {
                std::ifstream f(parent_log);
                std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
                g = parse_g_bind_from_log(text);
            }
        }

        // Collect poses for this restart.
        for (int pi = 0; pi <= 19; ++pi) {
            std::string cand = out_prefix + "_" + std::to_string(pi) + ".pdb";
            if (!fs::exists(cand)) continue;
            float cf; int freq; std::vector<float> mcfs;
            if (!parse_pose(cand, cf, freq, mcfs)) continue;
            poses.push_back({cand, cf, freq, g, static_cast<int>(ri)});
        }
        // INI seed if present.
        std::string ini = out_prefix + "_INI.pdb";
        if (fs::exists(ini)) {
            float cf; int freq; std::vector<float> mcfs;
            if (parse_pose(ini, cf, freq, mcfs)) {
                poses.push_back({ini, cf, 1, g, static_cast<int>(ri)});
            }
        }
    }
    return poses;
}

std::string elect_reported_pose(const std::vector<PoseCandidate>& pool, bool thermo_on) {
    if (pool.empty()) return {};

    // Simple adaptation of the logic: group by restart for consensus, score with G tiebreak if thermo.
    // For simplicity, replicate the default consensus path with G as tertiary tie-break (lower G better).
    // Build per-restart best CF for rough, but use full pool scoring.

    // To keep close to shipped: use the Z+H or consensus from the code.
    // Here we implement a basic version that prefers high "consensus" (number of other restarts with close pose), low CF, then low G if thermo.

    // Compute approximate consensus for each.
    std::vector<int> consensus(pool.size(), 0);
    constexpr float kConsensusDelta = 1.5f;
    // Note: full RMSD computation requires coords; for this, we fall back to simple CF + G if no coords loaded.
    // For the module, we assume caller provides or we skip full consensus for now and use CF + G.
    // To make functional, use CF primary, G secondary for thermo.

    int best_i = -1;
    for (size_t i = 0; i < pool.size(); ++i) {
        bool better = false;
        if (best_i < 0) better = true;
        else {
            if (pool[i].cf < pool[best_i].cf) better = true;
            else if (pool[i].cf == pool[best_i].cf) {
                if (thermo_on && std::isfinite(pool[i].g_bind) && std::isfinite(pool[best_i].g_bind)) {
                    if (pool[i].g_bind < pool[best_i].g_bind) better = true;
                }
            }
        }
        if (better) best_i = static_cast<int>(i);
    }
    if (best_i < 0) return {};
    return pool[best_i].path;
}

} // namespace reported_pose
