#pragma once

#include <string>
#include <vector>
#include <limits>

namespace reported_pose {

// Candidate pose with its properties from one restart.
struct PoseCandidate {
    std::string path;
    float cf = std::numeric_limits<float>::infinity();
    int freq = 1;
    float g_bind = std::numeric_limits<float>::quiet_NaN();  // lower is better for binding
    int restart_id = -1;
};

// Build the cross-restart pool of pose candidates from the restart prefixes.
// For each restart, parses its poses (_0.._N.pdb) and CF from REMARK, and assigns G_bind from the restart's [THERMO] log if available.
std::vector<PoseCandidate> build_cross_restart_pool(const std::vector<std::string>& prefixes);

// Elect the reported pose path using the shipped logic:
// freq-gated / consensus (Z+H or CF) , with optional G_bind as tie-break (min G) over full pool when thermo_on.
std::string elect_reported_pose(const std::vector<PoseCandidate>& pool, bool thermo_on);

// Helper to parse a single G_bind value from log text containing [THERMO] G_bind=...
float parse_g_bind_from_log(const std::string& log_text);

} // namespace reported_pose
