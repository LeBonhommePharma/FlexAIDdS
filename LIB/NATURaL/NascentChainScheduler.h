// NascentChainScheduler.h — Checkpoint scheduling for cotranslational docking
//
// Wraps natural::build_parallel_growth_schedule(tracks) into an iterator that emits
// one Checkpoint per chain length L_k spaced at a user-configurable interval. Pure
// scheduling logic — no GA dependency. The DualAssemblyRunner records outcomes back
// against each checkpoint as it advances.
//
// References: see docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md §5.1
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#pragma once

#include "NATURaLDualAssembly.h"

#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace natural {

// ─── Checkpoint ──────────────────────────────────────────────────────────────
struct Checkpoint {
    int               idx              = -1;     // 0-based checkpoint index
    int               L_k              = 0;      // chain length (aa or nt)
    double            t_arrival_s      = 0.0;    // real-time since initiation
    GrowthProcess     process          = GrowthProcess::Translation;
    DockingRolePolicy role_policy      = DockingRolePolicy::ProtofibrilAsTarget;
    int               track_index      = -1;
    std::string       track_name;
    bool              in_tunnel        = false;  // L_k ≤ tunnel_length
    bool              chaperone_shielded = false; // L_k - tunnel_length < 6
    // True for transcription tracks in eukaryotic cells (no direct protofibril encounter
    // during nuclear synthesis). Mirrors InVivoAssemblyTrack::direct_encounter_allowed.
    bool              direct_encounter_allowed = true;
};

// ─── CheckpointOutcome ───────────────────────────────────────────────────────
struct CheckpointOutcome {
    // Shannon entropy over backend-projected pose coordinates, not intrinsic
    // conformational entropy of the isolated molecule.
    double H_A_nats   = std::numeric_limits<double>::infinity();
    double H_B_nats   = std::numeric_limits<double>::infinity();
    double dG_A_kcal  = std::numeric_limits<double>::quiet_NaN();
    double dG_B_kcal  = std::numeric_limits<double>::quiet_NaN();
    // Pose-entropy role hint. This is a reproducible heuristic label, not a
    // proof that one molecule is intrinsically "the" receptor.
    // 'A' = Sim A pose ensemble lower entropy, 'B' = Sim B pose ensemble lower
    // entropy, '?' = deferred.
    char   tl_primary  = '?';
    double tl_weight   = 0.0;     // confidence ∈ [0,1]
    bool   tl_deferred = true;
    std::string tl_basis = "pose_entropy_heuristic";
    double p_elong     = std::numeric_limits<double>::quiet_NaN();
    double dG_elong    = std::numeric_limits<double>::quiet_NaN();
    bool   sim_c_evaluated = false;
    bool   sim_c_gated_in = false;
    int    protofibril_state_index = 0;
    bool   protofibril_structure_updated = false;
};

// ─── NascentChainScheduler ───────────────────────────────────────────────────
class NascentChainScheduler {
public:
    // tunnel_length is taken from each InVivoAssemblyTrack (process-native units —
    // aa for ribosome, nt for RNAP) rather than from a global default. The
    // chaperone_skip_residues parameter delays the first checkpoint past the NAC/SSB
    // binding region.
    NascentChainScheduler(std::vector<InVivoAssemblyTrack> tracks,
                          int                              checkpoint_interval,
                          int                              chaperone_skip_residues = 6);

    bool       has_next() const noexcept;
    Checkpoint next();                                        // advance the cursor
    void       record(const Checkpoint& ck, CheckpointOutcome out);

    const std::vector<std::pair<Checkpoint, CheckpointOutcome>>& history() const noexcept
    { return history_; }

    // Pre-computed schedule (read-only). Pre-computation is deterministic and
    // depends only on the input tracks + interval.
    const std::vector<Checkpoint>& schedule() const noexcept { return schedule_; }

    int  checkpoint_interval() const noexcept { return interval_; }
    int  cursor()              const noexcept { return cursor_; }
    void reset()                              { cursor_ = 0; history_.clear(); }

private:
    std::vector<InVivoAssemblyTrack> tracks_;
    int                              interval_;
    int                              chaperone_skip_;
    std::vector<Checkpoint>          schedule_;
    int                              cursor_ = 0;
    std::vector<std::pair<Checkpoint, CheckpointOutcome>> history_;
};

} // namespace natural
