// NascentChainScheduler.cpp — see header for design.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include "NascentChainScheduler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace natural {

NascentChainScheduler::NascentChainScheduler(
    std::vector<InVivoAssemblyTrack> tracks,
    int                              checkpoint_interval,
    int                              chaperone_skip_residues)
    : tracks_(std::move(tracks)),
      interval_(checkpoint_interval),
      chaperone_skip_(chaperone_skip_residues)
{
    if (interval_ <= 0)
        throw std::invalid_argument("NascentChainScheduler: checkpoint_interval must be > 0");

    for (int ti = 0; ti < static_cast<int>(tracks_.size()); ++ti) {
        const auto& tr = tracks_[ti];
        if (tr.chain_length <= 0) continue;
        if (!std::isfinite(tr.mean_elongation_rate) || tr.mean_elongation_rate <= 0.0)
            throw std::invalid_argument("NascentChainScheduler: positive elongation rate required");

        const double dwell = 1.0 / tr.mean_elongation_rate;
        const double init_wait = (std::isfinite(tr.initiation_rate) && tr.initiation_rate > 1e-12)
            ? 1.0 / tr.initiation_rate
            : 0.0;
        // tunnel measured in process-native units (aa for ribosome, nt for RNAP).
        const double tunnel_units = std::max(0.0, tr.tunnel_length);

        // First checkpoint at L = tunnel + 6 (first solvent-exposed residue past the
        // NAC/SSB region), then every interval_ units up to chain_length.
        const int first_L = static_cast<int>(std::floor(tunnel_units)) + chaperone_skip_;
        for (int L = first_L; L <= tr.chain_length; L += interval_) {
            Checkpoint ck;
            ck.idx              = static_cast<int>(schedule_.size());
            ck.L_k              = L;
            ck.t_arrival_s      = init_wait + L * dwell;
            ck.process          = tr.process;
            ck.role_policy      = tr.role_policy;
            ck.track_index      = ti;
            ck.track_name       = tr.name;
            ck.in_tunnel        = static_cast<double>(L) <= tunnel_units;
            ck.chaperone_shielded = (static_cast<double>(L) - tunnel_units) < chaperone_skip_;
            ck.direct_encounter_allowed = tr.direct_encounter_allowed;
            schedule_.push_back(std::move(ck));
        }
    }

    // Sort by (t_arrival, track_index) so the parallel timeline is monotone.
    std::stable_sort(schedule_.begin(), schedule_.end(),
                     [](const Checkpoint& a, const Checkpoint& b) {
                         if (a.t_arrival_s != b.t_arrival_s) return a.t_arrival_s < b.t_arrival_s;
                         return a.track_index < b.track_index;
                     });
    // Re-number idx in sorted order so consumers see a stable monotone sequence.
    for (int i = 0; i < static_cast<int>(schedule_.size()); ++i)
        schedule_[i].idx = i;
}

bool NascentChainScheduler::has_next() const noexcept {
    return cursor_ < static_cast<int>(schedule_.size());
}

Checkpoint NascentChainScheduler::next() {
    if (!has_next())
        throw std::out_of_range("NascentChainScheduler: schedule exhausted");
    return schedule_[cursor_++];
}

void NascentChainScheduler::record(const Checkpoint& ck, CheckpointOutcome out) {
    history_.emplace_back(ck, std::move(out));
}

} // namespace natural
