// DualAssemblyRunner.cpp — see header for design.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include "DualAssemblyRunner.h"

#include "../ShannonThermoStack/ShannonThermoStack.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace natural {

namespace {
constexpr int kShannonBins = 20;   // matches shannon_thermo::DEFAULT_HIST_BINS

std::vector<InVivoAssemblyTrack> make_tracks_from_config(const DualAssemblyConfig& cfg) {
    const int aa_len = static_cast<int>(cfg.sequence_fasta.size());
    const int nt_len = (cfg.transcript_nt > 0) ? cfg.transcript_nt : aa_len * 3;
    return make_human_protofibril_tracks(nt_len, aa_len, cfg.include_reciprocal_controls);
}
} // namespace

// ─── ctor ────────────────────────────────────────────────────────────────────
DualAssemblyRunner::DualAssemblyRunner(DualAssemblyConfig cfg,
                                       SimAFn             sim_a,
                                       SimBFn             sim_b,
                                       SimCFn             sim_c,
                                       TruncateFn         truncate,
                                       ProtofibrilAdvanceFn advance_protofibril)
    : cfg_(std::move(cfg)),
      scheduler_(make_tracks_from_config(cfg_), cfg_.checkpoint_interval),
      oracle_(cfg_.temperature_K, cfg_.acceptance_threshold),
      sim_a_(std::move(sim_a)),
      sim_b_(std::move(sim_b)),
      sim_c_(std::move(sim_c)),
      truncate_(std::move(truncate)),
      advance_protofibril_(std::move(advance_protofibril)),
      current_protofibril_pdb_(cfg_.protofibril_pdb)
{
    if (cfg_.protofibril_pdb.empty())
        throw std::invalid_argument("DualAssemblyRunner: protofibril_pdb is required");
    if (cfg_.sequence_fasta.empty())
        throw std::invalid_argument("DualAssemblyRunner: sequence_fasta is required");
    if (cfg_.sim_c_enabled && cfg_.monomer_pdb.empty())
        throw std::invalid_argument("DualAssemblyRunner: monomer_pdb required when sim_c_enabled");
    if (!sim_a_)
        throw std::invalid_argument("DualAssemblyRunner: sim_a callback is required");
    if (!truncate_)
        throw std::invalid_argument("DualAssemblyRunner: truncate callback is required");
    if (cfg_.sim_c_enabled && !sim_c_)
        throw std::invalid_argument("DualAssemblyRunner: sim_c callback required when sim_c_enabled");
}

// ─── Pose-entropy role discriminator ─────────────────────────────────────────
// Implements docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md §2.
//
// tl_primary convention:
//   'A' = Sim A's pose-coordinate ensemble is lower entropy.
//   'B' = Sim B's pose-coordinate ensemble is lower entropy.
//   '?' = deferred — neither pose ensemble has fallen below the soft threshold.
//
// This is intentionally a role hint derived from docking pose collapse. It is not
// intrinsic conformational entropy of the isolated chain/protofibril, so the CSV
// carries tl_basis=pose_entropy_heuristic.
void DualAssemblyRunner::assign_tl(double             H_A_nats,
                                    double             H_B_nats,
                                    char               prev_tl,
                                    CheckpointOutcome& out) noexcept
{
    using shannon_thermo::kHSC_soft_nats;
    using shannon_thermo::kHSC_hard_nats;

    const double H_lower = std::min(H_A_nats, H_B_nats);
    const bool   A_lower = (H_A_nats <= H_B_nats);

    if (!std::isfinite(H_lower)) {
        // Both H are infinite — neither sim contributed (most likely Sim B skipped
        // and Sim A failed). Retain previous assignment with zero weight.
        out.tl_primary  = prev_tl;
        out.tl_weight   = 0.0;
        out.tl_deferred = true;
        return;
    }

    if (H_lower < kHSC_hard_nats) {
        // Hard regime: the lower-H system is target with full confidence.
        out.tl_primary  = A_lower ? 'A' : 'B';
        out.tl_weight   = 1.0;
        out.tl_deferred = false;
        return;
    }

    if (H_lower >= kHSC_soft_nats) {
        // Deferred regime: inherit previous assignment with zero weight.
        out.tl_primary  = prev_tl;
        out.tl_weight   = 0.0;
        out.tl_deferred = true;
        return;
    }

    // Soft regime: w ∈ (0,1) linear in (soft − H_lower)/(soft − hard).
    const double w = (kHSC_soft_nats - H_lower) / (kHSC_soft_nats - kHSC_hard_nats);
    out.tl_primary  = A_lower ? 'A' : 'B';
    out.tl_weight   = std::clamp(w, 0.0, 1.0);
    out.tl_deferred = false;
}

// ─── Helpers ─────────────────────────────────────────────────────────────────
double DualAssemblyRunner::shannon_entropy_nats(const std::vector<double>& values) {
    if (values.empty()) return std::numeric_limits<double>::infinity();
    return shannon_thermo::compute_shannon_entropy(values, kShannonBins);
}

double DualAssemblyRunner::free_energy_kcal(const statmech::StatMechEngine& engine) {
    if (engine.size() == 0) return std::numeric_limits<double>::quiet_NaN();
    return engine.compute().free_energy;
}

// ─── CSV ─────────────────────────────────────────────────────────────────────
void DualAssemblyRunner::write_csv_header(const std::string& path) const {
    std::ofstream out(path, std::ios::trunc);
    if (!out) throw std::runtime_error("DualAssemblyRunner: cannot write " + path);
    out << "checkpoint_idx,L_k,process,track_name,role_policy,t_arrival_s,in_tunnel,"
           "H_A_nats,H_B_nats,dG_A_kcal,dG_B_kcal,"
           "tl_primary,tl_weight,tl_deferred,tl_basis,"
           "p_elong,dG_elong_kcal,sim_c_evaluated,sim_c_gated_in,"
           "protofibril_state_idx,protofibril_structure_updated\n";
}

void DualAssemblyRunner::write_csv_row(const std::string& path,
                                        const Checkpoint& ck,
                                        const CheckpointOutcome& out) const
{
    std::ofstream f(path, std::ios::app);
    if (!f) throw std::runtime_error("DualAssemblyRunner: cannot append to " + path);
    auto fmt = [](double x) -> std::string {
        if (!std::isfinite(x)) return "NaN";
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.6g", x);
        return buf;
    };
    const char* proc =
        (ck.process == GrowthProcess::Translation) ? "Translation" : "Transcription";
    const char* role =
        (ck.role_policy == DockingRolePolicy::ProtofibrilAsTarget)
            ? "ProtofibrilAsTarget" : "ReciprocalControl";
    f << ck.idx << ',' << ck.L_k << ',' << proc << ',' << ck.track_name << ','
      << role << ',' << fmt(ck.t_arrival_s) << ',' << (ck.in_tunnel ? "1" : "0") << ','
      << fmt(out.H_A_nats) << ',' << fmt(out.H_B_nats) << ','
      << fmt(out.dG_A_kcal) << ',' << fmt(out.dG_B_kcal) << ','
      << out.tl_primary << ',' << fmt(out.tl_weight) << ',' << (out.tl_deferred ? "1" : "0") << ','
      << out.tl_basis << ','
      << fmt(out.p_elong) << ',' << fmt(out.dG_elong) << ','
      << (out.sim_c_evaluated ? "1" : "0") << ','
      << (out.sim_c_gated_in ? "1" : "0") << ','
      << out.protofibril_state_index << ','
      << (out.protofibril_structure_updated ? "1" : "0") << '\n';
}

// ─── run ─────────────────────────────────────────────────────────────────────
std::vector<std::pair<Checkpoint, CheckpointOutcome>> DualAssemblyRunner::run() {
    write_csv_header(cfg_.output_csv);

    while (scheduler_.has_next()) {
        Checkpoint ck = scheduler_.next();
        CheckpointOutcome out;
        out.protofibril_state_index = protofibril_state_index_;
        TrackState& state = track_state_[ck.track_index];

        // Skip transcription tracks for direct docking (eukaryotic decoupling). The
        // transcription clock remains in the schedule as a synchronisation reference
        // but does not run Sim A.
        if (!ck.direct_encounter_allowed) {
            out.tl_deferred = true;
            scheduler_.record(ck, out);
            write_csv_row(cfg_.output_csv, ck, out);
            continue;
        }

        // Tunnel / chaperone gate.
        if (ck.in_tunnel || ck.chaperone_shielded) {
            out.tl_deferred = true;
            scheduler_.record(ck, out);
            write_csv_row(cfg_.output_csv, ck, out);
            continue;
        }

        // Generate the truncated nascent-chain PDB for this checkpoint.
        const std::string chain_pdb = truncate_(cfg_.sequence_fasta, ck.L_k,
                                                cfg_.nascent_pdb_dir);

        // ── Sim A ───────────────────────────────────────────────────────────
        // The reciprocal track must exercise the actual backend with swapped
        // target/ligand order; otherwise it is only a label in the CSV.
        const bool reciprocal =
            (ck.role_policy == DockingRolePolicy::ReciprocalControl);
        const std::string protofibril_pdb = current_protofibril_pdb_;
        GAResult a = reciprocal
            ? sim_a_(chain_pdb, protofibril_pdb, ck.L_k, cfg_.temperature_K)
            : sim_a_(protofibril_pdb, chain_pdb, ck.L_k, cfg_.temperature_K);
        out.H_A_nats   = shannon_entropy_nats(a.pose_rmsds);
        out.dG_A_kcal  = free_energy_kcal(a.engine);

        // ── Sim B (only if chain has collapsed at the prior checkpoint) ─────
        if (state.H_prev_chain_nats < shannon_thermo::kHSC_soft_nats && sim_b_ && !cfg_.monomer_pdb.empty()) {
            GAResult b = sim_b_(chain_pdb, cfg_.monomer_pdb, ck.L_k, cfg_.temperature_K);
            out.H_B_nats  = shannon_entropy_nats(b.pose_rmsds);
            out.dG_B_kcal = free_energy_kcal(b.engine);
        }

        // ── Sim C (every sim_c_interval checkpoints) ────────────────────────
        if (cfg_.sim_c_enabled && sim_c_) {
            ++state.checkpoints_since_sim_c;
            if (state.checkpoints_since_sim_c >= cfg_.sim_c_interval) {
                GAResult c = sim_c_(current_protofibril_pdb_, cfg_.monomer_pdb,
                                    cfg_.temperature_K);
                out.sim_c_evaluated = true;
                if (c.engine.size() > 0) {
                    auto dec   = oracle_.gate(c.engine, cfg_.monomer_conc_M);
                    out.p_elong = dec.p_elong;
                    out.dG_elong = dec.dG_elong;
                    out.sim_c_gated_in = dec.gated_in;
                    if (dec.gated_in) {
                        ++protofibril_state_index_;
                        out.protofibril_state_index = protofibril_state_index_;
                        if (advance_protofibril_) {
                            std::string next_protofibril = advance_protofibril_(
                                current_protofibril_pdb_,
                                cfg_.monomer_pdb,
                                protofibril_state_index_);
                            if (next_protofibril.empty())
                                throw std::runtime_error("DualAssemblyRunner: protofibril advance returned an empty path");
                            current_protofibril_pdb_ = std::move(next_protofibril);
                            out.protofibril_structure_updated = true;
                        }
                    }
                }
                state.checkpoints_since_sim_c = 0;
            }
        }

        // ── T/L discriminator ───────────────────────────────────────────────
        assign_tl(out.H_A_nats, out.H_B_nats, state.tl_prev, out);

        scheduler_.record(ck, out);
        write_csv_row(cfg_.output_csv, ck, out);

        state.H_prev_chain_nats = out.H_A_nats;
        if (!out.tl_deferred) state.tl_prev = out.tl_primary;
    }

    return scheduler_.history();
}

} // namespace natural
