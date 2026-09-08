// PoseHelixThermoRewrite.cpp — experimental 2014 NATURAL pose→helix ΔH/ΔS glue
//
// See PoseHelixThermoRewrite.h and docs/POSE_HELIX_THERMO_REWRITE.md.
// Apache-2.0 © 2026 Le Bonhomme Pharma
#include "PoseHelixThermoRewrite.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <utility>

namespace natural {

constexpr double kCalPerKcal = 1000.0;

namespace {

constexpr double kPi = 3.14159265358979323846;

bool finite3(double x, double y, double z) {
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

double dist_A(const ReceptorNtCoord& nt, const LigandAtomCoord& lig) {
    const double dx = nt.x - lig.x;
    const double dy = nt.y - lig.y;
    const double dz = nt.z - lig.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double vec_norm(double x, double y, double z) {
    return std::sqrt(x * x + y * y + z * z);
}

// Angle between ring normals in degrees, folded into [0, 90] so parallel and
// antiparallel stacks both count. Zero/degenerate normals → 0° (distance only).
double stack_angle_deg(const ReceptorNtCoord& nt, const LigandAtomCoord& lig) {
    const double n1 = vec_norm(nt.nx, nt.ny, nt.nz);
    const double n2 = vec_norm(lig.nx, lig.ny, lig.nz);
    if (n1 < 1e-12 || n2 < 1e-12) return 0.0;
    double cosang = (nt.nx * lig.nx + nt.ny * lig.ny + nt.nz * lig.nz) / (n1 * n2);
    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    const double acute = std::fabs(cosang); // parallel or antiparallel
    return std::acos(acute) * 180.0 / kPi;
}

bool in_inclusive(int x, int a, int b) {
    if (a > b) std::swap(a, b);
    return x >= a && x <= b;
}

bool nt_in_loop(const HelixSegment& h, int nt) {
    return in_inclusive(nt, h.loop_begin, h.loop_end);
}

bool nt_in_stem(const HelixSegment& h, int nt) {
    return in_inclusive(nt, h.stem_i_begin, h.stem_i_end) ||
           in_inclusive(nt, h.stem_j_begin, h.stem_j_end);
}

bool nt_in_helix(const HelixSegment& h, int nt) {
    return nt_in_loop(h, nt) || nt_in_stem(h, nt);
}

bool nt_in_any_helix(const std::vector<HelixSegment>& helices, int nt) {
    for (const auto& h : helices) {
        if (nt_in_helix(h, nt)) return true;
    }
    return false;
}

bool valid_helix(const HelixSegment& h) {
    return h.stem_i_end >= h.stem_i_begin ||
           h.stem_j_end >= h.stem_j_begin ||
           h.loop_end >= h.loop_begin;
}

struct HelixContactCount {
    int n_hbond = 0;
    int n_stack = 0;
    std::vector<int> touched;
};

HelixContactCount count_contacts(const HelixSegment& helix,
                                 const std::vector<ReceptorNtCoord>& rna_nts,
                                 const LigandPose& pose,
                                 const PoseHelixRewriteConfig& cfg) {
    HelixContactCount out;
    std::vector<int> hbond_nts;
    std::vector<int> stack_nts;
    for (const auto& nt : rna_nts) {
        if (nt.nt_index < 0) continue;
        if (!finite3(nt.x, nt.y, nt.z)) continue;
        const bool in_loop = nt_in_loop(helix, nt.nt_index);
        const bool in_stem = nt_in_stem(helix, nt.nt_index);
        if (!in_loop && !in_stem) continue;
        bool hit = false;
        for (const auto& lig : pose.atoms) {
            if (!finite3(lig.x, lig.y, lig.z)) continue;
            const double d = dist_A(nt, lig);
            if (in_loop && nt.is_hbond_partner && lig.is_hbond_partner &&
                d <= cfg.hbond_max_A) {
                hit = true;
                if (std::find(hbond_nts.begin(), hbond_nts.end(), nt.nt_index) ==
                    hbond_nts.end()) {
                    hbond_nts.push_back(nt.nt_index);
                }
            }
            if (in_stem && nt.is_aromatic && lig.is_aromatic &&
                d <= cfg.stack_max_A &&
                stack_angle_deg(nt, lig) <= cfg.stack_angle_max_deg) {
                hit = true;
                if (std::find(stack_nts.begin(), stack_nts.end(), nt.nt_index) ==
                    stack_nts.end()) {
                    stack_nts.push_back(nt.nt_index);
                }
            }
        }
        if (hit) out.touched.push_back(nt.nt_index);
    }
    out.n_hbond = static_cast<int>(hbond_nts.size());
    out.n_stack = static_cast<int>(stack_nts.size());
    std::sort(out.touched.begin(), out.touched.end());
    return out;
}

HelixThermoPatch make_patch(const std::string& label,
                            const HelixContactCount& c,
                            const PoseHelixRewriteConfig& cfg,
                            int pose_rank,
                            double weight) {
    HelixThermoPatch p;
    p.helix_label = label;
    p.delta_H_kcal = static_cast<double>(c.n_stack) * cfg.dH_per_stack_kcal;
    p.delta_S_cal_per_K = static_cast<double>(c.n_hbond) * cfg.dS_per_hbond_cal_K;
    p.delta_G_kcal = helix_thermo_delta_g_kcal(p.delta_H_kcal, p.delta_S_cal_per_K, cfg.T_K);
    p.touched_nt = c.touched;
    p.pose_rank = pose_rank;
    p.pose_boltzmann_w = weight;
    return p;
}

void add_touched(std::vector<int>& dst, const std::vector<int>& src) {
    for (int nt : src) {
        if (std::find(dst.begin(), dst.end(), nt) == dst.end()) dst.push_back(nt);
    }
    std::sort(dst.begin(), dst.end());
}

} // namespace

double helix_thermo_delta_g_kcal(double delta_H_kcal,
                                 double delta_S_cal_per_K,
                                 double T_K) {
    if (!std::isfinite(delta_H_kcal) || !std::isfinite(delta_S_cal_per_K) ||
        !std::isfinite(T_K)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return delta_H_kcal - T_K * (delta_S_cal_per_K / kCalPerKcal);
}

std::vector<double> boltzmann_weights_from_scores(const std::vector<double>& scores_kcal,
                                                   double T_K) {
    std::vector<double> w(scores_kcal.size(), 0.0);
    if (scores_kcal.empty()) return w;
    const double T = (std::isfinite(T_K) && T_K > 1e-12) ? T_K : 298.15;
    const double kT = kB_pose_helix_kcal * T;
    if (!(kT > 0.0) || !std::isfinite(kT)) {
        const double eq = 1.0 / static_cast<double>(scores_kcal.size());
        std::fill(w.begin(), w.end(), eq);
        return w;
    }

    std::vector<double> neg_betaE;
    neg_betaE.reserve(scores_kcal.size());
    for (double E : scores_kcal) {
        if (!std::isfinite(E)) {
            neg_betaE.push_back(-std::numeric_limits<double>::infinity());
            continue;
        }
        neg_betaE.push_back(-E / kT);
    }

    double m = -std::numeric_limits<double>::infinity();
    for (double v : neg_betaE) {
        if (v > m) m = v;
    }
    if (!std::isfinite(m)) {
        const double eq = 1.0 / static_cast<double>(scores_kcal.size());
        std::fill(w.begin(), w.end(), eq);
        return w;
    }

    double sum = 0.0;
    std::vector<double> e(neg_betaE.size(), 0.0);
    for (size_t i = 0; i < neg_betaE.size(); ++i) {
        if (!std::isfinite(neg_betaE[i])) {
            e[i] = 0.0;
            continue;
        }
        e[i] = std::exp(neg_betaE[i] - m);
        sum += e[i];
    }
    if (!(sum > 0.0) || !std::isfinite(sum)) {
        const double eq = 1.0 / static_cast<double>(scores_kcal.size());
        std::fill(w.begin(), w.end(), eq);
        return w;
    }
    for (size_t i = 0; i < w.size(); ++i) w[i] = e[i] / sum;
    return w;
}

PoseHelixRewriteResult rewrite_helices_from_poses(
    const std::vector<HelixSegment>& helices,
    const std::vector<ReceptorNtCoord>& rna_nts,
    const std::vector<LigandPose>& poses,
    const PoseHelixRewriteConfig& cfg)
{
    PoseHelixRewriteResult result;
    result.applied = true;

    if (helices.empty() || poses.empty() || cfg.top_k_poses <= 0) {
        return result;
    }

    std::vector<LigandPose> ranked = poses;
    std::sort(ranked.begin(), ranked.end(),
              [](const LigandPose& a, const LigandPose& b) {
                  if (a.rank != b.rank) return a.rank < b.rank;
                  return a.score_kcal < b.score_kcal;
              });
    const int k = std::min(cfg.top_k_poses, static_cast<int>(ranked.size()));
    ranked.resize(static_cast<size_t>(k));

    std::vector<double> scores;
    scores.reserve(ranked.size());
    for (const auto& p : ranked) scores.push_back(p.score_kcal);
    const std::vector<double> weights = boltzmann_weights_from_scores(scores, cfg.T_K);

    struct MixAcc {
        double dH = 0.0;
        double dS = 0.0;
        std::vector<int> touched;
        double w_sum = 0.0;
    };
    std::unordered_map<std::string, MixAcc> mix;

    for (size_t pi = 0; pi < ranked.size(); ++pi) {
        const LigandPose& pose = ranked[pi];
        const double w = (pi < weights.size()) ? weights[pi] : 0.0;

        for (const auto& helix : helices) {
            if (!valid_helix(helix)) continue;
            const HelixContactCount c = count_contacts(helix, rna_nts, pose, cfg);
            if (c.n_hbond == 0 && c.n_stack == 0) continue;
            const std::string label = helix.label.empty() ? std::string("helix") : helix.label;
            HelixThermoPatch patch = make_patch(label, c, cfg, pose.rank, w);
            result.patches.push_back(patch);
            MixAcc& acc = mix[label];
            acc.dH += w * patch.delta_H_kcal;
            acc.dS += w * patch.delta_S_cal_per_K;
            acc.w_sum += w;
            add_touched(acc.touched, patch.touched_nt);
        }

        if (!cfg.require_decision_helix_only) {
            HelixContactCount c;
            std::vector<int> hbond_nts;
            std::vector<int> stack_nts;
            for (const auto& nt : rna_nts) {
                if (nt.nt_index < 0) continue;
                if (nt_in_any_helix(helices, nt.nt_index)) continue;
                if (!finite3(nt.x, nt.y, nt.z)) continue;
                bool hit = false;
                for (const auto& lig : pose.atoms) {
                    if (!finite3(lig.x, lig.y, lig.z)) continue;
                    const double d = dist_A(nt, lig);
                    if (nt.is_hbond_partner && lig.is_hbond_partner && d <= cfg.hbond_max_A) {
                        hit = true;
                        if (std::find(hbond_nts.begin(), hbond_nts.end(), nt.nt_index) ==
                            hbond_nts.end()) {
                            hbond_nts.push_back(nt.nt_index);
                        }
                    }
                    if (nt.is_aromatic && lig.is_aromatic && d <= cfg.stack_max_A &&
                        stack_angle_deg(nt, lig) <= cfg.stack_angle_max_deg) {
                        hit = true;
                        if (std::find(stack_nts.begin(), stack_nts.end(), nt.nt_index) ==
                            stack_nts.end()) {
                            stack_nts.push_back(nt.nt_index);
                        }
                    }
                }
                if (hit) c.touched.push_back(nt.nt_index);
            }
            c.n_hbond = static_cast<int>(hbond_nts.size());
            c.n_stack = static_cast<int>(stack_nts.size());
            if (c.n_hbond > 0 || c.n_stack > 0) {
                HelixThermoPatch patch = make_patch("non_decision", c, cfg, pose.rank, w);
                result.patches.push_back(patch);
                MixAcc& acc = mix["non_decision"];
                acc.dH += w * patch.delta_H_kcal;
                acc.dS += w * patch.delta_S_cal_per_K;
                acc.w_sum += w;
                add_touched(acc.touched, patch.touched_nt);
            }
        }
    }

    for (auto& kv : mix) {
        HelixThermoPatch e;
        e.helix_label = kv.first;
        e.delta_H_kcal = kv.second.dH;
        e.delta_S_cal_per_K = kv.second.dS;
        e.delta_G_kcal = helix_thermo_delta_g_kcal(e.delta_H_kcal, e.delta_S_cal_per_K, cfg.T_K);
        e.touched_nt = kv.second.touched;
        e.pose_rank = -1;
        e.pose_boltzmann_w = kv.second.w_sum;
        result.ensemble.push_back(e);
    }
    std::sort(result.ensemble.begin(), result.ensemble.end(),
              [](const HelixThermoPatch& a, const HelixThermoPatch& b) {
                  return a.helix_label < b.helix_label;
              });
    return result;
}

} // namespace natural
