// tencom_ledger.h — ledger-only tENCoM eigenvalues (λ), inert on election
//
// Wave 4: write the real tENCoM spectrum into the thermodynamic REMARK / ledger
// behind FLEXAIDDS_LEDGER_TENCOM_LAMBDA (default OFF). λ must never enter CF
// ranking, clustering, BindingMode::compute_energy(), or output order.
//
// atom::eigen is an eigenvector grid, not λ — this header does not read it.
// Eigenvalues come only from tencm::TorsionalENM (the existing tENCoM channel).
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
#pragma once

#include "EnvFlags.h"
#include "tENCoM/tencm.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace flexaids {

/// Always true. A later change that feeds λ into fitness must fail tests that
/// assert this tag on every ledger record and scan ranking sources.
inline constexpr bool kTencomLambdaInertOnElection = true;

struct TencomLambdaLedger {
    bool flag_enabled = false;
    bool inert_on_election = kTencomLambdaInertOnElection;
    std::string status = "flag_off";  // flag_off | ok | unavailable
    std::vector<double> lambdas;      // real TorsionalENM eigenvalues; never fabricated
    int n = 0;
    double lambda_min = 0.0;
    double lambda_max = 0.0;
};

inline bool ledger_tencom_lambda_enabled() noexcept {
    return env_bool("FLEXAIDDS_LEDGER_TENCOM_LAMBDA", false);
}

inline void fill_lambda_stats(TencomLambdaLedger& rec) {
    rec.n = static_cast<int>(rec.lambdas.size());
    if (rec.lambdas.empty()) {
        rec.lambda_min = 0.0;
        rec.lambda_max = 0.0;
        return;
    }
    rec.lambda_min = rec.lambdas.front();
    rec.lambda_max = rec.lambdas.front();
    for (double lam : rec.lambdas) {
        rec.lambda_min = std::min(rec.lambda_min, lam);
        rec.lambda_max = std::max(rec.lambda_max, lam);
    }
}

/// Copy the real tENCoM spectrum. Does not invent eigenvalues.
inline TencomLambdaLedger spectrum_from_built_model(const tencm::TorsionalENM& model) {
    TencomLambdaLedger rec;
    rec.flag_enabled = ledger_tencom_lambda_enabled();
    rec.inert_on_election = kTencomLambdaInertOnElection;
    if (!model.is_built() || model.modes().empty()) {
        rec.status = "unavailable";
        return rec;
    }
    rec.lambdas.reserve(model.modes().size());
    for (const auto& mode : model.modes()) {
        if (!std::isfinite(mode.eigenvalue)) continue;
        rec.lambdas.push_back(mode.eigenvalue);
    }
    if (rec.lambdas.empty()) {
        rec.status = "unavailable";
        return rec;
    }
    fill_lambda_stats(rec);
    rec.status = "ok";
    return rec;
}

inline TencomLambdaLedger collect_tencom_lambda_from_ca(
    const std::vector<std::array<float, 3>>& ca_coords)
{
    TencomLambdaLedger rec;
    rec.inert_on_election = kTencomLambdaInertOnElection;
    rec.flag_enabled = ledger_tencom_lambda_enabled();
    if (!rec.flag_enabled) {
        rec.status = "flag_off";
        return rec;
    }
    if (ca_coords.size() < 3) {
        rec.status = "unavailable";
        return rec;
    }
    tencm::TorsionalENM model;
    model.build_from_ca(ca_coords);
    rec = spectrum_from_built_model(model);
    rec.flag_enabled = true;
    rec.inert_on_election = kTencomLambdaInertOnElection;
    return rec;
}

inline std::vector<std::array<float, 3>> extract_protein_ca_coords(
    const atom* atoms, const resid* residue, int res_cnt)
{
    std::vector<std::array<float, 3>> ca;
    if (atoms == nullptr || residue == nullptr || res_cnt < 1) return ca;
    for (int ri = 1; ri <= res_cnt; ++ri) {
        if (residue[ri].type != 0) continue;
        if (residue[ri].fatm == nullptr || residue[ri].latm == nullptr) continue;
        const int first = residue[ri].fatm[0];
        const int last = residue[ri].latm[0];
        for (int ai = first; ai <= last; ++ai) {
            const char* nm = atoms[ai].name;
            const bool is_ca =
                (nm[0] == ' ' && nm[1] == 'C' && nm[2] == 'A' && nm[3] == ' ') ||
                (nm[0] == 'C' && nm[1] == 'A' && nm[2] == ' ' && nm[3] == ' ');
            if (!is_ca) continue;
            ca.push_back({atoms[ai].coor[0], atoms[ai].coor[1], atoms[ai].coor[2]});
            break;
        }
    }
    return ca;
}

inline TencomLambdaLedger collect_tencom_lambda_from_atoms(
    const atom* atoms, const resid* residue, int res_cnt)
{
    TencomLambdaLedger rec;
    rec.inert_on_election = kTencomLambdaInertOnElection;
    rec.flag_enabled = ledger_tencom_lambda_enabled();
    if (!rec.flag_enabled) {
        rec.status = "flag_off";
        return rec;
    }
    return collect_tencom_lambda_from_ca(extract_protein_ca_coords(atoms, residue, res_cnt));
}

/// Empty when the flag is off (no REMARK). Otherwise one compact ledger line.
inline std::string format_tencom_lambda_remark(const TencomLambdaLedger& rec) {
    if (!rec.flag_enabled || rec.status == "flag_off") return {};
    char buf[320];
    const int n = std::snprintf(
        buf, sizeof(buf),
        "REMARK tencom_lambda_ledger inert_on_election=%d status=%s n=%d "
        "lambda_min=%.6g lambda_max=%.6g\n",
        rec.inert_on_election ? 1 : 0,
        rec.status.c_str(),
        rec.n,
        rec.lambda_min,
        rec.lambda_max);
    if (n <= 0) return {};
    return std::string(buf);
}

}  // namespace flexaids
