// pose_remarks.cpp — shared pose REMARK writer
//
// Implements LIB/pose_remarks.h. Gate FLEXAIDDS_UNIFIED_REMARKS, default OFF,
// read with flexaids::env_bool (empty string is OFF).
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "pose_remarks.h"
#include "EnvFlags.h"

#include <cstdarg>
#include <cstdio>

namespace flexaidds {
namespace remarks {

void PoseRemarkBuilder::line(const char* fmt, ...)
{
    if (!fmt) return;
    char tmp[2048];
    va_list ap;
    va_start(ap, fmt);
    const int n = vsnprintf(tmp, sizeof(tmp), fmt, ap);
    va_end(ap);
    if (n <= 0) return;
    if (static_cast<size_t>(n) < sizeof(tmp)) {
        buf_ += tmp;
        return;
    }
    // Rare: a single REMARK longer than tmp. Grow once rather than truncate.
    std::string grow(static_cast<size_t>(n) + 1, '\0');
    va_start(ap, fmt);
    vsnprintf(grow.data(), grow.size(), fmt, ap);
    va_end(ap);
    grow.resize(static_cast<size_t>(n));
    buf_ += grow;
}

void PoseRemarkBuilder::append_header()
{
    raw("REMARK optimized structure\n");
}

void PoseRemarkBuilder::append_cf_totals(double cf, double cf_app)
{
    line("REMARK CF=%8.5f\n", cf);
    line("REMARK CF.app=%8.5f\n", cf_app);
}

void PoseRemarkBuilder::append_cf_terms(const FA_Global* FA)
{
    append_cf_terms(FA, nullptr);
}

void PoseRemarkBuilder::append_cf_terms(const FA_Global* FA, const resid* residue)
{
    if (!FA || !FA->optres) return;
    for (int i = 0; i < FA->num_optres; ++i) {
        const OptRes& opt = FA->optres[i];
        if (residue) {
            const resid* res_ptr = &residue[opt.rnum];
            line("REMARK optimizable residue %s %c %d\n",
                 res_ptr->name, res_ptr->chn, res_ptr->number);
        }
        const cfstr* cf_ptr = &opt.cf;
        line("REMARK CF.com=%8.5f\n", cf_ptr->com);
        line("REMARK CF.sas=%8.5f\n", cf_ptr->sas);
        line("REMARK CF.wal=%8.5f\n", cf_ptr->wal);
        line("REMARK CF.con=%8.5f\n", cf_ptr->con);
        line("REMARK CF.gist=%8.5f\n", cf_ptr->gist);
        line("REMARK CF.hbond=%8.5f\n", cf_ptr->hbond);
    }
}

void PoseRemarkBuilder::append_residue_sas(const FA_Global* FA)
{
    if (!FA || !FA->optres) return;
    for (int i = 0; i < FA->num_optres; ++i) {
        line("REMARK Residue has an overall SAS of %.3f\n", FA->optres[i].cf.totsas);
    }
}

void PoseRemarkBuilder::append_torsions(const FA_Global* FA)
{
    if (!FA || !FA->opt_par) return;
    for (int i = 0; i < FA->npar; ++i) {
        line("REMARK [%8.3f]\n", FA->opt_par[i]);
    }
}

void PoseRemarkBuilder::emit_rmsd_pair(FA_Global* FA, atom* atoms, resid* residue,
                                       gridpoint* cleftgrid)
{
    // Local flag. No caller-supplied Hungarian, so a rank-N pose cannot inherit
    // rank-(N-1)'s flag. calc_rmsd is not pure: two calls, raw then sym, in
    // immediate succession, flag flipped between them. Do not memoize or reorder.
    if (!FA || FA->refstructure != 1) return;
    bool Hungarian = false;
    const double rmsd_raw = calc_rmsd(FA, atoms, residue, cleftgrid,
                                      FA->npar, FA->opt_par, Hungarian);
    line("REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n", rmsd_raw);
    line("REMARK rmsd_raw = %.5f\n", rmsd_raw);
    Hungarian = true;
    const double rmsd_sym = calc_rmsd(FA, atoms, residue, cleftgrid,
                                      FA->npar, FA->opt_par, Hungarian);
    line("REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n", rmsd_sym);
    line("REMARK rmsd_sym = %.5f\n", rmsd_sym);
}

void PoseRemarkBuilder::append_dsvib(const FA_Global* FA)
{
    // FA_Global does not carry the ic2cf aggregate cfstr. Emit the six-field
    // receipt with an explicit reason rather than omitting the block.
    (void)FA;
    append_dsvib(static_cast<const cfstr*>(nullptr));
}

void PoseRemarkBuilder::append_dsvib(const cfstr* scored)
{
    const int status = scored ? scored->dsvib_status : 0;
    const int ndof = scored ? scored->dsvib_ndof : 0;
    const int nflex = scored ? scored->dsvib_nflex : 0;
    line("REMARK DSVIB.status=%d basis=torsional "
         "n_torsion_dofs=%d n_flexbonds=%d "
         "units_S=kcal/mol/K units_G=kcal/mol\n",
         status, ndof, nflex);
    line("REMARK DSVIB.S_complex_minus_apo=%.8f\n",
         scored ? scored->dsvib_s_field : 0.0);
    line("REMARK DSVIB.S_ligand_free=%.8f\n",
         scored ? scored->dsvib_s_free : 0.0);
    line("REMARK DSVIB.S_apo=0.00000000 apo_treatment=frozen_receptor_dofs_cancel_exactly\n");
    line("REMARK DSVIB.dS_vib=%.8f\n", scored ? scored->dsvib_ds : 0.0);
    line("REMARK DSVIB.minus_T_dS_vib=%.8f weight=1.0_by_construction\n",
         scored ? scored->minus_T_dsvib : 0.0);
}

void PoseRemarkBuilder::append_inputs(const char* config_path, const char* ga_path)
{
    if (!config_path && !ga_path) return;
    line("REMARK inputs: %s & %s\n",
         config_path ? config_path : "",
         ga_path ? ga_path : "");
}

void PoseRemarkBuilder::append_interaction_contributions(const FA_Global* FA, atom* atoms)
{
    (void)FA;
    (void)atoms;
    // write_pdb.cpp:231 already dumps the +/- contribution block AFTER the
    // remark[] string. Composing it here would duplicate those 20 lines on
    // every pose. Method exists because the header lists it; it is a no-op
    // on purpose so the inheritance point at write_pdb stays the owner.
}

std::string compose_pose_remarks(Emitter e,
                                 FA_Global* FA,
                                 atom* atoms,
                                 resid* residue,
                                 gridpoint* cleftgrid,
                                 const char* config_path,
                                 const char* ga_path)
{
    PoseRemarkBuilder b(e);
    b.append_header();
    b.append_cf_terms(FA, residue);
    b.append_residue_sas(FA);
    b.append_torsions(FA);
    b.emit_rmsd_pair(FA, atoms, residue, cleftgrid);
    b.append_dsvib(FA);
    b.append_emitter_line();
    b.append_inputs(config_path, ga_path);
    b.append_interaction_contributions(FA, atoms);
    return b.str();
}

bool unified_remarks_enabled()
{
    // Default OFF. flexaids::env_bool: unset / empty / unparseable → fallback.
    return flexaids::env_bool("FLEXAIDDS_UNIFIED_REMARKS", false);
}

}  // namespace remarks
}  // namespace flexaidds
