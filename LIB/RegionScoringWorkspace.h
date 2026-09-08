#pragma once

#include "AtomCopyExtent.h"
#include "AtomOptResBinding.h"
#include "GAContext.h"
#include "Vcontacts.h"
#include "gaboom.h"
#include "ga_constants.h"
#include <algorithm>
#include <cstdlib>
#include <new>
#include <vector>

namespace flexaids {

// Scoring storage for one region. This does not establish complete ownership
// of every nested FA/atom/residue object: region GAs still run serially.
struct RegionScoringWorkspace {
    FA_Global fa;
    GB_Global gb;
    VC_Global vc;
    std::vector<atom> atoms_copy;
    std::vector<resid> residue_copy;
    std::vector<OptRes> optres;
    std::vector<int> contacts, calclist, ca_index, seed, scorable;
    std::vector<float> contributions;
    std::vector<atomsas> calc;
    std::vector<contactlist> contlist;
    std::vector<ptindex> ptorder;
    std::vector<vertex> centerpt, poly;
    std::vector<plane> cont;
    std::vector<edgevector> vedge;
    GAContext ga_ctx;

    RegionScoringWorkspace(const FA_Global& source_fa, const GB_Global& source_gb,
                           const VC_Global& source_vc, const atom* atoms,
                           const resid* residues)
        : fa(source_fa), gb(source_gb), vc(source_vc),
          atoms_copy(atoms, atoms + flexaid_one_based_copy_n(source_fa.atm_cnt)),
          residue_copy(residues, residues + source_fa.res_cnt + 1),
          contacts(CONTACTS_BUFFER_SIZE),
          calclist(source_fa.atm_cnt_real), ca_index(source_fa.atm_cnt_real, -1),
          seed(3 * source_fa.atm_cnt_real), scorable(source_fa.atm_cnt_real),
          contributions(source_fa.ntypes * source_fa.ntypes), calc(source_fa.atm_cnt_real),
          contlist(GA_CONTLIST_SIZE), ptorder(MAX_PT), centerpt(MAX_PT),
          poly(MAX_POLY), cont(MAX_PT), vedge(MAX_POLY) {
        // Never let an allocation failure transfer ownership of parent ca_rec.
        // Allocate this last, after every operation that could otherwise throw.
        vc.ca_rec = nullptr;
        if (source_fa.num_optres > 0)
            optres.assign(source_fa.optres, source_fa.optres + source_fa.num_optres);
        const AtomOptResBinding binding(
            std::span<const atom>(atoms + 1, source_fa.atm_cnt),
            std::span<const OptRes>(source_fa.optres, source_fa.num_optres));
        binding.bind(std::span<atom>(atoms_copy.data() + 1, source_fa.atm_cnt), optres);
        fa.contacts = contacts.data();
        fa.contributions = contributions.data();
        fa.optres = optres.data();
        vc.Calc = calc.data();
        vc.Calclist = calclist.data();
        vc.ca_index = ca_index.data();
        vc.seed = seed.data();
        vc.scorable_list = scorable.data();
        vc.scorable_cap = source_fa.atm_cnt_real;
        vc.n_scorable = 0;
        vc.fastpath_used = 0;
        vc.calc_count = 0;
        vc.numcarec = 0;
        vc.contlist = contlist.data();
        vc.ptorder = ptorder.data();
        vc.centerpt = centerpt.data();
        vc.poly = poly.data();
        vc.cont = cont.data();
        vc.vedge = vedge.data();
        // Vcontacts manages its thread-local indexed box cache. vindex==0
        // allocates/frees a box per scoring call; neither case is ours to free.
        vc.box = nullptr;
        // save_areas() grows ca_rec with realloc: vector storage is invalid here.
        vc.ca_recsize = std::max(source_vc.ca_recsize, 100);
        vc.ca_rec = static_cast<ca_struct*>(
            std::calloc(static_cast<std::size_t>(vc.ca_recsize), sizeof(ca_struct)));
        if (!vc.ca_rec) throw std::bad_alloc();
    }
    ~RegionScoringWorkspace() { std::free(vc.ca_rec); }
    RegionScoringWorkspace(const RegionScoringWorkspace&) = delete;
    RegionScoringWorkspace& operator=(const RegionScoringWorkspace&) = delete;
};

}  // namespace flexaids
