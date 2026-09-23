// tests/test_vcontacts_last_atom.cpp
// Regression test: index_protein() must index EVERY atom of the selected
// rotamers, using the engine's real 1-based atoms[] layout.
//
// atoms[] is 1-based (atoms[0] unused; SdfReader/read_pdb increment atm_cnt
// before storing), the ligand is the last residue, so its final atom is
// atoms[atm_cnt_real]; build_rotamers() appends rotamer copies above
// atm_cnt_real. Calc[] is 0-based with atm_cnt_real slots. index_protein()
// used to bound the atoms[] index by the Calc slot count (atmcnt ==
// atm_cnt_real, exclusive), which dropped the ligand's last atom -- and every
// atom of a residue sitting in a non-zero rotamer -- from every Vcontacts
// call: no contact, clash or SAS term for that atom.
//
// Isolated link, same as test_rigid_fastpath: Vcontacts.cpp + stubs.cpp +
// geometry.cpp. (test_rigid_fastpath uses 0-based indices precisely so the
// old bound "does not drop the last atom"; this test uses the engine layout.)
// Apache-2.0 (c) 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/Vcontacts.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

constexpr const char* kFastpathFlag = "FLEXAIDDS_RIGID_FASTPATH";

void clear_env()
{
    unsetenv(kFastpathFlag);
    unsetenv("FLEXAIDDS_HOIST_RECEPTOR_INDEX");
}

struct Contact {
    int atom = -1;
    double area = 0.0;
    double dist = 0.0;
    bool operator==(const Contact& o) const
    {
        return atom == o.atom && std::memcmp(&area, &o.area, sizeof area) == 0 &&
               std::memcmp(&dist, &o.dist, sizeof dist) == 0;
    }
};

// 40 rigid receptor atoms on a 5x4x2 lattice (5 A spacing) + a 2x3 ligand
// grid (4 A spacing) at z = 20 -- the geometry test_rigid_fastpath uses, but
// stored 1-based: receptor atoms[1..40], ligand atoms[41..46] == atm_cnt_real.
class VcontactsLastAtom : public ::testing::Test {
protected:
    static constexpr int kNRec = 40;
    static constexpr int kNLig = 6;
    static constexpr int kNReal = kNRec + kNLig;  // atm_cnt_real, last ligand atom index
    static constexpr float kSpace = 5.0f;
    static constexpr float kRadius = 1.7f;

    FA_Global fa_{};
    VC_Global vc_{};
    OptRes lig_opt_{};
    OptRes rec_opt_{};
    std::vector<atom> atoms_;
    std::vector<resid> residues_;
    std::vector<int> rec_fatm_;
    std::vector<int> rec_latm_;
    std::vector<int> lig_fatm_;
    std::vector<int> lig_latm_;
    std::vector<atomsas> calc_;
    std::vector<int> calclist_;
    std::vector<int> ca_index_;
    std::vector<int> seed_;
    std::vector<contactlist> contlist_;
    std::vector<ca_struct> ca_rec_;
    std::vector<ptindex> ptorder_;
    std::vector<vertex> centerpt_;
    std::vector<vertex> poly_;
    std::vector<plane> cont_;
    std::vector<edgevector> vedge_;
    std::vector<int> scorable_buf_;
    float lig_rest_[kNLig][3]{};

    void SetUp() override
    {
        clear_env();
        build_complex(/*with_rotamer=*/false);
        alloc_vc();
    }

    void TearDown() override
    {
        free_box();
        clear_env();
    }

    static void set_atom(atom& a, float x, float y, float z, int number, int ofres,
                         OptRes* optres)
    {
        a.coor[0] = x;
        a.coor[1] = y;
        a.coor[2] = z;
        a.radius = kRadius;
        a.pb_vdw_radius = 1.70;
        a.number = number;
        a.ofres = ofres;
        a.optres = optres;
        a.type = 1;
        std::strncpy(a.name, "C", 4);
        std::strncpy(a.element, "C", 2);
    }

    // with_rotamer: residue 1 also carries a rotamer copy (rot 1) stored at
    // atoms[kNReal+1 .. kNReal+kNRec], i.e. ABOVE atm_cnt_real, exactly where
    // build_rotamers() puts it (fatm[trot] = atm_cnt + 1). The copy is shifted
    // by +0.5 A in x so the two rotamers are distinguishable.
    void build_complex(bool with_rotamer)
    {
        const int n_rot_atoms = with_rotamer ? kNRec : 0;
        atoms_.assign(static_cast<size_t>(1 + kNReal + n_rot_atoms), atom{});
        residues_.assign(3, resid{});

        int ai = 1;  // atoms[0] is unused, as in the engine
        for (int z = 0; z < 2; ++z)
            for (int y = 0; y < 4; ++y)
                for (int x = 0; x < 5; ++x, ++ai)
                    set_atom(atoms_[static_cast<size_t>(ai)], x * kSpace, y * kSpace,
                             z * kSpace, ai, 1, with_rotamer ? &rec_opt_ : nullptr);
        ASSERT_EQ(ai, kNRec + 1);

        lig_opt_ = OptRes{};
        lig_opt_.rnum = 2;
        lig_opt_.type = 1;
        lig_opt_.tot = kNLig;
        for (int i = 0; i < kNLig; ++i) {
            const float x = 8.0f + static_cast<float>(i % 3) * 4.0f;
            const float y = 6.0f + static_cast<float>(i / 3) * 4.0f;
            lig_rest_[i][0] = x;
            lig_rest_[i][1] = y;
            lig_rest_[i][2] = 20.0f;
            set_atom(atoms_[static_cast<size_t>(kNRec + 1 + i)], x, y, 20.0f,
                     90001 + i, 2, &lig_opt_);
        }

        rec_fatm_ = {1};
        rec_latm_ = {kNRec};
        if (with_rotamer) {
            rec_opt_ = OptRes{};
            rec_opt_.rnum = 1;
            rec_opt_.type = 0;
            rec_opt_.tot = kNRec;
            for (int k = 0; k < kNRec; ++k) {
                atom& src = atoms_[static_cast<size_t>(1 + k)];
                atom& dst = atoms_[static_cast<size_t>(kNReal + 1 + k)];
                dst = src;
                dst.coor[0] += 0.5f;
            }
            rec_fatm_.push_back(kNReal + 1);
            rec_latm_.push_back(kNReal + kNRec);
        }
        lig_fatm_ = {kNRec + 1};
        lig_latm_ = {kNReal};

        residues_[1].type = 0;
        residues_[1].number = 1;
        residues_[1].rot = with_rotamer ? 1 : 0;
        residues_[1].trot = with_rotamer ? 1 : 0;
        residues_[1].fatm = rec_fatm_.data();
        residues_[1].latm = rec_latm_.data();
        residues_[1].bonded = nullptr;
        std::strncpy(residues_[1].name, "ALA", 3);
        residues_[1].chn = 'A';

        residues_[2].type = 1;
        residues_[2].number = 2;
        residues_[2].rot = 0;
        residues_[2].trot = 0;
        residues_[2].fatm = lig_fatm_.data();
        residues_[2].latm = lig_latm_.data();
        residues_[2].bonded = nullptr;
        std::strncpy(residues_[2].name, "LIG", 3);
        residues_[2].chn = 'L';

        fa_ = FA_Global{};
        fa_.atm_cnt = kNReal + n_rot_atoms;  // index of the last live atom
        fa_.atm_cnt_real = kNReal;           // Calc[] slot count
        fa_.res_cnt = 2;
        fa_.vindex = 0;
        fa_.num_optres = 1;
        fa_.optres = &lig_opt_;
        fa_.permeability = 1.0f;
        fa_.soft_wall_cutoff = 0.0f;
        fa_.intermolecular_clash_ratio = 0.0f;
        fa_.omit_buried = 0;
        fa_.vcontacts_planedef = 'X';
        fa_.globalmin[0] = -6.0f;
        fa_.globalmin[1] = -6.0f;
        fa_.globalmin[2] = -6.0f;
        fa_.globalmax[0] = 28.0f;
        fa_.globalmax[1] = 22.0f;
        fa_.globalmax[2] = 28.0f;
        fa_.maxwidth = 36.0f;
    }

    void alloc_vc()
    {
        const size_t n = static_cast<size_t>(kNReal);  // top.cpp: malloc(atm_cnt_real)
        calc_.assign(n, atomsas{});
        calclist_.assign(n, -1);
        ca_index_.assign(n, -1);
        seed_.assign(3 * n, -1);
        contlist_.assign(10000, contactlist{});
        ca_rec_.assign(4096, ca_struct{});
        ptorder_.assign(MAX_PT, ptindex{});
        centerpt_.assign(MAX_PT, vertex{});
        poly_.assign(MAX_POLY, vertex{});
        cont_.assign(MAX_PT, plane{});
        vedge_.assign(MAX_POLY, edgevector{});
        scorable_buf_.assign(n, 0);

        vc_ = VC_Global{};
        vc_.Calc = calc_.data();
        vc_.Calclist = calclist_.data();
        vc_.ca_index = ca_index_.data();
        vc_.seed = seed_.data();
        vc_.contlist = contlist_.data();
        vc_.ca_rec = ca_rec_.data();
        vc_.ca_recsize = static_cast<int>(ca_rec_.size());
        vc_.numcarec = 0;
        vc_.ptorder = ptorder_.data();
        vc_.centerpt = centerpt_.data();
        vc_.poly = poly_.data();
        vc_.cont = cont_.data();
        vc_.vedge = vedge_.data();
        vc_.planedef = 'X';
        vc_.recalc = 0;
        vc_.box = nullptr;
        vc_.calc_count = 0;
        vc_.scorable_list = scorable_buf_.data();
        vc_.n_scorable = 0;
        vc_.scorable_cap = kNReal;
        vc_.fastpath_used = 0;
    }

    void free_box()
    {
        if (vc_.box && !fa_.vindex) std::free(vc_.box);
        vc_.box = nullptr;
    }

    void move_ligand(float dx)
    {
        for (int i = 0; i < kNLig; ++i)
            atoms_[static_cast<size_t>(kNRec + 1 + i)].coor[0] = lig_rest_[i][0] + dx;
    }

    // Direct index_protein call, exactly as Vcontacts() makes it.
    int run_index_protein()
    {
        std::map<std::string, atomindex*> indexed;
        int dim = 0;
        int calc_count = -1;
        atomindex* box = index_protein(&fa_, atoms_.data(), residues_.data(), calc_.data(),
                                       calclist_.data(), &dim, fa_.atm_cnt_real, nullptr,
                                       indexed, &calc_count);
        if (box) std::free(box);
        return calc_count;
    }

    // 1-based atoms[] index behind each of the first n Calc slots.
    std::vector<int> indexed_atoms(int n) const
    {
        std::vector<int> out;
        for (int i = 0; i < n; ++i) {
            const atom* a = calc_[static_cast<size_t>(i)].atom;
            out.push_back(a ? static_cast<int>(a - atoms_.data()) : -1);
        }
        std::sort(out.begin(), out.end());
        return out;
    }

    int calc_slot_of(int atom_index, int n) const
    {
        for (int i = 0; i < n; ++i)
            if (calc_[static_cast<size_t>(i)].atom == &atoms_[static_cast<size_t>(atom_index)])
                return i;
        return -1;
    }

    std::vector<Contact> chain_of(int slot) const
    {
        std::vector<Contact> out;
        int idx = vc_.ca_index[slot];
        int guard = 0;
        while (idx != -1 && guard++ < vc_.ca_recsize) {
            const ca_struct& c = vc_.ca_rec[idx];
            out.push_back({c.atom, c.area, c.dist});
            idx = c.prev;
        }
        return out;
    }

    int run_vcontacts()
    {
        const int rc = Vcontacts(&fa_, atoms_.data(), residues_.data(), &vc_,
                                 /*clash_value=*/nullptr, /*non_scorable=*/false);
        return rc;
    }

    static std::vector<int> range(int lo, int hi)
    {
        std::vector<int> v;
        for (int k = lo; k <= hi; ++k) v.push_back(k);
        return v;
    }
};

// The ligand's last atom, atoms[atm_cnt_real], must get a Calc slot.
TEST_F(VcontactsLastAtom, IndexProteinIndexesLastOneBasedAtom)
{
    const int n = run_index_protein();
    ASSERT_EQ(n, kNReal) << "index_protein must fill all atm_cnt_real Calc slots";
    EXPECT_EQ(indexed_atoms(n), range(1, kNReal))
        << "Calc must hold atoms[1..atm_cnt_real]: atoms[0] is unused and "
           "atoms[atm_cnt_real] is the ligand's last atom";
    EXPECT_GE(calc_slot_of(kNReal, n), 0) << "atoms[atm_cnt_real] never entered Calc";
}

// Through Vcontacts(): the last ligand atom is scorable, has its own contact
// chain, and appears in its ligand neighbours' chains (4 A away).
TEST_F(VcontactsLastAtom, VcontactsGivesLastLigandAtomContacts)
{
    ASSERT_EQ(run_vcontacts(), 0);
    free_box();
    ASSERT_EQ(vc_.calc_count, kNReal);
    const int last = calc_slot_of(kNReal, vc_.calc_count);
    ASSERT_GE(last, 0) << "the ligand's last atom never entered Vcontacts";
    EXPECT_TRUE(calc_[static_cast<size_t>(last)].score);

    const std::vector<Contact> own = chain_of(last);
    EXPECT_FALSE(own.empty()) << "last ligand atom has no contact record";
    double area = 0.0;
    for (const Contact& c : own) area += c.area;
    EXPECT_GT(area, 0.0);

    bool listed_by_neighbour = false;
    for (int i = 0; i < vc_.calc_count; ++i) {
        if (i == last || !calc_[static_cast<size_t>(i)].score) continue;
        for (const Contact& c : chain_of(i))
            if (c.atom == last) listed_by_neighbour = true;
    }
    EXPECT_TRUE(listed_by_neighbour)
        << "no other ligand atom records a contact with the last ligand atom";
}

// A residue in a non-zero rotamer lives above atm_cnt_real (build_rotamers
// appends copies at atm_cnt+1); its atoms must be indexed, not the rot-0 ones.
TEST_F(VcontactsLastAtom, IndexProteinIndexesRotamerCopiesAboveAtmCntReal)
{
    build_complex(/*with_rotamer=*/true);
    alloc_vc();
    ASSERT_EQ(residues_[1].fatm[1], kNReal + 1);
    ASSERT_GT(fa_.atm_cnt, fa_.atm_cnt_real);

    const int n = run_index_protein();
    ASSERT_EQ(n, kNReal);
    std::vector<int> want = range(kNRec + 1, kNReal);             // ligand
    const std::vector<int> rot = range(kNReal + 1, kNReal + kNRec);  // rotamer 1 copy
    want.insert(want.end(), rot.begin(), rot.end());
    std::sort(want.begin(), want.end());
    EXPECT_EQ(indexed_atoms(n), want)
        << "Calc must hold the selected rotamer's atoms (atoms[" << kNReal + 1 << ".."
        << kNReal + kNRec << "]) plus the whole ligand";
}

// FLEXAIDDS_RIGID_FASTPATH must still engage once the last atom is indexed
// (its rebind guard used the same Calc-count bound) and must rebuild the
// same per-box index as the full path.
TEST_F(VcontactsLastAtom, RigidFastpathEngagesWithOneBasedLastAtom)
{
    auto box_slices = [&]() {
        std::vector<std::vector<int>> s;
        const int dim3 = vc_.dim * vc_.dim * vc_.dim;
        for (int b = 0; b < dim3; ++b) {
            std::vector<int> v;
            for (int k = 0; k < vc_.box[b].nument; ++k)
                v.push_back(vc_.Calclist[vc_.box[b].first + k]);
            std::sort(v.begin(), v.end());
            s.push_back(v);
        }
        return s;
    };

    // Full path at the moved pose.
    move_ligand(2.0f);
    ASSERT_EQ(run_vcontacts(), 0);
    ASSERT_EQ(vc_.calc_count, kNReal);
    const auto full = box_slices();
    const std::vector<Contact> full_last = chain_of(calc_slot_of(kNReal, kNReal));
    free_box();

    // Fresh workspace, fastpath ON: call 1 snapshots, call 2 (ligand moved,
    // same grid) must take the incremental path.
    alloc_vc();
    ASSERT_EQ(setenv(kFastpathFlag, "1", 1), 0);
    move_ligand(0.0f);
    ASSERT_EQ(run_vcontacts(), 0);
    free_box();
    move_ligand(2.0f);
    ASSERT_EQ(run_vcontacts(), 0);
    EXPECT_EQ(vc_.fastpath_used, 1)
        << "rigid fastpath refused to rebind the last ligand atom";
    ASSERT_EQ(vc_.calc_count, kNReal);
    EXPECT_EQ(box_slices(), full) << "fastpath per-box index differs from the full path";
    const int last = calc_slot_of(kNReal, kNReal);
    ASSERT_GE(last, 0);
    EXPECT_EQ(chain_of(last), full_last) << "last-atom contacts differ between paths";
    free_box();
}

}  // namespace
