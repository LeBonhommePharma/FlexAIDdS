// tests/test_rigid_fastpath.cpp
// Bit-identity of FLEXAIDDS_RIGID_FASTPATH (default OFF) vs the legacy
// index_protein + Vcontacts + calc_region path on a tiny synthetic complex.
// Isolated target: Vcontacts.cpp + stubs.cpp + geometry.cpp (no gaboom,
// no vcfunction energy matrix).
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/Vcontacts.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

constexpr const char* kFlag = "FLEXAIDDS_RIGID_FASTPATH";

// VC_Global may later grow scorable_list / n_scorable / scorable_cap /
// fastpath_used. These stay compiling if the header has not landed yet.
template <typename VC_T>
void wire_scorable_fields(VC_T& vc, int* buf, int cap)
{
    if constexpr (requires {
                      vc.scorable_list;
                      vc.n_scorable;
                      vc.scorable_cap;
                      vc.fastpath_used;
                  }) {
        vc.scorable_list = buf;
        vc.n_scorable = 0;
        vc.scorable_cap = cap;
        vc.fastpath_used = 0;
    }
}

template <typename VC_T>
int read_fastpath_used(const VC_T& vc)
{
    if constexpr (requires { vc.fastpath_used; }) {
        return vc.fastpath_used;
    }
    return -1;
}

void clear_fastpath_env()
{
    unsetenv(kFlag);
    unsetenv("FLEXAIDDS_HOIST_RECEPTOR_INDEX");
}

void set_rigid_fastpath(bool on)
{
    if (on) {
        ASSERT_EQ(setenv(kFlag, "1", 1), 0);
    } else {
        unsetenv(kFlag);
    }
}

bool bits_eq(double a, double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

std::string double_bits(double v)
{
    std::uint64_t u = 0;
    std::memcpy(&u, &v, sizeof(u));
    char buf[32];
    std::snprintf(buf, sizeof(buf), "0x%016llx", static_cast<unsigned long long>(u));
    return buf;
}

struct ContactTriple {
    int atom = -1;
    double area = 0.0;
    double dist = 0.0;

    bool operator==(const ContactTriple& o) const
    {
        return atom == o.atom && bits_eq(area, o.area) && bits_eq(dist, o.dist);
    }
};

struct AtomRecord {
    int boxnum = -1;
    double sas_field = 0.0;
    double sas_from_chain = 0.0;
    std::vector<ContactTriple> chain;

    bool operator==(const AtomRecord& o) const
    {
        return boxnum == o.boxnum && bits_eq(sas_field, o.sas_field) &&
               bits_eq(sas_from_chain, o.sas_from_chain) && chain == o.chain;
    }
};

struct PoseRecord {
    int vcontacts_rc = 0;
    int calc_count = 0;
    int dim = 0;
    int fastpath_used = -1;
    bool fastpath_fn = false;
    std::vector<int> boxnum;
    std::vector<int> calclist;
    std::vector<std::vector<int>> box_slices;  // sorted Calclist per box
    std::vector<AtomRecord> scorable;
};

struct IndexRecord {
    int dim = 0;
    int calc_count = 0;
    int fastpath_used = -1;
    bool fastpath_fn = false;
    std::vector<int> boxnum;
    std::vector<int> calclist;
    std::vector<std::vector<int>> box_slices;
};

void expect_pose_eq(const PoseRecord& a, const PoseRecord& b, const char* ctx)
{
    EXPECT_EQ(a.vcontacts_rc, b.vcontacts_rc) << ctx;
    EXPECT_EQ(a.calc_count, b.calc_count) << ctx;
    EXPECT_EQ(a.dim, b.dim) << ctx;
    ASSERT_EQ(a.boxnum.size(), b.boxnum.size()) << ctx;
    for (size_t i = 0; i < a.boxnum.size(); ++i) {
        EXPECT_EQ(a.boxnum[i], b.boxnum[i]) << ctx << " boxnum[" << i << "]";
    }
    ASSERT_EQ(a.calclist.size(), b.calclist.size()) << ctx;
    for (size_t i = 0; i < a.calclist.size(); ++i) {
        EXPECT_EQ(a.calclist[i], b.calclist[i]) << ctx << " Calclist[" << i << "]";
    }
    ASSERT_EQ(a.box_slices.size(), b.box_slices.size()) << ctx;
    for (size_t i = 0; i < a.box_slices.size(); ++i) {
        EXPECT_EQ(a.box_slices[i], b.box_slices[i]) << ctx << " box_slice[" << i << "]";
    }
    ASSERT_EQ(a.scorable.size(), b.scorable.size()) << ctx;
    for (size_t i = 0; i < a.scorable.size(); ++i) {
        EXPECT_EQ(a.scorable[i].boxnum, b.scorable[i].boxnum)
            << ctx << " scorable[" << i << "].boxnum";
        EXPECT_TRUE(bits_eq(a.scorable[i].sas_field, b.scorable[i].sas_field))
            << ctx << " scorable[" << i << "].SAS "
            << double_bits(a.scorable[i].sas_field) << " vs "
            << double_bits(b.scorable[i].sas_field);
        EXPECT_TRUE(bits_eq(a.scorable[i].sas_from_chain, b.scorable[i].sas_from_chain))
            << ctx << " scorable[" << i << "].sas_from_chain "
            << double_bits(a.scorable[i].sas_from_chain) << " vs "
            << double_bits(b.scorable[i].sas_from_chain);
        ASSERT_EQ(a.scorable[i].chain.size(), b.scorable[i].chain.size())
            << ctx << " scorable[" << i << "] ca_rec chain length";
        for (size_t k = 0; k < a.scorable[i].chain.size(); ++k) {
            EXPECT_EQ(a.scorable[i].chain[k].atom, b.scorable[i].chain[k].atom)
                << ctx << " scorable[" << i << "].chain[" << k << "].atom";
            EXPECT_TRUE(bits_eq(a.scorable[i].chain[k].area, b.scorable[i].chain[k].area))
                << ctx << " scorable[" << i << "].chain[" << k << "].area "
                << double_bits(a.scorable[i].chain[k].area) << " vs "
                << double_bits(b.scorable[i].chain[k].area);
            EXPECT_TRUE(bits_eq(a.scorable[i].chain[k].dist, b.scorable[i].chain[k].dist))
                << ctx << " scorable[" << i << "].chain[" << k << "].dist "
                << double_bits(a.scorable[i].chain[k].dist) << " vs "
                << double_bits(b.scorable[i].chain[k].dist);
        }
    }
}

void expect_index_eq(const IndexRecord& a, const IndexRecord& b, const char* ctx)
{
    EXPECT_EQ(a.dim, b.dim) << ctx;
    EXPECT_EQ(a.calc_count, b.calc_count) << ctx;
    ASSERT_EQ(a.boxnum.size(), b.boxnum.size()) << ctx;
    for (size_t i = 0; i < a.boxnum.size(); ++i) {
        EXPECT_EQ(a.boxnum[i], b.boxnum[i]) << ctx << " boxnum[" << i << "]";
    }
    ASSERT_EQ(a.calclist.size(), b.calclist.size()) << ctx;
    for (size_t i = 0; i < a.calclist.size(); ++i) {
        EXPECT_EQ(a.calclist[i], b.calclist[i]) << ctx << " Calclist[" << i << "]";
    }
    ASSERT_EQ(a.box_slices.size(), b.box_slices.size()) << ctx;
    for (size_t i = 0; i < a.box_slices.size(); ++i) {
        EXPECT_EQ(a.box_slices[i], b.box_slices[i]) << ctx << " box_slice[" << i << "]";
    }
}

std::vector<std::vector<int>> sorted_box_slices(const atomindex* box, int dim3,
                                                const int* calclist, int calc_count)
{
    std::vector<std::vector<int>> slices;
    if (!box || dim3 <= 0) return slices;
    slices.resize(static_cast<size_t>(dim3));
    for (int b = 0; b < dim3; ++b) {
        const int n = box[b].nument;
        const int first = box[b].first;
        if (n <= 0 || first < 0) continue;
        std::vector<int> sl;
        sl.reserve(static_cast<size_t>(n));
        for (int k = 0; k < n; ++k) {
            const int idx = first + k;
            if (idx >= 0 && idx < calc_count) sl.push_back(calclist[idx]);
        }
        std::sort(sl.begin(), sl.end());
        slices[static_cast<size_t>(b)] = std::move(sl);
    }
    return slices;
}

std::vector<ContactTriple> walk_ca_chain(const VC_Global& vc, int atom_i)
{
    std::vector<ContactTriple> out;
    if (!vc.ca_index || !vc.ca_rec || atom_i < 0) return out;
    int idx = vc.ca_index[atom_i];
    const int cap = vc.ca_recsize > 0 ? vc.ca_recsize : 0;
    int guard = 0;
    while (idx != -1 && guard++ < cap) {
        if (idx < 0 || idx >= cap) break;
        const ca_struct& c = vc.ca_rec[idx];
        out.push_back({c.atom, c.area, c.dist});
        idx = c.prev;
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// Fixture: 40 rigid receptor atoms on a 5×4×2 lattice + 6 scorable ligand
// atoms (optres != NULL). 0-based atom indices so index_protein's
// `atmi >= atmcnt` skip does not drop the last atom.
// ---------------------------------------------------------------------------
class RigidFastpathBitIdentity : public ::testing::Test {
protected:
    static constexpr int kNRec = 40;
    static constexpr int kNLig = 6;
    static constexpr int kNAtoms = kNRec + kNLig;
    static constexpr float kSpace = 5.0f;
    static constexpr float kRadius = 1.7f;

    // Pose 4 is an out-of-box / clashy translation (grid signature changes).
    static constexpr int kNPoses = 5;
    static constexpr float kPoseDx[kNPoses] = {0.0f, 2.0f, 0.0f, 4.0f, 80.0f};
    static constexpr float kPoseDy[kNPoses] = {0.0f, 0.0f, -1.5f, 1.0f, 0.0f};
    static constexpr float kPoseDz[kNPoses] = {0.0f, 0.0f, 0.8f, -0.5f, 0.0f};

    FA_Global fa_{};
    VC_Global vc_{};
    OptRes lig_opt_{};
    std::vector<atom> atoms_;
    std::vector<resid> residues_;
    std::vector<int> rec_fatm_{0};
    std::vector<int> rec_latm_{kNRec - 1};
    std::vector<int> lig_fatm_{kNRec};
    std::vector<int> lig_latm_{kNAtoms - 1};
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
        clear_fastpath_env();
        build_complex();
        alloc_vc();
    }

    void TearDown() override
    {
        free_vc_box();
        clear_fastpath_env();
    }

    void build_complex()
    {
        atoms_.assign(static_cast<size_t>(kNAtoms), atom{});
        residues_.assign(3, resid{});

        int ai = 0;
        for (int z = 0; z < 2; ++z) {
            for (int y = 0; y < 4; ++y) {
                for (int x = 0; x < 5; ++x, ++ai) {
                    atom& a = atoms_[static_cast<size_t>(ai)];
                    a.coor[0] = static_cast<float>(x) * kSpace;
                    a.coor[1] = static_cast<float>(y) * kSpace;
                    a.coor[2] = static_cast<float>(z) * kSpace;
                    a.radius = kRadius;
                    a.pb_vdw_radius = 1.70;
                    a.number = ai + 1;
                    a.ofres = 1;
                    a.optres = nullptr;
                    a.type = 1;
                    std::strncpy(a.name, "C", 4);
                    std::strncpy(a.element, "C", 2);
                }
            }
        }
        ASSERT_EQ(ai, kNRec);

        // 2×3 ligand grid at z=20 (far enough from the z=5 lattice that
        // voronoi_poly2 finishes) with 4 Å spacing so neighbour contacts exist.
        lig_opt_ = OptRes{};
        lig_opt_.rnum = 2;
        lig_opt_.type = 1;
        lig_opt_.tot = kNLig;

        for (int i = 0; i < kNLig; ++i) {
            atom& a = atoms_[static_cast<size_t>(kNRec + i)];
            a.coor[0] = 8.0f + static_cast<float>(i % 3) * 4.0f;
            a.coor[1] = 6.0f + static_cast<float>(i / 3) * 4.0f;
            a.coor[2] = 20.0f;
            lig_rest_[i][0] = a.coor[0];
            lig_rest_[i][1] = a.coor[1];
            lig_rest_[i][2] = a.coor[2];
            a.radius = kRadius;
            a.pb_vdw_radius = 1.70;
            a.number = kNRec + i + 1;
            a.ofres = 2;
            a.optres = &lig_opt_;
            a.type = 1;
            std::strncpy(a.name, "C", 4);
            std::strncpy(a.element, "C", 2);
        }

        residues_[1].type = 0;
        residues_[1].number = 1;
        residues_[1].rot = 0;
        residues_[1].trot = 1;
        residues_[1].fatm = rec_fatm_.data();
        residues_[1].latm = rec_latm_.data();
        residues_[1].bonded = nullptr;
        std::strncpy(residues_[1].name, "ALA", 3);
        residues_[1].chn = 'A';

        residues_[2].type = 1;
        residues_[2].number = 2;
        residues_[2].rot = 0;
        residues_[2].trot = 1;
        residues_[2].fatm = lig_fatm_.data();
        residues_[2].latm = lig_latm_.data();
        residues_[2].bonded = nullptr;
        std::strncpy(residues_[2].name, "LIG", 3);
        residues_[2].chn = 'L';

        fa_.atm_cnt = kNAtoms;
        fa_.atm_cnt_real = kNAtoms;
        fa_.res_cnt = 2;
        fa_.vindex = 0;
        fa_.num_optres = 1;
        fa_.optres = &lig_opt_;
        fa_.permeability = 1.0f;
        fa_.soft_wall_cutoff = 0.0f;
        fa_.intermolecular_clash_ratio = 0.0f;
        fa_.omit_buried = 0;
        fa_.vcontacts_planedef = 'X';
        // Envelope the in-box poses so only the +80 Å translation expands the grid.
        // Lattice [0,20]×[0,15]×[0,5]; ligand [8,16]×[6,10]×{20} at rest.
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
        calc_.assign(static_cast<size_t>(kNAtoms), atomsas{});
        calclist_.assign(static_cast<size_t>(kNAtoms), -1);
        ca_index_.assign(static_cast<size_t>(kNAtoms), -1);
        seed_.assign(static_cast<size_t>(3 * kNAtoms), -1);
        contlist_.assign(10000, contactlist{});
        ca_rec_.assign(4096, ca_struct{});
        ptorder_.assign(MAX_PT, ptindex{});
        centerpt_.assign(MAX_PT, vertex{});
        poly_.assign(MAX_POLY, vertex{});
        cont_.assign(MAX_PT, plane{});
        vedge_.assign(MAX_POLY, edgevector{});
        scorable_buf_.assign(static_cast<size_t>(kNAtoms), 0);

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
        vc_.recalc = 0;  // do not jitter coordinates on hull failure
        vc_.box = nullptr;
        wire_scorable_fields(vc_, scorable_buf_.data(), kNAtoms);
    }

    void free_vc_box()
    {
        if (vc_.box) {
            std::free(vc_.box);
            vc_.box = nullptr;
        }
    }

    void apply_pose(int pose)
    {
        ASSERT_GE(pose, 0);
        ASSERT_LT(pose, kNPoses);
        for (int i = 0; i < kNLig; ++i) {
            atom& a = atoms_[static_cast<size_t>(kNRec + i)];
            a.coor[0] = lig_rest_[i][0] + kPoseDx[pose];
            a.coor[1] = lig_rest_[i][1] + kPoseDy[pose];
            a.coor[2] = lig_rest_[i][2] + kPoseDz[pose];
        }
    }

    void restore_rest() { apply_pose(0); }

    PoseRecord capture_pose(int rc)
    {
        PoseRecord rec;
        rec.vcontacts_rc = rc;
        rec.calc_count = vc_.calc_count > 0 ? vc_.calc_count : fa_.atm_cnt_real;
        rec.dim = vc_.dim;
        rec.fastpath_used = read_fastpath_used(vc_);
        rec.fastpath_fn = vc_fastpath_active();
        const int n = rec.calc_count;
        rec.boxnum.resize(static_cast<size_t>(n), -1);
        rec.calclist.assign(calclist_.begin(), calclist_.begin() + n);
        for (int i = 0; i < n; ++i) rec.boxnum[static_cast<size_t>(i)] = calc_[static_cast<size_t>(i)].boxnum;

        const int dim3 = rec.dim * rec.dim * rec.dim;
        rec.box_slices = sorted_box_slices(vc_.box, dim3, calclist_.data(), n);

        rec.scorable.reserve(static_cast<size_t>(kNLig));
        for (int i = 0; i < n; ++i) {
            if (!calc_[static_cast<size_t>(i)].score) continue;
            AtomRecord ar;
            ar.boxnum = calc_[static_cast<size_t>(i)].boxnum;
            ar.sas_field = calc_[static_cast<size_t>(i)].SAS;
            ar.chain = walk_ca_chain(vc_, i);
            const float rad = calc_[static_cast<size_t>(i)].atom
                                  ? calc_[static_cast<size_t>(i)].atom->radius
                                  : kRadius;
            const double rado = static_cast<double>(rad) + static_cast<double>(Rw);
            double sas = 4.0 * static_cast<double>(PI) * rado * rado;
            for (const auto& t : ar.chain) sas -= t.area;
            ar.sas_from_chain = sas;
            rec.scorable.push_back(std::move(ar));
        }
        return rec;
    }

    PoseRecord eval_vcontacts(int pose)
    {
        apply_pose(pose);
        wire_scorable_fields(vc_, scorable_buf_.data(), kNAtoms);
        const int rc = Vcontacts(&fa_, atoms_.data(), residues_.data(), &vc_,
                                 /*clash_value=*/nullptr, /*non_scorable=*/false);
        PoseRecord rec = capture_pose(rc);
        free_vc_box();
        return rec;
    }

    std::vector<PoseRecord> eval_series()
    {
        std::vector<PoseRecord> out;
        out.reserve(static_cast<size_t>(kNPoses));
        for (int p = 0; p < kNPoses; ++p) out.push_back(eval_vcontacts(p));
        return out;
    }

    IndexRecord eval_index()
    {
        wire_scorable_fields(vc_, scorable_buf_.data(), kNAtoms);
        std::map<std::string, atomindex*> indexed;
        atomindex* prev = nullptr;
        int dim = 0;
        int calc_count = 0;
        atomindex* box = index_protein(&fa_, atoms_.data(), residues_.data(),
                                       calc_.data(), calclist_.data(), &dim,
                                       fa_.atm_cnt_real, prev, indexed, &calc_count);
        IndexRecord rec;
        rec.dim = dim;
        rec.calc_count = calc_count;
        rec.fastpath_used = read_fastpath_used(vc_);
        rec.fastpath_fn = vc_fastpath_active();
        rec.boxnum.resize(static_cast<size_t>(calc_count), -1);
        rec.calclist.assign(calclist_.begin(), calclist_.begin() + calc_count);
        for (int i = 0; i < calc_count; ++i) {
            rec.boxnum[static_cast<size_t>(i)] = calc_[static_cast<size_t>(i)].boxnum;
        }
        const int dim3 = dim * dim * dim;
        rec.box_slices = sorted_box_slices(box, dim3, calclist_.data(), calc_count);
        if (box) std::free(box);
        return rec;
    }

    void expect_fastpath_off()
    {
        EXPECT_FALSE(vc_fastpath_active())
            << "vc_fastpath_active() must be false when the flag is OFF";
        const int used = read_fastpath_used(vc_);
        if (used >= 0) {
            EXPECT_EQ(used, 0) << "fastpath_used must be 0 when the flag is OFF";
        }
    }
};

// Flag OFF, two evals of the same pose: no sticky state.
TEST_F(RigidFastpathBitIdentity, OffSecondCallMatchesFirstOff)
{
    set_rigid_fastpath(false);
    const PoseRecord a = eval_vcontacts(0);
    const PoseRecord b = eval_vcontacts(0);
    expect_pose_eq(a, b, "OFF call1 vs OFF call2 (pose 0)");
    EXPECT_GE(static_cast<int>(a.scorable.size()), 1);
    expect_fastpath_off();
}

// Gold OFF series vs ON series: boxnum, SAS, and ca_rec chain are bit-identical.
TEST_F(RigidFastpathBitIdentity, OnMatchesOffBitIdentical)
{
    set_rigid_fastpath(false);
    const std::vector<PoseRecord> off = eval_series();
    ASSERT_EQ(static_cast<int>(off.size()), kNPoses);

    restore_rest();
    set_rigid_fastpath(true);
    std::vector<PoseRecord> on;
    on.reserve(static_cast<size_t>(kNPoses));
    for (int p = 0; p < kNPoses; ++p) {
        on.push_back(eval_vcontacts(p));
    }
    ASSERT_EQ(on.size(), off.size());

    for (int p = 0; p < kNPoses; ++p) {
        const std::string ctx = "OFF vs ON pose " + std::to_string(p);
        expect_pose_eq(off[static_cast<size_t>(p)], on[static_cast<size_t>(p)], ctx.c_str());
    }

    // In-box poses should produce a real SAS/contact record, not just empty
    // bookkeeping. OOB (pose 4) may still produce ligand-ligand contacts.
    EXPECT_EQ(off[0].vcontacts_rc, 0) << "pose 0 Vcontacts rc (geometry/stub gap if nonzero)";
    EXPECT_FALSE(off[0].scorable.empty());
    bool any_chain = false;
    for (const auto& ar : off[0].scorable) {
        if (!ar.chain.empty()) any_chain = true;
    }
    EXPECT_TRUE(any_chain) << "expected at least one ca_rec contact on pose 0";
}

// After ON, unsetenv and re-run OFF: still matches the original OFF gold.
TEST_F(RigidFastpathBitIdentity, OffAfterOnStillMatchesOriginalOff)
{
    set_rigid_fastpath(false);
    const std::vector<PoseRecord> gold = eval_series();

    restore_rest();
    set_rigid_fastpath(true);
    (void)eval_series();

    restore_rest();
    set_rigid_fastpath(false);
    const std::vector<PoseRecord> off2 = eval_series();

    ASSERT_EQ(gold.size(), off2.size());
    for (int p = 0; p < kNPoses; ++p) {
        const std::string ctx = "gold OFF vs post-ON OFF pose " + std::to_string(p);
        expect_pose_eq(gold[static_cast<size_t>(p)], off2[static_cast<size_t>(p)], ctx.c_str());
    }
    expect_fastpath_off();
}

// index_protein-only bookkeeping: boxnum + per-box Calclist multiset.
TEST_F(RigidFastpathBitIdentity, IndexProteinOffVsOnBitIdentical)
{
    apply_pose(0);
    set_rigid_fastpath(false);
    const IndexRecord off1 = eval_index();
    const IndexRecord off2 = eval_index();
    expect_index_eq(off1, off2, "index OFF call1 vs OFF call2");
    EXPECT_EQ(off1.calc_count, kNAtoms);
    EXPECT_EQ(static_cast<int>(off1.boxnum.size()), kNAtoms);

    set_rigid_fastpath(true);
    const IndexRecord on1 = eval_index();  // snapshot / first ON
    const IndexRecord on2 = eval_index();  // incremental if the flag is wired
    expect_index_eq(off1, on1, "index OFF vs ON call1");
    expect_index_eq(off1, on2, "index OFF vs ON call2");

    // Move only scorable atoms; rigid boxnums stay put under ON (and should
    // recompute to the same cells under OFF).
    apply_pose(1);
    const IndexRecord on_moved = eval_index();
    ASSERT_EQ(on_moved.boxnum.size(), on2.boxnum.size());
    for (int i = 0; i < kNRec && i < on_moved.calc_count; ++i) {
        EXPECT_EQ(on_moved.boxnum[static_cast<size_t>(i)],
                  on2.boxnum[static_cast<size_t>(i)])
            << "rigid boxnum drifted under ON after ligand move, atom " << i;
    }

    set_rigid_fastpath(false);
    apply_pose(1);
    const IndexRecord off_moved = eval_index();
    expect_index_eq(off_moved, on_moved, "index OFF vs ON after ligand move");

    // OOB translation: still bit-identical, grid may grow.
    apply_pose(4);
    set_rigid_fastpath(false);
    const IndexRecord off_oob = eval_index();
    set_rigid_fastpath(true);
    (void)eval_index();  // first ON at this grid
    const IndexRecord on_oob = eval_index();
    expect_index_eq(off_oob, on_oob, "index OFF vs ON OOB pose");
}

// First ON eval may snapshot (legacy path). Call 2+ on a stable grid must
// report the incremental path once the flag is wired in index_protein.
TEST_F(RigidFastpathBitIdentity, FastpathActiveOnStableSecondOnEval)
{
    set_rigid_fastpath(true);
    const PoseRecord first = eval_vcontacts(0);
    (void)first;  // snapshot; do not require fastpath_used==1
    const PoseRecord second = eval_vcontacts(1);  // same grid, ligand moved
    EXPECT_TRUE(second.fastpath_fn)
        << "vc_fastpath_active() must be true on stable-grid ON call 2+";
    if (second.fastpath_used >= 0) {
        EXPECT_EQ(second.fastpath_used, 1)
            << "fastpath_used must be 1 on stable-grid ON call 2+";
    }

    restore_rest();
    set_rigid_fastpath(false);
    const PoseRecord off = eval_vcontacts(0);
    EXPECT_FALSE(off.fastpath_fn);
    if (off.fastpath_used >= 0) {
        EXPECT_EQ(off.fastpath_used, 0);
    }
}

// Default-OFF: never leave the flag set (TearDown also unsets).
TEST_F(RigidFastpathBitIdentity, FlagDoesNotLeak)
{
    set_rigid_fastpath(true);
    (void)eval_vcontacts(0);
    clear_fastpath_env();
    EXPECT_EQ(std::getenv(kFlag), nullptr);
}
