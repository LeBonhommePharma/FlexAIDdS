// bench_dsvib_cost.cpp — per-term cost of the thermodynamic dS_vib cycle.
//
// WHAT IS TIMED, and why each term is here:
//   cart_lig     build_from_ligand            Cartesian 3N ligand ANM. The term
//                                             ALREADY inside the objective via
//                                             ic2cf.cpp's tencom_weight*h_rep.
//                                             Timed as the reference point the
//                                             ~18 ms/pose figure came from.
//   tors_free    build_from_ligand_torsional  ligand dihedral basis, vacuum.
//   field_gather cell-grid neighbour query    receptor heavy atoms within the
//                                             ENM cutoff of the pose. The only
//                                             receptor-side per-pose work.
//   tors_field   ..._torsional_in_field       ligand dihedral basis, in the
//                                             static receptor field.
//   prot_ca      build_from_ca                the protein torsional model. Timed
//                                             ONCE per target, not per pose: the
//                                             receptor is rigid, so this is a
//                                             run-invariant cost if it were
//                                             needed at all.
//
// Reads the plain-text dumps produced beside this file. No FlexAID runtime, no
// input parsing in the timed region: every coordinate array is materialised
// before the clock starts, so the numbers are the ENM cost and nothing else.
#include "tencm.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace {

struct Case {
    std::string name;
    std::vector<std::array<float,3>> lig;      // heavy atoms
    std::vector<std::string>         lig_el;
    std::vector<std::array<float,3>> rec;      // receptor heavy atoms
    std::vector<std::array<float,3>> ca;       // receptor Calpha
    tencm::TorsionalENM::LigandTopology topo;
};

bool load(const std::string& path, Case& c) {
    std::ifstream f(path);
    if (!f) return false;
    int nl = 0, nr = 0, nc = 0, m = 0;
    f >> nl >> nr >> nc >> m;
    c.lig.resize(nl); c.lig_el.resize(nl);
    for (int i = 0; i < nl; ++i) {
        f >> c.lig[i][0] >> c.lig[i][1] >> c.lig[i][2] >> c.lig_el[i];
    }
    c.rec.resize(nr);
    for (int i = 0; i < nr; ++i) f >> c.rec[i][0] >> c.rec[i][1] >> c.rec[i][2];
    c.ca.resize(nc);
    for (int i = 0; i < nc; ++i) f >> c.ca[i][0] >> c.ca[i][1] >> c.ca[i][2];
    c.topo.rot_bonds.resize(m);
    for (int i = 0; i < m; ++i) f >> c.topo.rot_bonds[i].first >> c.topo.rot_bonds[i].second;
    c.topo.adjacency.assign(nl, {});
    for (int i = 0; i < nl; ++i) {
        int k = 0; f >> k;
        for (int j = 0; j < k; ++j) { int v; f >> v; c.topo.adjacency[i].push_back(v); }
    }
    return f.good() || f.eof();
}

// Same cell grid the objective uses, reimplemented here so the timed gather is
// the same algorithm and not an optimistic stand-in.
struct Grid {
    std::vector<std::array<float,3>> pts;
    std::vector<int> starts, items;
    float lo[3]{0,0,0}; float cell = 9.0f; int nx=0, ny=0, nz=0;
    std::size_t cell_of(const std::array<float,3>& p) const {
        int ix = (int)((p[0]-lo[0])/cell), iy = (int)((p[1]-lo[1])/cell), iz = (int)((p[2]-lo[2])/cell);
        ix = ix<0?0:(ix>=nx?nx-1:ix); iy = iy<0?0:(iy>=ny?ny-1:iy); iz = iz<0?0:(iz>=nz?nz-1:iz);
        return ((std::size_t)iz*ny + iy)*nx + ix;
    }
    void build(const std::vector<std::array<float,3>>& p, float cutoff) {
        pts = p; cell = cutoff;
        if (pts.empty()) return;
        float hi[3];
        for (int d=0;d<3;++d){ lo[d]=pts[0][d]; hi[d]=pts[0][d]; }
        for (const auto& q : pts) for (int d=0;d<3;++d){ if(q[d]<lo[d])lo[d]=q[d]; if(q[d]>hi[d])hi[d]=q[d]; }
        nx = std::max(1,(int)((hi[0]-lo[0])/cell)+1);
        ny = std::max(1,(int)((hi[1]-lo[1])/cell)+1);
        nz = std::max(1,(int)((hi[2]-lo[2])/cell)+1);
        const std::size_t nce = (std::size_t)nx*ny*nz;
        std::vector<int> cnt(nce+1,0);
        for (const auto& q : pts) ++cnt[cell_of(q)+1];
        for (std::size_t c=1;c<=nce;++c) cnt[c]+=cnt[c-1];
        starts = cnt; items.assign(pts.size(),0);
        std::vector<int> fill(nce,0);
        for (std::size_t i=0;i<pts.size();++i){
            const std::size_t c = cell_of(pts[i]);
            items[(std::size_t)starts[c]+fill[c]] = (int)i; ++fill[c];
        }
    }
    void gather(const std::vector<std::array<float,3>>& lig, float rc,
                std::vector<std::array<float,3>>& out, std::vector<char>& picked) const {
        out.clear(); picked.assign(pts.size(),0);
        const float rc2 = rc*rc;
        for (const auto& L : lig) {
            const int ix=(int)((L[0]-lo[0])/cell), iy=(int)((L[1]-lo[1])/cell), iz=(int)((L[2]-lo[2])/cell);
            for (int dz=-1;dz<=1;++dz) for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
                const int jx=ix+dx, jy=iy+dy, jz=iz+dz;
                if (jx<0||jy<0||jz<0||jx>=nx||jy>=ny||jz>=nz) continue;
                const std::size_t c = ((std::size_t)jz*ny + jy)*nx + jx;
                for (int t=starts[c]; t<starts[c+1]; ++t) {
                    const int pi = items[(std::size_t)t];
                    if (picked[(std::size_t)pi]) continue;
                    const auto& P = pts[(std::size_t)pi];
                    const float a=P[0]-L[0], b=P[1]-L[1], d=P[2]-L[2];
                    if (a*a+b*b+d*d > rc2) continue;
                    picked[(std::size_t)pi]=1; out.push_back(P);
                }
            }
        }
    }
};

double ms_since(const std::chrono::steady_clock::time_point& t0) {
    return std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t0).count();
}

// Sum of ln(lambda) over the vibrational subspace = ln det H. The
// calibration-free part of S_vib; printed so the -T*dS check can be done
// without trusting the absolute (uncalibrated) entropy.
double sum_log_eig(const tencm::TorsionalENM& m, int* n_modes) {
    const std::vector<double> ev = m.vibrational_eigenvalues();
    *n_modes = (int)ev.size();
    double s = 0.0;
    for (double e : ev) s += std::log(e);
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    // Last argument is the repetition count; everything before it is a case file.
    if (argc < 3) { std::fprintf(stderr, "usage: %s <case.txt>... <reps>\n", argv[0]); return 2; }
    const int reps  = std::atoi(argv[argc-1]);
    const int n_case = argc - 2;
    const float rc = tencm::DEFAULT_RC;

    std::printf("target,n_lig_heavy,n_rec_heavy,n_ca,M,n_field_nodes,reps,"
                "cart_lig_ms,tors_free_ms,field_gather_ms,tors_field_ms,prot_ca_ms_once,"
                "n_modes_free,n_modes_field,sumlog_free_nats,sumlog_field_nats\n");

    for (int a = 1; a <= n_case; ++a) {
        Case c;
        c.name = argv[a];
        const std::string base = c.name.substr(c.name.find_last_of('/') + 1);
        if (!load(c.name, c)) { std::fprintf(stderr, "load failed: %s\n", c.name.c_str()); continue; }

        // FlexAID atom[] array for the Cartesian reference path (coor + element
        // only, exactly what the post-hoc readers populate).
        std::vector<atom> latoms(c.lig.size());
        for (std::size_t i = 0; i < c.lig.size(); ++i) {
            latoms[i] = atom{};
            latoms[i].coor[0] = c.lig[i][0];
            latoms[i].coor[1] = c.lig[i][1];
            latoms[i].coor[2] = c.lig[i][2];
            std::strncpy(latoms[i].element, c.lig_el[i].c_str(), sizeof(latoms[i].element)-1);
        }

        Grid g; g.build(c.rec, rc);
        std::vector<std::array<float,3>> field; std::vector<char> picked;
        g.gather(c.lig, rc, field, picked);
        const int n_field = (int)field.size();

        // warm-up (page in Eigen workspaces; not timed)
        { tencm::TorsionalENM w; w.build_from_ligand_torsional(c.lig.data(), (int)c.lig.size(), c.topo, rc); }

        double t_cart = 0, t_free = 0, t_gath = 0, t_fld = 0;
        int nm_free = 0, nm_field = 0; double sl_free = 0, sl_field = 0;

        {   auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                tencm::TorsionalENM m;
                m.build_from_ligand(latoms.data(), 0, (int)latoms.size());
            }
            t_cart = ms_since(t0) / reps;
        }
        {   auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                tencm::TorsionalENM m;
                m.build_from_ligand_torsional(c.lig.data(), (int)c.lig.size(), c.topo, rc);
                if (r == reps-1 && m.is_built()) sl_free = sum_log_eig(m, &nm_free);
            }
            t_free = ms_since(t0) / reps;
        }
        {   auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) g.gather(c.lig, rc, field, picked);
            t_gath = ms_since(t0) / reps;
        }
        {   auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                tencm::TorsionalENM m;
                m.build_from_ligand_torsional_in_field(c.lig.data(), (int)c.lig.size(), c.topo,
                                                       field.data(), (int)field.size(), rc);
                if (r == reps-1 && m.is_built()) sl_field = sum_log_eig(m, &nm_field);
            }
            t_fld = ms_since(t0) / reps;
        }
        double t_prot = 0.0;
        if (c.ca.size() >= 3) {
            auto t0 = std::chrono::steady_clock::now();
            tencm::TorsionalENM p; p.build_from_ca(c.ca, rc);
            t_prot = ms_since(t0);
        }

        std::printf("%s,%zu,%zu,%zu,%zu,%d,%d,%.6f,%.6f,%.6f,%.6f,%.3f,%d,%d,%.6f,%.6f\n",
                    base.c_str(), c.lig.size(), c.rec.size(), c.ca.size(),
                    c.topo.rot_bonds.size(), n_field, reps,
                    t_cart, t_free, t_gath, t_fld, t_prot,
                    nm_free, nm_field, sl_free, sl_field);
        std::fflush(stdout);
    }
    return 0;
}