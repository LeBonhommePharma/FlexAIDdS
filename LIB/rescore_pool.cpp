// =============================================================================
// rescore_pool.cpp — offline pool rescoring driver (see rescore_pool.h)
//
// CF assembly mirrors ic2cf's post-vcfunction block exactly: vcfunction() fills
// per-DoF channel breakdowns into FA->optres[] (which holds ALL optimisable
// degrees of freedom), and the totals are the sum over those entries. On
// intra-ligand clash (error path) the wall penalty is reported with rclash=1,
// mirroring ic2cf's fail-closed behaviour.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
// =============================================================================

#include "rescore_pool.h"
#include "local_orient.h"
#include "flexaid.h"
#include "Vcontacts.h"

#include <algorithm>
#include <cstdlib>
#include <dirent.h>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace {

std::string basename_of(const std::string& path)
{
    const auto pos = path.find_last_of('/');
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

bool dir_exists(const std::string& p)
{
    struct stat st {};
    return ::stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

std::vector<std::string> list_pose_pdbs(const std::string& dir)
{
    std::vector<std::string> out;
    DIR* d = ::opendir(dir.c_str());
    if (!d) return out;
    while (const dirent* e = ::readdir(d)) {
        const std::string n = e->d_name;
        if (n.size() <= 4 || n.substr(n.size() - 4) != ".pdb") continue;
        if (n.find("elected_pose") != std::string::npos) continue;  // §10.1 duplicate
        if (n.find("_INI") != std::string::npos) continue;          // blinded seed
        if (n.find("_rmsd") != std::string::npos) continue;         // member copies
        out.push_back(dir + "/" + n);
    }
    ::closedir(d);
    std::sort(out.begin(), out.end());
    return out;
}

} // namespace

void rescore_pool_mode(FA_Global* FA, VC_Global* VC, atom* atoms,
                       resid* residue, gridpoint* /*cleftgrid*/)
{
    const char* root_env = std::getenv("FLEXAIDDS_RESCORE_POOL");
    if (!root_env || root_env[0] == '\0') return;

    const int lig_res = FA->res_cnt;
    if (lig_res < 1 || FA->atm_cnt <= 0) {
        std::fprintf(stderr, "[RESCORE] SKIP: no prepared atoms\n");
        return;
    }
    const int fa = residue[lig_res].fatm[0];
    const int la = residue[lig_res].latm[0];
    const int n_lig = la - fa + 1;
    if (n_lig <= 0) {
        std::fprintf(stderr, "[RESCORE] SKIP: empty ligand residue\n");
        return;
    }

    // Target label: FLEXAIDDS_RMSDST basename minus _ligand.sdf / .sdf.
    std::string target = "target";
    if (const char* r = std::getenv("FLEXAIDDS_RMSDST"); r && r[0]) {
        target = basename_of(r);
        const std::string suffix_lig = "_ligand.sdf";
        const std::string suffix_sdf = ".sdf";
        if (target.size() > suffix_lig.size() &&
            target.compare(target.size() - suffix_lig.size(),
                           suffix_lig.size(), suffix_lig) == 0)
            target.resize(target.size() - suffix_lig.size());
        else if (target.size() > suffix_sdf.size() &&
                 target.compare(target.size() - suffix_sdf.size(),
                                suffix_sdf.size(), suffix_sdf) == 0)
            target.resize(target.size() - suffix_sdf.size());
    }

    std::string dir = std::string(root_env) + "/" + target;
    if (!dir_exists(dir)) dir = root_env;
    auto files = list_pose_pdbs(dir);
    if (files.empty() && dir != root_env) {
        // Worst-case trap: an empty <root>/<label> shadowing a populated root
        // would silently score nothing. Fall back to the root itself.
        files = list_pose_pdbs(root_env);
        if (!files.empty())
            std::fprintf(stderr,
                "[RESCORE] note: %s held no poses; fell back to pool root\n",
                dir.c_str());
    }
    if (files.empty()) {
        std::fprintf(stderr, "[RESCORE] no pose files for target %s in %s\n",
                     target.c_str(), dir.c_str());
        return;
    }

    // Serial -> slot map. NOTE bounds: the ligand occupies the tail of the
    // atoms[] array and its last index EQUALS FA->atm_cnt (empirically la ==
    // atm_cnt on Astex redocks), so the inclusive upper bound is la.
    const int hi = la;
    std::unordered_map<int,int> serial_to_slot;
    serial_to_slot.reserve(static_cast<size_t>(hi + 1) * 2);
    for (int i = 0; i <= hi; ++i)
        serial_to_slot[atoms[i].number] = i;

    FILE* csv = nullptr;
    if (const char* o = std::getenv("FLEXAIDDS_RESCORE_OUT")) {
        if (o[0] != '\0') csv = std::fopen(o, "w");
        if (csv)
            std::fprintf(csv, "target,file,total,com,wal,sas,con,elec,hbond,"
                              "gist_desolv,metal_coord,entropy,pb_clash,rclash,"
                              "n_intraclashes,n_matched,n_skipped\n");
    }

    // ── WO-LOCALORIENT-1 setup (all inert unless FLEXAIDDS_LOCAL_ORIENT=1) ──
    const auto lo_cfg = flexaids::local_orient::config_from_env();
    const char* lo_dir = std::getenv("FLEXAIDDS_LOCAL_ORIENT_DIR");
    FILE* lo_out = nullptr;
    int lo_ran = 0, lo_improved = 0;
    if (lo_cfg.enabled) {
        if (const char* o = std::getenv("FLEXAIDDS_LOCAL_ORIENT_CSV")) {
            if (o[0] != '\0') lo_out = std::fopen(o, "w");
            if (lo_out)
                std::fprintf(lo_out,
                    "target,file,mode,cf_before,cf_after,cf_delta,evals,"
                    "applied_x,applied_y,applied_z,magnitude,centroid_shift,"
                    "max_bond_drift\n");
        }
        std::fprintf(stderr,
            "[LOCAL-ORIENT] enabled mode=%s budget=%d step0=%.3f min_step=%.3f "
            "out_dir=%s\n",
            (lo_cfg.mode == flexaids::local_orient::Mode::Orient ? "orient"
                                                                 : "jitter"),
            lo_cfg.budget, lo_cfg.step0, lo_cfg.min_step,
            (lo_dir && lo_dir[0]) ? lo_dir : "<none>");
    }

    int done = 0, refused = 0;
    if (std::getenv("FLEXAIDDS_RESCORE_DEBUG")) {
        for (int i = fa; i <= la; ++i)
            std::fprintf(stderr, "[RESCORE] ligslot %d number=%d name=%s\n",
                         i - fa, atoms[i].number, atoms[i].name);
    }
    for (const auto& f : files) {
        const int n_slots = hi + 1;
        std::vector<float> coor(static_cast<size_t>(n_slots) * 3, 0.0f);
        std::vector<unsigned char> mask(static_cast<size_t>(n_slots), 0);
        int matched = 0, skipped = 0;
        flexaids::load_complex_coor_from_pdb(f.c_str(), serial_to_slot,
                                             coor.data(), &matched, &skipped,
                                             mask.data());

        // Fail-closed where it matters: every LIGAND atom slot must be present
        // (vcfunction scores exactly these). Receptor-side stray records the
        // prep layer dropped (united-atom H) are tolerated and counted.
        bool lig_complete = true;
        for (int i = fa; i <= la; ++i) {
            if (!mask[static_cast<size_t>(i)]) { lig_complete = false; break; }
        }
        if (!lig_complete) {
            ++refused;
            std::fprintf(stderr,
                "[RESCORE] REFUSE %s: incomplete ligand coverage "
                "(matched %d records, skipped %d, ligand slots %d)\n",
                basename_of(f).c_str(), matched, skipped, n_lig);
            continue;
        }
        for (int i = 0; i <= hi; ++i) {
            if (mask[static_cast<size_t>(i)]) {
                atoms[i].coor[0] = coor[static_cast<size_t>(i) * 3 + 0];
                atoms[i].coor[1] = coor[static_cast<size_t>(i) * 3 + 1];
                atoms[i].coor[2] = coor[static_cast<size_t>(i) * 3 + 2];
            }
        }
        // ── WO-LOCALORIENT-1: optional local ORIENTATION refinement ─────────
        // Gate default OFF; when unset, lo.ran is false and not one coordinate
        // has moved, so the scoring below is bit-identical to prior behaviour.
        // When ON, atoms[] holds the refined pose and the CF assembled below is
        // the REFINED pose's CF — the receipt records both.
        const auto lo = flexaids::local_orient::refine(lo_cfg, FA, VC, atoms,
                                                       residue, fa, la);
        if (lo.ran) {
            ++lo_ran;
            if (lo.cf_after < lo.cf_before - 1e-9) ++lo_improved;
            std::fprintf(stderr,
                "[LOCAL-ORIENT] file=%s mode=%s cf_before=%.6f cf_after=%.6f "
                "delta=%.6f evals=%d magnitude=%.3f%s centroid_shift=%.6f "
                "max_bond_drift=%.6f\n",
                basename_of(f).c_str(),
                (lo_cfg.mode == flexaids::local_orient::Mode::Orient ? "orient"
                                                                     : "jitter"),
                lo.cf_before, lo.cf_after, lo.cf_after - lo.cf_before,
                lo.evals_used, lo.magnitude,
                (lo_cfg.mode == flexaids::local_orient::Mode::Orient ? "deg" : "A"),
                lo.centroid_shift, lo.max_bond_drift);
            // Fail LOUD, not silently: in orient mode a moved centroid or a
            // changed internal distance means the stage did something other
            // than a rigid rotation, and every downstream number would be
            // attributed to the wrong degree of freedom.
            if (lo_cfg.mode == flexaids::local_orient::Mode::Orient &&
                lo.centroid_shift > 1e-3) {
                std::fprintf(stderr,
                    "[LOCAL-ORIENT] FATAL: centroid moved %.6f A in orient mode "
                    "(must be 0)\n", lo.centroid_shift);
                std::exit(5);
            }
            if (lo.max_bond_drift > 1e-2) {
                std::fprintf(stderr,
                    "[LOCAL-ORIENT] FATAL: internal distance drifted %.6f A "
                    "(transform is not rigid; torsions must be frozen)\n",
                    lo.max_bond_drift);
                std::exit(5);
            }
            if (lo_out) {
                std::fprintf(lo_out,
                    "%s,%s,%s,%.6f,%.6f,%.6f,%d,%.4f,%.4f,%.4f,%.4f,%.6f,%.6f\n",
                    target.c_str(), basename_of(f).c_str(),
                    (lo_cfg.mode == flexaids::local_orient::Mode::Orient
                         ? "orient" : "jitter"),
                    lo.cf_before, lo.cf_after, lo.cf_after - lo.cf_before,
                    lo.evals_used, lo.applied[0], lo.applied[1], lo.applied[2],
                    lo.magnitude, lo.centroid_shift, lo.max_bond_drift);
            }
            // Refined pose is written ALONGSIDE the original, never over it.
            if (lo_dir && lo_dir[0]) {
                const std::string dst =
                    std::string(lo_dir) + "/" + basename_of(f);
                if (!flexaids::local_orient::write_refined_pdb(
                        f.c_str(), dst.c_str(), atoms, fa, la))
                    std::fprintf(stderr,
                        "[LOCAL-ORIENT] WARN: could not write %s\n", dst.c_str());
            }
        }

        std::vector<std::pair<int,int>> intraclashes;
        bool error = false;
        const double penalty =
            vcfunction(FA, VC, atoms, residue, intraclashes, &error);

        cfstr cf{};
        if (error) {
            cf.wal = penalty;
            cf.rclash = 1;
        } else {
            // Exact ic2cf post-vcfunction assembly: optres holds ALL optimisable
            // DoF (side chains AND ligand-associated channels); sum every channel.
            for (int j = 0; j < FA->num_optres; ++j) {
                cf.com         += FA->optres[j].cf.com;
                cf.wal         += FA->optres[j].cf.wal;
                cf.sas         += FA->optres[j].cf.sas;
                cf.con         += FA->optres[j].cf.con;
                cf.elec        += FA->optres[j].cf.elec;
                cf.gist_desolv += FA->optres[j].cf.gist_desolv;
                cf.metal_coord += FA->optres[j].cf.metal_coord;
                cf.hbond       += FA->optres[j].cf.hbond;
                cf.entropy     += FA->optres[j].cf.entropy;
                cf.pb_clash    += FA->optres[j].cf.pb_clash;
            }
        }

        const double total = get_cf_evalue(&cf, FA);
        ++done;
        std::fprintf(stderr,
            "[RESCORE] file=%s cf=%.6f breakdown=com:%.4f,wal:%.4f,sas:%.4f,"
            "con:%.4f,elec:%.4f,hbond:%.4f,gist_desolv:%.4f,metal_coord:%.4f,"
            "entropy:%.4f,pb_clash:%.4f,rclash=%d,intra=%zu\n",
            basename_of(f).c_str(), total,
            static_cast<double>(cf.com), static_cast<double>(cf.wal),
            static_cast<double>(cf.sas), static_cast<double>(cf.con),
            static_cast<double>(cf.elec), static_cast<double>(cf.hbond),
            static_cast<double>(cf.gist_desolv),
            static_cast<double>(cf.metal_coord),
            static_cast<double>(cf.entropy),
            static_cast<double>(cf.pb_clash),
            cf.rclash, intraclashes.size());
        if (csv) {
            std::fprintf(csv, "%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,"
                              "%.6f,%.6f,%.6f,%d,%zu,%d,%d\n",
                target.c_str(), basename_of(f).c_str(), total,
                static_cast<double>(cf.com), static_cast<double>(cf.wal),
                static_cast<double>(cf.sas), static_cast<double>(cf.con),
                static_cast<double>(cf.elec), static_cast<double>(cf.hbond),
                static_cast<double>(cf.gist_desolv),
                static_cast<double>(cf.metal_coord),
                static_cast<double>(cf.entropy),
                static_cast<double>(cf.pb_clash),
                cf.rclash, intraclashes.size(), matched, skipped);
        }
    }

    if (csv) std::fclose(csv);
    if (lo_out) std::fclose(lo_out);
    if (lo_cfg.enabled) {
        std::fprintf(stderr,
            "[LOCAL-ORIENT] done target=%s mode=%s refined=%d improved_cf=%d\n",
            target.c_str(),
            (lo_cfg.mode == flexaids::local_orient::Mode::Orient ? "orient"
                                                                 : "jitter"),
            lo_ran, lo_improved);
    }
    std::fprintf(stderr,
        "[RESCORE] done target=%s scored=%d refused=%d (GA skipped)\n",
        target.c_str(), done, refused);
}
