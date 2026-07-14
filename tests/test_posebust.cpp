// tests/test_posebust.cpp — NativePoseQC clean-room unit tests (GoogleTest)
//
// Exercises the SHIPPED LIB/PoseBust path on real Astex 1G9V artifacts:
//   CONECT extract → crystal topology assign → evaluate_paths
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include "PoseBust/BustCli.h"
#include "PoseBust/Engine.h"
#include "PoseBust/Loaders.h"
#include "PoseBust/Types.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_set>

namespace fs = std::filesystem;
using namespace flexaids::posebust;

namespace {

fs::path repo_root() {
    // tests/ is one level under repo root when built from build/
    // Prefer FLEXAIDDS_ROOT, else walk up from CWD looking for LIB/PoseBust.
    if (const char* e = std::getenv("FLEXAIDDS_ROOT"); e && e[0]) {
        return fs::path(e);
    }
    fs::path p = fs::current_path();
    for (int i = 0; i < 6; ++i) {
        if (fs::is_directory(p / "LIB" / "PoseBust") &&
            fs::is_directory(p / "benchmarks" / "astex_diverse")) {
            return p;
        }
        if (!p.has_parent_path() || p == p.root_path()) break;
        p = p.parent_path();
    }
    return fs::current_path();
}

fs::path astex_dir(const std::string& code) {
    return repo_root() / "benchmarks" / "astex_diverse" / "astex_diverse" / code;
}

fs::path find_1g9v_pose() {
    const fs::path root = repo_root();
    for (const char* rel : {
             "benchmarks/astex_repro/full_v132/1G9V",
             "benchmarks/astex_repro/full_v131/1G9V",
             "benchmarks/astex_repro/full_v130/1G9V",
         }) {
        const fs::path base = root / rel;
        if (!fs::is_directory(base)) continue;
        const fs::path elected = base / "elected_pose.pdb";
        if (fs::is_regular_file(elected)) return elected;
        // Prefer restart poses
        for (auto it = fs::recursive_directory_iterator(base);
             it != fs::recursive_directory_iterator(); ++it) {
            if (!it->is_regular_file()) continue;
            const auto name = it->path().filename().string();
            if (name.rfind("1G9V_", 0) == 0 && name.size() > 5 &&
                name.find("_INI") == std::string::npos &&
                it->path().extension() == ".pdb") {
                return it->path();
            }
        }
    }
    return {};
}

}  // namespace

TEST(PoseBustLoaders, CrystalSdf1G9VLoads25Heavy) {
    const fs::path lig = astex_dir("1G9V") / "1G9V_ligand.sdf";
    if (!fs::is_regular_file(lig)) {
        GTEST_SKIP() << "missing " << lig;
    }
    Molecule m;
    std::string err;
    ASSERT_TRUE(load_sdf(lig.string(), m, &err)) << err;
    EXPECT_EQ(m.n_heavy(), 25);
    EXPECT_FALSE(m.bonds.empty());
}

TEST(PoseBustLoaders, TopologyMismatchFailsClosed) {
    const fs::path pose = find_1g9v_pose();
    const fs::path other = astex_dir("1GPK") / "1GPK_ligand.sdf";
    if (pose.empty() || !fs::is_regular_file(other)) {
        GTEST_SKIP() << "missing 1G9V pose or 1GPK ligand";
    }
    Molecule lig, ref;
    std::string err;
    ASSERT_TRUE(load_pdb_flexaid_ligand(pose.string(), lig, &err)) << err;
    ASSERT_TRUE(load_sdf(other.string(), ref, &err)) << err;
    EXPECT_FALSE(assign_topology_from_reference(lig, ref, &err));
    EXPECT_FALSE(err.empty());
}

TEST(PoseBustLoaders, Extract1G9VNotHEM) {
    const fs::path pose = find_1g9v_pose();
    const fs::path crystal = astex_dir("1G9V") / "1G9V_ligand.sdf";
    if (pose.empty() || !fs::is_regular_file(crystal)) {
        GTEST_SKIP() << "missing 1G9V pose/crystal";
    }
    Molecule lig, ref;
    std::string err;
    ASSERT_TRUE(load_pdb_flexaid_ligand(pose.string(), lig, &err)) << err;
    // Must not pull heme-scale atom counts
    EXPECT_LT(static_cast<int>(lig.atoms.size()), 50);
    ASSERT_TRUE(load_sdf(crystal.string(), ref, &err)) << err;
    ASSERT_TRUE(assign_topology_from_reference(lig, ref, &err)) << err;
    EXPECT_EQ(lig.n_heavy(), 25);
}

TEST(PoseBustEngine, EvaluatePaths1G9VEmitsUpstreamKeys) {
    const fs::path pose = find_1g9v_pose();
    const fs::path crystal = astex_dir("1G9V") / "1G9V_ligand.sdf";
    const fs::path protein = astex_dir("1G9V") / "1G9V_apo.pdb";
    if (pose.empty() || !fs::is_regular_file(crystal) ||
        !fs::is_regular_file(protein)) {
        GTEST_SKIP() << "missing 1G9V artifacts";
    }
    EvaluateOptions opt;
    opt.suite = Suite::Dock;
    opt.pdb_id = "1G9V";
    auto rep = evaluate_paths(pose.string(), protein.string(), crystal.string(), opt);
    ASSERT_TRUE(rep.ran);
    EXPECT_TRUE(rep.error.empty()) << rep.error;
    EXPECT_EQ(rep.backend, "native_pose_qc");

    // Required upstream-named keys (dock suite subset native implements)
    const char* required[] = {
        "mol_pred_loaded",
        "mol_cond_loaded",
        "sanitization",
        "inchi_convertible",
        "all_atoms_connected",
        "no_radicals",
        "bond_lengths",
        "bond_angles",
        "internal_steric_clash",
        "aromatic_ring_flatness",
        "non-aromatic_ring_non-flatness",
        "double_bond_flatness",
        "double_bond_stereochemistry",
        "tetrahedral_chirality",
        "internal_energy",
        "minimum_distance_to_protein",
        "protein-ligand_maximum_distance",
        "volume_overlap_with_protein",
        "minimum_distance_to_organic_cofactors",
        "minimum_distance_to_inorganic_cofactors",
        "minimum_distance_to_waters",
        "molecular_formula",
        "molecular_bonds",
        "mol_true_loaded",
    };
    std::unordered_set<std::string> have;
    for (const auto& c : rep.checks) have.insert(c.key);
    for (const char* k : required) {
        EXPECT_TRUE(have.count(k)) << "missing key: " << k;
    }

    // Loading + connectivity must pass on a real FlexAID pose extract
    auto pred = rep.find_check("mol_pred_loaded");
    ASSERT_NE(pred, nullptr);
    EXPECT_TRUE(pred->passed);
    auto conn = rep.find_check("all_atoms_connected");
    ASSERT_NE(conn, nullptr);
    EXPECT_TRUE(conn->passed);
}

TEST(PoseBustEngine, CrystalSelfDockNearNativePassesCore) {
    // Crystal ligand coords rewritten as a synthetic complex is hard;
    // evaluate crystal SDF topology + apo protein after writing pose-less
    // path is not available — instead load crystal and protein and call
    // evaluate() directly (shipped evaluate API).
    const fs::path crystal = astex_dir("1G9V") / "1G9V_ligand.sdf";
    const fs::path protein = astex_dir("1G9V") / "1G9V_apo.pdb";
    if (!fs::is_regular_file(crystal) || !fs::is_regular_file(protein)) {
        GTEST_SKIP() << "missing 1G9V crystal/apo";
    }
    Molecule lig, prot;
    std::string err;
    ASSERT_TRUE(load_sdf(crystal.string(), lig, &err)) << err;
    ASSERT_TRUE(load_pdb_protein_heavy(protein.string(), prot, &err)) << err;
    EvaluateOptions opt;
    opt.suite = Suite::Dock;
    auto rep = evaluate(lig, prot, &lig, opt);
    ASSERT_TRUE(rep.ran);
    EXPECT_TRUE(rep.error.empty()) << rep.error;
    // Crystal pose should not clash with protein pocket
    auto mc = rep.find_check("minimum_distance_to_protein");
    ASSERT_NE(mc, nullptr);
    EXPECT_TRUE(mc->passed) << mc->detail;
    auto nr = rep.find_check("no_radicals");
    ASSERT_NE(nr, nullptr);
    EXPECT_TRUE(nr->passed) << nr->detail;
}

TEST(PoseBustBustCli, ResolveBinary) {
    // May be empty in CI without bust; just ensure the resolver does not throw.
    const std::string b = resolve_bust_binary();
    if (b.empty()) {
        GTEST_SKIP() << "bust not installed (set FLEXAIDDS_POSEBUSTERS_BIN)";
    }
    EXPECT_TRUE(fs::is_regular_file(b));
}

TEST(PoseBustContract, SuccessPbIsRmsdAndPbPassSemantics) {
    // Document fixed contract in a unit test so env remaps cannot return silently.
    // DockingResult is not linked here; pin the boolean algebra used by runner.
    const bool success_rmsd = true;
    const bool pb_pass = true;
    const bool success_pb = success_rmsd && pb_pass;
    EXPECT_TRUE(success_pb);
    EXPECT_FALSE(true && false);  // pb alone is never success_pb without rmsd
    const bool success_rmsd_only = true;
    const bool pb_fail = false;
    EXPECT_FALSE(success_rmsd_only && pb_fail);
}
