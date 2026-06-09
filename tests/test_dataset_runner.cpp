// =============================================================================
// test_dataset_runner.cpp — GTest suite for DatasetRunner
//
// Tests:
//   - PDB code list validity (correct count for each dataset)
//   - Statistical metric computation (Pearson, Spearman, Kendall on synthetic)
//   - RMSD computation on known coordinates
//   - BenchmarkSet enum parsing
//   - Report generation
//   - PDB HETATM parsing (with inline test data)
//   - Ligand extraction logic
//   - Excluded residue set
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// =============================================================================

// Include gtest FIRST to avoid macro pollution from flexaid.h (e.g. #define E)
#include <gtest/gtest.h>
#define private public
#include "DatasetRunner.h"
#undef private
#include <cmath>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;
using namespace dataset;

// =============================================================================
// PDB code list validity tests
// =============================================================================

TEST(DatasetRunnerCodes, AstexDiverse85Count) {
    auto codes = DatasetRunner::astex_diverse_codes();
    EXPECT_EQ(codes.size(), 85u);
}

TEST(DatasetRunnerCodes, AstexDiverse85NoDuplicates) {
    auto codes = DatasetRunner::astex_diverse_codes();
    std::set<std::string> unique(codes.begin(), codes.end());
    EXPECT_EQ(unique.size(), codes.size()) << "Astex Diverse list has duplicates";
}

TEST(DatasetRunnerCodes, AstexDiverse85ValidFormat) {
    auto codes = DatasetRunner::astex_diverse_codes();
    for (const auto& code : codes) {
        EXPECT_EQ(code.size(), 4u) << "Invalid PDB code length: " << code;
        // PDB codes: digit followed by 3 alphanumeric
        EXPECT_TRUE(std::isdigit(static_cast<unsigned char>(code[0])))
            << "First char should be digit: " << code;
        for (size_t i = 1; i < 4; ++i) {
            EXPECT_TRUE(std::isalnum(static_cast<unsigned char>(code[i])))
                << "Char " << i << " should be alphanumeric: " << code;
        }
    }
}

TEST(DatasetRunnerCodes, AstexDiverse85SpecificCodes) {
    auto codes = DatasetRunner::astex_diverse_codes();
    // Check first and last codes from the hardcoded list
    EXPECT_EQ(codes.front(), "1G9V");
    EXPECT_EQ(codes.back(), "2J62");
    // Check a few known codes are present
    auto has = [&](const std::string& c) {
        return std::find(codes.begin(), codes.end(), c) != codes.end();
    };
    EXPECT_TRUE(has("1Z95"));
    EXPECT_TRUE(has("1UNL"));
    EXPECT_TRUE(has("1HQ2"));
    EXPECT_TRUE(has("2BM2"));
}

TEST(DatasetRunnerCodes, CASF2016Count) {
    auto codes = DatasetRunner::casf2016_codes();
    EXPECT_EQ(codes.size(), 285u);
}

TEST(DatasetRunnerCodes, CASF2016NoDuplicates) {
    auto codes = DatasetRunner::casf2016_codes();
    std::set<std::string> unique(codes.begin(), codes.end());
    EXPECT_EQ(unique.size(), codes.size()) << "CASF-2016 list has duplicates";
}

TEST(DatasetRunnerCodes, DUDETargetCount) {
    auto targets = DatasetRunner::dude_targets();
    EXPECT_EQ(targets.size(), 102u);
}

TEST(DatasetRunnerCodes, DUDENoDuplicates) {
    auto targets = DatasetRunner::dude_targets();
    std::set<std::string> unique(targets.begin(), targets.end());
    EXPECT_EQ(unique.size(), targets.size()) << "DUD-E target list has duplicates";
}

TEST(DatasetRunnerCodes, HAP2Count) {
    auto codes = DatasetRunner::hap2_codes();
    EXPECT_EQ(codes.size(), 59u);
}

TEST(DatasetRunnerCodes, AstexNonNativeTargetCount) {
    auto targets = astex_nonnative_targets();
    // 65 targets (but we have a representative subset)
    EXPECT_GE(targets.size(), 30u);

    // Count total structures
    size_t total = 0;
    for (const auto& t : targets) {
        total += 1 + t.alternative_pdbs.size(); // native + alternatives
    }
    EXPECT_GE(total, 500u) << "Expected at least 500 total structures in Astex Non-Native";
}

// =============================================================================
// Statistical computation tests (synthetic data)
// =============================================================================

TEST(StatisticalMetrics, PearsonPerfectPositive) {
    // Perfect positive correlation
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {2.0, 4.0, 6.0, 8.0, 10.0};
    double r = compute_pearson_r(x, y);
    EXPECT_NEAR(r, 1.0, 1e-10);
}

TEST(StatisticalMetrics, PearsonPerfectNegative) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {10.0, 8.0, 6.0, 4.0, 2.0};
    double r = compute_pearson_r(x, y);
    EXPECT_NEAR(r, -1.0, 1e-10);
}

TEST(StatisticalMetrics, PearsonZero) {
    // Uncorrelated data
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {1.0, -1.0, 1.0, -1.0, 1.0};
    double r = compute_pearson_r(x, y);
    // Not exactly zero but close to zero
    EXPECT_LT(std::abs(r), 0.5);
}

TEST(StatisticalMetrics, PearsonKnownValue) {
    // Known Pearson r ≈ 0.8
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<double> y = {1.2, 2.5, 2.8, 4.1, 5.3, 5.8};
    double r = compute_pearson_r(x, y);
    EXPECT_NEAR(r, 0.995, 0.01);  // Nearly perfect correlation
}

TEST(StatisticalMetrics, SpearmanPerfect) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {10.0, 20.0, 30.0, 40.0, 50.0};
    double rho = compute_spearman_rho(x, y);
    EXPECT_NEAR(rho, 1.0, 1e-10);
}

TEST(StatisticalMetrics, SpearmanPerfectNegative) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {50.0, 40.0, 30.0, 20.0, 10.0};
    double rho = compute_spearman_rho(x, y);
    EXPECT_NEAR(rho, -1.0, 1e-10);
}

TEST(StatisticalMetrics, SpearmanMonotone) {
    // Monotone but non-linear: Spearman should be 1.0
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {1.0, 4.0, 9.0, 16.0, 25.0}; // y = x^2
    double rho = compute_spearman_rho(x, y);
    EXPECT_NEAR(rho, 1.0, 1e-10);
}

TEST(StatisticalMetrics, KendallPerfect) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0};
    double tau = compute_kendall_tau(x, y);
    EXPECT_NEAR(tau, 1.0, 1e-10);
}

TEST(StatisticalMetrics, KendallPerfectNegative) {
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {5.0, 4.0, 3.0, 2.0, 1.0};
    double tau = compute_kendall_tau(x, y);
    EXPECT_NEAR(tau, -1.0, 1e-10);
}

TEST(StatisticalMetrics, KendallWithTies) {
    // With ties: tau-b should handle them
    std::vector<double> x = {1.0, 2.0, 2.0, 3.0};
    std::vector<double> y = {1.0, 2.0, 2.0, 3.0};
    double tau = compute_kendall_tau(x, y);
    EXPECT_NEAR(tau, 1.0, 1e-10);
}

TEST(StatisticalMetrics, KendallKnownValue) {
    // Example from Wikipedia: x = (1,2,3,4,5), y = (3,4,1,2,5)
    // C = 6, D = 4, tau = (6-4)/10 = 0.2
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {3.0, 4.0, 1.0, 2.0, 5.0};
    double tau = compute_kendall_tau(x, y);
    EXPECT_NEAR(tau, 0.2, 1e-10);
}

TEST(StatisticalMetrics, EmptyInput) {
    std::vector<double> empty;
    EXPECT_EQ(compute_pearson_r(empty, empty), 0.0);
    EXPECT_EQ(compute_spearman_rho(empty, empty), 0.0);
    EXPECT_EQ(compute_kendall_tau(empty, empty), 0.0);
}

TEST(StatisticalMetrics, SingleElement) {
    std::vector<double> x = {1.0};
    std::vector<double> y = {2.0};
    EXPECT_EQ(compute_pearson_r(x, y), 0.0);
    EXPECT_EQ(compute_spearman_rho(x, y), 0.0);
    EXPECT_EQ(compute_kendall_tau(x, y), 0.0);
}

// =============================================================================
// RMSD computation tests
// =============================================================================

TEST(RMSDComputation, IdenticalCoords) {
    std::vector<float> coords = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};
    double rmsd = compute_rmsd(coords, coords);
    EXPECT_NEAR(rmsd, 0.0, 1e-10);
}

TEST(RMSDComputation, KnownRMSD) {
    // 2 atoms, one shifted by 1.0 in x, other identical
    // RMSD = sqrt((1^2 + 0 + 0 + 0 + 0 + 0) / 2) = sqrt(0.5) ≈ 0.707
    std::vector<float> a = {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f};
    std::vector<float> b = {1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f};
    double rmsd = compute_rmsd(a, b);
    EXPECT_NEAR(rmsd, std::sqrt(0.5), 1e-5);
}

TEST(RMSDComputation, UniformShift) {
    // All atoms shifted by same amount: RMSD = that amount
    std::vector<float> a = {0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f};
    std::vector<float> b = {2.0f, 0.0f, 0.0f, 3.0f, 1.0f, 1.0f};
    double rmsd = compute_rmsd(a, b);
    EXPECT_NEAR(rmsd, 2.0, 1e-5);
}

TEST(RMSDComputation, EmptyCoords) {
    std::vector<float> empty;
    double rmsd = compute_rmsd(empty, empty);
    EXPECT_GT(rmsd, 100.0); // Should return large value for invalid input
}

TEST(RMSDComputation, MismatchedSize) {
    std::vector<float> a = {0.0f, 0.0f, 0.0f};
    std::vector<float> b = {0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f};
    double rmsd = compute_rmsd(a, b);
    EXPECT_GT(rmsd, 100.0); // Should return large value for mismatched sizes
}

// =============================================================================
// BenchmarkSet enum parsing tests
// =============================================================================

TEST(BenchmarkSetParsing, ValidNames) {
    EXPECT_EQ(parse_benchmark_set("astex"), BenchmarkSet::ASTEX_DIVERSE);
    EXPECT_EQ(parse_benchmark_set("ASTEX"), BenchmarkSet::ASTEX_DIVERSE);
    EXPECT_EQ(parse_benchmark_set("astex_diverse"), BenchmarkSet::ASTEX_DIVERSE);
    EXPECT_EQ(parse_benchmark_set("astex_nonnative"), BenchmarkSet::ASTEX_NON_NATIVE);
    EXPECT_EQ(parse_benchmark_set("hap2"), BenchmarkSet::HAP2);
    EXPECT_EQ(parse_benchmark_set("casf2016"), BenchmarkSet::CASF_2016);
    EXPECT_EQ(parse_benchmark_set("posebusters"), BenchmarkSet::POSEBUSTERS);
    EXPECT_EQ(parse_benchmark_set("dude"), BenchmarkSet::DUD_E);
    EXPECT_EQ(parse_benchmark_set("bindingdb_itc"), BenchmarkSet::BINDINGDB_ITC);
    EXPECT_EQ(parse_benchmark_set("sampl6"), BenchmarkSet::SAMPL6_HG);
    EXPECT_EQ(parse_benchmark_set("sampl7"), BenchmarkSet::SAMPL7_HG);
    EXPECT_EQ(parse_benchmark_set("pdbbind"), BenchmarkSet::PDBBIND_REFINED);
}

TEST(BenchmarkSetParsing, InvalidNames) {
    EXPECT_FALSE(parse_benchmark_set("invalid").has_value());
    EXPECT_FALSE(parse_benchmark_set("").has_value());
    EXPECT_FALSE(parse_benchmark_set("xyz123").has_value());
}

TEST(BenchmarkSetParsing, BenchmarkSetName) {
    EXPECT_EQ(benchmark_set_name(BenchmarkSet::ASTEX_DIVERSE), "Astex Diverse");
    EXPECT_EQ(benchmark_set_name(BenchmarkSet::CASF_2016), "CASF-2016");
    EXPECT_EQ(benchmark_set_name(BenchmarkSet::DUD_E), "DUD-E");
    EXPECT_EQ(benchmark_set_name(BenchmarkSet::SAMPL6_HG), "SAMPL6 Host-Guest");
}

// =============================================================================
// DatasetEntry tests
// =============================================================================

TEST(DatasetEntry, HasAffinityFlags) {
    DatasetEntry entry;
    entry.experimental_affinity = -1.0f;
    EXPECT_FALSE(entry.has_affinity());

    entry.experimental_affinity = 6.5f;
    EXPECT_TRUE(entry.has_affinity());
}

TEST(DatasetEntry, HasEnthalpyFlags) {
    DatasetEntry entry;
    EXPECT_FALSE(entry.has_enthalpy());

    entry.experimental_dH = -7.5f;
    EXPECT_TRUE(entry.has_enthalpy());
}

TEST(DatasetEntry, HasEntropyFlags) {
    DatasetEntry entry;
    EXPECT_FALSE(entry.has_entropy());

    entry.experimental_TdS = -2.3f;
    EXPECT_TRUE(entry.has_entropy());
}

// =============================================================================
// PDB HETATM parsing tests (with inline test data)
// =============================================================================

TEST(PDBParsing, ParseHETATMRecords) {
    // Create a minimal PDB file for testing
    std::string test_dir = "/tmp/flexaidds_test_pdb";
    fs::create_directories(test_dir);
    std::string pdb_path = test_dir + "/test.pdb";

    {
        std::ofstream ofs(pdb_path);
        ofs << "HEADER    TEST PROTEIN\n";
        ofs << "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\n";
        ofs << "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C\n";
        ofs << "HETATM    3  C1  LIG B   1       5.000   6.000   7.000  1.00 20.00           C\n";
        ofs << "HETATM    4  N1  LIG B   1       6.000   7.000   8.000  1.00 20.00           N\n";
        ofs << "HETATM    5  O1  LIG B   1       7.000   8.000   9.000  1.00 20.00           O\n";
        ofs << "HETATM    6  C2  LIG B   1       8.000   9.000  10.000  1.00 20.00           C\n";
        ofs << "HETATM    7  O   HOH C   1      10.000  10.000  10.000  1.00  5.00           O\n";
        ofs << "END\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    auto atoms = runner.parse_pdb_hetatm(pdb_path);

    // Should have 5 HETATM records (4 ligand + 1 water)
    EXPECT_EQ(atoms.size(), 5u);

    // Check first ligand atom
    EXPECT_EQ(atoms[0].resName, "LIG");
    EXPECT_NEAR(atoms[0].x, 5.0f, 0.01f);
    EXPECT_NEAR(atoms[0].y, 6.0f, 0.01f);
    EXPECT_NEAR(atoms[0].z, 7.0f, 0.01f);
    EXPECT_EQ(atoms[0].element, "C");

    // Check water atom
    EXPECT_EQ(atoms[4].resName, "HOH");

    // Cleanup
    fs::remove_all(test_dir);
}

TEST(PDBParsing, ExtractLigandFromPDB) {
    std::string test_dir = "/tmp/flexaidds_test_extract";
    fs::create_directories(test_dir);
    std::string pdb_path = test_dir + "/test.pdb";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(pdb_path);
        ofs << "HEADER    TEST\n";
        ofs << "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\n";
        ofs << "ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C\n";
        // Ligand with 5 atoms
        ofs << "HETATM    3  C1  LIG B   1       5.000   6.000   7.000  1.00 20.00           C\n";
        ofs << "HETATM    4  N1  LIG B   1       6.000   7.000   8.000  1.00 20.00           N\n";
        ofs << "HETATM    5  O1  LIG B   1       7.000   8.000   9.000  1.00 20.00           O\n";
        ofs << "HETATM    6  C2  LIG B   1       8.000   9.000  10.000  1.00 20.00           C\n";
        ofs << "HETATM    7  N2  LIG B   1       9.000  10.000  11.000  1.00 20.00           N\n";
        // Water
        ofs << "HETATM    8  O   HOH C   1      20.000  20.000  20.000  1.00  5.00           O\n";
        // Ion
        ofs << "HETATM    9 ZN   ZN  D   1      25.000  25.000  25.000  1.00  5.00          ZN\n";
        ofs << "END\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    bool extracted = runner.extract_ligand(pdb_path, sdf_path);
    EXPECT_TRUE(extracted);
    EXPECT_TRUE(fs::exists(sdf_path));
    EXPECT_GT(fs::file_size(sdf_path), 0u);

    // Read SDF and verify header
    std::ifstream ifs(sdf_path);
    std::string line;
    std::getline(ifs, line); // molecule name
    EXPECT_EQ(line, "LIG");  // should be the ligand residue name

    // Cleanup
    fs::remove_all(test_dir);
}

// 1TW6 regression: the cognate ligand is the Smac AVPI tetrapeptide stored as
// ATOM records in a short chain; the only HETATM are blacklisted buffers/ions
// (BTB, EDO, ZN, LI).  HETATM-only extraction reports "only blacklisted
// cofactors" and skips the entry.  The peptide fallback must rescue it by
// harvesting the short standard-amino-acid chain and writing a multi-atom SDF.
TEST(PDBParsing, ExtractPeptideLigandFallback) {
    std::string test_dir = "/tmp/flexaidds_test_peptide";
    fs::create_directories(test_dir);
    std::string pdb_path = test_dir + "/test.pdb";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(pdb_path);
        ofs << "HEADER    TEST PEPTIDE LIGAND\n";
        // Receptor chain A: 60 standard residues so a clear receptor exists
        // (kMinRecRes = 50).  Coordinates are far from the peptide so they never
        // bridge into the ligand component.
        int serial = 1;
        for (int r = 1; r <= 60; ++r) {
            float bx = 200.0f + static_cast<float>(r);
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 10.00           C\n",
                serial++, r, bx, 0.0f, 0.0f);
            ofs << buf;
        }
        // Ligand chain C: AVPI tetrapeptide (real 1TW6 chain C coordinates).
        const char* avpi =
            "ATOM   1495  N   ALA C   1      86.680  60.766  15.403  1.00  5.51           N\n"
            "ATOM   1496  CA  ALA C   1      85.268  60.974  15.829  1.00  8.36           C\n"
            "ATOM   1497  C   ALA C   1      84.459  61.499  14.652  1.00  9.27           C\n"
            "ATOM   1498  O   ALA C   1      85.021  61.916  13.640  1.00  8.72           O\n"
            "ATOM   1499  CB  ALA C   1      85.198  61.933  17.006  1.00 10.14           C\n"
            "ATOM   1500  N   VAL C   2      83.140  61.448  14.785  1.00 10.44           N\n"
            "ATOM   1501  CA  VAL C   2      82.223  61.934  13.754  1.00  9.18           C\n"
            "ATOM   1502  C   VAL C   2      81.061  62.685  14.398  1.00  9.13           C\n"
            "ATOM   1503  O   VAL C   2      80.740  62.439  15.558  1.00  8.39           O\n"
            "ATOM   1504  CB  VAL C   2      81.663  60.772  12.895  1.00  9.29           C\n"
            "ATOM   1505  CG1 VAL C   2      82.786  60.057  12.156  1.00 10.46           C\n"
            "ATOM   1506  CG2 VAL C   2      80.853  59.787  13.748  1.00 10.38           C\n"
            "ATOM   1507  N   PRO C   3      80.418  63.597  13.669  1.00  8.40           N\n"
            "ATOM   1508  CA  PRO C   3      79.250  64.295  14.224  1.00  9.15           C\n"
            "ATOM   1509  C   PRO C   3      78.101  63.346  14.595  1.00  9.70           C\n"
            "ATOM   1510  O   PRO C   3      77.916  62.319  13.946  1.00  8.21           O\n"
            "ATOM   1511  CB  PRO C   3      78.840  65.251  13.101  1.00 10.26           C\n"
            "ATOM   1512  CG  PRO C   3      80.065  65.406  12.271  1.00 11.14           C\n"
            "ATOM   1513  CD  PRO C   3      80.741  64.070  12.312  1.00 11.16           C\n"
            "ATOM   1514  N   ILE C   4      77.342  63.698  15.628  1.00 10.20           N\n"
            "ATOM   1515  CA  ILE C   4      76.248  62.854  16.112  1.00 11.44           C\n"
            "ATOM   1516  C   ILE C   4      75.122  62.775  15.082  1.00 10.89           C\n"
            "ATOM   1517  O   ILE C   4      74.908  63.716  14.327  1.00  9.55           O\n"
            "ATOM   1518  CB  ILE C   4      75.716  63.384  17.468  1.00 10.54           C\n"
            "ATOM   1519  CG1 ILE C   4      76.766  63.159  18.562  1.00  9.08           C\n"
            "ATOM   1520  CG2 ILE C   4      74.383  62.711  17.841  1.00 11.57           C\n"
            "ATOM   1521  CD1 ILE C   4      76.443  63.841  19.889  1.00  9.00           C\n";
        ofs << avpi;
        // Blacklisted HETATM only: buffer + metal ions.
        ofs << "HETATM 1566  C1  BTB B 331      69.509  48.419  59.696  1.00 33.71           C\n";
        ofs << "HETATM 9001 ZN    ZN B 401      70.000  50.000  60.000  1.00  5.00          ZN\n";
        ofs << "HETATM 9002  O   HOH B 500      30.000  30.000  30.000  1.00  5.00           O\n";
        ofs << "END\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    bool extracted = runner.extract_ligand(pdb_path, sdf_path);
    EXPECT_TRUE(extracted);
    ASSERT_TRUE(fs::exists(sdf_path));
    EXPECT_GT(fs::file_size(sdf_path), 0u);

    // Title is the N-terminal residue name; counts line must report all 27
    // heavy atoms of the AVPI tetrapeptide.
    std::ifstream ifs(sdf_path);
    std::string l0, l1, l2, counts;
    std::getline(ifs, l0);      // title (N-term residue, e.g. ALA)
    std::getline(ifs, l1);
    std::getline(ifs, l2);
    std::getline(ifs, counts);  // V2000 counts line
    int natoms = std::stoi(counts.substr(0, 3));
    int nbonds = std::stoi(counts.substr(3, 3));
    EXPECT_EQ(natoms, 27);          // AVPI heavy-atom count
    EXPECT_GE(nbonds, 26);          // connected peptide (tree has natoms-1 bonds; PRO ring adds one)

    fs::remove_all(test_dir);
}

TEST(PDBParsing, ExtractLigandFromMMCIFWithChemCompBonds) {
    std::string test_dir = "/tmp/flexaidds_test_extract_mmcif";
    fs::create_directories(test_dir);
    std::string cif_path = test_dir + "/test.cif";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(cif_path);
        ofs << "data_test\n";
        ofs << "#\n";
        ofs << "loop_\n";
        ofs << "_atom_site.group_PDB\n";
        ofs << "_atom_site.id\n";
        ofs << "_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n";
        ofs << "_atom_site.label_alt_id\n";
        ofs << "_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n";
        ofs << "_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n";
        ofs << "_atom_site.Cartn_y\n";
        ofs << "_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n";
        ofs << "_atom_site.B_iso_or_equiv\n";
        ofs << "HETATM 1 C C1 . LIG A 1 0.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 2 C C2 . LIG A 1 5.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 3 C C3 . LIG A 1 10.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 4 C C4 . LIG A 1 15.000 0.000 0.000 1.00 10.00\n";
        ofs << "#\n";
        ofs << "loop_\n";
        ofs << "_chem_comp_bond.comp_id\n";
        ofs << "_chem_comp_bond.atom_id_1\n";
        ofs << "_chem_comp_bond.atom_id_2\n";
        ofs << "_chem_comp_bond.value_order\n";
        ofs << "_chem_comp_bond.pdbx_aromatic_flag\n";
        ofs << "LIG C1 C2 SING N\n";
        ofs << "LIG C2 C3 DOUB N\n";
        ofs << "LIG C3 C4 SING N\n";
        ofs << "#\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    bool extracted = runner.extract_ligand(cif_path, sdf_path);
    EXPECT_TRUE(extracted);
    EXPECT_TRUE(fs::exists(sdf_path));
    EXPECT_GT(fs::file_size(sdf_path), 0u);

    std::ifstream ifs(sdf_path);
    ASSERT_TRUE(ifs.good());
    std::string contents((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_NE(contents.find("FLEXAIDDS_LIGAND_EXTRACTOR_V4"), std::string::npos);

    std::istringstream iss(contents);
    std::string line1, line2, line3, counts;
    std::getline(iss, line1);
    std::getline(iss, line2);
    std::getline(iss, line3);
    std::getline(iss, counts);
    EXPECT_EQ(line1, "LIG");

    int atom_count = 0;
    int bond_count = 0;
    std::istringstream counts_stream(counts);
    counts_stream >> atom_count >> bond_count;
    EXPECT_EQ(atom_count, 4);
    EXPECT_EQ(bond_count, 3);

    fs::remove_all(test_dir);
}

// 1TW6 regression, mmCIF path: the peptide fallback must also fire when the
// structure is an mmCIF (the benchmark prefers the companion CIF for ligand
// identity).  Exercises parse_cif_hetatm_local(..., "ATOM").
TEST(PDBParsing, ExtractPeptideLigandFallbackMMCIF) {
    std::string test_dir = "/tmp/flexaidds_test_peptide_cif";
    fs::create_directories(test_dir);
    std::string cif_path = test_dir + "/test.cif";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(cif_path);
        ofs << "data_test\n#\nloop_\n";
        ofs << "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n";
        // Receptor chain A: 55 residues (> kMinRecRes) far from the peptide.
        int id = 1;
        for (int r = 1; r <= 55; ++r) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "ATOM %d C CA . ALA A %d %.3f 0.000 0.000 1.00 10.00\n",
                id++, r, 200.0 + r);
            ofs << buf;
        }
        // Ligand chain C: AVPI tetrapeptide (real 1TW6 chain C coordinates).
        struct A { const char* el; const char* nm; const char* res; int seq; double x, y, z; };
        const A avpi[] = {
            {"N","N","ALA",1,86.680,60.766,15.403},{"C","CA","ALA",1,85.268,60.974,15.829},
            {"C","C","ALA",1,84.459,61.499,14.652},{"O","O","ALA",1,85.021,61.916,13.640},
            {"C","CB","ALA",1,85.198,61.933,17.006},{"N","N","VAL",2,83.140,61.448,14.785},
            {"C","CA","VAL",2,82.223,61.934,13.754},{"C","C","VAL",2,81.061,62.685,14.398},
            {"O","O","VAL",2,80.740,62.439,15.558},{"C","CB","VAL",2,81.663,60.772,12.895},
            {"C","CG1","VAL",2,82.786,60.057,12.156},{"C","CG2","VAL",2,80.853,59.787,13.748},
            {"N","N","PRO",3,80.418,63.597,13.669},{"C","CA","PRO",3,79.250,64.295,14.224},
            {"C","C","PRO",3,78.101,63.346,14.595},{"O","O","PRO",3,77.916,62.319,13.946},
            {"C","CB","PRO",3,78.840,65.251,13.101},{"C","CG","PRO",3,80.065,65.406,12.271},
            {"C","CD","PRO",3,80.741,64.070,12.312},{"N","N","ILE",4,77.342,63.698,15.628},
            {"C","CA","ILE",4,76.248,62.854,16.112},{"C","C","ILE",4,75.122,62.775,15.082},
            {"O","O","ILE",4,74.908,63.716,14.327},{"C","CB","ILE",4,75.716,63.384,17.468},
            {"C","CG1","ILE",4,76.766,63.159,18.562},{"C","CG2","ILE",4,74.383,62.711,17.841},
            {"C","CD1","ILE",4,76.443,63.841,19.889},
        };
        for (const auto& a : avpi) {
            char buf[200];
            std::snprintf(buf, sizeof(buf),
                "ATOM %d %s %s . %s C %d %.3f %.3f %.3f 1.00 10.00\n",
                id++, a.el, a.nm, a.res, a.seq, a.x, a.y, a.z);
            ofs << buf;
        }
        // Blacklisted HETATM only.
        ofs << "HETATM " << id++ << " C C1 . BTB B 331 69.509 48.419 59.696 1.00 33.71\n";
        ofs << "HETATM " << id++ << " ZN ZN . ZN B 401 70.000 50.000 60.000 1.00 5.00\n";
        ofs << "#\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    bool extracted = runner.extract_ligand(cif_path, sdf_path);
    EXPECT_TRUE(extracted);
    ASSERT_TRUE(fs::exists(sdf_path));

    std::ifstream ifs(sdf_path);
    std::string l0, l1, l2, counts;
    std::getline(ifs, l0);
    std::getline(ifs, l1);
    std::getline(ifs, l2);
    std::getline(ifs, counts);
    int natoms = std::stoi(counts.substr(0, 3));
    EXPECT_EQ(natoms, 27);   // AVPI heavy atoms, harvested from the CIF ATOM block

    fs::remove_all(test_dir);
}

// Glycan deprioritisation: when an N-glycan (NAG/BMA…) is co-crystallised with a
// real small-molecule ligand, the linked sugar tree is the LARGEST connected
// component and would be mis-extracted as the docking target (1GPK: 4×NAG vs
// HUP).  extract_ligand() must drop the sugars and select the drug — but only
// when a non-sugar candidate exists (2GBP→BGC keeps the sugar when it is alone).
TEST(PDBParsing, ExtractLigandDeprioritisesGlycanWhenDrugPresent) {
    std::string test_dir = "/tmp/flexaidds_test_glycan_dedup";
    fs::create_directories(test_dir);
    std::string cif_path = test_dir + "/test.cif";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(cif_path);
        ofs << "data_test\n#\nloop_\n";
        ofs << "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n";
        // N-glycan tree: 8 NAG/BMA atoms, spatially contiguous (distance-bridged
        // into one large component) — larger than the 4-atom drug.
        ofs << "HETATM 1 C C1 . NAG A 1 0.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 2 C C2 . NAG A 1 1.500 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 3 C C3 . NAG A 1 3.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 4 O O3 . NAG A 1 4.500 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 5 C C1 . NAG A 2 6.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 6 C C2 . NAG A 2 7.500 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 7 C C1 . BMA A 3 9.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 8 O O5 . BMA A 3 10.500 0.000 0.000 1.00 10.00\n";
        // Real drug, far away (separate component), only 4 atoms.
        ofs << "HETATM 9 C C1 . DRG B 1 50.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 10 C C2 . DRG B 1 51.500 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 11 N N1 . DRG B 1 53.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 12 O O1 . DRG B 1 54.500 0.000 0.000 1.00 10.00\n";
        ofs << "#\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    ASSERT_TRUE(runner.extract_ligand(cif_path, sdf_path));
    std::ifstream ifs(sdf_path);
    std::string title;
    std::getline(ifs, title);
    EXPECT_EQ(title, "DRG");   // sugar tree dropped, drug selected

    fs::remove_all(test_dir);
}

// Inverse case: a sugar-binding protein where the monosaccharide IS the cognate
// ligand (2GBP→BGC).  With no non-sugar candidate present, the sugar must be
// kept rather than discarded.
TEST(PDBParsing, ExtractLigandKeepsSugarWhenItIsTheOnlyCandidate) {
    std::string test_dir = "/tmp/flexaidds_test_sugar_only";
    fs::create_directories(test_dir);
    std::string cif_path = test_dir + "/test.cif";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(cif_path);
        ofs << "data_test\n#\nloop_\n";
        ofs << "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n_atom_site.label_alt_id\n_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n";
        ofs << "HETATM 1 C C1 . BGC A 1 0.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 2 C C2 . BGC A 1 1.500 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 3 C C3 . BGC A 1 3.000 0.000 0.000 1.00 10.00\n";
        ofs << "HETATM 4 O O3 . BGC A 1 4.500 0.000 0.000 1.00 10.00\n";
        ofs << "#\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    ASSERT_TRUE(runner.extract_ligand(cif_path, sdf_path));
    std::ifstream ifs(sdf_path);
    std::string title;
    std::getline(ifs, title);
    EXPECT_EQ(title, "BGC");   // only candidate → sugar kept

    fs::remove_all(test_dir);
}

TEST(PDBParsing, PrepareEntryRegeneratesStaleLigandCache) {
    std::string test_dir = "/tmp/flexaidds_test_prepare_cache";
    std::string cache_dir = test_dir + "/cache";
    std::string entry_dir = cache_dir + "/Demo Dataset/TEST";
    fs::create_directories(entry_dir);

    std::string pdb_path = entry_dir + "/TEST.pdb";
    std::string sdf_path = entry_dir + "/TEST_ligand.sdf";

    {
        std::ofstream ofs(pdb_path);
        ofs << "HEADER    TEST CACHE ENTRY\n";
        for (int i = 0; i < 40; ++i) {
            ofs << "REMARK    CACHE PADDING LINE " << i << " ............................................................\n";
        }
        ofs << "HETATM    1  C1  LIG A   1       1.000   2.000   3.000  1.00 20.00           C\n";
        ofs << "HETATM    2  C2  LIG A   1       2.000   3.000   4.000  1.00 20.00           C\n";
        ofs << "HETATM    3  O1  LIG A   1       3.000   4.000   5.000  1.00 20.00           O\n";
        ofs << "HETATM    4  N1  LIG A   1       4.000   5.000   6.000  1.00 20.00           N\n";
        ofs << "END\n";
    }
    {
        std::ofstream ofs(sdf_path);
        ofs << "LIG\n";
        ofs << "  stale cache\n";
        ofs << "  missing extractor marker\n";
        ofs << "  0  0  0  0  0  0  0  0  0  0  0  0 V2000\n";
        ofs << "M  END\n";
        ofs << "$$$$\n";
    }

    auto source_time = fs::last_write_time(pdb_path);
    fs::last_write_time(sdf_path, source_time + std::chrono::hours(1));

    DatasetRunner runner(cache_dir);
    auto entry = runner.prepare_pdb_entry("TEST", "Demo Dataset", 0.0f, 0.0f, 0.0f);

    EXPECT_EQ(entry.pdb_id, "TEST");
    EXPECT_TRUE(fs::exists(entry.ligand_path));
    EXPECT_EQ(entry.ligand_path, sdf_path);

    std::ifstream ifs(sdf_path);
    ASSERT_TRUE(ifs.good());
    std::string contents((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_NE(contents.find("FLEXAIDDS_LIGAND_EXTRACTOR_V4"), std::string::npos);

    fs::remove_all(test_dir);
}

TEST(PDBParsing, PrepareEntryUsesCompanionCifForConvertedAtomLigands) {
    std::string test_dir = "/tmp/flexaidds_test_prepare_companion_cif";
    std::string cache_dir = test_dir + "/cache";
    std::string entry_dir = cache_dir + "/Demo Dataset/TEST";
    fs::create_directories(entry_dir);

    std::string pdb_path = entry_dir + "/TEST.pdb";
    std::string cif_path = entry_dir + "/TEST.cif";
    std::string sdf_path = entry_dir + "/TEST_ligand.sdf";

    {
        std::ofstream ofs(pdb_path);
        ofs << "HEADER    CONVERTED ATOM CACHE\n";
        for (int i = 0; i < 4; ++i) {
            ofs << "REMARK    padding padding padding padding padding padding padding\n";
        }
        ofs << "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N\n";
        ofs << "ATOM      2  C1  SAH B   0      10.000   0.000   0.000  1.00 20.00           C\n";
        ofs << "ATOM      3  C2  SAH B   0      11.500   0.000   0.000  1.00 20.00           C\n";
        ofs << "ATOM      4  O1  SAH B   0      13.000   0.000   0.000  1.00 20.00           O\n";
        ofs << "ATOM      5  C1  SKF C   0       1.000   2.000   3.000  1.00 20.00           C\n";
        ofs << "ATOM      6  C2  SKF C   0       2.400   2.000   3.000  1.00 20.00           C\n";
        ofs << "ATOM      7  N1  SKF C   0       3.800   2.000   3.000  1.00 20.00           N\n";
        ofs << "ATOM      8  O1  SKF C   0       5.200   2.000   3.000  1.00 20.00           O\n";
        ofs << "ATOM      9  C1  SKF D   0      21.000   2.000   3.000  1.00 20.00           C\n";
        ofs << "ATOM     10  C2  SKF D   0      22.400   2.000   3.000  1.00 20.00           C\n";
        ofs << "ATOM     11  N1  SKF D   0      23.800   2.000   3.000  1.00 20.00           N\n";
        ofs << "ATOM     12  O1  SKF D   0      25.200   2.000   3.000  1.00 20.00           O\n";
        ofs << "END\n";
    }

    {
        std::ofstream ofs(cif_path);
        ofs << "data_TEST\n";
        ofs << "loop_\n";
        ofs << "_atom_site.group_PDB\n";
        ofs << "_atom_site.id\n";
        ofs << "_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n";
        ofs << "_atom_site.label_alt_id\n";
        ofs << "_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n";
        ofs << "_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n";
        ofs << "_atom_site.Cartn_y\n";
        ofs << "_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n";
        ofs << "_atom_site.B_iso_or_equiv\n";
        ofs << "HETATM 1 C C1 . SAH B 0 10.000 0.000 0.000 1.00 20.00\n";
        ofs << "HETATM 2 C C2 . SAH B 0 11.500 0.000 0.000 1.00 20.00\n";
        ofs << "HETATM 3 O O1 . SAH B 0 13.000 0.000 0.000 1.00 20.00\n";
        ofs << "HETATM 4 C C1 . SKF C 0 1.000 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 5 C C2 . SKF C 0 2.400 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 6 N N1 . SKF C 0 3.800 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 7 O O1 . SKF C 0 5.200 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 8 C C1 . SKF D 0 21.000 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 9 C C2 . SKF D 0 22.400 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 10 N N1 . SKF D 0 23.800 2.000 3.000 1.00 20.00\n";
        ofs << "HETATM 11 O O1 . SKF D 0 25.200 2.000 3.000 1.00 20.00\n";
        ofs << "#\n";
    }

    {
        std::ofstream ofs(sdf_path);
        ofs << "SAH\n";
        ofs << "  stale cache\n";
        ofs << "  Extracted from structure HETATM records | FLEXAIDDS_LIGAND_EXTRACTOR_V4\n";
        ofs << "  3  2  0  0  0  0  0  0  0999 V2000\n";
        ofs << "   10.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n";
        ofs << "   11.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n";
        ofs << "   13.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n";
        ofs << "  1  2  1  0  0  0  0\n";
        ofs << "  2  3  1  0  0  0  0\n";
        ofs << "M  END\n$$$$\n";
    }
    fs::last_write_time(sdf_path, fs::last_write_time(cif_path) + std::chrono::hours(1));

    DatasetRunner runner(cache_dir);
    auto entry = runner.prepare_pdb_entry("TEST", "Demo Dataset", 0.0f, 0.0f, 0.0f);

    ASSERT_EQ(entry.ligand_path, sdf_path);
    std::ifstream ligand(sdf_path);
    ASSERT_TRUE(ligand.good());
    std::string title;
    std::getline(ligand, title);
    EXPECT_EQ(title, "SKF");

    ASSERT_TRUE(fs::exists(entry.receptor_path));
    std::ifstream apo(entry.receptor_path);
    ASSERT_TRUE(apo.good());
    std::string apo_contents((std::istreambuf_iterator<char>(apo)),
                             std::istreambuf_iterator<char>());
    EXPECT_EQ(apo_contents.find("SAH"), std::string::npos);
    EXPECT_EQ(apo_contents.find("SKF"), std::string::npos);
    EXPECT_NE(apo_contents.find("ALA"), std::string::npos);

    fs::remove_all(test_dir);
}

TEST(PDBParsing, ReadPdbSplitsConvertedNonpolymerResidues) {
    std::string test_dir = "/tmp/flexaidds_test_read_pdb_nag";
    fs::create_directories(test_dir);
    std::string pdb_path = test_dir + "/nag.pdb";

    {
        std::ofstream ofs(pdb_path);
        ofs << "ATOM      1  C1  NAG B   0       0.000   0.000   0.000  1.00 20.00           C\n";
        ofs << "ATOM      2  C2  NAG B   0       1.400   0.000   0.000  1.00 20.00           C\n";
        ofs << "ATOM      3  O5  NAG B   0       0.000   1.400   0.000  1.00 20.00           O\n";
        ofs << "ATOM      4  C1  NAG B   0       5.000   0.000   0.000  1.00 20.00           C\n";
        ofs << "ATOM      5  C2  NAG B   0       6.400   0.000   0.000  1.00 20.00           C\n";
        ofs << "ATOM      6  O5  NAG B   0       5.000   1.400   0.000  1.00 20.00           O\n";
        ofs << "END\n";
    }

    FA_Global FA;
    std::memset(static_cast<void*>(&FA), 0, sizeof(FA));
    FA.MIN_NUM_ATOM = 32;
    FA.MIN_NUM_RESIDUE = 16;
    FA.ntypes = 40;
    atom* atoms = nullptr;
    resid* residue = nullptr;

    read_pdb(&FA, &atoms, &residue, const_cast<char*>(pdb_path.c_str()));

    EXPECT_EQ(FA.res_cnt, 2);
    EXPECT_EQ(atoms[1].ofres, 1);
    EXPECT_EQ(atoms[4].ofres, 2);
    EXPECT_EQ(residue[1].type, 1);
    EXPECT_EQ(residue[2].type, 1);

    for (int r = 0; r <= FA.res_cnt; ++r) {
        free(residue[r].fatm);
        free(residue[r].latm);
        free(residue[r].bond);
    }
    free(residue);
    free(atoms);
    free(FA.num_atm);
    fs::remove_all(test_dir);
}

TEST(PDBParsing, DownloadStructureUsesCachedCIFBeforeNetwork) {
    std::string test_dir = "/tmp/flexaidds_test_download_cached_cif";
    std::string entry_dir = test_dir + "/2HKK";
    fs::create_directories(entry_dir);

    std::string cif_path = entry_dir + "/2HKK.cif";
    {
        std::ofstream ofs(cif_path);
        ofs << "data_2HKK\n";
        ofs << "loop_\n";
        ofs << "_atom_site.group_PDB\n";
        ofs << "_atom_site.id\n";
        ofs << "_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n";
        ofs << "_atom_site.label_alt_id\n";
        ofs << "_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n";
        ofs << "_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n";
        ofs << "_atom_site.Cartn_y\n";
        ofs << "_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n";
        ofs << "_atom_site.B_iso_or_equiv\n";
        for (int i = 0; i < 80; ++i) {
            ofs << "HETATM " << (i + 1) << " C C" << (i + 1)
                << " . LIG A 1 " << (0.5 * i) << " 0.000 0.000 1.00 10.00\n";
        }
        ofs << "#\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    std::string receptor_path;
    bool ok = runner.download_structure("2HKK", entry_dir, receptor_path);
    EXPECT_TRUE(ok);
    EXPECT_EQ(receptor_path, cif_path);
    EXPECT_TRUE(fs::exists(receptor_path));

    fs::remove_all(test_dir);
}

TEST(PDBParsing, DownloadStructureRejectsCachedHtmlPdbAndUsesCIF) {
    std::string test_dir = "/tmp/flexaidds_test_download_cached_html";
    std::string entry_dir = test_dir + "/2HKK";
    fs::create_directories(entry_dir);

    std::string pdb_path = entry_dir + "/2HKK.pdb";
    {
        std::ofstream ofs(pdb_path);
        ofs << "<!DOCTYPE html>\n";
        ofs << "<html><body>404</body></html>\n";
    }

    std::string cif_path = entry_dir + "/2HKK.cif";
    {
        std::ofstream ofs(cif_path);
        ofs << "data_2HKK\n";
        ofs << "loop_\n";
        ofs << "_atom_site.group_PDB\n";
        ofs << "_atom_site.id\n";
        ofs << "_atom_site.type_symbol\n";
        ofs << "_atom_site.auth_atom_id\n";
        ofs << "_atom_site.label_alt_id\n";
        ofs << "_atom_site.auth_comp_id\n";
        ofs << "_atom_site.auth_asym_id\n";
        ofs << "_atom_site.auth_seq_id\n";
        ofs << "_atom_site.Cartn_x\n";
        ofs << "_atom_site.Cartn_y\n";
        ofs << "_atom_site.Cartn_z\n";
        ofs << "_atom_site.occupancy\n";
        ofs << "_atom_site.B_iso_or_equiv\n";
        for (int i = 0; i < 80; ++i) {
            ofs << "HETATM " << (i + 1) << " C C" << (i + 1)
                << " . LIG A 1 " << (0.5 * i) << " 0.000 0.000 1.00 10.00\n";
        }
        ofs << "#\n";
    }

    DatasetRunner runner(test_dir + "/cache");
    std::string receptor_path;
    bool ok = runner.download_structure("2HKK", entry_dir, receptor_path);
    EXPECT_TRUE(ok);
    EXPECT_EQ(receptor_path, cif_path);
    EXPECT_FALSE(fs::exists(pdb_path)) << "HTML error page should be purged from cache";

    fs::remove_all(test_dir);
}

// =============================================================================
// Excluded residues tests
// =============================================================================

TEST(ExcludedResidues, WaterExcluded) {
    DatasetRunner runner("/tmp/flexaidds_test_excl/cache");
    // Water should be in the excluded set — verify via extract_ligand behavior
    std::string test_dir = "/tmp/flexaidds_test_excl";
    fs::create_directories(test_dir);
    std::string pdb_path = test_dir + "/water_only.pdb";
    std::string sdf_path = test_dir + "/ligand.sdf";

    {
        std::ofstream ofs(pdb_path);
        ofs << "HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00  5.00           O\n";
        ofs << "HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00  5.00           O\n";
        ofs << "END\n";
    }

    bool extracted = runner.extract_ligand(pdb_path, sdf_path);
    EXPECT_FALSE(extracted) << "Should not extract water as ligand";

    fs::remove_all(test_dir);
}

// =============================================================================
// Report generation tests
// =============================================================================

TEST(ReportGeneration, EmptyReport) {
    BenchmarkReport report;
    report.dataset_name = "Test Dataset";
    report.total_systems = 0;

    std::string test_dir = "/tmp/flexaidds_test_report";
    fs::create_directories(test_dir);

    DatasetRunner runner(test_dir + "/cache");
    runner.write_report(report, test_dir);

    // Should create markdown and CSV files
    EXPECT_TRUE(fs::exists(test_dir + "/test_dataset_report.md"));
    EXPECT_TRUE(fs::exists(test_dir + "/test_dataset_results.csv"));
    EXPECT_TRUE(fs::exists(test_dir + "/test_dataset_summary.csv"));

    fs::remove_all(test_dir);
}

TEST(ReportGeneration, WithResults) {
    BenchmarkReport report;
    report.dataset_name = "Astex Diverse";
    report.total_systems = 3;
    report.successful = 2;
    report.success_rate = 2.0 / 3.0;
    report.mean_rmsd = 1.5;
    report.median_rmsd = 1.2;
    report.pearson_r = 0.85;
    report.spearman_rho = 0.82;
    report.kendall_tau = 0.70;

    DockingResult r1;
    r1.pdb_id = "1ABC";
    r1.best_score = -8.5f;
    r1.rmsd_to_crystal = 0.9f;
    r1.predicted_dG = -8.5f;
    r1.predicted_dH = -6.0f;
    r1.predicted_TdS = -2.5f;
    r1.predicted_IEE = 3.2f;
    r1.num_poses = 15;
    r1.wall_time_s = 12.5;
    r1.success = true;

    DockingResult r2;
    r2.pdb_id = "2DEF";
    r2.best_score = -7.2f;
    r2.rmsd_to_crystal = 1.5f;
    r2.predicted_dG = -7.2f;
    r2.predicted_dH = -5.0f;
    r2.predicted_TdS = -2.2f;
    r2.predicted_IEE = 2.8f;
    r2.num_poses = 10;
    r2.wall_time_s = 15.0;
    r2.success = true;

    DockingResult r3;
    r3.pdb_id = "3GHI";
    r3.best_score = -5.0f;
    r3.rmsd_to_crystal = 3.5f;
    r3.predicted_dG = -5.0f;
    r3.predicted_dH = -3.0f;
    r3.predicted_TdS = -2.0f;
    r3.predicted_IEE = 4.1f;
    r3.num_poses = 5;
    r3.wall_time_s = 20.0;
    r3.success = false;
    report.results = {r1, r2, r3};

    std::string test_dir = "/tmp/flexaidds_test_report2";
    fs::create_directories(test_dir);

    DatasetRunner runner(test_dir + "/cache");
    runner.write_report(report, test_dir);

    // Verify CSV content
    std::ifstream csv(test_dir + "/astex_diverse_results.csv");
    EXPECT_TRUE(csv.good());

    std::string header;
    std::getline(csv, header);
    EXPECT_TRUE(header.find("pdb_id") != std::string::npos);
    EXPECT_TRUE(header.find("rmsd_to_crystal") != std::string::npos);

    // Read first data line
    std::string data_line;
    std::getline(csv, data_line);
    EXPECT_TRUE(data_line.find("1ABC") != std::string::npos);

    // Verify summary CSV
    std::ifstream summary(test_dir + "/astex_diverse_summary.csv");
    EXPECT_TRUE(summary.good());

    std::string summary_header;
    std::getline(summary, summary_header);
    EXPECT_TRUE(summary_header.find("pearson_r") != std::string::npos);

    fs::remove_all(test_dir);
}

// =============================================================================
// DatasetRunner construction and path tests
// =============================================================================

TEST(DatasetRunnerConstruction, DefaultCacheDir) {
    DatasetRunner runner;
    std::string cache = runner.cache_dir();
    EXPECT_FALSE(cache.empty());
    EXPECT_TRUE(cache.find("flexaidds") != std::string::npos ||
                cache.find("benchmarks") != std::string::npos);
}

TEST(DatasetRunnerConstruction, CustomCacheDir) {
    std::string custom_dir = "/tmp/flexaidds_test_custom_cache";
    DatasetRunner runner(custom_dir);
    EXPECT_EQ(runner.cache_dir(), custom_dir);

    // Should have created the directory
    EXPECT_TRUE(fs::exists(custom_dir));

    fs::remove_all(custom_dir);
}

// =============================================================================
// Additional statistical tests for edge cases
// =============================================================================

TEST(StatisticalMetrics, PearsonConstantX) {
    // All x values the same — undefined correlation
    std::vector<double> x = {5.0, 5.0, 5.0, 5.0};
    std::vector<double> y = {1.0, 2.0, 3.0, 4.0};
    double r = compute_pearson_r(x, y);
    EXPECT_NEAR(r, 0.0, 1e-10); // Should return 0 for degenerate case
}

TEST(StatisticalMetrics, PearsonLargeN) {
    // Test with larger dataset
    const int N = 1000;
    std::vector<double> x(N), y(N);
    for (int i = 0; i < N; ++i) {
        x[i] = static_cast<double>(i);
        y[i] = 2.0 * i + 1.0; // perfect linear
    }
    double r = compute_pearson_r(x, y);
    EXPECT_NEAR(r, 1.0, 1e-10);
}

TEST(StatisticalMetrics, SpearmanWithTies) {
    std::vector<double> x = {1.0, 2.0, 2.0, 4.0};
    std::vector<double> y = {1.0, 3.0, 3.0, 4.0};
    double rho = compute_spearman_rho(x, y);
    // With ties, Spearman should still be close to 1.0
    EXPECT_GT(rho, 0.9);
}

TEST(StatisticalMetrics, KendallTwoElements) {
    std::vector<double> x = {1.0, 2.0};
    std::vector<double> y = {1.0, 2.0};
    double tau = compute_kendall_tau(x, y);
    EXPECT_NEAR(tau, 1.0, 1e-10);
}

// =============================================================================
// Prepare from PDB list test
// =============================================================================

TEST(PrepareFromList, ParsePDBList) {
    std::string test_dir = "/tmp/flexaidds_test_pdblist";
    fs::create_directories(test_dir);

    // Create a PDB list file
    std::string list_path = test_dir + "/pdb_list.txt";
    {
        std::ofstream ofs(list_path);
        ofs << "# Comment line\n";
        ofs << "1UNL\n";
        ofs << "1HQ2 6.5\n";  // with affinity
        ofs << "  1Z95  \n";  // with whitespace
        ofs << "\n";           // empty line
        ofs << "2BM2\n";
    }

    DatasetRunner runner(test_dir + "/cache");

    // Don't actually download — just test the parsing logic
    // We can't easily test download without network, but we can verify
    // the function doesn't crash and returns entries
    // (In CI, these PDB downloads may not work)

    // Just verify the list file is readable
    std::ifstream ifs(list_path);
    EXPECT_TRUE(ifs.good());

    int line_count = 0;
    std::string line;
    while (std::getline(ifs, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty() || line[0] == '#') continue;
        line_count++;
    }
    EXPECT_EQ(line_count, 4); // 4 valid PDB entries

    fs::remove_all(test_dir);
}

// =============================================================================
// DockingResult/BenchmarkReport structure tests
// =============================================================================

TEST(DockingResult, DefaultValues) {
    DockingResult r;
    EXPECT_EQ(r.best_score, 0.0f);
    EXPECT_EQ(r.rmsd_to_crystal, 999.0f);
    EXPECT_FALSE(r.success);
    EXPECT_EQ(r.num_poses, 0);
}

TEST(BenchmarkReport, DefaultValues) {
    BenchmarkReport report;
    EXPECT_EQ(report.total_systems, 0);
    EXPECT_EQ(report.successful, 0);
    EXPECT_EQ(report.success_rate, 0.0);
    EXPECT_TRUE(report.results.empty());
}

// =============================================================================
// SAMPL6/7 data integrity tests
// =============================================================================

TEST(SAMPL6Data, ThermodynamicConsistency) {
    // ΔG = ΔH - TΔS, so ΔG ≈ ΔH - TΔS for each entry
    // This is a sanity check on the hardcoded data
    DatasetRunner runner("/tmp/flexaidds_test_sampl6/cache");
    auto entries = runner.prepare(BenchmarkSet::SAMPL6_HG);

    EXPECT_EQ(entries.size(), 27u); // 8 OA + 8 TEMOA + 11 CB8

    for (const auto& entry : entries) {
        if (entry.has_enthalpy() && entry.has_entropy()) {
            // ΔG = ΔH - TΔS
            float dG_from_affinity = -entry.experimental_affinity * 1.3636f;
            float dG_from_components = entry.experimental_dH - entry.experimental_TdS;
            // Allow some rounding tolerance
            EXPECT_NEAR(dG_from_affinity, dG_from_components, 0.5f)
                << "Thermodynamic inconsistency for " << entry.pdb_id
                << ": dG(aff)=" << dG_from_affinity
                << " dG(H-TS)=" << dG_from_components;
        }
    }

    fs::remove_all("/tmp/flexaidds_test_sampl6");
}

TEST(SAMPL7Data, EntryCount) {
    DatasetRunner runner("/tmp/flexaidds_test_sampl7/cache");
    auto entries = runner.prepare(BenchmarkSet::SAMPL7_HG);

    EXPECT_EQ(entries.size(), 30u);

    for (const auto& entry : entries) {
        EXPECT_EQ(entry.source, "SAMPL7-HG");
    }

    fs::remove_all("/tmp/flexaidds_test_sampl7");
}
