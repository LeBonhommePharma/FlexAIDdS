// tests/test_ion_handling.cpp
// Unit tests for HETATM / ion handling in assign_radii and assign_types.
// Verifies that metal ions and organic cofactors receive correct VdW radii
// and SYBYL atom types after the HETATM else-block was filled in.
// Apache-2.0 © 2026 Le Bonhomme Pharma

#include <gtest/gtest.h>
#include "../LIB/flexaid.h"
#include <cstring>
#include <cstdlib>

// Declarations of the functions under test
void assign_radii(atom* atoms, resid* residue, int atm_cnt);
void assign_types(FA_Global* FA, atom* atoms, resid* residue, char aminofile[]);

// ===========================================================================
// Helpers
// ===========================================================================

static FA_Global make_fa(int ntypes = 40) {
    FA_Global fa;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wnontrivial-memcall"
    memset(&fa, 0, sizeof(fa));
    #pragma clang diagnostic pop
    fa.ntypes = ntypes;
    fa.res_cnt = 1;  // one residue at index 1
    return fa;
}

// Build a minimal single-atom HETATM residue (type=1).
// Returns atom index 1 and residue index 1 inside the provided arrays.
static void make_ion_residue(atom* atoms, resid* residue,
                              const char* resname, const char* atomname)
{
    memset(&atoms[1], 0, sizeof(atom));
    memset(&residue[1], 0, sizeof(resid));

    strncpy(residue[1].name, resname, 3);  residue[1].name[3] = '\0';
    residue[1].type = 1;  // HETATM
    residue[1].fatm  = (int*)calloc(2, sizeof(int));
    residue[1].latm  = (int*)calloc(2, sizeof(int));
    residue[1].fatm[0] = 1;
    residue[1].latm[0] = 1;

    strncpy(atoms[1].name, atomname, 4);  atoms[1].name[4] = '\0';
    atoms[1].ofres  = 1;
    atoms[1].radius = 0.0f;
    atoms[1].type   = 39;  // dummy default (ntypes-1 for ntypes=40)
    atoms[1].optres = nullptr;
}

static void free_residue_arrays(resid* residue) {
    free(residue[1].fatm);
    free(residue[1].latm);
}

// ===========================================================================
// assign_radii — ion table
// ===========================================================================

struct IonRadiusCase { const char* resname; const char* atomname; float expected_r; };

class IonRadiusTest : public ::testing::TestWithParam<IonRadiusCase> {};

TEST_P(IonRadiusTest, CorrectRadius) {
    auto [resname, atomname, expected_r] = GetParam();

    atom   atoms[2];
    resid  residue[2];
    make_ion_residue(atoms, residue, resname, atomname);

    assign_radii(atoms, residue, /*atm_cnt=*/1);

    EXPECT_NEAR(atoms[1].radius, expected_r, 0.01f)
        << "Residue " << resname << " got wrong radius";

    free_residue_arrays(residue);
}

INSTANTIATE_TEST_SUITE_P(Ions, IonRadiusTest, ::testing::Values(
    IonRadiusCase{"MG ", " MG ", 1.73f},
    IonRadiusCase{"ZN ", " ZN ", 1.39f},
    IonRadiusCase{"CA ", " CA ", 1.74f},
    IonRadiusCase{"NA ", " NA ", 2.27f},
    IonRadiusCase{"K  ", " K  ", 2.75f},
    IonRadiusCase{"FE ", " FE ", 1.47f},
    IonRadiusCase{"FE2", " FE ", 1.47f},
    IonRadiusCase{"CU ", " CU ", 1.40f},
    IonRadiusCase{"MN ", " MN ", 1.61f},
    IonRadiusCase{"CL ", " CL ", 1.75f},
    IonRadiusCase{"BR ", " BR ", 1.85f}
));

// ===========================================================================
// assign_radii — organic cofactor element fallback
// ===========================================================================

TEST(OrganicCofactorRadius, CarbonAtom) {
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "HEM", " C1 ");
    assign_radii(atoms, residue, 1);
    EXPECT_NEAR(atoms[1].radius, 1.88f, 0.01f);
    free_residue_arrays(residue);
}

TEST(OrganicCofactorRadius, NitrogenAtom) {
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "FAD", " N1 ");
    assign_radii(atoms, residue, 1);
    EXPECT_NEAR(atoms[1].radius, 1.64f, 0.01f);
    free_residue_arrays(residue);
}

TEST(OrganicCofactorRadius, OxygenAtom) {
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "ATP", " O1 ");
    assign_radii(atoms, residue, 1);
    EXPECT_NEAR(atoms[1].radius, 1.42f, 0.01f);
    free_residue_arrays(residue);
}

TEST(OrganicCofactorRadius, SulfurAtom) {
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "COA", " S1 ");
    assign_radii(atoms, residue, 1);
    EXPECT_NEAR(atoms[1].radius, 1.77f, 0.01f);
    free_residue_arrays(residue);
}

TEST(OrganicCofactorRadius, PhosphorusAtom) {
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "ATP", " P  ");
    assign_radii(atoms, residue, 1);
    EXPECT_NEAR(atoms[1].radius, 1.80f, 0.01f);
    free_residue_arrays(residue);
}

TEST(OrganicCofactorRadius, UnknownFallback) {
    atom  atoms[2]; resid residue[2];
    // 'X' in position 1 of name → default 1.50 Å
    make_ion_residue(atoms, residue, "UNK", " X1 ");
    assign_radii(atoms, residue, 1);
    EXPECT_NEAR(atoms[1].radius, 1.50f, 0.01f);
    free_residue_arrays(residue);
}

// ===========================================================================
// assign_radii — ATOM record (protein) is unaffected
// ===========================================================================

TEST(ProteinRadiusUnchanged, CarbonAlpha) {
    atom  atoms[2]; resid residue[2];
    memset(&atoms[1], 0, sizeof(atom));
    memset(&residue[1], 0, sizeof(resid));
    strncpy(residue[1].name, "ALA", 3);  residue[1].name[3] = '\0';
    residue[1].type = 0;  // ATOM (protein)
    residue[1].fatm = (int*)calloc(2, sizeof(int));
    residue[1].latm = (int*)calloc(2, sizeof(int));
    residue[1].fatm[0] = 1; residue[1].latm[0] = 1;
    strncpy(atoms[1].name, " CA ", 4);  atoms[1].name[4] = '\0';
    atoms[1].ofres = 1; atoms[1].radius = 0.0f;

    assign_radii(atoms, residue, 1);

    // Cα is an aliphatic carbon → radius[2] = 1.88 Å (C4)
    EXPECT_NEAR(atoms[1].radius, 1.88f, 0.01f);

    free(residue[1].fatm); free(residue[1].latm);
}

// ===========================================================================
// assign_types — SYBYL ion type assignment
// ===========================================================================

// assign_types reads from AMINO.def — pass an empty (but valid) file
// so the AMINO.def loop is a no-op and only the ion-type loop runs.
#include <cstdio>

static std::string write_empty_aminodef() {
    char tmp[] = "/tmp/test_amino_XXXXXX.def";
    // mkstemps for .def suffix
    int fd = mkstemps(tmp, 4);
    if (fd >= 0) { close(fd); }
    return std::string(tmp);
}

struct IonTypeCase { const char* resname; const char* atomname; int expected_sybyl; };

class IonTypeTest : public ::testing::TestWithParam<IonTypeCase> {};

TEST_P(IonTypeTest, CorrectSybylType) {
    auto [resname, atomname, expected_sybyl] = GetParam();

    FA_Global fa = make_fa(/*ntypes=*/40);
    atom   atoms[2];
    resid  residue[2];
    make_ion_residue(atoms, residue, resname, atomname);

    std::string aminofile = write_empty_aminodef();
    assign_types(&fa, atoms, residue, const_cast<char*>(aminofile.c_str()));
    std::remove(aminofile.c_str());

    EXPECT_EQ(atoms[1].type, expected_sybyl)
        << "Residue " << resname << " got wrong SYBYL type";

    free_residue_arrays(residue);
}

INSTANTIATE_TEST_SUITE_P(IonTypes, IonTypeTest, ::testing::Values(
    IonTypeCase{"MG ", " MG ", 28},
    IonTypeCase{"ZN ", " ZN ", 35},
    IonTypeCase{"CA ", " CA ", 36},
    IonTypeCase{"FE ", " FE ", 37},
    IonTypeCase{"FE2", " FE ", 37},
    IonTypeCase{"MN ", " MN ", 31},
    IonTypeCase{"CU ", " CU ", 30},
    IonTypeCase{"NI ", " NI ", 34}
));

TEST(IonTypeTest, FallbackWhenNtypesSmall) {
    // If FA->ntypes < required SYBYL type, type should NOT be overwritten
    FA_Global fa = make_fa(/*ntypes=*/20);
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "ZN ", " ZN ");
    atoms[1].type = 19;  // dummy = ntypes-1

    std::string aminofile = write_empty_aminodef();
    assign_types(&fa, atoms, residue, const_cast<char*>(aminofile.c_str()));
    std::remove(aminofile.c_str());

    // SYBYL 35 > ntypes=20, so the ion block should skip assignment
    EXPECT_EQ(atoms[1].type, 19);  // unchanged
    free_residue_arrays(residue);
}

// Water oxygen must be typed O.3 (canonical VCT row 14), the same row as
// SER-OG / THR-OG1 / TYR-OH. This previously asserted type 1, which under the
// canonical numbering is C.1 — *sp carbon*, whose 1-13 (-198.3) and 1-14
// (-180.8) cells are the most attractive in MC_st0r5.2_6.dat. Every retained
// crystallographic water therefore radiated a large spurious attraction to
// ligand and protein oxygen. Row 40 (SOLVENT) is not an alternative: it is the
// reserved bulk-solvent/SAS pseudo-type and carries strongly repulsive entries.
TEST(IonTypeTest, WaterOxygenGetsO3) {
    FA_Global fa = make_fa(40);
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "HOH", " O  ");
    strncpy(atoms[1].element, "O", 1); atoms[1].element[1] = '\0';
    atoms[1].type = 39;

    std::string aminofile = write_empty_aminodef();
    assign_types(&fa, atoms, residue, const_cast<char*>(aminofile.c_str()));
    std::remove(aminofile.c_str());

    EXPECT_EQ(atoms[1].type, 14);  // O.3
    EXPECT_NE(atoms[1].type, 1);   // never C.1 — sp carbon
    free_residue_arrays(residue);
}

// A "WAT" water must be typed identically to an "HOH" water. Before the fix
// only the literal "HOH" was special-cased, so the two spellings took
// different paths through typing.
TEST(IonTypeTest, WaterResidueAliasesAgreeWithHOH) {
    for (const char* rname : {"WAT", "H2O", "DOD", "OHX"}) {
        FA_Global fa = make_fa(40);
        atom  atoms[2]; resid residue[2];
        make_ion_residue(atoms, residue, rname, " O  ");
        atoms[1].type = 39;

        std::string aminofile = write_empty_aminodef();
        assign_types(&fa, atoms, residue, const_cast<char*>(aminofile.c_str()));
        std::remove(aminofile.c_str());

        EXPECT_EQ(atoms[1].type, 14) << "water residue " << rname;
        free_residue_arrays(residue);
    }
}

// Water hydrogen (and deuterium, in neutron structures) must stay DUMMY —
// hydrogen is not scored against the VCT matrix.
TEST(IonTypeTest, WaterHydrogenStaysDummy) {
    FA_Global fa = make_fa(40);
    atom  atoms[2]; resid residue[2];
    make_ion_residue(atoms, residue, "HOH", " H1 ");
    atoms[1].type = 39;
    strncpy(atoms[1].element, "H", sizeof(atoms[1].element) - 1);

    std::string aminofile = write_empty_aminodef();
    assign_types(&fa, atoms, residue, const_cast<char*>(aminofile.c_str()));
    std::remove(aminofile.c_str());

    EXPECT_EQ(atoms[1].type, 39);  // DUMMY, untouched
    free_residue_arrays(residue);
}

// ===========================================================================
// tENCoM ion contact surface tests
//
// Build a synthetic receptor: N_PROT protein residues (Cα-only, linear chain
// along X-axis at 4 Å spacing) + 1 MG ion 6 Å from the last Cα.
// After assign_radii, MG gets radius=1.73 Å.  Surface distance =
//   6.0 − 1.73 (MG) − 1.88 (Cα) = 2.39 Å < 9 Å cutoff → spring is formed.
// ===========================================================================

#include "../LIB/tENCoM/tencm.h"
#include "../LIB/ion_utils.h"

namespace {

static constexpr int N_PROT = 10;  // need > 6 protein Cα so non-rigid modes exist
static constexpr float CA_SPACING = 4.0f;
static constexpr float MG_RADIUS  = 1.73f;
static constexpr float MG_DIST    = 6.0f;  // center-to-center from last Cα

// Helper: allocate fatm/latm arrays for one residue (rotamer 0 only)
static void alloc_res(resid* r, int first, int last) {
    r->fatm = (int*)calloc(2, sizeof(int));
    r->latm = (int*)calloc(2, sizeof(int));
    r->fatm[0] = first;
    r->latm[0] = last;
    r->rot  = 0;
}
static void free_res(resid* r) { free(r->fatm); free(r->latm); }

// Build a minimal atom/resid array with N_PROT protein Cα + 1 MG ion.
// atoms[1..N_PROT]  → protein residues (type=0), each with one " CA " atom
// atoms[N_PROT+1]   → MG ion (type=1, residue name "MG ")
// res_cnt = N_PROT + 1
static void build_synthetic(atom* atoms, resid* residues) {
    memset(atoms,   0, sizeof(atom)  * (N_PROT + 2));
    memset(residues,0, sizeof(resid) * (N_PROT + 2));

    // Protein Cα residues — zigzag ±1 Å in Y so pseudo-bond axes are non-collinear
    // (avoids degenerate linear chain where all cross products are zero)
    for (int i = 1; i <= N_PROT; ++i) {
        residues[i].type = 0;   // ATOM / protein
        strncpy(residues[i].name, "ALA", 3); residues[i].name[3] = '\0';
        alloc_res(&residues[i], i, i);

        strncpy(atoms[i].name, " CA ", 4); atoms[i].name[4] = '\0';
        atoms[i].ofres   = i;
        atoms[i].radius  = 1.88f;   // aliphatic C (pre-assigned)
        atoms[i].coor[0] = static_cast<float>(i - 1) * CA_SPACING;
        atoms[i].coor[1] = (i % 2 == 0) ? 1.0f : -1.0f;  // zigzag: non-linear
        atoms[i].coor[2] = 0.0f;
    }

    // MG ion (HETATM, residue N_PROT+1) — placed 6 Å from last Cα along X
    const int mg_idx = N_PROT + 1;
    const float last_ca_x = static_cast<float>(N_PROT - 1) * CA_SPACING;
    const float last_ca_y = (N_PROT % 2 == 0) ? 1.0f : -1.0f;

    residues[mg_idx].type = 1;  // HETATM
    strncpy(residues[mg_idx].name, "MG ", 3); residues[mg_idx].name[3] = '\0';
    alloc_res(&residues[mg_idx], mg_idx, mg_idx);

    strncpy(atoms[mg_idx].name, " MG ", 4); atoms[mg_idx].name[4] = '\0';
    atoms[mg_idx].ofres   = mg_idx;
    atoms[mg_idx].radius  = MG_RADIUS;
    atoms[mg_idx].coor[0] = last_ca_x + MG_DIST;
    atoms[mg_idx].coor[1] = last_ca_y;   // same Y as last Cα → center-to-center = MG_DIST
    atoms[mg_idx].coor[2] = 0.0f;
}

static void free_synthetic(resid* residues) {
    for (int i = 1; i <= N_PROT + 1; ++i) free_res(&residues[i]);
}

} // anonymous namespace

// After build(), n_protein_ca() must equal N_PROT and the model must be built.
TEST(TENCoMIonContact, IonNodePresent) {
    atom   atoms[N_PROT + 2];
    resid  residues[N_PROT + 2];
    build_synthetic(atoms, residues);

    tencm::TorsionalENM tenm;
    tenm.build(atoms, residues, N_PROT + 1);

    ASSERT_TRUE(tenm.is_built());
    EXPECT_EQ(tenm.n_protein_ca(), N_PROT);
    // Total nodes = protein Cα + 1 ion
    EXPECT_EQ(tenm.n_residues(), N_PROT + 1);

    free_synthetic(residues);
}

// tmcontsct() must contain at least one entry flagged as ion contact.
TEST(TENCoMIonContact, TmContSctHasIonEntry) {
    atom   atoms[N_PROT + 2];
    resid  residues[N_PROT + 2];
    build_synthetic(atoms, residues);

    tencm::TorsionalENM tenm;
    tenm.build(atoms, residues, N_PROT + 1);
    ASSERT_TRUE(tenm.is_built());

    bool found_ion = false;
    for (const auto& c : tenm.tmcontsct())
        found_ion |= c.is_ion;
    EXPECT_TRUE(found_ion) << "No ion contact found in tmcontsct()";

    free_synthetic(residues);
}

// Ion spring constant must be positive, finite, and scaled by area factor.
// For MG at d_surf = 6.0 - 1.73 - 1.88 = 2.39 Å:
//   area_scale = (1.73/1.88)^2 ≈ 0.848
//   k_raw = k0 * (9/2.39)^6 ≈ 3400 * k0
//   k_scaled ≈ 2880 * k0  — much larger than typical protein-protein spring
TEST(TENCoMIonContact, IonSpringKIsPositiveAndFinite) {
    atom   atoms[N_PROT + 2];
    resid  residues[N_PROT + 2];
    build_synthetic(atoms, residues);

    tencm::TorsionalENM tenm;
    tenm.build(atoms, residues, N_PROT + 1);
    ASSERT_TRUE(tenm.is_built());

    float max_ion_k = 0.0f;
    int   n_ion_contacts = 0;
    for (const auto& c : tenm.tmcontsct()) {
        if (!c.is_ion) continue;
        ++n_ion_contacts;
        EXPECT_GT(c.k_scaled, 0.0f)           << "Ion spring k_scaled must be positive";
        EXPECT_TRUE(std::isfinite(c.k_scaled)) << "Ion spring k_scaled must be finite";
        EXPECT_GT(c.d_surf, 0.0f)             << "Ion surface distance must be positive";
        if (c.k_scaled > max_ion_k) max_ion_k = c.k_scaled;
    }
    EXPECT_GT(n_ion_contacts, 0) << "At least one protein-ion contact expected";
    // The closest Cα (d_surf=2.39 Å) should produce a very stiff spring:
    //   k = k0*(9/2.39)^6*(1.73/1.88)^2 ≈ 2400; check the maximum is > 100.
    EXPECT_GT(max_ion_k, 100.0f) << "Maximum ion spring too weak for d_surf=2.39 Å";

    free_synthetic(residues);
}

// bfactors() must return exactly n_protein_ca() values (not n_residues()).
TEST(TENCoMIonContact, BFactorLengthEqualsProteinCaCount) {
    atom   atoms[N_PROT + 2];
    resid  residues[N_PROT + 2];
    build_synthetic(atoms, residues);

    tencm::TorsionalENM tenm;
    tenm.build(atoms, residues, N_PROT + 1);
    ASSERT_TRUE(tenm.is_built());

    auto bf = tenm.bfactors(300.0f);
    EXPECT_EQ(static_cast<int>(bf.size()), tenm.n_protein_ca());
    // All B-factors must be non-negative and finite
    for (float b : bf) {
        EXPECT_GE(b, 0.0f);
        EXPECT_TRUE(std::isfinite(b));
    }

    free_synthetic(residues);
}

// The ion acts as a restoring anchor: residues near the ion should have
// lower B-factors than residues far from it (ion stiffens the local network).
TEST(TENCoMIonContact, IonReducesBFactorNearIon) {
    atom   atoms[N_PROT + 2];
    resid  residues[N_PROT + 2];
    build_synthetic(atoms, residues);

    // Build with ion (MG at end)
    tencm::TorsionalENM tenm_with_ion;
    tenm_with_ion.build(atoms, residues, N_PROT + 1);

    // Build without ion (protein only)
    tencm::TorsionalENM tenm_no_ion;
    tenm_no_ion.build(atoms, residues, N_PROT);  // exclude ion residue

    ASSERT_TRUE(tenm_with_ion.is_built());
    ASSERT_TRUE(tenm_no_ion.is_built());

    auto bf_ion    = tenm_with_ion.bfactors(300.0f);
    auto bf_no_ion = tenm_no_ion.bfactors(300.0f);

    ASSERT_EQ(bf_ion.size(), bf_no_ion.size());

    // The last Cα (closest to Mg) should be stiffer with the ion present
    const int last = static_cast<int>(bf_ion.size()) - 1;
    {
        double mx_i = 0, mx_n = 0;
        for (float v : bf_ion) mx_i = std::max<double>(mx_i, v);
        for (float v : bf_no_ion) mx_n = std::max<double>(mx_n, v);
        const auto& mi = tenm_with_ion.modes();
        double lmin = mi.empty()?-1:mi.front().eigenvalue, lmax = mi.empty()?-1:mi.back().eigenvalue;
        std::cout << "PROBE bf_size=" << bf_ion.size() << " n_modes=" << mi.size()
                  << " n_bonds=" << tenm_with_ion.n_bonds()
                  << " lam_min=" << lmin << " lam_max=" << lmax
                  << " max_bf_ion=" << mx_i << " max_bf_noion=" << mx_n
                  << " bf_ion[last]=" << bf_ion[last] << std::endl;
    }
    EXPECT_LT(bf_ion[last], bf_no_ion[last])
        << "Residue adjacent to Mg should have lower B-factor with ion present";

    free_synthetic(residues);
}
