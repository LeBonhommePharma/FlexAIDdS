// test_pb_vdw_parity.cpp — Parity between the two PoseBusters vdW radius sources
//
// The pb_clash scan in vcfunction.cpp derives an atom's PoseBusters vdW radius
// from get_element(atoms[i].type). The atom struct also carries a precomputed
// atoms[i].pb_vdw_radius, populated in update_optres() from atoms[i].element
// (the structure file's element column). Vcontacts.cpp already reads the cached
// field; vcfunction.cpp reads it only under FLEXAIDDS_PB_VDW_CACHED.
//
// The two are NOT interchangeable, and the reason is the reader tables, not the
// PoseBusters table. `element` records the true element. `type` is an NRGDock
// VCT slot, and the readers deliberately map some elements onto a DIFFERENT
// element's slot (iodine -> the bromine row) or onto the DUMMY slot (anything
// the 40-type table does not represent). Going type -> element via get_element()
// therefore LOSES the true element for exactly those atoms.
//
// This test encodes the real reader tables (SdfReader::element_to_flexaid_type /
// element_radius and Mol2Reader::sybyl_to_flexaid_type / sybyl_radius, both read
// from source) and pins, per element, whether the two radius sources agree. It
// exists so the divergent set is a checked fact rather than an assumption.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
// SPDX-License-Identifier: Apache-2.0

#include <gtest/gtest.h>

#include <cstring>
#include <string>

#include "soft_wall.h"

extern "C++" const char* get_element(int type);

namespace {

// Strip the leading blank get_element() uses for one-letter elements, matching
// what the pb_clash scan does before the table lookup.
std::string type_element(int type) {
    const char* e = get_element(type);
    while (*e == ' ') ++e;
    return std::string(e);
}

// A ligand atom as the readers actually construct it: the element column they
// write, the VCT type they assign, and the NRG radius they assign (which is the
// fallback the type-derived lookup lands on when get_element() yields something
// outside the PoseBusters table).
struct ReaderAtom {
    const char* element;    // atoms[].element
    int         vct;        // atoms[].type
    double      nrg_radius; // atoms[].radius
};

// Mirrors SdfReader.cpp element_to_flexaid_type() + element_radius() exactly.
// The MOL2 path (sybyl_to_flexaid_type / sybyl_radius) agrees on every element
// below; only its sub-type discrimination differs, which never changes the
// ELEMENT that get_element() reports.
const ReaderAtom kReaderAtoms[] = {
    {"C",  3, 1.70}, {"N", 11, 1.55}, {"O", 14, 1.52}, {"S", 18, 1.80},
    {"P", 22, 1.80}, {"F", 23, 1.47}, {"Cl", 24, 1.75}, {"Br", 25, 1.85},
    {"I",  25, 1.98},                      // <-- iodine mapped onto the BROMINE row
    {"Se", 27, 1.90},
    {"Mg", 28, 1.73}, {"Sr", 29, 1.70}, {"Cu", 30, 1.70}, {"Mn", 31, 1.70},
    {"Hg", 32, 1.70}, {"Cd", 33, 1.70}, {"Ni", 34, 1.70}, {"Zn", 35, 1.22},
    {"Ca", 36, 2.31}, {"Fe", 37, 1.34}, {"Co", 38, 1.70},
    {"H",  39, 1.20},                      // <-- DUMMY slot
    {"Na", 39, 1.70}, {"K", 39, 1.70},     // <-- absent from the 40-type table
};

double from_element(const ReaderAtom& a) {
    return posebusters_vdw_radius(a.element, a.nrg_radius);
}
double from_type(const ReaderAtom& a) {
    return posebusters_vdw_radius(type_element(a.vct), a.nrg_radius);
}

bool is_expected_divergence(const char* el) {
    return std::strcmp(el, "I")  == 0 ||
           std::strcmp(el, "Na") == 0 ||
           std::strcmp(el, "K")  == 0;
}

}  // namespace

// Every element the readers map onto its OWN element's VCT row agrees on both
// paths, so switching those to the cached field is a pure no-op.
TEST(PbVdwParity, ElementPreservingTypesAgreeOnBothPaths) {
    for (const auto& a : kReaderAtoms) {
        if (is_expected_divergence(a.element)) continue;
        EXPECT_DOUBLE_EQ(from_type(a), from_element(a))
            << "element=" << a.element << " vct=" << a.vct
            << " get_element(vct)='" << type_element(a.vct) << "'";
    }
}

// Hydrogen lands in the DUMMY slot, so the type path falls back to the NRG
// radius — but the readers assign H a NRG radius of exactly 1.20, which is also
// the PoseBusters hydrogen radius. The two paths therefore COINCIDE for H.
// This corrects an earlier reading of this code that named type 39 as the
// divergence: the DUMMY slot only diverges when the fallback radius happens to
// differ from the PoseBusters value, which for H it does not.
TEST(PbVdwParity, DummySlotHydrogenCoincidesBecauseNrgRadiusIsAlso120) {
    const ReaderAtom h{"H", 39, 1.20};
    EXPECT_EQ(type_element(39), "Du");
    EXPECT_DOUBLE_EQ(posebusters_vdw_radius("H", 1.20), 1.20);
    EXPECT_DOUBLE_EQ(from_type(h), from_element(h));
}

// The real divergences. Iodine is mapped onto the bromine row on purpose (the
// iodine row is nearly untrained), and Na/K have no row at all. In all three the
// element column keeps the true element and the type column does not, so
// enabling FLEXAIDDS_PB_VDW_CACHED changes these atoms' clash cutoff.
TEST(PbVdwParity, IodineAndAlkaliDivergeAndAreTheOnlyOnes) {
    const ReaderAtom iodine{"I", 25, 1.98};
    EXPECT_EQ(type_element(25), "Br");
    EXPECT_DOUBLE_EQ(from_type(iodine), 1.90);      // bromine radius, via the type
    EXPECT_DOUBLE_EQ(from_element(iodine), 2.10);   // iodine radius, via the element

    const ReaderAtom na{"Na", 39, 1.70};
    EXPECT_DOUBLE_EQ(from_type(na), 1.70);          // DUMMY -> NRG fallback
    EXPECT_DOUBLE_EQ(from_element(na), 2.40);       // sodium radius, via the element

    const ReaderAtom k{"K", 39, 1.70};
    EXPECT_DOUBLE_EQ(from_type(k), 1.70);
    EXPECT_DOUBLE_EQ(from_element(k), 2.80);

    // And nothing else in the reader table diverges.
    for (const auto& a : kReaderAtoms) {
        const bool diverges = (from_type(a) != from_element(a));
        EXPECT_EQ(diverges, is_expected_divergence(a.element))
            << "element=" << a.element << " type=" << from_type(a)
            << " elem=" << from_element(a);
    }
}

// get_element() itself must stay a faithful 40-slot table: every VCT row the
// readers can assign has to report the element that row represents.
TEST(PbVdwParity, GetElementCoversAllFortyVctSlots) {
    const struct { int type; const char* elem; } kSlots[] = {
        {1,"C"},{2,"C"},{3,"C"},{4,"C"},{5,"C"},
        {6,"N"},{7,"N"},{8,"N"},{9,"N"},{10,"N"},{11,"N"},{12,"N"},
        {13,"O"},{14,"O"},{15,"O"},{16,"O"},
        {17,"S"},{18,"S"},{19,"S"},{20,"S"},{21,"S"},
        {22,"P"},{23,"F"},{24,"Cl"},{25,"Br"},{26,"I"},{27,"Se"},
        {28,"Mg"},{29,"Sr"},{30,"Cu"},{31,"Mn"},{32,"Hg"},{33,"Cd"},
        {34,"Ni"},{35,"Zn"},{36,"Ca"},{37,"Fe"},{38,"Co"},{39,"Du"},
    };
    for (const auto& s : kSlots)
        EXPECT_EQ(type_element(s.type), s.elem) << "vct slot " << s.type;

    // Slot 40 (SOLVENT) and out-of-range carry no element.
    for (int t : {0, 40, 99, -1})
        EXPECT_EQ(type_element(t), "") << "type " << t;
}

// The carve-out predicate must select exactly the documented metal set.
TEST(PbVdwParity, CoordinatingMetalSetIsExact) {
    for (const char* m : {"Zn","Fe","Mg","Ca","Mn","Co","Ni","Cu"})
        EXPECT_TRUE(posebusters_is_coordinating_metal(m)) << m;
    for (const char* nm : {"C","N","O","S","P","F","Cl","Br","I","Se","H","Na","K","Sr","Hg","Cd","Du",""})
        EXPECT_FALSE(posebusters_is_coordinating_metal(nm)) << nm;
}
