// test_pb_vdw_parity.cpp — Parity between the two PoseBusters vdW radius sources
//
// The pb_clash scan in vcfunction.cpp derives an atom's PoseBusters vdW radius
// from get_element(atoms[i].type). The atom struct also carries a precomputed
// atoms[i].pb_vdw_radius, populated in update_optres() from atoms[i].element
// (the PDB/MOL2/SDF element column). Vcontacts.cpp already reads the cached
// field; vcfunction.cpp does not.
//
// Before vcfunction.cpp can switch to the cached field the two sources must be
// shown to agree. This test enumerates every NRGDock type and pins exactly where
// they agree and where they do not, so the divergence is a checked fact rather
// than an assumption.
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

// The element column a correctly-parsed structure file carries for an atom of
// each NRGDock type. Types 1-38 are the element families get_element() encodes;
// type 39 ("Du") is FlexAID's dummy/hydrogen tag, whose structure files carry a
// real "H" element.
struct TypeElement { int type; const char* pdb_element; };

const TypeElement kTypeElements[] = {
    {1,"C"},{2,"C"},{3,"C"},{4,"C"},{5,"C"},
    {6,"N"},{7,"N"},{8,"N"},{9,"N"},{10,"N"},{11,"N"},{12,"N"},
    {13,"O"},{14,"O"},{15,"O"},{16,"O"},
    {17,"S"},{18,"S"},{19,"S"},{20,"S"},{21,"S"},
    {22,"P"},{23,"F"},{24,"Cl"},{25,"Br"},{26,"I"},{27,"Se"},
    {28,"Mg"},{29,"Sr"},{30,"Cu"},{31,"Mn"},{32,"Hg"},{33,"Cd"},
    {34,"Ni"},{35,"Zn"},{36,"Ca"},{37,"Fe"},{38,"Co"},
};

constexpr double kFallback = 1.234;  // sentinel: distinguishable from every table entry

}  // namespace

// Types 1-38: the type-derived and element-derived radii must be identical, so
// vcfunction.cpp may read the precomputed atoms[].pb_vdw_radius for them.
TEST(PbVdwParity, TypeAndElementAgreeForAllHeavyTypes) {
    for (const auto& te : kTypeElements) {
        const double from_type    = posebusters_vdw_radius(type_element(te.type), kFallback);
        const double from_element = posebusters_vdw_radius(te.pdb_element, kFallback);
        EXPECT_DOUBLE_EQ(from_type, from_element)
            << "type " << te.type << " get_element='" << type_element(te.type)
            << "' pdb_element='" << te.pdb_element << "'";
        EXPECT_NE(from_type, kFallback)
            << "type " << te.type << " unexpectedly fell through to the fallback radius";
    }
}

// Type 39 is the one real divergence, and it is why the switch to the cached
// field must be flag-gated rather than applied unconditionally: get_element(39)
// returns "Du", which is absent from the table and falls back to the NRG contact
// radius, while the element column says "H" and yields the PoseBusters 1.20.
TEST(PbVdwParity, Type39HydrogenDivergesAndIsTheOnlyDivergence) {
    EXPECT_EQ(type_element(39), "Du");

    const double from_type    = posebusters_vdw_radius(type_element(39), kFallback);
    const double from_element = posebusters_vdw_radius("H", kFallback);

    EXPECT_DOUBLE_EQ(from_type, kFallback);   // "Du" is not in the table
    EXPECT_DOUBLE_EQ(from_element, 1.20);     // PoseBusters hydrogen radius
    EXPECT_NE(from_type, from_element);
}

// Out-of-range types carry no element and fall back on both paths, so they are
// not a divergence source.
TEST(PbVdwParity, UnknownTypesFallBackOnBothPaths) {
    for (int t : {0, 40, 99, -1}) {
        EXPECT_DOUBLE_EQ(posebusters_vdw_radius(type_element(t), kFallback), kFallback)
            << "type " << t;
    }
}

// The carve-out predicate must select exactly the documented metal set.
TEST(PbVdwParity, CoordinatingMetalSetIsExact) {
    for (const char* m : {"Zn","Fe","Mg","Ca","Mn","Co","Ni","Cu"})
        EXPECT_TRUE(posebusters_is_coordinating_metal(m)) << m;
    for (const char* nm : {"C","N","O","S","P","F","Cl","Br","I","Se","H","Na","K","Sr","Hg","Cd","Du",""})
        EXPECT_FALSE(posebusters_is_coordinating_metal(nm)) << nm;
}
