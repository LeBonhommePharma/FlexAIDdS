// test_hbond_amine_roles.cpp — topology-based H-bond role classification
//
// Regression cover for the ligand sp3-nitrogen role defect: classify_hbond_donor
// returned false unconditionally for N_sp3, and classify_hbond_acceptor keyed on
// `partial_charge < 0.3f`. Ligand partial charges are identically 0 on the
// PDB-derived SDF path, so that test was always true and EVERY ligand sp3
// nitrogen entered scoring acceptor-only — a protonated amine had its donor and
// acceptor roles exactly inverted.
//
// The contract asserted here: amine substitution is resolved by heavy-atom
// topology, and the verdict is INDEPENDENT of partial charge. Every role
// assertion below is swept across a range of partial charges spanning the old
// 0.3f threshold; a classifier that reads charge cannot pass these.

#include <gtest/gtest.h>
#include <utility>
#include "../LIB/atom_typing_256.h"

namespace {

// Charges spanning the old `< 0.3f` acceptor threshold in both directions,
// including the identically-zero value the real SDF path supplies.
constexpr float kChargeSweep[] = {-0.8f, -0.3f, 0.0f, 0.29f, 0.3f, 0.5f, 1.0f};

atom256::HbondTopology topo(int n_heavy) {
    atom256::HbondTopology t;
    t.n_heavy_neighbors = n_heavy;
    t.known             = true;
    return t;
}

// SYBYL type numbers used below (1-indexed, see sybyl_to_base).
constexpr int SYBYL_N3 = 8;   // N.3
constexpr int SYBYL_N4 = 9;   // N.4
constexpr int SYBYL_O3 = 14;  // O.3

}  // namespace

// ── the headline defect ─────────────────────────────────────────────────────

TEST(HbondAmineRoles, ProtonatedPrimaryAmineIsDonor) {
    // R–NH3+ : one heavy substituent, three hydrogens. This is the dominant CNS
    // pharmacophore, and it was previously classified acceptor-only.
    for (float q : kChargeSweep) {
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_N4, q, /*n_hydrogens=*/3, topo(1));
        EXPECT_EQ(atom256::get_base(code), atom256::N_quat) << "q=" << q;
        EXPECT_TRUE(atom256::get_hbond_donor(code))
            << "protonated primary amine must donate, q=" << q;
        EXPECT_FALSE(atom256::get_hbond_acceptor(code))
            << "protonated amine has no free lone pair, q=" << q;
    }
}

TEST(HbondAmineRoles, ProtonatedAmineTypedNsp3IsDonor) {
    // Same chemistry arriving as N.3 with an explicit H count — the SDF path,
    // where no N.4 typing is available.
    for (float q : kChargeSweep) {
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_N3, q, /*n_hydrogens=*/2, topo(1));
        EXPECT_EQ(atom256::get_base(code), atom256::N_sp3) << "q=" << q;
        EXPECT_TRUE(atom256::get_hbond_donor(code)) << "q=" << q;
        // Neutral primary amine is amphoteric: N–H donates, lone pair accepts.
        EXPECT_TRUE(atom256::get_hbond_acceptor(code)) << "q=" << q;
    }
}

TEST(HbondAmineRoles, TertiaryAmineIsAcceptorOnly) {
    // R3N : three heavy substituents, no N–H. Acceptor, never donor.
    for (float q : kChargeSweep) {
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_N3, q, /*n_hydrogens=*/0, topo(3));
        EXPECT_EQ(atom256::get_base(code), atom256::N_sp3) << "q=" << q;
        EXPECT_FALSE(atom256::get_hbond_donor(code))
            << "tertiary amine has no labile H, q=" << q;
        EXPECT_TRUE(atom256::get_hbond_acceptor(code)) << "q=" << q;
    }
}

TEST(HbondAmineRoles, SecondaryAmineIsBothRoles) {
    for (float q : kChargeSweep) {
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_N3, q, /*n_hydrogens=*/1, topo(2));
        EXPECT_TRUE(atom256::get_hbond_donor(code)) << "q=" << q;
        EXPECT_TRUE(atom256::get_hbond_acceptor(code)) << "q=" << q;
    }
}

TEST(HbondAmineRoles, QuaternaryAmmoniumIsNeitherRole) {
    // R4N+ : four heavy substituents. No labile H, no free lone pair.
    for (float q : kChargeSweep) {
        for (int sybyl : {SYBYL_N3, SYBYL_N4}) {
            const uint8_t code = atom256::encode_from_sybyl(
                sybyl, q, /*n_hydrogens=*/0, topo(4));
            EXPECT_FALSE(atom256::get_hbond_donor(code))
                << "sybyl=" << sybyl << " q=" << q;
            EXPECT_FALSE(atom256::get_hbond_acceptor(code))
                << "sybyl=" << sybyl << " q=" << q;
        }
    }
}

TEST(HbondAmineRoles, ProtonatedTertiaryAmineIsDonorOnlyNotAmphoteric) {
    // R3NH+ : three heavy substituents plus one explicit H, formal charge NOT
    // known (the SDF path supplies none). Coordination is 4, so the lone pair
    // was spent forming the N-H — this is an ammonium, and it must NOT also
    // advertise an acceptor role.
    //
    // Regression guard: an earlier revision let the donor test accept three
    // independent forms of protonation evidence while the acceptor test
    // consulted only the perceived formal charge. This case came out donor AND
    // acceptor simultaneously, with a lone pair that does not exist. It is the
    // R-NH+ centre this whole fix exists for.
    for (float q : kChargeSweep) {
        atom256::HbondTopology t = topo(3);
        ASSERT_FALSE(t.charge_known);
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_N3, q, /*n_hydrogens=*/1, t);
        EXPECT_TRUE(atom256::get_hbond_donor(code))
            << "protonated tertiary amine must donate, q=" << q;
        EXPECT_FALSE(atom256::get_hbond_acceptor(code))
            << "lone pair is spent — must not also accept, q=" << q;
    }
}

TEST(HbondAmineRoles, ProtonatedPrimaryAndSecondaryAminesAreDonorOnly) {
    // Conjugate acids of 1° and 2° amines, evidenced by explicit H alone.
    // R-NH3+ is (1 heavy + 3 H) and R2NH2+ is (2 heavy + 2 H): both reach
    // coordination 4 and must lose the acceptor role.
    for (float q : kChargeSweep) {
        for (auto hc : {std::pair<int,int>{1,3}, std::pair<int,int>{2,2}}) {
            const uint8_t code = atom256::encode_from_sybyl(
                SYBYL_N3, q, /*n_hydrogens=*/hc.second, topo(hc.first));
            EXPECT_TRUE(atom256::get_hbond_donor(code))
                << "heavy=" << hc.first << " q=" << q;
            EXPECT_FALSE(atom256::get_hbond_acceptor(code))
                << "heavy=" << hc.first << " q=" << q;
        }
    }
}

TEST(HbondAmineRoles, NeutralAminesRemainAmphoteric) {
    // The protonation test must not misfire on legitimately amphoteric neutral
    // amines. Neutral sp3 N is 3-coordinate at every substitution level:
    // R-NH2 (1+2), R2NH (2+1), R3N (3+0). All keep the acceptor role.
    for (float q : kChargeSweep) {
        for (auto hc : {std::pair<int,int>{1,2}, std::pair<int,int>{2,1},
                        std::pair<int,int>{3,0}}) {
            const uint8_t code = atom256::encode_from_sybyl(
                SYBYL_N3, q, /*n_hydrogens=*/hc.second, topo(hc.first));
            EXPECT_TRUE(atom256::get_hbond_acceptor(code))
                << "neutral amine must accept, heavy=" << hc.first
                << " nH=" << hc.second << " q=" << q;
            EXPECT_EQ(atom256::get_hbond_donor(code), hc.second > 0)
                << "donates iff it has N-H, heavy=" << hc.first << " q=" << q;
        }
    }
}

TEST(HbondAmineRoles, DonorAndAcceptorEvidenceAreSymmetric) {
    // Structural invariant: for sp3 N, every (heavy, nH) combination the donor
    // test reads as protonated must also be excluded by the acceptor test.
    //
    // The invariant is stated as "amphoteric implies coordination < 4", NOT
    // "coordination == 3". Those differ when the H count is absent rather than
    // genuinely zero: heavy=1, nH=0 is an amine whose hydrogens were never
    // supplied, and the donor test infers N-H structurally from the low heavy
    // count while the coordination sum under-counts. Falling back to amphoteric
    // there is the conservative answer — it keeps the acceptor role the atom
    // always had and adds the donor role it was missing. Protonation simply
    // cannot be perceived when neither explicit H nor formal charge is present.
    for (int heavy = 1; heavy <= 4; ++heavy) {
        for (int nh = 0; nh <= 4; ++nh) {
            for (float q : kChargeSweep) {
                const bool d = atom256::classify_hbond_donor(
                    atom256::N_sp3, q, nh, topo(heavy));
                const bool a = atom256::classify_hbond_acceptor(
                    atom256::N_sp3, q, nh, topo(heavy));
                if (d && a)
                    EXPECT_LT(heavy + nh, 4)
                        << "amphoteric requires a free lone pair; heavy="
                        << heavy << " nH=" << nh << " q=" << q;
                if (heavy + nh >= 4)
                    EXPECT_FALSE(a)
                        << "no lone pair at coordination " << (heavy + nh)
                        << "; heavy=" << heavy << " nH=" << nh;
            }
        }
    }
}

// ── sp3 oxygen: hydroxyl vs ether ───────────────────────────────────────────

TEST(HbondAmineRoles, EtherOxygenIsAcceptorOnly) {
    // R–O–R : two heavy substituents, no O–H.
    for (float q : kChargeSweep) {
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_O3, q, /*n_hydrogens=*/0, topo(2));
        EXPECT_EQ(atom256::get_base(code), atom256::O_sp3) << "q=" << q;
        EXPECT_FALSE(atom256::get_hbond_donor(code))
            << "ether has no O-H, q=" << q;
        EXPECT_TRUE(atom256::get_hbond_acceptor(code)) << "q=" << q;
    }
}

TEST(HbondAmineRoles, HydroxylOxygenIsBothRoles) {
    // R–OH : one heavy substituent.
    for (float q : kChargeSweep) {
        const uint8_t code =
            atom256::encode_from_sybyl(SYBYL_O3, q, /*n_hydrogens=*/1, topo(1));
        EXPECT_TRUE(atom256::get_hbond_donor(code)) << "q=" << q;
        EXPECT_TRUE(atom256::get_hbond_acceptor(code)) << "q=" << q;
    }
}

// ── byte-identity guarantee when topology is unavailable ────────────────────

TEST(HbondAmineRoles, UnknownTopologyReproducesLegacyVerdict) {
    // A default-constructed HbondTopology has known=false. Every base type, H
    // count and charge must then yield exactly the pre-fix classification, so
    // call sites without connectivity (e.g. PDB receptor atoms with bond[0]==0)
    // are bit-for-bit unchanged.
    const atom256::HbondTopology unknown{};
    ASSERT_FALSE(unknown.known);

    for (uint8_t base = 0; base < atom256::BASE_TYPE_COUNT; ++base) {
        for (float q : kChargeSweep) {
            for (int nh = 0; nh <= 3; ++nh) {
                EXPECT_EQ(
                    atom256::classify_hbond_donor(base, q, nh, unknown),
                    atom256::classify_hbond_donor(base, q, nh))
                    << "base=" << int(base) << " q=" << q << " nH=" << nh;
                EXPECT_EQ(
                    atom256::classify_hbond_acceptor(base, q, nh, unknown),
                    atom256::classify_hbond_acceptor(base, q, nh))
                    << "base=" << int(base) << " q=" << q << " nH=" << nh;
            }
        }
    }
}

TEST(HbondAmineRoles, NonAmineBaseTypesUnaffectedByTopology) {
    // The topology overload must only change verdicts for N_sp3, N_quat and
    // O_sp3. Everything else delegates unchanged regardless of heavy count.
    for (uint8_t base = 0; base < atom256::BASE_TYPE_COUNT; ++base) {
        if (base == atom256::N_sp3 || base == atom256::N_quat ||
            base == atom256::O_sp3)
            continue;
        for (int heavy = 0; heavy <= 4; ++heavy) {
            for (float q : kChargeSweep) {
                for (int nh = 0; nh <= 3; ++nh) {
                    EXPECT_EQ(
                        atom256::classify_hbond_donor(base, q, nh, topo(heavy)),
                        atom256::classify_hbond_donor(base, q, nh))
                        << "base=" << int(base) << " heavy=" << heavy;
                    EXPECT_EQ(
                        atom256::classify_hbond_acceptor(base, q, nh,
                                                         topo(heavy)),
                        atom256::classify_hbond_acceptor(base, q, nh))
                        << "base=" << int(base) << " heavy=" << heavy;
                }
            }
        }
    }
}
