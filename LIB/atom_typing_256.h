// atom_typing_256.h — 8-bit atom type encoding for 256×256 soft contact matrix
//
// Extends FlexAID's 40-type SYBYL system to 256 types:
//   Bits 0–5: base type (64 classes, superset of SYBYL — no rare type collapse)
//   Bit    6: H-bond acceptor role
//   Bit    7: H-bond donor role
//
// All 40 SYBYL atom types map to distinct base types — no Solvent fallback.
// Includes sybyl_to_base() bridge from FlexAID's 40-type world, context-aware
// refinements for C_ar_hetadj and C_pi_bridging (NATURaL-critical for indole/
// tryptamine π-systems), and base_to_sybyl_parent() for the 256→40 projection.
#pragma once

#include <cstdint>
#include <cstring>
#include <array>
#include <cmath>

namespace atom256 {

// ─── base types (bits 0–5, 64 classes) ──────────────────────────────────────
enum BaseType : uint8_t {
    C_sp         =  0,   // C.1  — triple bond carbon
    C_sp2        =  1,   // C.2  — double bond carbon
    C_sp3        =  2,   // C.3  — alkane carbon
    C_ar         =  3,   // C.AR — aromatic carbon
    C_cat        =  4,   // C.CAT — cationic carbon
    N_sp         =  5,   // N.1
    N_sp2        =  6,   // N.2
    N_sp3        =  7,   // N.3
    N_quat       =  8,   // N.4  — quaternary ammonium
    N_ar         =  9,   // N.AR
    N_am         = 10,   // N.AM — amide nitrogen
    N_pl3        = 11,   // N.PL3 — planar sp3 nitrogen
    O_sp2        = 12,   // O.2  — carbonyl oxygen
    O_sp3        = 13,   // O.3  — hydroxyl / ether
    O_co2        = 14,   // O.CO2 — carboxylate
    O_ar         = 15,   // O.AR — aromatic oxygen
    S_sp2        = 16,   // S.2
    S_sp3        = 17,   // S.3
    S_oxide      = 18,   // S.O
    S_dioxide    = 19,   // S.O2
    S_ar         = 20,   // S.AR — aromatic sulfur
    P_sp3        = 21,   // P.3
    HAL_F        = 22,   // F
    HAL_Cl       = 23,   // CL
    HAL_Br       = 24,   // BR
    HAL_I        = 25,   // I
    C_ar_hetadj  = 26,   // aromatic C adjacent to heteroatom (indole C3a, C7a)
    C_pi_bridge  = 27,   // π-bridging carbon (tryptamine/indole bridge)
    Metal_Zn     = 28,   // ZN
    Metal_Ca     = 29,   // CA
    Metal_Fe     = 30,   // FE
    Solvent      = 31,   // SOL (true solvent probe only)
    // ─── extended types (32–63): no more Solvent fallback ────────────────
    HAL_Se       = 32,   // SE  — selenocysteine, distinct coordination
    Metal_Mg     = 33,   // MG  — 6 NRGRank interactions (O.2, O.3, O.CO2, P.3, …)
    Metal_Sr     = 34,   // SR
    Metal_Cu     = 35,   // CU  — distinct redox/coordination chemistry
    Metal_Mn     = 36,   // MN  — 1 NRGRank interaction (O.CO2)
    Metal_Hg     = 37,   // HG  — toxic heavy metal, distinct binding
    Metal_Cd     = 38,   // CD  — toxic heavy metal
    Metal_Ni     = 39,   // NI  — 1 NRGRank interaction (C.AR)
    Metal_Co     = 40,   // CO.OH — cobalt hydroxide
    Dummy        = 41,   // DUMMY — placeholder / unknown atom
    // 42–63: reserved for future types (RNA bases, lipid headgroups, PTMs, …)
    BASE_TYPE_COUNT = 64
};

// ─── encoding / decoding ────────────────────────────────────────────────────
// Layout: [D:1][A:1][B:6] = donor × acceptor × 64 base classes.

inline constexpr uint8_t encode_roles(uint8_t base_type, bool donor,
                                      bool acceptor) noexcept {
    return (static_cast<uint8_t>(donor) << 7) |
           (static_cast<uint8_t>(acceptor) << 6) |
           (base_type & 0x3F);
}

inline constexpr uint8_t encode(uint8_t base_type, bool donor,
                                bool acceptor) noexcept {
    return encode_roles(base_type, donor, acceptor);
}

// Catch stale v55 call sites that still try to encode [charge][H-bond].
inline constexpr uint8_t encode(uint8_t, uint8_t, bool) noexcept = delete;
inline constexpr uint8_t encode(uint8_t, int, bool) noexcept = delete;

inline constexpr uint8_t get_base(uint8_t code) noexcept { return code & 0x3F; }
inline constexpr bool    get_hbond_acceptor(uint8_t code) noexcept { return (code >> 6) & 0x01; }
inline constexpr bool    get_hbond_donor(uint8_t code) noexcept { return (code >> 7) & 0x01; }
inline constexpr bool    get_hbond(uint8_t code) noexcept {
    return get_hbond_donor(code) || get_hbond_acceptor(code);
}

// ─── H-bond role classification ─────────────────────────────────────────────
// Donor roles require attached H evidence. Runtime heavy-atom inputs provide
// this through explicit bonded H counts plus conservative implicit-H estimates
// computed before encode_from_sybyl(). This keeps the v57 donor/acceptor split
// without turning every polar heavy atom into a directional donor.

struct HbondRoles {
    bool donor;
    bool acceptor;
};

// ─── topology evidence for amine substitution ───────────────────────────────
// Partial charge cannot resolve amine substitution: ligand charges are
// identically 0 on the PDB-derived SDF path (PDB carries no partial charges by
// format definition), so `partial_charge < 0.3f` is unconditionally true and
// every ligand sp3 N entered scoring acceptor-only — a protonated amine, the
// dominant CNS pharmacophore, had its roles exactly inverted.
//
// The discriminator is heavy-atom substitution count, not charge and not
// n_hydrogens (which cannot tell "0 H because tertiary" from "0 H because the
// input had no H information"). `known` carries that distinction explicitly:
// when connectivity is absent the classifier falls back to the pre-existing
// charge-based verdict, byte for byte.
struct HbondTopology {
    int  n_heavy_neighbors = 0;      // bonded non-hydrogen count
    bool known             = false;  // connectivity actually available
    int  formal_charge     = 0;      // perceived formal charge (e)
    bool charge_known      = false;  // distinguishes "0" from "unknown"

    // Perceive formal charge from substitution where the topology allows it.
    // A 4-coordinate sp3 N is a quaternary/protonated ammonium (+1); nothing
    // else is inferable from heavy-atom count alone.
    constexpr bool is_quaternary_nitrogen() const noexcept {
        return known && n_heavy_neighbors >= 4;
    }
    constexpr bool is_cationic() const noexcept {
        return charge_known && formal_charge > 0;
    }
};

inline bool classify_hbond_donor(uint8_t base_type, float partial_charge,
                                 int n_hydrogens) noexcept {
    (void)partial_charge;
    switch (base_type) {
        case N_am:
            // Primary/secondary amide N (N-H present) is a good donor; tertiary
            // N.am (piperidine junction, N-methyl amide, ring-junction N) has no
            // labile H and must NOT be marked donor.
            return n_hydrogens > 0;
        case N_sp3:
            // Covers 1°/2°/3° amine. Tertiary amine (no N-H) dominates drug
            // scaffolds and n_hydrogens cannot distinguish it → conservative
            // acceptor-only default (no donor role).
            return false;
        case O_sp3:
            // Covers ether (no O-H) and hydroxyl. Ether dominates drug
            // scaffolds and n_hydrogens cannot distinguish it → conservative
            // acceptor-only default (no donor role).
            return false;
        case N_sp2:
        case N_ar:
        case N_pl3:
            return n_hydrogens > 0;
        case S_sp3:
            return n_hydrogens > 0;
        default:
            return false;
    }
}

inline bool classify_hbond_acceptor(uint8_t base_type, float partial_charge,
                                    int n_hydrogens) noexcept {
    switch (base_type) {
        case N_sp:
        case N_sp2:
        case N_sp3:
        case N_pl3:
            return partial_charge < 0.3f;
        case N_ar:
            // Heavy-atom-only aromatic N is treated as pyridine-like. Explicit
            // N-H evidence marks pyrrole-like donor and suppresses acceptance.
            return n_hydrogens == 0 && partial_charge < 0.3f;
        case N_am:
            return false;
        case O_sp2:
        case O_sp3:
        case O_co2:
        case O_ar:
            return true;
        case S_oxide:
        case S_dioxide:
            return true;
        case S_sp2:
        case S_sp3:
        case S_ar:
            return partial_charge < -0.3f;
        case HAL_F:
            return true;
        default:
            return false;
    }
}

// ─── topology-aware overloads ───────────────────────────────────────────────
// Identical to the 3-argument forms whenever `topo.known` is false, so call
// sites without connectivity keep their exact previous verdict.

inline bool classify_hbond_donor(uint8_t base_type, float partial_charge,
                                 int n_hydrogens,
                                 const HbondTopology& topo) noexcept {
    if (!topo.known)
        return classify_hbond_donor(base_type, partial_charge, n_hydrogens);

    switch (base_type) {
        case N_sp3: {
            // 4-coordinate → quaternary ammonium: no labile H, not a donor.
            if (topo.is_quaternary_nitrogen()) return false;
            // 1°/2° amine (≤2 heavy substituents) necessarily carries N–H.
            // 3° amine carries N–H only when protonated — evidenced by an
            // explicit H or a perceived positive formal charge.
            return topo.n_heavy_neighbors <= 2 || n_hydrogens > 0 ||
                   topo.is_cationic();
        }
        case N_quat:
            // SYBYL N.4 is by definition quaternary/protonated. It donates
            // through whatever H fill its fourth-plus valence.
            return topo.n_heavy_neighbors < 4 || n_hydrogens > 0;
        case O_sp3:
            // Hydroxyl (1 heavy substituent) donates; ether (2) does not.
            return topo.n_heavy_neighbors <= 1;
        default:
            return classify_hbond_donor(base_type, partial_charge, n_hydrogens);
    }
}

inline bool classify_hbond_acceptor(uint8_t base_type, float partial_charge,
                                    int n_hydrogens,
                                    const HbondTopology& topo) noexcept {
    if (!topo.known)
        return classify_hbond_acceptor(base_type, partial_charge, n_hydrogens);

    switch (base_type) {
        case N_sp3:
            // The lone pair is what accepts. Once the nitrogen is quaternised
            // or protonated it is spent — the atom is a donor, not an acceptor.
            return !topo.is_quaternary_nitrogen() && !topo.is_cationic();
        case N_quat:
            return false;
        default:
            return classify_hbond_acceptor(base_type, partial_charge,
                                           n_hydrogens);
    }
}

inline bool is_hbond_capable(uint8_t base_type, float partial_charge,
                              int n_hydrogens) noexcept {
    return classify_hbond_donor(base_type, partial_charge, n_hydrogens) ||
           classify_hbond_acceptor(base_type, partial_charge, n_hydrogens);
}

inline HbondRoles infer_hbond_roles(uint8_t base_type, float partial_charge,
                                    int n_hydrogens) noexcept {
    return {
        classify_hbond_donor(base_type, partial_charge, n_hydrogens),
        classify_hbond_acceptor(base_type, partial_charge, n_hydrogens),
    };
}

inline HbondRoles infer_hbond_roles(uint8_t base_type, float partial_charge,
                                    int n_hydrogens,
                                    const HbondTopology& topo) noexcept {
    return {
        classify_hbond_donor(base_type, partial_charge, n_hydrogens, topo),
        classify_hbond_acceptor(base_type, partial_charge, n_hydrogens, topo),
    };
}

// ─── SYBYL (1–40) ↔ base type (0–63) mapping ───────────────────────────────

// Forward mapping: SYBYL type → canonical base type (without context refinement)
// Every SYBYL type maps to a distinct base type — no Solvent fallback.
inline uint8_t sybyl_to_base(int sybyl_type) noexcept {
    // SYBYL types are 1-indexed (1–40)
    static constexpr uint8_t table[41] = {
        Solvent,     // 0: unused (placeholder)
        C_sp,        // 1: C.1
        C_sp2,       // 2: C.2
        C_sp3,       // 3: C.3
        C_ar,        // 4: C.AR
        C_cat,       // 5: C.CAT
        N_sp,        // 6: N.1
        N_sp2,       // 7: N.2
        N_sp3,       // 8: N.3
        N_quat,      // 9: N.4
        N_ar,        // 10: N.AR
        N_am,        // 11: N.AM
        N_pl3,       // 12: N.PL3
        O_sp2,       // 13: O.2
        O_sp3,       // 14: O.3
        O_co2,       // 15: O.CO2
        O_ar,        // 16: O.AR
        S_sp2,       // 17: S.2
        S_sp3,       // 18: S.3
        S_oxide,     // 19: S.O
        S_dioxide,   // 20: S.O2
        S_ar,        // 21: S.AR
        P_sp3,       // 22: P.3
        HAL_F,       // 23: F
        HAL_Cl,      // 24: CL
        HAL_Br,      // 25: BR
        HAL_I,       // 26: I
        HAL_Se,      // 27: SE
        Metal_Mg,    // 28: MG
        Metal_Sr,    // 29: SR
        Metal_Cu,    // 30: CU
        Metal_Mn,    // 31: MN
        Metal_Hg,    // 32: HG
        Metal_Cd,    // 33: CD
        Metal_Ni,    // 34: NI
        Metal_Zn,    // 35: ZN
        Metal_Ca,    // 36: CA
        Metal_Fe,    // 37: FE
        Metal_Co,    // 38: CO.OH
        Dummy,       // 39: DUMMY
        Solvent,     // 40: SOLVENT
    };
    if (sybyl_type < 0 || sybyl_type > 40) return Dummy;
    return table[sybyl_type];
}

// Reverse mapping: base type → SYBYL parent (1-indexed)
inline int base_to_sybyl_parent(uint8_t base_type) noexcept {
    static constexpr int table[BASE_TYPE_COUNT] = {
         1,  2,  3,  4,  5,        // 0–4:   C types
         6,  7,  8,  9, 10, 11, 12, // 5–11:  N types
        13, 14, 15, 16,            // 12–15: O types
        17, 18, 19, 20, 21,        // 16–20: S types
        22,                        // 21:    P.3
        23, 24, 25, 26,            // 22–25: halogens
         4,                        // 26:    C_ar_hetadj → C.AR
         2,                        // 27:    C_pi_bridge → C.2
        35, 36, 37,                // 28–30: Zn, Ca, Fe
        40,                        // 31:    Solvent
        // ── extended types (32–63) ──
        27,                        // 32:    HAL_Se → SE
        28,                        // 33:    Metal_Mg → MG
        29,                        // 34:    Metal_Sr → SR
        30,                        // 35:    Metal_Cu → CU
        31,                        // 36:    Metal_Mn → MN
        32,                        // 37:    Metal_Hg → HG
        33,                        // 38:    Metal_Cd → CD
        34,                        // 39:    Metal_Ni → NI
        38,                        // 40:    Metal_Co → CO.OH
        39,                        // 41:    Dummy → DUMMY
        // 42–63: reserved (map to DUMMY)
        39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
        39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
    };
    if (base_type >= BASE_TYPE_COUNT) return 39;
    return table[base_type];
}

// ─── context-aware refinement ───────────────────────────────────────────────
// Promotes C_ar to C_ar_hetadj or C_pi_bridge based on bonding environment.
// Call after initial sybyl_to_base() assignment when neighbor information is
// available.
//
// Parameters:
//   base:                initial base type from sybyl_to_base()
//   is_aromatic_carbon:  true if base == C_ar
//   has_heteroatom_neighbor: true if any bonded atom is N, O, or S
//   is_bridgehead:       true if atom is at ring junction (shared between two
//                        fused aromatic rings, e.g., indole C3a/C7a)

inline uint8_t refine_base_type(uint8_t base, bool is_aromatic_carbon,
                                 bool has_heteroatom_neighbor,
                                 bool is_bridgehead) noexcept {
    if (!is_aromatic_carbon || base != C_ar) return base;
    if (is_bridgehead) return C_pi_bridge;
    if (has_heteroatom_neighbor) return C_ar_hetadj;
    return base;
}

// ─── full encoding from SYBYL type + charge + H-bond info ──────────────────

inline uint8_t encode_from_sybyl(int sybyl_type, float partial_charge,
                                  int n_hydrogens,
                                  bool has_heteroatom_neighbor = false,
                                  bool is_bridgehead = false) noexcept {
    uint8_t base = sybyl_to_base(sybyl_type);
    bool aromatic_c = (sybyl_type == 4);  // C.AR
    base = refine_base_type(base, aromatic_c, has_heteroatom_neighbor,
                            is_bridgehead);
    const HbondRoles roles = infer_hbond_roles(base, partial_charge, n_hydrogens);
    return encode_roles(base, roles.donor, roles.acceptor);
}

// Topology-aware encoding. `topo.known == false` reproduces the 3-argument
// form exactly; SYBYL N.4 supplies its own formal charge, since that type is
// quaternary/protonated by definition.
inline uint8_t encode_from_sybyl(int sybyl_type, float partial_charge,
                                  int n_hydrogens,
                                  HbondTopology topo,
                                  bool has_heteroatom_neighbor = false,
                                  bool is_bridgehead = false) noexcept {
    uint8_t base = sybyl_to_base(sybyl_type);
    bool aromatic_c = (sybyl_type == 4);  // C.AR
    base = refine_base_type(base, aromatic_c, has_heteroatom_neighbor,
                            is_bridgehead);
    if (sybyl_type == 9) {  // N.4 — quaternary/protonated by SYBYL definition
        topo.formal_charge = 1;
        topo.charge_known  = true;
    } else if (topo.is_quaternary_nitrogen() && base == N_sp3) {
        topo.formal_charge = 1;
        topo.charge_known  = true;
    }
    const HbondRoles roles =
        infer_hbond_roles(base, partial_charge, n_hydrogens, topo);
    return encode_roles(base, roles.donor, roles.acceptor);
}

// ─── name table for debugging ───────────────────────────────────────────────

inline const char* base_type_name(uint8_t base) noexcept {
    static const char* names[BASE_TYPE_COUNT] = {
        "C.sp", "C.sp2", "C.sp3", "C.ar", "C.cat",
        "N.sp", "N.sp2", "N.sp3", "N.4", "N.ar", "N.am", "N.pl3",
        "O.sp2", "O.sp3", "O.co2", "O.ar",
        "S.sp2", "S.sp3", "S.O", "S.O2", "S.ar",
        "P.3",
        "F", "Cl", "Br", "I",
        "C.ar.het", "C.pi.br",
        "Zn", "Ca", "Fe",
        "SOL",
        // extended (32–41)
        "Se", "Mg", "Sr", "Cu", "Mn", "Hg", "Cd", "Ni", "Co.OH", "DUMMY",
        // reserved (42–63)
        "?42", "?43", "?44", "?45", "?46", "?47", "?48", "?49",
        "?50", "?51", "?52", "?53", "?54", "?55", "?56", "?57",
        "?58", "?59", "?60", "?61", "?62", "?63"
    };
    return (base < BASE_TYPE_COUNT) ? names[base] : "???";
}

} // namespace atom256
