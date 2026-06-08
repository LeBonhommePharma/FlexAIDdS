// SdfReader.cpp — Production-grade MDL V2000/V3000 MOL/SDF reader for FlexAIDdS
//
// Handles:
//   - V2000 standard fixed-width format
//   - V3000 extended format (space-delimited M  V30 records + line continuation)
//   - Windows (\r\n) and Unix (\n) line endings
//   - Mixed-case element symbols ("cl"→"Cl", "CL"→"Cl")
//   - M  CHG / M  RAD property blocks (formal charges)
//   - Multi-molecule SDF (reads first molecule only, stops at $$$$)
//   - Malformed counts lines with bond recovery
//   - Disconnected fragments (warns but continues)
//
// Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

#include "SdfReader.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cerrno>
#include <string>
#include <queue>
#include <vector>
#include <utility>

// ── helpers ──────────────────────────────────────────────────────────────────

static void strip_crlf(char* buf) {
    size_t len = strlen(buf);
    while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r'))
        buf[--len] = '\0';
}

// Capitalize first char, lowercase rest: "cl"→"Cl", "CL"→"Cl", "C"→"C"
static void normalize_element(char* elem) {
    if (!elem[0]) return;
    elem[0] = (char)toupper((unsigned char)elem[0]);
    for (int i = 1; elem[i]; ++i)
        elem[i] = (char)tolower((unsigned char)elem[i]);
}

/*
 * Element → canonical FlexAID VCT type number.
 *
 * These indices MUST match the canonical SYBYL type table in nrgrank_matrix.h
 * (kSybylTypes) and the receptor typing produced by assign_radii_types.cpp,
 * because atom.type is used directly as a row/column index into the VCT energy
 * matrix (MC_st0r5.2_6.dat / kEnergyMatrix). A wrong index scores a heteroatom
 * against the wrong element's energy row.
 *
 * SDF/MOL files carry no hybridization, so we fall back to the most common
 * generic perception for each element (sp3 for C/N/O/S). When richer typing is
 * available (MOL2 SYBYL strings), Mol2Reader::sybyl_to_flexaid_type is used.
 *
 * Canonical table (see nrgrank_matrix.h):
 *   1=C.1  2=C.2  3=C.3  4=C.AR  5=C.CAT
 *   6=N.1  7=N.2  8=N.3  9=N.4  10=N.AR 11=N.AM 12=N.PL3
 *  13=O.2 14=O.3 15=O.CO2 16=O.AR
 *  17=S.2 18=S.3 19=S.O  20=S.O2 21=S.AR
 *  22=P.3 23=F   24=CL   25=BR   26=I   27=SE
 *  28=MG 29=SR 30=CU 31=MN 32=HG 33=CD 34=NI 35=ZN 36=CA 37=FE 38=CO.OH
 *  39=DUMMY 40=SOLVENT
 */
static int element_to_flexaid_type(const char* elem) {
    if (!strcmp(elem, "C"))  return 3;   // C.3 (generic sp3 carbon)
    if (!strcmp(elem, "N"))  return 11;  // N.am (N.3/type-8 has zero matrix entries — N.am is the active generic N)
    if (!strcmp(elem, "O"))  return 14;  // O.3 (generic sp3 oxygen)
    if (!strcmp(elem, "S"))  return 18;  // S.3 (generic sp3 sulfur)
    if (!strcmp(elem, "P"))  return 22;  // P.3
    if (!strcmp(elem, "F"))  return 23;  // F
    if (!strcmp(elem, "Cl")) return 24;  // CL
    if (!strcmp(elem, "Br")) return 25;  // BR
    if (!strcmp(elem, "I"))  return 26;  // I
    if (!strcmp(elem, "Se")) return 27;  // SE
    if (!strcmp(elem, "Mg")) return 28;  // MG
    if (!strcmp(elem, "Sr")) return 29;  // SR
    if (!strcmp(elem, "Cu")) return 30;  // CU
    if (!strcmp(elem, "Mn")) return 31;  // MN
    if (!strcmp(elem, "Hg")) return 32;  // HG
    if (!strcmp(elem, "Cd")) return 33;  // CD
    if (!strcmp(elem, "Ni")) return 34;  // NI
    if (!strcmp(elem, "Zn")) return 35;  // ZN
    if (!strcmp(elem, "Ca")) return 36;  // CA
    if (!strcmp(elem, "Fe")) return 37;  // FE
    if (!strcmp(elem, "Co")) return 38;  // CO.OH
    return 39; // DUMMY (H and anything unknown — H is not scored)
}

static float element_radius(const char* elem) {
    if (!strcmp(elem, "C"))  return 1.70f;
    if (!strcmp(elem, "N"))  return 1.55f;
    if (!strcmp(elem, "O"))  return 1.52f;
    if (!strcmp(elem, "S"))  return 1.80f;
    if (!strcmp(elem, "P"))  return 1.80f;
    if (!strcmp(elem, "F"))  return 1.47f;
    if (!strcmp(elem, "Cl")) return 1.75f;
    if (!strcmp(elem, "Br")) return 1.85f;
    if (!strcmp(elem, "I"))  return 1.98f;
    if (!strcmp(elem, "H"))  return 1.20f;
    if (!strcmp(elem, "Fe")) return 1.34f;
    if (!strcmp(elem, "Zn")) return 1.22f;
    if (!strcmp(elem, "Ca")) return 2.31f;
    if (!strcmp(elem, "Mg")) return 1.73f;
    return 1.70f;
}

// ── temporary storage ────────────────────────────────────────────────────────

struct SdfAtom { float x, y, z; char elem[4]; float charge; };
struct SdfBond { int a1, a2, type; };

// ── V3000 parser ─────────────────────────────────────────────────────────────
//
// V3000 format uses extended connection-table blocks introduced with MDL Extended
// (ISIS/Base 2.5+).  Key grammar:
//   M  V30 BEGIN CTAB
//   M  V30 COUNTS na nb nSg n3D chiral
//   M  V30 BEGIN ATOM
//   M  V30 index element x y z mapping [CHG=n VAL=n ...]
//   M  V30 END ATOM
//   M  V30 BEGIN BOND
//   M  V30 index type atom1 atom2 [options]
//   M  V30 END BOND
//   M  V30 END CTAB
//   M  END
//
// Lines ending with '-' continue on the next line (V3000 §4.1.7).

static bool v30_read_logical_line(FILE* fp, char* buf, size_t bufsz) {
    buf[0] = '\0';
    char tmp[2048];
    while (fgets(tmp, (int)sizeof(tmp), fp)) {
        strip_crlf(tmp);
        size_t tlen = strlen(tmp);
        size_t curlen = strlen(buf);
        if (tlen > 0 && tmp[tlen-1] == '-') {
            // continuation: strip the '-' and accumulate
            tmp[tlen-1] = '\0';
            size_t space = bufsz - curlen - 1;
            strncat(buf, tmp, space);
        } else {
            size_t space = bufsz - curlen - 1;
            strncat(buf, tmp, space);
            return true;
        }
    }
    return buf[0] != '\0';
}

static int parse_sdf_v3000(FILE* fp,
                            std::vector<SdfAtom>& satoms,
                            std::vector<SdfBond>& sbonds)
{
    char buf[4096];
    bool in_atom = false, in_bond = false;

    while (v30_read_logical_line(fp, buf, sizeof(buf))) {
        if (strncmp(buf, "M  END", 6) == 0 ||
            strncmp(buf, "$$$$",   4) == 0) break;
        if (strncmp(buf, "M  V30 ", 7) != 0) continue;

        const char* rest = buf + 7;

        if (strncmp(rest, "BEGIN ATOM", 10) == 0) { in_atom = true;  continue; }
        if (strncmp(rest, "END ATOM",   8)  == 0) { in_atom = false; continue; }
        if (strncmp(rest, "BEGIN BOND", 10) == 0) { in_bond = true;  continue; }
        if (strncmp(rest, "END BOND",   8)  == 0) { in_bond = false; continue; }

        if (in_atom) {
            // M  V30 index element x y z mapping [key=val ...]
            int idx = 0;
            char raw_elem[32] = {};
            float x = 0, y = 0, z = 0;
            int mapnum = 0;
            if (sscanf(rest, "%d %31s %f %f %f %d",
                       &idx, raw_elem, &x, &y, &z, &mapnum) < 5) continue;

            // Strip brackets around query atoms like "[CH2]" or "[#6]"
            char* ep = raw_elem;
            if (*ep == '[') ep++;
            size_t elen = strlen(ep);
            if (elen > 0 && ep[elen-1] == ']') ep[elen-1] = '\0';

            SdfAtom a;
            a.x = x; a.y = y; a.z = z;
            strncpy(a.elem, ep, 3);
            a.elem[3] = '\0';
            normalize_element(a.elem);

            // Parse optional CHG= inline option
            a.charge = 0.0f;
            const char* chg = strstr(rest, "CHG=");
            if (chg) a.charge = (float)atof(chg + 4);

            satoms.push_back(a);
        }
        else if (in_bond) {
            // M  V30 index type atom1 atom2 [options]
            int idx = 0, btype = 0, a1 = 0, a2 = 0;
            if (sscanf(rest, "%d %d %d %d", &idx, &btype, &a1, &a2) == 4 &&
                a1 >= 1 && a2 >= 1 && btype >= 1) {
                sbonds.push_back({a1, a2, btype});
            }
        }
    }
    return (int)satoms.size();
}

// ── V2000 parser ─────────────────────────────────────────────────────────────

static int parse_sdf_v2000(FILE* fp, int natoms, int nbonds,
                            const char* sdf_file,
                            std::vector<SdfAtom>& satoms,
                            std::vector<SdfBond>& sbonds)
{
    char buf[512];
    satoms.reserve(natoms);
    sbonds.reserve(nbonds > 0 ? nbonds : 16);

    // ── atom block ──────────────────────────────────────────────────────────
    for (int i = 0; i < natoms; ++i) {
        if (!fgets(buf, sizeof(buf), fp)) {
            fprintf(stderr, "ERROR [SDF %s]: premature EOF in atom block (atom %d)\n",
                    sdf_file, i+1);
            return 0;
        }
        strip_crlf(buf);
        size_t blen = strlen(buf);

        SdfAtom a;
        a.charge = 0.0f;
        sscanf(buf, "%f %f %f", &a.x, &a.y, &a.z);

        // Element at columns 31-33 (0-indexed fixed-width)
        char sym[4] = {};
        int si = 0;
        for (int c = 31; c < 34 && c < (int)blen; ++c) {
            char ch = buf[c];
            if (ch != ' ' && ch != '\0' && ch != '\r' && ch != '\n')
                sym[si++] = ch;
        }
        sym[si] = '\0';
        if (si == 0) strncpy(sym, "Du", 3);  // empty element → dummy

        strncpy(a.elem, sym, 3);
        a.elem[3] = '\0';
        normalize_element(a.elem);

        satoms.push_back(a);
    }

    // ── bond block ──────────────────────────────────────────────────────────
    bool bonds_scanned_from_recovery = false;
    if (nbonds > 0) {
        for (int i = 0; i < nbonds; ++i) {
            if (!fgets(buf, sizeof(buf), fp)) {
                fprintf(stderr, "ERROR [SDF %s]: premature EOF in bond block (bond %d)\n",
                        sdf_file, i+1);
                return 0;
            }
            strip_crlf(buf);
            char f1[4], f2[4], f3[4];
            strncpy(f1, buf,   3); f1[3] = '\0';
            strncpy(f2, buf+3, 3); f2[3] = '\0';
            strncpy(f3, buf+6, 3); f3[3] = '\0';
            int a1 = atoi(f1), a2 = atoi(f2), btype = atoi(f3);
            if (a1 >= 1 && a1 <= natoms && a2 >= 1 && a2 <= natoms && btype >= 1)
                sbonds.push_back({a1, a2, btype});
        }
    } else {
        // Header claimed 0 bonds — scan remaining lines for bond records
        while (fgets(buf, sizeof(buf), fp)) {
            strip_crlf(buf);
            if (strncmp(buf, "M  END", 6) == 0 ||
                strncmp(buf, "$$$$",   4) == 0) {
                bonds_scanned_from_recovery = true;
                break;
            }
            if (strlen(buf) < 9) continue;
            bool has_decimal = false;
            for (int c = 0; c < 9; ++c)
                if (buf[c] == '.') { has_decimal = true; break; }
            if (has_decimal) continue;
            char f1[4], f2[4], f3[4];
            strncpy(f1, buf,   3); f1[3] = '\0';
            strncpy(f2, buf+3, 3); f2[3] = '\0';
            strncpy(f3, buf+6, 3); f3[3] = '\0';
            int a1 = atoi(f1), a2 = atoi(f2), btype = atoi(f3);
            if (a1 >= 1 && a1 <= natoms && a2 >= 1 && a2 <= natoms &&
                btype >= 1 && btype <= 8)
                sbonds.push_back({a1, a2, btype});
        }
        if (!sbonds.empty())
            fprintf(stderr, "WARN [SDF %s]: counts line claimed 0 bonds but %d recovered.\n",
                    sdf_file, (int)sbonds.size());
        // Already consumed through M  END / $$$$ — no further properties to parse
        if (bonds_scanned_from_recovery)
            return 1;
        return 1;
    }

    // ── properties block: M  CHG, M  RAD, M  END ────────────────────────────
    while (fgets(buf, sizeof(buf), fp)) {
        strip_crlf(buf);
        if (strncmp(buf, "M  END", 6) == 0 ||
            strncmp(buf, "$$$$",   4) == 0) break;

        // M  CHG  n  a1 v1 [a2 v2 ...]  — formal charges
        if (strncmp(buf, "M  CHG", 6) == 0) {
            int n = 0;
            const char* p = buf + 6;
            while (*p && !isdigit((unsigned char)*p)) ++p;
            n = (int)strtol(p, (char**)&p, 10);
            for (int k = 0; k < n; ++k) {
                while (*p && isspace((unsigned char)*p)) ++p;
                int aidx = (int)strtol(p, (char**)&p, 10);
                while (*p && isspace((unsigned char)*p)) ++p;
                int charge = (int)strtol(p, (char**)&p, 10);
                if (aidx >= 1 && aidx <= (int)satoms.size())
                    satoms[aidx-1].charge = (float)charge;
            }
        }
    }
    return 1;
}

// ── public API ───────────────────────────────────────────────────────────────

int read_sdf_ligand(FA_Global* FA, atom** atoms, resid** residue,
                    const char* sdf_file)
{
    FILE* fp = fopen(sdf_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR [SDF]: cannot open file: %s (%s)\n",
                sdf_file, strerror(errno));
        return 0;
    }

    printf("read_sdf_ligand: reading <%s>\n", sdf_file);

    char buf[512];
    char mol_name[64] = "LIG";

    // Line 1: molecule name
    if (fgets(buf, sizeof(buf), fp)) {
        strip_crlf(buf);
        const char* nm = buf;
        while (*nm == ' ' || *nm == '\t') ++nm;  // skip leading whitespace
        if (*nm) sscanf(nm, "%63s", mol_name);
    }
    // Line 2: program/timestamp (skip)
    if (!fgets(buf, sizeof(buf), fp)) {
        fprintf(stderr, "ERROR [SDF %s]: truncated at line 2\n", sdf_file);
        fclose(fp); return 0;
    }
    // Line 3: comment (skip)
    if (!fgets(buf, sizeof(buf), fp)) {
        fprintf(stderr, "ERROR [SDF %s]: truncated at line 3\n", sdf_file);
        fclose(fp); return 0;
    }

    // Line 4: counts line — determines V2000 vs V3000
    if (!fgets(buf, sizeof(buf), fp)) {
        fprintf(stderr, "ERROR [SDF %s]: missing counts line\n", sdf_file);
        fclose(fp); return 0;
    }
    strip_crlf(buf);

    const bool is_v3000 = (strstr(buf, "V3000") != nullptr);

    int natoms = 0, nbonds = 0;
    if (!is_v3000) {
        // V2000 fixed-width: cols 0-2 = natoms, cols 3-5 = nbonds
        char tmp[4];
        strncpy(tmp, buf,   3); tmp[3] = '\0'; natoms = atoi(tmp);
        strncpy(tmp, buf+3, 3); tmp[3] = '\0'; nbonds = atoi(tmp);

        if (natoms <= 0 || natoms > 9999) {
            fprintf(stderr, "ERROR [SDF %s]: invalid atom count %d in counts line\n",
                    sdf_file, natoms);
            fclose(fp); return 0;
        }
        if (nbonds < 0 || nbonds > 9999) {
            fprintf(stderr, "ERROR [SDF %s]: invalid bond count %d in counts line\n",
                    sdf_file, nbonds);
            fclose(fp); return 0;
        }
    }

    std::vector<SdfAtom> satoms;
    std::vector<SdfBond> sbonds;

    int ok;
    if (is_v3000)
        ok = parse_sdf_v3000(fp, satoms, sbonds);
    else
        ok = parse_sdf_v2000(fp, natoms, nbonds, sdf_file, satoms, sbonds);

    fclose(fp);

    if (!ok || satoms.empty()) {
        fprintf(stderr, "ERROR [SDF %s]: no atoms parsed\n", sdf_file);
        return 0;
    }

    natoms = (int)satoms.size();
    int nbonds_actual = (int)sbonds.size();
    printf("read_sdf_ligand: %d atoms, %d bonds (mol=%s)\n",
           natoms, nbonds_actual, mol_name);

    // ── Perceive hybridization from bond topology (C.ar / N.ar / O.co2) ───────
    // SDF/MOL files DO carry bond orders (MDL bond type 4 = aromatic). Pure
    // element typing never emits C.ar/N.ar/O.co2 — the only three VCT types that
    // actually differ from the degenerate sp2/sp3 rows in MC_st0r5.2_6.dat
    // (C.2↔C.3 and O.2↔O.3 are byte-identical). Derive these three from the
    // connection table so SDF scoring matches the obabel-converted MOL2 typing.
    std::vector<bool> is_aromatic(natoms, false);
    std::vector<bool> is_carboxylate_O(natoms, false);

    // Per-atom neighbour lists with bond order, 0-based (sbonds are 1-based).
    std::vector<std::vector<std::pair<int,int>>> nbr(natoms); // (neighbour, btype)
    for (const auto& sb : sbonds) {
        int i = sb.a1 - 1, j = sb.a2 - 1;
        if (i < 0 || i >= natoms || j < 0 || j >= natoms) continue;
        nbr[i].push_back({j, sb.type});
        nbr[j].push_back({i, sb.type});
        if (sb.type == 4) { is_aromatic[i] = true; is_aromatic[j] = true; }
    }

    // Carboxylate oxygens: O bonded to a C that has exactly two O neighbours,
    // one via a single bond and one via a double bond (COO⁻ / COOH pattern).
    for (int i = 0; i < natoms; ++i) {
        if (strcmp(satoms[i].elem, "O") != 0) continue;
        for (const auto& [cj, bt] : nbr[i]) {
            (void)bt;
            if (strcmp(satoms[cj].elem, "C") != 0) continue;
            int n_oxy = 0, n_single = 0, n_double = 0;
            for (const auto& [ok, obt] : nbr[cj]) {
                if (strcmp(satoms[ok].elem, "O") != 0) continue;
                ++n_oxy;
                if      (obt == 1) ++n_single;
                else if (obt == 2) ++n_double;
            }
            if (n_oxy == 2 && n_single == 1 && n_double == 1) {
                is_carboxylate_O[i] = true;
                break;
            }
        }
    }

    // ── Populate FA structures (same pattern as read_lig / Mol2Reader) ────────

    FA->optres = (OptRes*)malloc(FA->MIN_OPTRES * sizeof(OptRes));
    if (!FA->optres) { fprintf(stderr, "ERROR [SDF]: optres alloc failed\n"); return 0; }
    FA->MIN_OPTRES++;

    FA->num_het = 0;
    FA->num_het_atm = 0;

    FA->res_cnt++;
    if (FA->res_cnt >= FA->MIN_NUM_RESIDUE) {
        FA->MIN_NUM_RESIDUE = FA->res_cnt + 1;
        *residue = (resid*)realloc(*residue, FA->MIN_NUM_RESIDUE * sizeof(resid));
        if (!*residue) { fprintf(stderr, "ERROR [SDF]: residue realloc failed\n"); return 0; }
    }
    memset(&(*residue)[FA->res_cnt], 0, sizeof(resid));

    (*residue)[FA->res_cnt].fatm = (int*)malloc(sizeof(int));
    (*residue)[FA->res_cnt].latm = (int*)malloc(sizeof(int));
    (*residue)[FA->res_cnt].bond = (int*)malloc(FA->MIN_FLEX_BONDS * sizeof(int));
    if (!(*residue)[FA->res_cnt].fatm ||
        !(*residue)[FA->res_cnt].latm ||
        !(*residue)[FA->res_cnt].bond) {
        fprintf(stderr, "ERROR [SDF]: residue member alloc failed\n"); return 0;
    }
    memset((*residue)[FA->res_cnt].bond, 0, FA->MIN_FLEX_BONDS * sizeof(int));

    FA->num_het++;
    FA->het_res[FA->num_het] = FA->res_cnt;
    (*residue)[FA->res_cnt].bonded    = NULL;
    (*residue)[FA->res_cnt].shortpath = NULL;
    (*residue)[FA->res_cnt].shortflex = NULL;
    FA->resligand = &(*residue)[FA->res_cnt];
    (*residue)[FA->res_cnt].type   = 1;
    strncpy((*residue)[FA->res_cnt].name, mol_name, 3);
    (*residue)[FA->res_cnt].name[3] = '\0';
    (*residue)[FA->res_cnt].chn    = ' ';
    (*residue)[FA->res_cnt].number = 1;
    (*residue)[FA->res_cnt].rot    = 0;
    (*residue)[FA->res_cnt].fdih   = 0;

    // Map SDF 1-based atom index → internal FA index
    std::vector<int> idx_map(natoms + 1, 0);

    for (int ai = 0; ai < natoms; ++ai) {
        if (FA->atm_cnt + 1 >= FA->MIN_NUM_ATOM) {
            FA->MIN_NUM_ATOM += 50;
            *atoms = (atom*)realloc(*atoms, FA->MIN_NUM_ATOM * sizeof(atom));
            if (!*atoms) { fprintf(stderr, "ERROR [SDF]: atom realloc failed\n"); return 0; }
            memset(&(*atoms)[FA->MIN_NUM_ATOM - 50], 0, 50 * sizeof(atom));
        }

        // FlexAID uses 1-based atom indexing: atoms[0] is unused, atoms[1] is the
        // first real atom. Increment atm_cnt BEFORE using it as the storage index.
        FA->atm_cnt++;
        FA->atm_cnt_real++;
        FA->num_het_atm++;
        atom& a = (*atoms)[FA->atm_cnt];
        memset(&a, 0, sizeof(atom));

        int pdb_num = 90001 + ai;
        a.number = pdb_num;
        FA->num_atm[pdb_num] = FA->atm_cnt;
        idx_map[ai + 1] = FA->atm_cnt;

        a.coor[0] = a.coor_ori[0] = satoms[ai].x;
        a.coor[1] = a.coor_ori[1] = satoms[ai].y;
        a.coor[2] = a.coor_ori[2] = satoms[ai].z;
        a.coor_ref = NULL;

        snprintf(a.name, 5, "%-2s%d", satoms[ai].elem, (ai % 100));
        strncpy(a.element, satoms[ai].elem, 2);
        a.element[2] = '\0';

        // Element-only base type, then override with topology-perceived
        // hybridization for the three VCT types that matter (C.ar/N.ar/O.co2).
        int vct = element_to_flexaid_type(satoms[ai].elem);
        const char* el = satoms[ai].elem;
        const char* perceived = el;  // for diagnostics
        if      (!strcmp(el, "C") && is_aromatic[ai])       { vct = 4;  perceived = "C.ar"; }
        else if (!strcmp(el, "N") && is_aromatic[ai])       { vct = 10; perceived = "N.ar"; }
        else if (!strcmp(el, "O") && is_carboxylate_O[ai])  { vct = 15; perceived = "O.co2"; }

        a.type   = vct;
        a.radius = element_radius(satoms[ai].elem);
        a.charge = satoms[ai].charge;

        // [ATOM_TYPE] diagnostic — dump element→VCT integer mapping for each
        // ligand atom. SDF carries bond orders, so hybridization for the three
        // discriminating types (C.ar/N.ar/O.co2) is perceived from topology.
        if (getenv("FLEXAIDDS_DEBUG_TYPES")) {
            fprintf(stderr, "[ATOM_TYPE] idx=%d name=%s sybyl=%s vct=%d (SDF)\n",
                    ai, a.name, perceived, a.type);
        }
        a.ofres  = FA->res_cnt;
        a.recs   = 'f';
        a.bond[0] = 0;
        a.par    = NULL;
        a.cons   = NULL;
        a.optres = NULL;
        a.eigen  = NULL;

        if (ai == 0) (*residue)[FA->res_cnt].fatm[0] = FA->atm_cnt;
        (*residue)[FA->res_cnt].latm[0] = FA->atm_cnt;
    }

    // Populate bond arrays
    for (const auto& sb : sbonds) {
        if (sb.a1 < 1 || sb.a1 > natoms || sb.a2 < 1 || sb.a2 > natoms) continue;
        int fa1 = idx_map[sb.a1];
        int fa2 = idx_map[sb.a2];
        atom& ao = (*atoms)[fa1];
        atom& at = (*atoms)[fa2];
        if (ao.bond[0] < 6) { ao.bond[0]++; ao.bond[ao.bond[0]] = fa2; }
        if (at.bond[0] < 6) { at.bond[0]++; at.bond[at.bond[0]] = fa1; }
    }

    // ── Build IC reconstruction tree via BFS ──────────────────────────────────
    // Sets recs='m' and rec[0,1,2] so buildic()/buildcc() can reconstruct
    // 3D coordinates from internal coordinates during the GA.
    {
        int fa = (*residue)[FA->res_cnt].fatm[0];
        int la = (*residue)[FA->res_cnt].latm[0];
        int n  = la - fa + 1;

        std::vector<int>  parent(n, -1), grandpar(n, -1), grtgpar(n, -1);
        std::vector<bool> visited(n, false);
        std::queue<int>   q;
        q.push(fa);
        visited[0] = true;

        while (!q.empty()) {
            int cur = q.front(); q.pop();
            int ci  = cur - fa;
            for (int k = 1; k <= (*atoms)[cur].bond[0]; k++) {
                int nb = (*atoms)[cur].bond[k];
                if (nb < fa || nb > la) continue;
                int ni = nb - fa;
                if (!visited[ni]) {
                    visited[ni]  = true;
                    parent[ni]   = cur;
                    grandpar[ni] = parent[ci];
                    grtgpar[ni]  = grandpar[ci];
                    q.push(nb);
                }
            }
        }

        // Warn about disconnected atoms (missing bonds, salt forms, counterions)
        int unreached = 0;
        for (int i = 0; i < n; ++i) if (!visited[i]) ++unreached;
        if (unreached > 0)
            fprintf(stderr,
                "WARN [SDF %s]: %d atom(s) not reachable from atom 1 "
                "(disconnected fragment or missing bonds — check connectivity)\n",
                sdf_file, unreached);

        for (int ai = fa; ai <= la; ai++) {
            int li = ai - fa;
            atom& a  = (*atoms)[ai];
            a.recs   = 'm';
            a.rec[0] = (parent[li]   >= 0) ? parent[li]   : 0;
            a.rec[1] = (grandpar[li] >= 0) ? grandpar[li] : 0;
            a.rec[2] = (grtgpar[li]  >= 0) ? grtgpar[li]  : 0;
        }

        // Force GPA IC chain: first three atoms anchor rigid-body genes.
        if (n >= 2) {
            (*atoms)[fa+1].rec[0] = fa;
            (*atoms)[fa+1].rec[1] = 0;
            (*atoms)[fa+1].rec[2] = 0;
        }
        if (n >= 3) {
            (*atoms)[fa+2].rec[0] = fa+1;
            (*atoms)[fa+2].rec[1] = fa;
            (*atoms)[fa+2].rec[2] = 0;
        }

        buildic(FA, *atoms, *residue, FA->res_cnt);
    }

    // ── Build bonded matrix, shortest paths, shortflex ────────────────────────
    {
        int fa = (*residue)[FA->res_cnt].fatm[0];
        int la = (*residue)[FA->res_cnt].latm[0];
        int n  = la - fa + 1;
        int bondlist[MAX_ATM_HET];
        int neighbours[MAX_ATM_HET];
        int nbonded;
        for (int ai = fa; ai <= la; ai++) {
            nbonded = 0;
            bondedlist(*atoms, ai, FA->bloops, &nbonded, bondlist, neighbours);
            update_bonded(&(*residue)[FA->res_cnt], n, nbonded, bondlist, neighbours);
        }
        shortest_path(&(*residue)[FA->res_cnt], n, *atoms);
        assign_shortflex(&(*residue)[FA->res_cnt], n,
                         (*residue)[FA->res_cnt].fdih, *atoms);
    }

    FA->optres[0].rnum = FA->res_cnt;
    FA->optres[0].type = 1;
    FA->optres[0].tot  = FA->num_het_atm;
    FA->num_optres     = 1;

    printf("read_sdf_ligand: loaded %d atoms into FlexAID structures\n", natoms);
    return 1;
}
