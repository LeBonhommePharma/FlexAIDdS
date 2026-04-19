#include "SdfReader.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>

/*
 * MDL V2000 MOL/SDF file layout:
 *   Line 1:  molecule name
 *   Line 2:  program/timestamp
 *   Line 3:  comment
 *   Line 4:  counts line  "aaabbblll..."
 *            aaa = #atoms, bbb = #bonds (each 3 chars, right-justified)
 *   Atom block (one line per atom):
 *     xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
 *     positions: x(0-9) y(10-19) z(20-29) symbol(31-33) ...
 *   Bond block (one line per bond):
 *     111222tttsssxxxrrrccc
 *     111=first atom, 222=second atom, ttt=bond type (each 3 chars)
 *   Properties block (M  END terminates)
 *   $$$$ separates molecules in multi-molecule SDF
 */

static int element_to_flexaid_type(const char* elem) {
    if (!strcmp(elem, "C"))  return 1;
    if (!strcmp(elem, "N"))  return 4;
    if (!strcmp(elem, "O"))  return 10;
    if (!strcmp(elem, "S"))  return 16;
    if (!strcmp(elem, "P"))  return 20;
    if (!strcmp(elem, "F"))  return 13;
    if (!strcmp(elem, "Cl")) return 14;
    if (!strcmp(elem, "Br")) return 15;
    if (!strcmp(elem, "I"))  return 21;
    if (!strcmp(elem, "H"))  return 22;
    if (!strcmp(elem, "Fe")) return 30;
    if (!strcmp(elem, "Zn")) return 31;
    if (!strcmp(elem, "Ca")) return 32;
    if (!strcmp(elem, "Mg")) return 33;
    return 39; // dummy
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
    return 1.70f;
}

int read_sdf_ligand(FA_Global* FA, atom** atoms, resid** residue,
                    const char* sdf_file)
{
    FILE* fp = fopen(sdf_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open SDF file: %s\n", sdf_file);
        return 0;
    }

    printf("read_sdf_ligand: reading <%s>\n", sdf_file);

    char buf[256];
    char mol_name[64] = "LIG";

    // Line 1: molecule name
    if (fgets(buf, sizeof(buf), fp)) {
        size_t len = strlen(buf);
        if (len > 0 && buf[len - 1] == '\n') buf[len - 1] = '\0';
        if (strlen(buf) > 0) sscanf(buf, "%63s", mol_name);
    }
    // Line 2: program/timestamp (skip)
    if (!fgets(buf, sizeof(buf), fp)) { fclose(fp); return 0; }
    // Line 3: comment (skip)
    if (!fgets(buf, sizeof(buf), fp)) { fclose(fp); return 0; }

    // Line 4: counts
    if (!fgets(buf, sizeof(buf), fp)) {
        fprintf(stderr, "ERROR: premature EOF in SDF counts line\n");
        fclose(fp); return 0;
    }

    int natoms = 0, nbonds = 0;
    // V2000 counts: first 3 chars = natoms, next 3 = nbonds
    char tmp[4];
    strncpy(tmp, buf, 3); tmp[3] = '\0'; natoms = atoi(tmp);
    strncpy(tmp, buf + 3, 3); tmp[3] = '\0'; nbonds = atoi(tmp);

    if (natoms <= 0 || natoms > 9999) {
        fprintf(stderr, "ERROR: invalid atom count %d in SDF file\n", natoms);
        fclose(fp); return 0;
    }
    if (nbonds < 0 || nbonds > 9999) {
        fprintf(stderr, "ERROR: invalid bond count %d in SDF file\n", nbonds);
        fclose(fp); return 0;
    }

    /* ── Temporary atom storage ────────────────────────────── */
    struct SdfAtom { float x, y, z; char elem[4]; };
    struct SdfBond { int a1, a2, type; };

    std::vector<SdfAtom> satoms(natoms);
    std::vector<SdfBond> sbonds(nbonds);

    // Read atom block
    for (int i = 0; i < natoms; ++i) {
        if (!fgets(buf, sizeof(buf), fp)) {
            fprintf(stderr, "ERROR: premature EOF in SDF atom block at atom %d\n", i + 1);
            fclose(fp); return 0;
        }
        // V2000: x(0-9), y(10-19), z(20-29), symbol(31-33)
        sscanf(buf, "%f %f %f", &satoms[i].x, &satoms[i].y, &satoms[i].z);

        // Element symbol at column 31 (0-indexed), up to 3 chars
        char sym[4] = {};
        int si = 0;
        for (int c = 31; c < 34 && c < (int)strlen(buf); ++c) {
            if (buf[c] != ' ' && buf[c] != '\0')
                sym[si++] = buf[c];
        }
        sym[si] = '\0';
        strncpy(satoms[i].elem, sym, 3);
        satoms[i].elem[3] = '\0';
    }

    // Read bond block
    for (int i = 0; i < nbonds; ++i) {
        if (!fgets(buf, sizeof(buf), fp)) {
            fprintf(stderr, "ERROR: premature EOF in SDF bond block at bond %d\n", i + 1);
            fclose(fp); return 0;
        }
        // V2000: atom1(0-2), atom2(3-5), type(6-8) — 3-char right-justified integers
        char f1[4], f2[4], f3[4];
        strncpy(f1, buf, 3); f1[3] = '\0';
        strncpy(f2, buf + 3, 3); f2[3] = '\0';
        strncpy(f3, buf + 6, 3); f3[3] = '\0';
        sbonds[i].a1 = atoi(f1);
        sbonds[i].a2 = atoi(f2);
        sbonds[i].type = atoi(f3);
    }

    fclose(fp);

    printf("read_sdf_ligand: %d atoms, %d bonds\n", natoms, nbonds);

    /* ── Populate FA structures (same pattern as read_lig / Mol2Reader) ── */

    FA->optres = (OptRes*)malloc(FA->MIN_OPTRES * sizeof(OptRes));
    if (!FA->optres) { fprintf(stderr, "ERROR: optres alloc\n"); return 0; }
    FA->MIN_OPTRES++;

    FA->num_het = 0;
    FA->num_het_atm = 0;

    // New residue
    FA->res_cnt++;
    if (FA->res_cnt >= FA->MIN_NUM_RESIDUE) {
        FA->MIN_NUM_RESIDUE = FA->res_cnt + 1;
        *residue = (resid*)realloc(*residue, FA->MIN_NUM_RESIDUE * sizeof(resid));
        if (!*residue) { fprintf(stderr, "ERROR: residue realloc\n"); return 0; }
    }
    memset(&(*residue)[FA->res_cnt], 0, sizeof(resid));

    (*residue)[FA->res_cnt].fatm = (int*)malloc(sizeof(int));
    (*residue)[FA->res_cnt].latm = (int*)malloc(sizeof(int));
    (*residue)[FA->res_cnt].bond = (int*)malloc(FA->MIN_FLEX_BONDS * sizeof(int));
    if (!(*residue)[FA->res_cnt].fatm ||
        !(*residue)[FA->res_cnt].latm ||
        !(*residue)[FA->res_cnt].bond) {
        fprintf(stderr, "ERROR: residue member alloc\n"); return 0;
    }
    memset((*residue)[FA->res_cnt].bond, 0, FA->MIN_FLEX_BONDS * sizeof(int));

    FA->num_het++;
    FA->het_res[FA->num_het] = FA->res_cnt;
    (*residue)[FA->res_cnt].bonded = NULL;
    (*residue)[FA->res_cnt].shortpath = NULL;
    (*residue)[FA->res_cnt].shortflex = NULL;
    FA->resligand = &(*residue)[FA->res_cnt];
    (*residue)[FA->res_cnt].type = 1;
    strncpy((*residue)[FA->res_cnt].name, mol_name, 3);
    (*residue)[FA->res_cnt].name[3] = '\0';
    (*residue)[FA->res_cnt].chn = ' ';
    (*residue)[FA->res_cnt].number = 1;
    (*residue)[FA->res_cnt].rot = 0;
    (*residue)[FA->res_cnt].fdih = 0;

    // Map SDF 1-based atom index → internal FA index
    std::vector<int> idx_map(natoms + 1, 0); // 1-based

    for (int ai = 0; ai < natoms; ++ai) {
        FA->atm_cnt++;
        FA->atm_cnt_real++;
        FA->num_het_atm++;

        if (FA->atm_cnt >= FA->MIN_NUM_ATOM) {
            FA->MIN_NUM_ATOM += 50;
            *atoms = (atom*)realloc(*atoms, FA->MIN_NUM_ATOM * sizeof(atom));
            if (!*atoms) { fprintf(stderr, "ERROR: atom realloc\n"); return 0; }
            memset(&(*atoms)[FA->MIN_NUM_ATOM - 50], 0, 50 * sizeof(atom));
        }

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

        // Build atom name from element + index
        snprintf(a.name, 5, "%-2s%d", satoms[ai].elem, (ai % 100));
        strncpy(a.element, satoms[ai].elem, 2);
        a.element[2] = '\0';

        a.type   = element_to_flexaid_type(satoms[ai].elem);
        a.radius = element_radius(satoms[ai].elem);
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

    // Populate bonds
    for (int bi = 0; bi < nbonds; ++bi) {
        int i1 = sbonds[bi].a1;
        int i2 = sbonds[bi].a2;
        if (i1 < 1 || i1 > natoms || i2 < 1 || i2 > natoms) continue;

        int fa1 = idx_map[i1];
        int fa2 = idx_map[i2];

        atom& ao = (*atoms)[fa1];
        atom& at = (*atoms)[fa2];

        if (ao.bond[0] < 6) { ao.bond[0]++; ao.bond[ao.bond[0]] = fa2; }
        if (at.bond[0] < 6) { at.bond[0]++; at.bond[at.bond[0]] = fa1; }
    }

    // Finalise optres for the ligand (mirrors read_lig.cpp logic)
    FA->optres[0].rnum = FA->res_cnt;
    FA->optres[0].type = 1;
    FA->optres[0].tot  = FA->num_het_atm;
    FA->num_optres     = 1;

    printf("read_sdf_ligand: loaded %d atoms into FlexAID structures\n", natoms);
    return 1;
}
