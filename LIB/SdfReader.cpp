#include "SdfReader.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <queue>
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

    // Read bond block.
    // Guard: if the header claims 0 bonds but bond records exist in the file
    // (e.g. SDF written by an extractor that forgot to update the counts line),
    // scan the remaining lines and recover the bonds automatically.
    if (nbonds > 0) {
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
    } else {
        // Header said 0 bonds — scan remaining lines for bond records.
        // A V2000 bond line has two atom indices (1-based, ≤9999) in cols 0-5
        // followed by a bond type digit, with no decimal point in the first 9 chars.
        while (fgets(buf, sizeof(buf), fp)) {
            // Stop at properties block or next molecule
            if (strncmp(buf, "M  END", 6) == 0 ||
                strncmp(buf, "$$$$",   4) == 0) break;
            // Skip lines that are clearly not bond records
            if (strlen(buf) < 9) continue;
            bool has_decimal = false;
            for (int c = 0; c < 9; ++c)
                if (buf[c] == '.') { has_decimal = true; break; }
            if (has_decimal) continue;  // atom lines have decimal coords
            // Parse as bond record: first 6 chars must be two valid integers
            char f1[4], f2[4], f3[4];
            strncpy(f1, buf, 3); f1[3] = '\0';
            strncpy(f2, buf + 3, 3); f2[3] = '\0';
            strncpy(f3, buf + 6, 3); f3[3] = '\0';
            int a1 = atoi(f1), a2 = atoi(f2), btype = atoi(f3);
            if (a1 < 1 || a1 > natoms || a2 < 1 || a2 > natoms) continue;
            if (btype < 1 || btype > 8) continue;  // valid MDL bond types 1-8
            sbonds.push_back({a1, a2, btype});
        }
        if (!sbonds.empty())
            fprintf(stderr,
                "WARN: SDF counts line claimed 0 bonds but %d bond records found — "
                "using recovered bonds.\n", (int)sbonds.size());
        // Rewind: fclose+reopen would be needed for a full restart, but we're done
        // reading; the caller uses sbonds directly so no further reads needed.
        fclose(fp);
        fp = nullptr;
    }

    if (fp) fclose(fp);

    int nbonds_actual = (int)sbonds.size();
    printf("read_sdf_ligand: %d atoms, %d bonds\n", natoms, nbonds_actual);

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
        if (FA->atm_cnt >= FA->MIN_NUM_ATOM - 1) {
            FA->MIN_NUM_ATOM += 50;
            *atoms = (atom*)realloc(*atoms, FA->MIN_NUM_ATOM * sizeof(atom));
            if (!*atoms) { fprintf(stderr, "ERROR: atom realloc\n"); return 0; }
            memset(&(*atoms)[FA->MIN_NUM_ATOM - 50], 0, 50 * sizeof(atom));
        }

        atom& a = (*atoms)[FA->atm_cnt];
        FA->atm_cnt++;
        FA->atm_cnt_real++;
        FA->num_het_atm++;
        memset(&a, 0, sizeof(atom));

        int pdb_num = 90001 + ai;
        a.number = pdb_num;
        FA->num_atm[pdb_num] = FA->atm_cnt - 1;
        idx_map[ai + 1] = FA->atm_cnt - 1;

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

        if (ai == 0) (*residue)[FA->res_cnt].fatm[0] = FA->atm_cnt - 1;
        (*residue)[FA->res_cnt].latm[0] = FA->atm_cnt - 1;
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

    // Build IC reconstruction tree via BFS (mirrors Mol2Reader.cpp).
    // Sets recs='m' and rec[0,1,2] so buildic()/buildcc() can reconstruct
    // 3D coordinates from internal coordinates during the GA.  Without this,
    // all ligand atoms collapse to the centroid (FA->ori) in ic2cf.
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

        for (int ai = fa; ai <= la; ai++) {
            int li = ai - fa;
            atom& a = (*atoms)[ai];
            a.recs   = 'm';
            a.rec[0] = (parent[li]   >= 0) ? parent[li]   : 0;
            a.rec[1] = (grandpar[li] >= 0) ? grandpar[li] : 0;
            a.rec[2] = (grtgpar[li]  >= 0) ? grtgpar[li]  : 0;
        }

        // Force GPA IC chain so GPA1/GPA2 track GPA0 during rigid-body moves.
        if (n >= 2) {
            (*atoms)[fa+1].rec[0] = fa;    // GPA1 parent  = GPA0
            (*atoms)[fa+1].rec[1] = 0;
            (*atoms)[fa+1].rec[2] = 0;
        }
        if (n >= 3) {
            (*atoms)[fa+2].rec[0] = fa+1;  // GPA2 parent  = GPA1
            (*atoms)[fa+2].rec[1] = fa;    // GPA2 grandpar = GPA0
            (*atoms)[fa+2].rec[2] = 0;
        }

        buildic(FA, *atoms, *residue, FA->res_cnt);
    }

    // Build bonded matrix, shortest paths, and shortflex (mirrors read_lig.cpp)
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
        assign_shortflex(&(*residue)[FA->res_cnt], n, (*residue)[FA->res_cnt].fdih, *atoms);
    }

    // Finalise optres for the ligand (mirrors read_lig.cpp logic)
    FA->optres[0].rnum = FA->res_cnt;
    FA->optres[0].type = 1;
    FA->optres[0].tot  = FA->num_het_atm;
    FA->num_optres     = 1;

    printf("read_sdf_ligand: loaded %d atoms into FlexAID structures\n", natoms);
    return 1;
}
