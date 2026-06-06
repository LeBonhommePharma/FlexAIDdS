#include "Mol2Reader.h"
#include "flexaid.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>

/*
 * SYBYL atom type → FlexAID numeric type mapping.
 * FlexAID types are defined in the interaction matrix (MC_st0r5.2_6.dat).
 * This maps common SYBYL types to the closest FlexAID equivalents.
 */
static int sybyl_to_flexaid_type(const char* sybyl_type) {
    // Carbon types
    if (!strcmp(sybyl_type, "C.3"))   return 1;   // sp3 carbon
    if (!strcmp(sybyl_type, "C.2"))   return 2;   // sp2 carbon
    if (!strcmp(sybyl_type, "C.1"))   return 2;   // sp carbon (treat as sp2)
    if (!strcmp(sybyl_type, "C.ar"))  return 3;   // aromatic carbon
    if (!strcmp(sybyl_type, "C.cat")) return 2;   // carbocation

    // Nitrogen types
    if (!strcmp(sybyl_type, "N.3"))   return 4;   // sp3 nitrogen
    if (!strcmp(sybyl_type, "N.2"))   return 5;   // sp2 nitrogen
    if (!strcmp(sybyl_type, "N.1"))   return 5;   // sp nitrogen
    if (!strcmp(sybyl_type, "N.ar"))  return 6;   // aromatic nitrogen
    if (!strcmp(sybyl_type, "N.am"))  return 7;   // amide nitrogen
    if (!strcmp(sybyl_type, "N.pl3")) return 8;   // planar nitrogen
    if (!strcmp(sybyl_type, "N.4"))   return 9;   // quaternary nitrogen

    // Oxygen types
    if (!strcmp(sybyl_type, "O.3"))   return 10;  // sp3 oxygen
    if (!strcmp(sybyl_type, "O.2"))   return 11;  // sp2 oxygen
    if (!strcmp(sybyl_type, "O.co2")) return 12;  // carboxylate oxygen
    if (!strcmp(sybyl_type, "O.spc")) return 10;  // water oxygen
    if (!strcmp(sybyl_type, "O.t3p")) return 10;  // water oxygen

    // Sulfur types
    if (!strcmp(sybyl_type, "S.3"))   return 16;
    if (!strcmp(sybyl_type, "S.2"))   return 17;
    if (!strcmp(sybyl_type, "S.O"))   return 18;
    if (!strcmp(sybyl_type, "S.O2"))  return 19;

    // Phosphorus
    if (!strcmp(sybyl_type, "P.3"))   return 20;

    // Halogens
    if (!strcmp(sybyl_type, "F"))     return 13;
    if (!strcmp(sybyl_type, "Cl"))    return 14;
    if (!strcmp(sybyl_type, "Br"))    return 15;
    if (!strcmp(sybyl_type, "I"))     return 21;

    // Hydrogen
    if (!strcmp(sybyl_type, "H"))     return 22;
    if (!strcmp(sybyl_type, "H.spc")) return 22;
    if (!strcmp(sybyl_type, "H.t3p")) return 22;

    // Metals / misc
    if (!strcmp(sybyl_type, "Fe"))    return 30;
    if (!strcmp(sybyl_type, "Zn"))    return 31;
    if (!strcmp(sybyl_type, "Ca"))    return 32;
    if (!strcmp(sybyl_type, "Mg"))    return 33;

    return 39; // dummy type (unknown)
}

static float sybyl_radius(const char* sybyl_type) {
    // Multi-character element symbols must be checked before the switch
    if (!strcmp(sybyl_type, "Cl")) return 1.75f;
    if (!strcmp(sybyl_type, "Br")) return 1.85f;

    char elem = sybyl_type[0];
    switch (elem) {
        case 'C': return 1.70f;
        case 'N': return 1.55f;
        case 'O': return 1.52f;
        case 'S': return 1.80f;
        case 'P': return 1.80f;
        case 'F': return 1.47f;
        case 'H': return 1.20f;
        case 'I': return 1.98f;
        default:  return 1.70f;
    }
}

/* ── MOL2 section tags ─────────────────────────────────────────── */
enum Mol2Section { SEC_NONE, SEC_MOLECULE, SEC_ATOM, SEC_BOND, SEC_OTHER };

static Mol2Section parse_section_tag(const char* line) {
    if (strstr(line, "@<TRIPOS>MOLECULE")) return SEC_MOLECULE;
    if (strstr(line, "@<TRIPOS>ATOM"))     return SEC_ATOM;
    if (strstr(line, "@<TRIPOS>BOND"))     return SEC_BOND;
    if (line[0] == '@')                    return SEC_OTHER;
    return SEC_NONE;
}

/* ── main reader ───────────────────────────────────────────────── */

int read_mol2_ligand(FA_Global* FA, atom** atoms, resid** residue,
                     const char* mol2_file)
{
    FILE* fp = fopen(mol2_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: cannot open MOL2 file: %s\n", mol2_file);
        return 0;
    }

    printf("read_mol2_ligand: reading <%s>\n", mol2_file);

    /* Temporary storage */
    struct TmpAtom {
        int   id;
        char  name[8];
        float x, y, z;
        char  sybyl[16];
        int   res_id;
        char  res_name[8];
        float charge;
    };
    struct TmpBond {
        int origin, target, type;
    };

    std::vector<TmpAtom> tmp_atoms;
    std::vector<TmpBond> tmp_bonds;
    char mol_name[64] = "LIG";

    char buf[256];
    Mol2Section sec = SEC_NONE;
    int mol_line = 0; // line counter within MOLECULE section

    while (fgets(buf, sizeof(buf), fp)) {
        // Strip trailing newline
        size_t len = strlen(buf);
        while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r'))
            buf[--len] = '\0';

        // Check for section tag
        Mol2Section tag = parse_section_tag(buf);
        if (tag != SEC_NONE) { sec = tag; mol_line = 0; continue; }

        switch (sec) {
        case SEC_MOLECULE:
            mol_line++;
            if (mol_line == 1) {
                // molecule name
                sscanf(buf, "%63s", mol_name);
            }
            break;

        case SEC_ATOM: {
            TmpAtom a;
            memset(&a, 0, sizeof(a));
            // Fields: atom_id atom_name x y z atom_type [subst_id subst_name charge]
            int nf = sscanf(buf, "%d %7s %f %f %f %15s %d %7s %f",
                            &a.id, a.name, &a.x, &a.y, &a.z,
                            a.sybyl, &a.res_id, a.res_name, &a.charge);
            if (nf >= 6)
                tmp_atoms.push_back(a);
            break;
        }

        case SEC_BOND: {
            TmpBond b;
            char btype[8] = {};
            int nf = sscanf(buf, "%*d %d %d %7s", &b.origin, &b.target, btype);
            if (nf >= 2) {
                // "nc" = non-connected (Tripos placeholder), "du" = dummy bond — skip both
                if (!strcmp(btype, "nc") || !strcmp(btype, "du")) break;
                b.type = 1; // single bond default
                if (!strcmp(btype, "2") || !strcmp(btype, "do")) b.type = 2;
                if (!strcmp(btype, "3") || !strcmp(btype, "tr")) b.type = 3;
                if (!strcmp(btype, "ar"))                        b.type = 4;
                if (!strcmp(btype, "am"))                        b.type = 1; // amide→single
                tmp_bonds.push_back(b);
            }
            break;
        }

        default:
            break;
        }
    }
    fclose(fp);

    if (tmp_atoms.empty()) {
        fprintf(stderr, "ERROR: no atoms found in MOL2 file %s\n", mol2_file);
        return 0;
    }

    printf("read_mol2_ligand: %zu atoms, %zu bonds from %s\n",
           tmp_atoms.size(), tmp_bonds.size(), mol2_file);

    /* ── populate FA structures (mirrors read_lig logic) ────────── */

    // Set up optres
    FA->optres = (OptRes*)malloc(FA->MIN_OPTRES * sizeof(OptRes));
    if (!FA->optres) { fprintf(stderr, "ERROR: optres alloc\n"); return 0; }
    FA->MIN_OPTRES++;

    FA->num_het = 0;
    FA->num_het_atm = 0;

    // New residue for the ligand
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
    (*residue)[FA->res_cnt].type = 1; // ligand
    strncpy((*residue)[FA->res_cnt].name, mol_name, 3);
    (*residue)[FA->res_cnt].name[3] = '\0';
    (*residue)[FA->res_cnt].chn = ' ';
    (*residue)[FA->res_cnt].number = 1;
    (*residue)[FA->res_cnt].rot = 0;
    (*residue)[FA->res_cnt].fdih = 0;

    // Build mol2_id → internal index map
    std::map<int, int> id_map;
    [[maybe_unused]] int first_atm = FA->atm_cnt + 1;

    for (size_t ai = 0; ai < tmp_atoms.size(); ++ai) {
        // Grow atom array if needed
        if (FA->atm_cnt + 1 >= FA->MIN_NUM_ATOM) {
            FA->MIN_NUM_ATOM += 50;
            *atoms = (atom*)realloc(*atoms, FA->MIN_NUM_ATOM * sizeof(atom));
            if (!*atoms) { fprintf(stderr, "ERROR: atom realloc\n"); return 0; }
            memset(&(*atoms)[FA->MIN_NUM_ATOM - 50], 0, 50 * sizeof(atom));
        }

        // FlexAID uses 1-based atom indexing: atoms[0] is unused, atoms[1] is the
        // first real atom. Increment atm_cnt BEFORE using it as the storage index.
        FA->atm_cnt++;
        FA->atm_cnt_real++;
        FA->num_het_atm++;
        atom& a = (*atoms)[FA->atm_cnt];
        memset(&a, 0, sizeof(atom));

        // Assign a PDB-style number starting at 90001
        int pdb_num = 90001 + static_cast<int>(ai);
        a.number = pdb_num;
        FA->num_atm[pdb_num] = FA->atm_cnt;
        id_map[tmp_atoms[ai].id] = FA->atm_cnt;

        a.coor[0] = a.coor_ori[0] = tmp_atoms[ai].x;
        a.coor[1] = a.coor_ori[1] = tmp_atoms[ai].y;
        a.coor[2] = a.coor_ori[2] = tmp_atoms[ai].z;
        a.coor_ref = NULL;

        strncpy(a.name, tmp_atoms[ai].name, 4);
        a.name[4] = '\0';

        // Element from SYBYL type (first char(s) before '.')
        char elem[3] = {};
        elem[0] = tmp_atoms[ai].sybyl[0];
        if (islower(tmp_atoms[ai].sybyl[1]))
            elem[1] = tmp_atoms[ai].sybyl[1];
        strncpy(a.element, elem, 2);
        a.element[2] = '\0';

        a.type = sybyl_to_flexaid_type(tmp_atoms[ai].sybyl);
        a.radius = sybyl_radius(tmp_atoms[ai].sybyl);
        a.charge = tmp_atoms[ai].charge; // propagate MOL2 partial charge
        a.ofres = FA->res_cnt;
        a.recs = 'f'; // flexible by default
        a.bond[0] = 0;
        a.par = NULL;
        a.cons = NULL;
        a.optres = NULL;
        a.eigen = NULL;

        // Update residue first/last atom
        if (ai == 0) (*residue)[FA->res_cnt].fatm[0] = FA->atm_cnt;
        (*residue)[FA->res_cnt].latm[0] = FA->atm_cnt;
    }

    // Populate bond arrays from MOL2 bond table
    for (size_t bi = 0; bi < tmp_bonds.size(); ++bi) {
        auto it_o = id_map.find(tmp_bonds[bi].origin);
        auto it_t = id_map.find(tmp_bonds[bi].target);
        if (it_o == id_map.end() || it_t == id_map.end()) continue;

        int idx_o = it_o->second;
        int idx_t = it_t->second;

        atom& ao = (*atoms)[idx_o];
        atom& at = (*atoms)[idx_t];

        if (ao.bond[0] < 6) { ao.bond[0]++; ao.bond[ao.bond[0]] = idx_t; }
        if (at.bond[0] < 6) { at.bond[0]++; at.bond[at.bond[0]] = idx_o; }
    }

    // Build IC reconstruction tree via BFS from the first ligand atom.
    // Each atom needs rec[0,1,2] = internal indices of parent, grandparent,
    // great-grandparent in the spanning tree (0 = use FA->ori reference frame).
    // recs='m' marks atoms for IC-based Cartesian reconstruction in buildcc.
    {
        int fa = (*residue)[FA->res_cnt].fatm[0];
        int la = (*residue)[FA->res_cnt].latm[0];
        int n  = la - fa + 1;

        // BFS ancestors: parent_idx[i-fa] = internal index of parent in tree, -1 if root
        std::vector<int> parent(n, -1);
        std::vector<int> grandpar(n, -1);
        std::vector<int> grtgpar(n, -1);
        std::vector<bool> visited(n, false);

        std::queue<int> q;
        q.push(fa);
        visited[fa - fa] = true;

        while (!q.empty()) {
            int cur = q.front(); q.pop();
            int ci = cur - fa; // local index

            for (int k = 1; k <= (*atoms)[cur].bond[0]; k++) {
                int nb = (*atoms)[cur].bond[k]; // internal index of neighbour
                if (nb < fa || nb > la) continue;
                int ni = nb - fa;
                if (!visited[ni]) {
                    visited[ni] = true;
                    parent[ni]   = cur;
                    grandpar[ni] = parent[ci];   // may be -1
                    grtgpar[ni]  = grandpar[ci]; // may be -1
                    q.push(nb);
                }
            }
        }

        // Assign rec[] and recs for each ligand atom
        for (int ai = fa; ai <= la; ai++) {
            int li = ai - fa;
            atom& a = (*atoms)[ai];
            a.recs   = 'm';
            a.rec[0] = (parent[li]  >= 0) ? parent[li]  : 0;
            a.rec[1] = (grandpar[li] >= 0) ? grandpar[li] : 0;
            a.rec[2] = (grtgpar[li] >= 0) ? grtgpar[li] : 0;
        }

        // Force GPA IC chain: GPA1 parent=GPA0, GPA2 parent=GPA1.
        // The first three ligand atoms serve as GPA (global positioning) atoms
        // for rigid-body translation and rotation genes. ic2cf reconstructs
        // GPA0 from cleft-grid IC (relative to FA->ori), then reconstructs
        // GPA1/GPA2 from rotation genes relative to their parent GPA atom.
        // Without this override, BFS may assign ori-based rec[] to GPA1/GPA2
        // (e.g. when atom 0 and atom 1 are not directly bonded), which breaks
        // the GPA rotation reference frame when GPA0 translates.
        if (n >= 2) {
            (*atoms)[fa+1].rec[0] = fa;    // GPA1 parent  = GPA0
            (*atoms)[fa+1].rec[1] = 0;     // GPA1 grandpar = ori
            (*atoms)[fa+1].rec[2] = 0;
        }
        if (n >= 3) {
            (*atoms)[fa+2].rec[0] = fa+1;  // GPA2 parent   = GPA1
            (*atoms)[fa+2].rec[1] = fa;    // GPA2 grandpar = GPA0
            (*atoms)[fa+2].rec[2] = 0;
        }

        // Compute IC for all ligand atoms using buildic() which reads rec[]
        // and current coor[] to produce dis/ang/dih consistent with buildcc.
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

    printf("read_mol2_ligand: loaded %zu atoms into FlexAID structures\n",
           tmp_atoms.size());
    return 1;
}
