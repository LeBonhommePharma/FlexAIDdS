// direct_input.cpp — Direct PDB/MOL2 input pipeline (no .inp files)
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

#include "direct_input.h"
#include "fileio.h"
#include "Mol2Reader.h"
#include "CleftDetector.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <random>

// ── external prototypes (declared in flexaid.h) ──
// read_pdb, modify_pdb, calc_center, residue_conect, assign_types,
// assign_radii_types, read_emat, detect_cleft, write_cleft_spheres,
// generate_grid, calc_cleftic, ic_bounds, add2_optimiz_vec, update_optres,
// rna_structure

int setup_direct_input(FA_Global* FA, GB_Global* GB, VC_Global* VC,
                       atom** atoms, resid** residue,
                       rot** rotamer, gridpoint** cleftgrid,
                       const char* receptor_pdb, const char* ligand_mol2)
{
    int k;
    char emat[MAX_PATH__];
    char deftyp[MAX_PATH__];
    char tmpprotname[MAX_PATH__];

    // ─── 1. Interaction matrix ────────────────────────────────────────
    if (FA->dependencies_path[0] != '\0') {
        strcpy(emat, FA->dependencies_path);
    } else {
        strcpy(emat, FA->base_path);
    }
#ifdef _WIN32
    strcat(emat, "\\MC_st0r5.2_6.dat");
#else
    strcat(emat, "/MC_st0r5.2_6.dat");
#endif
    printf("interaction matrix is <%s>\n", emat);
    read_emat(FA, emat);

    // ─── 2. Determine if RNA ──────────────────────────────────────────
    if (rna_structure(const_cast<char*>(receptor_pdb))) {
        printf("target molecule is a RNA structure\n");
        FA->is_protein = 0;
    }

    // ─── 3. Type definitions ──────────────────────────────────────────
    if (FA->dependencies_path[0] != '\0') {
        strcpy(deftyp, FA->dependencies_path);
    } else {
        strcpy(deftyp, FA->base_path);
    }
    if (FA->is_protein) {
#ifdef _WIN32
        strcat(deftyp, "\\AMINO.def");
#else
        strcat(deftyp, "/AMINO.def");
#endif
    } else {
#ifdef _WIN32
        strcat(deftyp, "\\NUCLEOTIDES.def");
#else
        strcat(deftyp, "/NUCLEOTIDES.def");
#endif
    }
    printf("definition of types is <%s>\n", deftyp);

    // ─── 4. Read receptor PDB ─────────────────────────────────────────
    // modify_pdb cleans and reorders the PDB; write to a temp file.
    int random_num = static_cast<int>(std::random_device{}() % 900000 + 100000);
    const std::string tmpdir = std::filesystem::temp_directory_path().string();
    snprintf(tmpprotname, MAX_PATH__, "%s/flexaid_tmp_%d.pdb", tmpdir.c_str(), random_num);

    printf("read PDB file <%s>\n", receptor_pdb);
    modify_pdb(const_cast<char*>(receptor_pdb), tmpprotname,
               FA->exclude_het, FA->remove_water, FA->is_protein,
               FA->keep_ions, FA->keep_structural_waters, FA->structural_water_bfactor_max,
               FA->binding_site_water_radius, FA->binding_site_water_hbond_required,
               FA->reflig_file);
    read_pdb(FA, atoms, residue, tmpprotname);
    remove(tmpprotname);

    // Count real atoms per residue
    (*residue)[FA->res_cnt].latm[0] = FA->atm_cnt - 1;
    for (k = 1; k <= FA->res_cnt; k++) {
        FA->atm_cnt_real += (*residue)[k].latm[0] - (*residue)[k].fatm[0] + 1;
    }

    calc_center(FA, *atoms, *residue);

    if (FA->is_protein) {
        residue_conect(FA, *atoms, *residue, deftyp);
    }
    assign_types(FA, *atoms, *residue, deftyp);

    // ─── 5. Read ligand MOL2 ──────────────────────────────────────────
    printf("read ligand MOL2 <%s>\n", ligand_mol2);
    if (!read_mol2_ligand(FA, atoms, residue, ligand_mol2)) {
        fprintf(stderr, "ERROR: failed to read ligand MOL2: %s\n", ligand_mol2);
        return 1;
    }

    // ─── 6. Assign radii ─────────────────────────────────────────────
    assign_radii_types(FA, *atoms, *residue);
    printf("radii are now assigned\n");

    // Recalculate global bounds to include ligand atoms
    calc_center(FA, *atoms, *residue);

    // ─── 7. Auto-detect binding cleft ────────────────────────────────
    printf("AUTO binding-site detection (CleftDetector) ...\n");

    sphere* spheres = detect_cleft(*atoms, *residue,
                                   FA->atm_cnt_real, FA->res_cnt);
    if (spheres == NULL) {
        fprintf(stderr, "ERROR: AUTO cleft detection found no cavities.\n");
        return 1;
    }

    // Write detected spheres for inspection
    char auto_sph[MAX_PATH__];
    if (FA->temp_path[0] != '\0') {
        strcpy(auto_sph, FA->temp_path);
    } else {
        strcpy(auto_sph, ".");
    }
#ifdef _WIN32
    strcat(auto_sph, "\\auto_cleft.sph");
#else
    strcat(auto_sph, "/auto_cleft.sph");
#endif
    write_cleft_spheres(spheres, auto_sph);
    printf("detected cleft written to %s\n", auto_sph);

    // ─── 8. Generate grid ────────────────────────────────────────────
    strcpy(FA->rngopt, "locclf");
    *cleftgrid = generate_grid(FA, spheres, *atoms, *residue);
    calc_cleftic(FA, *cleftgrid);

    // Free spheres linked-list
    while (spheres != NULL) {
        sphere* prev = spheres->prev;
        free(spheres);
        spheres = prev;
    }

    // ─── 9. IC bounds ────────────────────────────────────────────────
    ic_bounds(FA, FA->rngopt);

    // ─── 10. Set up optimization vectors for the ligand ──────────────
    // In direct mode, there are no OPTIMZ lines from .inp, so we
    // set up the ligand residue as the optimization target.
    // The ligand is always the last residue added.
    {
        int opt[2];
        char chain = ' ';

        // Add the ligand het residue (opt[0]=residue number, opt[1]=0 for ligand)
        opt[0] = FA->resligand->number;
        opt[1] = 0;
        add2_optimiz_vec(FA, *atoms, *residue, opt, chain, "");

        // Side-chain and normal-mode finalization
        add2_optimiz_vec(FA, *atoms, *residue, opt, chain, "SC");
        add2_optimiz_vec(FA, *atoms, *residue, opt, chain, "NM");
    }

    // ─── 11. Update optres pointers in atom structs ──────────────────
    update_optres(*atoms, *residue, FA->atm_cnt, FA->optres, FA->num_optres);

    printf("Direct input pipeline complete: %d atoms, %d residues, %d grid points\n",
           FA->atm_cnt, FA->res_cnt, FA->num_grd);

    return 0;
}
