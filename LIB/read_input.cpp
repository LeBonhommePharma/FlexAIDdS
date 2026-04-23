#include "flexaid.h"
#include "fileio.h"
#include "CavityDetect/CavityDetect.h"
#include "CleftDetector.h"
#include "MIFGrid.h"
#include "RefLigSeed.h"
#include "CavityDetect/SpatialGrid.h"
#include "BindingResidues.h"
#include <vector>
#include <algorithm>
#include <random>

/*****************************************************************************
 * compute_mif_and_reflig — MIF computation, grid prioritization, RefLig seeding
 * Called after generate_grid() + calc_cleftic() for all detection modes.
 *****************************************************************************/
static void compute_mif_and_reflig(FA_Global* FA, atom* atoms,
                                    gridpoint** cleftgrid, const char* lig_file) {
	const bool need_mif = FA->mif_enabled || FA->grid_prio_percent < 100.0f;
	const bool need_reflig = FA->reflig_file[0] != '\0';
	const bool need_hetatm = !need_reflig && FA->reflig_hetatm_fallback && lig_file[0] != '\0';

	if (!need_mif && !need_reflig && !need_hetatm) return;

	// Build SpatialGrid for neighbor queries
	std::vector<atom> protein_atoms(atoms, atoms + FA->atm_cnt_real);
	cavity_detect::SpatialGrid sg;
	sg.build(protein_atoms);

	// ── MIF computation ──
	if (need_mif) {
		auto mif = mif::compute_mif(*cleftgrid, FA->num_grd,
		                             atoms, FA->atm_cnt_real, sg);
		free(FA->mif_energies); free(FA->mif_sorted); free(FA->mif_cdf);
		FA->mif_count = static_cast<int>(mif.sorted_indices.size());
		FA->mif_energies = static_cast<float*>(
		    malloc(mif.energies.size() * sizeof(float)));
		FA->mif_sorted = static_cast<int*>(
		    malloc(mif.sorted_indices.size() * sizeof(int)));
		std::copy_n(mif.energies.data(), mif.energies.size(), FA->mif_energies);
		std::copy_n(mif.sorted_indices.data(), mif.sorted_indices.size(), FA->mif_sorted);

		mif::build_sampling_cdf(mif, FA->mif_temperature);
		FA->mif_cdf = static_cast<double*>(
		    malloc(mif.cdf.size() * sizeof(double)));
		std::copy_n(mif.cdf.data(), mif.cdf.size(), FA->mif_cdf);

		// Grid prioritization: filter to top-K%
		if (FA->grid_prio_percent < 100.0f) {
			auto kept = mif::prioritize_grid(mif, FA->grid_prio_percent);
			gridpoint* new_grid = nullptr;
			int new_count = mif::rebuild_cleftgrid(
			    *cleftgrid, FA->num_grd, kept, &new_grid);
			if (new_grid && new_count > 0) {
				int old_count = FA->num_grd;
				free(*cleftgrid);
				*cleftgrid = new_grid;
				FA->num_grd = new_count;
				calc_cleftic(FA, *cleftgrid);
				printf("GRIDPRIO: kept %d/%d grid points (top %.0f%%)\n",
				       new_count - 1, old_count - 1, FA->grid_prio_percent);
			}
		}

		printf("MIF: computed for %d grid points (T=%.0fK)\n",
		       FA->mif_count, FA->mif_temperature);
	}

	// ── RefLig seeding (explicit file OR HETATM fallback) ──
	const char* reflig_path = need_reflig ? FA->reflig_file :
	                          need_hetatm ? lig_file : nullptr;
	if (reflig_path) {
		auto data = reflig::prepare_reflig_seed(
		    reflig_path, *cleftgrid, FA->num_grd, FA->reflig_k_nearest);
		free(FA->reflig_nearest_grid);
		FA->reflig_nearest_count = static_cast<int>(data.nearest_grid.size());
		FA->reflig_nearest_grid = static_cast<int*>(
		    malloc(data.nearest_grid.size() * sizeof(int)));
		std::copy_n(data.nearest_grid.data(), data.nearest_grid.size(),
		            FA->reflig_nearest_grid);

		printf("REFLIG: seeded from %s — centroid (%.1f, %.1f, %.1f), %d nearest points\n",
		       reflig_path, data.centroid[0], data.centroid[1], data.centroid[2],
		       FA->reflig_nearest_count);
	}
}

/*****************************************************************************
 * SUBROUTINE read_input reads input file.
 *****************************************************************************/
void read_input(FA_Global* FA,atom** atoms, resid** residue,rot** rotamer,gridpoint** cleftgrid,char* input_file){

	FILE *infile_ptr;          /* pointer to input file */
	char buffer[MAX_PATH__*2];   /* a line from the INPUT file */
	//char input_file[20];     /* input file name */
	char field[7];             /* field names on INPUT file */
	int  i,k;              /* dumb counter */

	//int  flag;               /* a simple flag */
	char pdb_name[MAX_PATH__];       /* 4 letter PDB filename */
	char lig_file[MAX_PATH__];       /* 4 letter PDB filename */
	char rmsd_file[MAX_PATH__];     /* 4 letter PDB filename */
	char clf_file[MAX_PATH__];      /* Cleft file that contains sphere coordinates */
	char normal_file[MAX_PATH__];   /* normal mode grid file */
	char eigen_file[MAX_PATH__];   /* normal mode grid file */

	char emat_forced[MAX_PATH__];
	char emat[MAX_PATH__];      /* interaction matrix file*/

	char deftyp_forced[MAX_PATH__];
	char deftyp[MAX_PATH__];     /* amino/nucleotide definition file*/

	char constraint_file[MAX_PATH__]; /* constraint file */

	char rotlib_file[MAX_PATH__];  /* rotamer library file */
	char rotobs_file[MAX_PATH__];  /* rotamer observations file */

	char rngopt[7] = "";
	char rngoptline[MAX_PATH__];

	//char anam[5];            /* temporary atom name */
	//char rnam[4];            /* temporary residue name */
	//int  type;               /* temporary atom type */
	//char num_char[7];        /* string to read integer */
	//int  rec[3];             /* temporary array for atom reconstruction data */
	int opt[2];              /* temporary array for reading optimization data*/
	char chain='-';
	char a[7],b[7]; //,mol_name[4]; 

	char gridfile[MAX_PATH__];
	//char sphfile[MAX_PATH__];
	char tmpprotname[MAX_PATH__];

	char optline[MAX_PAR][MAX_PATH__];
	char flexscline[MAX_PAR][MAX_PATH__];
	
	// opt lines counter
	int nopt=0;
	int nflexsc=0;
	
	sphere *spheres, * _sphere;
	
	rmsd_file[0]='\0';
	FA->state_path[0]='\0';
	FA->temp_path[0]='\0';
	constraint_file[0] = '\0';
	deftyp_forced[0] = '\0';
	emat_forced[0] = '\0';
	FA->dependencies_path[0] = '\0';
  
	spheres = NULL;

	infile_ptr=NULL;
	if (!OpenFile_B(input_file,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Could not read input file: %s\n",input_file);
		Terminate(8);
	}

	while (fgets(buffer, sizeof(buffer),infile_ptr)){
		size_t blen = strlen(buffer);
		if (blen > 0 && buffer[blen-1] == '\n')
			buffer[--blen] = '\0';

		if (blen < 6) continue;
		for(i=0;i<6;++i) field[i]=buffer[i];
		field[6]='\0';
		
		if(strcmp(field,"PDBNAM") == 0){strncpy(pdb_name,&buffer[7],MAX_PATH__-1);pdb_name[MAX_PATH__-1]='\0';}
		if(strcmp(field,"INPLIG") == 0){strncpy(lig_file,&buffer[7],MAX_PATH__-1);lig_file[MAX_PATH__-1]='\0';}
		if(strcmp(field,"METOPT") == 0){sscanf(buffer,"%s %2s",a,FA->metopt);}
		if(strcmp(field,"DEEFLX") == 0){FA->deelig_flex=1;}
		if(strcmp(field,"BPKENM") == 0){sscanf(buffer,"%s %2s",a,FA->bpkenm);}
		if(strcmp(field,"COMPLF") == 0){sscanf(buffer,"%s %3s",a,FA->complf);}
		if(strcmp(field,"VCTSCO") == 0){sscanf(buffer,"%s %5s",a,FA->vcontacts_self_consistency);}
		if(strcmp(field,"VCTPLA") == 0){sscanf(buffer,"%s %c",a,&FA->vcontacts_planedef);}
		if(strcmp(field,"NORMAR") == 0){FA->normalize_area=1;}
		if(strcmp(field,"USEACS") == 0){FA->useacs=1;}
		if(strcmp(field,"ACSWEI") == 0){sscanf(buffer,"%s %f",field,&FA->acsweight);}
		if(strcmp(field,"RNGOPT") == 0){strncpy(rngoptline,buffer,MAX_PATH__-1);rngoptline[MAX_PATH__-1]='\0';for(i=0;i<6;i++)rngopt[i]=buffer[7+i];rngopt[6]='\0';}
		if(strcmp(field,"OPTIMZ") == 0){if(nopt<MAX_PAR){strncpy(optline[nopt],buffer,MAX_PATH__-1);optline[nopt][MAX_PATH__-1]='\0';nopt++;}}
		if(strcmp(field,"FLEXSC") == 0){if(nflexsc<MAX_PAR){strncpy(flexscline[nflexsc],buffer,MAX_PATH__-1);flexscline[nflexsc][MAX_PATH__-1]='\0';nflexsc++;}}
		if(strcmp(field,"ROTOBS") == 0){FA->rotobs=1;}
		if(strcmp(field,"DEFTYP") == 0){strncpy(deftyp_forced,&buffer[7],MAX_PATH__-1);deftyp_forced[MAX_PATH__-1]='\0';}
		if(strcmp(field,"CLRMSD") == 0){sscanf(buffer,"%s %f",a,&FA->cluster_rmsd);}
		if(strcmp(field,"SUPCLU") == 0){FA->use_super_cluster=true;}
		if(strcmp(field,"TQCM__") == 0){FA->use_tqcm=true;}
		if(strcmp(field,"TQENS_") == 0){FA->use_tqens=true;}
		if(strcmp(field,"TQNN__") == 0){FA->use_tqnn=true;}
		if(strcmp(field,"MULTIM") == 0){FA->multi_model=true;}  // CCBM: multi-model ensemble docking
		if(strcmp(field,"ROTOUT") == 0){FA->rotout=1;}
		if(strcmp(field,"NMAMOD") == 0){sscanf(buffer,"%s %d",a,&FA->normal_modes);}
		if(strcmp(field,"NMAAMP") == 0){strncpy(normal_file,&buffer[7],MAX_PATH__-1);normal_file[MAX_PATH__-1]='\0';}
		if(strcmp(field,"NMAEIG") == 0){strncpy(eigen_file,&buffer[7],MAX_PATH__-1);eigen_file[MAX_PATH__-1]='\0';}
		if(strcmp(field,"RMSDST") == 0){strncpy(rmsd_file,&buffer[7],MAX_PATH__-1);rmsd_file[MAX_PATH__-1]='\0';}
		if(strcmp(field,"EXCHET") == 0){FA->exclude_het=1;}
		if(strcmp(field,"INCHOH") == 0){FA->remove_water=0;}
		if(strcmp(field,"NOINTR") == 0){FA->intramolecular=0;}
		if(strcmp(field,"OMITBU") == 0){FA->omit_buried=1;}
		if(strcmp(field,"VINDEX") == 0){FA->vindex=1;}
		if(strcmp(field,"HTPMOD") == 0){FA->htpmode=true;}
		if(strcmp(field,"PERMEA") == 0){sscanf(buffer,"%s %f",field,&FA->permeability);}
		if(strcmp(field,"INTRAF") == 0){sscanf(buffer,"%s %f",field,&FA->intrafraction);}
		if(strcmp(field,"VARDIS") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_angstron);}
		if(strcmp(field,"VARANG") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_angle);}
		if(strcmp(field,"VARDIH") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_dihedral);}
		if(strcmp(field,"VARFLX") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_flexible);}
		if(strcmp(field,"SLVPEN") == 0){sscanf(buffer,"%s %f",field,&FA->solventterm);}
		if(strcmp(field,"USEELC") == 0){FA->use_elec=1;}
		if(strcmp(field,"DIELEC") == 0){sscanf(buffer,"%s %f",field,&FA->dielectric);}
		if(strcmp(field,"USEGIS") == 0){FA->use_gist=1;}
		if(strcmp(field,"GISTDG") == 0){sscanf(buffer,"%s %s",field,FA->gist_dg_file);}
		if(strcmp(field,"GISTDN") == 0){sscanf(buffer,"%s %s",field,FA->gist_dens_file);}
		if(strcmp(field,"GISTWT") == 0){sscanf(buffer,"%s %lf",field,&FA->gist_weight);}
		if(strcmp(field,"GISTGC") == 0){sscanf(buffer,"%s %f",field,&FA->gist_dg_cutoff);}
		if(strcmp(field,"GISTRC") == 0){sscanf(buffer,"%s %f",field,&FA->gist_rho_cutoff);}
		if(strcmp(field,"GISTDV") == 0){sscanf(buffer,"%s %f",field,&FA->gist_divisor);}
		if(strcmp(field,"USEHBD") == 0){FA->use_hbond=1;}
		if(strcmp(field,"HBONDW") == 0){sscanf(buffer,"%s %lf",field,&FA->hbond_weight);}
		if(strcmp(field,"USEMTL") == 0){FA->use_metal_coord=1;}
		if(strcmp(field,"MTLWGT") == 0){sscanf(buffer,"%s %lf",field,&FA->metal_coord_weight);}
		if(strcmp(field,"MTLSIG") == 0){sscanf(buffer,"%s %lf",field,&FA->metal_coord_sigma);}
		if(strcmp(field,"MTLCNW") == 0){sscanf(buffer,"%s %lf",field,&FA->metal_coord_cn_weight);}
		if(strcmp(field,"OUTRNG") == 0){FA->output_range=1;}
		if(strcmp(field,"USEDEE") == 0){FA->useflexdee=1;}
		if(strcmp(field,"IMATRX") == 0){strncpy(emat_forced,&buffer[7],MAX_PATH__-1);emat_forced[MAX_PATH__-1]='\0';}
		if(strcmp(field,"DEECLA") == 0){sscanf(buffer,"%s %f",field,&FA->dee_clash);}
		if(strcmp(field,"ROTPER") == 0){sscanf(buffer,"%s %f",field,&FA->rotamer_permeability);}
		if(strcmp(field,"CONSTR") == 0){strncpy(constraint_file,&buffer[7],MAX_PATH__-1);constraint_file[MAX_PATH__-1]='\0';}
		if(strcmp(field,"MAXRES") == 0){sscanf(buffer,"%s %d",field,&FA->max_results);}
		if(strcmp(field,"SPACER") == 0){sscanf(buffer,"%s %f",field,&FA->spacer_length);}
		if(strcmp(field,"DEPSPA") == 0){strncpy(FA->dependencies_path,&buffer[7],MAX_PATH__-1);FA->dependencies_path[MAX_PATH__-1]='\0';}
		if(strcmp(field,"STATEP") == 0){strncpy(FA->state_path,&buffer[7],MAX_PATH__-1);FA->state_path[MAX_PATH__-1]='\0';}
		if(strcmp(field,"TEMPOP") == 0){strncpy(FA->temp_path,&buffer[7],MAX_PATH__-1);FA->temp_path[MAX_PATH__-1]='\0';}
		if(strcmp(field,"NRGSUI") == 0){FA->nrg_suite=1;}
		if(strcmp(field,"NRGOUT") == 0){sscanf(buffer,"%s %d",field,&FA->nrg_suite_timeout);}
		if(strcmp(field,"SCOLIG") == 0){FA->score_ligand_only=1;}
		if(strcmp(field,"SCOOUT") == 0){FA->output_scored_only=1;}
		if(strcmp(field,"TEMPER") == 0)
		{
			sscanf(buffer, "%s %u", field, &FA->temperature);
			if(FA->temperature > 0) 
			{
				FA->beta = (double) (1.0 / FA->temperature);
			}
			else FA->beta = 0;
		}
		if(strcmp(field,"CLUSTA") == 0)
		{
			if(FA->temperature > 0) sscanf(buffer, "%s %2s", field, FA->clustering_algorithm);
			else
			{
				fprintf(stdout,"Overriding the clustering algorithm to CF as the Temperature given in input parameter does not allow the consideration of conformational entropy.\n");
				strcpy(FA->clustering_algorithm, "CF");
			}
			
			if(strncmp(FA->clustering_algorithm,"FO",2) != 0 && strncmp(FA->clustering_algorithm,"DP",2) != 0 && strncmp(FA->clustering_algorithm,"CF",2) != 0)
			{
				fprintf(stderr,"ERROR: Invalid clustering algorithm given in input parameter.\n");
				Terminate(2);
			}
		}
	}
	
	CloseFile_B(&infile_ptr,"r");
  

	//////////////////////////////////////////////////
	////////// read input files afterwards ///////////
	//////////////////////////////////////////////////
  
	printf("dependencies path=%s\n", FA->dependencies_path);

	// default state path (controls pause-stop-abort)
	if(!strcmp(FA->state_path,"")){
		strncpy(FA->state_path,FA->base_path,MAX_PATH__-1);
		FA->state_path[MAX_PATH__-1]='\0';
	}

	if(!strcmp(FA->temp_path,"")){
		strncpy(FA->temp_path,FA->base_path,MAX_PATH__-1);
		FA->temp_path[MAX_PATH__-1]='\0';
	}

	// temporary pdb name
#ifdef _WIN32
	snprintf(tmpprotname,MAX_PATH__,"%s\\target.pdb",FA->temp_path);
#else
	snprintf(tmpprotname,MAX_PATH__,"%s/target.pdb",FA->temp_path);
#endif


	/************************************************************/
	/********          INTERACTION MATRIX              **********/
	/************************************************************/
	if(!strcmp(emat_forced,"")){
		// use default
		const char* emat_base = !strcmp(FA->dependencies_path,"") ? FA->base_path : FA->dependencies_path;
#ifdef _WIN32
		snprintf(emat,MAX_PATH__,"%s\\MC_st0r5.2_6.dat",emat_base);
#else
		snprintf(emat,MAX_PATH__,"%s/MC_st0r5.2_6.dat",emat_base);
#endif
	}else{
		// use forced matrix
		strncpy(emat,emat_forced,MAX_PATH__-1);emat[MAX_PATH__-1]='\0';
	}
	
	printf("interaction matrix is <%s>\n", emat);
	read_emat(FA,emat);

	printf("pdb target is <%s>\n", pdb_name);
	if(rna_structure(pdb_name)){
		printf("target molecule is a RNA structure\n");
		FA->is_protein = 0;
	}
	
	/************************************************************/
	/********          DEFINITION OF TYPES             **********/
	/************************************************************/
	if(!strcmp(deftyp_forced,"")){
		// use default definition
		const char* deftyp_base = !strcmp(FA->dependencies_path,"") ? FA->base_path : FA->dependencies_path;
		const char* deftyp_file = FA->is_protein ? "AMINO.def" : "NUCLEOTIDES.def";
#ifdef _WIN32
		snprintf(deftyp,MAX_PATH__,"%s\\%s",deftyp_base,deftyp_file);
#else
		snprintf(deftyp,MAX_PATH__,"%s/%s",deftyp_base,deftyp_file);
#endif
	}else{
		// use forced definition of types
		strncpy(deftyp,deftyp_forced,MAX_PATH__-1);deftyp[MAX_PATH__-1]='\0';
	}

	printf("definition of types is <%s>\n", deftyp);

	///////////////////////////////////////////////
  
	printf("read PDB file <%s>\n",pdb_name);

	// Create a copy of pdb_name that we can modify
	strcpy(tmpprotname, pdb_name);

	// Find the last directory separator
	char *filename = tmpprotname;
	char *last_separator = NULL;

	// Find the last directory separator
	char *slash = strrchr(tmpprotname, '/');
	if (slash != NULL) {
		last_separator = slash;
	}
#ifdef _WIN32
	char *backslash = strrchr(tmpprotname, '\\');
	if (backslash != NULL && (slash == NULL || backslash > slash)) {
		last_separator = backslash;
	}
#endif

	// If there's a separator, the filename starts after it
	if (last_separator != NULL) {
		filename = last_separator + 1;
	}

	// Find the first dot in the filename part
	char *dot = strchr(filename, '.');
	if (dot != NULL) {
		*dot = '\0'; // Remove everything after the first dot
	}

	// Generate random 6-digit number
	int random_num = static_cast<int>(std::random_device{}() % 900000 + 100000);

	// If we had a dot, restore the string terminator to its original position
	if (dot != NULL) {
		*dot = '.';
	}

	// Create the new filename with _tmp_random
	char *extension_pos = dot;
	char random_str[24]; // Buffer for "_tmp_XXXXXX.pdb"
	snprintf(random_str, sizeof(random_str), "_tmp_%d.pdb", random_num);

	if (extension_pos != NULL) {
		size_t avail = MAX_PATH__ - (size_t)(extension_pos - tmpprotname);
		strncpy(extension_pos, random_str, avail - 1);
		tmpprotname[MAX_PATH__ - 1] = '\0';
	} else {
		size_t cur = strlen(tmpprotname);
		strncpy(tmpprotname + cur, random_str, MAX_PATH__ - cur - 1);
		tmpprotname[MAX_PATH__ - 1] = '\0';
	}

	modify_pdb(pdb_name,tmpprotname,FA->exclude_het,FA->remove_water,FA->is_protein,
	           FA->keep_ions,FA->keep_structural_waters,FA->structural_water_bfactor_max);
	read_pdb(FA,atoms,residue,tmpprotname);
	remove(tmpprotname);

	(*residue)[FA->res_cnt].latm[0]=FA->atm_cnt;
	for(k=1;k<=FA->res_cnt;k++){
		FA->atm_cnt_real += (*residue)[k].latm[0]-(*residue)[k].fatm[0]+1;
	}
  
	calc_center(FA,*atoms,*residue);
  
	if(FA->is_protein){ residue_conect(FA,*atoms,*residue,deftyp); }

	assign_types(FA,*atoms,*residue,deftyp);
	
	//////////////////////////////////////////////

	printf("read ligand input file <%s>\n",lig_file);
	read_lig(FA,atoms,residue,lig_file);

	//////////////////////////////////////////////
	
	assign_radii_types(FA,(*atoms),(*residue));    
	printf("radii are now assigned\n");
    
	if(strcmp(rmsd_file,"")){
		printf("read rmsd structure <%s>: will match atom numbers\n",rmsd_file);

		int rmsd_atoms = read_rmsdst(FA,*atoms,*residue,rmsd_file);

		if(rmsd_atoms){
			FA->refstructure=1;
			printf("will use %d atoms to calculate RMSD\n", rmsd_atoms);
		}else{
			printf("no atoms is used to calculate RMSD\n");
		}
	}

	//////////////////////////////////////////////
  
	if(FA->normal_modes > 0){
		printf("read files related to NMA\n");
		read_normalgrid(FA,normal_file);
		read_eigen(FA,eigen_file);
		assign_eigen(FA,*atoms,*residue,FA->res_cnt,FA->normal_modes);
	}

  
	//////////////////////////////////////////////

	if((nflexsc || FA->autoflex_enabled) && FA->is_protein){

		if(nflexsc){
			read_flexscfile(FA,*residue,rotamer,flexscline,nflexsc,rotlib_file,rotobs_file);
		}else{
			// Auto-flex only: allocate rotamer array (read_flexscfile normally does this)
			(*rotamer) = (rot*)malloc(FA->MIN_ROTAMER_LIBRARY_SIZE*sizeof(rot));
			if(!(*rotamer)){
				fprintf(stderr,"ERROR: memory allocation error for rotamer\n");
				Terminate(2);
			}
			memset((*rotamer),0,FA->MIN_ROTAMER_LIBRARY_SIZE*sizeof(rot));
		}

		if (FA->rotobs) {

			// use rotamers found in bound conformations (hap2db)
			const char* rotobs_base = !strcmp(FA->dependencies_path,"") ? FA->base_path : FA->dependencies_path;
#ifdef _WIN32
			snprintf(rotobs_file,MAX_PATH__,"%s\\rotobs.lst",rotobs_base);
#else
			snprintf(rotobs_file,MAX_PATH__,"%s/rotobs.lst",rotobs_base);
#endif

			printf("read rotamer observations <%s>\n",rotobs_file);
			read_rotobs(FA,rotamer,rotobs_file);
		}else{

			// use penultimate rotamer library instances
			const char* rotlib_base = !strcmp(FA->dependencies_path,"") ? FA->base_path : FA->dependencies_path;
#ifdef _WIN32
			snprintf(rotlib_file,MAX_PATH__,"%s\\Lovell_LIB.dat",rotlib_base);
#else
			snprintf(rotlib_file,MAX_PATH__,"%s/Lovell_LIB.dat",rotlib_base);
#endif

			printf("read rotamer library <%s>\n",rotlib_file);
			read_rotlib(FA,rotamer,rotlib_file);
		}

		if(FA->nflxsc > 0 && FA->rotlibsize > 0){
			build_rotamers(FA,atoms,*residue,*rotamer);
			//build_close(FA,residue,atoms);
		}
	}

	//////////////////////////////////////////////

	if(strcmp(constraint_file,"")){
		printf("read constraint_file <%s>\n", constraint_file);
		read_constraints(FA,*atoms,*residue,constraint_file);
		assign_constraint_threshold(FA,*atoms,FA->constraints,FA->num_constraints);

	}

	///////////////////////////////////////////////

	if(!strcmp(rngopt,"LOCCEN")){
		strcpy(FA->rngopt,"loccen");
		
		_sphere = (sphere*)malloc(sizeof(sphere));
		if(_sphere == NULL){
			fprintf(stderr,"ERROR: memory allocation error for spheres (LOCCEN)\n");
			Terminate(2);
		}
        
		sscanf(rngoptline,"%s %s %f %f %f %f",
		       a,b,
		       &_sphere->center[0],&_sphere->center[1],&_sphere->center[2],
		       &_sphere->radius);
		_sphere->prev = NULL;
		
		// point to the new sphere created
		spheres = _sphere;
		
		(*cleftgrid) = generate_grid(FA,spheres,(*atoms),(*residue));
		calc_cleftic(FA,*cleftgrid);
		compute_mif_and_reflig(FA, *atoms, cleftgrid, lig_file);

	}else if(!strcmp(rngopt,"LOCCLF")){

		//RNGOPT LOCCLF filename.pdb
		strcpy(FA->rngopt,"locclf");
		strcpy(clf_file,&rngoptline[14]);

		printf("read binding-site grid <%s>\n",clf_file);
		spheres = read_spheres(clf_file);

		(*cleftgrid) = generate_grid(FA,spheres,(*atoms),(*residue));
		calc_cleftic(FA,*cleftgrid);
		compute_mif_and_reflig(FA, *atoms, cleftgrid, lig_file);

	}else if(!strcmp(rngopt,"LOCCDT")){

		// RNGOPT LOCCDT [cleft_id] [min_radius] [max_radius]
		// Use native CavityDetector (SURFNET + AVX-512/OpenMP/Metal) to
		// locate the largest binding cleft automatically.
		// Optional args: cleft_id (default 1), probe radii (default 1.4–4.0 Å).
		strcpy(FA->rngopt,"loccdt");

		int   cdt_cleft_id  = 1;
		float cdt_min_r     = 1.4f;
		float cdt_max_r     = 4.0f;
		sscanf(rngoptline,"%*s %d %f %f", &cdt_cleft_id, &cdt_min_r, &cdt_max_r);
		if(cdt_cleft_id < 1) cdt_cleft_id = 1;

		printf("LOCCDT: running native CavityDetector (cleft %d, probe %.2f–%.2f Å)\n",
		       cdt_cleft_id, cdt_min_r, cdt_max_r);

		{
			cavity_detect::CavityDetector detector;
			detector.load_from_fa(*atoms, *residue, FA->res_cnt);
			detector.detect(cdt_min_r, cdt_max_r);

			if(detector.clefts().empty()){
				fprintf(stderr,"ERROR: LOCCDT found no clefts — "
				        "check probe radii or fall back to LOCCLF/LOCCEN\n");
				Terminate(2);
			}

			// Write the detected cleft as a CLF sphere PDB for inspection
			detector.write_sphere_pdb("loccdt_cleft.pdb", cdt_cleft_id);
			const auto& cleft_list = detector.clefts();
			size_t cdt_idx = static_cast<size_t>(cdt_cleft_id - 1);
			printf("LOCCDT: cleft %d detected (%zu spheres) — "
			       "written to loccdt_cleft.pdb\n",
			       cdt_cleft_id,
			       cdt_idx < cleft_list.size() ? cleft_list[cdt_idx].spheres.size() : 0u);

			spheres = detector.to_flexaid_spheres(cdt_cleft_id);
		}

		if(!spheres){
			fprintf(stderr,"ERROR: LOCCDT cleft %d not found\n", cdt_cleft_id);
			Terminate(2);
		}

		// Generate grid from detected cavity spheres (was missing)
		(*cleftgrid) = generate_grid(FA, spheres, (*atoms), (*residue));
		calc_cleftic(FA, *cleftgrid);
		compute_mif_and_reflig(FA, *atoms, cleftgrid, lig_file);

	}else if(!strncmp(rngopt,"AUTO  ",4)){

		// RNGOPT AUTO — automatic cleft detection (SURFNET gap-sphere)
		strcpy(FA->rngopt,"locclf");
		printf("AUTO binding-site detection (CleftDetector) ...\n");

		spheres = detect_cleft(*atoms, *residue, FA->atm_cnt_real, FA->res_cnt);
		if(spheres == NULL){
			fprintf(stderr,"ERROR: AUTO cleft detection found no cavities.\n");
			Terminate(2);
		}

		// Optionally write detected spheres for inspection
		char auto_sph[MAX_PATH__];
#ifdef _WIN32
		snprintf(auto_sph, MAX_PATH__, "%s\\auto_cleft.sph", FA->temp_path);
#else
		snprintf(auto_sph, MAX_PATH__, "%s/auto_cleft.sph", FA->temp_path);
#endif
		write_cleft_spheres(spheres, auto_sph);

		(*cleftgrid) = generate_grid(FA,spheres,(*atoms),(*residue));
		calc_cleftic(FA,*cleftgrid);
		compute_mif_and_reflig(FA, *atoms, cleftgrid, lig_file);
	}
    
	//printf("IC bounds...\n");
	ic_bounds(FA,FA->rngopt);

	///////////////////////////////////////////////////////////////////////////////
    
	for(i=0;i<nopt;i++){
		sscanf(optline[i],"%s %d %s %d",a,&opt[0],a,&opt[1]);
		//printf("%d %d\n",opt[0],opt[1]);
		//getchar();
		//chain=buffer[11];
		chain=a[0];
		if(chain == '-'){chain = ' ';}
		//printf("Add2 optimiz vector...\n");
		add2_optimiz_vec(FA,*atoms,*residue,opt,chain,"");

	}

	// ── Auto-flex: add key binding residues as flexible side-chains ──
	if (FA->autoflex_enabled && FA->mif_energies && FA->mif_count > 0) {
		int nflxsc_before = FA->nflxsc;
		int n_added = binding_residues::add_key_residues_as_flexible(
			FA, *cleftgrid, *atoms, *residue, FA->autoflex_max);

		// Build rotamers for newly added flexible residues only
		if (n_added > 0 && FA->rotlibsize > 0) {
			int saved_nflxsc = FA->nflxsc;
			int saved_nflxsc_real = FA->nflxsc_real;
			flxsc* saved_flex_res = FA->flex_res;

			// Point flex_res to only the new entries so build_rotamers
			// processes only the auto-flexed residues
			FA->flex_res = &saved_flex_res[nflxsc_before];
			FA->nflxsc = n_added;

			build_rotamers(FA, atoms, *residue, *rotamer);

			// Accumulate nflxsc_real from both batches
			int new_nflxsc_real = FA->nflxsc_real;

			// Restore full flex_res array
			FA->flex_res = saved_flex_res;
			FA->nflxsc = saved_nflxsc;
			FA->nflxsc_real = saved_nflxsc_real + new_nflxsc_real;

			printf("AUTOFLEX: built rotamers for %d auto-flexed residues "
			       "(%d with rotamers, total flexible: %d)\n",
			       n_added, new_nflxsc_real, FA->nflxsc);
		}
	}

	add2_optimiz_vec(FA,*atoms,*residue,opt,chain,"SC");
	add2_optimiz_vec(FA,*atoms,*residue,opt,chain,"NM");
    
	//////////////////////////////////////////////
    
	if(FA->translational && strcmp(rngopt,"LOCCEN") && strcmp(rngopt,"LOCCLF")){
		fprintf(stderr,"ERROR: the binding-site is not defined\n");
		Terminate(2);
	}
    
	if(FA->output_range){
#ifdef _WIN32
		snprintf(gridfile,MAX_PATH__,"%s\\grid.sta.pdb",FA->temp_path);
#else
		snprintf(gridfile,MAX_PATH__,"%s/grid.sta.pdb",FA->temp_path);
#endif
		write_grid(FA,*cleftgrid,gridfile);
	}
    
	if(FA->translational && FA->num_grd==1){
		fprintf(stderr,"ERROR: the binding-site has no anchor points\n");
		Terminate(2);
	}
	
	// fill in optres pointer in atoms struct.
	update_optres(*atoms,*residue,FA->atm_cnt,FA->optres,FA->num_optres);
	
	if(FA->nrg_suite){
		if(FA->translational){
			for(i=1; i<FA->num_grd; i++){
				printf("Grid[%d]=%8.3f%8.3f%8.3f\n", i, (*cleftgrid)[i].coor[0], (*cleftgrid)[i].coor[1], (*cleftgrid)[i].coor[2]);
				if(i % 1000 == 0){
					fflush(stdout);
				}
				fflush(stdout);
			}
		}
	}
    
	/* FREE SPHERES LINKED-LIST */
	while(spheres != NULL){
		_sphere = spheres->prev;
		free(spheres);

		spheres = _sphere;
	}

	return;
}
