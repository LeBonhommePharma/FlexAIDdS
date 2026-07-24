#include "flexaid.h"
#include "fileio.h"

/******************************************************************************
 * SUBROUTINE assign_types assign atom types for protein and hetero group atoms
 ******************************************************************************/
void assign_types(FA_Global* FA,atom* atoms,resid* residue,char aminofile[]){

	FILE *infile_ptr;        /* pointer to input file */
	char buffer[81];         /* a line from the INPUT file */
	char bufnul[7];          /* dumb string to read input strings */
	int  i,j,k;                /* dumb counter */
	char anam[5];            /* temporary atom name */
	char rnam[4],rnamrev[4]; /* temporary residue name */
	int  type;               /* temporary atom type */
	//char field[7];
	char recs;
	int  rec[4];
	int rot;

	infile_ptr=NULL;
	if (!OpenFile_B(aminofile,"r",&infile_ptr))
		Terminate(8);

	while (fgets(buffer, sizeof(buffer),infile_ptr)){

		for(i=0;i<=5;i++){bufnul[i]=buffer[i];}
		bufnul[6]='\0';

		if(strcmp(bufnul,"RESIDU") == 0 || strcmp(bufnul,"NUCLEO") == 0){ 
			for(i=0;i<=2;i++){rnam[i]=buffer[i+7];}
			rnam[3]='\0';

			rnamrev[0]=rnam[2];
			rnamrev[1]=rnam[1];
			rnamrev[2]=rnam[0];
			rnamrev[3]='\0';
			//printf("%s\n",rnam);
		}

		if(strcmp(bufnul,"ATMTYP") == 0){ 
			bufnul[0]=buffer[10];
			bufnul[1]=buffer[11];
			bufnul[2]='\0';
			sscanf(bufnul,"%d",&type);
      
			for(i=0;i<=3;i++){anam[i]=buffer[i+12];}
			anam[4]='\0';
      
			if(type > FA->ntypes){
				printf("WARNING: res %s atom %s has atom type %d when %d types are defined\n",
				       rnam, anam, type, FA->ntypes);
				printf("WARNING: type %d is set to neutral (6)\n", type);

				type = 6;
			}


			recs=buffer[17];

			/*printf("recs:%c\n",recs);*/

			if(recs == 'm'){
				for(i=0;i<=3;i++){
					for(j=0;j<=2;j++){bufnul[j]=buffer[18+i*3+j];}
					bufnul[3]='\0';
					sscanf(bufnul,"%d",&rec[i]);
				}
				/*printf("%d %d %d\n",rec[0],rec[1],rec[2]);*/
			}
      
			/*PAUSE;*/
			for(k=1;k<=FA->res_cnt;k++){
				rot=residue[k].rot;
				for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){

					if(
						(strcmp(residue[atoms[i].ofres].name,rnam) == 0 ||
						 (!FA->is_protein && strcmp(residue[atoms[i].ofres].name,rnamrev) == 0)) &&
						strcmp(atoms[i].name,anam) == 0){
						
						atoms[i].type=type;
						atoms[i].recs=recs;
						
						if(atoms[i].recs == 'm'){
							for(j=0;j<=2;j++){
								atoms[i].rec[j]=rec[j]+
									residue[atoms[i].ofres].fatm[residue[atoms[i].ofres].rot]-1;
							}
							if(rec[3] != 0){
								atoms[i].rec[3]=rec[3]+
									residue[atoms[i].ofres].fatm[residue[atoms[i].ofres].rot]-1;
							}else{
								atoms[i].rec[3]=0;
							}
						}else{
							for(j=0;j<=3;j++){
								atoms[i].rec[j]=0;
							}
						}
					}
				}
			}
		}
	}
	CloseFile_B(&infile_ptr,"r");

	// ── Retained crystallographic water ────────────────────────────────────────
	// Canonical VCT row numbering (identical in AMINO.def, MC_st0r5.2_6.dat,
	// read_coor.cpp:canonical_vct_type_for_element, Mol2Reader.cpp:
	// sybyl_to_flexaid_type and top.cpp:sybyl_name_to_canonical_vct):
	//    1=C.1   2=C.2   3=C.3   4=C.ar  5=C.cat
	//    6=N.1   7=N.2   8=N.3   9=N.4  10=N.ar 11=N.am 12=N.pl3
	//   13=O.2  14=O.3  15=O.co2 16=O.ar
	//   17=S.2  18=S.3  19=S.O  20=S.O2 21=S.ar
	//   22=P.3  23=F    24=Cl   25=Br   26=I    27=Se
	//   28=Mg 29=Sr 30=Cu 31=Mn 32=Hg 33=Cd 34=Ni 35=Zn 36=Ca 37=Fe 38=Co.oh
	//   39=DUMMY  40=SOLVENT
	//
	// Water oxygen is an sp3 oxygen bearing two hydrogens, i.e. chemically an
	// O.3 (row 14) — the same row SER-OG / THR-OG1 / TYR-OH occupy above, and
	// the row read_coor.cpp already derives from the element for any retained
	// HETATM.  This loop used to overwrite that correct value with type 1
	// ("Set a Hydrophilic type"), which under the canonical numbering is C.1,
	// *sp carbon*.  Row 1's two strongest partners are 1-13 (C.1 x O.2) =
	// -198.3 and 1-14 (C.1 x O.3) = -180.8 — the most attractive cells in the
	// entire matrix — so every retained water radiated a large spurious
	// attraction to ligand and protein oxygen.  That is the source of the
	// CF.com blow-up documented in DatasetRunner.cpp (~-4269 from sub-Angstrom
	// ligand-O ... HOH-O contacts): the "C x O.3" dominant contact was water.
	//
	// Row 40 (SOLVENT) is NOT a valid alternative: it is a reserved pseudo-type
	// for the bulk-solvent / SAS desolvation channel, indexed as (ntypes-1) in
	// vcfunction.cpp, and it is neither empty nor neutral — it carries 20
	// non-zero, overwhelmingly repulsive entries (10-40 = +198.3, 24-40 =
	// +198.8).  Typing explicit waters as 40 would both double-count the
	// desolvation channel and make every structural water strongly repulsive.
	//
	// All water residue aliases are normalised here so a "WAT"/"H2O" water is
	// scored identically to an "HOH" water.  Water hydrogens (and deuterium in
	// neutron structures) stay DUMMY, as canonical_vct_type_for_element gives.
	{
		static const char* water_res[] = {"HOH","WAT","H2O","DOD","OHX",nullptr};
		for(k=1;k<=FA->res_cnt;k++){
			const char* rname = residue[k].name;
			bool is_water = false;
			for(int n=0; water_res[n]; ++n)
				if(!strncmp(rname, water_res[n], 3)){ is_water = true; break; }
			if(!is_water) continue;

			rot=residue[k].rot;
			for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){
				const char* el = atoms[i].element;
				if(el[0]=='H' && el[1]=='\0') continue;         // H / D stays DUMMY
				if(el[0]=='D' && el[1]=='\0') continue;
				if(14 <= FA->ntypes) atoms[i].type = 14;        // O.3
			}
		}
	}

	// Assign SYBYL atom types for metal ions by residue name (bounded by FA->ntypes)
	// SYBYL numbering from atom_typing_256.h: MG=28, CU=30, MN=31, HG=32, CD=33,
	// NI=34, ZN=35, CA=36, FE=37.  Fall back to dummy (ntypes-1) if matrix is smaller.
	{
		static const struct { const char* rnam; int sybyl; } ion_types[] = {
			{"MG ", 28}, {"CU ", 30}, {"CU1", 30}, {"CU2", 30},
			{"MN ", 31}, {"HG ", 32}, {"CD ", 33}, {"NI ", 34},
			{"ZN ", 35}, {"CA ", 36}, {"FE ", 37}, {"FE2", 37}, {"FE3", 37},
			{nullptr, 0}
		};
		for(k=1; k<=FA->res_cnt; k++){
			if(residue[k].type != 1) continue;  // HETATM only
			rot = residue[k].rot;
			const char* rname = residue[k].name;
			int sybyl = 0;
			for(int n = 0; ion_types[n].rnam; ++n)
				if(!strncmp(rname, ion_types[n].rnam, 3)){ sybyl = ion_types[n].sybyl; break; }
			if(sybyl > 0 && sybyl <= FA->ntypes)
				for(i = residue[k].fatm[rot]; i <= residue[k].latm[rot]; ++i)
					atoms[i].type = sybyl;
		}
	}

	return;
}
