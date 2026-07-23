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

	// Set a Hydrophilic type to water molecules
	for(k=1;k<=FA->res_cnt;k++){
		rot=residue[k].rot;
		for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){
			if(strcmp(residue[atoms[i].ofres].name,"HOH") == 0){atoms[i].type = 1;}
		}
	}

	// Assign atom types for metal ions by residue name (bounded by FA->ntypes).
	// These are *canonical VCT rows* (the table above), NOT the base-class enum
	// in atom_typing_256.h — that enum numbers metals differently (Metal_Mn=36,
	// Metal_Hg=37, …) and reconciling the two tables would mis-index the energy
	// matrix.  MG=28, CU=30, MN=31, HG=32, CD=33, NI=34, ZN=35, CA=36, FE=37.
	// Note rows 32 (Hg) and 33 (Cd) are near-dead in MC_st0r5.2_6.dat (0 and 1
	// non-zero entries): those ions are effectively steric-only.  Fall back to
	// dummy (ntypes-1) if the matrix is smaller.
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
