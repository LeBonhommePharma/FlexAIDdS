#include "flexaid.h"
#include "fileio.h"

#include <cctype>
#include <cstring>

namespace {

void pdb_element(const char* line, const char atom_name[5], char out[3])
{
	out[0] = '\0';
	out[1] = '\0';
	out[2] = '\0';

	const std::size_t len = std::strlen(line);
	char raw[3] = {'\0', '\0', '\0'};
	if (len > 76) raw[0] = line[76];
	if (len > 77) raw[1] = line[77];

	int n = 0;
	for (char c : raw) {
		if (c != '\0' && !std::isspace(static_cast<unsigned char>(c)) && n < 2)
			out[n++] = c;
	}

	// Element columns are optional in older PDBs. The PDB atom-name alignment
	// makes a leading blank the reliable one-letter-element case.
	if (n == 0) {
		int first = 0;
		while (first < 4 && (atom_name[first] == ' ' ||
		                     std::isdigit(static_cast<unsigned char>(atom_name[first]))))
			++first;
		if (first < 4) out[n++] = atom_name[first];
		if (first == 0 && first + 1 < 4 &&
		    std::isalpha(static_cast<unsigned char>(atom_name[first + 1])) &&
		    std::islower(static_cast<unsigned char>(atom_name[first + 1])))
			out[n++] = atom_name[first + 1];
	}

	if (n > 0) out[0] = static_cast<char>(
		std::toupper(static_cast<unsigned char>(out[0])));
	if (n > 1) out[1] = static_cast<char>(
		std::tolower(static_cast<unsigned char>(out[1])));
}

int canonical_vct_type_for_element(const char element[3], int ntypes)
{
	int type = ntypes > 1 ? ntypes - 1 : 1; // DUMMY fallback
	if      (!std::strcmp(element, "C"))  type = 3;
	else if (!std::strcmp(element, "N"))  type = 11;
	else if (!std::strcmp(element, "O"))  type = 14;
	else if (!std::strcmp(element, "S"))  type = 18;
	else if (!std::strcmp(element, "P"))  type = 22;
	else if (!std::strcmp(element, "F"))  type = 23;
	else if (!std::strcmp(element, "Cl")) type = 24;
	else if (!std::strcmp(element, "Br")) type = 25;
	// Iodine and selenium must alias exactly as the ligand readers do, otherwise
	// the same element scores differently depending on which side of the complex
	// it sits on. I → Br (row 26 has only 3 live entries); Se → S.3 (row 27 is
	// all-zero, and Se reaches the receptor as selenomethionine, a Met surrogate).
	// Keep in lockstep with Mol2Reader::sybyl_to_flexaid_type,
	// SdfReader::element_to_flexaid_type and top.cpp:sybyl_name_to_canonical_vct.
	else if (!std::strcmp(element, "I"))  type = 25;
	else if (!std::strcmp(element, "Se")) type = 18;
	else if (!std::strcmp(element, "Mg")) type = 28;
	else if (!std::strcmp(element, "Sr")) type = 29;
	else if (!std::strcmp(element, "Cu")) type = 30;
	else if (!std::strcmp(element, "Mn")) type = 31;
	else if (!std::strcmp(element, "Hg")) type = 32;
	else if (!std::strcmp(element, "Cd")) type = 33;
	else if (!std::strcmp(element, "Ni")) type = 34;
	else if (!std::strcmp(element, "Zn")) type = 35;
	else if (!std::strcmp(element, "Ca")) type = 36;
	else if (!std::strcmp(element, "Fe")) type = 37;
	else if (!std::strcmp(element, "Co")) type = 38;
	return (type >= 1 && type <= ntypes) ? type : (ntypes > 1 ? ntypes - 1 : 1);
}

} // namespace

static int standard_polymer_residue(const char rnam[4])
{
	static const char* aa[] = {
		"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
		"LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
		"SEC","PYL", NULL
	};
	for(int i=0; aa[i]; ++i)
		if(strncmp(rnam, aa[i], 3) == 0) return 1;

	static const char* nuc[] = {
		"  A","  C","  G","  T","  U",
		" DA"," DC"," DG"," DT"," DU",
		"A  ","C  ","G  ","T  ","U  ", NULL
	};
	for(int i=0; nuc[i]; ++i)
		if(strncmp(rnam, nuc[i], 3) == 0) return 1;
	return 0;
}

/***************************************************************************** 
 * SUBROUTINE read coor, gets a coordinates line from read_pdb and extracts
 * the coordinates, atom name, chain name, counts the number of ligands, 
 * residues, the initial and final atom of each residue and hetero group
 *****************************************************************************/
void read_coor(FA_Global* FA,atom** atoms,resid** residue,char line[], char res_numold[]){
	char  name[7];             /* 6 letter code of field name on PDB file, e.g. HETATM */
	char  coor_char[10];        /* string used to read the coordinates */
	char  num_char[6];         /* string to read the atom number field */
	char  atm_typ[5];          // temporary type
	char  res_new[4];          // temporary residue name 
	char  res_num[5];          // temporary residue number from the PDB annotation 
	//char  res_numold[5];       // another temporary residue number from the annot.
  
	int   i,j;                  /* dumb counters */

	/*

	  01234567890123456789012345678901234567890123456789012345678901234567890123456789
	  ATOM    239  CB  ALA A  46      17.761  -3.260 -12.974  1.00 20.41           C  
    
	*/
  
	for(i=0;i<=5;i++){name[i]=line[i];}
	name[6]='\0';
  
	for(i=0;i<=3;i++){atm_typ[i]=line[i+12];}
	atm_typ[4]='\0';

	if(line[16]==' ' || line[16]=='A'){
		FA->atm_cnt++;
    
		if(FA->atm_cnt==FA->MIN_NUM_ATOM){
			//printf("re-allocating memory for atoms\n");
			FA->MIN_NUM_ATOM*=2;

			(*atoms) = (atom*)realloc((*atoms),FA->MIN_NUM_ATOM*sizeof(atom));
      
			if(!(*atoms)){
				fprintf(stderr,"ERROR: memory allocation error for atoms.\n");
				Terminate(2);
			}
      
			memset(&(*atoms)[FA->MIN_NUM_ATOM/2],0,FA->MIN_NUM_ATOM/2*sizeof(atom));
			//printf("memory re-allocated for atoms\n");
		}
    
		// dummy atom type by default
		// dummy type is always second last (solvent is last)
		(*atoms)[FA->atm_cnt].type = FA->ntypes-1;

		(*atoms)[FA->atm_cnt].recs = 'r';
		(*atoms)[FA->atm_cnt].eigen = NULL;
		(*atoms)[FA->atm_cnt].ncons=0;
		(*atoms)[FA->atm_cnt].cons=NULL;
		(*atoms)[FA->atm_cnt].optres=NULL;
		(*atoms)[FA->atm_cnt].par=NULL;
		(*atoms)[FA->atm_cnt].graph=0;
		(*atoms)[FA->atm_cnt].coor_ref=NULL;
		(*atoms)[FA->atm_cnt].acs=-1.0f;		
		
		strncpy((*atoms)[FA->atm_cnt].name,atm_typ,sizeof((*atoms)[FA->atm_cnt].name)-1);
		(*atoms)[FA->atm_cnt].name[sizeof((*atoms)[FA->atm_cnt].name)-1]='\0';
		if(strcmp((*atoms)[FA->atm_cnt].name," OXT")==0){
			(*residue)[FA->res_cnt].ter = 1;
			//printf("Residue Ter: %d\n", (*residue)[FA->res_cnt].ter);
		}

		(*atoms)[FA->atm_cnt].isbb=0;
		
		if (!strcmp((*atoms)[FA->atm_cnt].name," CB ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," CA ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," N  ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," O  ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," C  ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," OXT")) 
		{ (*atoms)[FA->atm_cnt].isbb=1; }
		

		//(*atoms)[FA->atm_cnt].radius=assign_radius((*atoms)[FA->atm_cnt].name);
		
			char element[3];
			pdb_element(line, atm_typ, element);
			std::strncpy((*atoms)[FA->atm_cnt].element, element,
			             sizeof((*atoms)[FA->atm_cnt].element) - 1);
			(*atoms)[FA->atm_cnt].element[
				sizeof((*atoms)[FA->atm_cnt].element) - 1] = '\0';

			char type_field[3] = {'\0', '\0', '\0'};
			const std::size_t line_len = std::strlen(line);
			if (line_len > 76) type_field[0] = line[76];
			if (line_len > 77) type_field[1] = line[77];
			int type = std::atoi(type_field);
			if(type && type >= 1 && type <= FA->ntypes) {
				(*atoms)[FA->atm_cnt].type = type;
			} else {
				// Standard polymer atoms are refined by assign_types(). Retained
				// receptor HETATMs are not, so their initial type must already use
				// the canonical VCT matrix numbering.
				(*atoms)[FA->atm_cnt].type =
					canonical_vct_type_for_element(element, FA->ntypes);
			}
		
		for(j=0;j<=4;j++){num_char[j]=line[j+6];}
		num_char[5]='\0';
		sscanf(num_char,"%d",&i);
    
		(*atoms)[FA->atm_cnt].number=i;

		/* maps the PDB num into the internal counter */
		FA->num_atm[i]=FA->atm_cnt;

		for (j=0;j<=2;j++){
			for(i=0;i<=7;i++){
				coor_char[i]=line[30+i+j*8];
			}
			coor_char[8]='\0';
			sscanf(coor_char,"%f",&(*atoms)[FA->atm_cnt].coor[j]);
			sscanf(coor_char,"%f",&(*atoms)[FA->atm_cnt].coor_ori[j]);
		}
    
		for(i=0;i<=2;i++){res_new[i]=line[i+17];}
		res_new[3]='\0';
		for(i=0;i<=3;i++){res_num[i]=line[i+22];}
		res_num[4]='\0';
    
		int duplicate_nonpolymer_atom = 0;
		if(FA->res_cnt > 0 &&
		   strcmp(res_new,(*residue)[FA->res_cnt].name) == 0 &&
		   strcmp(res_num,res_numold) == 0 &&
		   line[21] == (*residue)[FA->res_cnt].chn &&
		   line[26] == (*residue)[FA->res_cnt].ins &&
		   !standard_polymer_residue(res_new)){
			for(int a=(*residue)[FA->res_cnt].fatm[0]; a<FA->atm_cnt; ++a){
				if(strcmp((*atoms)[a].name,atm_typ) == 0){
					duplicate_nonpolymer_atom = 1;
					break;
				}
			}
		}

		if(strcmp(res_new,(*residue)[FA->res_cnt].name) != 0   /* change of res name        */
		   || strcmp(res_num,res_numold) != 0           /* change of res number      */
		   || line[21] != (*residue)[FA->res_cnt].chn   /* change of chain           */
		   || line[26] != (*residue)[FA->res_cnt].ins   /* change of insertion code  */
		   || duplicate_nonpolymer_atom                 /* converted HETATM residue  */
			){
			FA->res_cnt++;

			if(FA->res_cnt==FA->MIN_NUM_RESIDUE){
				//printf("re-allocating memory for residue\n");
				FA->MIN_NUM_RESIDUE *= 2;
				(*residue) = (resid*)realloc((*residue),FA->MIN_NUM_RESIDUE*sizeof(resid));
				if(!(*residue)){
					fprintf(stderr,"ERROR: memory allocation error for residue.\n");
					Terminate(2);
				}
				memset(&(*residue)[FA->MIN_NUM_RESIDUE/2],0,FA->MIN_NUM_RESIDUE/2*sizeof(residue));
				//printf("memory re-allocated for residue\n");
			}

			// Guard against zero-size malloc when callers memset FA to 0
			// without setting top.cpp defaults (MIN_ROTAMER=1, MIN_FLEX_BONDS=5).
			// malloc(0) is implementation-defined and writing fatm[0]/latm[0]
			// then becomes a heap-buffer-overflow under ASan.
			if (FA->MIN_ROTAMER < 1) FA->MIN_ROTAMER = 1;
			if (FA->MIN_FLEX_BONDS < 1) FA->MIN_FLEX_BONDS = 5;

			(*residue)[FA->res_cnt].fatm = (int*)malloc(FA->MIN_ROTAMER*sizeof(int));
			(*residue)[FA->res_cnt].latm = (int*)malloc(FA->MIN_ROTAMER*sizeof(int));
			(*residue)[FA->res_cnt].bond = (int*)malloc(FA->MIN_FLEX_BONDS*sizeof(int));
			if(!(*residue)[FA->res_cnt].fatm ||
			   !(*residue)[FA->res_cnt].latm ||
			   !(*residue)[FA->res_cnt].bond){
				fprintf(stderr,"ERROR: memory allocation error for residue.fatm || .latm || .bond.\n");
				Terminate(2);
			}

			memset((*residue)[FA->res_cnt].fatm,0,FA->MIN_ROTAMER*sizeof(int));
			memset((*residue)[FA->res_cnt].latm,0,FA->MIN_ROTAMER*sizeof(int));
			memset((*residue)[FA->res_cnt].bond,0,FA->MIN_FLEX_BONDS*sizeof(int));

			(*residue)[FA->res_cnt].ter = 0;
			(*residue)[FA->res_cnt].rot=0;
			(*residue)[FA->res_cnt].fdih=0;
			(*residue)[FA->res_cnt].trot = 0;
			(*residue)[FA->res_cnt].gpa=NULL;
			(*residue)[FA->res_cnt].bonded=NULL;
			(*residue)[FA->res_cnt].fatm[0]=FA->atm_cnt;     
			(*residue)[FA->res_cnt-1].latm[0]=FA->atm_cnt-1;
			(*residue)[FA->res_cnt].chn=line[21];
			(*residue)[FA->res_cnt].ins=line[26];
      
			strcpy((*residue)[FA->res_cnt].name,res_new);
			strcpy(res_numold,res_num);

			if(strcmp(name,"ATOM  ")==0) {
				(*residue)[FA->res_cnt].type =
					standard_polymer_residue(res_new) ? 0 : 1;
			}else if(strcmp(name,"HETATM")==0) {
				(*residue)[FA->res_cnt].type=1;
			}
			
			sscanf(res_num,"%d",&(*residue)[FA->res_cnt].number);

			//printf("New residue: %d fatm[%d]=%d(%s-%s) latm[%d]=%d(%s-%s) :: %d\n", FA->res_cnt, FA->res_cnt, (*atoms)[(*residue)[FA->res_cnt].fatm[0]].number, (*residue)[FA->res_cnt].name, (*atoms)[(*residue)[FA->res_cnt].fatm[0]].name, FA->res_cnt-1, (*atoms)[(*residue)[FA->res_cnt-1].latm[0]].number, (*residue)[FA->res_cnt-1].name, (*atoms)[(*residue)[FA->res_cnt-1].latm[0]].name,(*residue)[FA->res_cnt].number);
            
		}

		(*atoms)[FA->atm_cnt].ofres=FA->res_cnt;

	}

	return;
}
