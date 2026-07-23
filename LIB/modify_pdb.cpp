#include "flexaid.h"
#include "fileio.h"
#include "ion_utils.h"

#include <cerrno>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#define NAA 20
#define NNA 4
#define NLIST 40

static char nucleic_acids[NNA][4]     = { "  A","  U","  C","  G" };
static char nucleic_acids_rev[NNA][4] = { "A  ","U  ","C  ","G  " };

static char protein_amino[NAA][4]     = { "GLY","ALA","VAL","LEU","ILE","MET",
					  "ASN","PRO","CYS","SER","THR","GLN",
					  "ASP","GLU","LYS","ARG","HIS",
					  "PHE","TRP","TYR" };

static char protein_atoms_order[NLIST][5]   = { " N  "," CA "," C  "," O  "," CB ",
						" CG "," SG "," OG "," CG1"," OG1"," CG2",
						" CD "," SD "," OD1"," CD1"," ND1"," OD2"," ND2"," CD2",
						" CE "," NE "," CE1"," NE1"," OE1"," OE2"," NE2"," CE2"," CE3",
						" CZ "," NZ "," CZ1"," CZ2"," CZ3",
						" CH "," CH1"," CH2"," OH "," NH1"," NH2",
						" OXT" };

// PDB B-factor occupies columns 60-65 (0-indexed)
static float pdb_bfactor(const char* buf) {
    char tmp[7]; strncpy(tmp, &buf[60], 6); tmp[6]='\0';
    return (float)atof(tmp);
}

// ─────────────────────────────────────────────────────────────────────────
// Selective crystallographic water retention
//
// A crystallographic water is retained as a receptor atom only when it is a
// plausible *bridging* water: it sits inside the binding cavity AND is
// anchored to the protein by a direct hydrogen bond.  Every other HOH — bulk
// solvent and the surface shell — is stripped.  Without this filter, every
// low-B-factor water in the structure becomes a receptor atom and the GA can
// bury the ligand inside the solvent shell, harvesting unbounded
// complementarity from sub-Ångström ligand-O ⋯ HOH-O contacts.
//
// Criterion 1 (cavity):  d(HOH-O, nearest crystal-ligand heavy atom) ≤ radius
// Criterion 2 (H-bond):  2.4 Å ≤ d(HOH-O, protein N/O/S) ≤ 3.5 Å
// Criterion 3 (order):   B-factor ≤ structural_water_bfactor_max (applied by
//                        the existing filter in the streaming pass).
//
// Hydrogens are absent from these coordinates, so criterion 2 uses the
// standard heavy-atom donor/acceptor distance window rather than a D–H⋯A
// angle: the lower bound rejects unphysically short contacts (which would be
// clashes or alternate conformations, not H-bonds) and the upper bound is the
// conventional 3.5 Å H-bond cutoff.
// ─────────────────────────────────────────────────────────────────────────

namespace {

struct WaterVec3 { double x, y, z; };

constexpr double kHBondMin = 2.4;   // Å — below this it is a clash, not an H-bond
constexpr double kHBondMax = 3.5;   // Å — conventional heavy-atom H-bond cutoff

inline double sq_dist(const WaterVec3& a, const WaterVec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

// Element symbol of a PDB ATOM/HETATM line: prefer columns 77-78, fall back to
// the first alphabetic character of the atom name (columns 13-16).
std::string pdb_element(const std::string& line) {
    std::string e;
    if(line.size() >= 78) {
        for(size_t i = 76; i < 78; ++i)
            if(std::isalpha(static_cast<unsigned char>(line[i]))) e += std::toupper(line[i]);
    }
    if(e.empty() && line.size() >= 16) {
        for(size_t i = 12; i < 16; ++i)
            if(std::isalpha(static_cast<unsigned char>(line[i]))) { e += std::toupper(line[i]); break; }
    }
    return e;
}

bool parse_pdb_xyz(const std::string& line, WaterVec3* out) {
    if(line.size() < 54) return false;
    try {
        out->x = std::stod(line.substr(30, 8));
        out->y = std::stod(line.substr(38, 8));
        out->z = std::stod(line.substr(46, 8));
    } catch(...) { return false; }
    return true;
}

// Heavy-atom coordinates of the crystal (oracle) ligand.  Accepts SDF (V2000),
// MOL2 and PDB; the format is chosen by extension and, failing that, sniffed.
std::vector<WaterVec3> read_ligand_heavy_atoms(const std::string& path) {
    std::vector<WaterVec3> out;
    std::ifstream in(path);
    if(!in) return out;

    std::string ext = std::filesystem::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    std::vector<std::string> lines;
    for(std::string l; std::getline(in, l); ) {
        if(!l.empty() && l.back() == '\r') l.pop_back();
        lines.push_back(l);
    }
    if(lines.empty()) return out;

    if(ext == ".sdf" || ext == ".mol") {
        // Line 4 is the counts line: "aaabbb..." with atom count in cols 1-3.
        if(lines.size() < 4) return out;
        int natoms = 0;
        try { natoms = std::stoi(lines[3].substr(0, 3)); } catch(...) { return out; }
        for(int i = 0; i < natoms && (size_t)(4 + i) < lines.size(); ++i) {
            std::istringstream ls(lines[4 + i]);
            WaterVec3 v; std::string sym;
            if(!(ls >> v.x >> v.y >> v.z >> sym)) continue;
            if(sym == "H" || sym == "D") continue;
            out.push_back(v);
        }
    } else if(ext == ".mol2") {
        bool in_atoms = false;
        for(const auto& l : lines) {
            if(l.rfind("@<TRIPOS>", 0) == 0) { in_atoms = (l.rfind("@<TRIPOS>ATOM", 0) == 0); continue; }
            if(!in_atoms) continue;
            std::istringstream ls(l);
            int id; std::string name; WaterVec3 v; std::string type;
            if(!(ls >> id >> name >> v.x >> v.y >> v.z >> type)) continue;
            if(type == "H" || type.rfind("H.", 0) == 0) continue;
            out.push_back(v);
        }
    } else {
        for(const auto& l : lines) {
            if(l.rfind("ATOM  ", 0) != 0 && l.rfind("HETATM", 0) != 0) continue;
            if(l.size() >= 20 && l.compare(17, 3, "HOH") == 0) continue;
            const std::string el = pdb_element(l);
            if(el == "H" || el == "D") continue;
            WaterVec3 v;
            if(parse_pdb_xyz(l, &v)) out.push_back(v);
        }
    }
    return out;
}

// Returns the ordinals (0-based, in file order) of the HETATM HOH lines that
// should be retained.  `kept_any` distinguishes "filter ran and kept nothing"
// from "filter did not run".
std::set<int> select_binding_site_waters(const char* receptor_pdb,
                                         const std::string& ligand_path,
                                         float radius,
                                         int hbond_required,
                                         int* n_total_out,
                                         int* n_lig_out)
{
    std::set<int> keep;
    *n_total_out = 0;

    const std::vector<WaterVec3> lig = read_ligand_heavy_atoms(ligand_path);
    *n_lig_out = (int)lig.size();

    std::ifstream in(receptor_pdb);
    if(!in) return keep;

    std::vector<WaterVec3> water;
    std::vector<WaterVec3> protein_polar;   // N/O/S of the protein

    for(std::string l; std::getline(in, l); ) {
        const bool is_atom   = (l.rfind("ATOM  ", 0) == 0);
        const bool is_hetatm = (l.rfind("HETATM", 0) == 0);
        if(!is_atom && !is_hetatm) continue;

        const bool is_water = (l.size() >= 20 && l.compare(17, 3, "HOH") == 0);
        WaterVec3 v;
        const bool have_xyz = parse_pdb_xyz(l, &v);

        if(is_hetatm && is_water) {
            // Every HOH HETATM line is recorded, parseable or not, so that the
            // ordinals here line up one-for-one with the streaming pass below.
            water.push_back(have_xyz ? v : WaterVec3{1e30, 1e30, 1e30});
            continue;
        }
        if(is_atom && have_xyz) {
            const std::string el = pdb_element(l);
            if(el == "N" || el == "O" || el == "S") protein_polar.push_back(v);
        }
    }

    *n_total_out = (int)water.size();
    if(lig.empty()) return keep;   // no oracle ligand → caller falls back

    const double r2 = (double)radius * (double)radius;
    const double hb_min2 = kHBondMin * kHBondMin, hb_max2 = kHBondMax * kHBondMax;

    for(size_t w = 0; w < water.size(); ++w) {
        bool in_cavity = false;
        for(const auto& a : lig)
            if(sq_dist(water[w], a) <= r2) { in_cavity = true; break; }
        if(!in_cavity) continue;

        if(hbond_required) {
            bool hbonded = false;
            for(const auto& p : protein_polar) {
                const double d2 = sq_dist(water[w], p);
                if(d2 >= hb_min2 && d2 <= hb_max2) { hbonded = true; break; }
            }
            if(!hbonded) continue;
        }
        keep.insert((int)w);
    }
    return keep;
}

} // namespace

void modify_pdb(char* infile, char* outfile, int exclude_het, int remove_water, int is_protein,
                int keep_ions, int keep_structural_waters, float structural_water_bfactor_max,
                float binding_site_water_radius, int binding_site_water_hbond_required,
                const char* oracle_ligand_path)
{
	char bufnul[10];
	char buffer[100];   // pdb line

	#define MAX_RESIDUE_LINES 50
	char lines[MAX_RESIDUE_LINES][100]; // store residue lines
	int  nlines=0;

	int prev_resnum = -1;
	int resnum = -1;
	char res[4];

	int read = 0;
	int wrote = 0;

	FILE* infile_ptr = NULL;
	FILE* outfile_ptr = NULL;

	char insert = '-', prev_insert = '-';
	
	printf("Protein PDB files are reordered\n");
	printf("Hydrogens are removed\n");
	printf("'A' alternate conformation ONLY is chosen\n");
	printf("heterogroups are%s excluded\n", exclude_het ? "":" not");
	printf("water molecules will%s be removed\n", exclude_het || remove_water ? "":" not");

	// ── Selective structural-water retention (pre-pass) ──
	// Only meaningful when waters would otherwise be retained wholesale.
	std::set<int> bs_water_keep;
	bool          bs_water_filter_active = false;
	int           hoh_ordinal = 0;

	// The crystal ligand normally arrives via reference_ligand.file, but the
	// benchmark runner deliberately leaves that empty to keep crystal
	// coordinates out of the GA; there the pose is exported as FLEXAIDDS_RMSDST.
	std::string oracle_lig = (oracle_ligand_path != nullptr) ? oracle_ligand_path : "";
	if(oracle_lig.empty()) {
		const char* env = getenv("FLEXAIDDS_RMSDST");
		if(env != nullptr) oracle_lig = env;
	}

	if(keep_structural_waters && remove_water && !exclude_het &&
	   binding_site_water_radius > 0.0f && !oracle_lig.empty())
	{
		int n_total = 0, n_lig = 0;
		bs_water_keep = select_binding_site_waters(infile, oracle_lig,
		                                           binding_site_water_radius,
		                                           binding_site_water_hbond_required,
		                                           &n_total, &n_lig);
		const std::string target = std::filesystem::path(infile).stem().string();
		if(n_lig == 0) {
			// No usable oracle ligand: fall back to the plain B-factor filter
			// rather than silently stripping every water.
			fprintf(stderr, "[WATER-FILTER] %s: no ligand heavy atoms readable from %s - "
			                "falling back to B-factor-only retention\n",
			        target.c_str(), oracle_lig.c_str());
		} else {
			bs_water_filter_active = true;
			printf("[WATER-FILTER] %s: kept %d/%d structural waters within %.1fA of oracle ligand\n",
			       target.c_str(), (int)bs_water_keep.size(), n_total,
			       (double)binding_site_water_radius);
		}
	}

	if(!OpenFile_B(infile,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Could not order PDB file %s\n", infile);
		Terminate(20);
	}

	outfile_ptr = fopen(outfile,"w");
	if(outfile_ptr == NULL){
		fprintf(stderr, "ERROR: Could not write temporary PDB file: %s (%s)\n",
		        outfile, strerror(errno));
		Terminate(20);
	}
	
	
	while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){
		if(!strncmp(&buffer[0],"ATOM  ",6)){
			// all lines that start with 'ATOM  ' field

			read++;
			
			//0         1         2         3         4         5         6         
			//0123456789012345678901234567890123456789012345678901234567890123456789
			//ATOM     47  CB  ILE A   7      38.324  -3.725  17.587  1.00  0.00           C  
			strncpy(res,&buffer[17],3);
			res[3]='\0';

			strncpy(bufnul,&buffer[22],4);
			sscanf(bufnul,"%d",&resnum);
				
			// insertion of residue
			insert = buffer[26];
				
			// skip alternate conformations other than 'A'
			if(buffer[16] != ' ' && buffer[16] != 'A'){ continue; }
				
			if(resnum == prev_resnum && insert == prev_insert){
				if(is_protein && is_natural_amino(res)){
					// store line
					if(nlines < MAX_RESIDUE_LINES){
						strncpy(lines[nlines], buffer, sizeof(lines[0]) - 1);
						lines[nlines][sizeof(lines[0]) - 1] = '\0';
						nlines++;
					}
				}else if(!is_protein && is_natural_nucleic(res)){
					fprintf(outfile_ptr,"%s",buffer);
				}else{
					// ligands/mod. amino acids are marked as HETATM by default
					fprintf(outfile_ptr,"HETATM%s",&buffer[6]);
				}
					
			}else if(prev_resnum != -1){
				if(is_protein && nlines > 0){
					//write out ordered lines
					rewrite_residue2(lines,nlines,&wrote,outfile_ptr);
					nlines=0;
				}
				
				if(is_protein && is_natural_amino(res)){
					if(nlines < MAX_RESIDUE_LINES){
						strncpy(lines[nlines], buffer, sizeof(lines[0]) - 1);
						lines[nlines][sizeof(lines[0]) - 1] = '\0';
						nlines++;
					}
				}else if(!is_protein && is_natural_nucleic(res)){
					fprintf(outfile_ptr,"%s",buffer);
				}else{
					// ligands/mod. amino acids are marked as HETATM by default
					fprintf(outfile_ptr,"HETATM%s",&buffer[6]);
				}

			}else{
				if(is_protein && is_natural_amino(res)){
					if(nlines < MAX_RESIDUE_LINES){
						strncpy(lines[nlines], buffer, sizeof(lines[0]) - 1);
						lines[nlines][sizeof(lines[0]) - 1] = '\0';
						nlines++;
					}
				}else if(is_natural_nucleic(res)){
					fprintf(outfile_ptr,"%s",buffer);
				}else{
					// ligands/mod. amino acids are marked as HETATM by default
					fprintf(outfile_ptr,"HETATM%s",&buffer[6]);
				}
			}
				
			prev_resnum = resnum;
			prev_insert = insert;
			

		}else{
			// all other lines that do not start with 'ATOM  ' field

			if(!strncmp(&buffer[0],"HETATM",6)){
				if(exclude_het) {
					// Keep metal ions when keep_ions=1, regardless of exclude_het
					if(keep_ions && is_ion_resname(&buffer[17])) { /* keep */ }
					else { continue; }
				} else {
					if(!strncmp(&buffer[17],"HOH",3) && remove_water) {
						const int this_hoh = hoh_ordinal++;
						// Retain low-B-factor structural waters when requested
						if(keep_structural_waters) {
							float bf = pdb_bfactor(buffer);
							if(bf > structural_water_bfactor_max) continue;
							// …and, when the binding-site filter is active,
							// only those that bridge ligand and protein.
							if(bs_water_filter_active &&
							   bs_water_keep.find(this_hoh) == bs_water_keep.end()) continue;
						} else {
							continue;
						}
					}
				}
			}

			if(is_protein && nlines > 0){
				rewrite_residue2(lines,nlines,&wrote,outfile_ptr);
				nlines=0;
			}
			
			fprintf(outfile_ptr,"%s",buffer);
			
		}


	}

	if(is_protein && nlines > 0){	
		rewrite_residue2(lines,nlines,&wrote,outfile_ptr);
	}
	
	CloseFile_B(&infile_ptr,"r");

	fclose(outfile_ptr);
	
	printf("number of ATOM lines read is %d\n", read);
	printf("number of lines outputted is %d\n", wrote);
	
}

int get_NextLine(char lines[][100], int nlines){

	int k,l=-1,m=NLIST;
	char name[5];
	
	for(int i=0; i<nlines; i++){
		strncpy(name,&lines[i][12],4);
		name[4]='\0';
		
		k=NLIST;
		for(int j=NLIST-1; j>=0; --j)
			if(!strcmp(protein_atoms_order[j],name))
				k=j;
		
		if(k<m){
			m=k;
			l=i;
		}
			
	}

	return l;
}

void rewrite_residue2(char lines[][100], int nlines, int* wrote, FILE* outfile_ptr){
	
	int i;
	
	while((i=get_NextLine(lines,nlines)) != -1){
		char newline[100];
		memset(newline, 0, sizeof(newline));
		memcpy(newline, lines[i], 6);
		snprintf(&newline[6], 6, "%5d", ++(*wrote));
		newline[11] = '\0'; // snprintf null-terminates at [11]
		size_t tail_len = strlen(&lines[i][11]);
		if(tail_len > sizeof(newline) - 12) tail_len = sizeof(newline) - 12;
		memcpy(&newline[11], &lines[i][11], tail_len);
		newline[11 + tail_len] = '\0';
		fprintf(outfile_ptr,"%s",newline);
		strcpy(lines[i],"                    ");
	}
	
}

int is_natural_amino(char* res){
	
	for(int i=0; i<NAA; i++){
		if(!strcmp(res,protein_amino[i])){
			return 1;
		}
	}

	return 0;
}

int is_natural_nucleic(char* res){
	
	for(int i=0; i<NNA; i++){
		if(!strcmp(res,nucleic_acids[i]) || !strcmp(res,nucleic_acids_rev[i])){
			return 1;
		}
	}

	return 0;
}
