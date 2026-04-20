// CifReader.cpp — PDBx/mmCIF parser for FlexAIDdS
//
// Parses the _atom_site loop from mmCIF files. Column order is auto-detected
// from the header, so this works with any compliant mmCIF file regardless
// of column arrangement.
//
// Key _atom_site columns used:
//   group_PDB      → ATOM / HETATM
//   type_symbol     → element (C, N, O, ...)
//   label_atom_id   → atom name (CA, CB, ...)
//   label_comp_id   → residue name (ALA, LIG, ...)
//   label_asym_id   → chain ID
//   label_seq_id    → residue sequence number
//   Cartn_x/y/z    → coordinates
//   occupancy       → occupancy
//   B_iso_or_equiv  → B-factor
//   auth_seq_id     → author residue number (preferred for FlexAID)
//   auth_comp_id    → author residue name
//   auth_asym_id    → author chain ID
//   auth_atom_id    → author atom name
//   pdbx_PDB_ins_code → insertion code
//   pdbx_formal_charge → formal charge
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.

#include "CifReader.h"
#include "fileio.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>

// Column indices in the _atom_site loop (set to -1 if not present)
struct CifAtomSiteColumns {
    int group_PDB       = -1;  // ATOM / HETATM
    int id              = -1;  // atom serial number
    int type_symbol     = -1;  // element symbol
    int label_atom_id   = -1;  // atom name
    int label_alt_id    = -1;  // alt location
    int label_comp_id   = -1;  // residue name
    int label_asym_id   = -1;  // chain ID
    int label_seq_id    = -1;  // residue number
    int Cartn_x         = -1;
    int Cartn_y         = -1;
    int Cartn_z         = -1;
    int occupancy       = -1;
    int B_iso_or_equiv  = -1;
    int auth_seq_id     = -1;  // author residue number
    int auth_comp_id    = -1;  // author residue name
    int auth_asym_id    = -1;  // author chain ID
    int auth_atom_id    = -1;  // author atom name
    int pdbx_PDB_ins_code    = -1;
    int pdbx_formal_charge   = -1;
    int pdbx_PDB_model_num   = -1;

    int num_columns = 0;
};

// Tokenize a whitespace-delimited line, respecting single/double quotes
static std::vector<std::string> tokenize_cif_line(const std::string& line) {
    std::vector<std::string> tokens;
    size_t i = 0;
    size_t n = line.size();
    while (i < n) {
        // Skip whitespace
        while (i < n && std::isspace(static_cast<unsigned char>(line[i]))) i++;
        if (i >= n) break;

        std::string tok;
        if (line[i] == '\'' || line[i] == '"') {
            // Quoted value
            char quote = line[i++];
            while (i < n && line[i] != quote) tok += line[i++];
            if (i < n) i++; // skip closing quote
        } else {
            // Unquoted value
            while (i < n && !std::isspace(static_cast<unsigned char>(line[i])))
                tok += line[i++];
        }
        tokens.push_back(tok);
    }
    return tokens;
}

// Map column name to index in CifAtomSiteColumns
static void map_column(CifAtomSiteColumns& cols, const std::string& name, int idx) {
    // Strip "_atom_site." prefix
    std::string attr = name;
    if (attr.find("_atom_site.") == 0)
        attr = attr.substr(11);

    if (attr == "group_PDB")          cols.group_PDB = idx;
    else if (attr == "id")            cols.id = idx;
    else if (attr == "type_symbol")   cols.type_symbol = idx;
    else if (attr == "label_atom_id") cols.label_atom_id = idx;
    else if (attr == "label_alt_id")  cols.label_alt_id = idx;
    else if (attr == "label_comp_id") cols.label_comp_id = idx;
    else if (attr == "label_asym_id") cols.label_asym_id = idx;
    else if (attr == "label_seq_id")  cols.label_seq_id = idx;
    else if (attr == "Cartn_x")       cols.Cartn_x = idx;
    else if (attr == "Cartn_y")       cols.Cartn_y = idx;
    else if (attr == "Cartn_z")       cols.Cartn_z = idx;
    else if (attr == "occupancy")     cols.occupancy = idx;
    else if (attr == "B_iso_or_equiv") cols.B_iso_or_equiv = idx;
    else if (attr == "auth_seq_id")   cols.auth_seq_id = idx;
    else if (attr == "auth_comp_id")  cols.auth_comp_id = idx;
    else if (attr == "auth_asym_id")  cols.auth_asym_id = idx;
    else if (attr == "auth_atom_id")  cols.auth_atom_id = idx;
    else if (attr == "pdbx_PDB_ins_code")    cols.pdbx_PDB_ins_code = idx;
    else if (attr == "pdbx_formal_charge")   cols.pdbx_formal_charge = idx;
    else if (attr == "pdbx_PDB_model_num")   cols.pdbx_PDB_model_num = idx;
}

// Get string value from tokens, return "?" if column not mapped
static std::string get_val(const std::vector<std::string>& toks, int col_idx) {
    if (col_idx < 0 || col_idx >= static_cast<int>(toks.size())) return "?";
    return toks[col_idx];
}

// Get float value, return 0.0 if missing
static float get_float(const std::vector<std::string>& toks, int col_idx) {
    std::string v = get_val(toks, col_idx);
    if (v == "?" || v == ".") return 0.0f;
    return static_cast<float>(std::atof(v.c_str()));
}

// Get int value, return 0 if missing
static int get_int(const std::vector<std::string>& toks, int col_idx) {
    std::string v = get_val(toks, col_idx);
    if (v == "?" || v == ".") return 0;
    return std::atoi(v.c_str());
}

// Element symbol → default FlexAID type (same mapping as SdfReader)
#if 0  // currently unused — kept for future CIF element→type mapping
static int element_to_type(const char* elem) {
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
    if (!strcmp(elem, "Na")) return 34;
    if (!strcmp(elem, "K"))  return 35;
    return 39; // unknown → dummy
}
#endif  // element_to_type

#if 0  // currently unused — kept for future CIF element→radius mapping
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
#endif  // element_radius

// Core CIF parser — reads _atom_site loop, filters by group (ATOM/HETATM/both)
// and populates FlexAID atom/resid arrays.
static int parse_cif_atom_site(
    FA_Global* FA, atom** atoms_ptr, resid** residue_ptr,
    const char* cif_file,
    bool read_atom,     // include ATOM records
    bool read_hetatm,   // include HETATM records
    int  model_num)     // which model (1 = first; 0 = all)
{
    FILE* fp = fopen(cif_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open CIF file: %s\n", cif_file);
        return 0;
    }

    // Phase 1: Find the _atom_site loop and parse column headers
    CifAtomSiteColumns cols;
    char buf[2048];
    bool in_atom_site_loop = false;
    bool reading_headers = false;
    int col_count = 0;

    // Scan for "loop_" followed by "_atom_site.*" headers
    while (fgets(buf, sizeof(buf), fp)) {
        std::string line(buf);
        // Trim trailing whitespace
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();

        if (line == "loop_") {
            in_atom_site_loop = false;
            reading_headers = true;
            col_count = 0;
            continue;
        }

        if (reading_headers) {
            if (line.find("_atom_site.") == 0) {
                in_atom_site_loop = true;
                map_column(cols, line, col_count);
                col_count++;
                continue;
            } else if (in_atom_site_loop) {
                // First non-header line after _atom_site headers = data
                reading_headers = false;
                cols.num_columns = col_count;
                // Don't consume this line — process it as data below
                break;
            } else {
                // Not an _atom_site loop — reset
                reading_headers = false;
                continue;
            }
        }
    }

    if (!in_atom_site_loop || cols.num_columns == 0) {
        fprintf(stderr, "ERROR: No _atom_site loop found in CIF file: %s\n", cif_file);
        fclose(fp);
        return 0;
    }

    // Validate required columns
    if (cols.Cartn_x < 0 || cols.Cartn_y < 0 || cols.Cartn_z < 0) {
        fprintf(stderr, "ERROR: CIF file missing Cartn_x/y/z columns: %s\n", cif_file);
        fclose(fp);
        return 0;
    }

    // Phase 2: Parse data rows
    // We already have the first data line in buf from the header scan
    // Re-read from the position after headers

    // Collect all atom data first, then build arrays
    struct CifAtom {
        std::string group;    // ATOM or HETATM
        int    serial;
        std::string element;
        std::string atom_name;
        char   alt_loc;
        std::string res_name;
        std::string chain;
        int    res_seq;
        char   ins_code;
        float  x, y, z;
        float  occupancy;
        float  bfactor;
        int    charge;
        int    model;
    };

    std::vector<CifAtom> cif_atoms;

    // Process the current line (first data line) and continue
    auto process_line = [&](const std::string& line) {
        if (line.empty() || line[0] == '#' || line[0] == '_' || line == "loop_")
            return false;  // end of data block

        auto toks = tokenize_cif_line(line);
        if (toks.empty() || static_cast<int>(toks.size()) < cols.num_columns)
            return false;  // malformed or end of block

        // Check if this is still atom data (first token should be ATOM or HETATM)
        std::string group = get_val(toks, cols.group_PDB);
        if (group != "ATOM" && group != "HETATM")
            return false;

        // Filter by group
        if (group == "ATOM" && !read_atom) return true;
        if (group == "HETATM" && !read_hetatm) return true;

        // Filter by model number
        int mdl = get_int(toks, cols.pdbx_PDB_model_num);
        if (model_num > 0 && mdl > 0 && mdl != model_num) return true;

        // Filter alt locations — keep '.' and 'A' only
        std::string alt = get_val(toks, cols.label_alt_id);
        if (alt != "." && alt != "?" && alt != "A") return true;

        CifAtom a;
        a.group = group;
        a.serial = get_int(toks, cols.id);

        // Prefer auth_* fields over label_* (matches PDB convention)
        std::string elem = get_val(toks, cols.type_symbol);
        a.element = (elem != "?" && elem != ".") ? elem : "X";

        std::string aname = get_val(toks, cols.auth_atom_id);
        if (aname == "?" || aname == ".") aname = get_val(toks, cols.label_atom_id);
        a.atom_name = aname;

        a.alt_loc = (alt == "." || alt == "?") ? ' ' : alt[0];

        std::string rname = get_val(toks, cols.auth_comp_id);
        if (rname == "?" || rname == ".") rname = get_val(toks, cols.label_comp_id);
        a.res_name = rname;

        std::string chn = get_val(toks, cols.auth_asym_id);
        if (chn == "?" || chn == ".") chn = get_val(toks, cols.label_asym_id);
        a.chain = chn;

        int rseq = get_int(toks, cols.auth_seq_id);
        if (rseq == 0) rseq = get_int(toks, cols.label_seq_id);
        a.res_seq = rseq;

        std::string ins = get_val(toks, cols.pdbx_PDB_ins_code);
        a.ins_code = (ins == "?" || ins == ".") ? ' ' : ins[0];

        a.x = get_float(toks, cols.Cartn_x);
        a.y = get_float(toks, cols.Cartn_y);
        a.z = get_float(toks, cols.Cartn_z);
        a.occupancy = get_float(toks, cols.occupancy);
        a.bfactor   = get_float(toks, cols.B_iso_or_equiv);
        a.charge    = get_int(toks, cols.pdbx_formal_charge);
        a.model     = mdl;

        cif_atoms.push_back(a);
        return true;
    };

    // Process the first data line already in buf
    {
        std::string line(buf);
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        process_line(line);
    }

    // Continue reading
    while (fgets(buf, sizeof(buf), fp)) {
        std::string line(buf);
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        if (!process_line(line)) break;
    }

    fclose(fp);

    if (cif_atoms.empty()) {
        fprintf(stderr, "ERROR: No atoms found in CIF file: %s\n", cif_file);
        return 0;
    }

    printf("CIF: parsed %zu atoms from %s\n", cif_atoms.size(), cif_file);

    // Phase 3: Convert to PDB-style temp file and let read_pdb handle it
    // This reuses the existing PDB reading infrastructure which handles
    // residue connectivity, type assignment, and all the FlexAID-specific
    // bookkeeping. Writing a temp PDB is the most robust approach.
    char tmp_pdb[MAX_PATH__];
    snprintf(tmp_pdb, MAX_PATH__, "/tmp/flexaid_cif_%d.pdb", static_cast<int>(std::random_device{}() % 900000 + 100000));

    FILE* out = fopen(tmp_pdb, "w");
    if (!out) {
        fprintf(stderr, "ERROR: Cannot create temp PDB for CIF conversion\n");
        return 0;
    }

    for (const auto& a : cif_atoms) {
        // Format as PDB ATOM/HETATM record (fixed-width 80 chars)
        // Columns: 1-6 record, 7-11 serial, 13-16 name, 17 altloc,
        //          18-20 resname, 22 chain, 23-26 resseq, 27 icode,
        //          31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 bfac,
        //          77-78 element

        // Atom name formatting: if 4 chars use cols 13-16, else right-justify in 13-16
        char name_buf[5] = "    ";
        if (a.atom_name.size() >= 4) {
            for (int k = 0; k < 4 && k < (int)a.atom_name.size(); k++)
                name_buf[k] = a.atom_name[k];
        } else {
            // Right-justify: leading space for 1-3 char names
            int offset = (a.element.size() > 1) ? 0 : 1;
            for (int k = 0; k < (int)a.atom_name.size(); k++)
                name_buf[offset + k] = a.atom_name[k];
        }

        char chain_c = a.chain.empty() ? ' ' : a.chain[0];
        char ins_c   = a.ins_code;

        // Truncate residue name to 3 chars
        char rname[4] = "   ";
        for (int k = 0; k < 3 && k < (int)a.res_name.size(); k++)
            rname[k] = a.res_name[k];

        fprintf(out, "%-6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                a.group.c_str(),
                a.serial % 100000,
                name_buf,
                a.alt_loc,
                rname,
                chain_c,
                a.res_seq % 10000,
                ins_c,
                a.x, a.y, a.z,
                a.occupancy > 0.0f ? a.occupancy : 1.0f,
                a.bfactor,
                a.element.c_str());
    }
    fprintf(out, "END\n");
    fclose(out);

    // Now use the existing PDB reader on the temp file
    int result;
    if (read_atom && !read_hetatm) {
        // Receptor mode — use read_pdb
        read_pdb(FA, atoms_ptr, residue_ptr, tmp_pdb);
        result = (FA->atm_cnt > 0) ? 1 : 0;
    } else {
        // Ligand mode — need to handle differently
        // For now, treat as MOL2-style ligand reading via the PDB temp
        read_pdb(FA, atoms_ptr, residue_ptr, tmp_pdb);
        result = (FA->atm_cnt > 0) ? 1 : 0;
    }

    remove(tmp_pdb);
    return result;
}

// Public API: read receptor (ATOM records)
int read_cif_receptor(FA_Global* FA, atom** atoms, resid** residue,
                      const char* cif_file) {
    return parse_cif_atom_site(FA, atoms, residue, cif_file,
                               /*read_atom=*/true, /*read_hetatm=*/false,
                               /*model_num=*/1);
}

// Public API: read ligand (HETATM records)
int read_cif_ligand(FA_Global* FA, atom** atoms, resid** residue,
                    const char* cif_file) {
    return parse_cif_atom_site(FA, atoms, residue, cif_file,
                               /*read_atom=*/false, /*read_hetatm=*/true,
                               /*model_num=*/1);
}


// =============================================================================
// CCBM: Multi-model PDB reader
// =============================================================================
// Reads a PDB file with MODEL/ENDMDL records and extracts coordinates for
// each model into FA->model_coords. The first model is loaded fully via
// read_pdb(), then subsequent models have their coordinates extracted.

int read_multi_model_pdb(FA_Global* FA, atom** atoms, resid** residue,
                         const char* pdb_file) {
    // Phase 1: Count models and collect per-model ATOM coordinate lines
    FILE* fp = fopen(pdb_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open PDB file for multi-model: %s\n", pdb_file);
        return 0;
    }

    // Collect all coordinates per model
    struct ModelCoord {
        int    model_num;
        std::vector<std::array<float,3>> coords;
    };
    std::vector<ModelCoord> models;
    int current_model = 1;  // default model if no MODEL records
    bool has_model_records = false;
    [[maybe_unused]] bool first_model_started = false;

    char buf[256];
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "MODEL", 5) == 0) {
            has_model_records = true;
            sscanf(buf + 5, "%d", &current_model);
            models.push_back({current_model, {}});
            first_model_started = true;
            continue;
        }
        if (strncmp(buf, "ENDMDL", 6) == 0) {
            continue;
        }
        if (strncmp(buf, "ATOM", 4) == 0 || strncmp(buf, "HETATM", 6) == 0) {
            // Parse coordinates from columns 31-54 (PDB standard)
            if (strlen(buf) < 54) continue;
            float x, y, z;
            char xbuf[9], ybuf[9], zbuf[9];
            strncpy(xbuf, buf + 30, 8); xbuf[8] = '\0';
            strncpy(ybuf, buf + 38, 8); ybuf[8] = '\0';
            strncpy(zbuf, buf + 46, 8); zbuf[8] = '\0';
            x = static_cast<float>(atof(xbuf));
            y = static_cast<float>(atof(ybuf));
            z = static_cast<float>(atof(zbuf));

            if (!has_model_records && models.empty()) {
                models.push_back({1, {}});
            }
            if (!models.empty()) {
                models.back().coords.push_back({x, y, z});
            }
        }
    }
    fclose(fp);

    if (models.empty()) {
        // No atoms found — fall back to single-model
        fprintf(stderr, "WARNING: No ATOM records found in %s, attempting single-model read\n", pdb_file);
        FA->n_models = 1;
        FA->model_coords.clear();
        FA->model_strain.clear();
        FA->model_strain.push_back(0.0);
        return 1;
    }

    // Phase 2: Load first model via standard read_pdb for full topology
    // Create a temp PDB with only the first model's atoms
    char tmp_pdb[MAX_PATH__];
    snprintf(tmp_pdb, MAX_PATH__, "/tmp/flexaid_mm_%d.pdb", static_cast<int>(std::random_device{}() % 900000 + 100000));
    FILE* out = fopen(tmp_pdb, "w");
    if (!out) {
        fprintf(stderr, "ERROR: Cannot create temp PDB for multi-model\n");
        return 0;
    }

    // Re-read and write only the first model
    fp = fopen(pdb_file, "r");
    if (!fp) { fclose(out); return 0; }

    int write_model = -1;
    bool writing = !has_model_records;  // if no MODEL records, write everything
    while (fgets(buf, sizeof(buf), fp)) {
        if (strncmp(buf, "MODEL", 5) == 0) {
            int m = 0;
            sscanf(buf + 5, "%d", &m);
            if (write_model < 0) {
                write_model = m;
                writing = true;
            } else {
                writing = false;
            }
            continue;
        }
        if (strncmp(buf, "ENDMDL", 6) == 0) {
            if (writing) {
                writing = false;
            }
            continue;
        }
        if (writing) {
            fputs(buf, out);
        }
    }
    fprintf(out, "END\n");
    fclose(fp);
    fclose(out);

    // Load the first model with full FlexAID infrastructure
    read_pdb(FA, atoms, residue, tmp_pdb);
    remove(tmp_pdb);

    // Phase 3: Build model_coords arrays
    int n_atoms_ref = static_cast<int>(models[0].coords.size());
    FA->n_models = static_cast<int>(models.size());
    FA->model_coords.resize(FA->n_models);
    FA->model_strain.resize(FA->n_models, 0.0);  // default: all strains = 0

    for (int m = 0; m < FA->n_models; ++m) {
        int n_atoms_m = static_cast<int>(models[m].coords.size());
        int n_copy = std::min(n_atoms_ref, n_atoms_m);
        FA->model_coords[m].resize(n_copy * 3, 0.0f);
        for (int a = 0; a < n_copy; ++a) {
            FA->model_coords[m][a * 3 + 0] = models[m].coords[a][0];
            FA->model_coords[m][a * 3 + 1] = models[m].coords[a][1];
            FA->model_coords[m][a * 3 + 2] = models[m].coords[a][2];
        }
    }

    printf("CCBM: parsed %d models from %s (%d atoms each)\n",
           FA->n_models, pdb_file, n_atoms_ref);

    return FA->n_models;
}


// =============================================================================
// CCBM: Multi-model CIF reader
// =============================================================================
// Reads a CIF file and groups atoms by pdbx_PDB_model_num into separate
// coordinate arrays in FA->model_coords.

int read_multi_model_cif(FA_Global* FA, atom** atoms, resid** residue,
                         const char* cif_file) {
    // Phase 1: Parse all atoms with their model numbers
    FILE* fp = fopen(cif_file, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open CIF file for multi-model: %s\n", cif_file);
        return 0;
    }

    CifAtomSiteColumns cols;
    char buf[2048];
    bool in_atom_site_loop = false;
    bool reading_headers = false;
    int col_count = 0;

    while (fgets(buf, sizeof(buf), fp)) {
        std::string line(buf);
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();

        if (line == "loop_") {
            in_atom_site_loop = false;
            reading_headers = true;
            col_count = 0;
            continue;
        }

        if (reading_headers) {
            if (line.find("_atom_site.") == 0) {
                in_atom_site_loop = true;
                map_column(cols, line, col_count);
                col_count++;
                continue;
            } else if (in_atom_site_loop) {
                reading_headers = false;
                cols.num_columns = col_count;
                break;
            } else {
                reading_headers = false;
                continue;
            }
        }
    }

    if (!in_atom_site_loop || cols.num_columns == 0) {
        fprintf(stderr, "ERROR: No _atom_site loop in CIF: %s\n", cif_file);
        fclose(fp);
        return 0;
    }

    // Collect atoms grouped by model number
    struct CifModelAtom {
        int    model;
        float  x, y, z;
    };
    std::vector<CifModelAtom> all_atoms;
    std::map<int, int> model_map;  // model_num -> sequential index

    auto process = [&](const std::string& line) -> bool {
        if (line.empty() || line[0] == '#' || line[0] == '_' || line == "loop_")
            return false;
        auto toks = tokenize_cif_line(line);
        if (toks.empty() || static_cast<int>(toks.size()) < cols.num_columns)
            return false;
        std::string group = get_val(toks, cols.group_PDB);
        if (group != "ATOM" && group != "HETATM")
            return false;

        // Filter alt locations
        std::string alt = get_val(toks, cols.label_alt_id);
        if (alt != "." && alt != "?" && alt != "A") return true;

        CifModelAtom a;
        a.model = get_int(toks, cols.pdbx_PDB_model_num);
        if (a.model <= 0) a.model = 1;
        a.x = get_float(toks, cols.Cartn_x);
        a.y = get_float(toks, cols.Cartn_y);
        a.z = get_float(toks, cols.Cartn_z);
        all_atoms.push_back(a);

        if (model_map.find(a.model) == model_map.end()) {
            int idx = static_cast<int>(model_map.size());
            model_map[a.model] = idx;
        }
        return true;
    };

    // Process first data line in buf
    {
        std::string line(buf);
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        process(line);
    }

    while (fgets(buf, sizeof(buf), fp)) {
        std::string line(buf);
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        if (!process(line)) break;
    }
    fclose(fp);

    int n_models = static_cast<int>(model_map.size());
    if (n_models == 0) {
        fprintf(stderr, "WARNING: No atoms found in CIF multi-model: %s\n", cif_file);
        FA->n_models = 1;
        return 0;
    }

    // Load first model via standard reader for topology
    [[maybe_unused]] int first_model = model_map.begin()->first;
    read_cif_receptor(FA, atoms, residue, cif_file);

    // Group coordinates by model
    std::vector<std::vector<std::array<float,3>>> model_coords(n_models);
    for (const auto& a : all_atoms) {
        int idx = model_map[a.model];
        model_coords[idx].push_back({a.x, a.y, a.z});
    }

    // Build FA->model_coords
    int n_atoms_ref = model_coords.empty() ? 0 : static_cast<int>(model_coords[0].size());
    FA->n_models = n_models;
    FA->model_coords.resize(n_models);
    FA->model_strain.resize(n_models, 0.0);

    for (int m = 0; m < n_models; ++m) {
        int n_copy = std::min(n_atoms_ref, static_cast<int>(model_coords[m].size()));
        FA->model_coords[m].resize(n_copy * 3, 0.0f);
        for (int a = 0; a < n_copy; ++a) {
            FA->model_coords[m][a * 3 + 0] = model_coords[m][a][0];
            FA->model_coords[m][a * 3 + 1] = model_coords[m][a][1];
            FA->model_coords[m][a * 3 + 2] = model_coords[m][a][2];
        }
    }

    printf("CCBM: parsed %d models from CIF %s (%d atoms each)\n",
           n_models, cif_file, n_atoms_ref);

    return n_models;
}
