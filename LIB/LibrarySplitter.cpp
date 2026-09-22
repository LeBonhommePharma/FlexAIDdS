// LibrarySplitter.cpp — Multi-molecule library + multi-model receptor splitting
//
// Handles all input library formats:
//   - Multi-molecule SDF ($$$$ separators)
//   - SMILES files (.smi, .smiles — one per line)
//   - Directories of MOL2/SDF/PDB files
//   - Multi-model PDB (MODEL/ENDMDL) for NMR ensembles, cryo-EM conformers
//   - Multi-model CIF (pdbx_PDB_model_num) for mmCIF ensembles
//   - Single files (passthrough)
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.

#include "LibrarySplitter.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <algorithm>
#include <cctype>
#include "temp_dir.h"
#include "fs_safe.h"

namespace library {

namespace fs = std::filesystem;

// Count $$$$ separators in an SDF file
static int count_sdf_molecules(const std::string& path) {
    std::ifstream in(path);
    if (!in) return 0;
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() >= 4 && line.substr(0, 4) == "$$$$")
            count++;
    }
    return std::max(count, 1); // at least 1 if file non-empty
}

// Count MODEL records in a PDB file (NMR/cryo-EM ensembles)
static int count_pdb_models(const std::string& path) {
    std::ifstream in(path);
    if (!in) return 1;
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() >= 5 && line.substr(0, 5) == "MODEL")
            count++;
    }
    return std::max(count, 1);
}

// Count non-empty lines in a SMILES file
static int count_smiles_lines(const std::string& path) {
    std::ifstream in(path);
    if (!in) return 0;
    int count = 0;
    std::string line;
    while (std::getline(in, line)) {
        // Trim
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        // Skip empty lines and comments
        if (!line.empty() && line[0] != '#')
            count++;
    }
    return count;
}

// Get lowercase extension
static std::string get_ext(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return "";
    std::string ext = path.substr(dot);
    for (auto& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext;
}

int detect_library_size(const std::string& path) {
    if (fs::is_directory(path)) {
        int count = 0;
        for (const auto& entry : fs::directory_iterator(path)) {
            std::string ext = get_ext(entry.path().string());
            if (ext == ".mol2" || ext == ".sdf" || ext == ".mol" || ext == ".pdb")
                count++;
        }
        return count;
    }

    std::string ext = get_ext(path);

    if (ext == ".sdf" || ext == ".mol")
        return count_sdf_molecules(path);

    if (ext == ".smi" || ext == ".smiles" || ext == ".txt")
        return count_smiles_lines(path);

    if (ext == ".pdb" || ext == ".ent")
        return count_pdb_models(path);

    return 1;
}

LibraryInfo split_library(const std::string& path) {
    LibraryInfo lib;
    lib.total = 0;

    // === Directory of files ===
    if (fs::is_directory(path)) {
        std::vector<std::string> files;
        for (const auto& entry : fs::directory_iterator(path)) {
            std::string ext = get_ext(entry.path().string());
            if (ext == ".mol2" || ext == ".sdf" || ext == ".mol" || ext == ".pdb")
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());

        for (const auto& f : files) {
            LigandEntry e;
            e.path = f;
            e.name = fs::path(f).stem().string();
            e.format = get_ext(f).substr(1); // remove dot
            e.is_temp = false;
            lib.ligands.push_back(e);
        }
        lib.total = static_cast<int>(lib.ligands.size());
        return lib;
    }

    std::string ext = get_ext(path);

    // === Multi-molecule SDF ===
    if (ext == ".sdf" || ext == ".mol") {
        int n_mols = count_sdf_molecules(path);

        if (n_mols <= 1) {
            // Single molecule — passthrough
            LigandEntry e;
            e.path = path;
            e.name = fs::path(path).stem().string();
            e.format = "sdf";
            e.is_temp = false;
            lib.ligands.push_back(e);
            lib.total = 1;
            return lib;
        }

        // Multi-molecule — split at $$$$
        lib.temp_dir = flexaids::temp_dir() + "/flexaid_lib_" +
                       std::to_string(std::random_device{}() % 900000 + 100000);
        // Sits directly downstream of the resolver that aborted nine cells on
        // 2026-09-20. If temp_dir() degraded to "." in a read-only cwd this
        // used to abort; if it merely failed and we carried on, every
        // LigandEntry below would name a file that was never written and the
        // campaign would dock a silently short library. Refuse loudly and
        // return an empty LibraryInfo -- callers size their result vector
        // from lib.total, so zero surfaces as "no ligands", not as success.
        {
            std::string why;
            if (!flexaids::fs_safe::ensure_dir(lib.temp_dir, &why)) {
                std::fprintf(stderr,
                    "[SPLIT] ERROR: cannot create SDF split directory '%s' (%s). "
                    "Set FLEXAIDDS_TMPDIR to a writable, durable path. "
                    "Returning an empty library rather than ligand entries "
                    "pointing at files that were never written.\n",
                    lib.temp_dir.c_str(), why.c_str());
                lib.temp_dir.clear();   // nothing to clean up
                lib.ligands.clear();
                lib.total = 0;
                return lib;
            }
        }

        std::ifstream in(path);
        int mol_idx = 0;
        std::string line;
        std::ostringstream mol_buf;

        while (std::getline(in, line)) {
            mol_buf << line << "\n";
            if (line.size() >= 4 && line.substr(0, 4) == "$$$$") {
                // Write this molecule
                std::string mol_path = lib.temp_dir + "/lig_" +
                                       std::to_string(mol_idx) + ".sdf";
                std::ofstream out(mol_path);
                out << mol_buf.str();
                out.close();

                // Extract molecule name from first line
                std::string mol_content = mol_buf.str();
                std::string mol_name = "lig_" + std::to_string(mol_idx);
                auto nl = mol_content.find('\n');
                if (nl != std::string::npos) {
                    std::string first_line = mol_content.substr(0, nl);
                    while (!first_line.empty() && std::isspace(static_cast<unsigned char>(first_line.back())))
                        first_line.pop_back();
                    if (!first_line.empty())
                        mol_name = first_line;
                }

                LigandEntry e;
                e.path = mol_path;
                e.name = mol_name;
                e.format = "sdf";
                e.is_temp = true;
                lib.ligands.push_back(e);

                mol_buf.str("");
                mol_buf.clear();
                mol_idx++;
            }
        }
        lib.total = static_cast<int>(lib.ligands.size());
        return lib;
    }

    // === SMILES file ===
    if (ext == ".smi" || ext == ".smiles" || ext == ".txt") {
        std::ifstream in(path);
        int idx = 0;
        std::string line;

        while (std::getline(in, line)) {
            while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
                line.pop_back();
            if (line.empty() || line[0] == '#') continue;

            // Format: SMILES [whitespace NAME]
            std::string smiles, name;
            auto tab = line.find('\t');
            if (tab == std::string::npos) tab = line.find(' ');

            if (tab != std::string::npos) {
                smiles = line.substr(0, tab);
                name = line.substr(tab + 1);
                while (!name.empty() && std::isspace(static_cast<unsigned char>(name.front())))
                    name.erase(name.begin());
            } else {
                smiles = line;
                name = "mol_" + std::to_string(idx);
            }

            LigandEntry e;
            e.path = smiles;  // path IS the SMILES string
            e.name = name;
            e.format = "smiles";
            e.is_temp = false;
            lib.ligands.push_back(e);
            idx++;
        }
        lib.total = static_cast<int>(lib.ligands.size());
        return lib;
    }

    // === Multi-model PDB (NMR ensemble / cryo-EM conformers) ===
    if (ext == ".pdb" || ext == ".ent") {
        int n_models = count_pdb_models(path);

        if (n_models <= 1) {
            LigandEntry e;
            e.path = path;
            e.name = fs::path(path).stem().string();
            e.format = "pdb";
            e.is_temp = false;
            lib.ligands.push_back(e);
            lib.total = 1;
            return lib;
        }

        // Multi-model — split at MODEL/ENDMDL boundaries
        lib.temp_dir = flexaids::temp_dir() + "/flexaid_models_" +
                       std::to_string(std::random_device{}() % 900000 + 100000);
        // Same reasoning as the SDF branch above: a multi-model receptor that
        // silently splits into zero files would dock against nothing.
        {
            std::string why;
            if (!flexaids::fs_safe::ensure_dir(lib.temp_dir, &why)) {
                std::fprintf(stderr,
                    "[SPLIT] ERROR: cannot create model split directory '%s' (%s). "
                    "Set FLEXAIDDS_TMPDIR to a writable, durable path. "
                    "Returning an empty library.\n",
                    lib.temp_dir.c_str(), why.c_str());
                lib.temp_dir.clear();
                lib.ligands.clear();
                lib.total = 0;
                return lib;
            }
        }

        std::ifstream in(path);
        std::string line;
        std::ostringstream model_buf;
        std::string header_buf;  // CRYST1, REMARK, etc. before first MODEL
        int model_idx = 0;
        bool in_model = false;
        bool seen_first_model = false;

        while (std::getline(in, line)) {
            if (line.size() >= 5 && line.substr(0, 5) == "MODEL") {
                in_model = true;
                seen_first_model = true;
                model_buf.str("");
                model_buf.clear();
                model_buf << header_buf;  // prepend shared header
                continue;
            }

            if (!seen_first_model) {
                // Accumulate header lines (CRYST1, REMARK, SCALE, etc.)
                header_buf += line + "\n";
                continue;
            }

            if (line.size() >= 6 && line.substr(0, 6) == "ENDMDL") {
                model_buf << "END\n";

                std::string model_path = lib.temp_dir + "/model_" +
                                         std::to_string(model_idx + 1) + ".pdb";
                std::ofstream out(model_path);
                out << model_buf.str();
                out.close();

                LigandEntry e;
                e.path = model_path;
                e.name = "model_" + std::to_string(model_idx + 1);
                e.format = "pdb";
                e.is_temp = true;
                lib.ligands.push_back(e);

                model_idx++;
                in_model = false;
                continue;
            }

            if (in_model) {
                model_buf << line << "\n";
            }
        }

        lib.total = static_cast<int>(lib.ligands.size());
        printf("NMR/cryo-EM ensemble: %d models extracted from %s\n",
               lib.total, path.c_str());
        return lib;
    }

    // === Single file passthrough ===
    LigandEntry e;
    e.path = path;
    e.name = fs::path(path).stem().string();
    e.format = ext.empty() ? "unknown" : ext.substr(1);
    e.is_temp = false;
    lib.ligands.push_back(e);
    lib.total = 1;
    return lib;
}

void cleanup_library(const LibraryInfo& lib) {
    if (lib.temp_dir.empty()) return;
    std::error_code ec;
    fs::remove_all(lib.temp_dir, ec);
}

} // namespace library
