// RefLigSeed.h — Reference ligand seeding for GA initialization
// Reads a reference ligand PDB/MOL2 to bias initial population toward
// a known binding mode. Header-only, Apache-2.0 © 2026 Le Bonhomme Pharma

#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <numeric>

#include "flexaid.h"

namespace reflig {

// ── Lightweight coordinate extraction ───────────────────────────────────────
// Parses atom coordinates from PDB or MOL2 without needing FA_Global.

struct RefLigAtom {
    float x, y, z;
};

struct RefLigData {
    std::vector<RefLigAtom> atoms;
    float centroid[3];
    std::vector<int> nearest_grid;  // K nearest grid point indices (1-based)
};

// ── Parse PDB atom coordinates ──────────────────────────────────────────────
inline std::vector<RefLigAtom> parse_pdb_coords(const std::string& path) {
    std::vector<RefLigAtom> atoms;
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return atoms;

    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0) {
            if (strlen(line) >= 54) {
                RefLigAtom a;
                a.x = static_cast<float>(atof(std::string(line + 30, 8).c_str()));
                a.y = static_cast<float>(atof(std::string(line + 38, 8).c_str()));
                a.z = static_cast<float>(atof(std::string(line + 46, 8).c_str()));
                atoms.push_back(a);
            }
        }
    }
    fclose(f);
    return atoms;
}

// ── Parse MOL2 atom coordinates ─────────────────────────────────────────────
inline std::vector<RefLigAtom> parse_mol2_coords(const std::string& path) {
    std::vector<RefLigAtom> atoms;
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return atoms;

    char line[256];
    bool in_atom_block = false;

    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "@<TRIPOS>ATOM", 13) == 0) {
            in_atom_block = true;
            continue;
        }
        if (line[0] == '@') {
            in_atom_block = false;
            continue;
        }
        if (in_atom_block) {
            RefLigAtom a;
            // MOL2 ATOM: id name x y z type [res_id res_name charge]
            int id;
            char name[32];
            if (sscanf(line, "%d %31s %f %f %f", &id, name, &a.x, &a.y, &a.z) >= 5) {
                atoms.push_back(a);
            }
        }
    }
    fclose(f);
    return atoms;
}

// ── Parse SDF/MOL atom coordinates ─────────────────────────────────────────
inline std::vector<RefLigAtom> parse_sdf_coords(const std::string& path) {
    std::vector<RefLigAtom> atoms;
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return atoms;

    char line[256];
    int line_no = 0;
    int atom_count = -1;
    while (fgets(line, sizeof(line), f)) {
        ++line_no;
        if (line_no == 4) {
            if (std::strlen(line) >= 3) {
                atom_count = std::atoi(std::string(line, 3).c_str());
            }
            if (atom_count <= 0) break;
            for (int i = 0; i < atom_count && fgets(line, sizeof(line), f); ++i) {
                if (std::strlen(line) < 34) continue;
                char elem_buf[4] = {0, 0, 0, 0};
                std::strncpy(elem_buf, line + 31, 3);
                std::string elem(elem_buf);
                elem.erase(std::remove_if(elem.begin(), elem.end(), ::isspace), elem.end());
                if (elem == "H" || elem == "D" || elem == "Du") continue;
                RefLigAtom a;
                a.x = static_cast<float>(std::atof(std::string(line, 0, 10).c_str()));
                a.y = static_cast<float>(std::atof(std::string(line, 10, 10).c_str()));
                a.z = static_cast<float>(std::atof(std::string(line, 20, 10).c_str()));
                atoms.push_back(a);
            }
            break;
        }
    }
    fclose(f);
    return atoms;
}

// ── Auto-detect file type and parse ─────────────────────────────────────────
inline std::vector<RefLigAtom> parse_reflig(const std::string& path) {
    // Check extension
    std::string lower_path = path;
    std::transform(lower_path.begin(), lower_path.end(), lower_path.begin(), ::tolower);

    if (lower_path.size() >= 4 && lower_path.substr(lower_path.size() - 4) == ".pdb") {
        return parse_pdb_coords(path);
    }
    if (lower_path.size() >= 5 && lower_path.substr(lower_path.size() - 5) == ".mol2") {
        return parse_mol2_coords(path);
    }
    if ((lower_path.size() >= 4 && lower_path.substr(lower_path.size() - 4) == ".sdf") ||
        (lower_path.size() >= 4 && lower_path.substr(lower_path.size() - 4) == ".mol")) {
        return parse_sdf_coords(path);
    }
    // Try PDB first, then MOL2, then SDF/MOL.
    auto atoms = parse_pdb_coords(path);
    if (atoms.empty()) atoms = parse_mol2_coords(path);
    if (atoms.empty()) atoms = parse_sdf_coords(path);
    return atoms;
}

// ── Compute centroid ────────────────────────────────────────────────────────
inline void compute_centroid(const std::vector<RefLigAtom>& atoms, float out[3]) {
    out[0] = out[1] = out[2] = 0.0f;
    if (atoms.empty()) return;
    for (const auto& a : atoms) {
        out[0] += a.x;
        out[1] += a.y;
        out[2] += a.z;
    }
    float n = static_cast<float>(atoms.size());
    out[0] /= n;
    out[1] /= n;
    out[2] /= n;
}

// ── Find K nearest grid points to centroid ──────────────────────────────────
inline std::vector<int> find_nearest_grid_points(
    const float centroid[3],
    const gridpoint* cleftgrid, int num_grd,
    int k)
{
    if (num_grd <= 1) return {};

    // Compute distances from centroid to each grid point (skip index 0)
    struct IndexDist {
        int index;
        float dist_sq;
    };
    std::vector<IndexDist> distances;
    distances.reserve(static_cast<std::size_t>(num_grd - 1));

    for (int i = 1; i < num_grd; ++i) {
        float dx = cleftgrid[i].coor[0] - centroid[0];
        float dy = cleftgrid[i].coor[1] - centroid[1];
        float dz = cleftgrid[i].coor[2] - centroid[2];
        distances.push_back({i, dx*dx + dy*dy + dz*dz});
    }

    // Partial sort to find K smallest
    int actual_k = std::min(k, static_cast<int>(distances.size()));
    std::partial_sort(distances.begin(),
                      distances.begin() + actual_k,
                      distances.end(),
                      [](const IndexDist& a, const IndexDist& b) {
                          return a.dist_sq < b.dist_sq;
                      });

    std::vector<int> result;
    result.reserve(static_cast<std::size_t>(actual_k));
    for (int i = 0; i < actual_k; ++i)
        result.push_back(distances[static_cast<std::size_t>(i)].index);
    return result;
}

// ── Prepare reference ligand seeding data ───────────────────────────────────
inline RefLigData prepare_reflig_seed(
    const std::string& reflig_file,
    const gridpoint* cleftgrid, int num_grd,
    int k_nearest = 10)
{
    RefLigData data;
    data.atoms = parse_reflig(reflig_file);
    if (data.atoms.empty()) {
        data.centroid[0] = data.centroid[1] = data.centroid[2] = 0.0f;
        return data;
    }

    compute_centroid(data.atoms, data.centroid);
    data.nearest_grid = find_nearest_grid_points(
        data.centroid, cleftgrid, num_grd, k_nearest);

    return data;
}

} // namespace reflig
