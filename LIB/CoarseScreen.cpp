// CoarseScreen.cpp — NRGRank coarse-grained screening for FlexAIDdS
//
// Full C++20 translation of:
//   - process_target.py: build_index_cubes(), get_cf_list(), load_ligand_test_dots(),
//     clean_binding_site_grid(), get_clash_per_dot()
//   - rank_molecules.py: center_coords(), rotate_ligand(), get_cf(), get_cf_main(),
//     get_cf_with_clash(), get_cf_main_clash()
//
// Reference:
//   DesCôteaux T, Mailhot O, Najmanovich RJ. "NRGRank: Coarse-grained
//   structurally-informed ultra-massive virtual screening."
//   bioRxiv 2025.02.17.638675.
//
// SPDX-License-Identifier: Apache-2.0

#include "CoarseScreen.h"
#include "nrgrank_matrix.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// ─── SIMD includes (optional AVX2) ────────────────────────────────────
#if defined(FLEXAIDS_USE_AVX2) || defined(FLEXAIDS_USE_AVX512)
#include <immintrin.h>
#endif

namespace nrgrank {

// ═══════════════════════════════════════════════════════════════════════
//  IndexCubeGrid — spatial hash mapping atoms to cubes
// ═══════════════════════════════════════════════════════════════════════

void IndexCubeGrid::build(const std::vector<TargetAtom>& atoms, float cell_width) {
    if (atoms.empty()) return;

    cell_width_ = cell_width;

    // If cell_width not specified, compute from max radius + water
    if (cell_width_ <= 0.f) {
        float max_rad = 0.f;
        for (auto& a : atoms)
            max_rad = std::max(max_rad, a.radius);
        cell_width_ = 2.f * (max_rad + kWaterRadius);
    }

    // Compute bounding box with padding
    Vec3 lo{atoms[0].pos}, hi{atoms[0].pos};
    for (auto& a : atoms) {
        lo.x = std::min(lo.x, a.pos.x);
        lo.y = std::min(lo.y, a.pos.y);
        lo.z = std::min(lo.z, a.pos.z);
        hi.x = std::max(hi.x, a.pos.x);
        hi.y = std::max(hi.y, a.pos.y);
        hi.z = std::max(hi.z, a.pos.z);
    }

    // Pad by cell_width (cw_factor=1 in Python)
    min_xyz_ = {lo.x - cell_width_, lo.y - cell_width_, lo.z - cell_width_};
    max_xyz_ = {hi.x + cell_width_, hi.y + cell_width_, hi.z + cell_width_};

    // Grid dimensions: ((max - min) / cell_width).astype(int32) + 1
    nx_ = static_cast<int>((max_xyz_.x - min_xyz_.x) / cell_width_) + 1;
    ny_ = static_cast<int>((max_xyz_.y - min_xyz_.y) / cell_width_) + 1;
    nz_ = static_cast<int>((max_xyz_.z - min_xyz_.z) / cell_width_) + 1;

    // First pass: count atoms per cell to find max_per_cell
    const int n_cells = nx_ * ny_ * nz_;
    std::vector<int> counts(n_cells, 0);

    for (size_t i = 0; i < atoms.size(); ++i) {
        int ix = static_cast<int>((atoms[i].pos.x - min_xyz_.x) / cell_width_);
        int iy = static_cast<int>((atoms[i].pos.y - min_xyz_.y) / cell_width_);
        int iz = static_cast<int>((atoms[i].pos.z - min_xyz_.z) / cell_width_);
        ix = std::clamp(ix, 0, nx_-1);
        iy = std::clamp(iy, 0, ny_-1);
        iz = std::clamp(iz, 0, nz_-1);
        counts[(ix * ny_ + iy) * nz_ + iz]++;
    }

    max_per_cell_ = *std::max_element(counts.begin(), counts.end());
    if (max_per_cell_ == 0) max_per_cell_ = 1;

    // Allocate grid filled with placeholder
    const size_t total = static_cast<size_t>(n_cells) * max_per_cell_;
    grid_.assign(total, kPlaceholder);

    // Second pass: fill grid
    std::fill(counts.begin(), counts.end(), 0);
    for (size_t i = 0; i < atoms.size(); ++i) {
        int ix = static_cast<int>((atoms[i].pos.x - min_xyz_.x) / cell_width_);
        int iy = static_cast<int>((atoms[i].pos.y - min_xyz_.y) / cell_width_);
        int iz = static_cast<int>((atoms[i].pos.z - min_xyz_.z) / cell_width_);
        ix = std::clamp(ix, 0, nx_-1);
        iy = std::clamp(iy, 0, ny_-1);
        iz = std::clamp(iz, 0, nz_-1);
        int flat = (ix * ny_ + iy) * nz_ + iz;
        int slot = counts[flat]++;
        grid_[static_cast<size_t>(flat) * max_per_cell_ + slot] = static_cast<int>(i);
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  File I/O
// ═══════════════════════════════════════════════════════════════════════

std::vector<TargetAtom> parse_target_mol2(const std::string& path) {
    std::vector<TargetAtom> atoms;
    std::ifstream fin(path);
    if (!fin.is_open()) return atoms;

    bool in_atoms = false;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        // Trim leading whitespace
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;

        if (line.find("@<TRIPOS>ATOM") != std::string::npos) {
            in_atoms = true;
            continue;
        }
        if (in_atoms && line[first] == '@') {
            break; // end of atom block
        }

        if (!in_atoms) continue;

        // Parse: atom_id name x y z type [res_id res_name charge]
        std::istringstream iss(line);
        int atom_id;
        std::string name;
        float x, y, z;
        std::string sybyl_type;

        if (!(iss >> atom_id >> name >> x >> y >> z >> sybyl_type)) continue;

        // Skip hydrogens (first element before '.' is 'H')
        std::string elem = sybyl_type.substr(0, sybyl_type.find('.'));
        if (elem == "H" || elem == "h") continue;

        TargetAtom ta;
        ta.pos = {x, y, z};
        ta.type = sybyl_type_lookup(sybyl_type, ta.radius);
        atoms.push_back(ta);
    }
    return atoms;
}

std::vector<BindingSiteSphere> parse_binding_site_pdb(const std::string& path) {
    std::vector<BindingSiteSphere> spheres;
    std::ifstream fin(path);
    if (!fin.is_open()) return spheres;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.size() < 66) continue;
        if (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM") continue;

        BindingSiteSphere s;
        try {
            s.center.x = std::stof(line.substr(30, 8));
            s.center.y = std::stof(line.substr(38, 8));
            s.center.z = std::stof(line.substr(46, 8));
            s.radius   = std::stof(line.substr(60, 6));
        } catch (...) {
            continue;
        }
        spheres.push_back(s);
    }
    return spheres;
}

// ═══════════════════════════════════════════════════════════════════════
//  CoarseScreener
// ═══════════════════════════════════════════════════════════════════════

CoarseScreener::CoarseScreener()  = default;
CoarseScreener::~CoarseScreener() = default;

bool CoarseScreener::load_target_mol2(const std::string& path) {
    target_atoms_ = parse_target_mol2(path);
    return !target_atoms_.empty();
}

bool CoarseScreener::load_binding_site_pdb(const std::string& path) {
    binding_spheres_ = parse_binding_site_pdb(path);
    return !binding_spheres_.empty();
}

void CoarseScreener::prepare_target(const std::vector<TargetAtom>& atoms,
                                     const std::vector<BindingSiteSphere>& spheres) {
    target_atoms_    = atoms;
    binding_spheres_ = spheres;
    prepared_        = false;

    if (target_atoms_.empty() || binding_spheres_.empty()) return;

    build_grid();
    precompute_cf();
    generate_anchors();
    if (config_.use_clash) {
        build_clash_grid();
    }
    clean_anchors();

    prepared_ = true;
}

// ───── Grid construction ───────────────────────────────────────────────

void CoarseScreener::build_grid() {
    grid_.build(target_atoms_, config_.cell_width);
}

// ───── CF precomputation (get_cf_list kernel) ──────────────────────────
// This is the core hot loop from process_target.py, translated faithfully.
// 6 nested loops: x,y,z grid × atom_types × 3³ neighbors × atoms_in_cell

void CoarseScreener::precompute_cf() {
    const int gx = grid_.nx();
    const int gy = grid_.ny();
    const int gz = grid_.nz();
    cf_nx_ = gx;
    cf_ny_ = gy;
    cf_nz_ = gz;

    const int n_types = NUM_ATOM_TYPES; // 39
    const size_t total = static_cast<size_t>(gx) * gy * gz * n_types;
    cf_grid_.assign(total, 0.f);

    const auto& ematf = EnergyMatrixF::instance();
    const int max_per_cell = grid_.max_atoms_per_cell();

    // OpenMP parallelise over x-slabs
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) collapse(2)
    #endif
    for (int x = 0; x < gx; ++x) {
        for (int y = 0; y < gy; ++y) {
            for (int z = 0; z < gz; ++z) {
                // For each atom type (1..39)
                for (int at_idx = 0; at_idx < n_types; ++at_idx) {
                    const int atom_type = at_idx + 1; // 1-based
                    float cf = 0.f;

                    // 27-neighbor scan
                    for (int di = -1; di <= 1; ++di) {
                        const int ni = x + di;
                        for (int dj = -1; dj <= 1; ++dj) {
                            const int nj = y + dj;
                            for (int dk = -1; dk <= 1; ++dk) {
                                const int nk = z + dk;

                                // Python uses strict inequality: 0 < i < len
                                // This means boundary cubes (index 0 and max-1)
                                // are excluded from neighbor lookup
                                if (ni <= 0 || ni >= gx ||
                                    nj <= 0 || nj >= gy ||
                                    nk <= 0 || nk >= gz)
                                    continue;

                                const int* cell_data = grid_.cell(ni, nj, nk);
                                if (cell_data[0] == IndexCubeGrid::kPlaceholder)
                                    continue;

                                // Sum energy for all atoms in this neighbor cell
                                for (int s = 0; s < max_per_cell; ++s) {
                                    const int neighbour = cell_data[s];
                                    if (neighbour == IndexCubeGrid::kPlaceholder)
                                        break;
                                    const int type_2 = target_atoms_[neighbour].type;
                                    cf += ematf.data[atom_type][type_2];
                                }
                            }
                        }
                    }

                    // Store: cf_grid_[(((x*ny)+y)*nz+z)*39 + at_idx]
                    const size_t idx = static_cast<size_t>(((x * cf_ny_) + y) * cf_nz_ + z) * n_types + at_idx;
                    cf_grid_[idx] = cf;
                }
            }
        }
    }
}

// ───── Anchor point generation (load_ligand_test_dots) ─────────────────

void CoarseScreener::generate_anchors() {
    anchor_points_.clear();
    if (binding_spheres_.empty()) return;

    const float sep = config_.test_dot_separation;

    // Bounding box around binding site spheres (accounting for sphere radii)
    float xlo = binding_spheres_[0].center.x - binding_spheres_[0].radius;
    float xhi = binding_spheres_[0].center.x + binding_spheres_[0].radius;
    float ylo = binding_spheres_[0].center.y - binding_spheres_[0].radius;
    float yhi = binding_spheres_[0].center.y + binding_spheres_[0].radius;
    float zlo = binding_spheres_[0].center.z - binding_spheres_[0].radius;
    float zhi = binding_spheres_[0].center.z + binding_spheres_[0].radius;

    for (size_t i = 1; i < binding_spheres_.size(); ++i) {
        const auto& s = binding_spheres_[i];
        xlo = std::min(xlo, s.center.x - s.radius);
        xhi = std::max(xhi, s.center.x + s.radius);
        ylo = std::min(ylo, s.center.y - s.radius);
        yhi = std::max(yhi, s.center.y + s.radius);
        zlo = std::min(zlo, s.center.z - s.radius);
        zhi = std::max(zhi, s.center.z + s.radius);
    }

    // Generate grid of test dots, keeping only those inside a binding sphere
    for (float dx = xlo; dx < xhi; dx += sep) {
        for (float dy = ylo; dy < yhi; dy += sep) {
            for (float dz = zlo; dz < zhi; dz += sep) {
                // Round to 3 decimals like Python
                float px = std::round(dx * 1000.f) / 1000.f;
                float py = std::round(dy * 1000.f) / 1000.f;
                float pz = std::round(dz * 1000.f) / 1000.f;

                Vec3 point{px, py, pz};

                // Check if inside any binding sphere
                for (const auto& s : binding_spheres_) {
                    Vec3 diff = point - s.center;
                    float dist = diff.length();
                    if (dist < s.radius) {
                        anchor_points_.push_back(point);
                        break;
                    }
                }
            }
        }
    }
}

// ───── Clean anchors (clean_binding_site_grid) ─────────────────────────
// Remove anchor points that are within 2.0 Å of any target atom

void CoarseScreener::clean_anchors() {
    if (anchor_points_.empty()) return;

    const float clash_dist = config_.clash_distance; // 2.0 Å
    const float clash_dist_sq = clash_dist * clash_dist;
    const int max_per_cell = grid_.max_atoms_per_cell();

    std::vector<bool> remove(anchor_points_.size(), false);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64)
    #endif
    for (size_t a = 0; a < anchor_points_.size(); ++a) {
        const Vec3& point = anchor_points_[a];
        int gx, gy, gz;
        grid_.world_to_grid(point, gx, gy, gz);

        bool found_clash = false;
        for (int di = -1; di <= 1 && !found_clash; ++di) {
            for (int dj = -1; dj <= 1 && !found_clash; ++dj) {
                for (int dk = -1; dk <= 1 && !found_clash; ++dk) {
                    int ni = gx + di, nj = gy + dj, nk = gz + dk;
                    if (!grid_.in_bounds(ni, nj, nk)) continue;

                    const int* cell_data = grid_.cell(ni, nj, nk);
                    if (cell_data[0] == IndexCubeGrid::kPlaceholder) continue;

                    for (int s = 0; s < max_per_cell; ++s) {
                        int idx = cell_data[s];
                        if (idx == IndexCubeGrid::kPlaceholder) break;

                        Vec3 diff = target_atoms_[idx].pos - point;
                        float dist_sq = diff.length_sq();
                        if (dist_sq <= clash_dist_sq) {
                            found_clash = true;
                            break;
                        }
                    }
                }
            }
        }
        remove[a] = found_clash;
    }

    // Compact
    std::vector<Vec3> cleaned;
    cleaned.reserve(anchor_points_.size());
    for (size_t i = 0; i < anchor_points_.size(); ++i) {
        if (!remove[i]) cleaned.push_back(anchor_points_[i]);
    }
    anchor_points_ = std::move(cleaned);
}

// ───── Clash grid (get_clash_per_dot) ──────────────────────────────────

void CoarseScreener::build_clash_grid() {
    if (binding_spheres_.empty()) return;

    const float spacing = config_.clash_dot_distance;
    clash_spacing_ = spacing;
    const float padding = config_.bd_site_padding;
    const float clash_dist = config_.clash_distance;

    // Bounding box with padding (same as make_binding_site_cuboid)
    float xlo = binding_spheres_[0].center.x - binding_spheres_[0].radius;
    float xhi = binding_spheres_[0].center.x + binding_spheres_[0].radius;
    float ylo = binding_spheres_[0].center.y - binding_spheres_[0].radius;
    float yhi = binding_spheres_[0].center.y + binding_spheres_[0].radius;
    float zlo = binding_spheres_[0].center.z - binding_spheres_[0].radius;
    float zhi = binding_spheres_[0].center.z + binding_spheres_[0].radius;

    for (size_t i = 1; i < binding_spheres_.size(); ++i) {
        const auto& s = binding_spheres_[i];
        xlo = std::min(xlo, s.center.x - s.radius);
        xhi = std::max(xhi, s.center.x + s.radius);
        ylo = std::min(ylo, s.center.y - s.radius);
        yhi = std::max(yhi, s.center.y + s.radius);
        zlo = std::min(zlo, s.center.z - s.radius);
        zhi = std::max(zhi, s.center.z + s.radius);
    }

    xlo -= padding; xhi += padding;
    ylo -= padding; yhi += padding;
    zlo -= padding; zhi += padding;

    clash_origin_ = {xlo, ylo, zlo};

    // Number of cells along each axis
    clash_nx_ = static_cast<int>((xhi - xlo) / spacing) + 1;
    clash_ny_ = static_cast<int>((yhi - ylo) / spacing) + 1;
    clash_nz_ = static_cast<int>((zhi - zlo) / spacing) + 1;

    const size_t total = static_cast<size_t>(clash_nx_) * clash_ny_ * clash_nz_;
    clash_grid_.assign(total, 0);

    const int max_per_cell = grid_.max_atoms_per_cell();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) collapse(2)
    #endif
    for (int a = 0; a < clash_nx_; ++a) {
        for (int b = 0; b < clash_ny_; ++b) {
            for (int c = 0; c < clash_nz_; ++c) {
                float x_val = xlo + a * spacing;
                float y_val = ylo + b * spacing;
                float z_val = zlo + c * spacing;
                Vec3 point{x_val, y_val, z_val};

                int gx, gy, gz;
                grid_.world_to_grid(point, gx, gy, gz);

                bool clash = false;
                for (int di = -1; di <= 1 && !clash; ++di) {
                    for (int dj = -1; dj <= 1 && !clash; ++dj) {
                        for (int dk = -1; dk <= 1 && !clash; ++dk) {
                            int ni = gx + di, nj = gy + dj, nk = gz + dk;
                            if (ni < 0 || ni >= grid_.nx() ||
                                nj < 0 || nj >= grid_.ny() ||
                                nk < 0 || nk >= grid_.nz())
                                continue;

                            const int* cell_data = grid_.cell(ni, nj, nk);
                            if (cell_data[0] == IndexCubeGrid::kPlaceholder) continue;

                            for (int s = 0; s < max_per_cell; ++s) {
                                int idx = cell_data[s];
                                if (idx == IndexCubeGrid::kPlaceholder) break;
                                Vec3 diff = target_atoms_[idx].pos - point;
                                if (diff.length() <= clash_dist) {
                                    clash = true;
                                    break;
                                }
                            }
                        }
                    }
                }

                clash_grid_[static_cast<size_t>((a * clash_ny_) + b) * clash_nz_ + c] =
                    clash ? 1 : 0;
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  Rotation generation
// ═══════════════════════════════════════════════════════════════════════

std::vector<std::array<float,9>> CoarseScreener::generate_rotations(int per_axis) {
    const float step = 360.f / per_axis;
    const float deg2rad = static_cast<float>(M_PI) / 180.f;

    // Pre-generate all rotation matrices, then deduplicate
    std::vector<std::array<float,9>> matrices;
    matrices.reserve(static_cast<size_t>(per_axis) * per_axis * per_axis);

    for (int rx = 0; rx < per_axis; ++rx) {
        float ax = rx * step * deg2rad;
        float cx = std::cos(ax), sx = std::sin(ax);
        for (int ry = 0; ry < per_axis; ++ry) {
            float ay = ry * step * deg2rad;
            float cy = std::cos(ay), sy = std::sin(ay);
            for (int rz = 0; rz < per_axis; ++rz) {
                float az = rz * step * deg2rad;
                float cz = std::cos(az), sz = std::sin(az);

                // R = Rx * Ry * Rz (same order as Python)
                // Rx = [[1,0,0],[0,cx,-sx],[0,sx,cx]]
                // Ry = [[cy,0,sy],[0,1,0],[-sy,0,cy]]
                // Rz = [[cz,-sz,0],[sz,cz,0],[0,0,1]]
                std::array<float,9> m{};
                m[0] = cy*cz;
                m[1] = -cy*sz;
                m[2] = sy;
                m[3] = sx*sy*cz + cx*sz;
                m[4] = -sx*sy*sz + cx*cz;
                m[5] = -sx*cy;
                m[6] = -cx*sy*cz + sx*sz;
                m[7] = cx*sy*sz + sx*cz;
                m[8] = cx*cy;

                matrices.push_back(m);
            }
        }
    }

    // Deduplicate (matching Python's np.unique on axis=0)
    // Sort and remove near-duplicates (within 1e-4 tolerance)
    auto approx_less = [](const std::array<float,9>& a, const std::array<float,9>& b) {
        for (int i = 0; i < 9; ++i) {
            if (a[i] < b[i] - 1e-4f) return true;
            if (a[i] > b[i] + 1e-4f) return false;
        }
        return false;
    };
    std::sort(matrices.begin(), matrices.end(), approx_less);

    auto approx_eq = [](const std::array<float,9>& a, const std::array<float,9>& b) {
        for (int i = 0; i < 9; ++i)
            if (std::fabs(a[i] - b[i]) > 1e-4f) return false;
        return true;
    };
    matrices.erase(std::unique(matrices.begin(), matrices.end(), approx_eq),
                   matrices.end());

    return matrices;
}

// ═══════════════════════════════════════════════════════════════════════
//  Pose scoring — get_cf() from rank_molecules.py
// ═══════════════════════════════════════════════════════════════════════

float CoarseScreener::score_pose(const float* coords, const int* atom_types,
                                  int num_atoms, const Vec3& anchor) const {
    const float cw    = grid_.cell_width();
    const Vec3  mxyz  = grid_.min_xyz();
    const float inv_cw = 1.f / cw;

    // Translate ligand to anchor and compute grid indices for each atom
    // Check bounds first
    for (int i = 0; i < num_atoms; ++i) {
        float px = coords[i*3+0] + anchor.x;
        float py = coords[i*3+1] + anchor.y;
        float pz = coords[i*3+2] + anchor.z;

        int ix = static_cast<int>((px - mxyz.x) * inv_cw);
        int iy = static_cast<int>((py - mxyz.y) * inv_cw);
        int iz = static_cast<int>((pz - mxyz.z) * inv_cw);

        if (ix < 0 || iy < 0 || iz < 0 ||
            ix >= cf_nx_ || iy >= cf_ny_ || iz >= cf_nz_) {
            return config_.default_cf;
        }
    }

    // Sum CF
    float cf = 0.f;
    for (int i = 0; i < num_atoms; ++i) {
        float px = coords[i*3+0] + anchor.x;
        float py = coords[i*3+1] + anchor.y;
        float pz = coords[i*3+2] + anchor.z;

        int ix = static_cast<int>((px - mxyz.x) * inv_cw);
        int iy = static_cast<int>((py - mxyz.y) * inv_cw);
        int iz = static_cast<int>((pz - mxyz.z) * inv_cw);

        int at_idx = atom_types[i] - 1; // 0-based for cf_grid_
        if (at_idx < 0 || at_idx >= NUM_ATOM_TYPES) {
            return config_.default_cf;
        }

        float temp_cf = cf_at(ix, iy, iz, at_idx);
        if (temp_cf == config_.default_cf) {
            return config_.default_cf;
        }
        cf += temp_cf;
    }
    return cf;
}

float CoarseScreener::score_pose_with_clash(const float* coords, const int* atom_types,
                                             int num_atoms, const Vec3& anchor) const {
    // First check clash grid
    const float inv_spacing = 1.f / clash_spacing_;
    const Vec3& corigin = clash_origin_;

    for (int i = 0; i < num_atoms; ++i) {
        float px = coords[i*3+0] + anchor.x;
        float py = coords[i*3+1] + anchor.y;
        float pz = coords[i*3+2] + anchor.z;

        // Round to nearest clash grid cell (matching Python's np.round)
        int cx = static_cast<int>(std::round((px - corigin.x) * inv_spacing));
        int cy = static_cast<int>(std::round((py - corigin.y) * inv_spacing));
        int cz = static_cast<int>(std::round((pz - corigin.z) * inv_spacing));

        // Clamp to boundary (Python: arr[arr == size] -= 1)
        if (cx == clash_nx_) cx--;
        if (cy == clash_ny_) cy--;
        if (cz == clash_nz_) cz--;

        if (cx < 0 || cy < 0 || cz < 0 ||
            cx >= clash_nx_ || cy >= clash_ny_ || cz >= clash_nz_) {
            return config_.default_cf;
        }

        if (clash_grid_[static_cast<size_t>((cx * clash_ny_) + cy) * clash_nz_ + cz]) {
            return config_.default_cf;
        }
    }

    // No clash — compute CF score
    return score_pose(coords, atom_types, num_atoms, anchor);
}

// ═══════════════════════════════════════════════════════════════════════
//  Screening — main ligand ranking loop
// ═══════════════════════════════════════════════════════════════════════

ScreenResult CoarseScreener::screen_one(const ScreenLigand& ligand) const {
    return screen_one_impl(ligand, /*use_omp=*/true);
}

ScreenResult CoarseScreener::screen_one_impl(const ScreenLigand& ligand,
                                              bool use_omp) const {
    ScreenResult result;
    result.name  = ligand.name;
    result.score = config_.default_cf;

    if (!prepared_ || ligand.atoms.empty() || anchor_points_.empty())
        return result;

    const int n_atoms = static_cast<int>(ligand.atoms.size());

    // Center coordinates (subtract centroid) — matches center_coords()
    float cx = 0.f, cy = 0.f, cz = 0.f;
    for (const auto& a : ligand.atoms) {
        cx += a.pos.x; cy += a.pos.y; cz += a.pos.z;
    }
    cx /= n_atoms; cy /= n_atoms; cz /= n_atoms;

    std::vector<float> centered(n_atoms * 3);
    std::vector<int>   types(n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        centered[i*3+0] = ligand.atoms[i].pos.x - cx;
        centered[i*3+1] = ligand.atoms[i].pos.y - cy;
        centered[i*3+2] = ligand.atoms[i].pos.z - cz;
        types[i]         = ligand.atoms[i].type;
    }

    // Generate rotation matrices — 9³=729, deduplicated
    auto rotations = generate_rotations(config_.rotations_per_axis);
    const int n_rot = static_cast<int>(rotations.size());

    // Pre-rotate all orientations: [n_rot][n_atoms*3]
    std::vector<float> all_rotated(static_cast<size_t>(n_rot) * n_atoms * 3);
    for (int r = 0; r < n_rot; ++r) {
        const auto& m = rotations[r];
        float* out = all_rotated.data() + static_cast<size_t>(r) * n_atoms * 3;
        for (int a = 0; a < n_atoms; ++a) {
            float ix = centered[a*3+0];
            float iy = centered[a*3+1];
            float iz = centered[a*3+2];
            out[a*3+0] = m[0]*ix + m[1]*iy + m[2]*iz;
            out[a*3+1] = m[3]*ix + m[4]*iy + m[5]*iz;
            out[a*3+2] = m[6]*ix + m[7]*iy + m[8]*iz;
        }
    }

    float best_cf = config_.default_cf;
    int best_rot = -1, best_anc = -1;

    const bool use_clash = config_.use_clash && !clash_grid_.empty();
    const int n_anchors = static_cast<int>(anchor_points_.size());

    // Inner scoring loop — optionally parallelised with OpenMP
    auto score_fn = [&](const float* pose, int n, const Vec3& anc) -> float {
        return use_clash
            ? score_pose_with_clash(pose, types.data(), n, anc)
            : score_pose(pose, types.data(), n, anc);
    };

    #ifdef _OPENMP
    if (use_omp) {
        #pragma omp parallel
        {
            float local_best = config_.default_cf;
            int local_best_rot = -1, local_best_anc = -1;

            #pragma omp for schedule(dynamic, 4) nowait
            for (int ai = 0; ai < n_anchors; ++ai) {
                const Vec3& anchor = anchor_points_[ai];
                for (int ri = 0; ri < n_rot; ++ri) {
                    const float* pose = all_rotated.data() +
                        static_cast<size_t>(ri) * n_atoms * 3;
                    float cf = score_fn(pose, n_atoms, anchor);
                    if (cf < local_best ||
                        (cf == local_best && ai < local_best_anc) ||
                        (cf == local_best && ai == local_best_anc && ri < local_best_rot)) {
                        local_best     = cf;
                        local_best_rot = ri;
                        local_best_anc = ai;
                    }
                }
            }

            #pragma omp critical
            {
                if (local_best < best_cf ||
                    (local_best == best_cf && local_best_anc < best_anc) ||
                    (local_best == best_cf && local_best_anc == best_anc && local_best_rot < best_rot)) {
                    best_cf  = local_best;
                    best_rot = local_best_rot;
                    best_anc = local_best_anc;
                }
            }
        }
    } else
    #else
    (void)use_omp;
    #endif
    {
        // Serial path
        for (int ai = 0; ai < n_anchors; ++ai) {
            const Vec3& anchor = anchor_points_[ai];
            for (int ri = 0; ri < n_rot; ++ri) {
                const float* pose = all_rotated.data() +
                    static_cast<size_t>(ri) * n_atoms * 3;
                float cf = score_fn(pose, n_atoms, anchor);
                if (cf < best_cf ||
                    (cf == best_cf && ai < best_anc) ||
                    (cf == best_cf && ai == best_anc && ri < best_rot)) {
                    best_cf  = cf;
                    best_rot = ri;
                    best_anc = ai;
                }
            }
        }
    }

    result.score         = best_cf;
    result.best_rotation = best_rot;
    result.best_anchor   = best_anc;
    if (best_anc >= 0)
        result.best_position = anchor_points_[best_anc];

    return result;
}

std::vector<ScreenResult> CoarseScreener::screen(
        const std::vector<ScreenLigand>& ligands) const {

    std::vector<ScreenResult> results(ligands.size());

    // For large libraries, parallelize over ligands (each ligand scored sequentially).
    // For small libraries, the inner OpenMP in screen_one() handles it.
    const bool parallelize_outer = ligands.size() > 4;

    if (parallelize_outer) {
        // When parallelising over ligands, each ligand is scored serially
        // to avoid nested OMP overhead. We inline the scoring loop.
        #ifdef _OPENMP
        const int n = static_cast<int>(ligands.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < n; ++i) {
            results[i] = screen_one_impl(ligands[i], /*use_omp=*/false);
        }
        #else
        for (size_t i = 0; i < ligands.size(); ++i)
            results[i] = screen_one_impl(ligands[i], /*use_omp=*/false);
        #endif
    } else {
        for (size_t i = 0; i < ligands.size(); ++i)
            results[i] = screen_one_impl(ligands[i], /*use_omp=*/true);
    }

    // Sort by score (best = most negative first)
    std::sort(results.begin(), results.end(),
              [](const ScreenResult& a, const ScreenResult& b) {
                  return a.score < b.score;
              });

    return results;
}

// ═══════════════════════════════════════════════════════════════════════
//  Ligand file loading
// ═══════════════════════════════════════════════════════════════════════

std::vector<ScreenLigand> CoarseScreener::load_ligands_mol2(const std::string& path) {
    std::vector<ScreenLigand> ligands;
    std::ifstream fin(path);
    if (!fin.is_open()) return ligands;

    std::string line;
    ScreenLigand current;
    bool in_molecule = false;
    bool in_atoms    = false;
    bool got_name    = false;

    while (std::getline(fin, line)) {
        if (line.find("@<TRIPOS>MOLECULE") != std::string::npos) {
            // Save previous molecule if any
            if (!current.atoms.empty()) {
                ligands.push_back(std::move(current));
                current = ScreenLigand{};
            }
            in_molecule = true;
            in_atoms    = false;
            got_name    = false;
            continue;
        }

        if (in_molecule && !got_name) {
            // Next non-empty line is the molecule name
            size_t first = line.find_first_not_of(" \t\r\n");
            if (first != std::string::npos) {
                current.name = line.substr(first);
                // Trim trailing whitespace
                size_t last = current.name.find_last_not_of(" \t\r\n");
                if (last != std::string::npos)
                    current.name = current.name.substr(0, last + 1);
                got_name = true;
            }
            continue;
        }

        if (line.find("@<TRIPOS>ATOM") != std::string::npos) {
            in_atoms = true;
            in_molecule = false;
            continue;
        }

        if (in_atoms && !line.empty() && line[0] == '@') {
            in_atoms = false;
            continue;
        }

        if (in_atoms) {
            std::istringstream iss(line);
            int atom_id;
            std::string name;
            float x, y, z;
            std::string sybyl_type;

            if (!(iss >> atom_id >> name >> x >> y >> z >> sybyl_type)) continue;

            // Skip hydrogens
            std::string elem = sybyl_type.substr(0, sybyl_type.find('.'));
            if (elem == "H" || elem == "h") continue;

            LigandAtom la;
            la.pos = {x, y, z};
            float dummy_rad;
            la.type = sybyl_type_lookup(sybyl_type, dummy_rad);
            current.atoms.push_back(la);
        }
    }

    // Don't forget last molecule
    if (!current.atoms.empty())
        ligands.push_back(std::move(current));

    return ligands;
}

std::vector<ScreenLigand> CoarseScreener::load_ligands_sdf(const std::string& path) {
    std::vector<ScreenLigand> ligands;
    std::ifstream fin(path);
    if (!fin.is_open()) return ligands;

    std::string line;
    while (std::getline(fin, line)) {
        ScreenLigand lig;

        // Line 1: molecule name
        lig.name = line;
        size_t last = lig.name.find_last_not_of(" \t\r\n");
        if (last != std::string::npos)
            lig.name = lig.name.substr(0, last + 1);

        // Line 2: program/timestamp
        if (!std::getline(fin, line)) break;
        // Line 3: comment
        if (!std::getline(fin, line)) break;
        // Line 4: counts
        if (!std::getline(fin, line)) break;

        int n_atoms = 0, n_bonds = 0;
        std::istringstream iss(line);
        iss >> n_atoms >> n_bonds;

        // Atom block
        for (int i = 0; i < n_atoms; ++i) {
            if (!std::getline(fin, line)) break;
            std::istringstream aiss(line);
            float x, y, z;
            std::string symbol;
            if (!(aiss >> x >> y >> z >> symbol)) continue;

            if (symbol == "H" || symbol == "h") continue;

            LigandAtom la;
            la.pos = {x, y, z};

            // SDF has element symbols, not SYBYL types.
            // Use basic element→SYBYL mapping
            float dummy_rad;
            std::string sybyl;
            if      (symbol == "C")  sybyl = "C.3";
            else if (symbol == "N")  sybyl = "N.3";
            else if (symbol == "O")  sybyl = "O.3";
            else if (symbol == "S")  sybyl = "S.3";
            else if (symbol == "P")  sybyl = "P.3";
            else if (symbol == "F")  sybyl = "F";
            else if (symbol == "Cl") sybyl = "CL";
            else if (symbol == "Br") sybyl = "BR";
            else if (symbol == "I")  sybyl = "I";
            // ── Transition-metal / metalloid cofactor atoms ──────────────────
            // Metal-containing ligands (organometallics, metallodrugs, Se in
            // selenoproteins/selenocysteine analogs) previously fell through to
            // "DUMMY" and were typed as inert spheres, zeroing their pairwise
            // interaction-matrix terms.  Map each to its dedicated SYBYL token so
            // sybyl_type_lookup() returns the correct radius and emat row.
            // Token spellings must match kSybylTypes[] in nrgrank_matrix.h
            // (matched case-insensitively, spaces stripped).  Co only exists
            // there as the "CO.OH" cobalt entry; "NI"/"FE"/"ZN"/... map 1:1.
            else if (symbol == "Fe") sybyl = "FE";
            else if (symbol == "Zn") sybyl = "ZN";
            else if (symbol == "Mg") sybyl = "MG";
            else if (symbol == "Ca") sybyl = "CA";
            else if (symbol == "Mn") sybyl = "MN";
            else if (symbol == "Cu") sybyl = "CU";
            else if (symbol == "Co") sybyl = "CO.OH";
            else if (symbol == "Ni") sybyl = "NI";
            else if (symbol == "Se") sybyl = "SE";
            else                     sybyl = "DUMMY";

            la.type = sybyl_type_lookup(sybyl, dummy_rad);
            lig.atoms.push_back(la);
        }

        // Skip rest until "$$$$"
        while (std::getline(fin, line)) {
            if (line.find("$$$$") != std::string::npos) break;
        }

        if (!lig.atoms.empty())
            ligands.push_back(std::move(lig));
    }

    return ligands;
}

} // namespace nrgrank
