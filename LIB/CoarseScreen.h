// CoarseScreen.h — coarse-grained cube screening for FlexAIDdS
//
// Reimplemented from the published method (DesCôteaux, Mailhot, Najmanovich,
// bioRxiv 2025.02.17.638675v2, doi 10.1101/2025.02.17.638675;
// supplementary data Zenodo 10.5281/zenodo.16861024), clean-room per
// docs/licensing/clean-room-policy.md. Apache-2.0. Not a translation of
// NRGlab/NRGRank source (that tree is GPL-3.0 and is not a dependency).
//
// Pipeline:
//   1. Build index-cube grid (cell_width=6.56 Å) from target atoms
//   2. Precompute CF energy at every cube × atom-type (27-neighbor kernel)
//   3. Generate anchor points (1.5 Å spacing), remove clashing ones
//   4. For each ligand: generate 729 rotations × anchor translations,
//      look up precomputed CF, return best score
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace nrgrank {

// ───── Data structures ─────────────────────────────────────────────────

/// 3D coordinate
struct Vec3 {
    float x = 0.f, y = 0.f, z = 0.f;

    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(float s)       const { return {x*s, y*s, z*s}; }
    float dot(const Vec3& o)      const { return x*o.x + y*o.y + z*o.z; }
    float length_sq()             const { return x*x + y*y + z*z; }
    float length()                const { return std::sqrt(length_sq()); }
};

/// Target atom (from MOL2)
struct TargetAtom {
    Vec3  pos;
    int   type   = 0;   // 1–39 SYBYL type
    float radius = 0.f;
};

/// Binding-site sphere (from PDB cleft definition)
struct BindingSiteSphere {
    Vec3  center;
    float radius = 0.f;
};

/// Ligand atom for coarse screening
struct LigandAtom {
    Vec3  pos;
    int   type = 0;  // 1–39 SYBYL type
};

/// A single ligand to screen
struct ScreenLigand {
    std::string name;
    std::vector<LigandAtom> atoms;
};

/// Result for one screened ligand
struct ScreenResult {
    std::string name;
    float       score         = 0.f;
    int         best_rotation = -1;
    int         best_anchor   = -1;
    /// Centroid position of best pose (anchor + ligand centroid offset)
    Vec3        best_position{};
};

// ───── Index-cube 3D spatial grid ──────────────────────────────────────

/// Spatial grid that maps atoms into cubes for O(1) neighbor lookup.
/// Mirrors NRGRank's build_index_cubes().
class IndexCubeGrid {
public:
    // Default parameters matching NRGRank
    static constexpr float kWaterRadius   = 1.4f;
    static constexpr float kCellWidth     = 6.56f;
    static constexpr int   kPlaceholder   = -1;

    /// Build from target atoms. If cell_width <= 0, computes automatically.
    void build(const std::vector<TargetAtom>& atoms, float cell_width = kCellWidth);

    /// Grid dimensions
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    int max_atoms_per_cell() const { return max_per_cell_; }

    Vec3 min_xyz() const { return min_xyz_; }
    Vec3 max_xyz() const { return max_xyz_; }
    float cell_width() const { return cell_width_; }

    /// Get atom indices in cell (ix, iy, iz). Returns pointer to array of
    /// max_per_cell_ ints, terminated by kPlaceholder.
    const int* cell(int ix, int iy, int iz) const {
        return grid_.data() + static_cast<size_t>(((ix * ny_) + iy) * nz_ + iz) * max_per_cell_;
    }

    /// Convert world position to grid indices (truncating, like Python astype(int32))
    void world_to_grid(const Vec3& pos, int& ix, int& iy, int& iz) const {
        ix = static_cast<int>((pos.x - min_xyz_.x) / cell_width_);
        iy = static_cast<int>((pos.y - min_xyz_.y) / cell_width_);
        iz = static_cast<int>((pos.z - min_xyz_.z) / cell_width_);
    }

    /// Check if grid index is in bounds
    bool in_bounds(int ix, int iy, int iz) const {
        return ix >= 0 && ix < nx_ && iy >= 0 && iy < ny_ && iz >= 0 && iz < nz_;
    }

private:
    std::vector<int> grid_;   // flat [nx][ny][nz][max_per_cell]
    Vec3  min_xyz_{}, max_xyz_{};
    float cell_width_ = kCellWidth;
    int   nx_ = 0, ny_ = 0, nz_ = 0;
    int   max_per_cell_ = 0;
};

// ───── Configuration ───────────────────────────────────────────────────

struct CoarseScreenConfig {
    // Target preprocessing
    float cell_width          = 6.56f;
    float water_radius        = 1.4f;
    float test_dot_separation = 1.5f;
    float clash_distance      = 2.0f;
    float bd_site_padding     = 2.0f;
    float clash_dot_distance  = 0.25f;

    // Ligand screening
    int   rotations_per_axis  = 9;   // 9³ = 729 orientations
    bool  use_clash           = true;
    float default_cf          = 1e8f;

    // Output
    int   top_n               = 100;  // return top-N for full docking
};

// ───── Main screener class ─────────────────────────────────────────────

/// Full NRGRank coarse-grained screener.
///
/// Usage:
///   CoarseScreener cs;
///   cs.set_config(cfg);
///   cs.prepare_target(target_atoms, binding_site_spheres);
///   auto results = cs.screen(ligands);
///
class CoarseScreener {
public:
    CoarseScreener();
    ~CoarseScreener();

    /// Set screening parameters
    void set_config(const CoarseScreenConfig& cfg) { config_ = cfg; }
    const CoarseScreenConfig& config() const { return config_; }

    // ── Target preparation (paper §2 cube/anchor construction) ──

    /// Load target from MOL2 file
    bool load_target_mol2(const std::string& path);

    /// Load binding site from PDB cleft file
    bool load_binding_site_pdb(const std::string& path);

    /// Prepare target from pre-loaded data
    void prepare_target(const std::vector<TargetAtom>& atoms,
                        const std::vector<BindingSiteSphere>& spheres);

    /// Check if target is prepared and ready for screening
    bool is_prepared() const { return prepared_; }

    /// Number of anchor points (after clash removal)
    size_t num_anchors() const { return anchor_points_.size(); }

    /// CF grid dimensions
    int cf_nx() const { return grid_.nx(); }
    int cf_ny() const { return grid_.ny(); }
    int cf_nz() const { return grid_.nz(); }

    // ── Ligand screening (rank_molecules.py equivalent) ──

    /// Screen a batch of ligands. Returns results sorted by score (best first).
    /// Thread-safe: uses OpenMP internally, but do not call concurrently.
    std::vector<ScreenResult> screen(const std::vector<ScreenLigand>& ligands) const;

    /// Screen a single ligand (convenience wrapper).
    ScreenResult screen_one(const ScreenLigand& ligand) const;

    // ── MOL2/SDF file loading utilities ──

    /// Load ligands from multi-mol2 file
    static std::vector<ScreenLigand> load_ligands_mol2(const std::string& path);

    /// Load ligands from SDF file
    static std::vector<ScreenLigand> load_ligands_sdf(const std::string& path);

private:
    CoarseScreenConfig config_;

    // Target data
    std::vector<TargetAtom>          target_atoms_;
    std::vector<BindingSiteSphere>   binding_spheres_;
    IndexCubeGrid                    grid_;
    bool                             prepared_ = false;

    // Precomputed CF grid: [nx][ny][nz][NUM_ATOM_TYPES]
    // Flat storage, indexed as cf_grid_[(((x*ny)+y)*nz+z)*39 + (type-1)]
    std::vector<float> cf_grid_;
    int cf_nx_ = 0, cf_ny_ = 0, cf_nz_ = 0;

    // Anchor points (cleaned binding site dots)
    std::vector<Vec3> anchor_points_;

    // Clash grid: [cx][cy][cz] — true if clashing
    std::vector<uint8_t> clash_grid_;
    int clash_nx_ = 0, clash_ny_ = 0, clash_nz_ = 0;
    Vec3 clash_origin_{};
    float clash_spacing_ = 0.25f;

    // ── Internal methods ──

    /// Build index cube grid
    void build_grid();

    /// Precompute CF energies (get_cf_list kernel)
    void precompute_cf();

    /// Generate anchor points from binding site
    void generate_anchors();

    /// Remove clashing anchor points (clean_binding_site_grid)
    void clean_anchors();

    /// Build clash detection grid (get_clash_per_dot)
    void build_clash_grid();

public:
    /// Generate rotation matrices (729 total for 9 per axis)
    static std::vector<std::array<float,9>> generate_rotations(int per_axis);

private:

    /// Score a single ligand orientation at a single anchor
    float score_pose(const float* centered_coords, const int* atom_types,
                     int num_atoms, const Vec3& anchor) const;

    /// Score with clash detection
    float score_pose_with_clash(const float* centered_coords, const int* atom_types,
                                int num_atoms, const Vec3& anchor) const;

    /// Access precomputed CF value
    float cf_at(int ix, int iy, int iz, int atom_type_minus1) const {
        return cf_grid_[static_cast<size_t>(((ix * cf_ny_) + iy) * cf_nz_ + iz) * 39 + atom_type_minus1];
    }

    /// Internal implementation with optional OMP parallelism
    ScreenResult screen_one_impl(const ScreenLigand& ligand, bool use_omp) const;
};

// ───── File I/O helpers ────────────────────────────────────────────────

/// Parse target atoms from MOL2 file (skipping hydrogens)
std::vector<TargetAtom> parse_target_mol2(const std::string& path);

/// Parse binding site spheres from PDB cleft file
std::vector<BindingSiteSphere> parse_binding_site_pdb(const std::string& path);

} // namespace nrgrank
