// tencm.h — Torsional Elastic Network Contact Model (TENCM) for backbone flexibility
//
// Implements the torsional ENM of Delarue & Sanejouand (J. Mol. Biol. 2002) and
// Yang, Song & Cui (Biophys J. 2009):
//   – Spring network over Cα contacts within r_cutoff
//   – Torsional DOFs: one pseudo-torsion per Cα–Cα bond (i → i+1)
//   – Hessian H_kl = Σ_{contacts(i,j)} k_ij (ΔJ_ki · ΔJ_li)
//   – Normal modes via symmetric Jacobi diagonalisation
//   – Mode-weighted Boltzmann sampling for backbone perturbation during GA
//
// Metal ion support (tmcontsct):
//   – Ions (HETATM, matched by is_ion_resname) are appended as rigid pseudo-nodes
//     after all protein Cα nodes (index ≥ n_protein_ca_)
//   – Protein–ion contacts use surface-to-surface distance:
//       d_surf = r_center − r_ion − R_CA_EFF   (R_CA_EFF = 1.88 Å)
//   – Spring constant scaled by ion contact surface area:
//       k_ij = k0 * (rc/d_surf)^6 * (r_ion/R_CA_EFF)²
//   – Ion nodes have zero Jacobian (rigid; no torsional DOF)
//   – tmcontsct() exposes per-contact surface-scaled contributions for scoring
//
// Used by FlexAIDdS to generate protein backbone flexibility without rebuilding
// the full rotamer library every GA generation.
#pragma once

#include <vector>
#include <array>
#include <utility>   // std::pair, used by LigandTopology below
#include <span>
#include <concepts>
#include <cmath>
#include <memory>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <cstring>

#include "../flexaid.h"

namespace tencm {

// ─── constants ───────────────────────────────────────────────────────────────
inline constexpr float kB_kcal    = 0.001987206f; // kcal mol⁻¹ K⁻¹
inline constexpr float DEFAULT_RC = 9.0f;          // Å contact cutoff
inline constexpr float DEFAULT_K0 = 1.0f;          // spring constant (kcal mol⁻¹ Å⁻²)
inline constexpr int   N_MODES    = 20;            // low-frequency modes kept

// ─── data structures ─────────────────────────────────────────────────────────

// Elastic contact between Cα i and Cα j (i < j)
struct Contact {
    int   i, j;
    float k;     // spring constant k_ij = k0 * (rc/r0)^6
    float r0;    // equilibrium Cα–Cα distance (or surface-to-surface for ions)
};

// Per-contact surface-area-weighted spring score (tmcontsct).
// For protein–protein contacts: k_scaled == k (area_scale = 1).
// For protein–ion contacts: k_scaled = k0*(rc/d_surf)^6 * (r_ion/R_CA_EFF)²
//   where d_surf = r_center − r_ion − R_CA_EFF (surface-to-surface distance).
struct TmContSct {
    int   i, j;        // node indices (same as Contact; j >= n_protein_ca_ for ions)
    float k_scaled;    // area-scaled spring constant (kcal mol⁻¹ Å⁻²)
    float d_surf;      // equilibrium distance used (Å)
    bool  is_ion;      // true when one node is a metal ion pseudo-node
};

// Pseudo-torsion DOF: rotation about bond between Cα_k and Cα_{k+1}
struct PseudoBond {
    int   k;          // 0-based index: bond connects residue k to k+1
    float axis[3];    // unit vector (r_{k+1} - r_k) / |...|
    float pivot[3];   // midpoint (r_k + r_{k+1}) / 2
};

// One normal mode: eigenvalue (stiffness) + eigenvector over torsion DOFs
struct NormalMode {
    double eigenvalue;                // kcal mol⁻¹ rad⁻²
    std::vector<double> eigenvector; // length = n_bonds
};

// A perturbed backbone conformation
struct Conformer {
    std::vector<float> delta_theta;           // torsion perturbations (rad)
    std::vector<std::array<float,3>> ca;      // perturbed Cα positions
    float strain_energy;                      // ½ δθᵀ H δθ (kcal/mol)
};

// ─── concepts ────────────────────────────────────────────────────────────────
template<typename T>
concept FloatLike = std::floating_point<T>;

// ─── main engine ─────────────────────────────────────────────────────────────
class TorsionalENM {
public:
    // Build network from parsed atom/residue arrays.
    // Selects Cα atoms from protein residues.
    void build(const atom*  atoms,
               const resid* residue,
               int          res_cnt,
               float        cutoff = DEFAULT_RC,
               float        k0     = DEFAULT_K0);

    // Sample one perturbed backbone conformation.
    // temperature: Kelvin; rng: seeded generator passed in for reproducibility.
    Conformer sample(float temperature, std::mt19937& rng) const;

    // Apply a Conformer back to the FA atoms array (in-place coordinate update).
    void apply(const Conformer& conf,
               atom*            atoms,
               const resid*     residue) const;

    // Predicted Cα B-factors at given T (Å²)
    std::vector<float> bfactors(float temperature) const;

    // Build directly from Cα coordinates (no FA atom/resid dependency).
    // ca_coords: sequential Cα positions (x,y,z) along the chain.
    void build_from_ca(const std::vector<std::array<float,3>>& ca_coords,
                       float cutoff = DEFAULT_RC,
                       float k0     = DEFAULT_K0);

    // Build a Cartesian 3N×3N Anisotropic Network Model (ANM) Hessian over the
    // ligand HEAVY atoms in the half-open index range [lig_start, lig_end) of
    // atoms[].  Unlike build()/build_from_ca() (which assemble a *torsional*
    // pseudo-bond Hessian over Cα backbone DOFs), this assembles the classic
    // Cartesian ANM super-element Hessian and diagonalises it directly into
    // modes_.  Eigenvalues (model-scale stiffness, λ ≥ 0; ~6 rigid-body modes
    // near zero) are exposed via .modes(); eigenvectors are left empty since the
    // Level-3 H(ω) diagnostic consumes eigenvalues only.
    //
    // Contact potential matches build_from_ca: step-function within `cutoff`
    // with spring k_ij = k0 * (cutoff / r0)^6.  Hydrogens (element "H") are
    // excluded; if fewer than 3 heavy atoms are found the model is not built
    // (is_built() == false) and modes() is empty.
    void build_from_ligand(const atom* atoms,
                           int   lig_start,
                           int   lig_end,
                           float cutoff = 7.0f,
                           float k0     = DEFAULT_K0);

    // Build the ligand ENM Hessian in INTERNAL (dihedral) coordinates -- the same
    // basis build()/build_from_ca() use for the protein.
    //
    // WHY THIS EXISTS. The thermodynamic cycle
    //     dS_vib(bind) = S_vib(complex) - S_vib(apo receptor) - S_vib(free ligand)
    // is only a differential if all three terms live in ONE basis. Until this
    // function, the two protein terms were Torsional (dihedral DOFs, no rigid-body
    // null space) while the ligand term came from build_from_ligand and was
    // Cartesian 3N (six rigid-body modes stripped by a cutoff). Subtracting across
    // bases is not a validity-preserving operation: the bases differ not only in
    // dimension but in WHICH MODES EXIST. The same incoherence reached the
    // objective, where ic2cf.cpp adds tencom_weight * cf.h_rep with h_rep taken
    // from the Cartesian path.
    //
    // DOF DEFINITION, matching the GA's own. FlexAID already represents the ligand
    // in internal coordinates: buildcc() reconstructs Cartesians from each atom's
    // rec[0..2] frame plus dis/ang/dih, and a ligand dihedral gene sets
    // atoms[FA->map_par[i].atm].dih directly (FOPTICS.cpp:563), with dependent
    // atoms following via .shift (:571). So one torsional DOF is taken per DISTINCT
    // rotatable bond rec[1]--rec[0] carried by a reconstruction-flagged
    // (recs == 'm') ligand atom. That is deliberately the GA's own DOF set rather
    // than an independent perception, so the entropy is computed over the same
    // coordinates the search moves in. Callers should assert
    // n_torsion_dofs() == FA->nflexbonds; a mismatch means the two perceptions have
    // diverged and the number is not comparable to a search result.
    //
    // Hessian: H_kl = sum_contacts k_ij (u_ij . dJ_k)(u_ij . dJ_l), identical in
    // form to the protein assembly, with the Jacobian J_k[a] = u_k x (r_a - r_pivot)
    // for atoms downstream of DOF k's bond and 0 otherwise.
    //
    // DEGENERATE CASES ARE REAL PHYSICS, NOT FAILURES, and are reported rather than
    // silently zeroed. Measured on Astex-84: torsional DOF count has median 4 and
    // range 0-11, against 3N-6 median 63 for the Cartesian path. Four ligands
    // (1GPK 1Q41 1U4D 1W1P) have ZERO rotatable bonds and five (1HNN 1P2Y 1XOZ
    // 1YV3 2GBP) have exactly one. A rigid ligand genuinely has no internal
    // vibrational entropy to lose, so the thermodynamic S_vib is 0 (an empty sum)
    // and well-defined -- but a 32-bin SHANNON entropy of the log-frequency
    // spectrum is UNDEFINED over an empty spectrum and identically 0 over a single
    // mode. Those are different functionals and only the first is valid at low DOF
    // count. is_built() is false when fewer than 2 DOFs exist, so a caller cannot
    // mistake a degenerate spectrum for a computed one.
    void build_from_ligand_torsional(const atom* atoms,
                                     int   lig_start,
                                     int   lig_end,
                                     float cutoff = 7.0f,
                                     float k0     = DEFAULT_K0);

    /// Ligand topology in NODE-INDEX space (0..n-1 over the coordinate array
    /// passed alongside it), so a caller that has no FlexAID atom[] can still
    /// build the torsional model.
    ///
    /// WHY THIS FORM EXISTS. Every post-hoc consumer of the ligand ENM reads a
    /// finished pose file and has coordinates and elements only -- recs, rec[0..2]
    /// and bond[] do not survive into a PDB. Handing the atom[] overload such an
    /// array yields ZERO DOFs and a refusal for every ligand, indistinguishable
    /// from "all ligands are rigid". The engine therefore emits its own DOF list
    /// as a sidecar at dock time (top.cpp, <name>_ligtopo.json) and this overload
    /// consumes it, so the entropy is computed over the SEARCH's coordinates
    /// rather than a re-perception of the bond graph.
    struct LigandTopology {
        /// Rotatable bonds as (pivot_node, distal_node). Must come from the
        /// engine's flexbond range, NOT from map_par typ==2: typ==2 also covers
        /// the ligand's rigid-body placement pseudo-dihedrals, which carry
        /// bnd = -1 and no second reference atom. Measured on 1JD0: typ==2 gives
        /// 3 where nflexbonds is 1.
        std::vector<std::pair<int,int>> rot_bonds;
        /// Undirected adjacency over node indices, used to find the atoms
        /// downstream of each rotated bond.
        std::vector<std::vector<int>>   adjacency;
    };

    /// Same Hessian, same spring, same basis as the atom[] overload -- only the
    /// source of the topology differs. Refuses below 2 DOFs for the same reason.
    void build_from_ligand_torsional(const std::array<float,3>* xyz,
                                     int   n_nodes,
                                     const LigandTopology& topo,
                                     float cutoff = 7.0f,
                                     float k0     = DEFAULT_K0);

private:
    /// Shared core: assembles and diagonalises the torsional Hessian from a node
    /// coordinate array plus explicit topology. BOTH public overloads delegate
    /// here, so the two entry points cannot drift into different physics -- the
    /// same consolidation this file already applied to the rigid-mode cutoff
    /// after the CLI and the objective disagreed by 56% on one ligand.
    void assemble_torsional_(const std::array<float,3>* xyz,
                             int n_nodes,
                             const LigandTopology& topo,
                             float cutoff,
                             float k0);
public:

    /// Number of internal-coordinate DOFs the torsional ligand path found.
    /// Zero on every other build path. Compare against FA->nflexbonds.
    int n_torsion_dofs() const noexcept { return n_torsion_dofs_; }

    // Getters
    int n_residues()    const noexcept { return static_cast<int>(ca_.size()); }
    int n_protein_ca()  const noexcept { return n_protein_ca_; }
    int n_bonds()       const noexcept { return static_cast<int>(bonds_.size()); }
    const std::vector<NormalMode>&  modes()      const noexcept { return modes_; }
    bool is_built()                              const noexcept { return built_; }
    const std::vector<std::array<float,3>>& ca_positions() const noexcept { return ca_; }

    /// Which coordinate system the assembled Hessian lives in.
    ///
    /// Torsional (build/build_from_ca): dihedral DOFs. Global translation and
    ///   rotation are not representable, so there is NO rigid-body null space.
    /// Cartesian (build_from_ligand): 3N DOFs, so the spectrum carries SIX
    ///   rigid-body modes whose eigenvalues are numerically zero.
    enum class Basis { Torsional, Cartesian };
    Basis basis() const noexcept { return basis_; }

    /// Eigenvalues of the VIBRATIONAL subspace, rigid-body modes removed.
    ///
    /// WHY THIS EXISTS. Filtering with `eigenvalue > 0.0` does NOT remove the
    /// rigid-body subspace on the Cartesian path: those eigenvalues land at
    /// ~1e-12 with round-off SIGNS, so roughly half survive as positive and
    /// enter the spectrum. MEASURED on this ligand path (BU72, 32 heavy atoms,
    /// 96 modes): H_pop = 1.896 with them against 2.950 without -- a 56%
    /// change -- because log(1e-13) = -30 becomes the lower bin edge and
    /// stretches the log-frequency span from 9.65 to 39.0, compressing every
    /// real mode into a quarter of the 32 bins. The contamination is in the
    /// BINNING GEOMETRY, not the histogram mass. Worse, the survivor count is
    /// unstable: a 1e-6 A coordinate jitter flips it between 2 and 4 of 6, so
    /// the old value carried a component set by floating-point sign noise and
    /// was not reproducible across BLAS builds or thread counts.
    ///
    /// The cutoff is RELATIVE to the stiffest mode, max(1e-10, lam_max*1e-8),
    /// matching the two call sites that already got this right
    /// (ligand_tencom_pose.cpp:85, DatasetRunner.cpp:974) so every consumer
    /// agrees on one rule. The margin is enormous: measured separation between
    /// the largest rigid eigenvalue (7e-12) and the softest real mode (3.2) is
    /// ~6.6e11x, so the threshold sits in eleven decades of empty space.
    ///
    /// Negative eigenvalues fail the comparison rather than reaching log(),
    /// so no NaN can propagate into H(omega).
    ///
    /// `n_dropped` receives how many modes were removed; `n_expected` receives
    /// 6 on the Cartesian path and 0 on the torsional one. A MISMATCH means the
    /// contact graph fragmented (each disconnected component contributes its
    /// own rigid modes) and is reported to stderr rather than silently absorbed
    /// into the entropy.
    std::vector<double> vibrational_eigenvalues(int* n_dropped  = nullptr,
                                                int* n_expected = nullptr) const;

    /// Per-contact surface-area-weighted spring scores.
    /// Ion contacts (is_ion==true) have k_scaled = k0*(rc/d_surf)^6 * (r_ion/R_CA_EFF)²,
    /// enabling callers to weight vibrational entropy by ion contact surface area.
    const std::vector<TmContSct>& tmcontsct() const noexcept { return tmcontsct_; }

private:
    // Coordinate system of the assembled Hessian. Defaults to Torsional and is
    // set explicitly by EVERY builder, so a reused object cannot carry a stale
    // basis from a previous build of the other kind.
    Basis basis_ = Basis::Torsional;

    // Number of internal-coordinate DOFs found by build_from_ligand_torsional.
    // Left at 0 by every other builder so a caller reading it after the wrong
    // build path gets 0 rather than a stale count from a previous object use.
    int n_torsion_dofs_ = 0;

    // Internal node coordinate store:
    //   indices [0, n_protein_ca_)  → protein Cα atoms
    //   indices [n_protein_ca_, N)  → metal ion pseudo-nodes (rigid)
    std::vector<std::array<float,3>> ca_;
    // Map: sequential node index → atom index in FA atoms[]
    std::vector<int> ca_atom_idx_;
    // Map: sequential node index → residue index (1-based FA convention)
    std::vector<int> res_idx_;

    // Number of protein Cα nodes; ion nodes are appended after this boundary.
    int n_protein_ca_ = 0;
    // VdW radius for each ion node (indexed by node_index − n_protein_ca_).
    std::vector<float> ion_radii_;

    std::vector<Contact>    contacts_;
    std::vector<TmContSct>  tmcontsct_;  // surface-area-weighted per-contact scores
    std::vector<PseudoBond> bonds_;
    std::vector<NormalMode> modes_;

    // Hessian stored as dense symmetric matrix (n_bonds × n_bonds)
    std::vector<double> H_;

    // Cached Jacobian matrix: J_cached_[k*N + i] = {Jx, Jy, Jz}
    // Pre-computed once during assemble_hessian(), reused by bfactors()/sample().
    std::vector<std::array<float,3>> J_cached_;
    bool jac_cached_ = false;

    bool  built_  = false;
    float cutoff_ = DEFAULT_RC;
    float k0_     = DEFAULT_K0;

    // Build steps
    void extract_ca(const atom* atoms, const resid* residue, int res_cnt);
    void build_contacts();
    void build_bonds();
    void assemble_hessian();
    void diagonalize();     // Jacobi iteration on H_

    // Jacobian: ∂r_{atom_i}/∂θ_{bond_k}.
    // Returns {0,0,0} when atom_i is upstream of bond_k.
    std::array<float,3> jac(int bond_k, int atom_i) const noexcept;

    // Jacobi sweep helper
    static void jacobi_rotate(std::vector<double>& A,
                               std::vector<double>& V,
                               int n, int p, int q) noexcept;
};

// ─── free function: quick fluctuation amplitude at given T ───────────────────
// Returns the rms displacement (Å) of residue i predicted by the model.
float residue_rms_fluctuation(const TorsionalENM& tencm,
                              int residue_idx,
                              float temperature);

}  // namespace tencm
