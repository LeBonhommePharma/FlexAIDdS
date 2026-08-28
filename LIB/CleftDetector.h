#ifndef CLEFT_DETECTOR_H
#define CLEFT_DETECTOR_H

/*
 * CleftDetector — automatic binding-site detection for FlexAID
 *
 * Implements a SURFNET-style gap-sphere algorithm (the same geometric
 * principle used by GetCleft):
 *   1. For every pair of protein surface atoms within a distance cutoff,
 *      place a probe sphere midway between them. Ligand / HETATM residues
 *      (residue.type == 1) are excluded so a loaded cognate ligand cannot
 *      occupy the pocket during SURFNET (required for --redock).
 *   2. Shrink the probe until no other protein atom penetrates it
 *      (or discard if radius falls below a minimum).
 *   3. Cluster surviving spheres by spatial proximity (single-linkage)
 *      and keep EVERY cluster meeting min_cluster_size (all genuine pockets,
 *      not just the largest — the cognate site is often a smaller cavity).
 *      Downstream site-confinement (top.cpp) trims the grid to the cognate
 *      pocket, so handing it all pockets guarantees the right one is present.
 *   4. Return a linked list of sphere_struct* ready for generate_grid().
 *
 * The implementation is header + .cpp so it can be compiled as part of
 * the FlexAID executable without any new dependencies beyond what
 * CMakeLists.txt already pulls in (Eigen3 optional, OpenMP optional).
 */

#include "flexaid.h"
#include <vector>
#include <string>

struct CleftDetectorParams {
    float max_pair_dist;     // max distance between atom pair for probe placement (A)
    float probe_radius_max;  // initial probe sphere radius (A)
    float probe_radius_min;  // minimum acceptable probe radius (A)
    float probe_shrink_step; // radius decrement per iteration (A)
    float cluster_cutoff;    // single-linkage clustering distance (A)
    int   min_cluster_size;  // discard clusters smaller than this
    // Layer 2: keep at most this many ligandable clusters (volume×enclosure).
    // Default 5 — larger than multi-cleft fan-out; 0 = keep all min_cluster_size.
    int   top_k_clefts;
    // Optional spatial pre-filter: when oracle_radius > 0, SURFNET only processes
    // atoms within this radius of oracle_center.  Eliminates O(N^3) blowup on
    // multimeric receptors (e.g. 1OF6 octamer, 20826 atoms): a 15 A filter reduces
    // the working set to ~200-400 atoms, yielding a ~10,000x speedup.
    float oracle_center[3];  // centroid of oracle binding site (A)
    float oracle_radius;     // pre-filter radius (A); 0.0 = disabled
};

// Default parameters matching typical GetCleft behaviour
inline CleftDetectorParams default_cleft_params() {
    CleftDetectorParams p;
    p.max_pair_dist    = 12.0f;
    p.probe_radius_max =  5.0f;
    p.probe_radius_min =  1.5f;
    p.probe_shrink_step=  0.1f;
    p.cluster_cutoff   =  4.0f;
    p.min_cluster_size =  10;
    p.top_k_clefts     =  5;   // ligandable top-K (ensemble layer 2)
    p.oracle_center[0] = 0.0f; p.oracle_center[1] = 0.0f; p.oracle_center[2] = 0.0f;
    p.oracle_radius    = 0.0f;  // disabled by default
    return p;
}

/*  detect_cleft
 *
 *  atoms    – protein atom array (already read by read_pdb)
 *  residue  – residue array
 *  atm_cnt  – total number of atoms
 *  res_cnt  – total number of residues
 *  params   – tuning knobs (use default_cleft_params() for sane defaults)
 *
 *  Returns a linked list of sphere_struct* identical to what
 *  read_spheres() produces, so it plugs straight into generate_grid().
 *  Caller owns the memory (free with free_sphere_list).
 */
sphere* detect_cleft(const atom* atoms, const resid* residue,
                     int atm_cnt, int res_cnt,
                     const CleftDetectorParams& params = default_cleft_params());

/*  write_cleft_spheres
 *
 *  Writes detected spheres to a PDB-format sphere file
 *  (same format read_spheres expects), useful for caching / inspection.
 */
void write_cleft_spheres(const sphere* spheres, const char* filename);

/*  free_sphere_list – frees the linked list returned by detect_cleft */
void free_sphere_list(sphere* head);

// ─── Task 8: Cleft Annotation & Flexible Residue Selection (PREPROCESSING) ──
// This module turns cleft geometry + optional active-site info into
// human-readable annotation and a recommended set of flexible residues.
// It is strictly preprocessing — it must never change scoring behaviour.
//
// External annotations (Pfam, UniProt, etc.) must come from user-supplied
// files or a separate preprocessing script. No web lookups are allowed here.

struct CleftAnnotation {
    int         cleft_id = -1;
    std::string class_label;                    // "orthosteric", "allosteric", "unknown"
    double      confidence = 0.0;               // 0.0 – 1.0
    double      volume_A3 = 0.0;
    double      distance_to_active_site_A = -1.0;

    std::vector<int>         nearby_active_site_residues;
    std::vector<int>         recommended_flexible_residues;
    std::vector<std::string> evidence;          // free-text reasons
};

// Flexible residue selector rules (Task 8):
// - Include residues within 'distance_shell_A' of any cleft sphere or ligand atom.
// - Optionally include 'active_site_residues' if provided.
// - Never include Gly/Ala unless a backbone flexibility module is active (we ignore for now).
// - Respect 'user_fixed_residues' (these are excluded).
// - Force-include anything in 'user_forced_flexible' (unless invalid).
// - Deduplicate and sort deterministically (chain then residue number).
std::vector<int> select_flexible_residues(
    const atom* atoms,
    const resid* residue,
    int atm_cnt,
    int res_cnt,
    const std::vector<int>& cleft_sphere_residues,   // residues near detected cleft
    double distance_shell_A,
    const std::vector<int>& active_site_residues = {},
    const std::vector<int>& user_fixed_residues = {},
    const std::vector<int>& user_forced_flexible = {}
);

#endif // CLEFT_DETECTOR_H
