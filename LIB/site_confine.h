// site_confine.h — Cognate-site grid confinement decision
//
// Direct-mode / --redock SITE-CONFINE (top.cpp) trims the SURFNET grid to
// points near the cognate ligand centroid. A previous MIN_SITE_GRID=500
// floor silently fell back to the full protein grid when the pocket was
// sparse (1STP biotin: keep=188 → 20 Å poses). Sparse-but-nonempty pockets
// must stay confined.
//
// Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
#pragma once

namespace flexaids {

/// Rebuild the cleft grid from `keep_count` cognate-site points?
/// `total_grid_pts` is the searchable count (typically FA->num_grd - 1,
/// excluding the reflig reference slot). `hard_floor` is a crash guard
/// for degenerate 1–few-point grids; MIN_SITE_GRID is a warning-only
/// density target and must not trigger full-grid fallback.
inline bool site_confine_should_rebuild(int keep_count, int total_grid_pts,
                                        int hard_floor = 25) {
    return keep_count >= hard_floor && keep_count < total_grid_pts;
}

} // namespace flexaids
