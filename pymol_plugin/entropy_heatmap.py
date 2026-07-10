"""Entropy heatmap visualization for FlexAID∆S binding modes.

Computes spatial entropy density from pose ensembles and renders
it as a pseudoatom-based heatmap in PyMOL (Phase 3, deliverable 3.1).

Usage:
    PyMOL> flexaids_load_results /path/to/output
    PyMOL> flexaids_entropy_heatmap 1
    PyMOL> flexaids_entropy_heatmap 1, renderer=cgo
"""

from __future__ import annotations

import math
from typing import List, Optional, Tuple

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

try:
    from pymol import cmd
except ImportError as exc:
    raise ImportError("PyMOL not available") from exc

# flexaidds is only required for type/runtime mode objects; pure helpers
# (_read_pose_coords, entropy color, grid math) work without it.
try:
    from flexaidds import BindingModeResult  # noqa: F401
except ImportError:
    BindingModeResult = object  # type: ignore[misc, assignment]


def _read_pose_coords(
    pdb_path: str,
    ligand_only: bool = True,
) -> List[Tuple[float, float, float]]:
    """Extract heavy-atom coordinates from a PDB file.

    Args:
        pdb_path: Path to PDB file.
        ligand_only: If True, only read HETATM records (ligand atoms).
            This dramatically reduces coordinate count and improves
            heatmap performance.  Set False to include protein atoms.

    Returns a list of (x, y, z) tuples for non-hydrogen atoms.
    """
    coords = []
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if ligand_only:
                    if not line.startswith("HETATM"):
                        continue
                else:
                    if not (line.startswith("ATOM") or line.startswith("HETATM")):
                        continue
                element = line[76:78].strip() if len(line) >= 78 else ""
                atom_name = line[12:16].strip()
                if element == "H" or (not element and atom_name.startswith("H")):
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
    except (OSError, ValueError):
        pass
    # Fall back to all atoms if no HETATM found
    if not coords and ligand_only:
        return _read_pose_coords(pdb_path, ligand_only=False)
    return coords


def _compute_grid_bounds(
    all_coords: List[Tuple[float, float, float]],
    padding: float = 3.0,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Compute bounding box around all coordinates with padding."""
    xs = [c[0] for c in all_coords]
    ys = [c[1] for c in all_coords]
    zs = [c[2] for c in all_coords]
    return (
        (min(xs) - padding, min(ys) - padding, min(zs) - padding),
        (max(xs) + padding, max(ys) + padding, max(zs) + padding),
    )


def _compute_spatial_entropy_numpy(
    pose_coords_list: List[List[Tuple[float, float, float]]],
    bounds: Tuple[Tuple[float, float, float], Tuple[float, float, float]],
    grid_spacing: float,
    sigma: float,
    max_points: int = 30,
) -> List[Tuple[float, float, float, float]]:
    """NumPy-accelerated spatial entropy computation.

    Vectorises the Gaussian density evaluation using broadcasting,
    giving ~10-100x speedup over pure Python on typical grids.
    """
    (x_min, y_min, z_min), (x_max, y_max, z_max) = bounds

    nx = min(max(1, int((x_max - x_min) / grid_spacing) + 1), max_points)
    ny = min(max(1, int((y_max - y_min) / grid_spacing) + 1), max_points)
    nz = min(max(1, int((z_max - z_min) / grid_spacing) + 1), max_points)

    gx = np.linspace(x_min, x_min + (nx - 1) * grid_spacing, nx)
    gy = np.linspace(y_min, y_min + (ny - 1) * grid_spacing, ny)
    gz = np.linspace(z_min, z_min + (nz - 1) * grid_spacing, nz)

    # Grid points: (N_grid, 3) where N_grid = nx*ny*nz
    grid_x, grid_y, grid_z = np.meshgrid(gx, gy, gz, indexing="ij")
    grid_pts = np.stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()], axis=1)
    n_grid = grid_pts.shape[0]

    inv_2sigma2 = 1.0 / (2.0 * sigma * sigma)
    r2_cutoff = (4.0 * sigma) ** 2
    n_poses = len(pose_coords_list)

    # Accumulate per-pose densities: shape (n_grid, n_poses)
    densities = np.zeros((n_grid, n_poses), dtype=np.float64)

    for p_idx, pose_coords in enumerate(pose_coords_list):
        atoms = np.asarray(pose_coords, dtype=np.float64)  # (n_atoms, 3)
        # Process in chunks to limit memory: each chunk evaluates
        # (n_grid, chunk_size) distances
        chunk_size = max(1, min(len(atoms), 500))
        for start in range(0, len(atoms), chunk_size):
            atom_chunk = atoms[start:start + chunk_size]  # (chunk, 3)
            # Broadcast: (n_grid, 1, 3) - (1, chunk, 3) -> (n_grid, chunk, 3)
            diff = grid_pts[:, np.newaxis, :] - atom_chunk[np.newaxis, :, :]
            r2 = np.sum(diff * diff, axis=2)  # (n_grid, chunk)
            mask = r2 < r2_cutoff
            contrib = np.where(mask, np.exp(-r2 * inv_2sigma2), 0.0)
            densities[:, p_idx] += contrib.sum(axis=1)

    # Shannon entropy per grid point
    totals = densities.sum(axis=1)  # (n_grid,)
    valid = totals > 0.1

    results = []
    if not np.any(valid):
        return results

    # Normalise to probabilities
    probs = densities[valid] / totals[valid, np.newaxis]  # (n_valid, n_poses)
    # -p * ln(p), with 0*ln(0) = 0
    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(probs > 1e-15, np.log(probs), 0.0)
    entropy = -np.sum(probs * log_p, axis=1)

    max_s = math.log(n_poses) if n_poses > 1 else 1.0
    entropy_norm = entropy / max_s if max_s > 0 else entropy

    valid_pts = grid_pts[valid]
    for i in range(len(valid_pts)):
        results.append((
            float(valid_pts[i, 0]),
            float(valid_pts[i, 1]),
            float(valid_pts[i, 2]),
            float(entropy_norm[i]),
        ))

    return results


def _compute_spatial_entropy_pure(
    pose_coords_list: List[List[Tuple[float, float, float]]],
    bounds: Tuple[Tuple[float, float, float], Tuple[float, float, float]],
    grid_spacing: float,
    sigma: float,
    max_points: int = 30,
) -> List[Tuple[float, float, float, float]]:
    """Pure-Python spatial entropy computation (fallback when NumPy unavailable)."""
    (x_min, y_min, z_min), (x_max, y_max, z_max) = bounds

    nx = min(max(1, int((x_max - x_min) / grid_spacing) + 1), max_points)
    ny = min(max(1, int((y_max - y_min) / grid_spacing) + 1), max_points)
    nz = min(max(1, int((z_max - z_min) / grid_spacing) + 1), max_points)

    inv_2sigma2 = 1.0 / (2.0 * sigma * sigma)
    n_poses = len(pose_coords_list)
    r2_cutoff = (4.0 * sigma) ** 2

    results = []

    for ix in range(nx):
        gx = x_min + ix * grid_spacing
        for iy in range(ny):
            gy = y_min + iy * grid_spacing
            for iz in range(nz):
                gz = z_min + iz * grid_spacing

                pose_densities = []
                for pose_coords in pose_coords_list:
                    density = 0.0
                    for ax, ay, az in pose_coords:
                        dx = gx - ax
                        dy = gy - ay
                        dz = gz - az
                        r2 = dx * dx + dy * dy + dz * dz
                        if r2 < r2_cutoff:
                            density += math.exp(-r2 * inv_2sigma2)
                    pose_densities.append(density)

                total = sum(pose_densities)
                if total < 1e-12:
                    continue

                entropy = 0.0
                for d in pose_densities:
                    p = d / total
                    if p > 1e-15:
                        entropy -= p * math.log(p)

                max_s = math.log(n_poses) if n_poses > 1 else 1.0
                entropy_norm = entropy / max_s if max_s > 0 else 0.0

                if total > 0.1:
                    results.append((gx, gy, gz, entropy_norm))

    return results


def _compute_spatial_entropy(
    mode: BindingModeResult,
    grid_spacing: float = 2.0,
    sigma: float = 2.0,
) -> List[Tuple[float, float, float, float]]:
    """Compute spatial entropy density on a 3D grid from pose ensemble.

    Each grid point accumulates a Gaussian-weighted atom count from all
    poses.  The local density is then converted to a Shannon-style
    entropy measure:  S_i = -sum_j p_j ln(p_j)  where p_j is the
    fraction of total density contributed by pose j at grid point i.

    Uses NumPy vectorisation when available for ~10-100x speedup.

    Args:
        mode: BindingModeResult with pose PDB paths.
        grid_spacing: Distance between grid points in Angstroms.
        sigma: Gaussian smoothing width in Angstroms.

    Returns:
        List of (x, y, z, entropy) tuples for grid points with non-zero
        entropy.
    """
    if not mode.poses:
        return []

    pose_coords_list = []
    for pose in mode.poses:
        coords = _read_pose_coords(str(pose.path))
        if coords:
            pose_coords_list.append(coords)

    if not pose_coords_list:
        return []

    all_coords = [c for pose_coords in pose_coords_list for c in pose_coords]
    if not all_coords:
        return []

    bounds = _compute_grid_bounds(all_coords)

    if _HAS_NUMPY:
        return _compute_spatial_entropy_numpy(
            pose_coords_list, bounds, grid_spacing, sigma,
        )
    return _compute_spatial_entropy_pure(
        pose_coords_list, bounds, grid_spacing, sigma,
    )


def _entropy_color(entropy_norm: float) -> Tuple[float, float, float]:
    """Map normalised entropy [0,1] to a burgundy-white-purple colormap.

    0.0 (low entropy / ordered) -> burgundy red (0.502, 0.0, 0.125)
    0.5 (medium)                -> white
    1.0 (high entropy / disordered) -> purple blue (0.294, 0.0, 0.510)
    """
    t = max(0.0, min(1.0, entropy_norm))
    if t < 0.5:
        s = t / 0.5
        # burgundy red → white
        return (0.502 + s * (1.0 - 0.502), s, 0.125 + s * (1.0 - 0.125))
    else:
        s = (t - 0.5) / 0.5
        # white → purple blue
        return (1.0 + s * (0.294 - 1.0), 1.0 - s, 1.0 + s * (0.510 - 1.0))


def _render_pseudoatom(
    grid_points: List[Tuple[float, float, float, float]],
    obj_name: str,
    mode_id: int,
    sphere_scale: float,
    transparency: float,
) -> None:
    """Render heatmap using pseudoatom spheres (legacy method)."""
    for i, (x, y, z, s_norm) in enumerate(grid_points):
        r, g, b = _entropy_color(s_norm)
        color_name = f"_ent_m{mode_id}_{i}"
        cmd.set_color(color_name, [r, g, b])
        cmd.pseudoatom(
            obj_name,
            pos=[x, y, z],
            name=f"E{i}",
            b=s_norm,
        )

    cmd.show("spheres", obj_name)
    cmd.set("sphere_scale", sphere_scale, obj_name)
    cmd.set("sphere_transparency", transparency, obj_name)

    for i, (x, y, z, s_norm) in enumerate(grid_points):
        color_name = f"_ent_m{mode_id}_{i}"
        cmd.color(color_name, f"{obj_name} and name E{i}")


def _render_cgo(
    grid_points: List[Tuple[float, float, float, float]],
    obj_name: str,
    sphere_scale: float,
    transparency: float,
) -> None:
    """Render heatmap using PyMOL CGO spheres for better performance.

    CGO (Compiled Graphics Objects) bypass the molecular object overhead,
    rendering faster with large point counts and supporting transparency
    natively.
    """
    # CGO constants
    COLOR = 6.0    # 0x6
    SPHERE = 7.0   # 0x7
    ALPHA = 25.0   # 0x19

    cgo_list: list = [ALPHA, 1.0 - transparency]

    for x, y, z, s_norm in grid_points:
        r, g, b = _entropy_color(s_norm)
        cgo_list.extend([COLOR, r, g, b])
        cgo_list.extend([SPHERE, x, y, z, sphere_scale])

    cmd.load_cgo(cgo_list, obj_name)


def render_entropy_heatmap(
    mode_id: int,
    grid_spacing: float = 2.0,
    sigma: float = 2.0,
    sphere_scale: float = 0.4,
    transparency: float = 0.3,
    renderer: str = "cgo",
) -> None:
    """Compute and render an entropy heatmap for a loaded binding mode.

    Requires results to be loaded first via ``flexaids_load_results``.

    Args:
        mode_id: Numeric binding-mode identifier.
        grid_spacing: Grid spacing in Angstroms (default 2.0).
        sigma: Gaussian smoothing width in Angstroms (default 2.0).
        sphere_scale: PyMOL sphere scale for heatmap points.
        transparency: Sphere transparency (0=opaque, 1=invisible).
        renderer: 'cgo' (default, faster) or 'pseudoatom' (legacy).

    Example:
        PyMOL> flexaids_entropy_heatmap 1
        PyMOL> flexaids_entropy_heatmap 1, renderer=pseudoatom
        PyMOL> flexaids_entropy_heatmap 1, grid_spacing=1.5, sigma=1.5
    """
    from . import results_adapter

    mode = results_adapter._get_mode(int(mode_id))
    if mode is None:
        print("ERROR: No loaded result set or mode not found. "
              "Use 'flexaids_load_results' first.")
        return

    if not mode.poses:
        print(f"ERROR: Mode {mode_id} has no poses.")
        return

    grid_spacing = float(grid_spacing)
    sigma = float(sigma)
    renderer = str(renderer).strip().lower()

    print(f"Computing entropy heatmap for mode {mode_id} "
          f"(grid={grid_spacing}A, sigma={sigma}A, renderer={renderer})...")

    grid_points = _compute_spatial_entropy(mode, grid_spacing, sigma)

    if not grid_points:
        print(f"WARNING: No entropy density computed for mode {mode_id}.")
        print("  Recovery tips:")
        print("    • Need ≥2 poses with readable ATOM/HETATM coordinates")
        print("    • Try larger sigma / coarser grid_spacing (e.g. sigma=3, grid=3)")
        print("    • Confirm PDBs are not empty or ligand-only without coords")
        print("    • Check that mode_id matches a loaded mode (flexaids_status)")
        return

    obj_name = f"entropy_map_mode{mode_id}"

    # Remove previous heatmap if it exists
    try:
        cmd.delete(obj_name)
    except Exception:
        pass

    if renderer == "pseudoatom":
        _render_pseudoatom(grid_points, obj_name, mode_id, sphere_scale, transparency)
    else:
        _render_cgo(grid_points, obj_name, sphere_scale, transparency)

    n_high = sum(1 for _, _, _, s in grid_points if s > 0.5)
    print(
        f"Entropy heatmap rendered: {len(grid_points)} grid points "
        f"({n_high} high-entropy). Object: '{obj_name}'"
    )
    print("  Burgundy = ordered (low entropy), Purple = disordered (high entropy)")
