"""Binding mode animation for FlexAID∆S (Phase 3, deliverable 3.3).

Generates smooth coordinate interpolation between two binding mode
representative poses and renders as a PyMOL movie.

Usage:
    PyMOL> flexaids_load_results /path/to/output
    PyMOL> flexaids_animate 1 2
    PyMOL> flexaids_animate 1 2, n_frames=100
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
    try:
        from pymol import stored as _pymol_stored
    except ImportError:
        _pymol_stored = getattr(cmd, "stored", None)
except ImportError as exc:
    raise ImportError("PyMOL not available") from exc


def _kabsch_align(
    coords_mobile: List[List[float]],
    coords_target: List[List[float]],
) -> List[List[float]]:
    """Kabsch optimal rotation to align coords_mobile onto coords_target.

    Minimises RMSD between the two coordinate sets using SVD-based
    superposition.  Requires NumPy.  Returns the aligned mobile coordinates.

    If NumPy is unavailable or the sets are empty, returns coords_mobile
    unchanged.
    """
    if not _HAS_NUMPY:
        return coords_mobile
    n = min(len(coords_mobile), len(coords_target))
    if n < 3:
        return coords_mobile

    P = np.asarray(coords_mobile[:n], dtype=np.float64)
    Q = np.asarray(coords_target[:n], dtype=np.float64)

    # Centre both sets
    centroid_P = P.mean(axis=0)
    centroid_Q = Q.mean(axis=0)
    P_c = P - centroid_P
    Q_c = Q - centroid_Q

    # Cross-covariance matrix
    H = P_c.T @ Q_c

    U, S, Vt = np.linalg.svd(H)

    # Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, 1.0 if d >= 0 else -1.0])

    R = Vt.T @ sign_matrix @ U.T

    # Apply rotation and translation to all mobile coords
    all_P = np.asarray(coords_mobile, dtype=np.float64)
    aligned = (all_P - centroid_P) @ R.T + centroid_Q

    return aligned.tolist()


def _read_atom_coords(pdb_path: str) -> List[List[float]]:
    """Read all ATOM/HETATM coordinates from a PDB file.

    Returns a list of [x, y, z] lists preserving atom order.
    """
    coords = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return coords


def _interpolate_coords(
    coords1: List[List[float]],
    coords2: List[List[float]],
    t: float,
) -> List[List[float]]:
    """Linear interpolation between two coordinate sets.

    Args:
        coords1: Starting coordinates (list of [x, y, z]).
        coords2: Ending coordinates (list of [x, y, z]).
        t: Interpolation parameter in [0, 1].

    Returns:
        Interpolated coordinates.
    """
    n = min(len(coords1), len(coords2))
    if n == 0:
        return []
    if _HAS_NUMPY:
        a = np.asarray(coords1[:n], dtype=np.float64)
        b = np.asarray(coords2[:n], dtype=np.float64)
        return ((1.0 - t) * a + t * b).tolist()
    result = []
    for i in range(n):
        result.append([
            (1.0 - t) * coords1[i][0] + t * coords2[i][0],
            (1.0 - t) * coords1[i][1] + t * coords2[i][1],
            (1.0 - t) * coords1[i][2] + t * coords2[i][2],
        ])
    return result


def animate_binding_modes(
    mode_id_1: int,
    mode_id_2: int,
    n_frames: int = 50,
    movie_name: str = "",
    align: bool = True,
) -> None:
    """Create a smooth animation between two binding mode representative poses.

    Loads the best pose from each mode, optionally performs Kabsch RMSD
    alignment, then does linear Cartesian interpolation and sets up a
    PyMOL movie with the interpolated frames.

    Args:
        mode_id_1: First binding mode ID.
        mode_id_2: Second binding mode ID.
        n_frames: Number of interpolation frames (default 50).
        movie_name: If non-empty, save movie to this filename (e.g. "morph.mpg").
        align: If True (default), Kabsch-align mode 2 onto mode 1 before
               interpolation.  Requires NumPy.

    Example:
        PyMOL> flexaids_animate 1 2
        PyMOL> flexaids_animate 1 3, n_frames=100, align=0
    """
    from . import results_adapter

    mode1 = results_adapter._get_mode(int(mode_id_1))
    mode2 = results_adapter._get_mode(int(mode_id_2))

    if mode1 is None or mode2 is None:
        print("ERROR: Mode(s) not found. Use 'flexaids_load_results' first.")
        return

    best1 = mode1.best_pose()
    best2 = mode2.best_pose()

    if best1 is None or best2 is None:
        print("ERROR: No representative pose available for one or both modes.")
        return

    n_frames = int(n_frames)
    if n_frames < 2:
        n_frames = 2

    coords1 = _read_atom_coords(str(best1.path))
    coords2 = _read_atom_coords(str(best2.path))

    if not coords1 or not coords2:
        print("ERROR: Could not read coordinates from PDB files.")
        return

    n_atoms = min(len(coords1), len(coords2))
    if len(coords1) != len(coords2):
        print(
            f"WARNING: Atom count mismatch ({len(coords1)} vs {len(coords2)}). "
            f"Using first {n_atoms} atoms."
        )

    # Kabsch alignment: superpose coords2 onto coords1
    align = bool(int(align)) if isinstance(align, str) else bool(align)
    if align and _HAS_NUMPY:
        coords2 = _kabsch_align(coords2, coords1)
        print("  Kabsch alignment applied.")

    obj_name = f"morph_m{mode_id_1}_m{mode_id_2}"
    # Remove prior morph with the same name
    try:
        cmd.delete(obj_name)
    except Exception:
        pass

    # Load the first pose as the base object
    cmd.load(str(best1.path), obj_name, state=1)

    # Resolve PyMOL stored namespace (cmd.stored or pymol.stored)
    stored_ns = _pymol_stored if _pymol_stored is not None else getattr(cmd, "stored", None)
    if stored_ns is None:
        print("ERROR: PyMOL stored namespace unavailable; cannot build morph frames.")
        return

    # Create interpolated states
    for frame_idx in range(1, n_frames):
        t = frame_idx / (n_frames - 1)
        interp = _interpolate_coords(coords1, coords2, t)

        # Load the base PDB into a new state, then update coordinates
        cmd.load(str(best1.path), obj_name, state=frame_idx + 1)

        stored_key = f"_morph_coords_{frame_idx}"
        try:
            stored_ns.__dict__[stored_key] = iter(interp)
        except Exception:
            setattr(stored_ns, stored_key, iter(interp))

        cmd.alter_state(
            frame_idx + 1,
            obj_name,
            f"(x, y, z) = next(stored.{stored_key})",
        )

    # Set up the movie
    cmd.mset(f"1 -{n_frames}")
    cmd.show("sticks", f"{obj_name} and organic")
    cmd.show("cartoon", obj_name)
    cmd.zoom(obj_name)

    if movie_name:
        movie_name = str(movie_name)
        print(f"Saving movie to {movie_name}...")
        cmd.movie.produce(movie_name)
        print(f"Movie saved: {movie_name}")

    print(
        f"Animation created: '{obj_name}' ({n_frames} frames, "
        f"mode {mode_id_1} -> mode {mode_id_2}). "
        "Use PyMOL movie controls to play."
    )
