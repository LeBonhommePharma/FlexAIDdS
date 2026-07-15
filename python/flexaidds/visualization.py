"""PyMOL integration for FlexAID∆S pose visualization.

Provides helpers for rendering docked poses, binding modes, and thermodynamic
heat maps within PyMOL.  All public functions gracefully degrade to a warning
when PyMOL is not installed.
"""

import warnings
from typing import List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .docking import BindingMode, BindingPopulation, Pose
    from .io import PDBStructure

try:
    import pymol
    from pymol import cmd as _cmd
    _PYMOL_AVAILABLE = True
except ImportError:
    _PYMOL_AVAILABLE = False


def _require_pymol() -> None:
    if not _PYMOL_AVAILABLE:
        raise ImportError(
            "PyMOL is required for visualization. "
            "Install it with: pip install pymol-open-source"
        )


def load_binding_mode(
    mode: "BindingMode",
    pdb_paths: List[str],
    mode_name: str = "binding_mode",
    *,
    color: str = "cyan",
    show: str = "sticks",
) -> None:
    """Load all poses of a binding mode into PyMOL.

    Args:
        mode: BindingMode whose poses will be visualized.
        pdb_paths: List of PDB file paths for the poses in this mode.
        mode_name: Object name prefix in PyMOL.
        color: PyMOL color string for the ligand poses.
        show: PyMOL representation ('sticks', 'lines', 'spheres', etc.).
    """
    _require_pymol()
    for i, pdb_path in enumerate(pdb_paths):
        obj = f"{mode_name}_{i:03d}"
        _cmd.load(pdb_path, obj)
        _cmd.color(color, obj)
        _cmd.show(show, obj)


def load_population(
    population: "BindingPopulation",
    pdb_dir: str,
    *,
    palette: Optional[List[str]] = None,
) -> None:
    """Load a full BindingPopulation into PyMOL with per-mode colors.

    Args:
        population: BindingPopulation to visualize.
        pdb_dir: Directory containing output PDB files.
        palette: List of PyMOL color names, one per mode.  Cycles if fewer
                 colors than modes.
    """
    _require_pymol()
    from pathlib import Path

    # Burgundy-to-purple spectrum palette
    _cmd.set_color("_bp_0", [0.502, 0.000, 0.125])  # burgundy red
    _cmd.set_color("_bp_1", [0.470, 0.000, 0.200])
    _cmd.set_color("_bp_2", [0.435, 0.000, 0.270])
    _cmd.set_color("_bp_3", [0.400, 0.000, 0.340])
    _cmd.set_color("_bp_4", [0.365, 0.000, 0.400])
    _cmd.set_color("_bp_5", [0.340, 0.000, 0.440])
    _cmd.set_color("_bp_6", [0.294, 0.000, 0.510])  # purple blue
    default_palette = ["_bp_0", "_bp_1", "_bp_2", "_bp_3",
                       "_bp_4", "_bp_5", "_bp_6"]
    palette = palette or default_palette
    pdb_dir_path = Path(pdb_dir)

    # Canonical naming aligns with pymol_plugin: flexaids_mode{N}_pose{R}
    prefix = "flexaids"
    for mode_idx, mode in enumerate(population):
        color = palette[mode_idx % len(palette)]
        mode_id = mode_idx + 1
        pattern = f"*_{mode_id}_*.pdb"
        mode_pdbs = sorted(pdb_dir_path.glob(pattern))
        if not mode_pdbs:
            # Fallback: any PDBs with mode token
            mode_pdbs = sorted(pdb_dir_path.glob(f"*mode*{mode_id}*.pdb"))
        if not mode_pdbs:
            warnings.warn(
                f"No PDB files found for binding mode {mode_id} "
                f"(pattern: {pattern})", RuntimeWarning
            )
            continue
        for pose_rank, pdb_path in enumerate(mode_pdbs, start=1):
            obj = f"{prefix}_mode{mode_id}_pose{pose_rank}"
            _cmd.load(str(pdb_path), obj)
            _cmd.color(color, obj)
            _cmd.show("sticks", f"{obj} and organic")
            _cmd.show("cartoon", f"{obj} and polymer")
            try:
                _cmd.group(f"{prefix}_mode{mode_id}", obj)
            except Exception:
                pass


def _burgundy_purple_rgb(frac: float):
    """Interpolate from burgundy red to purple blue.

    frac = 0.0 → burgundy red (0.502, 0.0, 0.125)
    frac = 1.0 → purple blue (0.294, 0.0, 0.510)
    """
    frac = max(0.0, min(1.0, frac))
    r = 0.502 + frac * (0.294 - 0.502)
    g = 0.0
    b = 0.125 + frac * (0.510 - 0.125)
    return [r, g, b]


def color_by_energy(
    obj_name: str,
    energies: List[float],
    *,
    low_color: str = "blue",
    high_color: str = "red",
) -> None:
    """Color PyMOL object states by energy using a burgundy-to-purple gradient.

    Args:
        obj_name: PyMOL object name.
        energies: Energy value for each state (kcal/mol).
        low_color: Ignored (kept for API compat). Uses burgundy red.
        high_color: Ignored (kept for API compat). Uses purple blue.
    """
    _require_pymol()
    if not energies:
        return

    e_min = min(energies)
    e_max = max(energies)
    e_range = e_max - e_min or 1.0

    for state_idx, energy in enumerate(energies, start=1):
        frac = (energy - e_min) / e_range  # 0 = most favorable, 1 = least
        # Burgundy red → purple blue spectrum
        rgb = _burgundy_purple_rgb(frac)
        _cmd.set_color(f"_energy_{state_idx}", rgb)
        _cmd.color(f"_energy_{state_idx}", f"{obj_name} and state {state_idx}")


def color_by_boltzmann_weight(
    obj_name: str,
    weights: List[float],
) -> None:
    """Color PyMOL object states by Boltzmann weight (burgundy-to-purple).

    High probability (weight ~ 1.0) → burgundy red.
    Low probability (weight ~ 0.0) → purple blue.

    Args:
        obj_name: PyMOL object name.
        weights: Boltzmann weight for each state (normalized, 0-1).
    """
    _require_pymol()
    if not weights:
        return

    for state_idx, w in enumerate(weights, start=1):
        # Invert: high weight = burgundy (frac=0), low weight = purple (frac=1)
        frac = 1.0 - max(0.0, min(1.0, w))
        rgb = _burgundy_purple_rgb(frac)
        _cmd.set_color(f"_bw_{state_idx}", rgb)
        _cmd.color(f"_bw_{state_idx}", f"{obj_name} and state {state_idx}")


def show_cleft_spheres(sphere_pdb: str, *, color: str = "yellow",
                       transparency: float = 0.5) -> None:
    """Visualize a FlexAID cleft sphere file as transparent spheres in PyMOL.

    Args:
        sphere_pdb: Path to sphere PDB written by CavityDetector.
        color: Sphere color.
        transparency: 0.0 = opaque, 1.0 = invisible.
    """
    _require_pymol()
    obj = "cleft_spheres"
    _cmd.load(sphere_pdb, obj)
    _cmd.show("spheres", obj)
    _cmd.color(color, obj)
    _cmd.set("sphere_transparency", transparency, obj)


def setup_publication_view(
    receptor_obj: str = "receptor",
    ligand_obj: str = "ligand",
    *,
    bg_color: str = "white",
    ray_trace: bool = False,
    output_png: Optional[str] = None,
) -> None:
    """Apply publication-quality PyMOL settings and optionally render.

    Args:
        receptor_obj: PyMOL object name for receptor.
        ligand_obj: PyMOL object name (or selection) for ligand poses.
        bg_color: Background color.
        ray_trace: Whether to ray-trace before saving.
        output_png: If given, save rendered image to this path.
    """
    _require_pymol()
    _cmd.bg_color(bg_color)
    _cmd.set("ray_shadows", 0)
    _cmd.set("ambient", 0.4)
    _cmd.set("specular", 0.0)
    _cmd.set("depth_cue", 1)

    if receptor_obj in _cmd.get_object_list():
        _cmd.show("cartoon", receptor_obj)
        _cmd.color("grey80", receptor_obj)
        _cmd.set("cartoon_transparency", 0.3, receptor_obj)

    _cmd.orient()
    _cmd.zoom("all", buffer=2.0)

    if output_png:
        if ray_trace:
            _cmd.ray(1200, 900)
        _cmd.png(output_png, dpi=300)
