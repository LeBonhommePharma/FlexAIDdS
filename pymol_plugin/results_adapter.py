"""Read-only PyMOL adapter for FlexAID∆S docking results.

Loads result directories through ``flexaidds.load_results`` and creates
grouped PyMOL objects using the shared :mod:`pymol_plugin.session` state.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

try:
    from pymol import cmd
except ImportError as exc:
    raise ImportError("PyMOL not available") from exc

try:
    from flexaidds import BindingModeResult, DockingResult, load_results
except ImportError as exc:
    raise ImportError(
        "flexaidds Python package is required for the read-only PyMOL adapter"
    ) from exc

from .session import SESSION, group_name, object_name


def _safe_int(value, name: str = "value") -> Optional[int]:
    """Convert to int with a user-friendly error message."""
    try:
        return int(value)
    except (ValueError, TypeError):
        print(f"ERROR: '{value}' is not a valid integer for {name}.")
        return None


def _sync_compat_aliases() -> None:
    """Keep legacy module globals in sync with SESSION.

    Older callers (GUI, tests) still read/write ``_loaded_result`` etc.
    """
    global _loaded_result, _loaded_objects, _loaded_prefix
    _loaded_result = SESSION.result
    _loaded_objects = SESSION.objects
    _loaded_prefix = SESSION.prefix


# Legacy module-level aliases (kept for GUI / external scripts)
_loaded_result: Optional[DockingResult] = None
_loaded_objects: Dict[int, List[str]] = {}
_loaded_prefix: str = "flexaids"
_sync_compat_aliases()


def _get_mode(mode_id: int) -> Optional[BindingModeResult]:
    # Prefer SESSION; fall back to legacy alias if an external caller
    # assigned _loaded_result without going through load_docking_results.
    if SESSION.result is None and _loaded_result is not None:
        SESSION.result = _loaded_result
    return SESSION.get_mode(mode_id)


def _delete_previous_objects(prefix: Optional[str] = None) -> None:
    """Remove previously created FlexAID∆S groups/objects from PyMOL."""
    pref = prefix or SESSION.prefix
    for names in list(SESSION.objects.values()):
        for name in names:
            try:
                cmd.delete(name)
            except Exception:
                pass
    for mode_id in list(SESSION.objects.keys()):
        try:
            cmd.delete(group_name(pref, mode_id))
        except Exception:
            pass
    try:
        for obj in list(cmd.get_object_list("all") or []):
            s = str(obj)
            if s.startswith(f"{pref}_mode") or s.startswith(f"{pref}_"):
                try:
                    cmd.delete(obj)
                except Exception:
                    pass
    except Exception:
        pass


def unload_results(delete_objects: int = 1) -> None:
    """Clear the loaded docking result session.

    Safe unload: deletes only ``{prefix}_*`` objects plus heatmaps/morphs
    created by the plugin.

    Args:
        delete_objects: 1 (default) also deletes PyMOL objects created by
            the last load; 0 only clears the in-memory session.
    """
    pref = SESSION.prefix or "flexaids"
    if int(delete_objects):
        _delete_previous_objects(pref)
        # Heatmaps, morphs, contact maps
        try:
            for obj in list(cmd.get_object_list("all") or []):
                s = str(obj)
                if (
                    s.startswith("entropy_map_mode")
                    or s.startswith(f"{pref}_contacts_")
                    or s.startswith(f"{pref}_morph")
                    or s.startswith(f"{pref}_cleft")
                    or "flexaids_morph" in s
                ):
                    try:
                        cmd.delete(obj)
                    except Exception:
                        pass
        except Exception:
            pass
    SESSION.clear()
    _sync_compat_aliases()
    try:
        from . import visualization
        visualization.clear_session_view()
    except Exception:
        pass
    print("FlexAID∆S results unloaded (prefix-scoped cleanup).")


def load_docking_results(
    results_dir: str,
    prefix: str = "flexaids",
    temperature: float = 300.0,
) -> None:
    """Load FlexAID∆S results — alias of the unified ``flexaids_load`` path.

    Deprecated as a distinct entrypoint; prefer ``flexaids_load``.
    Kept for back-compat with scripts and GUI code that still call it.

    Args:
        results_dir: Directory containing docking result PDB files.
        prefix: Prefix used to create PyMOL object and group names.
        temperature: Fallback temperature (K) if REMARKs omit T.

    Example:
        PyMOL> flexaids_load_results /path/to/output
        PyMOL> flexaids_load /path/to/output [, prefix] [, temperature]
    """
    from . import visualization

    # Delegate fully to the canonical loader (same object model + readiness card)
    visualization.load_binding_modes(
        results_dir,
        prefix=prefix,
        temperature=temperature,
    )
    _sync_compat_aliases()


def show_binding_mode(mode_id: int, show_all: int = 0) -> None:
    """Show one loaded binding mode.

    Args:
        mode_id: Numeric binding-mode identifier.
        show_all: 1 to show all poses in the mode, 0 for best pose only.
    """
    mid = _safe_int(mode_id, "mode_id")
    if mid is None:
        return
    # Re-sync if GUI/tests assigned aliases directly
    if SESSION.result is None and _loaded_result is not None:
        SESSION.result = _loaded_result
        SESSION.objects = _loaded_objects or {}
        SESSION.prefix = _loaded_prefix or "flexaids"

    mode = _get_mode(mid)
    if mode is None:
        print("ERROR: No loaded result set or mode not found.")
        return

    object_names = SESSION.object_names_for(mode.mode_id)
    if not object_names and _loaded_objects:
        object_names = list(_loaded_objects.get(mode.mode_id, []))
    if not object_names:
        print(f"ERROR: No PyMOL objects loaded for mode {mode.mode_id}.")
        return

    best_pose = mode.best_pose()
    best_name = None
    if best_pose is not None:
        best_name = object_name(SESSION.prefix, mode.mode_id, best_pose.pose_rank)

    for names in (SESSION.objects or _loaded_objects).values():
        for name in names:
            cmd.disable(name)

    if int(show_all):
        for name in object_names:
            cmd.enable(name)
    elif best_name is not None and best_name in object_names:
        cmd.enable(best_name)
    else:
        cmd.enable(object_names[0])

    gname = group_name(SESSION.prefix, mode.mode_id)
    try:
        cmd.zoom(gname)
    except Exception:
        cmd.zoom(" ".join(object_names))

    print(
        f"Mode {mode.mode_id}: n_poses={mode.n_poses}, "
        f"free_energy={mode.free_energy}, best_cf={mode.best_cf}"
    )


def color_mode_by_score(mode_id: int, metric: str = "cf") -> None:
    """Color poses within one mode using a score gradient.

    Args:
        mode_id: Numeric binding-mode identifier.
        metric: 'cf' (CF/contact-function scoring proxy) or 'free_energy'.
            Lower values are colored burgundy red.
    """
    mid = _safe_int(mode_id, "mode_id")
    if mid is None:
        return
    if SESSION.result is None and _loaded_result is not None:
        SESSION.result = _loaded_result
        SESSION.objects = _loaded_objects or {}
        SESSION.prefix = _loaded_prefix or "flexaids"

    mode = _get_mode(mid)
    if mode is None:
        print("ERROR: No loaded result set or mode not found.")
        return

    object_names = SESSION.object_names_for(mode.mode_id) or list(
        (_loaded_objects or {}).get(mode.mode_id, [])
    )
    if not object_names:
        print(f"ERROR: No PyMOL objects loaded for mode {mode.mode_id}.")
        return

    metric = str(metric).strip().lower()
    pose_by_name: Dict[str, object] = {}
    for pose in mode.poses:
        pose_by_name[object_name(SESSION.prefix, mode.mode_id, pose.pose_rank)] = pose

    values = []
    ordered_names = []
    for obj_name in object_names:
        pose = pose_by_name.get(obj_name)
        if pose is None:
            continue
        if metric == "free_energy":
            value = pose.free_energy if pose.free_energy is not None else mode.free_energy
        else:
            value = pose.cf if pose.cf is not None else mode.best_cf
        values.append(value)
        ordered_names.append(obj_name)

    finite = [v for v in values if v is not None]
    if not finite:
        print(f"ERROR: No numeric values available for metric '{metric}'.")
        return

    vmin = min(finite)
    vmax = max(finite)
    vrange = (vmax - vmin) if vmax > vmin else 1.0

    for obj_name, value in zip(ordered_names, values):
        if value is None:
            continue
        t = max(0.0, min(1.0, (value - vmin) / vrange))
        color_name = f"{SESSION.prefix}_{metric}_{obj_name}"
        r = 0.502 + t * (0.294 - 0.502)
        g = 0.0
        b = 0.125 + t * (0.510 - 0.125)
        cmd.set_color(color_name, [r, g, b])
        cmd.color(color_name, obj_name)
        cmd.enable(obj_name)

    label = "CF scoring proxy" if metric == "cf" else metric
    print(
        f"Colored mode {mode.mode_id} by {label} "
        "(burgundy=lower/favourable, purple=higher)."
    )


def show_mode_details(mode_id: int) -> None:
    """Print thermodynamic summary for one loaded mode."""
    mid = _safe_int(mode_id, "mode_id")
    if mid is None:
        return
    if SESSION.result is None and _loaded_result is not None:
        SESSION.result = _loaded_result

    mode = _get_mode(mid)
    if mode is None:
        print("ERROR: No loaded result set or mode not found.")
        return

    temperature = mode.temperature if mode.temperature is not None else (
        SESSION.result.temperature if SESSION.result is not None else SESSION.temperature_K
    )
    entropy_term = None
    if mode.entropy is not None and temperature is not None:
        entropy_term = mode.entropy * temperature

    print(f"Mode {mode.mode_id} details (ensemble thermodynamic ledger)")
    print(f"  rank:          {mode.rank}")
    print(f"  poses:         {mode.n_poses}")
    print(f"  best_cf:       {mode.best_cf}  (CF/contact-function scoring proxy)")
    print(f"  free_energy F: {mode.free_energy}  kcal/mol")
    print(f"  enthalpy H:    {mode.enthalpy}  kcal/mol")
    print(f"  entropy S:     {mode.entropy}  kcal/(mol·K)")
    print(f"  temperature:   {temperature}  K")
    print(f"  T·S:           {entropy_term}  kcal/mol")
    if mode.heat_capacity is not None:
        print(f"  heat capacity: {mode.heat_capacity}  kcal/(mol·K²)")
