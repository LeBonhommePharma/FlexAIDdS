"""PyMOL visualization functions for FlexAID∆S binding modes.

These functions can be called from PyMOL command line:
    PyMOL> flexaids_load /path/to/output
    PyMOL> flexaids_show_ensemble mode1
    PyMOL> flexaids_color_boltzmann mode1
    PyMOL> flexaids_thermo mode1

This module delegates result-directory parsing to the canonical
``flexaidds.load_results()`` API and shares state via
:mod:`pymol_plugin.session` so entropy heatmaps / animation / ITC work
regardless of which load command was used.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

try:
    from pymol import cmd
    import pymol  # noqa: F401
except ImportError as exc:
    raise ImportError("PyMOL not available") from exc

try:
    from flexaidds import BindingModeResult, DockingResult, PoseResult, load_results
    from flexaidds.thermodynamics import kB_kcal as _kB_kcal
except ImportError as exc:
    raise ImportError(
        "flexaidds Python package is required for PyMOL result loading"
    ) from exc

from .session import SESSION, group_name, mode_label, object_name


@dataclass
class _ModeRecord:
    """Thermodynamic and structural data for one FlexAID binding mode."""

    mode_id: int
    pdb_objects: List[str] = field(default_factory=list)
    poses: List[PoseResult] = field(default_factory=list)
    cf_values: List[float] = field(default_factory=list)
    best_cf: Optional[float] = None
    total_cf: Optional[float] = None
    frequency: int = 0
    free_energy: Optional[float] = None
    enthalpy: Optional[float] = None
    entropy: Optional[float] = None
    heat_capacity: Optional[float] = None
    boltzmann_weights: List[float] = field(default_factory=list)


# Module-level mirrors kept for back-compat with any external scripts
_loaded_modes: Dict[str, _ModeRecord] = {}
_loaded_result: Optional[DockingResult] = None
_output_dir: Optional[Path] = None
_temperature_K: float = 300.0


def _score_value(pose: PoseResult) -> Optional[float]:
    """Prefer CF scoring proxy, then cf_app, then free energy."""
    if pose.cf is not None:
        return pose.cf
    if pose.cf_app is not None:
        return pose.cf_app
    if pose.free_energy is not None:
        return pose.free_energy
    return None


def _compute_boltzmann_weights(values: List[float], temperature_K: float) -> List[float]:
    """Log-sum-exp stable Boltzmann weights for a list of energies."""
    if not values:
        return []
    if temperature_K <= 0.0:
        temperature_K = 300.0
    beta = 1.0 / (_kB_kcal * temperature_K)
    neg_beta_e = [-beta * value for value in values]
    max_val = max(neg_beta_e)
    shifted = [math.exp(v - max_val) for v in neg_beta_e]
    total = sum(shifted)
    if total <= 0.0:
        return []
    return [value / total for value in shifted]


def _make_mode_record(mode: BindingModeResult) -> _ModeRecord:
    cf_values = [value for value in (_score_value(pose) for pose in mode.poses) if value is not None]
    temperature = mode.temperature if mode.temperature is not None else SESSION.temperature_K
    weights = _compute_boltzmann_weights(cf_values, temperature) if cf_values else []
    total_cf = sum(cf_values) if cf_values else None
    frequency = mode.frequency if mode.frequency is not None else mode.n_poses
    return _ModeRecord(
        mode_id=mode.mode_id,
        poses=list(mode.poses),
        cf_values=cf_values,
        best_cf=mode.best_cf,
        total_cf=total_cf,
        frequency=frequency,
        free_energy=mode.free_energy,
        enthalpy=mode.enthalpy,
        entropy=mode.entropy,
        heat_capacity=mode.heat_capacity,
        boltzmann_weights=weights,
    )


def _mode_name(mode_id: int) -> str:
    return mode_label(mode_id)


def _find_mode(mode_name: str) -> Optional[_ModeRecord]:
    rec = _loaded_modes.get(mode_name)
    if rec is not None:
        return rec
    # Accept bare integers: "1" → mode1
    if mode_name.isdigit():
        return _loaded_modes.get(_mode_name(int(mode_name)))
    return None


def _mirror_session_locals() -> None:
    global _loaded_modes, _loaded_result, _output_dir, _temperature_K
    _loaded_modes = SESSION.mode_records
    _loaded_result = SESSION.result
    _output_dir = SESSION.output_dir
    _temperature_K = SESSION.temperature_K


def clear_session_view() -> None:
    """Clear visualization bookkeeping (called on unload)."""
    global _loaded_modes, _loaded_result, _output_dir
    SESSION.mode_records.clear()
    _loaded_modes = SESSION.mode_records
    _loaded_result = None
    _output_dir = None


def sync_from_session() -> None:
    """Rebuild visualization mode records from the shared SESSION.

    Called after ``flexaids_load_results`` so ensemble/thermo commands work.
    """
    global _loaded_modes, _loaded_result, _output_dir, _temperature_K

    result = SESSION.result
    _loaded_result = result
    _output_dir = SESSION.output_dir
    _temperature_K = SESSION.temperature_K
    SESSION.mode_records.clear()
    _loaded_modes = SESSION.mode_records

    if result is None:
        return

    for mode in result.binding_modes:
        mode_name = _mode_name(mode.mode_id)
        record = _make_mode_record(mode)
        # Prefer already-created object names from results_adapter
        names = SESSION.object_names_for(mode.mode_id)
        if names:
            record.pdb_objects = list(names)
        else:
            for pose in mode.poses:
                record.pdb_objects.append(
                    object_name(SESSION.prefix, mode.mode_id, pose.pose_rank)
                )
        SESSION.mode_records[mode_name] = record

    _loaded_modes = SESSION.mode_records


def _print_readiness_card(result: DockingResult, temperature: float) -> None:
    """Print a one-line readiness summary after load."""
    n_modes = result.n_modes
    n_poses = sum(m.n_poses for m in result.binding_modes)
    has_cf = any(m.best_cf is not None for m in result.binding_modes)
    sources = {
        (m.metadata or {}).get("ledger_source", "missing")
        for m in result.binding_modes
    }
    if "engine_remark" in sources:
        fhs = "yes"
    elif "ensemble_estimate_from_cf" in sources:
        fhs = "recomputed"
    elif any(m.free_energy is not None for m in result.binding_modes):
        fhs = "yes"
    else:
        fhs = "missing"
    rmsd_vals = []
    for m in result.binding_modes:
        for p in m.poses:
            if p.rmsd_sym is not None:
                rmsd_vals.append(p.rmsd_sym)
            elif p.rmsd_raw is not None:
                rmsd_vals.append(p.rmsd_raw)
    rmsd_str = f"{min(rmsd_vals):.2f}–{max(rmsd_vals):.2f}" if rmsd_vals else "n/a"
    print(
        f"  readiness: modes={n_modes} poses={n_poses} "
        f"CF={'yes' if has_cf else 'no'} F/H/S={fhs} "
        f"T={temperature:.1f} RMSD={rmsd_str}"
    )
    if fhs == "recomputed":
        print(
            "  note: F/H/S are ensemble estimates from CF proxy "
            "(not full vib/solvent ledger)."
        )


def load_binding_modes(
    output_dir: str,
    prefix: str = "flexaids",
    temperature: float = 300.0,
) -> None:
    """Load FlexAID∆S docking results (canonical plugin entrypoint).

    Usage::

        flexaids_load <dir> [, prefix] [, temperature]

    Back-compat: if the second argument is numeric and no third arg is given,
    it is treated as temperature (legacy ``flexaids_load dir, T``).

    Objects: ``{prefix}_mode{N}_pose{R}`` grouped as ``{prefix}_mode{N}``.
    ``flexaids_load_results`` is an alias of this function.
    """
    global _loaded_modes, _loaded_result, _output_dir, _temperature_K

    output_path = Path(str(output_dir).strip())
    if not output_path.exists():
        print(f"ERROR: Directory not found: {output_dir}")
        return

    # PyMOL string args + back-compat: load(dir, T) vs load(dir, prefix [, T])
    # If second positional arg is purely numeric, treat it as temperature.
    prefix_s = str(prefix).strip() if prefix is not None else "flexaids"
    temperature_f = float(temperature)
    if re.fullmatch(r"[-+]?\d+(?:\.\d*)?", prefix_s):
        temperature_f = float(prefix_s)
        prefix_s = "flexaids"
    if not prefix_s:
        prefix_s = "flexaids"

    try:
        result = load_results(output_path)
    except Exception as exc:
        print(f"ERROR: Could not load docking results: {exc}")
        return

    if not result.binding_modes:
        print(f"ERROR: No binding modes found in {output_path.resolve()}.")
        return

    # Clear previous PyMOL objects via results_adapter cleanup
    try:
        from . import results_adapter
        results_adapter._delete_previous_objects(SESSION.prefix or prefix_s)
    except Exception:
        pass

    SESSION.clear()
    SESSION.result = result
    SESSION.prefix = prefix_s
    SESSION.temperature_K = temperature_f
    SESSION.output_dir = getattr(result, "source_dir", None) or output_path.resolve()
    if result.temperature is not None:
        SESSION.temperature_K = float(result.temperature)

    SESSION.mode_records.clear()
    SESSION.objects = {}

    for mode in result.binding_modes:
        mode_name = _mode_name(mode.mode_id)
        record = _make_mode_record(mode)
        gname = group_name(prefix_s, mode.mode_id)
        object_names: List[str] = []
        best_pose = mode.best_pose()

        for pose in mode.poses:
            obj_name = object_name(prefix_s, mode.mode_id, pose.pose_rank)
            try:
                cmd.load(str(pose.path), obj_name)
            except Exception as exc:
                print(f"WARNING: Failed to load {pose.path.name}: {exc}")
                continue
            cmd.group(gname, obj_name)
            cmd.hide("everything", obj_name)
            cmd.show("sticks", f"{obj_name} and organic")
            cmd.show("cartoon", f"{obj_name} and polymer")
            cmd.show("lines", f"{obj_name} and organic")
            if best_pose is None or pose.path != best_pose.path:
                cmd.disable(obj_name)
            object_names.append(obj_name)
            record.pdb_objects.append(obj_name)

        SESSION.objects[mode.mode_id] = object_names
        SESSION.mode_records[mode_name] = record

    _mirror_session_locals()

    # Sync results_adapter legacy aliases
    try:
        from . import results_adapter
        results_adapter._sync_compat_aliases()
    except Exception:
        pass

    n_modes = len(SESSION.mode_records)
    n_poses = sum(len(rec.pdb_objects) for rec in SESSION.mode_records.values())
    print(
        f"Loaded {n_modes} binding modes ({n_poses} PDB objects) from "
        f"{SESSION.output_dir} (prefix='{prefix_s}')"
    )
    _print_readiness_card(result, SESSION.temperature_K)
    print("Use 'flexaids_show_mode N' / 'flexaids_show_ensemble modeN' to visualize.")


def show_pose_ensemble(mode_name: str, show_all: bool = True) -> None:
    """Display all poses belonging to a binding mode."""
    if not SESSION.mode_records and SESSION.result is not None:
        sync_from_session()
    if not SESSION.mode_records:
        print("ERROR: No modes loaded. Use 'flexaids_load' or 'flexaids_load_results' first.")
        return

    # PyMOL often passes strings; coerce show_all
    if isinstance(show_all, str):
        show_all = show_all.strip().lower() not in ("0", "false", "no")

    rec = _find_mode(str(mode_name).strip())
    if rec is None:
        available = ", ".join(sorted(SESSION.mode_records))
        print(f"ERROR: Mode '{mode_name}' not found. Available: {available}")
        return

    if not rec.pdb_objects:
        print(f"ERROR: No PDB objects for {mode_name}.")
        return

    # Disable other modes first for a clean view
    for other in SESSION.mode_records.values():
        if other.mode_id == rec.mode_id:
            continue
        for obj in other.pdb_objects:
            cmd.disable(obj)

    if show_all:
        for obj in rec.pdb_objects:
            cmd.enable(obj)
            cmd.show("cartoon", f"{obj} and polymer")
            cmd.show("sticks", f"{obj} and organic")
    else:
        if rec.boltzmann_weights and len(rec.boltzmann_weights) == len(rec.pdb_objects):
            rep_idx = rec.boltzmann_weights.index(max(rec.boltzmann_weights))
        else:
            rep_idx = 0

        for index, obj in enumerate(rec.pdb_objects):
            if index == rep_idx:
                cmd.enable(obj)
                cmd.show("cartoon", f"{obj} and polymer")
                cmd.show("sticks", f"{obj} and organic")
            else:
                cmd.disable(obj)

    try:
        cmd.zoom(group_name(SESSION.prefix, rec.mode_id))
    except Exception:
        cmd.zoom(" ".join(rec.pdb_objects))
    label = "all poses" if show_all else "representative pose"
    print(f"Showing {label} for {_mode_name(rec.mode_id)} ({len(rec.pdb_objects)} PDB objects).")


def _burgundy_purple_rgb(t: float):
    """Interpolate burgundy red → purple blue.

    t = 0.0 → burgundy red (0.502, 0.0, 0.125)
    t = 1.0 → purple blue  (0.294, 0.0, 0.510)
    """
    t = max(0.0, min(1.0, t))
    r = 0.502 + t * (0.294 - 0.502)
    g = 0.0
    b = 0.125 + t * (0.510 - 0.125)
    return [r, g, b]


def color_by_boltzmann_weight(mode_name: str) -> None:
    """Color poses by Boltzmann weight (burgundy = high probability, purple = low)."""
    if not SESSION.mode_records and SESSION.result is not None:
        sync_from_session()
    if not SESSION.mode_records:
        print("ERROR: No modes loaded. Use 'flexaids_load' or 'flexaids_load_results' first.")
        return

    rec = _find_mode(str(mode_name).strip())
    if rec is None:
        available = ", ".join(sorted(SESSION.mode_records))
        print(f"ERROR: Mode '{mode_name}' not found. Available: {available}")
        return

    if not rec.pdb_objects:
        print(f"ERROR: No poses for {mode_name}.")
        return

    weights = rec.boltzmann_weights
    if not weights or len(weights) != len(rec.pdb_objects):
        n = len(rec.pdb_objects)
        weights = [1.0 / n] * n

    w_min = min(weights)
    w_max = max(weights)
    w_range = w_max - w_min if w_max > w_min else 1.0

    for index, (obj, weight) in enumerate(zip(rec.pdb_objects, weights)):
        t = (weight - w_min) / w_range
        # High weight = burgundy red (t=1 → frac=0), low weight = purple blue
        frac = 1.0 - t
        color_name = f"flexaids_bw_{_mode_name(rec.mode_id)}_{index}"
        cmd.set_color(color_name, _burgundy_purple_rgb(frac))
        cmd.color(color_name, obj)
        cmd.enable(obj)

    print(
        f"Colored {len(rec.pdb_objects)} poses for {_mode_name(rec.mode_id)} by Boltzmann weight "
        "(burgundy=high, purple=low). "
        "NOTE: weights use CF/contact-function scoring proxy energies, "
        "not ensemble free energy F — do not read as true thermodynamic p_i."
    )
    # Optional opacity by weight (high weight → more opaque)
    try:
        for obj, weight in zip(rec.pdb_objects, weights):
            # transparency: 0 = opaque, 1 = invisible
            t_alpha = max(0.0, min(0.85, 1.0 - float(weight)))
            cmd.set("stick_transparency", t_alpha, obj)
    except Exception:
        pass


def show_thermodynamics(mode_name: str) -> None:
    """Print thermodynamic properties of a binding mode to PyMOL console."""
    if not SESSION.mode_records and SESSION.result is not None:
        sync_from_session()
    if not SESSION.mode_records:
        print("ERROR: No modes loaded. Use 'flexaids_load' or 'flexaids_load_results' first.")
        return

    rec = _find_mode(str(mode_name).strip())
    if rec is None:
        available = ", ".join(sorted(SESSION.mode_records))
        print(f"ERROR: Mode '{mode_name}' not found. Available: {available}")
        return

    temperature = SESSION.temperature_K
    ledger_source = "unknown"
    if SESSION.result is not None:
        for mode in SESSION.result.binding_modes:
            if mode.mode_id == rec.mode_id:
                if mode.temperature is not None:
                    temperature = mode.temperature
                ledger_source = (mode.metadata or {}).get("ledger_source", "engine_remark")
                break
        else:
            if SESSION.result.temperature is not None:
                temperature = SESSION.result.temperature

    entropy_term = (rec.entropy * temperature) if rec.entropy is not None else None
    mname = _mode_name(rec.mode_id)

    print(f"\nThermodynamic ledger for {mname} (T = {temperature:.1f} K):")
    print(f"  ΔG / F (Free Energy): {rec.free_energy:10.4f} kcal/mol" if rec.free_energy is not None else "  ΔG / F (Free Energy): N/A")
    print(f"  ΔH / H (Enthalpy):    {rec.enthalpy:10.4f} kcal/mol" if rec.enthalpy is not None else "  ΔH / H (Enthalpy):    N/A")
    print(f"  S (Entropy):          {rec.entropy:10.6f} kcal/(mol·K)" if rec.entropy is not None else "  S (Entropy):          N/A")
    print(f"  T·S (Entropy term):   {entropy_term:10.4f} kcal/mol" if entropy_term is not None else "  T·S (Entropy term):   N/A")
    print(f"  Heat Capacity (Cv):   {rec.heat_capacity:10.4f} kcal/(mol·K²)" if rec.heat_capacity is not None else "  Heat Capacity (Cv):   N/A")
    print(f"  Best CF (proxy):      {rec.best_cf:10.5f}" if rec.best_cf is not None else "  Best CF (proxy):      N/A")
    print(f"  # Poses / frequency:  {rec.frequency:10d}")
    print(f"  ledger_source:        {ledger_source}")
    if ledger_source == "ensemble_estimate_from_cf":
        print("  note: ensemble estimate from CF proxy (not full vib/solvent ledger)")
    if rec.heat_capacity is None or rec.entropy is None:
        print("  flag: incomplete ledger (missing Cv and/or S; Gate 6 language)")
    print()


def export_to_nrgsuite(output_dir: str, nrgsuite_file: str) -> None:
    """Export binding modes to NRGSuite-compatible format."""
    if not SESSION.mode_records:
        load_binding_modes(output_dir)
        if not SESSION.mode_records:
            print(f"ERROR: Could not load any modes from {output_dir}")
            return

    out_path = Path(nrgsuite_file)
    try:
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(
                "# FlexAID∆S → NRGSuite export\n"
                "# mode_id\tbest_cf\tfree_energy_kcal_mol\t"
                "enthalpy_kcal_mol\tentropy_kcal_mol_K\tn_poses\n"
            )
            for mode_name in sorted(
                SESSION.mode_records, key=lambda name: SESSION.mode_records[name].mode_id
            ):
                rec = SESSION.mode_records[mode_name]
                if None not in (rec.free_energy, rec.enthalpy, rec.entropy):
                    fh.write(
                        f"{rec.mode_id}\t{rec.best_cf:.5f}\t{rec.free_energy:.4f}\t"
                        f"{rec.enthalpy:.4f}\t{rec.entropy:.6f}\t{rec.frequency}\n"
                    )
                else:
                    best_cf = f"{rec.best_cf:.5f}" if rec.best_cf is not None else "N/A"
                    fh.write(f"{rec.mode_id}\t{best_cf}\tN/A\tN/A\tN/A\t{rec.frequency}\n")
    except OSError as exc:
        print(f"ERROR: Could not write NRGSuite file: {exc}")
        return

    print(f"Exported {len(SESSION.mode_records)} binding modes to {out_path}")
