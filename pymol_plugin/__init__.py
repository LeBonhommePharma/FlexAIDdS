"""FlexAID∆S PyMOL Plugin.

Visualization and analysis interface for molecular docking results.

Features:
- Binding mode cluster rendering
- Thermodynamic property display (ensemble ledger)
- Pose ensemble visualization with Boltzmann weighting
- Integration with NRGSuite workflow
- Read-only loading of docking result ensembles through the flexaidds Python API
- Entropy heatmap visualization (spatial entropy density)
- Interactive docking workflow from PyMOL
- Binding mode animation (coordinate interpolation morph)
- ITC-style thermogram comparison plots

Installation:
    1. PyMOL > Plugin Manager > Install New Plugin
    2. Select this directory or ZIP file
    3. Restart PyMOL

Usage:
    PyMOL> Plugin > FlexAID∆S
    PyMOL> flexaids_help
"""

from __future__ import annotations

try:
    from pymol import cmd

    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    import warnings

    warnings.warn("PyMOL not available. Plugin functionality disabled.", ImportWarning)

__version__ = "2.1.0"
__author__ = "Louis-Philippe Morency"
__email__ = "lp@thebonhomme.com"

_COMMANDS_REGISTERED = False


def __init_plugin__(app=None):
    """PyMOL plugin initialization (called automatically by PyMOL)."""
    if not PYMOL_AVAILABLE:
        print("FlexAID∆S Plugin: PyMOL not available")
        return

    try:
        from pymol.plugins import addmenuitemqt

        addmenuitemqt("FlexAID∆S", run_plugin_gui)
    except Exception as exc:
        print(f"FlexAID∆S Plugin: could not register menu item: {exc}")

    _register_commands()


def run_plugin_gui():
    """Launch the main FlexAID∆S GUI panel."""
    if not PYMOL_AVAILABLE:
        print("ERROR: PyMOL not available")
        return

    try:
        from .gui import FlexAIDSPanel
    except ImportError as exc:
        print(f"ERROR: Could not load FlexAID∆S GUI: {exc}")
        return

    dialog = FlexAIDSPanel()
    dialog.show()
    # Keep a reference so the dialog is not garbage-collected
    run_plugin_gui._dialog = dialog  # type: ignore[attr-defined]


def flexaids_help(*_args, **_kwargs) -> None:
    """Print available FlexAID∆S PyMOL commands."""
    print(
        f"""
FlexAID∆S PyMOL Plugin v{__version__}
====================================
Load / inspect
  flexaids_load <dir> [, temperature]
  flexaids_load_results <dir> [, prefix]
  flexaids_unload [, delete_objects]
  flexaids_show_mode <mode_id> [, show_all]
  flexaids_show_ensemble <modeN|id> [, show_all]
  flexaids_mode_details <mode_id>
  flexaids_thermo <modeN|id>

Color / viz
  flexaids_color_boltzmann <modeN|id>
  flexaids_color_mode <mode_id> [, metric=cf|free_energy]
  flexaids_entropy_heatmap <mode_id> [, grid_spacing, sigma, renderer]

Analysis
  flexaids_animate <mode_a>, <mode_b> [, n_frames, align]
  flexaids_itc_plot [, output_png]
  flexaids_itc_compare <itc.csv> [, output_png]

Docking
  flexaids_dock <receptor_obj>, <ligand.mol2> [, site_selection, temperature]
  flexaids_dock_cancel

Notes
  - Ranking uses the CF/contact-function scoring proxy during GA search.
  - Free energy / entropy shown here are ensemble thermodynamic ledger values
    from StatMech when present in PDB REMARKs.
"""
    )


def _register_commands() -> None:
    """Register PyMOL cmd.extend entry points (idempotent, lazy imports)."""
    global _COMMANDS_REGISTERED
    if not PYMOL_AVAILABLE or _COMMANDS_REGISTERED:
        return

    try:
        from .visualization import (
            load_binding_modes,
            show_pose_ensemble,
            color_by_boltzmann_weight,
            show_thermodynamics,
        )
        from .results_adapter import (
            load_docking_results,
            show_binding_mode,
            color_mode_by_score,
            show_mode_details,
            unload_results,
        )
        from .entropy_heatmap import render_entropy_heatmap
        from .mode_animation import animate_binding_modes
        from .itc_comparison import (
            plot_enthalpy_entropy_compensation,
            plot_free_energy_comparison,
        )
        from .interactive_docking import dock_interactive, dock_cancel
    except ImportError as exc:
        print(f"FlexAID∆S Plugin: command registration incomplete: {exc}")
        print("  Install the flexaidds package (pip install -e ./python) and retry.")
        return

    cmd.extend("flexaids_load", load_binding_modes)
    cmd.extend("flexaids_show_ensemble", show_pose_ensemble)
    cmd.extend("flexaids_color_boltzmann", color_by_boltzmann_weight)
    cmd.extend("flexaids_thermo", show_thermodynamics)
    cmd.extend("flexaids_load_results", load_docking_results)
    cmd.extend("flexaids_show_mode", show_binding_mode)
    cmd.extend("flexaids_color_mode", color_mode_by_score)
    cmd.extend("flexaids_mode_details", show_mode_details)
    cmd.extend("flexaids_unload", unload_results)
    cmd.extend("flexaids_entropy_heatmap", render_entropy_heatmap)
    cmd.extend("flexaids_animate", animate_binding_modes)
    cmd.extend("flexaids_itc_plot", plot_enthalpy_entropy_compensation)
    cmd.extend("flexaids_itc_compare", plot_free_energy_comparison)
    cmd.extend("flexaids_dock", dock_interactive)
    cmd.extend("flexaids_dock_cancel", dock_cancel)
    cmd.extend("flexaids_help", flexaids_help)

    _COMMANDS_REGISTERED = True
    print(f"FlexAID∆S Plugin v{__version__} ready. Type 'flexaids_help' for commands.")


# Auto-register when imported inside PyMOL (Plugin Manager also calls __init_plugin__)
if PYMOL_AVAILABLE:
    _register_commands()
