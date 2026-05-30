"""
JSON-driven reporting and visualization helpers for FlexAID∆S (Task 9).

All functions in this module consume only `DockingResult` / `BindingModeResult`
objects (or their JSON equivalents). They must never modify scoring, ranking,
or docking behaviour.

Experimental fields are explicitly labelled in all generated output.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .models import DockingResult, BindingModeResult


def _fmt(val: Optional[float], digits: int = 2) -> str:
    if val is None:
        return "N/A"
    return f"{val:.{digits}f}"


def generate_pymol_script(
    result: Union[DockingResult, Dict[str, Any]],
    output_path: Union[str, Path],
    *,
    receptor_name: str = "receptor",
    ligand_name: str = "ligand",
    show_experimental: bool = True,
) -> str:
    """
    Generate a self-contained PyMOL script that visualises the docking result.

    The script:
    - Loads receptor and ligand
    - Shows the top binding modes
    - Labels G_total, H_eff, -T*S, rank
    - Clearly marks experimental fields
    - Does NOT draw fake interactions

    Returns the script as a string (also written to output_path).
    """
    if isinstance(result, dict):
        # Allow raw dict / JSON-loaded data
        modes = result.get("binding_modes", [])
        source = result.get("source_dir", ".")
    else:
        modes = result.binding_modes
        source = str(result.source_dir)

    lines: List[str] = [
        "# FlexAIDdS Thermodynamic Report - PyMOL Script",
        "# Generated from JSON results only (Task 9)",
        "# Do not edit manually - re-generate from results if needed",
        "",
        f"receptor = '{receptor_name}'",
        f"ligand   = '{ligand_name}'",
        "",
        "cmd.reinitialize()",
        f"cmd.load('{source}/{receptor_name}.pdb', receptor)",
        f"cmd.load('{source}/{ligand_name}.pdb', ligand)",
        "",
        "cmd.hide('everything')",
        "cmd.show('cartoon', receptor)",
        "cmd.show('sticks', ligand)",
        "",
        "# Color by free energy (best = green, worst = red) - simple ramp",
        "cmd.spectrum('b', 'green_red', ligand)",
        "",
        "print('=== FlexAIDdS Thermodynamic Summary ===')",
        "",
    ]

    for i, mode in enumerate(modes[:5]):  # top 5 for clarity
        thermo = getattr(mode, "thermodynamics", None) or (mode.get("thermodynamics") if isinstance(mode, dict) else None)

        if thermo:
            g_total = thermo.get("G_total_kcal_mol", mode.get("free_energy"))
            h_eff = thermo.get("H_eff_kcal_mol")
            mts = thermo.get("minus_T_S_config_kcal_mol")
            is_exp = thermo.get("components_complete", True) is False
        else:
            g_total = getattr(mode, "free_energy", None)
            h_eff = getattr(mode, "enthalpy", None)
            mts = getattr(mode, "entropy", None)
            is_exp = False

        label = (
            f"Mode {getattr(mode, 'rank', i+1)} | "
            f"G={_fmt(g_total)} | "
            f"H={_fmt(h_eff)} | "
            f"-T*S={_fmt(mts)} kcal/mol"
        )
        if is_exp and show_experimental:
            label += " [EXPERIMENTAL]"

        lines.append(f"print('{label}')")

    # Add experimental warning at the end
    if show_experimental:
        lines.extend([
            "",
            "print('')",
            "print('WARNING: Fields marked [EXPERIMENTAL] are diagnostic only.')",
            "print('They have not been calibrated against experimental data.')",
        ])

    script = "\n".join(lines)

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(script)

    return script


def generate_markdown_report(
    result: Union[DockingResult, Dict[str, Any]],
    output_path: Union[str, Path],
    *,
    title: str = "FlexAIDdS Thermodynamic Report",
) -> str:
    """
    Generate a human-readable Markdown thermodynamic report from JSON results.
    All experimental fields are clearly labelled.
    """
    if isinstance(result, dict):
        modes = result.get("binding_modes", [])
    else:
        modes = result.binding_modes

    lines = [
        f"# {title}",
        "",
        "> **Important**: This report is generated purely from docking result JSON.",
        "> No new calculations that affect ranking were performed.",
        "",
        "| Rank | G_total (kcal/mol) | H_eff | -T*S | Experimental? | Notes |",
        "|------|---------------------|--------|------|---------------|-------|",
    ]

    for mode in modes:
        if isinstance(mode, dict):
            rank = mode.get("rank")
            thermo = mode.get("thermodynamics") or {}
            g = thermo.get("G_total_kcal_mol") or mode.get("free_energy")
            h = thermo.get("H_eff_kcal_mol") or mode.get("enthalpy")
            mts = thermo.get("minus_T_S_config_kcal_mol")
            is_exp = not thermo.get("components_complete", True)
        else:
            rank = mode.rank
            thermo = getattr(mode, "thermodynamics", None) or {}
            g = getattr(mode, "free_energy", None)
            h = getattr(mode, "enthalpy", None)
            mts = None
            is_exp = False

        exp_str = "⚠️ Yes" if is_exp else "No"
        note = " (see diagnostics)" if is_exp else ""

        lines.append(
            f"| {rank} | {_fmt(g, 2)} | {_fmt(h, 2)} | {_fmt(mts, 2)} | {exp_str} | {note} |"
        )

    lines.extend([
        "",
        "## Legend",
        "- G_total = G_config + G_vib + G_natural + G_other",
        "- All values in kcal/mol unless noted.",
        "- Fields marked 'Experimental' are diagnostic only and should not be interpreted as physical affinities.",
    ])

    report = "\n".join(lines)
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report)

    return report


def generate_temperature_scan_plot(
    scan_points: List[Dict[str, float]],
    output_path: Union[str, Path],
    title: str = "Temperature Dependence",
) -> Optional[str]:
    """
    Generate a temperature scan plot (G, H, S, Cv vs T) from scan data.

    Requires matplotlib. Returns the path to the saved figure or None if unavailable.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    if not scan_points:
        return None

    Ts = [p.get("T_K") for p in scan_points if p.get("T_K") is not None]
    if not Ts:
        return None

    Gs = [p.get("G_kcal_mol") for p in scan_points]
    Hs = [p.get("H_kcal_mol") for p in scan_points]
    Ss = [p.get("S_kcal_mol_K") for p in scan_points]
    Cvs = [p.get("Cv_kcal_mol_K") for p in scan_points]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(title)

    axes[0, 0].plot(Ts, Gs, marker="o", label="G (kcal/mol)")
    axes[0, 0].set_ylabel("G (kcal/mol)")
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(Ts, Hs, marker="s", label="H (kcal/mol)", color="orange")
    axes[0, 1].set_ylabel("H (kcal/mol)")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(Ts, Ss, marker="^", label="S (kcal/mol/K)", color="green")
    axes[1, 0].set_xlabel("Temperature (K)")
    axes[1, 0].set_ylabel("S (kcal/mol/K)")
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(Ts, Cvs, marker="d", label="Cv (kcal/mol/K)", color="red")
    axes[1, 1].set_xlabel("Temperature (K)")
    axes[1, 1].set_ylabel("Cv (kcal/mol/K)")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return str(out)


# Convenience entry point for CLI / scripts
def write_all_reports(
    result: Union[DockingResult, Dict[str, Any]],
    output_dir: Union[str, Path],
    prefix: str = "flexaidds_report",
) -> Dict[str, str]:
    """
    Write PyMOL script + Markdown report (and temperature plot if scan data is present)
    into the given directory.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    written = {}

    pymol_path = out_dir / f"{prefix}.pml"
    written["pymol"] = generate_pymol_script(result, pymol_path)

    md_path = out_dir / f"{prefix}.md"
    written["markdown"] = generate_markdown_report(result, md_path)

    # Optional: temperature scan plot if the result contains scan data
    if isinstance(result, dict) and "temperature_scan" in result:
        scan = result["temperature_scan"]
        plot_path = out_dir / f"{prefix}_temp_scan.png"
        out = generate_temperature_scan_plot(scan, plot_path)
        if out:
            written["temperature_plot"] = out

    return written
