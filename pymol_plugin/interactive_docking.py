"""Interactive docking workflow for FlexAID∆S (Phase 3, deliverable 3.2).

Drives the **modern** FlexAID CLI from PyMOL::

    FlexAID <receptor.pdb> <ligand.mol2> -c dock_config.json -o <prefix>

Usage:
    PyMOL> flexaids_dock receptor_obj, ligand.mol2
    PyMOL> flexaids_dock receptor_obj, ligand.mol2, site_selection=sele
"""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
import threading
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    from pymol import cmd, stored
except ImportError as exc:
    raise ImportError("PyMOL not available") from exc

try:
    from flexaidds import load_results
except ImportError as exc:
    raise ImportError(
        "flexaidds Python package is required for interactive docking"
    ) from exc


class DockingProgressCallback:
    """Tracks docking progress and updates PyMOL status."""

    def __init__(self) -> None:
        self.generation: int = 0
        self.best_cf: float = float("inf")
        self.running: bool = False
        self.cancelled: bool = False
        self.work_dir: Optional[str] = None
        self.proc = None  # subprocess.Popen handle when streaming

    def on_log_line(self, line: str) -> None:
        """Parse progress from engine log lines."""
        # Match generation lines commonly printed by the GA
        m = re.search(
            r"(?:Gen(?:eration)?|GA)\s*[:=]?\s*(\d+).*?(?:CF|best|E)\s*[:=]?\s*([-+]?\d+\.?\d*)",
            line,
            re.IGNORECASE,
        )
        if m:
            self.generation = int(m.group(1))
            self.best_cf = float(m.group(2))
            print(
                f"  Gen {self.generation}: Best CF = {self.best_cf:.4f}  "
                f"(CF/contact-function scoring proxy)"
            )
            return
        # Surface other interesting status lines briefly
        low = line.lower()
        if any(k in low for k in ("binding mode", "cluster", "free energy", "metal", "done")):
            print(f"  {line[:160]}")

    def on_generation(self, gen_num: int, best_cf: float, mean_entropy: float = 0.0) -> None:
        """Update progress after each GA generation."""
        self.generation = gen_num
        self.best_cf = best_cf
        print(
            f"  Gen {gen_num}: Best CF = {best_cf:.4f}, "
            f"Mean S = {mean_entropy:.6f}"
        )

    def cancel(self) -> None:
        """Request cancellation of the running docking job."""
        self.cancelled = True
        if self.proc is not None:
            try:
                self.proc.terminate()
            except Exception:
                pass
        print("Docking cancellation requested...")


# Module-level handle for the active background docking job
_active_docking: Optional[DockingProgressCallback] = None


def _get_selection_center(selection: str) -> Optional[tuple]:
    """Get the center of mass of a PyMOL selection.

    Returns (x, y, z) tuple or None if selection is empty/invalid.
    """
    try:
        model = cmd.get_model(selection)
        if not model.atom:
            return None
        xs = [a.coord[0] for a in model.atom]
        ys = [a.coord[1] for a in model.atom]
        zs = [a.coord[2] for a in model.atom]
        return (
            sum(xs) / len(xs),
            sum(ys) / len(ys),
            sum(zs) / len(zs),
        )
    except Exception:
        return None


def _get_selection_radius(selection: str) -> float:
    """Estimate the binding site radius from a PyMOL selection.

    Returns the maximum distance from the center of mass to any
    selected atom, plus a 2A padding.
    """
    try:
        model = cmd.get_model(selection)
        if not model.atom:
            return 10.0
        xs = [a.coord[0] for a in model.atom]
        ys = [a.coord[1] for a in model.atom]
        zs = [a.coord[2] for a in model.atom]
        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        cz = sum(zs) / len(zs)
        max_r = 0.0
        for a in model.atom:
            dx = a.coord[0] - cx
            dy = a.coord[1] - cy
            dz = a.coord[2] - cz
            r = (dx * dx + dy * dy + dz * dz) ** 0.5
            if r > max_r:
                max_r = r
        return max_r + 2.0
    except Exception:
        return 10.0


def _save_receptor_pdb(obj_name: str, output_path: str) -> bool:
    """Save a PyMOL object as a PDB file."""
    try:
        cmd.save(output_path, obj_name)
        return True
    except Exception as exc:
        print(f"ERROR: Could not save receptor: {exc}")
        return False


def _session_work_dir() -> Path:
    """User-visible session folder under FLEXAIDDS_RESULTS or temp."""
    base = os.environ.get("FLEXAIDDS_RESULTS") or os.environ.get("FLEXAIDDS_ICLOUD")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if base:
        d = Path(base) / "pymol_sessions" / f"dock_{stamp}"
        d.mkdir(parents=True, exist_ok=True)
        return d
    return Path(tempfile.mkdtemp(prefix="flexaids_dock_"))


def _write_modern_json(
    config_path: Path,
    temperature: int,
    n_results: int,
    center: tuple,
    radius: float,
) -> None:
    """Write a minimal modern JSON config for interactive docking."""
    cfg = {
        "thermodynamics": {
            "temperature": int(temperature),
            "clustering_algorithm": "CF",
            "cluster_rmsd": 2.0,
        },
        "output": {"max_results": int(n_results)},
        "optimization": {"grid_spacing": 0.375},
        "ga": {
            "num_chromosomes": 500,
            "num_generations": 800,
            "fitness_model": "SMFREE",
        },
        # Site metadata for provenance (engine may auto-detect cleft)
        "pymol_site": {
            "center": [float(center[0]), float(center[1]), float(center[2])],
            "radius_A": float(radius),
        },
    }
    config_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")


def _write_site_sphere_pdb(path: Path, center: tuple, radius: float) -> None:
    """Write a single-sphere cleft PDB (B-factor = radius) for optional inspection."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("REMARK FlexAID∆S PyMOL site sphere\n")
        fh.write(
            f"HETATM    1  C   SPH A   1    "
            f"{center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}"
            f"  1.00{radius:6.2f}           C\n"
        )
        fh.write("END\n")


def _run_docking_worker(
    receptor_pdb: str,
    ligand_file: str,
    config_json: str,
    work_dir: str,
    timeout: int,
    callback: DockingProgressCallback,
) -> None:
    """Background worker that runs the modern docking CLI."""
    global _active_docking
    cleanup_work_dir = False  # keep session folder visible by default

    try:
        from flexaidds.docking import Docking

        # Minimal .inp shell so Docking can parse temperature/paths;
        # actual invocation uses modern receptor+ligand CLI.
        shell_inp = Path(work_dir) / "dock_shell.inp"
        shell_inp.write_text(
            f"PDBNAM {receptor_pdb}\n"
            f"INPLIG {ligand_file}\n"
            f"TEMPER {300}\n"
            f"NRGOUT 10\n",
            encoding="utf-8",
        )

        docking = Docking(str(shell_inp))

        def _progress(line: str) -> None:
            if callback.cancelled:
                raise RuntimeError("Docking cancelled by user")
            callback.on_log_line(line)

        population = docking.run(
            timeout=timeout,
            receptor=receptor_pdb,
            ligand=ligand_file,
            config_json=config_json,
            output_prefix="flexaid_out",
            progress_callback=_progress,
        )

        callback.running = False

        if callback.cancelled:
            print("Docking cancelled by user.")
            return

        n_modes = len(population)
        print(f"Docking complete: {n_modes} binding mode(s) found.")
        print(f"  Session folder: {work_dir}")

        for mode_idx, mode in enumerate(population):
            thermo = mode.get_thermodynamics()
            print(
                f"  Mode {mode_idx + 1}: F={thermo.free_energy:.2f} kcal/mol, "
                f"H={thermo.mean_energy:.2f}, S={thermo.entropy:.6f}, "
                f"n_poses={mode.n_poses}"
            )

        from . import results_adapter
        results_adapter.load_docking_results(work_dir, prefix="dock")
        print("Results loaded into PyMOL with prefix 'dock'.")

    except FileNotFoundError:
        print(
            "ERROR: FlexAID binary not found. Set FLEXAIDDS_BINARY or build with:\n"
            "  cmake --build build --target FlexAID"
        )
        cleanup_work_dir = True
    except RuntimeError as exc:
        if callback.cancelled:
            print("Docking cancelled by user.")
        else:
            print(f"ERROR: Docking failed: {exc}")
    except Exception as exc:
        print(f"ERROR: Unexpected error: {exc}")
    finally:
        callback.running = False
        _active_docking = None
        if cleanup_work_dir and work_dir:
            shutil.rmtree(work_dir, ignore_errors=True)


def dock_interactive(
    receptor_obj: str,
    ligand_file: str,
    site_selection: str = "sele",
    temperature: int = 300,
    timeout: int = 300,
    n_results: int = 10,
    async_mode: bool = True,
) -> None:
    """Run FlexAID∆S docking from PyMOL using the modern binary CLI.

    Workflow:
    1. Export receptor from PyMOL object to session folder PDB
    2. Determine binding site from PyMOL selection (COM + radius)
    3. Write modern JSON config (+ site sphere for inspection)
    4. Invoke ``FLEXAIDDS_BINARY`` with ``receptor ligand -c json -o prefix``
    5. Stream log progress; auto-load results into PyMOL

    By default runs asynchronously in a background thread so PyMOL
    remains responsive.  Use ``flexaids_dock_cancel`` to abort.

    Args:
        receptor_obj: Name of the receptor PyMOL object.
        ligand_file: Path to ligand file (MOL2 or SDF).
        site_selection: PyMOL selection defining the binding site
                       (default: 'sele').
        temperature: Docking temperature in Kelvin (default: 300).
        timeout: Maximum docking time in seconds (default: 300).
        n_results: Number of output poses to generate (default: 10).
        async_mode: If True (default), run in a background thread.

    Example:
        PyMOL> flexaids_dock receptor, /path/to/ligand.mol2
        PyMOL> flexaids_dock receptor, ligand.mol2, site_selection=active_site
        PyMOL> flexaids_dock receptor, ligand.mol2, async_mode=0
    """
    global _active_docking

    receptor_obj = str(receptor_obj).strip()
    ligand_file = str(ligand_file).strip()
    site_selection = str(site_selection).strip()
    temperature = int(temperature)
    timeout = int(timeout)
    n_results = int(n_results)
    async_mode = bool(int(async_mode)) if isinstance(async_mode, str) else bool(async_mode)

    if _active_docking is not None and _active_docking.running:
        print("ERROR: A docking job is already running. "
              "Use 'flexaids_dock_cancel' to abort it first.")
        return

    # Validate receptor object exists in PyMOL
    if receptor_obj not in cmd.get_object_list():
        print(f"ERROR: Receptor object '{receptor_obj}' not found in PyMOL.")
        available = ", ".join(cmd.get_object_list())
        if available:
            print(f"  Available objects: {available}")
        return

    # Validate ligand file exists
    ligand_path = Path(ligand_file)
    if not ligand_path.is_file():
        print(f"ERROR: Ligand file not found: {ligand_file}")
        return

    # Get binding site center
    center = _get_selection_center(site_selection)
    if center is None:
        print(
            f"WARNING: Selection '{site_selection}' is empty or invalid. "
            "Using receptor center as binding site."
        )
        center = _get_selection_center(receptor_obj)
        if center is None:
            print("ERROR: Could not determine binding site center.")
            return

    radius = _get_selection_radius(site_selection)

    work_dir = _session_work_dir()
    receptor_pdb = str(work_dir / "receptor.pdb")
    ligand_copy = work_dir / ligand_path.name
    try:
        shutil.copy2(str(ligand_path.resolve()), str(ligand_copy))
    except Exception:
        ligand_copy = ligand_path  # use original path
    config_json = work_dir / "dock_config.json"
    sphere_pdb = work_dir / "site_sphere.pdb"

    print("Preparing docking (modern CLI)...")
    print(f"  Receptor: {receptor_obj}")
    print(f"  Ligand: {ligand_file}")
    print(f"  Binding site center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    print(f"  Binding site radius: {radius:.1f} A")
    print(f"  Temperature: {temperature} K")
    print(f"  Session folder: {work_dir}")

    if not _save_receptor_pdb(receptor_obj, receptor_pdb):
        if str(work_dir).startswith(tempfile.gettempdir()):
            shutil.rmtree(work_dir, ignore_errors=True)
        return

    _write_modern_json(config_json, temperature, n_results, center, radius)
    _write_site_sphere_pdb(sphere_pdb, center, radius)

    callback = DockingProgressCallback()
    callback.running = True
    callback.work_dir = str(work_dir)
    _active_docking = callback

    print(f"Starting FlexAID∆S docking (timeout: {timeout}s, "
          f"async={async_mode})...")
    print(f"  Binary: $FLEXAIDDS_BINARY or PATH (FlexAID / FlexAIDdS)")

    args = (
        receptor_pdb,
        str(ligand_copy.resolve()),
        str(config_json),
        str(work_dir),
        timeout,
        callback,
    )

    if async_mode:
        t = threading.Thread(
            target=_run_docking_worker,
            args=args,
            daemon=True,
        )
        t.start()
        print("Docking running in background. PyMOL remains usable.")
        print("  Use 'flexaids_dock_cancel' to abort.")
    else:
        _run_docking_worker(*args)


def dock_cancel() -> None:
    """Cancel a running background docking job.

    Example:
        PyMOL> flexaids_dock_cancel
    """
    global _active_docking
    if _active_docking is not None and _active_docking.running:
        _active_docking.cancel()
    else:
        print("No active docking job to cancel.")
