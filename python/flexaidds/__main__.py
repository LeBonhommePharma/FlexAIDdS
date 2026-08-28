"""Command-line entry point for the ``flexaidds`` package.

Invoked as::

    python -m flexaidds <results_dir> [--json]

Scans *results_dir* for FlexAID∆S docking output PDB files and prints a
human-readable summary of the binding modes to stdout.  Pass ``--json`` to
emit a machine-readable JSON payload instead.

Exit codes:
    0 – success
    1 – unhandled error (propagated as an exception)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

from .__version__ import __version__
from .results import load_results
from . import _tui as tui
from .__version__ import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m flexaidds",
        description="Inspect FlexAID∆S docking result directories from Python.",
    )
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument("results_dir", nargs="?", type=Path, default=None,
                        help="Directory containing docking result PDB files")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON instead of a human summary.",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable coloured output (also honours NO_COLOR / FLEXAIDDS_NO_COLOR).",
    )
    parser.add_argument(
        "--csv",
        metavar="PATH",
        default=None,
        help="Write binding-mode summary to a CSV file at PATH.",
    )
    parser.add_argument(
        "--top",
        metavar="N",
        type=int,
        default=None,
        help="Show only the top N binding modes in the summary table (default: all).",
    )
    parser.add_argument(
        "--best-only",
        "--best-mode",
        action="store_true",
        help="Print *only* the rank-1 BindingMode (ranked by the engine-emitted soft_beta_G election objective, falling back to the legacy ensemble transform) + key fields/PDB pointer.",
    )
    parser.add_argument(
        "--check-update",
        action="store_true",
        help="Check for newer versions of FlexAID∆S on GitHub.",
    )
    parser.add_argument(
        "--self-update",
        action="store_true",
        help="Check for and install the latest version via pip.",
    )
    return parser


def _fmt(value: object, width: int = 10, precision: int = 3) -> str:
    if value is None:
        return "-".center(width)
    if isinstance(value, float):
        return f"{value:{width}.{precision}f}"
    return str(value).rjust(width)


def _unit_suffix(modes) -> str:
    """Return a unit label only when provenance authorizes physical units.

    Unauthorized ensembles are proxy diagnostics in the declared input domain,
    so they must not be labelled kcal/mol.
    """
    if modes and all(
        m.scientific_provenance.allows_canonical_claims() for m in modes
    ):
        return "kcal/mol"
    return "proxy"


def _print_table(result, top_n: Optional[int]) -> None:
    modes = result.binding_modes
    if top_n is not None and top_n > 0:
        modes = modes[:top_n]

    unit = _unit_suffix(modes)
    header = (
        f"{'Mode':>5}  {'Rank':>5}  {'N_poses':>7}  "
        f"{f'F ({unit})':>14}  {f'H ({unit})':>14}  "
        f"{f'S ({unit}/K)':>16}  {'Best CF':>10}"
    )
    separator = "-" * len(header)
    print(separator)
    print(header)
    print(separator)
    for mode in modes:
        print(
            f"{mode.mode_id:>5}  {mode.rank:>5}  {mode.n_poses:>7}  "
            f"{_fmt(mode.free_energy, 14)}  {_fmt(mode.enthalpy, 14)}  "
            f"{_fmt(mode.entropy, 16, 6)}  {_fmt(mode.best_cf, 10)}"
        )
    print(separator)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.no_color or args.json or args.csv is not None:
        tui.set_enabled(False)

    # Handle update flags first
    if args.check_update or args.self_update:
        from .updater import check_for_updates, update_pip

        info = check_for_updates()
        if info is None:
            print("Could not reach GitHub API to check for updates.")
            return 1

        print(f"Current version: {info.current_version}")
        print(f"Latest version:  {info.latest_version}")

        if not info.update_available:
            print("You are up to date.")
            return 0

        print(f"Update available: {info.release_url}")
        if args.self_update:
            print("Installing update via pip...")
            version = info.latest_version.lstrip("v")
            rc = update_pip(version)
            if rc == 0:
                print("Update installed successfully.")
            else:
                print(f"pip exited with code {rc}")
            return rc

        return 0

    if args.results_dir is None:
        parser.print_help()
        return 0

    result = load_results(args.results_dir)

    if getattr(args, "best_only", False) or getattr(args, "best_mode", False):  # --best-only / --best-mode
        top = result.top_mode()
        if top is None:
            print(tui.warn("No binding modes found."))
            return 1
        # The election objective (soft_beta_G) is what the engine ranked on;
        # free_energy is a legacy ensemble transform whose units/interpretation
        # are gated by provenance, so it is not described as a free energy here.
        print(
            "Best BindingMode (engine election objective soft_beta_G, "
            f"falling back to the legacy ensemble transform; T={result.temperature} K):"
        )
        print("  " + tui.kv("mode_id", top.mode_id) + " "
              + tui.kv("rank", top.rank) + " " + tui.kv("n_poses", top.n_poses))
        print(f"  claim_validity={top.claim_validity.value}")
        print(
            "  " + tui.kv("soft_beta_G", top.soft_beta_G, tui.dG()) + " "
            + tui.kv("proxy_free_energy", top.proxy_free_energy, tui.dG()) + " "
            + tui.kv("free_energy", top.free_energy, tui.dG())
        )
        print("  " + tui.kv("enthalpy", top.enthalpy, tui.dH()) + " "
              + tui.kv("entropy", top.entropy, tui.dS()))
        print("  " + tui.kv("temperature", top.temperature, tui.T()) + " "
              + tui.kv("best_cf", top.best_cf, tui.mint()))
        # Suggest the artifact path (common layout)
        print(f"  (Look for corresponding *_mode_{top.mode_id}_*.pdb or rank 1 pose in {result.source_dir} subdirs for full REMARK thermo + coords)")
        return 0

    if args.json:
        print(result.to_json(sort_keys=True))
        return 0

    if args.csv is not None:
        result.to_csv(args.csv)
        print(f"Wrote {result.n_modes} binding mode(s) to {args.csv}")
        return 0

    print(tui.brand() + "  " + tui.muted() + "binding-mode summary" + tui.reset())
    print("  " + tui.kv("results directory", result.source_dir))
    print("  " + tui.kv("binding modes", result.n_modes, tui.mint()))
    if result.temperature is not None:
        print("  " + tui.kv("temperature", f"{result.temperature} K", tui.T()))
    print()

    if not result.binding_modes:
        return 0

    _print_table(result, args.top)

    top = result.top_mode()
    if top is not None:
        print(
            f"\nTop mode: mode_id={top.mode_id}, "
            f"free_energy={top.free_energy}, best_cf={top.best_cf}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
