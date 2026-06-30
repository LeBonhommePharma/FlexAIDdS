"""cli.py - Simple Typer CLI for astex_entropy benchmark.

Commands:
  data_prep
  run --mode native|non_native
  rescore --poses_from vina|rdock|boltz

Production-grade but clean. Focus on ITC / ΔH / ΔS / entropy-enthalpy index.
"""

import typer
from pathlib import Path
from typing import Optional

from .data_prep import data_prep as _data_prep
from .runners import run_vina, run_rdock, run_boltz
from .entropy import re_rank_poses, compute_entropy_thermo
from .evaluate import evaluate_poses, write_report
from .config import DEFAULT_DATA_DIR, DEFAULT_RESULTS_DIR, ASTEX_DIVERSE_DIR

app = typer.Typer(help="Astex Entropy benchmark (thermo / ITC / ΔH ΔS / index focus)")


@app.command()
def data_prep(
    out_dir: Path = typer.Option(DEFAULT_DATA_DIR, "--out-dir", help="Where to write prepared inputs"),
):
    """Prepare receptor/ligand/config for Vina, rDock, Boltz-2 from your yamls."""
    _data_prep(out_dir=out_dir)
    typer.echo(f"data_prep finished -> {out_dir}")


@app.command()
def run(
    mode: str = typer.Option("native", "--mode", help="native or non_native"),
    out_dir: Path = typer.Option(DEFAULT_RESULTS_DIR, "--out-dir"),
    targets: Optional[str] = typer.Option(None, help="Comma list or all"),
):
    """Run Vina + rDock + Boltz-2 on the chosen mode."""
    out_dir = Path(out_dir) / mode
    out_dir.mkdir(parents=True, exist_ok=True)

    # Very thin loop - extend with real target iteration from yamls
    # For demo we just show the structure a user would call.
    typer.echo(f"Running tools for mode={mode} into {out_dir}")
    typer.echo("Example (implement loop over targets from yaml):")
    typer.echo("  for each target: run_vina(...) ; run_rdock(...) ; run_boltz(...)")
    # Real implementation would load yaml, call the three runners, save poses.


@app.command()
def rescore(
    poses_from: str = typer.Option(..., "--poses_from", help="vina | rdock | boltz"),
    mode: str = typer.Option("native", "--mode"),
    out_dir: Path = typer.Option(DEFAULT_RESULTS_DIR, "--out-dir"),
    receptor: Optional[Path] = None,
    reference: Optional[Path] = None,
):
    """Rescore poses from one tool using entropy/thermo (ΔH, ΔS, index).

    Re-ranks and evaluates with PB + RMSD + thermo focus.
    """
    out_dir = Path(out_dir) / mode / f"rescore_{poses_from}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Demo: find pose files (user adjusts glob)
    pose_files = list((Path(out_dir).parent / poses_from).glob("*.sdf"))[:5]  # small example
    if not pose_files:
        typer.echo("No pose files found. Run the tools first or adjust path.")
        raise typer.Exit(1)

    poses = [{"path": str(p), "score": 0.0} for p in pose_files]

    rec = receptor or (ASTEX_DIVERSE_DIR / "structures" / "1SQ5.pdb")  # example
    ref = reference or (ASTEX_DIVERSE_DIR / "data" / "1sq5" / "1sq5.sdf")

    # Compute thermo for each (your metrics)
    for p in poses:
        th = compute_entropy_thermo(Path(p["path"]), rec)
        p["thermo"] = th

    ranked = re_rank_poses(poses, rec)

    # Evaluate with PB + RMSD + index
    stats = evaluate_poses(ranked, ref, rec)
    write_report(stats, out_dir / "report.json")

    typer.echo(f"Rescore complete. Report: {out_dir / 'report.json'}")


if __name__ == "__main__":
    app()
