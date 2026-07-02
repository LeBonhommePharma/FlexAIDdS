from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from .config import load_config
from .entropy import rescore_poses
from .prep import prepare_data
from .tools import run_pose_generators


app = typer.Typer(help="Astex entropy benchmark: Vina/rDock/Boltz-2 plus FlexAIDdS entropy metrics.")


ModeOpt = Annotated[str, typer.Option("--mode", help="Benchmark mode: native or non_native.")]
ConfigOpt = Annotated[Path | None, typer.Option("--config", help="Path to astex_entropy config.yaml.")]


def _abort(exc: Exception) -> None:
    typer.secho(str(exc), fg=typer.colors.RED, err=True)
    raise typer.Exit(1)


def _check_mode(mode: str, *, allow_all: bool = False) -> str:
    allowed = {"native", "non_native"} | ({"all"} if allow_all else set())
    if mode not in allowed:
        raise typer.BadParameter(f"mode must be one of: {', '.join(sorted(allowed))}")
    return mode


@app.command("data_prep")
def data_prep(
    mode: Annotated[str, typer.Option("--mode", help="native, non_native, or all.")] = "all",
    config: ConfigOpt = None,
    max_targets: Annotated[int | None, typer.Option("--max-targets", help="Limit targets/pairs for smoke runs.")] = None,
    download_missing: Annotated[bool, typer.Option("--download-missing", help="Download missing non-native PDBs from RCSB.")] = False,
    force: Annotated[bool, typer.Option("--force", help="Regenerate prepared inputs.")] = False,
) -> None:
    try:
        cfg = load_config(config)
        mode = _check_mode(mode, allow_all=True)
        counts = prepare_data(cfg, mode, max_targets=max_targets, download_missing=download_missing, force=force)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        _abort(exc)
    typer.echo(f"Prepared manifests in {Path(cfg['work_dir']) / 'manifests'}")
    for key, value in counts.items():
        typer.echo(f"{key}: {value} targets")


@app.command("run")
def run(
    mode: ModeOpt,
    config: ConfigOpt = None,
    tools: Annotated[str, typer.Option("--tools", help="Comma list: vina,rdock,boltz.")] = "vina,rdock,boltz",
    dry_run: Annotated[bool, typer.Option("--dry-run", help="Prepare commands/files but do not execute docking tools.")] = False,
    skip_missing_tools: Annotated[bool, typer.Option("--skip-missing-tools", help="Skip unavailable executables instead of failing.")] = False,
) -> None:
    try:
        cfg = load_config(config)
        mode = _check_mode(mode)
        selected = [item.strip() for item in tools.split(",") if item.strip()]
        counts = run_pose_generators(cfg, mode, selected, dry_run=dry_run, skip_missing_tools=skip_missing_tools)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        _abort(exc)
    poses_dir = Path(cfg["work_dir"]) / "poses"
    if dry_run:
        typer.echo(f"Dry-run complete; pose CSVs unchanged under {poses_dir}")
    else:
        typer.echo(f"Pose CSVs written under {poses_dir}")
    for tool, count in counts.items():
        typer.echo(f"{tool}: {count} poses")


@app.command("rescore")
def rescore(
    poses_from: Annotated[str, typer.Option("--poses_from", help="Pose source: vina, rdock, or boltz.")],
    mode: ModeOpt = "native",
    config: ConfigOpt = None,
) -> None:
    try:
        cfg = load_config(config)
        mode = _check_mode(mode)
        if poses_from not in {"vina", "rdock", "boltz"}:
            raise typer.BadParameter("poses_from must be one of: vina, rdock, boltz")
        out_csv = rescore_poses(
            cfg,
            mode=mode,
            poses_from=poses_from,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        _abort(exc)
    typer.echo(f"Rescored poses: {out_csv}")
