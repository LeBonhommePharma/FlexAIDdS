from __future__ import annotations

import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Annotated, Any

import typer

try:
    from .config import load_config
except ImportError:
    from benchmarks.astex_entropy.config import load_config


VALID_MODES = {"native", "non_native"}
VALID_TOOLS = {"flexaidds", "vina", "rdock", "boltz"}
DEFAULT_TOOLS = "flexaidds,vina,rdock,boltz"


def _split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def _validate_modes(value: str) -> list[str]:
    if value == "all":
        return ["native", "non_native"]
    modes = _split_csv(value)
    bad = sorted(set(modes) - VALID_MODES)
    if bad:
        raise typer.BadParameter(f"unknown mode(s): {', '.join(bad)}")
    return modes


def _validate_tools(value: str) -> list[str]:
    tools = _split_csv(value)
    bad = sorted(set(tools) - VALID_TOOLS)
    if bad:
        raise typer.BadParameter(f"unknown tool(s): {', '.join(bad)}")
    return tools


def _run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _resolve_executable(value: str) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if path.exists():
        return path
    resolved = shutil.which(value)
    return Path(resolved) if resolved else None


def _require_executable(label: str, value: str, *, hint: str) -> str:
    resolved = _resolve_executable(value)
    if resolved is None:
        raise RuntimeError(f"Required {label} missing: {value}. {hint}")
    if resolved.is_file() and not os.access(resolved, os.X_OK):
        raise RuntimeError(f"Required {label} is not executable: {resolved}. {hint}")
    return str(resolved)


def preflight_required_validators(cfg: dict[str, Any]) -> dict[str, str]:
    entropy_cfg = cfg.get("entropy", {})
    tools_cfg = cfg.get("tools", {})
    posebusters_command = str(entropy_cfg.get("posebusters_command", "") or "").strip()
    if not posebusters_command:
        raise RuntimeError("PoseBusters is required: entropy.posebusters_command is empty.")

    return {
        "posebusters": _require_executable(
            "PoseBusters validator",
            str(tools_cfg.get("posebusters", "bust")),
            hint="Install PoseBusters or fix tools.posebusters in config.yaml.",
        ),
        "tencom_eigen": _require_executable(
            "tENCoM/Eigen validator",
            str(entropy_cfg.get("tencom_binary", "")),
            hint="Build FlexAIDdS with Eigen/tENCoM support so tencom_entropy_diff exists.",
        ),
        "posebusters_command": posebusters_command,
    }


def preflight_dock_session_guard(
    cfg: dict[str, Any],
    *,
    dry_run: bool,
    workers: int = 1,
) -> dict[str, Any]:
    """Sol #9 multi-session guard: hold, mkdir lock, disk floor, binary pin.

    Import path is repo-relative so unit tests can load the same helper used by
    bash launchers. Dry-run still enforces hold + disk + workers but skips the
    exclusive lock so offline analysis cannot strand a lock.
    """
    import importlib.util
    import sys

    repo_root = Path(cfg.get("repo_root") or Path(__file__).resolve().parents[2])
    guard_path = repo_root / "scripts" / "dock_session_guard.py"
    if not guard_path.is_file():
        raise RuntimeError(f"dock_session_guard.py missing at {guard_path}")

    name = "dock_session_guard"
    spec = importlib.util.spec_from_file_location(name, guard_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {guard_path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)

    binary = Path(
        str(
            cfg.get("entropy", {}).get("flexaidds_binary")
            or cfg.get("tools", {}).get("flexaidds", {}).get("benchmark_datasets")
            or ""
        )
    )
    out_dir = Path(cfg["work_dir"])
    result = mod.preflight_dock(
        out_dir=out_dir,
        binary=binary if binary.is_file() else None,
        workers=int(workers),
        acquire_lock=not dry_run,
        copy_binary=binary.is_file() and not dry_run,
        repo_root=repo_root,
        owner="astex_entropy.orchestrate",
        note="astex_entropy orchestrator",
    )
    if not result.ok:
        raise RuntimeError("; ".join(result.messages))
    return {
        "messages": list(result.messages),
        "lock_dir": result.lock_dir,
        "free_gib": result.free_gib,
        "binary_pin": result.binary_pin or {},
    }


def _write_summary(summary: dict[str, Any], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "orchestrator_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Astex Entropy Orchestrator Summary",
        "",
        f"- Run ID: `{summary['run_id']}`",
        f"- Work dir: `{summary['work_dir']}`",
        f"- Modes: `{', '.join(summary['modes'])}`",
        f"- Tools: `{', '.join(summary['tools'])}`",
        f"- Dry run: `{summary['dry_run']}`",
        f"- PoseBusters: `{summary['preflight']['posebusters']}`",
        f"- tENCoM/Eigen: `{summary['preflight']['tencom_eigen']}`",
        "",
        "## Outputs",
        "",
    ]
    for mode, mode_summary in summary["results"].items():
        lines.append(f"### {mode}")
        lines.append("")
        lines.append(f"- Prepared targets: `{mode_summary.get('prepared_targets', 0)}`")
        for tool, count in mode_summary.get("pose_counts", {}).items():
            lines.append(f"- {tool} poses: `{count}`")
        for tool, path in mode_summary.get("rescored", {}).items():
            lines.append(f"- {tool} rescored CSV: `{path}`")
        if mode_summary.get("errors"):
            lines.append("- Errors:")
            for error in mode_summary["errors"]:
                lines.append(f"  - `{error['stage']}` / `{error.get('tool', mode)}`: {error['message']}")
        lines.append("")
    (out_dir / "orchestrator_summary.md").write_text("\n".join(lines))


def orchestrate(
    *,
    modes: list[str],
    tools: list[str],
    config_path: Path | None,
    max_targets: int | None,
    download_missing: bool,
    force: bool,
    dry_run: bool,
    skip_missing_tools: bool,
    skip_rescore: bool,
    continue_on_error: bool,
) -> Path:
    try:
        from .entropy import rescore_poses
        from .prep import prepare_data
        from .tools import run_pose_generators
    except ImportError:
        from benchmarks.astex_entropy.entropy import rescore_poses
        from benchmarks.astex_entropy.prep import prepare_data
        from benchmarks.astex_entropy.tools import run_pose_generators

    cfg = load_config(config_path)
    # Multi-session dual-dock refusal before any prep/tool work.
    flex_workers = int(cfg.get("tools", {}).get("flexaidds", {}).get("threads", 1) or 1)
    session_guard = preflight_dock_session_guard(
        cfg, dry_run=dry_run, workers=flex_workers
    )
    preflight = preflight_required_validators(cfg)
    preflight["session_guard"] = session_guard
    run_id = _run_id()
    summary: dict[str, Any] = {
        "run_id": run_id,
        "work_dir": cfg["work_dir"],
        "config_path": cfg["config_path"],
        "modes": modes,
        "tools": tools,
        "dry_run": dry_run,
        "skip_rescore": skip_rescore,
        "preflight": preflight,
        "results": {},
    }
    out_dir = Path(cfg["work_dir"]) / "orchestrator_runs" / run_id

    for mode in modes:
        mode_summary: dict[str, Any] = {
            "prepared_targets": 0,
            "pose_counts": {},
            "rescored": {},
            "errors": [],
        }
        summary["results"][mode] = mode_summary

        try:
            prep_counts = prepare_data(
                cfg,
                mode,
                max_targets=max_targets,
                download_missing=download_missing,
                force=force,
            )
            mode_summary["prepared_targets"] = int(prep_counts.get(mode, 0))
        except Exception as exc:
            mode_summary["errors"].append({"stage": "data_prep", "message": str(exc)})
            if not continue_on_error:
                _write_summary(summary, out_dir)
                raise
            continue

        try:
            pose_counts = run_pose_generators(
                cfg,
                mode,
                tools,
                dry_run=dry_run,
                skip_missing_tools=skip_missing_tools,
            )
            mode_summary["pose_counts"] = {tool: int(count) for tool, count in pose_counts.items()}
        except Exception as exc:
            mode_summary["errors"].append({"stage": "run", "message": str(exc)})
            if not continue_on_error:
                _write_summary(summary, out_dir)
                raise
            continue

        if dry_run or skip_rescore:
            continue

        for tool in tools:
            if mode_summary["pose_counts"].get(tool, 0) <= 0:
                continue
            try:
                out_csv = rescore_poses(cfg, mode=mode, poses_from=tool)
                mode_summary["rescored"][tool] = str(out_csv)
            except Exception as exc:
                mode_summary["errors"].append({"stage": "rescore", "tool": tool, "message": str(exc)})
                if not continue_on_error:
                    _write_summary(summary, out_dir)
                    raise

    _write_summary(summary, out_dir)
    return out_dir


def main(
    mode: Annotated[
        str,
        typer.Option("--mode", help="native, non_native, all, or comma list."),
    ] = "native",
    tools: Annotated[
        str,
        typer.Option("--tools", help="Comma list: flexaidds,vina,rdock,boltz."),
    ] = DEFAULT_TOOLS,
    config: Annotated[
        Path | None,
        typer.Option("--config", help="Path to astex_entropy config.yaml."),
    ] = None,
    max_targets: Annotated[
        int | None,
        typer.Option("--max-targets", help="Limit targets/pairs for smoke runs."),
    ] = None,
    download_missing: Annotated[
        bool,
        typer.Option("--download-missing", help="Download missing non-native PDB files from RCSB."),
    ] = False,
    force: Annotated[
        bool,
        typer.Option("--force", help="Regenerate prepared inputs and pass configured force behavior to tools."),
    ] = False,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Prepare commands/files but do not execute docking or rescoring."),
    ] = False,
    skip_missing_tools: Annotated[
        bool,
        typer.Option("--skip-missing-tools", help="Skip unavailable tool executables instead of failing."),
    ] = False,
    skip_rescore: Annotated[
        bool,
        typer.Option("--skip-rescore", help="Run pose generators but do not run entropy/PoseBusters rescoring."),
    ] = False,
    continue_on_error: Annotated[
        bool,
        typer.Option("--continue-on-error", help="Keep later modes/tools running after a failure."),
    ] = False,
) -> None:
    try:
        out_dir = orchestrate(
            modes=_validate_modes(mode),
            tools=_validate_tools(tools),
            config_path=config,
            max_targets=max_targets,
            download_missing=download_missing,
            force=force,
            dry_run=dry_run,
            skip_missing_tools=skip_missing_tools,
            skip_rescore=skip_rescore,
            continue_on_error=continue_on_error,
        )
    except Exception as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    typer.echo(f"Orchestrator summary: {out_dir}")


if __name__ == "__main__":
    typer.run(main)
