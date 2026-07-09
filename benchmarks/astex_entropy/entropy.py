from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from .io_utils import (
    _local_scratch_root,
    _materialize_local,
    _safe_exists,
    _safe_mtime,
    _safe_read_csv,
    _safe_read_text,
    _safe_unlink,
    _safe_write_csv,
    _safe_write_text,
    read_poses,
    run_command,
)
from .models import PoseRecord
from .validation import rmsd_diagnostics, rmsd_to_reference, run_posebusters


K_B_KCAL = 0.001987206


def _ensure_flexaidds_python(cfg: dict[str, Any]) -> None:
    python_path = cfg.get("entropy", {}).get("python_path")
    if python_path and Path(python_path).exists() and str(python_path) not in sys.path:
        sys.path.insert(0, str(python_path))


def _compute_shannon_entropy(values: list[float], bins: int, cfg: dict[str, Any]) -> float:
    _ensure_flexaidds_python(cfg)
    try:
        from flexaidds.tencm import compute_shannon_entropy  # type: ignore

        return float(compute_shannon_entropy(values, bins))
    except Exception:
        if not values:
            return 0.0
        lo = min(values)
        hi = max(values)
        if math.isclose(lo, hi):
            return 0.0
        width = (hi - lo) / bins
        counts = [0] * bins
        for value in values:
            idx = min(bins - 1, max(0, int((value - lo) / width)))
            counts[idx] += 1
        total = sum(counts)
        return -sum((c / total) * math.log(c / total) for c in counts if c)


def _ensemble_free_energy(values: list[float], temperature: float, cfg: dict[str, Any]) -> float:
    _ensure_flexaidds_python(cfg)
    try:
        from flexaidds.thermodynamics import StatMechEngine  # type: ignore

        engine = StatMechEngine(temperature)
        engine.add_samples(values)
        return float(engine.compute().free_energy)
    except Exception:
        beta = 1.0 / (K_B_KCAL * temperature)
        e_min = min(values)
        log_z = math.log(sum(math.exp(-beta * (value - e_min)) for value in values)) - beta * e_min
        return -K_B_KCAL * temperature * log_z


def _score_as_energy(record: PoseRecord) -> float | None:
    try:
        score = float(record.raw_score)
    except (TypeError, ValueError):
        return None
    if record.score_direction.lower() in {"higher", "max", "higher_better"}:
        return -score
    return score


def _analysis_target_id(target_id: str) -> str:
    return target_id.split("__clf", 1)[0]


def _raw_score_sort_key(record: PoseRecord) -> float:
    try:
        score = float(record.raw_score)
    except (TypeError, ValueError):
        return float("inf")
    if record.score_direction.lower() in {"higher", "max", "higher_better"}:
        return -score
    return score


def limit_poses_per_target(records: list[PoseRecord], max_poses: int) -> list[PoseRecord]:
    if max_poses <= 0:
        return list(records)
    grouped: dict[str, list[PoseRecord]] = defaultdict(list)
    for record in records:
        grouped[_analysis_target_id(record.target_id)].append(record)
    limited: list[PoseRecord] = []
    for group in grouped.values():
        ranked = sorted(group, key=_raw_score_sort_key)
        limited.extend(ranked[:max_poses])
    return limited


def _boltzmann_probabilities(energies: list[float], temperature: float) -> list[float]:
    beta = 1.0 / (K_B_KCAL * temperature)
    e_min = min(energies)
    weights = [math.exp(-beta * (energy - e_min)) for energy in energies]
    total = sum(weights)
    if total <= 0:
        return [1.0 / len(energies)] * len(energies)
    return [weight / total for weight in weights]


def _pose_ligand_pdb(pose_sdf: Path, out_pdb: Path, obabel: str) -> Path:
    if not _safe_exists(out_pdb) or _safe_mtime(pose_sdf) > _safe_mtime(out_pdb):
        run_command([obabel, "-isdf", str(pose_sdf), "-opdb", "-O", str(out_pdb)], log_path=out_pdb.with_suffix(".obabel.log"))
    return out_pdb


def _make_complex_pdb(receptor_pdb: Path, ligand_pdb: Path, complex_pdb: Path) -> Path:
    if _safe_exists(complex_pdb) and _safe_mtime(complex_pdb) > _safe_mtime(ligand_pdb):
        return complex_pdb
    receptor_lines = [
        line for line in _safe_read_text(receptor_pdb).splitlines()
        if not line.startswith("END")
    ]
    ligand_lines = []
    for line in _safe_read_text(ligand_pdb).splitlines():
        if line.startswith(("ATOM", "HETATM")):
            ligand_lines.append("HETATM" + line[6:])
    _safe_write_text(complex_pdb, "\n".join(receptor_lines + ligand_lines + ["END"]) + "\n")
    return complex_pdb


def _parse_delta_f_vib(outdir: Path) -> float:
    pattern = re.compile(
        r"(?:DELTA_F_VIB|delta_F_vib_star)\s*=?\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        re.IGNORECASE,
    )
    paths: list[Path] = []
    try:
        paths = sorted(outdir.glob("*.pdb"), key=_safe_mtime, reverse=True)
    except (TimeoutError, OSError):
        paths = []
    log_path = outdir / "tencom.log"
    if _safe_exists(log_path):
        paths.append(log_path)
    for path in paths:
        for line in _safe_read_text(path).splitlines():
            match = pattern.search(line)
            if match:
                try:
                    return float(match.group(1))
                except ValueError:
                    pass
    return 0.0


def _run_tencom_for_pose(
    record: PoseRecord,
    cfg: dict[str, Any],
    *,
    require_tencom: bool,
) -> tuple[float, str]:
    tencom = Path(str(cfg["entropy"].get("tencom_binary", "")))
    if not tencom.exists():
        if require_tencom:
            raise RuntimeError(f"tENCoM binary not found: {tencom}")
        return 0.0, "tencom_missing"
    pose_sdf = _materialize_local(Path(record.pose_sdf))
    pose_dir = _local_scratch_root() / "tencom" / record.target_id / record.pose_id
    pose_dir.mkdir(parents=True, exist_ok=True)
    ligand_pdb = _pose_ligand_pdb(pose_sdf, pose_dir / "ligand_pose.pdb", cfg["tools"]["obabel"])
    receptor_pdb = _materialize_local(Path(record.receptor_pdb))
    complex_pdb = _make_complex_pdb(receptor_pdb, ligand_pdb, pose_dir / "complex.pdb")
    log_path = pose_dir / "tencom.log"
    if not _safe_exists(log_path):
        timeout = int(
            cfg.get("entropy", {}).get(
                "tencom_timeout_seconds",
                cfg.get("runtime", {}).get("command_timeout_seconds", 7200),
            )
            or 7200
        )
        result = run_command(
            [
                str(tencom),
                "--ref", str(receptor_pdb),
                str(complex_pdb),
                "--temp", str(float(cfg["entropy"].get("temperature_K", 298.15))),
                "--outdir", str(pose_dir),
                "--quiet",
            ],
            log_path=log_path,
            check=False,
            timeout=timeout,
        )
        if result.returncode != 0:
            if require_tencom:
                raise RuntimeError(f"tENCoM failed for {record.pose_id}: {result.stdout[-2000:]}")
            return 0.0, "tencom_failed"
    return _parse_delta_f_vib(pose_dir), "ok"


def _cross_rescore_label(poses_from: str, scorer_tool: str | None) -> str:
    if not scorer_tool or scorer_tool == poses_from:
        return poses_from
    return f"{poses_from}__scored_by_{scorer_tool}"


def _write_report(
    df: pd.DataFrame,
    out_dir: Path,
    poses_from: str,
    mode: str,
    *,
    scorer_tool: str | None = None,
) -> None:
    target_col = "analysis_target_id" if "analysis_target_id" in df.columns else "target_id"
    label = _cross_rescore_label(poses_from, scorer_tool)
    target_success = 0
    target_success_rate = 0.0
    if len(df):
        best = df.sort_values([target_col, "rank_entropy"]).groupby(target_col).head(1)
        target_success = int(best["success_pb"].sum())
        target_success_rate = float(best["success_pb"].mean())
    summary = {
        "mode": mode,
        "poses_from": poses_from,
        "scorer_tool": scorer_tool or poses_from,
        "label": label,
        "n_poses": int(len(df)),
        "n_targets": int(df[target_col].nunique()) if len(df) else 0,
        "success_pb_poses": int(df["success_pb"].sum()) if len(df) else 0,
        "success_pb_pose_rate": float(df["success_pb"].mean()) if len(df) else 0.0,
        "success_pb_targets": target_success,
        "success_pb_target_rate": target_success_rate,
    }
    lines = [
        f"# Astex Entropy Rescore: {mode} / {label}",
        "",
        f"- Pose source: {poses_from}",
        f"- Scoring function: {summary['scorer_tool']}",
        f"- Poses: {summary['n_poses']}",
        f"- Targets: {summary['n_targets']}",
        f"- Pose-level successes: {summary['success_pb_poses']} ({summary['success_pb_pose_rate']:.3f})",
        f"- Target-level successes (best entropy rank): {summary['success_pb_targets']} ({summary['success_pb_target_rate']:.3f})",
        "",
        "A pose counts as successful only when RMSD <= 2.0 A and PoseBusters passes.",
        "Ranking uses G_bind = H_proxy + TdS_shannon - TdS_vib, matching the FlexAIDdS thermodynamic sign convention.",
        "",
    ]
    if len(df):
        top = df.sort_values([target_col, "rank_entropy"]).groupby(target_col).head(1)
        lines.extend([
            "## Top Entropy-Ranked Pose Per Target",
            "",
            "| target_id | cavity_target_id | pose_id | G_bind | RMSD_A | success_pb |",
            "|---|---|---:|---:|---:|---:|",
        ])
        for row in top.itertuples():
            target_id = getattr(row, target_col)
            lines.append(
                f"| {target_id} | {row.target_id} | {row.pose_id} | {row.G_bind:.4f} | "
                f"{row.rmsd_A if pd.notna(row.rmsd_A) else 'NA'} | {bool(row.success_pb)} |"
            )
        lines.append("")
    (out_dir / "report.md").write_text("\n".join(lines))

    try:
        cache_dir = out_dir / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
        os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        if len(df):
            plt.figure(figsize=(9, 5))
            try:
                import seaborn as sns

                sns.scatterplot(data=df, x="G_bind", y="rmsd_A", hue="success_pb", style="target_id", legend=False)
            except ImportError:
                colors = ["#1f77b4" if bool(value) else "#d62728" for value in df["success_pb"]]
                plt.scatter(df["G_bind"], df["rmsd_A"], c=colors, edgecolors="black", linewidths=0.4)
            plt.axhline(2.0, color="black", linewidth=1, linestyle="--")
            plt.xlabel("G_bind")
            plt.ylabel("RMSD (A)")
            plt.tight_layout()
            plt.savefig(out_dir / "gbind_vs_rmsd.png", dpi=160)
            plt.close()
    except Exception as exc:
        (out_dir / "plot_error.txt").write_text(f"{type(exc).__name__}: {exc}\n")


def _rescore_records(
    records: list[PoseRecord],
    cfg: dict[str, Any],
    *,
    mode: str,
    poses_from: str,
    scorer_tool: str | None = None,
) -> pd.DataFrame:
    temperature = float(cfg["entropy"].get("temperature_K", 298.15))
    bins = int(cfg["entropy"].get("shannon_bins", 20))
    rmsd_cutoff = float(cfg["entropy"].get("rmsd_cutoff_A", 2.0))

    grouped: dict[str, list[PoseRecord]] = defaultdict(list)
    for record in records:
        grouped[_analysis_target_id(record.target_id)].append(record)

    rows: list[dict[str, Any]] = []
    for target_id, target_records in grouped.items():
        energies: list[float] = []
        for record in target_records:
            energy = _score_as_energy(record)
            if energy is None or not math.isfinite(energy):
                raise RuntimeError(
                    f"Pose {record.pose_id} from {poses_from}/{target_id} has no finite raw score; "
                    "refusing to drop it from entropy rescoring."
                )
            energies.append(energy)

        probabilities = _boltzmann_probabilities(energies, temperature)
        ensemble_shannon = _compute_shannon_entropy(energies, bins, cfg)
        ensemble_f = _ensemble_free_energy(energies, temperature, cfg)

        for record, energy, probability in zip(target_records, energies, probabilities):
            # success_pb uses whole-ligand RMSD only; fragment MCS is diagnostic.
            rmsd_info = rmsd_diagnostics(
                record.pose_sdf,
                record.reference_sdf,
                obabel=cfg["tools"].get("obabel"),
            )
            rmsd = rmsd_info.get("rmsd_A")
            pb = run_posebusters(
                record.pose_sdf,
                record.reference_sdf,
                record.receptor_pdb,
                rmsd=rmsd,
                config=cfg,
                require_posebusters=True,
            )
            tds_vib, tencom_status = _run_tencom_for_pose(record, cfg, require_tencom=True)
            surprisal = -math.log(max(probability, 1e-300))
            tds_shannon = K_B_KCAL * temperature * surprisal
            g_bind = energy + tds_shannon - tds_vib
            success_pb = bool(
                rmsd_info.get("whole_ligand")
                and rmsd is not None
                and rmsd <= rmsd_cutoff
                and pb["all_passed"]
            )
            rows.append(
                {
                    **record.to_dict(),
                    "analysis_target_id": target_id,
                    "poses_from": poses_from,
                    "scorer_tool": scorer_tool or poses_from,
                    "H_vct_proxy": energy,
                    "ensemble_free_energy": ensemble_f,
                    "ensemble_shannon_nats": ensemble_shannon,
                    "pose_probability": probability,
                    "shannon_energy_collapse": surprisal,
                    "TdS_shannon": tds_shannon,
                    "TdS_vib": tds_vib,
                    "G_bind": g_bind,
                    "rmsd_A": rmsd,
                    "fragment_mcs_rmsd_A": rmsd_info.get("fragment_mcs_rmsd_A"),
                    "whole_ligand_rmsd": bool(rmsd_info.get("whole_ligand")),
                    "pose_heavy_atoms": rmsd_info.get("pose_heavy_atoms"),
                    "ref_heavy_atoms": rmsd_info.get("ref_heavy_atoms"),
                    "posebusters_all_passed": bool(pb["all_passed"]),
                    "success_pb": success_pb,
                    "posebusters_failures": ";".join(pb.get("failures", [])),
                    "tencom_status": tencom_status,
                }
            )

    df = pd.DataFrame(rows)
    if len(df):
        rank_group = "analysis_target_id" if "analysis_target_id" in df.columns else "target_id"
        df["rank_entropy"] = df.groupby(rank_group)["G_bind"].rank(method="first", ascending=True).astype(int)
        df["rank_raw"] = df.groupby(rank_group)["H_vct_proxy"].rank(method="first", ascending=True).astype(int)
    return df


def rescore_poses(
    cfg: dict[str, Any],
    *,
    mode: str,
    poses_from: str,
    scorer_tool: str | None = None,
) -> Path:
    pose_csv = Path(cfg["work_dir"]) / "poses" / f"{mode}_{poses_from}_poses.csv"
    records = read_poses(pose_csv)
    if not records:
        raise RuntimeError(f"No poses found: {pose_csv}")

    max_poses = int(cfg.get("entropy", {}).get("max_poses_per_target", 0) or 0)
    records = limit_poses_per_target(records, max_poses)

    scorer = scorer_tool or poses_from
    if scorer != poses_from:
        from .cross_score import cross_score_records

        records = cross_score_records(
            records,
            cfg,
            mode=mode,
            poses_from=poses_from,
            scorer_tool=scorer,
        )

    label = _cross_rescore_label(poses_from, scorer if scorer != poses_from else None)
    out_dir = Path(cfg["work_dir"]) / "rescored" / mode / label
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "rescored_poses.csv"
    checkpoint_csv = out_dir / "rescored_poses.checkpoint.csv"

    grouped: dict[str, list[PoseRecord]] = defaultdict(list)
    for record in records:
        grouped[_analysis_target_id(record.target_id)].append(record)

    frames: list[pd.DataFrame] = []
    completed_targets: set[str] = set()
    prior = _safe_read_csv(checkpoint_csv)
    if prior is not None and len(prior):
        target_col = "analysis_target_id" if "analysis_target_id" in prior.columns else "target_id"
        completed_targets = set(prior[target_col].astype(str).unique())
        frames.append(prior)

    for target_id in sorted(grouped):
        if target_id in completed_targets:
            continue
        target_df = _rescore_records(
            grouped[target_id],
            cfg,
            mode=mode,
            poses_from=poses_from,
            scorer_tool=scorer if scorer != poses_from else None,
        )
        frames.append(target_df)
        partial = pd.concat(frames, ignore_index=True) if frames else target_df
        _safe_write_csv(partial, checkpoint_csv)

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    _safe_write_csv(df, out_csv)
    if _safe_exists(checkpoint_csv):
        _safe_unlink(checkpoint_csv)
    _write_report(
        df,
        out_dir,
        poses_from,
        mode,
        scorer_tool=scorer if scorer != poses_from else None,
    )
    return out_csv


def _write_cross_matrix_summary(
    cfg: dict[str, Any],
    *,
    mode: str,
    matrix: dict[str, dict[str, Any]],
) -> Path:
    out_dir = Path(cfg["work_dir"]) / "rescored" / mode
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / "cross_score_matrix.json"
    out_json.write_text(json.dumps(matrix, indent=2, sort_keys=True) + "\n")

    lines = [
        f"# Cross-Score Matrix: {mode}",
        "",
        "Rows are pose sources; columns are scoring functions. "
        "Off-diagonal cells run cross-scoring before Shannon/tENCoM rescoring.",
        "",
        "| poses_from \\\\ scorer | " + " | ".join(matrix["scorers"]) + " |",
        "|---|" + "|".join(["---:"] * len(matrix["scorers"])) + "|",
    ]
    for poses_from in matrix["poses_sources"]:
        cells = []
        for scorer in matrix["scorers"]:
            entry = matrix["cells"].get(poses_from, {}).get(scorer, {})
            if entry.get("status") == "ok":
                cells.append(f"{entry.get('success_pb', 0)}/{entry.get('n_poses', 0)}")
            else:
                cells.append(entry.get("status", "missing"))
        lines.append(f"| {poses_from} | " + " | ".join(cells) + " |")
    lines.append("")
    (out_dir / "cross_score_matrix.md").write_text("\n".join(lines))
    return out_json


def cross_rescore_matrix(
    cfg: dict[str, Any],
    *,
    mode: str,
    poses_sources: list[str],
    scorers: list[str],
    continue_on_error: bool = False,
) -> dict[str, Any]:
    matrix: dict[str, Any] = {
        "mode": mode,
        "poses_sources": poses_sources,
        "scorers": scorers,
        "cells": {},
    }
    for poses_from in poses_sources:
        matrix["cells"][poses_from] = {}
        pose_csv = Path(cfg["work_dir"]) / "poses" / f"{mode}_{poses_from}_poses.csv"
        if not pose_csv.exists() or not read_poses(pose_csv):
            for scorer in scorers:
                matrix["cells"][poses_from][scorer] = {
                    "status": "no_poses",
                    "pose_csv": str(pose_csv),
                }
            continue
        for scorer in scorers:
            try:
                out_csv = rescore_poses(
                    cfg,
                    mode=mode,
                    poses_from=poses_from,
                    scorer_tool=scorer,
                )
                df = pd.read_csv(out_csv)
                matrix["cells"][poses_from][scorer] = {
                    "status": "ok",
                    "rescored_csv": str(out_csv),
                    "n_poses": int(len(df)),
                    "success_pb": int(df["success_pb"].sum()) if len(df) and "success_pb" in df.columns else 0,
                    "success_pb_rate": float(df["success_pb"].mean()) if len(df) and "success_pb" in df.columns else 0.0,
                }
            except Exception as exc:
                matrix["cells"][poses_from][scorer] = {
                    "status": "error",
                    "message": str(exc),
                }
                if not continue_on_error:
                    _write_cross_matrix_summary(cfg, mode=mode, matrix=matrix)
                    raise
    _write_cross_matrix_summary(cfg, mode=mode, matrix=matrix)
    return matrix
