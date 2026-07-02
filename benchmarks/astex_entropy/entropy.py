from __future__ import annotations

import csv
import math
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from .io_utils import read_poses, run_command
from .models import PoseRecord
from .validation import rmsd_to_reference, run_posebusters


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


def _boltzmann_probabilities(energies: list[float], temperature: float) -> list[float]:
    beta = 1.0 / (K_B_KCAL * temperature)
    e_min = min(energies)
    weights = [math.exp(-beta * (energy - e_min)) for energy in energies]
    total = sum(weights)
    if total <= 0:
        return [1.0 / len(energies)] * len(energies)
    return [weight / total for weight in weights]


def _pose_ligand_pdb(pose_sdf: Path, out_pdb: Path, obabel: str) -> Path:
    if not out_pdb.exists() or pose_sdf.stat().st_mtime > out_pdb.stat().st_mtime:
        run_command([obabel, "-isdf", str(pose_sdf), "-opdb", "-O", str(out_pdb)], log_path=out_pdb.with_suffix(".obabel.log"))
    return out_pdb


def _make_complex_pdb(receptor_pdb: Path, ligand_pdb: Path, complex_pdb: Path) -> Path:
    if complex_pdb.exists() and complex_pdb.stat().st_mtime > ligand_pdb.stat().st_mtime:
        return complex_pdb
    receptor_lines = [
        line for line in receptor_pdb.read_text(errors="ignore").splitlines()
        if not line.startswith("END")
    ]
    ligand_lines = []
    for line in ligand_pdb.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            ligand_lines.append("HETATM" + line[6:])
    complex_pdb.write_text("\n".join(receptor_lines + ligand_lines + ["END"]) + "\n")
    return complex_pdb


def _parse_delta_f_vib(outdir: Path) -> float:
    pdbs = sorted(outdir.glob("*.pdb"), key=lambda p: p.stat().st_mtime, reverse=True)
    for pdb in pdbs:
        for line in pdb.read_text(errors="ignore").splitlines():
            if line.startswith("REMARK") and "DELTA_F_VIB=" in line:
                try:
                    return float(line.split("DELTA_F_VIB=", 1)[1].split()[0])
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
    pose_sdf = Path(record.pose_sdf)
    pose_dir = pose_sdf.parent / "tencom" / record.pose_id
    pose_dir.mkdir(parents=True, exist_ok=True)
    ligand_pdb = _pose_ligand_pdb(pose_sdf, pose_dir / "ligand_pose.pdb", cfg["tools"]["obabel"])
    complex_pdb = _make_complex_pdb(Path(record.receptor_pdb), ligand_pdb, pose_dir / "complex.pdb")
    log_path = pose_dir / "tencom.log"
    if not log_path.exists():
        result = run_command(
            [
                str(tencom),
                "--ref", record.receptor_pdb,
                str(complex_pdb),
                "--temp", str(float(cfg["entropy"].get("temperature_K", 298.15))),
                "--outdir", str(pose_dir),
                "--quiet",
            ],
            log_path=log_path,
            check=False,
        )
        if result.returncode != 0:
            if require_tencom:
                raise RuntimeError(f"tENCoM failed for {record.pose_id}: {result.stdout[-2000:]}")
            return 0.0, "tencom_failed"
    return _parse_delta_f_vib(pose_dir), "ok"


def _write_report(df: pd.DataFrame, out_dir: Path, poses_from: str, mode: str) -> None:
    summary = {
        "mode": mode,
        "poses_from": poses_from,
        "n_poses": int(len(df)),
        "n_targets": int(df["target_id"].nunique()) if len(df) else 0,
        "success_pb": int(df["success_pb"].sum()) if len(df) else 0,
        "success_pb_rate": float(df["success_pb"].mean()) if len(df) else 0.0,
    }
    lines = [
        f"# Astex Entropy Rescore: {mode} / {poses_from}",
        "",
        f"- Poses: {summary['n_poses']}",
        f"- Targets: {summary['n_targets']}",
        f"- PoseBusters successes: {summary['success_pb']} ({summary['success_pb_rate']:.3f})",
        "",
        "A pose counts as successful only when RMSD <= 2.0 A and PoseBusters passes.",
        "Ranking uses G_bind = H_proxy + TdS_shannon - TdS_vib, matching the FlexAIDdS thermodynamic sign convention.",
        "",
    ]
    if len(df):
        top = df.sort_values(["target_id", "rank_entropy"]).groupby("target_id").head(1)
        lines.extend([
            "## Top Entropy-Ranked Pose Per Target",
            "",
            "| target_id | pose_id | G_bind | RMSD_A | success_pb |",
            "|---|---:|---:|---:|---:|",
        ])
        for row in top.itertuples():
            lines.append(
                f"| {row.target_id} | {row.pose_id} | {row.G_bind:.4f} | "
                f"{row.rmsd_A if pd.notna(row.rmsd_A) else 'NA'} | {bool(row.success_pb)} |"
            )
        lines.append("")
    (out_dir / "report.md").write_text("\n".join(lines))

    try:
        cache_dir = out_dir / ".cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
        os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
        import matplotlib.pyplot as plt
        import seaborn as sns

        if len(df):
            plt.figure(figsize=(9, 5))
            sns.scatterplot(data=df, x="G_bind", y="rmsd_A", hue="success_pb", style="target_id", legend=False)
            plt.axhline(2.0, color="black", linewidth=1, linestyle="--")
            plt.tight_layout()
            plt.savefig(out_dir / "gbind_vs_rmsd.png", dpi=160)
            plt.close()
    except Exception:
        pass


def rescore_poses(
    cfg: dict[str, Any],
    *,
    mode: str,
    poses_from: str,
) -> Path:
    pose_csv = Path(cfg["work_dir"]) / "poses" / f"{mode}_{poses_from}_poses.csv"
    records = read_poses(pose_csv)
    if not records:
        raise RuntimeError(f"No poses found: {pose_csv}")

    temperature = float(cfg["entropy"].get("temperature_K", 298.15))
    bins = int(cfg["entropy"].get("shannon_bins", 20))
    rmsd_cutoff = float(cfg["entropy"].get("rmsd_cutoff_A", 2.0))

    grouped: dict[str, list[PoseRecord]] = defaultdict(list)
    for record in records:
        grouped[record.target_id].append(record)

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
            rmsd = rmsd_to_reference(
                record.pose_sdf,
                record.reference_sdf,
                obabel=cfg["tools"].get("obabel"),
            )
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
            success_pb = bool(rmsd is not None and rmsd <= rmsd_cutoff and pb["all_passed"])
            rows.append(
                {
                    **record.to_dict(),
                    "H_vct_proxy": energy,
                    "ensemble_free_energy": ensemble_f,
                    "ensemble_shannon_nats": ensemble_shannon,
                    "pose_probability": probability,
                    "shannon_energy_collapse": surprisal,
                    "TdS_shannon": tds_shannon,
                    "TdS_vib": tds_vib,
                    "G_bind": g_bind,
                    "rmsd_A": rmsd,
                    "posebusters_all_passed": bool(pb["all_passed"]),
                    "success_pb": success_pb,
                    "posebusters_failures": ";".join(pb.get("failures", [])),
                    "tencom_status": tencom_status,
                }
            )

    df = pd.DataFrame(rows)
    if len(df):
        df["rank_entropy"] = df.groupby("target_id")["G_bind"].rank(method="first", ascending=True).astype(int)
        df["rank_raw"] = df.groupby("target_id")["H_vct_proxy"].rank(method="first", ascending=True).astype(int)
    out_dir = Path(cfg["work_dir"]) / "rescored" / mode / poses_from
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "rescored_poses.csv"
    df.to_csv(out_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    _write_report(df, out_dir, poses_from, mode)
    return out_csv
