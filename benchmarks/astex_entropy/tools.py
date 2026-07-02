from __future__ import annotations

import json
import math
import os
import shutil
from pathlib import Path
from typing import Any, Iterable

import yaml
from rdkit import Chem
from rdkit import RDLogger

from .io_utils import run_command, write_poses
from .models import PoseRecord, TargetRecord

RDLogger.DisableLog("rdApp.error")


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _require_executable(name: str, *, skip_missing: bool) -> str | None:
    if not str(name).strip():
        resolved = None
    else:
        path = Path(name)
        resolved = str(path) if path.is_file() and os.access(path, os.X_OK) else shutil.which(name)
    if resolved:
        return resolved
    if skip_missing:
        return None
    raise RuntimeError(f"Required executable not found or not executable: {name}")


def _read_mol(path: str | Path, obabel: str | None = None) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    mol = supplier[0] if supplier and len(supplier) else None
    if mol is None:
        obabel = obabel or shutil.which("obabel")
        repaired = Path(path).with_suffix(".rdkit.sdf")
        if obabel:
            run_command(
                [obabel, "-isdf", str(path), "-osdf", "-O", str(repaired)],
                log_path=repaired.with_suffix(".obabel.log"),
                check=False,
            )
            supplier = Chem.SDMolSupplier(str(repaired), removeHs=False)
            mol = supplier[0] if supplier and len(supplier) else None
    if mol is None:
        raise RuntimeError(f"Could not read molecule: {path}")
    return mol


def _box_from_ligand(
    path: str | Path,
    padding: float,
    obabel: str | None = None,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    mol = _read_mol(path, obabel)
    conf = mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        xs.append(pos.x)
        ys.append(pos.y)
        zs.append(pos.z)
    if not xs:
        raise RuntimeError(f"No heavy atoms found in ligand: {path}")
    center = (
        (min(xs) + max(xs)) / 2.0,
        (min(ys) + max(ys)) / 2.0,
        (min(zs) + max(zs)) / 2.0,
    )
    size = (
        max(12.0, max(xs) - min(xs) + 2 * padding),
        max(12.0, max(ys) - min(ys) + 2 * padding),
        max(12.0, max(zs) - min(zs) + 2 * padding),
    )
    return center, size


def _ligand_smiles(path: str | Path, obabel: str | None = None) -> str:
    mol = _read_mol(path, obabel)
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def _protein_sequences_from_pdb(path: str | Path) -> dict[str, str]:
    seen: set[tuple[str, str]] = set()
    seqs: dict[str, list[str]] = {}
    with Path(path).open(errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = (line[21].strip() or "A")
            resname = line[17:20].strip().upper()
            resseq = line[22:27].strip()
            key = (chain, resseq)
            if key in seen:
                continue
            seen.add(key)
            seqs.setdefault(chain, []).append(AA3_TO_1.get(resname, "X"))
    return {chain: "".join(seq) for chain, seq in seqs.items() if seq}


def _write_vina_config(target: TargetRecord, cfg: dict[str, Any], target_dir: Path) -> None:
    vina_dir = target_dir / "vina"
    vina_dir.mkdir(parents=True, exist_ok=True)
    padding = float(cfg["tools"]["vina"].get("box_padding_A", 8.0))
    center, size = _box_from_ligand(target.reference_sdf, padding, cfg["tools"].get("obabel"))
    text = "\n".join(
        [
            f"center_x = {center[0]:.3f}",
            f"center_y = {center[1]:.3f}",
            f"center_z = {center[2]:.3f}",
            f"size_x = {size[0]:.3f}",
            f"size_y = {size[1]:.3f}",
            f"size_z = {size[2]:.3f}",
            f"exhaustiveness = {int(cfg['tools']['vina'].get('exhaustiveness', 16))}",
            f"num_modes = {int(cfg['tools']['vina'].get('num_modes', 20))}",
            f"energy_range = {float(cfg['tools']['vina'].get('energy_range', 8))}",
            "",
        ]
    )
    (vina_dir / "vina_config.txt").write_text(text)


def _write_rdock_prm(target: TargetRecord, target_dir: Path, reference_sdf: Path) -> None:
    rdock_dir = target_dir / "rdock"
    rdock_dir.mkdir(parents=True, exist_ok=True)
    prm = f"""RBT_PARAMETER_FILE_V1.00
TITLE {target.target_id}
RECEPTOR_FILE receptor.mol2
RECEPTOR_FLEX 3.0

SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {reference_sdf.name}
    RADIUS 6.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
END_SECTION

SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
"""
    (rdock_dir / "receptor.prm").write_text(prm)


def _write_boltz_input(target: TargetRecord, cfg: dict[str, Any], target_dir: Path) -> None:
    boltz_dir = target_dir / "boltz"
    boltz_dir.mkdir(parents=True, exist_ok=True)
    seqs = _protein_sequences_from_pdb(target.receptor_pdb)
    sequences: list[dict[str, Any]] = []
    for chain, sequence in sorted(seqs.items()):
        sequences.append({"protein": {"id": chain, "sequence": sequence}})
    sequences.append({"ligand": {"id": "L", "smiles": _ligand_smiles(target.ligand_sdf, cfg["tools"].get("obabel"))}})
    data = {
        "version": 1,
        "sequences": sequences,
        "properties": [{"affinity": {"binder": "L"}}],
    }
    (boltz_dir / "input.yaml").write_text(yaml.safe_dump(data, sort_keys=False))


def prepare_tool_inputs(target: TargetRecord, cfg: dict[str, Any], target_dir: Path, *, force: bool = False) -> None:
    obabel = cfg["tools"]["obabel"]
    vina_dir = target_dir / "vina"
    rdock_dir = target_dir / "rdock"
    boltz_dir = target_dir / "boltz"
    vina_dir.mkdir(parents=True, exist_ok=True)
    rdock_dir.mkdir(parents=True, exist_ok=True)
    boltz_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = vina_dir / "receptor.pdbqt"
    ligand_pdbqt = vina_dir / "ligand.pdbqt"
    if force or not receptor_pdbqt.exists():
        run_command([obabel, "-ipdb", target.receptor_pdb, "-opdbqt", "-O", str(receptor_pdbqt), "-xr", "-xc"], log_path=vina_dir / "obabel_receptor.log")
    if force or not ligand_pdbqt.exists():
        run_command([obabel, "-isdf", target.ligand_sdf, "-opdbqt", "-O", str(ligand_pdbqt), "-h"], log_path=vina_dir / "obabel_ligand.log")
    _write_vina_config(target, cfg, target_dir)

    receptor_mol2 = rdock_dir / "receptor.mol2"
    if force or not receptor_mol2.exists():
        run_command([obabel, "-ipdb", target.receptor_pdb, "-omol2", "-O", str(receptor_mol2)], log_path=rdock_dir / "obabel_receptor.log")
    rdock_ligand_sdf = rdock_dir / "ligand.sdf"
    rdock_reference_sdf = rdock_dir / "reference.sdf"
    if force or not rdock_ligand_sdf.exists():
        run_command([obabel, "-isdf", target.ligand_sdf, "-osdf", "-O", str(rdock_ligand_sdf)], log_path=rdock_dir / "obabel_ligand.log")
    if force or not rdock_reference_sdf.exists():
        run_command([obabel, "-isdf", target.reference_sdf, "-osdf", "-O", str(rdock_reference_sdf)], log_path=rdock_dir / "obabel_reference.log")
    _write_rdock_prm(target, target_dir, rdock_reference_sdf)
    _write_boltz_input(target, cfg, target_dir)


def _split_sdf(path: Path, out_dir: Path, prefix: str, scores: list[float | None], tool: str, target: TargetRecord, source: Path) -> list[PoseRecord]:
    out_dir.mkdir(parents=True, exist_ok=True)
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    records: list[PoseRecord] = []
    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        pose_id = f"{prefix}_{idx + 1:03d}"
        out_path = out_dir / f"{pose_id}.sdf"
        writer = Chem.SDWriter(str(out_path))
        if idx < len(scores) and scores[idx] is not None:
            mol.SetProp("raw_score", f"{scores[idx]:.6g}")
        mol.SetProp("_Name", pose_id)
        writer.write(mol)
        writer.close()
        raw_score = "" if idx >= len(scores) or scores[idx] is None else f"{scores[idx]:.6g}"
        records.append(
            PoseRecord(
                target_id=target.target_id,
                mode=target.mode,
                tool=tool,
                pose_id=pose_id,
                pose_sdf=str(out_path),
                receptor_pdb=target.receptor_pdb,
                reference_sdf=target.reference_sdf,
                raw_score=raw_score,
                score_direction=str(tool),
                source_file=str(source),
            )
        )
    return records


def _parse_vina_scores(log_path: Path) -> list[float | None]:
    scores: list[float | None] = []
    if not log_path.exists():
        return scores
    for line in log_path.read_text(errors="ignore").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                scores.append(float(parts[1]))
            except ValueError:
                pass
    return scores


def _parse_rdock_scores(sdf_path: Path) -> list[float | None]:
    scores: list[float | None] = []
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    for mol in supplier:
        score: float | None = None
        if mol is not None:
            for prop in ("SCORE", "SCORE.INTER", "Rbt.Score"):
                if mol.HasProp(prop):
                    try:
                        score = float(mol.GetProp(prop))
                        break
                    except ValueError:
                        pass
        scores.append(score)
    return scores


def _tool_score_direction(cfg: dict[str, Any], tool: str) -> str:
    return str(cfg["tools"].get(tool, {}).get("score_direction", "lower"))


def _command_timeout(cfg: dict[str, Any]) -> int | None:
    value = cfg.get("runtime", {}).get("command_timeout_seconds", 7200)
    if value in (None, "", 0, "0"):
        return None
    return int(value)


def _fix_direction(records: Iterable[PoseRecord], cfg: dict[str, Any], tool: str) -> list[PoseRecord]:
    direction = _tool_score_direction(cfg, tool)
    return [
        PoseRecord(**{**record.to_dict(), "score_direction": direction})
        for record in records
    ]


def run_vina(target: TargetRecord, cfg: dict[str, Any], target_dir: Path, *, dry_run: bool, skip_missing: bool) -> list[PoseRecord]:
    exe = _require_executable(str(cfg["tools"]["vina"]["executable"]), skip_missing=skip_missing)
    if exe is None:
        return []
    vina_dir = target_dir / "vina"
    out_pdbqt = vina_dir / "vina_out.pdbqt"
    log_path = vina_dir / "vina.stdout.log"
    if dry_run:
        return []
    args = [
            exe,
            "--receptor", str(vina_dir / "receptor.pdbqt"),
            "--ligand", str(vina_dir / "ligand.pdbqt"),
            "--config", str(vina_dir / "vina_config.txt"),
            "--out", str(out_pdbqt),
    ]
    if cfg["tools"]["vina"].get("cpu"):
        args.extend(["--cpu", str(cfg["tools"]["vina"]["cpu"])])
    run_command(
        args,
        log_path=vina_dir / "vina.stdout.log",
        timeout=_command_timeout(cfg),
    )
    poses_sdf = vina_dir / "vina_poses.sdf"
    run_command(
        [cfg["tools"]["obabel"], "-ipdbqt", str(out_pdbqt), "-osdf", "-O", str(poses_sdf)],
        log_path=vina_dir / "obabel_poses.log",
        timeout=_command_timeout(cfg),
    )
    if not poses_sdf.exists():
        return []
    scores = _parse_vina_scores(log_path)
    records = _split_sdf(poses_sdf, target_dir / "poses" / "vina", f"{target.target_id}_vina", scores, "vina", target, poses_sdf)
    return _fix_direction(records, cfg, "vina")


def _rdock_env(cfg: dict[str, Any]) -> dict[str, str]:
    rbdock_path = Path(str(cfg["tools"]["rdock"]["rbdock"]))
    if rbdock_path.parent.name != "bin":
        return {}
    prefix = rbdock_path.parent.parent
    env: dict[str, str] = {}
    for candidate in (prefix / "share" / "rdock", prefix / "rDock", prefix):
        if (candidate / "data" / "RbtElements.dat").exists():
            env["RBT_ROOT"] = str(candidate)
            break
    lib = prefix / "lib"
    if lib.exists():
        existing = os.environ.get("DYLD_LIBRARY_PATH", "")
        env["DYLD_LIBRARY_PATH"] = str(lib) if not existing else f"{lib}:{existing}"
    return env


def run_rdock(target: TargetRecord, cfg: dict[str, Any], target_dir: Path, *, dry_run: bool, skip_missing: bool) -> list[PoseRecord]:
    rbcavity = _require_executable(str(cfg["tools"]["rdock"]["rbcavity"]), skip_missing=skip_missing)
    rbdock = _require_executable(str(cfg["tools"]["rdock"]["rbdock"]), skip_missing=skip_missing)
    if rbcavity is None or rbdock is None:
        return []
    rdock_dir = target_dir / "rdock"
    prm = rdock_dir / "receptor.prm"
    out_prefix = "rdock_out"
    out_sdf = rdock_dir / "rdock_out.sd"
    if dry_run:
        return []
    env = _rdock_env(cfg)
    run_command(
        [rbcavity, "-r", prm.name, "-was"],
        cwd=rdock_dir,
        log_path=rdock_dir / "rbcavity.log",
        env=env,
        timeout=_command_timeout(cfg),
    )
    ligand_sdf = rdock_dir / "ligand.sdf"
    if not ligand_sdf.exists():
        ligand_sdf = Path(target.ligand_sdf)
    ligand_arg = ligand_sdf.name if ligand_sdf.parent == rdock_dir else str(ligand_sdf)
    run_command(
        [
            rbdock,
            "-i", ligand_arg,
            "-o", out_prefix,
            "-r", prm.name,
            "-p", "dock.prm",
            "-n", str(int(cfg["tools"]["rdock"].get("nrun", 20))),
        ],
        cwd=rdock_dir,
        log_path=rdock_dir / "rbdock.log",
        env=env,
        timeout=_command_timeout(cfg),
    )
    if not out_sdf.exists():
        return []
    records = _split_sdf(out_sdf, target_dir / "poses" / "rdock", f"{target.target_id}_rdock", _parse_rdock_scores(out_sdf), "rdock", target, out_sdf)
    return _fix_direction(records, cfg, "rdock")


def _collect_boltz_scores(out_dir: Path) -> list[float | None]:
    def walk(obj: Any) -> float | None:
        if isinstance(obj, dict):
            if "affinity_pred_value" in obj:
                try:
                    return float(obj["affinity_pred_value"])
                except (TypeError, ValueError):
                    pass
            for key in ("confidence_score", "ptm", "iptm"):
                if key in obj:
                    try:
                        return -float(obj[key])
                    except (TypeError, ValueError):
                        pass
            for value in obj.values():
                found = walk(value)
                if found is not None:
                    return found
        elif isinstance(obj, list):
            for value in obj:
                found = walk(value)
                if found is not None:
                    return found
        return None

    scores: list[float | None] = []
    for path in sorted(out_dir.rglob("*confidence*.json")):
        try:
            score = walk(json.loads(path.read_text(errors="ignore")))
        except json.JSONDecodeError:
            continue
        if score is not None:
            scores.append(score)
    return scores


def _extract_boltz_ligand_pdb(structure_pdb: Path, ligand_pdb: Path) -> Path | None:
    lines = structure_pdb.read_text(errors="ignore").splitlines()
    ligand_lines = [
        line for line in lines
        if line.startswith("HETATM") and (len(line) <= 21 or line[21].strip() in {"L", ""})
    ]
    if not ligand_lines:
        ligand_lines = [line for line in lines if line.startswith("HETATM")]
    if not ligand_lines:
        return None
    ligand_pdb.write_text("\n".join(ligand_lines + ["END"]) + "\n")
    return ligand_pdb


def _boltz_structure_to_ligand_sdf(structure_path: Path, out_sdf: Path, obabel: str) -> Path | None:
    work_pdb = out_sdf.with_suffix(".complex.pdb")
    if structure_path.suffix.lower() == ".pdb":
        work_pdb = structure_path
    else:
        fmt = structure_path.suffix.lower().lstrip(".")
        run_command(
            [obabel, f"-i{fmt}", str(structure_path), "-opdb", "-O", str(work_pdb)],
            log_path=out_sdf.with_suffix(".complex_obabel.log"),
            check=False,
        )
        if not work_pdb.exists():
            return None
    ligand_pdb = out_sdf.with_suffix(".ligand.pdb")
    if _extract_boltz_ligand_pdb(work_pdb, ligand_pdb) is None:
        return None
    run_command([obabel, "-ipdb", str(ligand_pdb), "-osdf", "-O", str(out_sdf)], log_path=out_sdf.with_suffix(".obabel.log"))
    return out_sdf


def run_boltz(target: TargetRecord, cfg: dict[str, Any], target_dir: Path, *, dry_run: bool, skip_missing: bool) -> list[PoseRecord]:
    exe = _require_executable(str(cfg["tools"]["boltz"]["executable"]), skip_missing=skip_missing)
    if exe is None:
        return []
    boltz_dir = target_dir / "boltz"
    out_dir = boltz_dir / "prediction"
    input_yaml = boltz_dir / "input.yaml"
    boltz_cache = Path(str(cfg["tools"]["boltz"].get("cache_dir", Path(cfg["work_dir"]) / "cache" / "boltz")))
    boltz_cache.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return []
    args = [exe, "predict", str(input_yaml), "--out_dir", str(out_dir), "--cache", str(boltz_cache)]
    if bool(cfg["tools"]["boltz"].get("use_msa_server", True)):
        args.append("--use_msa_server")
    if bool(cfg["tools"]["boltz"].get("use_potentials", True)):
        args.append("--use_potentials")
    if cfg["tools"]["boltz"].get("accelerator"):
        args.extend(["--accelerator", str(cfg["tools"]["boltz"]["accelerator"])])
    if cfg["tools"]["boltz"].get("recycling_steps"):
        args.extend(["--recycling_steps", str(cfg["tools"]["boltz"]["recycling_steps"])])
    if cfg["tools"]["boltz"].get("sampling_steps"):
        args.extend(["--sampling_steps", str(cfg["tools"]["boltz"]["sampling_steps"])])
    if cfg["tools"]["boltz"].get("diffusion_samples"):
        args.extend(["--diffusion_samples", str(cfg["tools"]["boltz"]["diffusion_samples"])])
    if cfg["tools"]["boltz"].get("devices"):
        args.extend(["--devices", str(cfg["tools"]["boltz"]["devices"])])
    if cfg["tools"]["boltz"].get("num_workers"):
        args.extend(["--num_workers", str(cfg["tools"]["boltz"]["num_workers"])])
    if cfg["tools"]["boltz"].get("preprocessing_threads"):
        args.extend(["--preprocessing-threads", str(cfg["tools"]["boltz"]["preprocessing_threads"])])
    if cfg["tools"]["boltz"].get("max_parallel_samples"):
        args.extend(["--max_parallel_samples", str(cfg["tools"]["boltz"]["max_parallel_samples"])])
    if cfg["tools"]["boltz"].get("model"):
        args.extend(["--model", str(cfg["tools"]["boltz"]["model"])])
    if cfg["tools"]["boltz"].get("output_format"):
        args.extend(["--output_format", str(cfg["tools"]["boltz"]["output_format"])])
    numba_cache = Path(cfg["work_dir"]) / "cache" / "numba"
    numba_cache.mkdir(parents=True, exist_ok=True)
    run_command(
        args,
        log_path=boltz_dir / "boltz.log",
        env={"NUMBA_CACHE_DIR": str(numba_cache), "BOLTZ_CACHE": str(boltz_cache)},
        timeout=_command_timeout(cfg),
    )

    candidates = sorted(out_dir.rglob("*.sdf")) + sorted(out_dir.rglob("*.mol2"))
    converted: list[Path] = []
    for idx, candidate in enumerate(candidates[:20]):
        out_sdf = boltz_dir / f"boltz_pose_source_{idx + 1:03d}.sdf"
        if candidate.suffix.lower() == ".sdf":
            shutil.copy2(candidate, out_sdf)
        else:
            run_command([cfg["tools"]["obabel"], str(candidate), "-osdf", "-O", str(out_sdf)], log_path=out_sdf.with_suffix(".obabel.log"))
        converted.append(out_sdf)
    if not converted:
        structures = sorted(out_dir.rglob("*.pdb")) + sorted(out_dir.rglob("*.cif")) + sorted(out_dir.rglob("*.mmcif"))
        for idx, structure in enumerate(structures[:20]):
            out_sdf = boltz_dir / f"boltz_pose_source_{idx + 1:03d}.sdf"
            if _boltz_structure_to_ligand_sdf(structure, out_sdf, cfg["tools"]["obabel"]) is None:
                continue
            converted.append(out_sdf)
    if not converted:
        return []

    records: list[PoseRecord] = []
    scores = _collect_boltz_scores(out_dir)
    for source_idx, sdf in enumerate(converted):
        source_scores = [scores[source_idx]] if source_idx < len(scores) else []
        source_records = _split_sdf(
            sdf,
            target_dir / "poses" / "boltz",
            f"{target.target_id}_boltz_{source_idx + 1:03d}",
            source_scores,
            "boltz",
            target,
            sdf,
        )
        records.extend(source_records)
    return _fix_direction(records, cfg, "boltz")


RUNNERS = {
    "vina": run_vina,
    "rdock": run_rdock,
    "boltz": run_boltz,
}


def run_pose_generators(
    cfg: dict[str, Any],
    mode: str,
    tools: list[str],
    *,
    dry_run: bool = False,
    skip_missing_tools: bool = False,
) -> dict[str, int]:
    from .io_utils import read_targets

    manifest = Path(cfg["work_dir"]) / "manifests" / f"{mode}_manifest.csv"
    targets = read_targets(manifest)
    if not targets:
        raise RuntimeError(f"No targets found in manifest: {manifest}. Run data_prep first.")
    if not tools:
        raise ValueError("No tools selected; use --tools vina,rdock,boltz or a non-empty subset.")
    unknown = sorted(set(tools) - set(RUNNERS))
    if unknown:
        raise ValueError(f"Unknown tool(s): {', '.join(unknown)}")
    counts: dict[str, int] = {tool: 0 for tool in tools}
    collected: dict[str, list[PoseRecord]] = {tool: [] for tool in tools}
    poses_dir = Path(cfg["work_dir"]) / "poses"
    if not dry_run:
        for tool in tools:
            write_poses([], poses_dir / f"{mode}_{tool}_poses.csv")
        write_poses([], poses_dir / f"{mode}_all_poses.csv")
    for target in targets:
        target_dir = Path(cfg["work_dir"]) / "prepared" / mode / target.target_dir_name
        prepare_tool_inputs(target, cfg, target_dir)
        for tool in tools:
            records = RUNNERS[tool](target, cfg, target_dir, dry_run=dry_run, skip_missing=skip_missing_tools)
            if records:
                collected[tool].extend(records)
                pose_csv = poses_dir / f"{mode}_{tool}_poses.csv"
                all_csv = poses_dir / f"{mode}_all_poses.csv"
                write_poses(collected[tool], pose_csv)
                all_records = [record for tool_records in collected.values() for record in tool_records]
                write_poses(all_records, all_csv)
                counts[tool] += len(records)
    return counts
