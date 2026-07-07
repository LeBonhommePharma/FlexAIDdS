from __future__ import annotations

import json
import math
import os
import csv
import shutil
from pathlib import Path
from typing import Any, Iterable

import yaml
from rdkit import Chem
from rdkit import RDLogger

from .io_utils import run_command, write_poses
from .models import PoseRecord, TargetRecord
from .sdf_utils import read_first_sdf_mol, smiles_from_sdf

RDLogger.DisableLog("rdApp.error")


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}
NON_LIGAND_HETATM = {
    "HOH", "WAT", "DOD", "MG", "MN", "ZN", "CA", "FE", "FE2", "CU", "NA", "K",
    "CO", "NI", "CL", "BR", "I",
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
    return read_first_sdf_mol(path, obabel, allow_unsanitized=True)


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
    return smiles_from_sdf(path, obabel)


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
                score_direction="lower",
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


def _copy_target(record: TargetRecord, **updates: str) -> TargetRecord:
    data = record.to_dict()
    data.update({key: str(value) for key, value in updates.items()})
    return TargetRecord(**data)


def _rank_from_cleft_name(path: Path) -> int:
    stem = path.stem
    for marker in ("_clf_", "_sph_"):
        if marker in stem:
            try:
                return int(stem.rsplit(marker, 1)[1])
            except ValueError:
                pass
    return 10_000


def _getcleft_anchor_from_complex(pdb_path: Path, target_id: str) -> str | None:
    try:
        from benchmarks.astex_diverse.prepare_oracle_clefts import (
            LIGAND_RESNAME,
            PEPTIDE_LIGAND,
            build_anchor_string,
            find_ligand_anchor,
            find_peptide_anchor,
        )
    except Exception:
        LIGAND_RESNAME = {}
        PEPTIDE_LIGAND = {}
        build_anchor_string = None
        find_ligand_anchor = None
        find_peptide_anchor = None

    code = target_id.split("_x_", 1)[0].split("__", 1)[0].upper()
    if code in PEPTIDE_LIGAND and find_peptide_anchor and build_anchor_string:
        resname, resnum, chain = PEPTIDE_LIGAND[code]
        found = find_peptide_anchor(pdb_path, resname, resnum, chain)
        if found:
            return build_anchor_string(resname, found[0], found[1])

    if code in LIGAND_RESNAME and find_ligand_anchor and build_anchor_string:
        resname = LIGAND_RESNAME[code]
        found = find_ligand_anchor(pdb_path, resname)
        if found:
            return build_anchor_string(resname, found[0], found[1])

    with pdb_path.open(errors="ignore") as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper() if len(line) >= 20 else ""
            if resname in NON_LIGAND_HETATM:
                continue
            chain = line[21].strip() or "-"
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            if build_anchor_string:
                return build_anchor_string(resname, resnum, chain)
            if len(resname) < 3:
                resname = "-" * (3 - len(resname)) + resname
            return f"{resname}{resnum}{chain}-"
    return None


def _getcleft_executable(cfg: dict[str, Any], *, skip_missing: bool) -> str | None:
    flex_cfg = cfg["tools"].get("flexaidds", {})
    cavity_cfg = flex_cfg.get("cavity", {}) if isinstance(flex_cfg.get("cavity", {}), dict) else {}
    configured = cavity_cfg.get("getcleft") or flex_cfg.get("getcleft") or "/Users/lp.more/Projects/Get_Cleft/Get_Cleft"
    return _require_executable(str(configured), skip_missing=skip_missing)


def _cavity_detector_kind(flex_cfg: dict[str, Any], cavity_cfg: dict[str, Any]) -> str:
    return str(cavity_cfg.get("detector", flex_cfg.get("cavity_detector", "native"))).lower()


def _native_cavity_detector_executable(cfg: dict[str, Any], *, skip_missing: bool) -> str | None:
    flex_cfg = cfg["tools"].get("flexaidds", {})
    cavity_cfg = flex_cfg.get("cavity", {}) if isinstance(flex_cfg.get("cavity", {}), dict) else {}
    configured = (
        cavity_cfg.get("cavity_detect_cli")
        or flex_cfg.get("cavity_detect_cli")
        or f"{cfg['repo_root']}/build_lto/cavity_detect_cli"
    )
    return _require_executable(str(configured), skip_missing=skip_missing)


def _cavity_radius_args(cavity_cfg: dict[str, Any]) -> list[str]:
    min_radius = float(cavity_cfg.get("min_radius_A", 1.4))
    max_radius = float(cavity_cfg.get("max_radius_A", 4.0))
    return ["--min-radius", str(min_radius), "--max-radius", str(max_radius)]


def _cavity_site_cutoff_args(cavity_cfg: dict[str, Any]) -> list[str]:
    cutoff = float(cavity_cfg.get("site_cutoff_A", 15.0))
    return ["--site-cutoff", str(cutoff)]


def _run_native_cavity_occupied(
    record: TargetRecord,
    cfg: dict[str, Any],
    mode: str,
    *,
    skip_missing: bool,
) -> Path:
    receptor = Path(record.receptor_pdb)
    if not receptor.exists():
        raise RuntimeError(f"{record.target_id}: receptor missing for native CavityDetector run: {receptor}")

    ligand = Path(record.cavity_ligand_sdf or record.reference_sdf)
    if not ligand.exists():
        raise RuntimeError(f"{record.target_id}: original ligand missing for occupied-cavity selection: {ligand}")

    exe = _native_cavity_detector_executable(cfg, skip_missing=skip_missing)
    if exe is None:
        raise RuntimeError(f"{record.target_id}: native cavity_detect_cli is required for occupied-cavity generation")

    flex_cfg = cfg["tools"].get("flexaidds", {})
    cavity_cfg = flex_cfg.get("cavity", {}) if isinstance(flex_cfg.get("cavity", {}), dict) else {}
    out_dir = Path(cfg["work_dir"]) / "cavities" / mode / record.target_dir_name / "native_occupied"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{record.target_dir_name}_occupied_sph_1.pdb"
    if out_path.exists():
        return out_path

    args = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_path),
        *_cavity_radius_args(cavity_cfg),
        *_cavity_site_cutoff_args(cavity_cfg),
    ]
    run_command(args, cwd=cfg["repo_root"], log_path=out_dir / "cavity_detect_cli.log")
    if not out_path.exists():
        raise RuntimeError(f"{record.target_id}: native CavityDetector produced no occupied cavity at {out_path}")
    return out_path


def _run_native_cavity_multi(
    record: TargetRecord,
    cfg: dict[str, Any],
    mode: str,
    *,
    top_cavities: int,
    skip_missing: bool,
) -> list[Path]:
    receptor = Path(record.receptor_pdb)
    if not receptor.exists():
        raise RuntimeError(f"{record.target_id}: receptor missing for native CavityDetector run: {receptor}")

    exe = _native_cavity_detector_executable(cfg, skip_missing=skip_missing)
    if exe is None:
        raise RuntimeError(f"{record.target_id}: native cavity_detect_cli is required for multi-cavity generation")

    flex_cfg = cfg["tools"].get("flexaidds", {})
    cavity_cfg = flex_cfg.get("cavity", {}) if isinstance(flex_cfg.get("cavity", {}), dict) else {}
    out_dir = Path(cfg["work_dir"]) / "cavities" / mode / record.target_dir_name / "native_multi"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = record.target_dir_name
    existing = sorted(out_dir.glob(f"{prefix}_sph_*.pdb"), key=_rank_from_cleft_name)
    if len(existing) >= top_cavities:
        return existing[:top_cavities]
    for stale in out_dir.glob(f"{prefix}_sph_*.pdb"):
        stale.unlink()
    args = [
        exe,
        "--receptor", str(receptor),
        "--out-dir", str(out_dir),
        "--prefix", prefix,
        "--top", str(top_cavities),
        *_cavity_radius_args(cavity_cfg),
    ]
    run_command(args, cwd=cfg["repo_root"], log_path=out_dir / "cavity_detect_cli.log")
    matches = sorted(out_dir.glob(f"{prefix}_sph_*.pdb"), key=_rank_from_cleft_name)
    if not matches:
        raise RuntimeError(f"{record.target_id}: native CavityDetector produced no ranked cavities in {out_dir}")
    return matches[:top_cavities]


def _run_getcleft_occupied(
    record: TargetRecord,
    cfg: dict[str, Any],
    mode: str,
    *,
    skip_missing: bool,
) -> Path:
    source = Path(record.cavity_source_pdb or record.source_complex)
    if not source.exists():
        raise RuntimeError(f"{record.target_id}: no cavity source PDB for occupied-cavity GetCleft run")

    exe = _getcleft_executable(cfg, skip_missing=skip_missing)
    if exe is None:
        raise RuntimeError(f"{record.target_id}: GetCleft is required for occupied-cavity generation")

    anchor = _getcleft_anchor_from_complex(source, record.target_id)
    if not anchor:
        raise RuntimeError(f"{record.target_id}: could not resolve original-ligand anchor in {source}")

    out_dir = Path(cfg["work_dir"]) / "cavities" / mode / record.target_dir_name / "occupied"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = record.target_dir_name
    existing = sorted(out_dir.glob(f"{prefix}*_clf_*.pdb"), key=_rank_from_cleft_name)
    if existing:
        return existing[0]
    run_command([exe, "-p", str(source), "-o", prefix, "-a", anchor], cwd=out_dir, log_path=out_dir / "getcleft.log")
    matches = sorted(out_dir.glob(f"{prefix}*_clf_*.pdb"), key=_rank_from_cleft_name)
    if not matches:
        raise RuntimeError(f"{record.target_id}: GetCleft produced no occupied cavity in {out_dir}")
    return matches[0]


def _run_getcleft_multi(
    record: TargetRecord,
    cfg: dict[str, Any],
    mode: str,
    *,
    top_cavities: int,
    skip_missing: bool,
) -> list[Path]:
    receptor = Path(record.receptor_pdb)
    if not receptor.exists():
        raise RuntimeError(f"{record.target_id}: receptor missing for multi-cavity GetCleft run: {receptor}")

    exe = _getcleft_executable(cfg, skip_missing=skip_missing)
    if exe is None:
        raise RuntimeError(f"{record.target_id}: GetCleft is required for multi-cavity generation")

    out_dir = Path(cfg["work_dir"]) / "cavities" / mode / record.target_dir_name / "multi"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = record.target_dir_name
    existing = sorted(out_dir.glob(f"{prefix}_sph_*.pdb"), key=_rank_from_cleft_name)
    if len(existing) >= top_cavities:
        return existing[:top_cavities]
    for stale in out_dir.glob(f"{prefix}_*"):
        stale.unlink()
    run_command(
        [exe, "-p", str(receptor), "-t", str(top_cavities), "-s", "-o", prefix],
        cwd=out_dir,
        log_path=out_dir / "getcleft.log",
    )
    matches = sorted(out_dir.glob(f"{prefix}_sph_*.pdb"), key=_rank_from_cleft_name)
    if not matches:
        raise RuntimeError(f"{record.target_id}: GetCleft produced no ranked cavities in {out_dir}")
    return matches[:top_cavities]


def _expand_flexaidds_cavity_targets(
    records: list[TargetRecord],
    cfg: dict[str, Any],
    mode: str,
    *,
    skip_missing: bool,
) -> list[TargetRecord]:
    flex_cfg = cfg["tools"].get("flexaidds", {})
    cavity_cfg = flex_cfg.get("cavity", {}) if isinstance(flex_cfg.get("cavity", {}), dict) else {}
    protocol = str(cavity_cfg.get("protocol", flex_cfg.get("cavity_protocol", "main_occupied_cavity"))).lower()
    detector = _cavity_detector_kind(flex_cfg, cavity_cfg)
    if protocol in {"", "none", "off", "legacy"}:
        return records

    if protocol in {"main", "main_occupied", "main_occupied_cavity", "occupied"}:
        expanded: list[TargetRecord] = []
        use_precomputed = bool(cavity_cfg.get("use_precomputed", False))
        for record in records:
            cleft = Path(record.pocket_pdb) if use_precomputed and record.pocket_pdb else None
            if not cleft or not cleft.exists():
                if detector in {"native", "native_cavity_detector", "cavity_detect", "cavity_detector", "flexaidds"}:
                    cleft = _run_native_cavity_occupied(record, cfg, mode, skip_missing=skip_missing)
                elif detector in {"getcleft", "get_cleft", "external_getcleft"}:
                    cleft = _run_getcleft_occupied(record, cfg, mode, skip_missing=skip_missing)
                else:
                    raise RuntimeError(f"Unknown FlexAIDdS cavity detector: {detector}")
            expanded.append(_copy_target(record, cleft_sphere_pdb=str(cleft)))
        return expanded

    if protocol in {"multi", "multi_cavity", "all_cavities"}:
        top_cavities = int(cavity_cfg.get("top_cavities", flex_cfg.get("top_cavities", 3)))
        expanded = []
        for record in records:
            if detector in {"native", "native_cavity_detector", "cavity_detect", "cavity_detector", "flexaidds"}:
                clefts = _run_native_cavity_multi(
                    record, cfg, mode, top_cavities=top_cavities, skip_missing=skip_missing
                )
            elif detector in {"getcleft", "get_cleft", "external_getcleft"}:
                clefts = _run_getcleft_multi(
                    record, cfg, mode, top_cavities=top_cavities, skip_missing=skip_missing
                )
            else:
                raise RuntimeError(f"Unknown FlexAIDdS cavity detector: {detector}")
            for rank, cleft in enumerate(clefts, 1):
                expanded.append(
                    _copy_target(
                        record,
                        target_id=f"{record.target_id}__clf{rank}",
                        cleft_sphere_pdb=str(cleft),
                        notes=f"{record.notes};source_target_id={record.target_id};cavity_rank={rank}".strip(";"),
                    )
                )
        return expanded

    raise RuntimeError(f"Unknown FlexAIDdS cavity protocol: {protocol}")


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
    boltz_cfg = cfg["tools"]["boltz"]
    wrapper = str(boltz_cfg.get("wrapper", "") or "")
    if wrapper:
        python_exe = _require_executable(str(boltz_cfg.get("python", "")), skip_missing=skip_missing)
        if python_exe is None:
            return []
        if not Path(wrapper).exists():
            if skip_missing:
                return []
            raise RuntimeError(f"Boltz wrapper not found: {wrapper}")
        base_args = [python_exe, wrapper]
    else:
        exe = _require_executable(str(boltz_cfg["executable"]), skip_missing=skip_missing)
        if exe is None:
            return []
        base_args = [exe]

    boltz_dir = target_dir / "boltz"
    out_dir = boltz_dir / "prediction"
    input_yaml = boltz_dir / "input.yaml"
    boltz_cache = Path(str(boltz_cfg.get("cache_dir", "~/.boltz"))).expanduser().resolve()
    boltz_cache.mkdir(parents=True, exist_ok=True)
    if dry_run:
        return []

    args = [*base_args, "predict", str(input_yaml), "--out_dir", str(out_dir), "--cache", str(boltz_cache)]
    if bool(boltz_cfg.get("use_msa_server", True)):
        args.append("--use_msa_server")
    if bool(boltz_cfg.get("use_potentials", True)):
        args.append("--use_potentials")
    if boltz_cfg.get("accelerator"):
        args.extend(["--accelerator", str(boltz_cfg["accelerator"])])
    if boltz_cfg.get("recycling_steps"):
        args.extend(["--recycling_steps", str(boltz_cfg["recycling_steps"])])
    if boltz_cfg.get("sampling_steps"):
        args.extend(["--sampling_steps", str(boltz_cfg["sampling_steps"])])
    if boltz_cfg.get("diffusion_samples"):
        args.extend(["--diffusion_samples", str(boltz_cfg["diffusion_samples"])])
    if boltz_cfg.get("sampling_steps_affinity"):
        args.extend(["--sampling_steps_affinity", str(boltz_cfg["sampling_steps_affinity"])])
    if boltz_cfg.get("diffusion_samples_affinity"):
        args.extend(["--diffusion_samples_affinity", str(boltz_cfg["diffusion_samples_affinity"])])
    if boltz_cfg.get("devices"):
        args.extend(["--devices", str(boltz_cfg["devices"])])
    if boltz_cfg.get("num_workers"):
        args.extend(["--num_workers", str(boltz_cfg["num_workers"])])
    if boltz_cfg.get("preprocessing_threads"):
        args.extend(["--preprocessing-threads", str(boltz_cfg["preprocessing_threads"])])
    if boltz_cfg.get("max_parallel_samples"):
        args.extend(["--max_parallel_samples", str(boltz_cfg["max_parallel_samples"])])
    if boltz_cfg.get("model"):
        args.extend(["--model", str(boltz_cfg["model"])])
    if boltz_cfg.get("output_format"):
        args.extend(["--output_format", str(boltz_cfg["output_format"])])
    if boltz_cfg.get("seed") not in (None, ""):
        args.extend(["--seed", str(boltz_cfg["seed"])])
    numba_cache = Path(cfg["work_dir"]) / "cache" / "numba"
    numba_cache.mkdir(parents=True, exist_ok=True)
    run_command(
        args,
        log_path=boltz_dir / "boltz.log",
        env={
            "NUMBA_CACHE_DIR": str(numba_cache),
            "BOLTZ_CACHE": str(boltz_cache),
            "XDG_CACHE_HOME": str(Path(cfg["work_dir"]) / "cache"),
        },
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


def _flexaidds_pair_json(records: list[TargetRecord], cfg: dict[str, Any], mode: str) -> Path:
    flex_cfg = cfg["tools"]["flexaidds"]
    include_oracle = bool(flex_cfg.get("include_oracle_site", False))
    cavity_cfg = flex_cfg.get("cavity", {}) if isinstance(flex_cfg.get("cavity", {}), dict) else {}
    pairs: list[dict[str, Any]] = []
    for idx, record in enumerate(records):
        pair = {
            "index": idx,
            "receptor_id": record.target_dir_name,
            "ligand_id": record.ligand_source_id or record.target_id,
            "receptor_pdb": record.receptor_pdb,
            "ligand_sdf": record.ligand_sdf,
            "rmsd_ref_sdf": record.reference_sdf,
        }
        if record.cleft_sphere_pdb:
            pair["cleft_sphere_file"] = record.cleft_sphere_pdb
            pair["source_target_id"] = record.target_id.split("__clf", 1)[0]
        if include_oracle and record.pocket_pdb:
            pair["oracle_site_pdb"] = record.pocket_pdb
        pairs.append(pair)

    payload = {
        "schema_version": 1,
        "name": f"astex_entropy_{mode}_flexaidds",
        "description": "Generated by benchmarks.astex_entropy for FlexAIDdS head-to-head benchmarking.",
        "oracle_mode": include_oracle,
        "cavity_protocol": cavity_cfg.get("protocol", flex_cfg.get("cavity_protocol", "")),
        "cavity_detector": _cavity_detector_kind(flex_cfg, cavity_cfg),
        "pairs": pairs,
    }
    out_path = Path(cfg["work_dir"]) / "manifests" / f"{mode}_flexaidds_crossdock.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return out_path


def _flexaidds_output_dir(cfg: dict[str, Any], mode: str) -> Path:
    return Path(cfg["work_dir"]) / "flexaidds" / mode


def _read_result_row(result_csv: Path) -> dict[str, str]:
    if not result_csv.exists():
        return {}
    with result_csv.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows[0] if rows else {}


def _parse_flexaidds_pose_score(pdb_path: Path, result_row: dict[str, str]) -> str:
    for line in pdb_path.read_text(errors="ignore").splitlines()[:80]:
        if line.startswith("REMARK CF="):
            return line.split("=", 1)[1].strip()
    return result_row.get("best_score", "")


def _extract_flexaidds_ligand_pdb(complex_pdb: Path, ligand_pdb: Path) -> bool:
    ligand_lines: list[str] = []
    ligand_serials: set[str] = set()
    conect_lines: list[str] = []
    for line in complex_pdb.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        if rec == "HETATM":
            resname = line[17:20].strip().upper() if len(line) >= 20 else ""
            if resname in NON_LIGAND_HETATM:
                continue
            ligand_lines.append(line)
            ligand_serials.add(line[6:11].strip())
        elif rec == "CONECT":
            conect_lines.append(line)
    if not ligand_lines:
        return False
    kept_conect = [
        line for line in conect_lines
        if any(line[i:i + 5].strip() in ligand_serials for i in range(6, len(line), 5))
    ]
    ligand_pdb.write_text("\n".join(ligand_lines + kept_conect + ["END"]) + "\n")
    return True


def _numeric_pose_key(path: Path) -> tuple[int, str]:
    stem = path.stem
    try:
        return int(stem.rsplit("_", 1)[1]), stem
    except (IndexError, ValueError):
        return 10_000, stem


def _collect_flexaidds_records(
    targets: list[TargetRecord],
    cfg: dict[str, Any],
    mode: str,
) -> list[PoseRecord]:
    output_root = _flexaidds_output_dir(cfg, mode)
    records: list[PoseRecord] = []
    for target in targets:
        result_dir = output_root / target.target_dir_name
        if not result_dir.exists():
            continue
        result_row = _read_result_row(result_dir / "result.csv")
        pose_pdbs = [
            path for path in result_dir.glob(f"{target.target_dir_name}_*.pdb")
            if not path.name.endswith("_INI.pdb")
        ]
        for idx, pdb_path in enumerate(sorted(pose_pdbs, key=_numeric_pose_key), start=1):
            pose_id = f"{target.target_dir_name}_flexaidds_{idx:03d}"
            pose_dir = Path(cfg["work_dir"]) / "prepared" / mode / target.target_dir_name / "poses" / "flexaidds"
            pose_dir.mkdir(parents=True, exist_ok=True)
            ligand_pdb = pose_dir / f"{pose_id}.ligand.pdb"
            pose_sdf = pose_dir / f"{pose_id}.sdf"
            if not _extract_flexaidds_ligand_pdb(pdb_path, ligand_pdb):
                continue
            run_command(
                [cfg["tools"]["obabel"], "-ipdb", str(ligand_pdb), "-osdf", "-O", str(pose_sdf)],
                log_path=pose_sdf.with_suffix(".obabel.log"),
                check=False,
            )
            if not pose_sdf.exists():
                continue
            records.append(
                PoseRecord(
                    target_id=target.target_id,
                    mode=target.mode,
                    tool="flexaidds",
                    pose_id=pose_id,
                    pose_sdf=str(pose_sdf),
                    receptor_pdb=target.receptor_pdb,
                    reference_sdf=target.reference_sdf,
                    raw_score=_parse_flexaidds_pose_score(pdb_path, result_row),
                    score_direction="lower",
                    source_file=str(pdb_path),
                )
            )
    return _fix_direction(records, cfg, "flexaidds")


def run_flexaidds_batch(
    targets: list[TargetRecord],
    cfg: dict[str, Any],
    mode: str,
    *,
    dry_run: bool,
    skip_missing: bool,
) -> list[PoseRecord]:
    exe = _require_executable(str(cfg["tools"]["flexaidds"]["benchmark_datasets"]), skip_missing=skip_missing)
    if exe is None:
        return []
    targets = _expand_flexaidds_cavity_targets(targets, cfg, mode, skip_missing=skip_missing)
    flex_cfg = cfg["tools"]["flexaidds"]
    output_dir = _flexaidds_output_dir(cfg, mode)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(str(flex_cfg.get("cache_dir", Path(cfg["work_dir"]) / "cache" / "flexaidds")))
    cache_dir.mkdir(parents=True, exist_ok=True)
    pair_json = _flexaidds_pair_json(targets, cfg, mode)
    args = [
        exe,
        "--benchmark", f"crossdock_json:{pair_json}",
        "--output", str(output_dir),
        "--threads", str(int(flex_cfg.get("threads", 1))),
        "--omp-threads", str(int(flex_cfg.get("omp_threads", 6))),
        "--cache", str(cache_dir),
        "--mode", str(flex_cfg.get("mode", "autonomous")),
        "--ga-generations", str(int(flex_cfg.get("ga_generations", 500))),
        "--ga-population", str(int(flex_cfg.get("ga_population", 1000))),
        "--job-timeout-seconds", str(int(flex_cfg.get("job_timeout_seconds", 7200))),
    ]
    if flex_cfg.get("grid_spacing"):
        args.extend(["--grid-spacing", str(flex_cfg["grid_spacing"])])
    if bool(flex_cfg.get("force", False)):
        args.append("--force")

    env = {
        "OMP_NUM_THREADS": str(int(flex_cfg.get("omp_threads", 6))),
        "FLEXAIDDS_RESTARTS": str(int(flex_cfg.get("restarts", 1))),
        "FLEXAIDDS_PARALLEL_RESTARTS": str(int(flex_cfg.get("parallel_restarts", 0))),
    }
    flexaidds_binary = cfg.get("entropy", {}).get("flexaidds_binary")
    if flexaidds_binary:
        binary = _require_executable(str(flexaidds_binary), skip_missing=skip_missing)
        if binary is not None:
            env["FLEXAIDDS_BINARY"] = binary

    command_path = output_dir / "flexaidds_command.txt"
    command_lines = [f"{key}={value}" for key, value in sorted(env.items())]
    command_lines.append(" ".join(str(arg) for arg in args))
    command_path.write_text("\n".join(command_lines) + "\n")
    if dry_run:
        return []
    run_command(args, cwd=cfg["repo_root"], log_path=output_dir / "flexaidds_benchmark.log", env=env, timeout=None)
    return _collect_flexaidds_records(targets, cfg, mode)


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
        raise ValueError("No tools selected; use --tools flexaidds,vina,rdock,boltz or a non-empty subset.")
    valid_tools = set(RUNNERS) | {"flexaidds"}
    unknown = sorted(set(tools) - valid_tools)
    if unknown:
        raise ValueError(f"Unknown tool(s): {', '.join(unknown)}")

    counts: dict[str, int] = {tool: 0 for tool in tools}
    collected: dict[str, list[PoseRecord]] = {tool: [] for tool in tools}
    poses_dir = Path(cfg["work_dir"]) / "poses"
    if not dry_run:
        for tool in tools:
            write_poses([], poses_dir / f"{mode}_{tool}_poses.csv")
        write_poses([], poses_dir / f"{mode}_all_poses.csv")

    if "flexaidds" in tools:
        records = run_flexaidds_batch(targets, cfg, mode, dry_run=dry_run, skip_missing=skip_missing_tools)
        collected["flexaidds"].extend(records)
        counts["flexaidds"] = len(records)
        if not dry_run:
            write_poses(records, poses_dir / f"{mode}_flexaidds_poses.csv")

    if any(tool != "flexaidds" for tool in tools):
        for target in targets:
            target_dir = Path(cfg["work_dir"]) / "prepared" / mode / target.target_dir_name
            prepare_tool_inputs(target, cfg, target_dir)
            for tool in tools:
                if tool == "flexaidds":
                    continue
                records = RUNNERS[tool](target, cfg, target_dir, dry_run=dry_run, skip_missing=skip_missing_tools)
                if records:
                    collected[tool].extend(records)
                    counts[tool] += len(records)
                    if not dry_run:
                        pose_csv = poses_dir / f"{mode}_{tool}_poses.csv"
                        all_csv = poses_dir / f"{mode}_all_poses.csv"
                        write_poses(collected[tool], pose_csv)
                        all_records = [record for tool_records in collected.values() for record in tool_records]
                        write_poses(all_records, all_csv)

    if not dry_run:
        all_records = [record for tool_records in collected.values() for record in tool_records]
        write_poses(all_records, poses_dir / f"{mode}_all_poses.csv")
    return counts
