from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import shlex
from typing import Any

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdMolAlign

from .io_utils import run_command
from .sdf_utils import normalise_v2000_counts_file, read_first_sdf_mol

RDLogger.DisableLog("rdApp.error")


def _read_first_sdf(path: str | Path, obabel: str | None = None) -> Chem.Mol:
    return read_first_sdf_mol(path, obabel, allow_unsanitized=True)


def _remove_hs_for_rmsd(mol: Chem.Mol) -> Chem.Mol:
    try:
        return Chem.RemoveHs(mol, sanitize=False)
    except TypeError:
        return Chem.RemoveHs(mol)
    except Exception:
        return mol


def _ordered_heavy_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float | None:
    if pose.GetNumConformers() == 0 or ref.GetNumConformers() == 0:
        return None
    pose_conf = pose.GetConformer()
    ref_conf = ref.GetConformer()
    pose_atoms = [atom for atom in pose.GetAtoms() if atom.GetAtomicNum() > 1]
    ref_atoms = [atom for atom in ref.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(pose_atoms) != len(ref_atoms):
        return None
    total = 0.0
    for pose_atom, ref_atom in zip(pose_atoms, ref_atoms, strict=True):
        if pose_atom.GetAtomicNum() != ref_atom.GetAtomicNum():
            return None
        p = pose_conf.GetAtomPosition(pose_atom.GetIdx())
        r = ref_conf.GetAtomPosition(ref_atom.GetIdx())
        total += (p.x - r.x) ** 2 + (p.y - r.y) ** 2 + (p.z - r.z) ** 2
    return (total / len(pose_atoms)) ** 0.5 if pose_atoms else None


def _normalise_sdf_for_posebusters(path: str | Path, suffix: str, obabel: str | None = None) -> Path:
    source = Path(path)
    out_path = source.with_suffix(suffix)
    executable = obabel or shutil.which("obabel")
    if not executable:
        return source
    run_command(
        [executable, "-isdf", str(source), "-osdf", "-O", str(out_path)],
        log_path=out_path.with_suffix(".obabel.log"),
        check=False,
    )
    if out_path.exists():
        normalise_v2000_counts_file(out_path)
    return out_path if out_path.exists() else source


def rmsd_to_reference(
    pose_sdf: str | Path,
    reference_sdf: str | Path,
    *,
    obabel: str | None = None,
) -> float | None:
    try:
        pose = _remove_hs_for_rmsd(_read_first_sdf(pose_sdf, obabel))
        ref = _remove_hs_for_rmsd(_read_first_sdf(reference_sdf, obabel))
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            return None
        try:
            return float(rdMolAlign.GetBestRMS(pose, ref))
        except Exception:
            try:
                return float(rdMolAlign.CalcRMS(pose, ref))
            except Exception:
                return _ordered_heavy_rmsd(pose, ref)
    except Exception:
        return None


def run_posebusters(
    pose_sdf: str | Path,
    reference_sdf: str | Path,
    receptor_pdb: str | Path,
    *,
    rmsd: float | None,
    config: dict[str, Any],
    require_posebusters: bool = True,
) -> dict[str, Any]:
    entropy_cfg = config.get("entropy", {})
    tools_cfg = config.get("tools", {})
    command_template = str(entropy_cfg.get("posebusters_command", "") or "")
    executable = shutil.which(str(tools_cfg.get("posebusters", "bust"))) or str(tools_cfg.get("posebusters", "bust"))
    if not command_template or not executable:
        if require_posebusters:
            raise RuntimeError("PoseBusters command is required for success_pb; set tools.posebusters and entropy.posebusters_command.")
        return {
            "success_pb": False,
            "all_passed": False,
            "failures": ["posebusters_command_missing"],
            "has_posebusters": False,
        }

    pose_path = Path(pose_sdf)
    obabel = str(tools_cfg.get("obabel", "")) or None
    pose_input = _normalise_sdf_for_posebusters(pose_sdf, ".posebusters_pose.sdf", obabel)
    reference_input = _normalise_sdf_for_posebusters(reference_sdf, ".posebusters_reference.sdf", obabel)
    out_csv = pose_path.with_suffix(".posebusters.csv")
    out_json = pose_path.with_suffix(".posebusters.json")
    log_path = pose_path.with_suffix(".posebusters.log")
    context = {
        "posebusters": executable,
        "pose_sdf": str(pose_input),
        "reference_sdf": str(reference_input),
        "receptor_pdb": str(receptor_pdb),
        "posebusters_config": str(entropy_cfg.get("posebusters_config", "redock")),
        "out_csv": str(out_csv),
        "out_json": str(out_json),
    }
    args = [part.format_map(context) for part in shlex.split(command_template)]
    try:
        result = run_command(args, log_path=log_path, check=False)
    except Exception as exc:
        if require_posebusters:
            raise RuntimeError(f"PoseBusters command failed to launch for {pose_sdf}: {exc}") from exc
        return {
            "success_pb": False,
            "all_passed": False,
            "failures": ["posebusters_launch_failed"],
            "has_posebusters": False,
        }
    if result.returncode != 0:
        if require_posebusters:
            raise RuntimeError(f"PoseBusters failed for {pose_sdf}; see {log_path}")
        return {
            "success_pb": False,
            "all_passed": False,
            "failures": ["posebusters_failed"],
            "has_posebusters": False,
        }

    if out_csv.exists():
        parsed = _read_posebusters_output(out_csv)
    elif out_json.exists():
        parsed = _read_posebusters_output(out_json)
    else:
        parsed = _read_posebusters_output(log_path)
    all_passed = parsed["all_passed"]
    failures = parsed["failures"]
    success = bool(rmsd is not None and rmsd <= 2.0 and all_passed)
    return {
        "success_pb": success,
        "all_passed": all_passed,
        "failures": failures,
        "has_posebusters": True,
    }


def _read_posebusters_output(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"all_passed": False, "failures": ["posebusters_output_missing"]}
    text = path.read_text(errors="ignore").strip()
    if not text:
        return {"all_passed": False, "failures": ["posebusters_output_empty"]}
    if path.suffix.lower() == ".csv" or text.splitlines()[0].count(",") >= 1:
        return _parse_posebusters_csv(text)
    try:
        return _parse_posebusters_json(json.loads(text))
    except json.JSONDecodeError:
        return {"all_passed": False, "failures": ["posebusters_output_unparseable"]}


def _parse_posebusters_csv(text: str) -> dict[str, Any]:
    rows = list(csv.DictReader(text.splitlines()))
    if not rows:
        return {"all_passed": False, "failures": ["posebusters_csv_empty"]}
    row = rows[0]
    failures: list[str] = []
    bool_values: list[bool] = []
    metadata_keys = {"file", "molecule", "position", "mol_pred", "mol_true", "mol_cond"}
    for key, value in row.items():
        if str(key).strip().lower() in metadata_keys:
            continue
        lowered = str(value).strip().lower()
        if lowered in {"true", "1", "yes", "pass", "passed"}:
            bool_values.append(True)
        elif lowered in {"false", "0", "no", "fail", "failed"}:
            bool_values.append(False)
            failures.append(key)
    return {"all_passed": bool(bool_values) and all(bool_values), "failures": failures}


def _parse_posebusters_json(data: Any) -> dict[str, Any]:
    failures: list[str] = []
    values: list[bool] = []

    def walk(obj: Any, path: str = "") -> None:
        if isinstance(obj, bool):
            values.append(obj)
            if not obj:
                failures.append(path or "posebusters_check")
        elif isinstance(obj, dict):
            for key, value in obj.items():
                if str(key).lower() in {"file", "path", "molecule", "mol_pred", "mol_true", "mol_cond"}:
                    continue
                walk(value, f"{path}.{key}" if path else str(key))
        elif isinstance(obj, list):
            for idx, value in enumerate(obj):
                walk(value, f"{path}[{idx}]")

    walk(data)
    return {"all_passed": bool(values) and all(values), "failures": failures}
