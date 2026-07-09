from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import shlex
from typing import Any

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

from .io_utils import (
    _local_scratch_root,
    _materialize_local,
    _safe_exists,
    _safe_read_text,
    _safe_stat_size,
    _safe_unlink,
    _safe_write_text,
    run_command,
)
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


def _heavy_formula(mol: Chem.Mol) -> tuple[tuple[str, int], ...]:
    counts: dict[str, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1:
            continue
        symbol = atom.GetSymbol()
        counts[symbol] = counts.get(symbol, 0) + 1
    return tuple(sorted(counts.items()))


def _whole_ligand_compatible(pose: Chem.Mol, ref: Chem.Mol) -> bool:
    """Require complete whole-ligand topology for success-bearing RMSD."""
    if pose.GetNumAtoms() != ref.GetNumAtoms():
        return False
    if pose.GetNumAtoms() == 0:
        return False
    return _heavy_formula(pose) == _heavy_formula(ref)


def _mcs_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float | None:
    """Diagnostic-only MCS/fragment RMSD. Never use for success_pb."""
    if pose.GetNumConformers() == 0 or ref.GetNumConformers() == 0:
        return None
    mcs = rdFMCS.FindMCS([pose, ref], timeout=30)
    if mcs.numAtoms <= 0:
        return None
    pattern = Chem.MolFromSmarts(mcs.smartsString)
    if pattern is None:
        return None
    pose_match = pose.GetSubstructMatch(pattern)
    ref_match = ref.GetSubstructMatch(pattern)
    if not pose_match or not ref_match or len(pose_match) != len(ref_match):
        return None
    atom_map = list(zip(pose_match, ref_match, strict=True))
    return float(AllChem.AlignMol(pose, ref, atomMap=atom_map))


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


def _posebusters_scratch(pose_path: Path) -> Path:
    scratch = _local_scratch_root() / "posebusters" / pose_path.stem
    scratch.mkdir(parents=True, exist_ok=True)
    return scratch


def _file_fingerprint(path: Path) -> dict[str, Any]:
    try:
        st = path.stat()
        return {
            "path": str(path),
            "mtime_ns": int(getattr(st, "st_mtime_ns", int(st.st_mtime * 1e9))),
            "size": int(st.st_size),
        }
    except OSError:
        return {"path": str(path), "mtime_ns": -1, "size": -1}


def _executable_fingerprint(executable: str) -> dict[str, Any]:
    path = Path(executable)
    if path.is_file():
        return _file_fingerprint(path)
    resolved = shutil.which(executable)
    if resolved:
        return _file_fingerprint(Path(resolved))
    return {"path": executable, "mtime_ns": -1, "size": -1}


def _normalise_sdf_for_posebusters(path: str | Path, suffix: str, obabel: str | None = None) -> Path:
    source = _materialize_local(Path(path))
    out_path = _posebusters_scratch(source) / f"{source.stem}{suffix}"
    executable = obabel or shutil.which("obabel")
    if not executable:
        return source

    source_fp = _file_fingerprint(source)
    meta_path = out_path.with_suffix(out_path.suffix + ".meta.json")
    reuse = False
    if _safe_exists(out_path) and _safe_stat_size(out_path) > 0 and _safe_exists(meta_path):
        try:
            meta = json.loads(_safe_read_text(meta_path))
            reuse = meta.get("source") == source_fp and meta.get("obabel") == executable
        except (json.JSONDecodeError, OSError, TypeError):
            reuse = False
    if not reuse:
        if _safe_exists(out_path):
            _safe_unlink(out_path)
        run_command(
            [executable, "-isdf", str(source), "-osdf", "-O", str(out_path)],
            log_path=out_path.with_suffix(".obabel.log"),
            check=False,
        )
        if _safe_exists(out_path):
            normalise_v2000_counts_file(out_path)
            _safe_write_text(
                meta_path,
                json.dumps({"source": source_fp, "obabel": executable}, indent=2) + "\n",
            )
    return out_path if _safe_exists(out_path) and _safe_stat_size(out_path) > 0 else source


def rmsd_diagnostics(
    pose_sdf: str | Path,
    reference_sdf: str | Path,
    *,
    obabel: str | None = None,
) -> dict[str, Any]:
    """
    Compute whole-ligand RMSD (success-bearing) and diagnostic fragment MCS RMSD.

    success_pb must use only `rmsd_A` when `whole_ligand` is True. Fragment MCS
    values are diagnostics only — truncated ligands can produce false sub-2 A
    fragment scores (Astex Diverse success definition requires complete ligand).
    """
    result: dict[str, Any] = {
        "rmsd_A": None,
        "fragment_mcs_rmsd_A": None,
        "whole_ligand": False,
        "pose_heavy_atoms": None,
        "ref_heavy_atoms": None,
        "formula_match": False,
    }
    try:
        pose = _remove_hs_for_rmsd(_read_first_sdf(_materialize_local(Path(pose_sdf)), obabel))
        ref = _remove_hs_for_rmsd(_read_first_sdf(_materialize_local(Path(reference_sdf)), obabel))
    except Exception:
        return result

    result["pose_heavy_atoms"] = int(pose.GetNumAtoms())
    result["ref_heavy_atoms"] = int(ref.GetNumAtoms())
    result["formula_match"] = _heavy_formula(pose) == _heavy_formula(ref)
    result["whole_ligand"] = _whole_ligand_compatible(pose, ref)

    # Diagnostic MCS always attempted (does not feed success_pb).
    try:
        mcs_value = _mcs_rmsd(pose, ref)
        if mcs_value is not None and math.isfinite(float(mcs_value)):
            result["fragment_mcs_rmsd_A"] = float(mcs_value)
    except Exception:
        pass

    if not result["whole_ligand"]:
        return result

    for calculator in (
        lambda: float(rdMolAlign.GetBestRMS(pose, ref)),
        lambda: float(rdMolAlign.CalcRMS(pose, ref)),
        lambda: _ordered_heavy_rmsd(pose, ref),
    ):
        try:
            value = calculator()
        except Exception:
            value = None
        if value is not None and math.isfinite(float(value)):
            result["rmsd_A"] = float(value)
            return result
    return result


def rmsd_to_reference(
    pose_sdf: str | Path,
    reference_sdf: str | Path,
    *,
    obabel: str | None = None,
) -> float | None:
    """Whole-ligand RMSD only. Returns None for incomplete/truncated poses."""
    return rmsd_diagnostics(pose_sdf, reference_sdf, obabel=obabel).get("rmsd_A")


def _posebusters_cache_meta_path(scratch: Path) -> Path:
    return scratch / "posebusters.cache.json"


def _build_posebusters_cache_meta(
    *,
    pose_sdf: Path,
    reference_sdf: Path,
    receptor_pdb: Path,
    pose_input: Path,
    reference_input: Path,
    executable: str,
    command_template: str,
    args: list[str],
    config_name: str,
) -> dict[str, Any]:
    return {
        "version": 1,
        "pose_sdf": _file_fingerprint(pose_sdf),
        "reference_sdf": _file_fingerprint(reference_sdf),
        "receptor_pdb": _file_fingerprint(receptor_pdb),
        "pose_input": _file_fingerprint(pose_input),
        "reference_input": _file_fingerprint(reference_input),
        "posebusters_executable": _executable_fingerprint(executable),
        "command_template": command_template,
        "command_args": args,
        "posebusters_config": config_name,
    }


def _posebusters_cache_is_valid(
    meta_path: Path,
    expected: dict[str, Any],
    out_csv: Path,
    out_json: Path,
) -> bool:
    if not _safe_exists(meta_path):
        return False
    if not (
        (_safe_exists(out_csv) and _safe_stat_size(out_csv) > 0)
        or (_safe_exists(out_json) and _safe_stat_size(out_json) > 0)
    ):
        return False
    try:
        meta = json.loads(_safe_read_text(meta_path))
    except (json.JSONDecodeError, OSError, TypeError):
        return False
    for key in (
        "pose_sdf",
        "reference_sdf",
        "receptor_pdb",
        "pose_input",
        "reference_input",
        "posebusters_executable",
        "command_template",
        "command_args",
        "posebusters_config",
    ):
        if meta.get(key) != expected.get(key):
            return False
    return True


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

    pose_path = _materialize_local(Path(pose_sdf))
    reference_path = _materialize_local(Path(reference_sdf))
    receptor_path = _materialize_local(Path(receptor_pdb))
    obabel = str(tools_cfg.get("obabel", "")) or None
    pose_input = _normalise_sdf_for_posebusters(pose_path, ".posebusters_pose.sdf", obabel)
    reference_input = _normalise_sdf_for_posebusters(reference_path, ".posebusters_reference.sdf", obabel)
    scratch = _posebusters_scratch(pose_path)
    out_csv = scratch / "posebusters.csv"
    out_json = scratch / "posebusters.json"
    log_path = scratch / "posebusters.log"
    meta_path = _posebusters_cache_meta_path(scratch)
    config_name = str(entropy_cfg.get("posebusters_config", "redock"))

    for stale in (out_csv, out_json):
        if _safe_exists(stale) and _safe_stat_size(stale) == 0:
            _safe_unlink(stale)

    context = {
        "posebusters": executable,
        "pose_sdf": str(pose_input),
        "reference_sdf": str(reference_input),
        "receptor_pdb": str(receptor_path),
        "posebusters_config": config_name,
        "out_csv": str(out_csv),
        "out_json": str(out_json),
    }
    args = [part.format_map(context) for part in shlex.split(command_template)]
    expected_meta = _build_posebusters_cache_meta(
        pose_sdf=pose_path,
        reference_sdf=reference_path,
        receptor_pdb=receptor_path,
        pose_input=pose_input,
        reference_input=reference_input,
        executable=executable,
        command_template=command_template,
        args=args,
        config_name=config_name,
    )

    if _posebusters_cache_is_valid(meta_path, expected_meta, out_csv, out_json):
        cached = out_csv if _safe_stat_size(out_csv) > 0 else out_json
        parsed = _read_posebusters_output(cached)
        if "posebusters_output_unreadable" not in parsed["failures"]:
            all_passed = parsed["all_passed"]
            failures = parsed["failures"]
            success = bool(rmsd is not None and rmsd <= 2.0 and all_passed)
            return {
                "success_pb": success,
                "all_passed": all_passed,
                "failures": failures,
                "has_posebusters": True,
                "cache_hit": True,
            }

    # Invalidate stale caches when inputs/command changed.
    for path in (out_csv, out_json, meta_path):
        if _safe_exists(path):
            _safe_unlink(path)

    timeout = int(entropy_cfg.get("posebusters_timeout_seconds", 300))
    try:
        result = run_command(args, log_path=log_path, check=False, timeout=timeout)
    except Exception as exc:
        _safe_write_text(log_path, f"PoseBusters command error for {pose_sdf}: {exc}\n")
        return {
            "success_pb": False,
            "all_passed": False,
            "failures": ["posebusters_command_error"],
            "has_posebusters": True,
            "cache_hit": False,
        }
    if result.returncode != 0:
        return {
            "success_pb": False,
            "all_passed": False,
            "failures": ["posebusters_failed"],
            "has_posebusters": True,
            "cache_hit": False,
        }

    if _safe_exists(out_csv) and _safe_stat_size(out_csv) > 0:
        parsed = _read_posebusters_output(out_csv)
    elif _safe_exists(out_json) and _safe_stat_size(out_json) > 0:
        parsed = _read_posebusters_output(out_json)
    else:
        parsed = _read_posebusters_output(log_path)

    if _safe_stat_size(out_csv) > 0 or _safe_stat_size(out_json) > 0:
        _safe_write_text(meta_path, json.dumps(expected_meta, indent=2) + "\n")

    all_passed = parsed["all_passed"]
    failures = parsed["failures"]
    success = bool(rmsd is not None and rmsd <= 2.0 and all_passed)
    return {
        "success_pb": success,
        "all_passed": all_passed,
        "failures": failures,
        "has_posebusters": True,
        "cache_hit": False,
    }


def _read_posebusters_output(path: Path) -> dict[str, Any]:
    if not _safe_exists(path):
        return {"all_passed": False, "failures": ["posebusters_output_missing"]}
    text = _safe_read_text(path).strip()
    if not text:
        return {"all_passed": False, "failures": ["posebusters_output_empty", "posebusters_output_unreadable"]}
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
