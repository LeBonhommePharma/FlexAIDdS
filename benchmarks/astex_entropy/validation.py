from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import shutil
import shlex
from typing import Any, Sequence

from rdkit import Chem
from rdkit import RDLogger

from .io_utils import run_command
from .sdf_utils import normalise_v2000_counts_file, read_first_sdf_mol

RDLogger.DisableLog("rdApp.error")

# Hard caps: unbounded substructure / backtracking can hang on bad graphs.
_MAX_GRAPH_MAPS = 64
_MAX_DFS_NODES = 50_000


def _read_first_sdf(path: str | Path, obabel: str | None = None) -> Chem.Mol:
    return read_first_sdf_mol(path, obabel, allow_unsanitized=True)


def _remove_hs_for_rmsd(mol: Chem.Mol) -> Chem.Mol:
    try:
        return Chem.RemoveHs(mol, sanitize=False)
    except TypeError:
        return Chem.RemoveHs(mol)
    except Exception:
        return mol


def _atom_signature(atom: Chem.Atom) -> tuple[int, int]:
    """Element + heavy-atom degree (stable without full sanitization)."""
    try:
        heavy_deg = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)
    except Exception:
        heavy_deg = atom.GetDegree()
    return (int(atom.GetAtomicNum()), int(heavy_deg))


def _element_sequence(mol: Chem.Mol) -> list[int]:
    return [int(a.GetAtomicNum()) for a in mol.GetAtoms()]


def _bond_order_key(mol: Chem.Mol, i: int, j: int) -> int | None:
    bond = mol.GetBondBetweenAtoms(i, j)
    if bond is None:
        return None
    try:
        return int(round(bond.GetBondTypeAsDouble() * 2.0))
    except Exception:
        return 2  # single


def _direct_rmsd_for_map(
    pose: Chem.Mol,
    ref: Chem.Mol,
    mapping: Sequence[int],
) -> float | None:
    """Unaligned (direct-coordinate) RMSD for a full pose->ref atom map."""
    if pose.GetNumConformers() == 0 or ref.GetNumConformers() == 0:
        return None
    n = pose.GetNumAtoms()
    if n == 0 or len(mapping) != n:
        return None
    if len(set(int(x) for x in mapping)) != n:
        return None
    pose_conf = pose.GetConformer()
    ref_conf = ref.GetConformer()
    total = 0.0
    for pose_i, ref_i in enumerate(mapping):
        ri = int(ref_i)
        if ri < 0 or ri >= ref.GetNumAtoms():
            return None
        if pose.GetAtomWithIdx(pose_i).GetAtomicNum() != ref.GetAtomWithIdx(ri).GetAtomicNum():
            return None
        p = pose_conf.GetAtomPosition(pose_i)
        r = ref_conf.GetAtomPosition(ri)
        total += (p.x - r.x) ** 2 + (p.y - r.y) ** 2 + (p.z - r.z) ** 2
    return math.sqrt(total / n)


def _ordered_identity_map(n: int) -> list[int]:
    return list(range(n))


def _maps_via_substruct(pose: Chem.Mol, ref: Chem.Mol) -> list[list[int]]:
    """Full-graph maps from RDKit isomorphism (capped)."""
    n = pose.GetNumAtoms()
    maps: list[list[int]] = []
    try:
        # uniquify=True keeps the first representative of each unique atom set;
        # still enumerate a few for symmetry with maxMatches.
        matches = ref.GetSubstructMatches(
            pose, uniquify=False, maxMatches=_MAX_GRAPH_MAPS
        )
        for m in matches:
            if len(m) == n and len(set(m)) == n:
                maps.append(list(m))
                if len(maps) >= _MAX_GRAPH_MAPS:
                    break
    except Exception:
        return maps
    return maps


def _maps_via_backtrack(pose: Chem.Mol, ref: Chem.Mol) -> list[list[int]]:
    """
    Bounded DFS for chemically valid full-graph bijections.

    Uses element + heavy degree + bond-order consistency. Aborts after
    _MAX_DFS_NODES expansions so graph-mismatch cases cannot hang.
    """
    n = pose.GetNumAtoms()
    pose_sig = [_atom_signature(a) for a in pose.GetAtoms()]
    ref_sig = [_atom_signature(a) for a in ref.GetAtoms()]

    # Quick reject: multiset of signatures must match.
    if sorted(pose_sig) != sorted(ref_sig):
        return []

    candidates: list[list[int]] = []
    for i in range(n):
        cands = [j for j in range(n) if pose_sig[i] == ref_sig[j]]
        if not cands:
            return []
        candidates.append(cands)

    # Fail fast when branching factor is combinatorial for large N.
    product = 1
    for c in candidates:
        product *= len(c)
        if product > 1_000_000 and n > 12:
            # Fall back to ordered map only if elements align; else empty.
            return []

    order = sorted(range(n), key=lambda i: len(candidates[i]))
    assignment = [-1] * n
    used = [False] * n
    maps: list[list[int]] = []
    nodes = 0

    def bonds_ok(pose_i: int, ref_j: int) -> bool:
        for bond in pose.GetAtomWithIdx(pose_i).GetBonds():
            other = bond.GetOtherAtomIdx(pose_i)
            if assignment[other] < 0:
                continue
            ref_other = assignment[other]
            po = _bond_order_key(pose, pose_i, other)
            ro = _bond_order_key(ref, ref_j, ref_other)
            if po is None or ro is None or po != ro:
                return False
        return True

    def dfs(depth: int) -> None:
        nonlocal nodes
        if len(maps) >= _MAX_GRAPH_MAPS:
            return
        nodes += 1
        if nodes > _MAX_DFS_NODES:
            return
        if depth == n:
            maps.append(assignment.copy())
            return
        pose_i = order[depth]
        for ref_j in candidates[pose_i]:
            if used[ref_j]:
                continue
            if not bonds_ok(pose_i, ref_j):
                continue
            used[ref_j] = True
            assignment[pose_i] = ref_j
            dfs(depth + 1)
            assignment[pose_i] = -1
            used[ref_j] = False
            if nodes > _MAX_DFS_NODES or len(maps) >= _MAX_GRAPH_MAPS:
                return

    dfs(0)
    return maps


def _enumerate_full_graph_maps(pose: Chem.Mol, ref: Chem.Mol) -> list[list[int]]:
    """Chemically valid full-graph pose→ref maps; never unbounded."""
    n = pose.GetNumAtoms()
    if n == 0 or n != ref.GetNumAtoms():
        return []

    # Fast path: identical element order — always include ordered map first.
    maps: list[list[int]] = []
    if _element_sequence(pose) == _element_sequence(ref):
        maps.append(_ordered_identity_map(n))

    # RDKit isomorphism (symmetry-aware). Cap matches.
    for m in _maps_via_substruct(pose, ref):
        if m not in maps:
            maps.append(m)
        if len(maps) >= _MAX_GRAPH_MAPS:
            return maps

    # If we already have ≥1 full map, skip expensive backtracking unless
    # symmetry search returned nothing beyond ordered (still OK for RMSD min).
    if maps:
        return maps

    # No ordered map and no RDKit map: try bounded signature DFS.
    for m in _maps_via_backtrack(pose, ref):
        if m not in maps:
            maps.append(m)
        if len(maps) >= _MAX_GRAPH_MAPS:
            break
    return maps


def direct_graph_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float | None:
    """
    Symmetry-aware whole-ligand RMSD without Kabsch/alignment.

    Enumerates chemically valid full-graph atom mappings and returns the
    minimum direct-coordinate RMSD. Returns None when no full-graph map exists
    (count mismatch, graph mismatch, or missing conformers).
    """
    if pose is None or ref is None:
        return None
    if pose.GetNumAtoms() == 0 or ref.GetNumAtoms() == 0:
        return None
    if pose.GetNumAtoms() != ref.GetNumAtoms():
        return None
    if pose.GetNumConformers() == 0 or ref.GetNumConformers() == 0:
        return None

    maps = _enumerate_full_graph_maps(pose, ref)
    if not maps:
        return None

    best: float | None = None
    for m in maps:
        val = _direct_rmsd_for_map(pose, ref, m)
        if val is None or not math.isfinite(val):
            continue
        if best is None or val < best:
            best = val
    return best


def _ordered_heavy_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float | None:
    """Legacy ordered heavy-atom helper (diagnostics / tests)."""
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
    """
    Direct-coordinate, chemically valid full-graph RMSD (no alignment).

    Does **not** use RDKit GetBestRMS / Kabsch. A pure translation reports
    non-zero RMSD. Graph-incompatible atom sets return None.
    """
    try:
        pose = _remove_hs_for_rmsd(_read_first_sdf(pose_sdf, obabel))
        ref = _remove_hs_for_rmsd(_read_first_sdf(reference_sdf, obabel))
        return direct_graph_rmsd(pose, ref)
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
