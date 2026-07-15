"""Benchmark data path resolution for DatasetRunner."""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Tuple

from .astex_targets import (
    ASTEX_NONNATIVE_BY_NAME,
    lookup_nonnative_family,
    parse_crossdock_entry_id,
)

_RMSD_REMARK_RE = re.compile(
    r"REMARK\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+RMSD\s+to\s+ref",
    re.IGNORECASE,
)
_CF_REMARK_RE = re.compile(r"REMARK\s+CF=([+-]?(?:\d+(?:\.\d*)?|\.\d+))", re.IGNORECASE)
_CF_APP_REMARK_RE = re.compile(
    r"REMARK\s+CF\.app=([+-]?(?:\d+(?:\.\d*)?|\.\d+))", re.IGNORECASE
)


def _first_existing(paths: List[Path]) -> Optional[Path]:
    for p in paths:
        if p.is_file():
            return p
    return None


def _structure_probes(data_dir: Path, pdb_id: str, state: str = "") -> List[Path]:
    """Candidate receptor PDB paths for a PDB-code directory."""
    tid = pdb_id.upper()
    tdir = data_dir / tid
    probes = [
        tdir / f"{tid}_{state}.pdb" if state else tdir / f"{tid}.pdb",
        tdir / f"{tid}.pdb",
        tdir / f"{tid}_holo.pdb",
        tdir / f"{tid}_protein.pdb",
        tdir / "receptor.pdb",
        data_dir / f"{tid}.pdb",
    ]
    if state and state != "holo":
        probes.insert(0, tdir / f"{tid}_{state}.pdb")
    return probes


def _ligand_probes(data_dir: Path, pdb_id: str) -> List[Path]:
    tid = pdb_id.upper()
    tdir = data_dir / tid
    return [
        tdir / f"{tid}_ligand.sdf",
        tdir / f"{tid}_ligand.mol2",
        tdir / f"{tid}.sdf",
        tdir / f"{tid}.mol2",
    ]


def find_structure_pdb(data_dir: Path, pdb_id: str, state: str = "") -> Optional[Path]:
    """Locate a receptor/structure PDB under ``data_dir``."""
    return _first_existing(_structure_probes(data_dir, pdb_id, state))


def find_ligand_file(data_dir: Path, pdb_id: str) -> Optional[Path]:
    """Locate the reference ligand for a PDB-code directory."""
    return _first_existing(_ligand_probes(data_dir, pdb_id))


def resolve_astex_nonnative_paths(
    data_dir: Path,
    entry_id: str,
    state: str,
) -> Tuple[Optional[Path], List[Path]]:
    """Resolve receptor + ligand paths for astex_nonnative entries.

    Supports:
    - Catalog ids: ``1G9V_1EVE`` with state ``crossdock``
    - YAML family ids: ``ACE`` with state ``holo`` / ``apo`` / ``alternative``
    """
    cross = parse_crossdock_entry_id(entry_id)
    if cross is not None or state == "crossdock":
        if cross is None:
            return None, []
        native_pdb, receptor_pdb = cross
        receptor = find_structure_pdb(data_dir, receptor_pdb)
        ligand = find_ligand_file(data_dir, native_pdb)
        return receptor, [ligand] if ligand else []

    family = lookup_nonnative_family(entry_id)
    if family is None:
        return None, []

    native = family.native_pdb
    ligand = find_ligand_file(data_dir, native)
    ligands = [ligand] if ligand else []

    if state == "holo":
        receptor = find_structure_pdb(data_dir, native, state="holo")
        return receptor, ligands

    if state == "apo":
        receptor = find_structure_pdb(data_dir, native, state="apo")
        if receptor is None:
            for alt in family.alternative_pdbs:
                receptor = find_structure_pdb(data_dir, alt, state="apo")
                if receptor is not None:
                    break
        return receptor, ligands

    if state == "alternative":
        for alt in family.alternative_pdbs:
            if alt.upper() == native.upper():
                continue
            alt_dir = data_dir / alt.upper()
            if not alt_dir.is_dir():
                continue
            receptor = find_structure_pdb(data_dir, alt)
            if receptor is not None:
                return receptor, ligands

    return None, ligands


def resolve_benchmark_paths(
    slug: str,
    data_dir: Path,
    entry_id: str,
    state: str,
) -> Tuple[Optional[Path], List[Path]]:
    """Dataset-aware receptor/ligand resolution."""
    if slug == "astex_nonnative":
        return resolve_astex_nonnative_paths(data_dir, entry_id, state)

    # Default: PDB-id targets (astex_diverse, hap2, casf, etc.)
    receptor = find_structure_pdb(data_dir, entry_id, state=state)
    ligand = find_ligand_file(data_dir, entry_id)
    ligands: List[Path] = []
    if ligand is not None:
        ligands.append(ligand)
    else:
        target_dir = data_dir / entry_id.upper()
        if not target_dir.is_dir():
            target_dir = data_dir / entry_id
        if target_dir.is_dir():
            ligands = (
                list(target_dir.glob("*.mol2"))
                + list(target_dir.glob("*.sdf"))
            )
    return receptor, ligands


def validate_exec_path(path: str | Path) -> None:
    """Refuse paths that must never reach a shell or corrupt provenance.

    Raises ValueError on empty path, embedded NUL, newlines, or other C0
    control characters (TAB is allowed). Mirrors ``LIB/shell_exec.h``.
    """
    s = str(path)
    if not s:
        raise ValueError("path is empty")
    if '\x00' in s:
        raise ValueError("path contains NUL byte")
    for ch in s:
        if ch in ('\n', '\r'):
            raise ValueError("path contains newline/CR")
        o = ord(ch)
        if o < 0x20 and ch != '\t':
            raise ValueError("path contains control character")
        if o == 0x7F:
            raise ValueError("path contains control character")


def is_safe_exec_path(path: str | Path) -> bool:
    try:
        validate_exec_path(path)
        return True
    except ValueError:
        return False

