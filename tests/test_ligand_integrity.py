#!/usr/bin/env python3
"""Unit tests for scripts/ligand_integrity.py + validate_ligand_integrity CLI."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    # Ensure scripts/ is importable for sibling modules
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def lig():
    return _load("ligand_integrity", SCRIPTS / "ligand_integrity.py")


@pytest.fixture(scope="module")
def cli():
    return _load("validate_ligand_integrity", SCRIPTS / "validate_ligand_integrity.py")


def _lig_ref_18() -> str:
    """Minimal LIG_ref-style PDB (18 heavy atoms + CONECT), inspired by 1P62 counts."""
    lines = []
    # 18 carbons in a chain with 1.5 Å bonds
    for i in range(18):
        ser = 90000 + i
        x = float(i) * 1.5
        lines.append(
            f"HETATM{ser:5d} C    LIG  9999    {x:8.3f}   0.000   0.000  1.00  1.00           C"
        )
    for i in range(18):
        ser = 90000 + i
        if i == 0:
            lines.append(f"CONECT{ser:5d}{ser+1:5d}")
        elif i == 17:
            lines.append(f"CONECT{ser:5d}{ser-1:5d}")
        else:
            lines.append(f"CONECT{ser:5d}{ser-1:5d}{ser+1:5d}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def _pose_missing_one() -> str:
    """17 heavy atoms (missing 90017) — 1P62-style emission bug fixture."""
    lines = []
    for i in range(17):
        ser = 90000 + i
        x = float(i) * 1.5
        lines.append(
            f"HETATM{ser:5d} C    GEO     1    {x:8.3f}   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def _split_ligand_1t40_style() -> str:
    """Two islands with CONECT claiming a bond across a huge gap (1T40 torn ligand)."""
    lines = [
        "HETATM90000 C    LIG  9999      0.000   0.000   0.000  1.00  1.00           C",
        "HETATM90001 C    LIG  9999      1.500   0.000   0.000  1.00  1.00           C",
        "HETATM90002 C    LIG  9999     20.000   0.000   0.000  1.00  1.00           C",
        "HETATM90003 C    LIG  9999     21.500   0.000   0.000  1.00  1.00           C",
        "CONECT9000090001",
        "CONECT900019000090002",  # false bond 90001-90002 ~18.5 Å
        "CONECT900029000190003",
        "CONECT9000390002",
        "END",
    ]
    return "\n".join(lines) + "\n"


def test_count_heavy_lig_ref(lig, tmp_path: Path):
    p = tmp_path / "LIG_ref.pdb"
    p.write_text(_lig_ref_18())
    n, atoms = lig.count_heavy_atoms_pdb(p, ligand_only=False)
    assert n == 18
    assert len(atoms) == 18


def test_match_ok(lig, tmp_path: Path):
    ref = tmp_path / "LIG_ref.pdb"
    pose = tmp_path / "pose.pdb"
    text = _lig_ref_18()
    ref.write_text(text)
    pose.write_text(text)
    res = lig.validate_ligand_integrity(ref, pose, max_bond=3.0)
    assert res.ok
    assert res.exit_code == lig.EXIT_OK
    assert res.ref_heavy == 18
    assert res.pose_heavy == 18


def test_count_mismatch_1p62_style(lig, tmp_path: Path):
    ref = tmp_path / "LIG_ref.pdb"
    pose = tmp_path / "1P62_0.pdb"
    ref.write_text(_lig_ref_18())
    pose.write_text(_pose_missing_one())
    res = lig.validate_ligand_integrity(ref, pose, check_bonds=False)
    assert not res.ok
    assert res.exit_code == lig.EXIT_COUNT_MISMATCH
    assert res.ref_heavy == 18
    assert res.pose_heavy == 17
    assert 90017 in res.missing_serials_in_pose


def test_bond_breach_1t40_style(lig, tmp_path: Path):
    ref = tmp_path / "LIG_ref.pdb"
    pose = tmp_path / "torn.pdb"
    torn = _split_ligand_1t40_style()
    ref.write_text(torn)
    pose.write_text(torn)
    res = lig.validate_ligand_integrity(ref, pose, max_bond=3.0, check_bonds=True)
    assert not res.ok
    assert res.exit_code == lig.EXIT_BOND_TOO_LONG
    assert res.max_bond_a is not None and res.max_bond_a > 10.0
    assert any(b.distance > 3.0 for b in res.bond_breaches)


def test_work_dir_prep_only(lig, tmp_path: Path):
    work = tmp_path / "1P62"
    work.mkdir()
    (work / "LIG_ref.pdb").write_text(_lig_ref_18())
    res = lig.validate_work_dir(work, require_ini=False)
    assert res.ok
    assert any("no INI yet" in m for m in res.messages)


def test_work_dir_require_ini_missing(lig, tmp_path: Path):
    work = tmp_path / "1P62"
    work.mkdir()
    (work / "LIG_ref.pdb").write_text(_lig_ref_18())
    res = lig.validate_work_dir(work, require_ini=True)
    assert not res.ok
    assert res.exit_code == lig.EXIT_MISSING_POSE


def test_work_dir_with_ini(lig, tmp_path: Path):
    work = tmp_path / "1P62"
    work.mkdir()
    text = _lig_ref_18()
    (work / "LIG_ref.pdb").write_text(text)
    (work / "1P62_INI.pdb").write_text(text)
    res = lig.validate_work_dir(work, require_ini=True)
    assert res.ok
    assert res.pose_heavy == 18


def test_cli_exit_codes(cli, tmp_path: Path):
    ref = tmp_path / "LIG_ref.pdb"
    pose = tmp_path / "pose.pdb"
    ref.write_text(_lig_ref_18())
    pose.write_text(_pose_missing_one())
    rc = cli.main(["--ref", str(ref), "--pose", str(pose), "--no-bonds", "-q"])
    assert rc == 12  # COUNT_MISMATCH


def test_cli_json(cli, tmp_path: Path):
    ref = tmp_path / "LIG_ref.pdb"
    ref.write_text(_lig_ref_18())
    out = tmp_path / "out.json"
    rc = cli.main(["--ref", str(ref), "--json", str(out), "-q"])
    assert rc == 0
    data = json.loads(out.read_text())
    assert data["ok"] is True
    assert data["ref_heavy"] == 18


def test_is_heavy_element(lig):
    assert lig.is_heavy_element("C") is True
    assert lig.is_heavy_element("N") is True
    assert lig.is_heavy_element("CL") is True
    assert lig.is_heavy_element("H") is False
    assert lig.is_heavy_element("HE") is True
    assert lig.is_heavy_element("") is False
