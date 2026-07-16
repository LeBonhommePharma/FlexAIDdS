#!/usr/bin/env python3
"""Unit tests for scripts/clean_target_apo.py."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"


def _load():
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    path = SCRIPTS / "clean_target_apo.py"
    spec = importlib.util.spec_from_file_location("clean_target_apo", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["clean_target_apo"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cta():
    return _load()


TINY_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.400   2.400   0.000  1.00  0.00           O
HETATM    5 MG    MG A 401      5.000   5.000   5.000  1.00 30.00          MG
HETATM    6  O   HOH A 402      6.000   6.000   6.000  1.00 20.00           O
HETATM    7  O   HOH A 403      7.000   7.000   7.000  1.00 20.00           O
HETATM    8 FE    FE A 404      8.000   8.000   8.000  1.00 30.00          FE
HETATM    9 FE   HEM A 500      9.000   9.000   9.000  1.00 20.00          FE
HETATM   10  NA  HEM A 500      9.500   9.500   9.500  1.00 20.00           N
CONECT    5    6
CONECT    6    5
CONECT    1    2
END
"""


def test_strip_water_and_metals_default(cta):
    out, rep = cta.clean_apo_pdb(TINY_PDB, keep_hoh=False, keep_metals=False)
    assert rep.waters_removed == 2
    assert rep.metals_removed >= 2  # MG + FE ion (HEM Fe may stay with res HEM)
    assert "HOH" not in out
    assert "MG A 401" not in out and "MG    MG" not in out
    # Protein kept
    assert "ALA A   1" in out
    # HEM kept by default (cofactor)
    assert "HEM" in out
    # Orphan CONECT to removed MG/water dropped
    assert "CONECT    5" not in out
    # Protein CONECT kept
    assert "CONECT    1" in out
    assert rep.atoms_out < rep.atoms_in


def test_keep_hoh(cta):
    out, rep = cta.clean_apo_pdb(TINY_PDB, keep_hoh=True, keep_metals=False)
    assert rep.waters_removed == 0
    assert "HOH" in out
    assert rep.metals_removed >= 1


def test_keep_metals(cta):
    out, rep = cta.clean_apo_pdb(TINY_PDB, keep_hoh=False, keep_metals=True)
    assert rep.metals_removed == 0
    assert "MG" in out
    assert rep.waters_removed == 2


def test_file_roundtrip(cta, tmp_path: Path):
    src = tmp_path / "in.pdb"
    dst = tmp_path / "out.pdb"
    src.write_text(TINY_PDB)
    rep = cta.clean_apo_file(src, dst, keep_hoh=False, keep_metals=False)
    text = dst.read_text()
    assert rep.waters_removed == 2
    assert "HOH" not in text
    assert "ALA" in text


def test_cli_dry_run(cta, tmp_path: Path, monkeypatch):
    src = tmp_path / "in.pdb"
    src.write_text(TINY_PDB)
    rc = cta.main([str(src), "--in-place", "--dry-run", "-q"])
    assert rc == 0
    # Unchanged
    assert src.read_text() == TINY_PDB


def test_env_keep_hoh(cta, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_KEEP_HOH", "1")
    monkeypatch.delenv("FLEXAIDDS_KEEP_METALS", raising=False)
    hoh, metals = cta.resolve_keep_flags()
    assert hoh is True
    assert metals is False
