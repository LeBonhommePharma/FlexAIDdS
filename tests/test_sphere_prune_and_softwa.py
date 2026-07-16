#!/usr/bin/env python3
"""Tests for LOCCLF sphere prune + classic SOFTWA emission (CF competitiveness prep)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"


def _load_gen():
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))
    path = SCRIPTS / "generate_flexaid_inp.py"
    spec = importlib.util.spec_from_file_location("generate_flexaid_inp", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["generate_flexaid_inp"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def gen():
    return _load_gen()


def test_matrix_pin_is_repo_validated(gen):
    assert gen.MATRIX_MD5_PIN == "9dc93717dfed0698006d88dd6a9627bc"


def test_soft_wall_default_and_config_keyword(gen, tmp_path, monkeypatch):
    monkeypatch.delenv("FLEXAIDDS_SOFT_WALL", raising=False)
    assert abs(gen.resolve_soft_wall_cutoff() - 0.40) < 1e-9
    monkeypatch.setenv("FLEXAIDDS_SOFT_WALL", "0.0")
    assert gen.resolve_soft_wall_cutoff() == 0.0
    monkeypatch.delenv("FLEXAIDDS_SOFT_WALL", raising=False)

    cfg = tmp_path / "CONFIG.inp"
    gen.write_config(
        cfg,
        target_pdb=tmp_path / "t.pdb",
        ligand_inp=tmp_path / "l.inp",
        spheres_pdb=tmp_path / "s.pdb",
        matrix_path=tmp_path / "m.dat",
        depspa=tmp_path,
        statep=tmp_path / "state",
        tempop=tmp_path / "tmp",
        rmsd_ref=None,
        temper=0,
        n_flex=0,
        soft_wall_cutoff=0.40,
    )
    text = cfg.read_text()
    assert "SOFTWA 0.40" in text
    assert "TEMPER 0" in text
    assert "CLUSTA CF" in text


def test_prune_spheres_near_ligand(gen, tmp_path):
    lig = tmp_path / "LIG_ref.pdb"
    lig.write_text(
        "HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O   LIG A   1       1.200   0.000   0.000  1.00  0.00           O\n"
    )
    sph = tmp_path / "sph.pdb"
    # near (kept), mid (kept if max large), far (pruned by distance)
    lines = ["REMARK spheres\n"]
    for i, (x, y, z) in enumerate(
        [
            (0.1, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (3.0, 0.0, 0.0),
            (20.0, 0.0, 0.0),
            (30.0, 0.0, 0.0),
        ],
        start=1,
    ):
        lines.append(
            f"ATOM  {i:5d}  SPH SURF    1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  1.50           S\n"
        )
    sph.write_text("".join(lines))
    info = gen.prune_spheres_near_ligand(
        sph, lig, max_dist=5.0, max_spheres=2
    )
    assert info["pruned"] is True
    assert info["n_in"] == 5
    assert info["n_out"] == 2
    out_lines = [
        ln for ln in sph.read_text().splitlines() if ln.startswith("ATOM")
    ]
    assert len(out_lines) == 2
    # nearest-first: 0.1 and 2.0 kept
    assert "0.100" in out_lines[0] or "  0.100" in out_lines[0]


def test_arm_c_requires_oracle_flag(gen):
    assert gen.ARM_SPEC["C"]["requires_oracle_pass"] is True
    assert gen.ARM_SPEC["B0"].get("science_control") is False
