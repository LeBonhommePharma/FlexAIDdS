"""Graph-identity guards for benchmarks/astex_entropy/validation.py.

Ordered identity must require full bond-order graph equality.
One-directional subgraph matches must be rejected.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

from benchmarks.astex_entropy.validation import (  # noqa: E402
    _enumerate_full_graph_maps,
    _full_graph_bond_order_equal,
    direct_graph_rmsd,
)


def _mol_xyz(smiles: str, conf_xyz: list[tuple[float, float, float]]) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    mol = Chem.AddHs(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(conf_xyz):
        conf.SetAtomPosition(i, (x, y, z))
    mol.AddConformer(conf, assignId=True)
    return mol


def test_ordered_identity_requires_full_graph_bond_order():
    # Same heavy elements in order but different bond orders: C-C vs C=C
    ethane = Chem.MolFromSmiles("CC")
    ethene = Chem.MolFromSmiles("C=C")
    assert ethane is not None and ethene is not None
    # Force same atom count without H for simplicity
    assert ethane.GetNumAtoms() == ethene.GetNumAtoms() == 2
    conf_a = Chem.Conformer(2)
    conf_a.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf_a.SetAtomPosition(1, (1.5, 0.0, 0.0))
    conf_b = Chem.Conformer(2)
    conf_b.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf_b.SetAtomPosition(1, (1.3, 0.0, 0.0))
    ethane.AddConformer(conf_a, assignId=True)
    ethene.AddConformer(conf_b, assignId=True)

    identity = [0, 1]
    assert not _full_graph_bond_order_equal(ethane, ethene, identity)
    maps = _enumerate_full_graph_maps(ethane, ethene)
    assert maps == []
    assert direct_graph_rmsd(ethane, ethene) is None


def test_rejects_one_directional_subgraph_match():
    # Propane as "pose" vs ethane as "ref" — different atom counts → no map.
    pose = Chem.MolFromSmiles("CCC")
    ref = Chem.MolFromSmiles("CC")
    assert pose is not None and ref is not None
    maps = _enumerate_full_graph_maps(pose, ref)
    assert maps == []


def test_full_graph_identity_allows_rmsd():
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    conf0 = Chem.Conformer(mol.GetNumAtoms())
    conf1 = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf0.SetAtomPosition(i, (float(i), 0.0, 0.0))
        conf1.SetAtomPosition(i, (float(i) + 0.1, 0.0, 0.0))
    mol.AddConformer(conf0, assignId=True)
    # Clone for ref with same graph, shifted coords
    ref = Chem.Mol(mol)
    ref.RemoveAllConformers()
    ref.AddConformer(conf1, assignId=True)
    maps = _enumerate_full_graph_maps(mol, ref)
    assert maps
    assert any(_full_graph_bond_order_equal(mol, ref, m) for m in maps)
    rmsd = direct_graph_rmsd(mol, ref)
    assert rmsd is not None and math.isfinite(rmsd)
    assert abs(rmsd - 0.1) < 1e-6
