"""Skeleton tests for Grand Canonical Partition Function (Ξ) Python exposure (P2).

Pure-Python path only for now (no C++ required).
Covers:
- _PyGrandPartitionFunction construction, add_ligand (float + StatMechEngine), errors
- log_Xi, p(empty), p(bound), mean_occupancy, variance
- selectivity (apparent vs intrinsic)
- rank()
- LigandSpec validation
- compute_grand_partition helper
- DockingResult grand field roundtrips (to/from json)
- load_results + grand.json sidecar (mock)
- Numerical parity with hand calculations (from test_grand_partition.cpp cases)

Future: when C++ wired, add @requires_core tests + exact cross roundtrip asserts.
See GPF_IMPLEMENTATION_PLAN.md for full validation gates.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from flexaidds import (
    GrandPartitionFunction,
    _PyGrandPartitionFunction,
    LigandSpec,
    LigandRank,
    compute_grand_partition,
    HAS_GRAND_BINDINGS,
    HAS_CORE_BINDINGS,
    StatMechEngine,
    load_results,
    DockingResult,
)
from flexaidds.models import LigandSpec as LigandSpecFromModels


def test_has_grand_bindings_flag_is_false_until_cpp_wired():
    """Initially False; becomes True only after explicit _core grand exposure."""
    assert HAS_GRAND_BINDINGS is False
    # Pure fallback always importable
    assert _PyGrandPartitionFunction is not None
    assert GrandPartitionFunction is _PyGrandPartitionFunction  # alias for now


def test_empty_xi():
    gpf = GrandPartitionFunction(300.0)
    assert gpf.num_ligands() == 0
    assert abs(gpf.log_Xi() - 0.0) < 1e-12
    assert abs(gpf.empty_probability() - 1.0) < 1e-12
    assert abs(gpf.mean_occupancy() - 0.0) < 1e-12


def test_invalid_temperature():
    with pytest.raises(ValueError):
        GrandPartitionFunction(0.0)
    with pytest.raises(ValueError):
        GrandPartitionFunction(-10.0)


def test_single_ligand():
    gpf = GrandPartitionFunction(300.0)
    log_Z_A = 10.0
    gpf.add_ligand("ligandA", log_Z_A)

    assert gpf.num_ligands() == 1
    assert gpf.has_ligand("ligandA")

    expected_log_xi = math.log(1.0 + math.exp(10.0))
    assert abs(gpf.log_Xi() - expected_log_xi) < 1e-10

    expected_pA = math.exp(10.0) / (1.0 + math.exp(10.0))
    assert abs(gpf.binding_probability("ligandA") - expected_pA) < 1e-8
    assert abs(gpf.empty_probability() - (1.0 - expected_pA)) < 1e-8

    # F_bound = -kT ln Z
    kT = 0.001987206 * 300.0
    assert abs(gpf.F_bound("ligandA") - (-kT * log_Z_A)) < 1e-10


def test_two_ligands_competitive():
    gpf = GrandPartitionFunction(300.0)
    gpf.add_ligand("A", 10.0)
    gpf.add_ligand("B", 8.0)

    xi = 1.0 + math.exp(10.0) + math.exp(8.0)
    assert abs(math.exp(gpf.log_Xi()) - xi) < xi * 1e-10

    # selectivity A/B (apparent, here conc=1M so = intrinsic)
    assert abs(gpf.selectivity("A", "B") - math.exp(2.0)) < 1e-10
    assert abs(gpf.log_selectivity("A", "B") - 2.0) < 1e-12
    assert abs(gpf.log_intrinsic_selectivity("A", "B") - 2.0) < 1e-12

    pA = gpf.binding_probability("A")
    pB = gpf.binding_probability("B")
    pE = gpf.empty_probability()
    assert abs(pA + pB + pE - 1.0) < 1e-10

    ranks = gpf.rank()
    assert len(ranks) == 2
    assert isinstance(ranks[0], LigandRank)
    assert ranks[0].name == "A"  # better dG
    assert ranks[0].p_bound == pytest.approx(pA, rel=1e-9)


def test_equal_ligands_probs():
    gpf = GrandPartitionFunction(300.0)
    gpf.add_ligand("X", 5.0)
    gpf.add_ligand("Y", 5.0)
    gpf.add_ligand("Z", 5.0)

    pX = gpf.binding_probability("X")
    pY = gpf.binding_probability("Y")
    pZ = gpf.binding_probability("Z")
    assert abs(pX - pY) < 1e-12
    assert abs(pY - pZ) < 1e-12
    assert abs(gpf.empty_probability() + 3 * pX - 1.0) < 1e-10


def test_concentration_affects_apparent_not_intrinsic():
    gpf = GrandPartitionFunction(300.0)
    # Z_A = e^10 , Z_B = e^8
    # c_A=1uM=1e-6 M , c_B=1M
    gpf.add_ligand("A", 10.0, concentration_M=1e-6)
    gpf.add_ligand("B", 8.0, concentration_M=1.0)

    log_intr = gpf.log_intrinsic_selectivity("A", "B")
    assert abs(log_intr - 2.0) < 1e-12   # pure Z ratio

    log_app = gpf.log_selectivity("A", "B")
    # log_app = (log_cA + 10) - (log_cB + 8) = log(1e-6) + 2 ≈ -13.8155 + 2
    expected_log_app = math.log(1e-6) + 2.0
    assert abs(log_app - expected_log_app) < 1e-10

    # apparent selectivity << 1 because low conc A
    sel = gpf.selectivity("A", "B")
    assert sel < 1.0
    assert sel > 0.0


def test_add_ligand_from_statmech_engine():
    engine = StatMechEngine(300.0)
    for e in [-12.0, -11.5, -13.0]:
        engine.add_sample(e)
    thermo = engine.compute()
    assert thermo.log_Z is not None

    gpf = GrandPartitionFunction(300.0)
    gpf.add_ligand("from_engine", engine, concentration_M=0.001)
    assert gpf.has_ligand("from_engine")
    assert abs(gpf.binding_probability("from_engine")) > 0.0


def test_ligand_spec_and_compute_helper():
    specA = LigandSpec(name="A", concentration_M=1.0)
    specB = LigandSpec(name="B", concentration_M=1e-3, ligand_id="B001")

    assert isinstance(specA, LigandSpecFromModels)  # same class

    logz = {"A": 12.0, "B": 9.5}
    gpf = compute_grand_partition(
        [specA, specB],
        temperature_K=298.0,
        log_Z_map=logz,
    )
    assert gpf.temperature == 298.0
    assert gpf.num_ligands() == 2
    p_sum = sum(gpf.probabilities().values())
    assert abs(p_sum - 1.0) < 1e-10

    # error if no logz
    with pytest.raises(ValueError):
        compute_grand_partition([specA])


def test_ligand_spec_validation():
    with pytest.raises(ValueError):
        LigandSpec("bad", concentration_M=0.0)
    with pytest.raises(ValueError):
        LigandSpec("bad2", concentration_M=2000.0)


def test_grand_docking_result_fields_and_json_roundtrip(tmp_path: Path):
    # Create a minimal DockingResult with grand populated
    res = DockingResult(
        source_dir=tmp_path,
        binding_modes=[],
        grand_log_xi=math.log(1 + math.exp(10.0)),
        ligand_occupancies={"A": 0.99995, "__empty__": 5e-5},
        selectivities={"A/B": 7.389},
        empty_probability=5e-5,
        mean_occupancy=0.99995,
    )
    assert res.grand_log_xi is not None
    assert "A" in res.ligand_occupancies

    j = res.to_json()
    assert "grand_log_xi" in j
    loaded = DockingResult.from_json(j)
    assert abs(loaded.grand_log_xi - res.grand_log_xi) < 1e-9
    assert loaded.ligand_occupancies.get("A") == pytest.approx(0.99995, rel=1e-5)


def test_load_results_with_grand_sidecar(tmp_path: Path):
    # Write a dummy pose PDB
    pdb = tmp_path / "mode_1_pose_1.pdb"
    pdb.write_text(
        "REMARK binding_mode = 1\nREMARK pose_rank = 1\nREMARK CF = -10.0\n"
        "ATOM      1  C   LIG A   1       0   0   0  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    # Write grand sidecar
    side = {
        "grand_log_xi": 10.2,
        "ligand_occupancies": {"L1": 0.8, "__empty__": 0.2},
        "empty_probability": 0.2,
    }
    (tmp_path / "grand.json").write_text(json.dumps(side), encoding="utf-8")

    loaded = load_results(tmp_path)
    assert loaded.grand_log_xi == pytest.approx(10.2)
    assert loaded.ligand_occupancies.get("L1") == 0.8
    assert loaded.metadata.get("grand_sidecar") is True


def test_grand_remarks_parsed_into_metadata(tmp_path: Path):
    """REMARKs with grand keys are captured (via io aliases + remarks dict)."""
    pdb = tmp_path / "bm_grand.pdb"
    pdb.write_text(
        "REMARK grand_log_xi = 5.123\nREMARK p_empty = 0.01\n"
        "ATOM      1  C   LIG A   1       0   0   0  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    res = load_results(tmp_path)
    # The top level grand not auto-lifted from per-pose remarks here, but
    # remarks on poses contain the keys (future lifting possible).
    pose0 = res.binding_modes[0].poses[0]
    assert "grand_log_xi" in pose0.remarks or "grand_log_xi" in str(pose0.remarks)
    # Non-breaking: legacy path still works
    assert res.n_modes == 1


# Placeholder for C++/Py roundtrip gate (to be filled when bindings wired)
def test_grand_pure_py_roundtrip_plan():
    """Documented plan: once C++ GrandPartitionFunction bound,
    instantiate both, feed same log_Z + concs, assert all quantities
    (log_Xi, p, selectivities, rank dG) match within 1e-9 rel.
    Also test load_results on a multi-ligand GPF-augmented output dir.
    """
    # This test always passes; it serves as executable reminder.
    assert not HAS_GRAND_BINDINGS or HAS_CORE_BINDINGS  # either way ok for skeleton
    # TODO(P3/P4): add exact match tests + fixtures with known analytical solutions.
