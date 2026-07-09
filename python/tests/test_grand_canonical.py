"""Pure-Python μVT helpers — mirror of GrandPartitionFunction semantics."""

from __future__ import annotations

import math

import pytest

from flexaidds.grand_canonical import (
    CompetitiveSite,
    set_concentration,
    kB_kcal,
)


def test_empty_xi():
    s = CompetitiveSite(300.0)
    assert s.log_Xi() == pytest.approx(0.0)
    assert s.empty_probability() == pytest.approx(1.0)
    assert s.mean_N() == pytest.approx(0.0)


def test_set_concentration_api():
    s = CompetitiveSite(300.0)
    s.add("fentanyl", log_Z=10.0, c_M=1e-9)
    s.add("naloxone", log_Z=8.0, c_M=1e-9)
    p0 = s.binding_probability("naloxone")
    set_concentration(s, {"naloxone": 1e-4})
    assert s.binding_probability("naloxone") > p0


def test_set_concentration_sequence():
    s = CompetitiveSite(300.0)
    s.add("A", 5.0, 1.0)
    s.add("B", 5.0, 1.0)
    names = ["A", "B"]
    set_concentration(s, [1e-6, 1e-3], names=names)
    assert s.ligands["A"].concentration_M == pytest.approx(1e-6)
    assert s.ligands["B"].concentration_M == pytest.approx(1e-3)


def test_mor_naloxone_toy_nvt_vs_muvt():
    """Synthetic MOR competitive binding — NVT F independent of c; μVT p depends on c."""
    T = 310.0
    # Stand-in log_Z from "NVT" docking
    log_Z_fen = 15.0
    log_Z_nal = 12.0

    # NVT: stronger fentanyl (more negative F)
    F_fen = -kB_kcal * T * log_Z_fen
    F_nal = -kB_kcal * T * log_Z_nal
    assert F_fen < F_nal

    site = CompetitiveSite(T)
    site.add("fentanyl", log_Z_fen, 1e-9)
    site.add("naloxone", log_Z_nal, 1e-6)
    assert site.binding_probability("fentanyl") + site.binding_probability(
        "naloxone"
    ) + site.empty_probability() == pytest.approx(1.0)

    p_n_lo = site.binding_probability("naloxone")
    set_concentration(site, {"naloxone": 1e-3})
    assert site.binding_probability("naloxone") > p_n_lo
    assert 0.0 <= site.ligand_entropy_collapse() <= 1.0


def test_mixing_entropy_half_sat():
    s = CompetitiveSite(300.0)
    s.add("A", 0.0, 1.0)  # half-sat
    assert s.mixing_entropy() == pytest.approx(kB_kcal * math.log(2.0))


def test_occupancy_curve_monotonic():
    s = CompetitiveSite(300.0)
    s.add("A", 5.0, 1e-9)
    curve = s.occupancy_vs_concentration("A", [1e-9, 1e-6, 1e-3, 1.0])
    assert len(curve) == 4
    for i in range(1, len(curve)):
        assert curve[i].p_species >= curve[i - 1].p_species
