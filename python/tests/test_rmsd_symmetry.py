"""Symmetry-corrected RMSD, matching FlexAID's Hungarian branch.

``FlexAID/LIB/calc_rmsd.c::calc_Hungarian_RMSD`` builds one cost matrix per
atom TYPE and solves the assignment problem inside each type.  Without it, a
symmetric ligand is penalised for an equivalent relabelling: a benzene ring
rotated by one vertex is the *same pose* and the *same molecule*, but paired
positionally every atom appears displaced.

That is not a cosmetic difference.  The Astex 2.0 A redocking criterion is
defined on the symmetry-corrected value, so an uncorrected metric fails poses
that are correct and reports a docking-power number that is too low.

The pairing is constrained to run within an element: carbon may only be
reassigned to carbon.  A cross-element swap could only ever lower the RMSD by
matching chemically different atoms, which would flatter the result.
"""

import numpy as np
import pytest

from flexaidds.benchmark import _symmetry_permutation, compute_rmsd


def _square_ring(rotation: int = 0):
    """Four equivalent carbons on a unit square, optionally relabelled."""
    xyz = np.array(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
        dtype=np.float64,
    )
    return np.roll(xyz, rotation, axis=0)


def test_relabelled_symmetric_ring_is_zero_rmsd():
    """The headline case: same atoms, same places, different labels."""
    ref = _square_ring(0)
    pred = _square_ring(1)  # identical set of positions, rotated labelling
    elements = ["C"] * 4

    positional = compute_rmsd(pred, ref)
    corrected = compute_rmsd(pred, ref, elements)

    assert positional > 1.0, "control: positional pairing must see a difference"
    assert corrected == pytest.approx(0.0, abs=1e-9)


def test_a_genuinely_displaced_pose_is_not_rescued():
    """Symmetry correction must not turn a real miss into a hit.

    The whole risk of this change is that it flatters results.  Translating
    every atom by 3 A is a real 3 A error and no relabelling can remove it.
    """
    ref = _square_ring(0)
    pred = ref + np.array([3.0, 0.0, 0.0])
    elements = ["C"] * 4

    assert compute_rmsd(pred, ref, elements) == pytest.approx(3.0, abs=1e-9)


def test_elements_never_swap_across_types():
    """Carbon must not be paired with oxygen even when that would score better.

    Here the O sits exactly on a C position and vice versa, so a type-blind
    assignment would report 0.0.  FlexAID's per-type matrices forbid it, and
    so must this.
    """
    ref = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=np.float64)
    pred = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float64)
    elements = ["C", "O"]

    corrected = compute_rmsd(pred, ref, elements)
    assert corrected == pytest.approx(2.0, abs=1e-9), (
        "a C/O swap was allowed -- the assignment is not type-constrained"
    )


def test_permutation_is_a_permutation():
    """Every atom used exactly once: no atom may be matched twice."""
    ref = _square_ring(0)
    pred = _square_ring(2)
    perm = _symmetry_permutation(pred, ref, ["C"] * 4)

    assert sorted(perm.tolist()) == [0, 1, 2, 3]


def test_correction_never_increases_rmsd():
    """The assignment is a minimisation, so it is bounded by the positional value."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        ref = rng.normal(size=(6, 3))
        pred = rng.normal(size=(6, 3))
        elements = ["C", "C", "C", "N", "N", "O"]
        assert compute_rmsd(pred, ref, elements) <= compute_rmsd(pred, ref) + 1e-12


def test_single_atom_of_a_type_is_untouched():
    """One atom of an element has no partner to swap with."""
    ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    pred = np.array([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0]], dtype=np.float64)
    assert compute_rmsd(pred, ref, ["C", "O"]) == pytest.approx(0.5, abs=1e-9)


def test_element_count_mismatch_raises():
    """A wrong-length element list is a caller bug, not a silent fallback."""
    ref = _square_ring(0)
    with pytest.raises(ValueError, match="elements has"):
        compute_rmsd(ref, ref, ["C", "C"])
