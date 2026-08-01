"""The ligand is one residue in the pose file, not every non-water HETATM.

FlexAID writes the receptor alongside the docked ligand.  36 of the 85 Astex
Diverse receptors carry a non-water HETATM cofactor -- ions in most cases, a
172-atom HEM in four.  Selecting "HETATM and not HOH" therefore collects the
cofactor as if it were part of the ligand.

On 1mq6 that yields 37 atoms against a 36-atom reference, with a calcium at
index 0 shifting every subsequent atom by one position.  The old caller then
truncated both arrays to the shorter length and computed an RMSD anyway --
producing 6.49 A where the true best-of-10 is 7.18 A.  Optimistic, plausible,
and indistinguishable from a measurement.

Both halves are pinned here: correct selection, and refusal to guess.
"""

import numpy as np
import pytest

from flexaidds.benchmark import extract_ligand_coords_from_pdb


def _het(serial, name, resname, chain, resseq, xyz, element):
    x, y, z = xyz
    return (
        f"HETATM{serial:>5} {name:<4}{resname:>3} {chain}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}\n"
    )


def _pose_with_cofactor(tmp_path, n_ligand=36):
    """A pose file shaped like 1mq6's: one Ca ion, then the ligand."""
    p = tmp_path / "flexaid_0.pdb"
    body = _het(1, "CA", "CA", "A", 1, (0.0, 0.0, 0.0), "Ca")
    for i in range(n_ligand):
        body += _het(2 + i, f"C {i}", "XLD", "A", 2, (float(i), 0.0, 0.0), "C")
    body += "END\n"
    p.write_text(body)
    return p


def test_cofactor_is_not_counted_as_ligand(tmp_path):
    """The regression: 37 atoms in, 36 out, and the Ca is not among them."""
    pose = _pose_with_cofactor(tmp_path)

    # Without a count the extractor used to concatenate every non-water
    # HETATM residue and return 37.  That union WAS the bug, so it now
    # raises instead of handing back a silently contaminated array.
    with pytest.raises(ValueError, match="no expected_n_atoms"):
        extract_ligand_coords_from_pdb(pose)

    coords = extract_ligand_coords_from_pdb(pose, expected_n_atoms=36)
    assert coords.shape == (36, 3)
    # the calcium sat at the origin; the ligand starts at x=0 too, so check
    # the *last* atom instead -- it must be the 36th ligand atom, not the 35th
    assert coords[-1][0] == pytest.approx(35.0)


def test_correspondence_is_not_shifted(tmp_path):
    """The consequence the atom count only hints at.

    With the cofactor included, atom i of the pose is atom i-1 of the ligand.
    Comparing against a reference then measures the wrong pairs.
    """
    pose = _pose_with_cofactor(tmp_path)
    ref = np.array([[float(i), 0.0, 0.0] for i in range(36)])

    good = extract_ligand_coords_from_pdb(pose, expected_n_atoms=36)
    assert np.allclose(good, ref)  # exact correspondence, zero RMSD

    # The shifted correspondence is no longer reachable: the countless path
    # refuses the union rather than returning 37 atoms for a 36-atom ligand.
    with pytest.raises(ValueError, match="Refusing to concatenate"):
        extract_ligand_coords_from_pdb(pose)


def test_refuses_to_guess_when_no_residue_matches(tmp_path):
    """A mismatch is a bug in selection, not a condition to work around."""
    pose = _pose_with_cofactor(tmp_path)
    with pytest.raises(ValueError, match="Cannot identify the ligand"):
        extract_ligand_coords_from_pdb(pose, expected_n_atoms=99)


def test_refuses_to_guess_when_two_residues_match(tmp_path):
    """Ambiguity must fail loudly rather than take the first."""
    p = tmp_path / "two.pdb"
    body = ""
    for r, resname in enumerate(("LIG", "XLD")):
        for i in range(5):
            body += _het(
                1 + r * 5 + i, f"C {i}", resname, "A", r + 1, (float(i), 0.0, 0.0), "C"
            )
    p.write_text(body + "END\n")
    with pytest.raises(ValueError, match="found 2 candidates"):
        extract_ligand_coords_from_pdb(p, expected_n_atoms=5)


def test_a_large_cofactor_does_not_win(tmp_path):
    """Why 'pick the biggest residue' is not the fix.

    1g9v carries a 172-atom HEM against a 25-atom ligand.  Size selects the
    cofactor; the reference atom count selects the ligand.
    """
    p = tmp_path / "hem.pdb"
    body = ""
    for i in range(172):
        body += _het(1 + i, f"C {i}", "HEM", "A", 1, (float(i), 1.0, 0.0), "C")
    for i in range(25):
        body += _het(200 + i, f"C {i}", "LIG", "A", 2, (float(i), 0.0, 0.0), "C")
    p.write_text(body + "END\n")

    coords = extract_ligand_coords_from_pdb(p, expected_n_atoms=25)
    assert coords.shape == (25, 3)
    assert coords[0][1] == pytest.approx(0.0)  # the ligand's y, not HEM's


def test_single_residue_needs_no_count(tmp_path):
    """A ligand-only PDB is unambiguous, so the countless path still works.

    The union was only ever wrong when there was more than one candidate.
    Refusing outright would have broken every ligand-only reference file,
    so one residue is returned as before.
    """
    p = tmp_path / "one.pdb"
    p.write_text(
        "".join(
            f"HETATM{i:5d}  C   LIG A   1    {float(i):8.3f}{0.0:8.3f}{0.0:8.3f}"
            f"  1.00  0.00           C\n"
            for i in range(5)
        )
    )
    coords = extract_ligand_coords_from_pdb(p)
    assert coords.shape == (5, 3)
    assert coords[-1][0] == pytest.approx(4.0)


def test_reference_loader_refuses_contaminated_pdb(tmp_path):
    """The REFERENCE side, which #363 did not cover.

    _reference_ligand_coords is where expected_n_atoms comes FROM, so it
    cannot pass one.  A reference PDB carrying a cofactor therefore used to
    define a contaminated count, which then selected a contaminated-sized
    residue from every pose compared against it.
    """
    from flexaidds.dataset_runner.runner import _reference_ligand_coords

    ref = _pose_with_cofactor(tmp_path)
    with pytest.raises(ValueError, match="Refusing to concatenate"):
        _reference_ligand_coords(ref)
