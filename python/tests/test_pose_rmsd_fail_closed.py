"""Pose RMSD must fail closed on atom-count mismatch (no prefix truncation).

Science audit finding 1.1: when docked pose and crystal reference heavy-atom
counts differ, silently computing RMSD on ``pred[:n]`` / ``ref[:n]`` produces
meaningless rates. Sentinel ``-1.0`` is the fail-closed contract so docking
power cannot promote a truncated RMSD into a success.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from flexaidds.dataset_runner.runner import _pose_rmsd_vs_reference


def _write_hetatm_pdb(path: Path, coords: list[tuple[float, float, float]]) -> None:
    """Write HETATM lines with fixed PDB columns (compatible with extract_ligand_coords)."""
    lines = []
    for i, (x, y, z) in enumerate(coords, start=1):
        # Columns: 1-6 record, 7-11 serial, 13-16 name, 18-20 res, 22 chain,
        # 23-26 resseq, 31-38 x, 39-46 y, 47-54 z, 77-78 element
        name = f"C{i}"[:4]
        lines.append(
            f"HETATM{i:5d} {name:>4s} LIG A{1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
        )
    lines.append("END\n")
    path.write_text("".join(lines))


def test_matching_atom_counts_return_finite_rmsd(tmp_path: Path) -> None:
    pose = tmp_path / "pose.pdb"
    _write_hetatm_pdb(pose, [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    ref = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64
    )
    rmsd = _pose_rmsd_vs_reference(pose, ref)
    assert rmsd >= 0.0
    assert rmsd < 1e-6


def test_atom_count_mismatch_returns_sentinel_not_truncated_rmsd(
    tmp_path: Path,
) -> None:
    # Pose has 4 heavy atoms; reference has 3 — old code would truncate to 3
    # and return a finite (wrong) RMSD.
    pose = tmp_path / "pose_extra.pdb"
    _write_hetatm_pdb(
        pose,
        [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 5.0),  # extra atom
        ],
    )
    ref = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64
    )
    rmsd = _pose_rmsd_vs_reference(pose, ref)
    assert rmsd == -1.0


def test_ref_longer_than_pose_also_fails_closed(tmp_path: Path) -> None:
    pose = tmp_path / "pose_short.pdb"
    _write_hetatm_pdb(pose, [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    ref = np.zeros((5, 3), dtype=np.float64)
    assert _pose_rmsd_vs_reference(pose, ref) == -1.0
