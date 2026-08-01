"""Regression: the dataset_runner PDB parser must capture the entropy the
engine actually writes.

The C++ writer emits ``REMARK entropy = <value>`` (lowercase, space-equals) at
three sites (LIB/BindingMode.cpp:714,856 and LIB/cluster.cpp:548). The parser's
explicit ``elif "ENTROPY:"`` token (uppercase-colon) is never emitted, and the
fallback regex block rescued RMSD and CF but had no entropy pattern -- so
``entropy_correction`` silently kept its 0.0 initialiser on every pose. That is
a parse failure that looks like a measurement (see issue #350).

These tests pin the writer's spelling to the parser's vocabulary so the
mismatch cannot silently return.
"""

from pathlib import Path

from flexaidds.dataset_runner.runner import DatasetRunner


def _write_pose(work_dir: Path, name: str, body: str) -> None:
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / name).write_text(body)


def test_entropy_remark_written_spelling_is_captured(tmp_path: Path) -> None:
    """``REMARK entropy = <v>`` (the spelling the engine emits) is parsed."""
    _write_pose(
        tmp_path,
        "flexaid_0.pdb",
        "REMARK CF=-5.07911\n"
        "REMARK enthalpy = -5.079114\n"
        "REMARK entropy = 0.00000780\n"
        "REMARK 3.15582 RMSD to ref. structure\n",
    )

    poses = DatasetRunner._parse_flexaid_output(
        tmp_path, target_id="1G9V", ligand_id="lig", structural_state="holo"
    )

    assert len(poses) == 1
    pose = poses[0]
    # The defect: this used to be 0.0 because nothing read `entropy = `.
    assert pose.entropy_correction == 0.00000780
    # enthalpy_score stays sourced from CF= (unchanged), not the enthalpy line.
    assert pose.enthalpy_score == -5.07911
    # total_score = enthalpy_score - entropy_correction, now with a real ΔS.
    assert pose.total_score == -5.07911 - 0.00000780


def test_entropy_defaults_to_zero_when_writer_omits_it(tmp_path: Path) -> None:
    """No entropy REMARK -> entropy_correction keeps its 0.0 default (no crash)."""
    _write_pose(
        tmp_path,
        "flexaid_0.pdb",
        "REMARK CF=-5.07911\n"
        "REMARK 3.15582 RMSD to ref. structure\n",
    )

    poses = DatasetRunner._parse_flexaid_output(
        tmp_path, target_id="1G9V", ligand_id="lig", structural_state="holo"
    )

    assert len(poses) == 1
    assert poses[0].entropy_correction == 0.0


def test_legacy_entropy_token_still_parses(tmp_path: Path) -> None:
    """The explicit ``ENTROPY:`` token remains honoured for back-compat."""
    _write_pose(
        tmp_path,
        "flexaid_0.pdb",
        "REMARK CF_SCORE: -5.07911\n"
        "REMARK ENTROPY: 0.0123\n",
    )

    poses = DatasetRunner._parse_flexaid_output(
        tmp_path, target_id="1G9V", ligand_id="lig", structural_state="holo"
    )

    assert len(poses) == 1
    assert poses[0].entropy_correction == 0.0123
