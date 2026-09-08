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
import pytest

from flexaidds.dataset_runner.runner import DatasetRunner
from flexaidds.dataset_runner.metrics import docking_power, entropy_rescue_rate


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
        "REMARK temperature = 300.00\n"
        "REMARK 3.15582 RMSD to ref. structure\n",
    )

    poses = DatasetRunner._parse_flexaid_output(
        tmp_path, target_id="1G9V", ligand_id="lig", structural_state="holo"
    )

    assert len(poses) == 1
    pose = poses[0]
    # The defect: this used to be 0.0 because nothing read `entropy = `.
    assert pose.entropy == 0.00000780
    assert pose.temperature == 300.0
    assert pose.entropy_correction == 300.0 * 0.00000780
    # enthalpy_score stays sourced from CF= (unchanged), not the enthalpy line.
    assert pose.enthalpy_score == -5.07911
    # The separate reranking proxy is CF - T*S, using the emitted temperature.
    assert pose.total_score == -5.07911 - 300.0 * 0.00000780
    assert pose.ensemble_mean_energy == -5.079114


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


def test_legacy_entropy_token_requires_emitted_temperature(tmp_path: Path) -> None:
    """An old spelling supplies S, not permission to guess its temperature."""
    _write_pose(
        tmp_path,
        "flexaid_0.pdb",
        "REMARK CF_SCORE: -5.07911\n"
        "REMARK ENTROPY: 0.0123\n",
    )

    with pytest.raises(ValueError, match="without emitted temperature"):
        DatasetRunner._parse_flexaid_output(
            tmp_path, target_id="1G9V", ligand_id="lig", structural_state="holo")


def _parse(path, objective="cf_minus_ts"):
    return DatasetRunner._parse_flexaid_output(
        path, "target", "lig", "holo", ranking_objective=objective)


@pytest.mark.parametrize("temperature", [100.0, 300.0, 600.0])
def test_correction_depends_on_emitted_temperature(tmp_path, temperature):
    _write_pose(tmp_path, "flexaid_0.pdb",
                f"REMARK CF=-10\nREMARK entropy = 0.01\nREMARK temperature = {temperature}\n")
    pose = _parse(tmp_path)[0]
    assert pose.entropy_correction == temperature * 0.01
    assert pose.total_score == -10 - temperature * 0.01


def test_dimensional_counterexample_reverses_ranking_but_not_generator(tmp_path):
    # Arbitrary scalar counterexample, not a molecular accuracy claim.
    _write_pose(tmp_path, "flexaid_0.pdb",
                "REMARK CF=-11\nREMARK entropy = 0.001\nREMARK temperature = 300.00\n"
                "REMARK enthalpy = -11\nREMARK free_energy = -11.3\n"
                "REMARK soft_beta_G = -15\nREMARK binding_mode = 0\nREMARK RMSD: 5\n")
    _write_pose(tmp_path, "flexaid_1.pdb",
                "REMARK CF=-10\nREMARK entropy = 0.010\nREMARK temperature = 300.00\n"
                "REMARK enthalpy = -10\nREMARK free_energy = -13\n"
                "REMARK soft_beta_G = -14\nREMARK binding_mode = 1\nREMARK RMSD: 1\n")
    corrected, legacy = _parse(tmp_path), _parse(tmp_path, "legacy_cf_minus_s")
    assert docking_power(corrected) == 1.0
    assert docking_power(legacy) == 0.0
    assert entropy_rescue_rate(corrected, rank_threshold=1) == 1.0
    assert entropy_rescue_rate(legacy, rank_threshold=1) == 0.0
    assert docking_power(corrected, ranking="generator") == 0.0
    assert docking_power(legacy, ranking="generator") == 0.0
    assert corrected[1].ensemble_free_energy == -13
    assert corrected[1].enthalpy_score == -10
    assert corrected[1].generator_score == -14
    assert corrected[1].pose_rank == 2


def test_zero_entropy_is_score_identical_without_assumed_temperature(tmp_path):
    _write_pose(tmp_path, "flexaid_0.pdb", "REMARK CF=-10\nREMARK entropy = 0.00000000\n")
    corrected, legacy = _parse(tmp_path)[0], _parse(tmp_path, "legacy_cf_minus_s")[0]
    assert corrected.entropy_correction == legacy.entropy_correction == 0.0
    assert corrected.total_score == legacy.total_score == -10.0
    assert corrected.temperature is None


@pytest.mark.parametrize("temperature", ["0", "-1", "nan", "inf"])
def test_invalid_temperature_refuses_entire_election(tmp_path, temperature):
    _write_pose(tmp_path, "flexaid_0.pdb",
                f"REMARK CF=-10\nREMARK entropy = 0.01\nREMARK temperature = {temperature}\n")
    with pytest.raises(ValueError, match="Invalid emitted temperature"):
        _parse(tmp_path)


def test_explicit_zero_score_is_preserved_as_emitted_metadata(tmp_path):
    _write_pose(tmp_path, "flexaid_0.pdb",
                "REMARK CF=-10\nREMARK entropy = 0.01\nREMARK temperature = 300\n"
                "REMARK TOTAL_SCORE: 0\nREMARK free_energy = -80\nREMARK enthalpy = -77\n")
    pose = _parse(tmp_path)[0]
    assert pose.total_score == -13
    assert pose.emitted_total_score == 0
    assert pose.ensemble_mean_energy == -77
    assert pose.ensemble_free_energy == -80


def test_generator_election_uses_numeric_filename_not_cluster_identity(tmp_path):
    # DensityPeak sorts clusters before writing _i.pdb but retains original IDs.
    for output_index, cluster_id, rmsd in [(10, 0, 5), (0, 9, 1), (2, 1, 5)]:
        _write_pose(tmp_path, f"flexaid_{output_index}.pdb",
                    f"REMARK CF=-10\nREMARK binding_mode = {cluster_id}\nREMARK RMSD: {rmsd}\n")
    poses = _parse(tmp_path)
    assert [pose.pose_rank for pose in poses] == [1, 3, 11]
    assert [pose.binding_mode_id for pose in poses] == [9, 1, 0]
    assert docking_power(poses, ranking="generator") == 1.0


def test_missing_generator_top1_is_not_replaced_by_later_output(tmp_path):
    _write_pose(tmp_path, "flexaid_1.pdb",
                "REMARK CF=-10\nREMARK binding_mode = 0\nREMARK RMSD: 1\n")
    poses = _parse(tmp_path)
    assert poses[0].pose_rank == 2
    assert docking_power(poses, ranking="generator") == 0.0
    assert docking_power(poses, ranking="generator", top_n=3) == 1.0
