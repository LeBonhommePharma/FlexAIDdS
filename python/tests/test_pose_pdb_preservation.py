"""Pose PDBs must survive the docking loop, and land where CI can reach them.

The engine writes ``flexaid_*.pdb`` into a ``tempfile.TemporaryDirectory``
created per ligand (``runner.py`` ``_run_flexaid``).  That directory is removed
when the ligand's iteration ends -- before the target finishes, let alone the
job -- so every coordinate the engine produces is destroyed on the same pass
that extracts its scalars.

Two properties are pinned here, and the second is the one that is easy to lose:

1. With ``FLEXAIDDS_KEEP_POSES`` set, the PDBs exist after the temp dir is gone.
2. They are under ``results_dir`` -- the directory ``benchmark-tier1.yml``
   actually uploads.  A copy that persists in a sibling directory passes any
   "did the files survive" check and still uploads nothing.

Pure Python: no compiled bindings and no FlexAID binary required.
"""

import tempfile
from pathlib import Path

import pytest

from flexaidds.dataset_runner.runner import DatasetRunner


PDB_BODY = (
    "REMARK CF=-5.08119\n"
    "REMARK enthalpy = -5.079114\n"
    "REMARK entropy = 0.00000780\n"
    "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00\n"
    "END\n"
)


def _runner(tmp_path: Path, keep: bool, monkeypatch) -> DatasetRunner:
    if keep:
        monkeypatch.setenv("FLEXAIDDS_KEEP_POSES", "1")
    else:
        monkeypatch.delenv("FLEXAIDDS_KEEP_POSES", raising=False)
    return DatasetRunner(results_dir=tmp_path / "results", dry_run=True)


def _engine_output(work_dir: Path, names=("flexaid_0.pdb", "flexaid_1.pdb")) -> None:
    for name in names:
        (work_dir / name).write_text(PDB_BODY)


def test_poses_survive_the_temp_dir(tmp_path, monkeypatch):
    """The whole point: coordinates outlive the directory they were written in."""
    runner = _runner(tmp_path, keep=True, monkeypatch=monkeypatch)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        _engine_output(work)
        runner._preserve_pose_pdbs(work, "1gpk", "1gpk_ligand", "holo")
    # temp dir is now gone -- exactly the state the real loop leaves behind
    assert not work.exists()

    dest = runner.pose_pdb_dir("1gpk", "1gpk_ligand", "holo")
    assert sorted(p.name for p in dest.glob("*.pdb")) == [
        "flexaid_0.pdb", "flexaid_1.pdb"
    ]
    assert dest.joinpath("flexaid_0.pdb").read_text() == PDB_BODY


def test_destination_is_inside_results_dir(tmp_path, monkeypatch):
    """Reachability, not just persistence.

    ``benchmark-tier1.yml`` uploads ``${{ env.RESULTS_DIR }}/`` and nothing
    else.  If this relation breaks, the files still exist and CI still gets an
    empty artifact -- which looks like the fix worked.
    """
    runner = _runner(tmp_path, keep=True, monkeypatch=monkeypatch)
    dest = runner.pose_pdb_dir("1gpk", "1gpk_ligand", "holo")
    assert runner.results_dir in dest.parents


def test_disabled_by_default(tmp_path, monkeypatch):
    """Default OFF: an existing run's behaviour and disk usage are unchanged."""
    runner = _runner(tmp_path, keep=False, monkeypatch=monkeypatch)
    assert runner.keep_poses is False


def test_flag_is_read_from_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_KEEP_POSES", "1")
    assert DatasetRunner(results_dir=tmp_path / "r", dry_run=True).keep_poses is True
    # "0" is an explicit disable, not merely a set variable
    monkeypatch.setenv("FLEXAIDDS_KEEP_POSES", "0")
    assert DatasetRunner(results_dir=tmp_path / "r", dry_run=True).keep_poses is False


def test_no_pdbs_is_not_an_error(tmp_path, monkeypatch):
    """A crashed entry produces no output; that must not fail the run."""
    runner = _runner(tmp_path, keep=True, monkeypatch=monkeypatch)
    with tempfile.TemporaryDirectory() as tmp:
        assert runner._preserve_pose_pdbs(Path(tmp), "t", "l", "holo") == []
    assert not runner.pose_pdb_dir("t", "l", "holo").exists()


def test_entries_do_not_collide(tmp_path, monkeypatch):
    """Same ligand filename under two targets must not overwrite each other."""
    runner = _runner(tmp_path, keep=True, monkeypatch=monkeypatch)
    for target in ("1gpk", "1hnn"):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            _engine_output(work, names=("flexaid_0.pdb",))
            runner._preserve_pose_pdbs(work, target, "ligand", "holo")
    assert runner.pose_pdb_dir("1gpk", "ligand", "holo").joinpath("flexaid_0.pdb").exists()
    assert runner.pose_pdb_dir("1hnn", "ligand", "holo").joinpath("flexaid_0.pdb").exists()


def test_structural_state_separates_destinations(tmp_path, monkeypatch):
    """holo and apo runs of one target/ligand pair are different measurements."""
    runner = _runner(tmp_path, keep=True, monkeypatch=monkeypatch)
    holo = runner.pose_pdb_dir("1gpk", "ligand", "holo")
    apo = runner.pose_pdb_dir("1gpk", "ligand", "apo")
    assert holo != apo


@pytest.mark.parametrize("hostile", ["../escape", "a/b", "with space"])
def test_ids_cannot_escape_results_dir(tmp_path, monkeypatch, hostile):
    """Target ids reach this from dataset YAML; a path separator must not escape."""
    runner = _runner(tmp_path, keep=True, monkeypatch=monkeypatch)
    dest = runner.pose_pdb_dir(hostile, hostile, "holo")
    assert runner.results_dir in dest.parents


def test_docking_loop_actually_preserves(tmp_path, monkeypatch):
    """The wiring, not just the helper.

    Every other test here calls ``_preserve_pose_pdbs`` directly, so deleting
    the call site inside ``_run_flexaid`` leaves them all green -- a test you
    know how to pass.  This one drives the real loop with a stubbed engine and
    fails if the call site is removed, which is the only version of the check
    that could have caught the defect in the first place.
    """
    from flexaidds.dataset_runner import runner as runner_mod

    class _Result:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, **kwargs):
        # The engine writes into its cwd -- the per-ligand temp dir.
        Path(kwargs["cwd"]).joinpath("flexaid_0.pdb").write_text(PDB_BODY)
        return _Result()

    monkeypatch.setattr(runner_mod.subprocess, "run", fake_run)
    monkeypatch.setenv("FLEXAIDDS_KEEP_POSES", "1")

    receptor = tmp_path / "rec.pdb"
    receptor.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
    ligand = tmp_path / "lig.mol2"
    ligand.write_text(
        "@<TRIPOS>MOLECULE\n"
        "lig\n"
        " 1 0 0 0 0\n"
        "SMALL\nNO_CHARGES\n"
        "@<TRIPOS>ATOM\n"
        "      1 C1          0.0000    0.0000    0.0000 C.3     1  LIG     0.0000\n"
    )

    runner = DatasetRunner(results_dir=tmp_path / "results", dry_run=False)
    runner._run_flexaid("1gpk", receptor, [ligand], "holo", True)

    dest = runner.pose_pdb_dir("1gpk", "lig", "holo")
    assert dest.joinpath("flexaid_0.pdb").exists(), (
        "the docking loop did not preserve the pose PDB -- call site missing?"
    )
