"""FlexAID's ``flexaid_INI.pdb`` starting structure must never be a pose.

On a self-docking benchmark the input placement IS the crystal pose, so this
file scores ~0 A against the reference (measured 0.0320 A on 1gpk).  Globbing it
into the pose list both depresses ``mean_rmsd`` by a constant amount in every
cell and creates a latent false positive for ``docking_power``.
"""

from pathlib import Path

from flexaidds.dataset_runner.runner import _is_initial_pose_file


class TestIsInitialPoseFile:
    def test_flags_the_initial_structure(self):
        assert _is_initial_pose_file(Path("flexaid_INI.pdb"))

    def test_flags_case_and_prefix_variants(self):
        # FlexAID has shipped both capitalisations of the binary's output prefix.
        assert _is_initial_pose_file(Path("FlexAID_INI.pdb"))
        assert _is_initial_pose_file(Path("flexaid_ini.pdb"))

    def test_does_not_flag_real_poses(self):
        for name in ("flexaid_0.pdb", "flexaid_9.pdb", "flexaid_10.pdb"):
            assert not _is_initial_pose_file(Path(name)), name

    def test_does_not_flag_names_merely_starting_with_ini(self):
        # Guard the suffix match: only a trailing ``_INI`` is the start structure.
        assert not _is_initial_pose_file(Path("flexaid_INIT.pdb"))

    def test_matches_on_full_path(self):
        assert _is_initial_pose_file(Path("/tmp/poses/1gpk/flexaid_INI.pdb"))
        assert not _is_initial_pose_file(Path("/tmp/poses/1gpk/flexaid_3.pdb"))
