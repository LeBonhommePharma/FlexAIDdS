"""Liveness gate: a binary that never executes must count as a crash.

The #326 liveness gate scores a run INCONCLUSIVE when ``flexaid_crashes > 0``.
Before this fix, the counter incremented only on a non-zero *return code* — but
when the binary is missing or not executable, ``subprocess.run`` raises
``FileNotFoundError`` / ``PermissionError`` before any return code exists, the
exception is swallowed per-entry by ``_process_one_item``, and the entry looks
like "executed, produced 0 poses" (a productivity signal) rather than "the
engine never ran" (a liveness signal).

That matters because productivity cover is not guaranteed: it is switched off
per-dataset by ``FLEXAIDDS_BENCH_ALLOW_EMPTY`` and is vacuous for a dataset that
declares zero baselines. Liveness is the only gate whose subject is the engine
itself, so it must see the case it is named for.
"""

from pathlib import Path

from flexaidds.dataset_runner.runner import DatasetRunner


def _runner(tmp_path: Path, binary: str) -> DatasetRunner:
    return DatasetRunner(results_dir=tmp_path / "results", binary=binary)


def test_missing_binary_is_recorded_as_a_crash(tmp_path):
    # An absolute path that does not exist -> subprocess.run raises
    # FileNotFoundError (an OSError) before producing a return code.
    runner = _runner(tmp_path, binary=str(tmp_path / "does_not_exist_FlexAID"))
    receptor = tmp_path / "rec.pdb"
    receptor.write_text("")
    ligand = tmp_path / "lig.mol2"
    ligand.write_text("")

    poses = runner._run_flexaid(
        "T1", receptor, [ligand], "holo", with_entropy=False
    )

    assert poses == []                        # no synthetic poses invented
    assert runner._flexaid_crashes == 1       # liveness gate will fire
    # -1 sentinel marks "did not execute", distinct from any real exit code.
    assert runner._entry_exit_codes["T1/lig"] == -1


def test_non_executable_binary_is_recorded_as_a_crash(tmp_path):
    # A file that exists but is not executable -> PermissionError (also OSError).
    binary = tmp_path / "FlexAID_not_exec"
    binary.write_text("#!/bin/sh\necho hi\n")
    binary.chmod(0o644)  # readable, NOT executable
    runner = _runner(tmp_path, binary=str(binary))
    receptor = tmp_path / "rec.pdb"
    receptor.write_text("")
    ligand = tmp_path / "lig.mol2"
    ligand.write_text("")

    poses = runner._run_flexaid(
        "T2", receptor, [ligand], "holo", with_entropy=False
    )

    assert poses == []
    assert runner._flexaid_crashes == 1
    assert runner._entry_exit_codes["T2/lig"] == -1


def test_every_missing_ligand_entry_counts(tmp_path):
    # Two ligands, both hitting a missing binary -> two independent crashes,
    # so the counter reflects the true number of dead invocations.
    runner = _runner(tmp_path, binary=str(tmp_path / "nope_FlexAID"))
    receptor = tmp_path / "rec.pdb"
    receptor.write_text("")
    ligs = [tmp_path / "a.mol2", tmp_path / "b.mol2"]
    for lg in ligs:
        lg.write_text("")

    poses = runner._run_flexaid(
        "T3", receptor, ligs, "holo", with_entropy=False
    )

    assert poses == []
    assert runner._flexaid_crashes == 2
