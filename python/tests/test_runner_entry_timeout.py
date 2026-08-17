"""Per-entry FlexAID timeout is configurable (WO-1).

``runner.py`` used to hardcode ``subprocess.run(..., timeout=3600)``. One slow
Astex target exceeded that cap, the liveness gate fired (a timeout counts as
"the engine produced no result"), and every PR that drew that target exited 3 —
including docs-only and .gitignore-only PRs. The cap is now resolved as
``--entry-timeout-seconds`` > ``FLEXAIDDS_ENTRY_TIMEOUT_SECONDS`` > 3600, with
the DEFAULT UNCHANGED so nothing moves for existing runs.
"""

from pathlib import Path

import pytest

from flexaidds.dataset_runner.runner import DatasetRunner


def _runner(tmp_path: Path, binary: str, **kw) -> DatasetRunner:
    return DatasetRunner(results_dir=tmp_path / "results", binary=binary, **kw)


def _slow_binary(tmp_path: Path) -> str:
    binary = tmp_path / "flexaid_slow"
    binary.write_text("#!/bin/sh\nsleep 30\n")
    binary.chmod(0o755)
    return str(binary)


def test_default_cap_is_3600(tmp_path, monkeypatch):
    """With nothing set, the cap is exactly the old hardcoded value."""
    monkeypatch.delenv("FLEXAIDDS_ENTRY_TIMEOUT_SECONDS", raising=False)
    runner = _runner(tmp_path, binary="FlexAID")
    assert runner.entry_timeout_seconds == 3600


def test_env_override_sets_cap(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_ENTRY_TIMEOUT_SECONDS", "5400")
    runner = _runner(tmp_path, binary="FlexAID")
    assert runner.entry_timeout_seconds == 5400


def test_explicit_arg_beats_env(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_ENTRY_TIMEOUT_SECONDS", "5400")
    runner = _runner(tmp_path, binary="FlexAID", entry_timeout_seconds=120)
    assert runner.entry_timeout_seconds == 120


def test_non_positive_cap_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="entry_timeout_seconds"):
        _runner(tmp_path, binary="FlexAID", entry_timeout_seconds=0)


def test_slow_entry_times_out_and_names_the_cap(tmp_path):
    """WO-1 acceptance: with a small cap, a deliberately slow entry times out,
    counts as a liveness crash, and the record carries the cap that killed it
    (so the gate's reason can say TIMED OUT at N s — WO-4)."""
    runner = _runner(
        tmp_path, binary=_slow_binary(tmp_path), entry_timeout_seconds=2
    )
    receptor = tmp_path / "rec.pdb"
    receptor.write_text("")
    ligand = tmp_path / "lig.mol2"
    ligand.write_text("")

    poses = runner._run_flexaid("T1", receptor, [ligand], "holo", with_entropy=False)

    assert poses == []
    assert runner._flexaid_crashes == 1
    assert runner._entry_exit_codes["T1/lig"] is None
    assert runner._entry_timeouts["T1/lig"] == 2


def test_fast_entry_under_small_cap_is_clean(tmp_path):
    """The cap must only kill entries that exceed it — a fast entry stays
    liveness-clean even with a tiny cap."""
    import shutil

    true_bin = shutil.which("true")
    if not true_bin:
        pytest.skip("no `true` binary available to stand in for a fast exec")
    runner = _runner(tmp_path, binary=true_bin, entry_timeout_seconds=30)
    receptor = tmp_path / "rec.pdb"
    receptor.write_text("")
    ligand = tmp_path / "lig.mol2"
    ligand.write_text("")

    runner._run_flexaid("T2", receptor, [ligand], "holo", with_entropy=False)

    assert runner._flexaid_crashes == 0
    assert runner._entry_timeouts == {}
