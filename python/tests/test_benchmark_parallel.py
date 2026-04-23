"""Tests for parallel benchmark orchestration and checkpoint support."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from flexaidds.benchmark import (
    BenchmarkSystem,
    BenchmarkResult,
    MethodResult,
    SystemBenchmarkResult,
    auto_workers,
    run_benchmark,
    _save_checkpoint,
    _load_checkpoint,
    _run_single_system,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_system(sid: str) -> BenchmarkSystem:
    """Create a minimal BenchmarkSystem for testing."""
    return BenchmarkSystem(
        system_id=sid,
        protein_pdb_path=Path("/tmp/receptor.pdb"),
        protein_sequence="ACGT",
        ligand_mol2_path=Path("/tmp/ligand.mol2"),
        ligand_smiles="C",
        reference_pose_pdb_path=Path("/tmp/ref.pdb"),
    )


@pytest.fixture
def three_systems():
    return [_make_system(f"SYS_{i}") for i in range(3)]


# ---------------------------------------------------------------------------
# auto_workers
# ---------------------------------------------------------------------------


def test_auto_workers_returns_positive():
    n = auto_workers()
    assert isinstance(n, int)
    assert 1 <= n <= 16


def test_auto_workers_fallback():
    """When sysconf is unavailable, auto_workers still returns a valid count."""
    with patch("os.sysconf", side_effect=AttributeError):
        n = auto_workers()
        assert 1 <= n <= 16


# ---------------------------------------------------------------------------
# _run_single_system
# ---------------------------------------------------------------------------


def test_run_single_system_skips_on_error(three_systems):
    """_run_single_system returns empty result when docking fails."""
    with patch("flexaidds.benchmark.run_flexaidds", side_effect=RuntimeError("boom")):
        result = _run_single_system(
            three_systems[0],
            methods=("flexaidds",),
            flexaidds_binary=None,
            timeout_per_system=10,
            boltz2_predict_affinity=False,
        )
    assert result.flexaidds_result is None
    assert result.boltz2_result is None
    assert result.system.system_id == "SYS_0"


def test_run_single_system_success(three_systems):
    """_run_single_system captures MethodResult when docking succeeds."""
    fake_result = MethodResult(
        method="flexaidds", system_id="SYS_0",
        best_pose_rmsd_angstrom=1.5, n_poses=10, wall_time_seconds=5.0,
    )
    with patch("flexaidds.benchmark.run_flexaidds", return_value=fake_result):
        result = _run_single_system(
            three_systems[0],
            methods=("flexaidds",),
            flexaidds_binary=None,
            timeout_per_system=10,
            boltz2_predict_affinity=False,
        )
    assert result.flexaidds_result is not None
    assert result.flexaidds_result.best_pose_rmsd_angstrom == 1.5


# ---------------------------------------------------------------------------
# Checkpoint save / load
# ---------------------------------------------------------------------------


def test_checkpoint_round_trip(tmp_path, three_systems):
    """Checkpoint save + load recovers completed systems."""
    ckpt = tmp_path / "checkpoint.json"

    results_by_idx = {
        0: SystemBenchmarkResult(
            system=three_systems[0],
            flexaidds_result=MethodResult(
                method="flexaidds", system_id="SYS_0",
                best_pose_rmsd_angstrom=1.2, wall_time_seconds=3.0,
            ),
        ),
        2: SystemBenchmarkResult(
            system=three_systems[2],
            flexaidds_result=MethodResult(
                method="flexaidds", system_id="SYS_2",
                best_pose_rmsd_angstrom=0.8, wall_time_seconds=2.0,
            ),
        ),
    }

    _save_checkpoint(ckpt, results_by_idx, three_systems)
    assert ckpt.exists()

    loaded = _load_checkpoint(ckpt, three_systems)
    assert 0 in loaded
    assert 2 in loaded
    assert 1 not in loaded
    assert loaded[0].system.system_id == "SYS_0"
    assert loaded[2].system.system_id == "SYS_2"


def test_checkpoint_missing_file(tmp_path, three_systems):
    """_load_checkpoint returns empty dict for nonexistent file."""
    result = _load_checkpoint(tmp_path / "nope.json", three_systems)
    assert result == {}


def test_checkpoint_corrupt_file(tmp_path, three_systems):
    """_load_checkpoint returns empty dict for corrupt JSON."""
    bad = tmp_path / "bad.json"
    bad.write_text("{invalid json")
    result = _load_checkpoint(bad, three_systems)
    assert result == {}


# ---------------------------------------------------------------------------
# run_benchmark — parallel mode
# ---------------------------------------------------------------------------


def test_parallel_benchmark_with_mock(three_systems):
    """run_benchmark with max_workers=2 processes systems in parallel."""
    fake_result = MethodResult(
        method="flexaidds", system_id="",
        best_pose_rmsd_angstrom=1.0, wall_time_seconds=0.1,
    )

    def mock_flexaidds(system, *, binary=None, timeout=3600, temperature=300.0):
        return MethodResult(
            method="flexaidds", system_id=system.system_id,
            best_pose_rmsd_angstrom=1.0, wall_time_seconds=0.01,
        )

    with patch("flexaidds.benchmark.run_flexaidds", side_effect=mock_flexaidds):
        result = run_benchmark(
            three_systems,
            methods=("flexaidds",),
            max_workers=2,
        )

    assert result.n_systems == 3
    ids = [sr.system.system_id for sr in result.systems]
    assert set(ids) == {"SYS_0", "SYS_1", "SYS_2"}
    for sr in result.systems:
        assert sr.flexaidds_result is not None


def test_parallel_error_skip(three_systems):
    """on_error='skip' leaves result as None for failing systems in parallel."""
    call_count = 0

    def sometimes_fail(system, *, binary=None, timeout=3600, temperature=300.0):
        nonlocal call_count
        call_count += 1
        if system.system_id == "SYS_1":
            raise RuntimeError("intentional failure")
        return MethodResult(
            method="flexaidds", system_id=system.system_id,
            wall_time_seconds=0.01,
        )

    with patch("flexaidds.benchmark.run_flexaidds", side_effect=sometimes_fail):
        result = run_benchmark(
            three_systems,
            methods=("flexaidds",),
            max_workers=2,
            on_error="skip",
        )

    assert result.n_systems == 3
    sys1 = next(sr for sr in result.systems if sr.system.system_id == "SYS_1")
    assert sys1.flexaidds_result is None


def test_progress_callback_parallel(three_systems):
    """Progress callback fires once per system in parallel mode."""
    calls = []

    def cb(system_id, completed, total):
        calls.append((system_id, completed, total))

    def mock_fa(system, *, binary=None, timeout=3600, temperature=300.0):
        return MethodResult(
            method="flexaidds", system_id=system.system_id,
            wall_time_seconds=0.01,
        )

    with patch("flexaidds.benchmark.run_flexaidds", side_effect=mock_fa):
        run_benchmark(
            three_systems,
            methods=("flexaidds",),
            max_workers=2,
            progress_callback=cb,
        )

    assert len(calls) == 3
    assert all(t == 3 for _, _, t in calls)


# ---------------------------------------------------------------------------
# run_benchmark — sequential with checkpoint resume
# ---------------------------------------------------------------------------


def test_sequential_checkpoint_resume(tmp_path, three_systems):
    """Sequential run with checkpoint resumes from previously saved state."""
    ckpt = tmp_path / "ckpt.json"

    # Pre-populate checkpoint with SYS_0 done
    pre = BenchmarkResult(systems=(
        SystemBenchmarkResult(
            system=three_systems[0],
            flexaidds_result=MethodResult(
                method="flexaidds", system_id="SYS_0",
                best_pose_rmsd_angstrom=0.5, wall_time_seconds=1.0,
            ),
        ),
    ))
    pre.to_json(ckpt)

    call_ids = []

    def mock_fa(system, *, binary=None, timeout=3600, temperature=300.0):
        call_ids.append(system.system_id)
        return MethodResult(
            method="flexaidds", system_id=system.system_id,
            wall_time_seconds=0.01,
        )

    with patch("flexaidds.benchmark.run_flexaidds", side_effect=mock_fa):
        result = run_benchmark(
            three_systems,
            methods=("flexaidds",),
            max_workers=None,
            checkpoint_path=ckpt,
        )

    # SYS_0 should NOT have been re-run
    assert "SYS_0" not in call_ids
    assert "SYS_1" in call_ids
    assert "SYS_2" in call_ids
    assert result.n_systems == 3


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------


def test_sequential_default(three_systems):
    """run_benchmark without max_workers behaves as before (sequential)."""
    order = []

    def mock_fa(system, *, binary=None, timeout=3600, temperature=300.0):
        order.append(system.system_id)
        return MethodResult(
            method="flexaidds", system_id=system.system_id,
            wall_time_seconds=0.01,
        )

    with patch("flexaidds.benchmark.run_flexaidds", side_effect=mock_fa):
        result = run_benchmark(
            three_systems,
            methods=("flexaidds",),
        )

    # Sequential means deterministic order
    assert order == ["SYS_0", "SYS_1", "SYS_2"]
    assert result.n_systems == 3
