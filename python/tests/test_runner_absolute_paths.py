"""Paths handed to the engine must be absolute at construction -- except the
binary, where "absolute" is the wrong rule.

FlexAID is invoked with ``cwd=<per-entry temp dir>`` (``_run_flexaid``), so a
relative directory means one thing to the runner and another to the child. The
data roots therefore get ``.resolve()``.

The binary is NOT a directory and NOT reliably a path: ``subprocess`` looks a
bare name up on ``PATH`` and ignores ``cwd``, so ``.resolve()`` there would
convert a working lookup into a nonexistent absolute path. The separator is what
says which of the two you have. These tests pin both halves, because the obvious
uniform fix passes the first three and breaks the fourth.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from flexaidds.dataset_runner.runner import DatasetConfig, DatasetRunner


def _runner(tmp_path, **kw):
    return DatasetRunner(results_dir=str(tmp_path / "results"), **kw)


def test_cache_dir_default_is_absolute(tmp_path, monkeypatch):
    """The fallback literal is relative -- the default is the trap."""
    monkeypatch.delenv("FLEXAIDDS_BENCHMARK_DATA", raising=False)
    monkeypatch.chdir(tmp_path)
    assert _runner(tmp_path).cache_dir.is_absolute()


def test_cache_dir_from_env_is_absolute(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_BENCHMARK_DATA", "benchmark_data")
    monkeypatch.chdir(tmp_path)
    r = _runner(tmp_path)
    assert r.cache_dir.is_absolute()
    assert r.cache_dir == (tmp_path / "benchmark_data").resolve()


def test_data_dir_from_yaml_is_absolute(tmp_path, monkeypatch):
    """A relative data_dir: in a config is a claim about every caller."""
    monkeypatch.chdir(tmp_path)
    yaml_path = tmp_path / "toy.yaml"
    yaml_path.write_text(
        "slug: toy\nname: Toy\ndocking_mode: self_docking\n"
        "data_dir: benchmarks/toy\ntargets: [t1]\n"
    )
    cfg = DatasetConfig.from_yaml(yaml_path)
    assert cfg.data_dir is not None
    assert cfg.data_dir.is_absolute()
    assert cfg.data_dir == (tmp_path / "benchmarks" / "toy").resolve()


def test_bare_binary_name_is_NOT_absolutised(tmp_path):
    """The case the obvious fix destroys.

    ``Path("FlexAID").resolve()`` yields ``<cwd>/FlexAID`` -- a path that does
    not exist -- turning a working PATH lookup into a hard FileNotFoundError.
    A bare name must survive untouched unless PATH can actually resolve it.
    """
    r = _runner(tmp_path, binary="definitely-not-on-path-xyzzy")
    assert r.binary == "definitely-not-on-path-xyzzy"
    assert not os.path.isabs(r.binary)


def test_bare_name_on_PATH_still_resolves(tmp_path):
    """A bare name that PATH *can* resolve becomes absolute, via PATH -- not cwd."""
    import shutil

    if shutil.which("sh") is None:  # pragma: no cover - POSIX runners have sh
        pytest.skip("no 'sh' on PATH")
    r = _runner(tmp_path, binary="sh")
    assert os.path.isabs(r.binary)
    assert Path(r.binary).exists()
    assert Path(r.binary).name == "sh"


def test_relative_binary_with_separator_is_absolutised(tmp_path, monkeypatch):
    """Contains a separator -> subprocess would resolve it against the tmp dir."""
    build = tmp_path / "build"
    build.mkdir()
    exe = build / "FlexAID"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(0o755)
    monkeypatch.chdir(tmp_path)

    r = _runner(tmp_path, binary=os.path.join("build", "FlexAID"))
    assert os.path.isabs(r.binary)
    assert Path(r.binary) == exe.resolve()
