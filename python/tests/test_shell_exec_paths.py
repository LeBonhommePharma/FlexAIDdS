"""Unit tests for dataset_runner path safety helpers (mirrors LIB/shell_exec.h)."""

from __future__ import annotations

import pytest

from flexaidds.dataset_runner.data_paths import is_safe_exec_path, validate_exec_path


def test_accept_normal_paths():
    validate_exec_path("/tmp/data/MC_st0r5.2_6.dat")
    validate_exec_path("relative/path with spaces.pdb")
    assert is_safe_exec_path("/tmp/ligand.sdf")


def test_accept_tab():
    validate_exec_path("name\twith\ttab")


def test_reject_empty():
    with pytest.raises(ValueError, match="empty"):
        validate_exec_path("")
    assert not is_safe_exec_path("")


def test_reject_nul():
    with pytest.raises(ValueError, match="NUL"):
        validate_exec_path("evil\x00payload")


def test_reject_newline():
    with pytest.raises(ValueError, match="newline"):
        validate_exec_path("a\nb")
    with pytest.raises(ValueError, match="newline"):
        validate_exec_path("a\rb")


def test_reject_control():
    with pytest.raises(ValueError, match="control"):
        validate_exec_path("a\x01b")
