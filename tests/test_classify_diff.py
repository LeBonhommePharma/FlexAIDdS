"""Path buckets: engine-critical vs swarm/audit pack vs hygiene."""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import check_repo_hygiene as hygiene  # noqa: E402


def test_classify_engine_vs_pack():
    globs = ["LIB/gaboom.h", "LIB/RngSeed.h", "LIB/statmech.cpp"]
    buckets = hygiene.classify_paths(
        [
            "LIB/gaboom.h",
            "docs/swarm/2026-08-13/ASTEX84_FROZEN_POSE_BENCHMARK.csv",
            "tests/test_gaboom.cpp",
            "scripts/check_repo_hygiene.py",
        ],
        globs,
    )
    assert buckets["engine_critical"] == ["LIB/gaboom.h"]
    assert buckets["science_pack"] == [
        "docs/swarm/2026-08-13/ASTEX84_FROZEN_POSE_BENCHMARK.csv"
    ]
    assert buckets["tests"] == ["tests/test_gaboom.cpp"]
    assert buckets["ci_hygiene"] == ["scripts/check_repo_hygiene.py"]


def test_pack_bundling_flags_mix(tmp_path: Path, monkeypatch):
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "SCIENCE_CRITICAL_PATHS.txt").write_text(
        "LIB/gaboom.h\n", encoding="utf-8"
    )
    monkeypatch.setattr(
        hygiene,
        "_changed_files",
        lambda base="origin/main": (
            ["LIB/gaboom.h", "docs/swarm/README.md"],
            None,
        ),
    )
    errors = hygiene.check_science_pack_bundling(tmp_path)
    assert errors
    assert "split the PR" in errors[0]


def test_pack_bundling_ok_when_split(tmp_path: Path, monkeypatch):
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "SCIENCE_CRITICAL_PATHS.txt").write_text(
        "LIB/gaboom.h\n", encoding="utf-8"
    )
    monkeypatch.setattr(
        hygiene,
        "_changed_files",
        lambda base="origin/main": (["LIB/gaboom.h", "tests/test_gaboom.cpp"], None),
    )
    assert hygiene.check_science_pack_bundling(tmp_path) == []
