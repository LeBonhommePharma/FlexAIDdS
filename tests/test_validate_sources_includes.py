"""Quoted #include from a compiled TU counts as a source reference."""
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import validate_sources as vs  # noqa: E402


def test_include_from_compiled_tu_is_not_orphan(tmp_path: Path):
    lib = tmp_path / "LIB"
    lib.mkdir()
    (lib / "used.cpp").write_text('#include "EnvFlags.h"\nvoid f() {}\n', encoding="utf-8")
    (lib / "EnvFlags.h").write_text("#pragma once\n", encoding="utf-8")
    (lib / "orphan.h").write_text("#pragma once\n", encoding="utf-8")
    (tmp_path / "CMakeLists.txt").write_text("add_library(x LIB/used.cpp)\n", encoding="utf-8")

    result = vs.validate_sources(root=tmp_path, strict=True)
    orphans = {p.name for p in result.orphans}
    assert "EnvFlags.h" not in orphans
    assert "orphan.h" in orphans


def test_transitive_header_include_is_not_orphan(tmp_path: Path):
    lib = tmp_path / "LIB"
    lib.mkdir()
    (lib / "used.cpp").write_text('#include "Mid.h"\n', encoding="utf-8")
    (lib / "Mid.h").write_text('#pragma once\n#include "Leaf.h"\n', encoding="utf-8")
    (lib / "Leaf.h").write_text("#pragma once\n", encoding="utf-8")
    (tmp_path / "CMakeLists.txt").write_text("add_library(x LIB/used.cpp)\n", encoding="utf-8")

    result = vs.validate_sources(root=tmp_path, strict=True)
    names = {p.name for p in result.orphans}
    assert "Mid.h" not in names
    assert "Leaf.h" not in names


def test_angle_bracket_include_does_not_count(tmp_path: Path):
    lib = tmp_path / "LIB"
    lib.mkdir()
    (lib / "used.cpp").write_text("#include <EnvFlags.h>\n", encoding="utf-8")
    (lib / "EnvFlags.h").write_text("#pragma once\n", encoding="utf-8")
    (tmp_path / "CMakeLists.txt").write_text("add_library(x LIB/used.cpp)\n", encoding="utf-8")

    result = vs.validate_sources(root=tmp_path, strict=True)
    assert any(p.name == "EnvFlags.h" for p in result.orphans)
