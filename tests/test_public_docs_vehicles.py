#!/usr/bin/env python3
"""Public-vehicle honesty: landing README, install how-to, user guide, changelog, atom types.

Drives the shipped docs and the shipped config parser / typing functions.
Does not re-implement mapping tables or invent JSON keys.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PUBLIC_VEHICLES = (
    ROOT / "README.md",
    ROOT / "docs" / "INSTALL.md",
    ROOT / "docs" / "USERGUIDE.md",
    ROOT / "VERSION.md",
    ROOT / "docs" / "ATOM_TYPES.md",
)

REMAP_STRINGS = ("N.2→N.ar", "N.3→N.am", "C.1→C.2", "I→BR")

# Live VCT rows returned by the mapping *functions* (not comments).
MOL2_TOP_REMAPS = (("N.2", 10), ("N.3", 11), ("C.1", 2), ("I", 25))
SDF_IODINE_REMAP = ("I", 25)

JSON_FENCE_RE = re.compile(r"```json\s*\n(.*?)```", re.S)
APPLIED_KEY_RE = re.compile(
    r"j(?:bool|int|dbl|flt|str)\(\s*config\s*,\s*\"([A-Za-z0-9_]+)\"\s*,\s*\"([A-Za-z0-9_]+)\""
)
APPLIED_INDEX_RE = re.compile(
    r"config\[\"([A-Za-z0-9_]+)\"\]\[\"([A-Za-z0-9_]+)\"\]"
)
STRCMP_RETURN_RE = re.compile(
    r"strcmp\([^,]+,\s*\"([^\"]+)\"\)\)\s*return\s+(\d+)"
)
DEAD_JSON_KEYS = frozenset({"ring_conformers", "chirality"})
KELVIN_DEFAULT_RE = re.compile(
    r"(?:entropy at 300\s*K|300\s*K entropy|flexibility at 300\s*K|"
    r"Temperature in K(?:elvin)?(?:\s|\||$)|"
    r"defaults enable full flexibility at 300\s*K)",
)


def _applied_parser_keys(parser_src: str) -> set[tuple[str, str]]:
    keys = set(APPLIED_KEY_RE.findall(parser_src))
    keys.update(APPLIED_INDEX_RE.findall(parser_src))
    return keys


def _function_body(source: str, name: str) -> str:
    match = re.search(rf"\b{re.escape(name)}\s*\([^;]*?\)\s*\{{", source)
    assert match, f"missing function {name}"
    start = match.end() - 1
    depth = 0
    for index in range(start, len(source)):
        char = source[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return source[start : index + 1]
    raise AssertionError(f"unclosed function {name}")


def _strip_cpp_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//.*?$", "", text, flags=re.M)


def _mapping_returns(path: Path, func_name: str) -> dict[str, int]:
    body = _strip_cpp_comments(_function_body(path.read_text(encoding="utf-8"), func_name))
    found = {key: int(val) for key, val in STRCMP_RETURN_RE.findall(body)}
    assert found, f"no strcmp/return pairs in {path.name}:{func_name}"
    return found


def test_landing_readme_install_links_share_one_howto() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    header = re.search(r"\*\*\[Installation\]\(([^)]+)\)\*\*", readme)
    assert header, "GitHub landing README is missing the header Installation link"
    header_target = header.group(1)

    quick = re.search(
        r"## Quick Start.*?\[`([^`]+)`\]\(([^)]+)\)",
        readme,
        flags=re.S,
    )
    assert quick, "Quick Start is missing an install how-to path"
    label, quick_target = quick.group(1), quick.group(2)

    assert header_target == quick_target, (
        f"header Installation={header_target!r} vs Quick Start={quick_target!r}"
    )
    assert label == header_target, (
        f"Quick Start label {label!r} must match the how-to path {header_target!r}"
    )
    assert "INSTALLATION.md" not in header_target

    howto = ROOT / header_target
    assert howto.is_file(), f"install how-to does not exist: {header_target}"
    text = howto.read_text(encoding="utf-8")
    assert "Python package" in text
    assert "separate" in text.lower()
    assert (
        "does **not** contain the docking engine" in text
        or "does not contain the docking engine" in text
    )
    assert "CMake" in text and "Homebrew" in text and "Docker" in text
    assert "pip" in text.lower()


def test_user_guide_json_examples_use_only_applied_parser_keys() -> None:
    guide = (ROOT / "docs" / "USERGUIDE.md").read_text(encoding="utf-8")
    parser_src = (ROOT / "LIB" / "config_parser.cpp").read_text(encoding="utf-8")
    applied = _applied_parser_keys(parser_src)
    assert ("ga", "num_chromosomes") in applied
    assert ("thermodynamics", "temperature") in applied
    assert ("flexibility", "ring_conformers") not in applied
    assert ("flexibility", "chirality") not in applied

    fences = JSON_FENCE_RE.findall(guide)
    assert fences, "user guide has no fenced JSON examples"
    unknown: list[str] = []
    for block in fences:
        payload = json.loads(block)
        assert isinstance(payload, dict), "JSON examples must be objects"
        for section, body in payload.items():
            assert isinstance(body, dict), f"section {section!r} must be an object"
            for key in body:
                assert key not in DEAD_JSON_KEYS, (
                    f"user guide advertises ignored JSON key {section}.{key}"
                )
                if (section, key) not in applied:
                    unknown.append(f"{section}.{key}")
    assert unknown == [], (
        "user guide JSON keys not applied by LIB/config_parser.cpp: "
        + ", ".join(unknown)
    )


def test_user_guide_does_not_call_default_temperature_kelvin() -> None:
    text = (ROOT / "docs" / "USERGUIDE.md").read_text(encoding="utf-8")
    match = KELVIN_DEFAULT_RE.search(text)
    assert match is None, (
        "user guide still describes default docking temperature as Kelvin / "
        f"300 K entropy: {match.group(0)!r}"
    )
    assert "300 K" not in text
    assert "score-scale" in text.lower()


def test_four_remaps_match_mapping_functions_and_readme_table() -> None:
    mol2 = _mapping_returns(ROOT / "LIB" / "Mol2Reader.cpp", "sybyl_to_flexaid_type")
    top = _mapping_returns(ROOT / "LIB" / "top_helpers.cpp", "sybyl_name_to_canonical_vct")
    sdf = _mapping_returns(ROOT / "LIB" / "SdfReader.cpp", "element_to_flexaid_type")

    for key, expected in MOL2_TOP_REMAPS:
        assert mol2[key] == expected, f"Mol2Reader {key} -> {mol2.get(key)}"
        assert top[key] == expected, f"sybyl_name_to_canonical_vct {key} -> {top.get(key)}"
    assert sdf[SDF_IODINE_REMAP[0]] == SDF_IODINE_REMAP[1]

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    atom_types = (ROOT / "docs" / "ATOM_TYPES.md").read_text(encoding="utf-8")
    for remap in REMAP_STRINGS:
        assert remap in readme, f"README remap table missing {remap}"
        assert remap in atom_types, f"atom-types page missing {remap}"
    assert "MC_st0r5.2_6.dat" in atom_types
    assert "40" in atom_types


def test_python_package_version_remains_2_0_3() -> None:
    version_py = (ROOT / "python" / "flexaidds" / "__version__.py").read_text(
        encoding="utf-8"
    )
    changelog = (ROOT / "VERSION.md").read_text(encoding="utf-8")
    assert '__version__ = "2.0.3"' in version_py
    assert "**Current Version**: 2.0.3" in changelog
    assert "Unreleased (post-v2.0.3)" in changelog


def test_public_vehicles_exist() -> None:
    for path in PUBLIC_VEHICLES:
        assert path.is_file(), f"missing public vehicle {path}"
