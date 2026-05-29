#!/usr/bin/env python3
"""Validate the flexaid-docking skill package.

Checks:
  1. Required skill/package files and directories exist.
  2. SKILL.md has YAML frontmatter with exact `name` and `description`.
  3. SKILL.md documents all required FlexAID/FlexAIDdS aliases.
  4. Every XML file under the skill is well-formed UTF-8 with one root.
  5. The XML manifest mirrors the canonical metadata and aliases.
  6. XML IDs and aliases are unique.
  7. Every local Markdown link and declared file reference resolves on disk.

Exit code is 0 on success, 1 on any failure. Designed to run with no
dependencies beyond the Python 3 standard library.
"""
from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parent.parent
SKILL_MD = SKILL_DIR / "SKILL.md"
MANIFEST_XML = SKILL_DIR / "references" / "skill-manifest.xml"

REQUIRED_NAME = "flexaid-docking"
REQUIRED_DESCRIPTION = (
    "Use this skill for FlexAID, FlexAIDdS, and FlexAID∆S docking workflows, "
    "including safe repo review, implementation planning, XML/package "
    "validation, and docking/thermodynamic-roadmap task decomposition."
)
REQUIRED_ALIASES = (
    "/FlexAid docking",
    "/FlexAidDS",
    "FlexAIDdS",
    "FlexAID∆S",
)
REQUIRED_LAYOUT = (
    SKILL_MD,
    SKILL_DIR / "references",
    SKILL_DIR / "references" / "flexaid-docking-guidance.md",
    MANIFEST_XML,
    SKILL_DIR / "scripts",
    SKILL_DIR / "scripts" / "validate_skill.py",
    SKILL_DIR / "assets",
    REPO_ROOT / "AGENTS.md",
)


class ValidationError(Exception):
    """Raised for any single validation failure."""


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise ValidationError(f"missing file: {_display_path(path)}") from exc
    except UnicodeDecodeError as exc:
        raise ValidationError(
            f"file is not valid UTF-8: {_display_path(path)}"
        ) from exc


def _normalize_ws(text: str) -> str:
    return " ".join(text.split())


def _parse_frontmatter(text: str) -> dict[str, str]:
    """Minimal YAML-ish frontmatter parser: top-level `key: value` only."""
    if not text.startswith("---"):
        raise ValidationError("SKILL.md does not start with a `---` frontmatter fence")
    parts = text.split("---", 2)
    if len(parts) < 3:
        raise ValidationError("SKILL.md frontmatter has no closing `---` fence")
    body = parts[1]
    meta: dict[str, str] = {}
    for raw_line in body.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValidationError(
                f"frontmatter line missing colon: {raw_line!r}"
            )
        key, _, value = line.partition(":")
        meta[key.strip()] = value.strip()
    return meta


def check_skill_layout() -> list[str]:
    issues: list[str] = []
    for path in REQUIRED_LAYOUT:
        if not path.exists():
            issues.append(f"required skill path is missing: {_display_path(path)}")
    return issues


def check_skill_md() -> list[str]:
    issues: list[str] = []
    text = _read_text(SKILL_MD)
    meta = _parse_frontmatter(text)

    if "name" not in meta:
        issues.append("SKILL.md frontmatter missing `name`")
    elif meta["name"] != REQUIRED_NAME:
        issues.append(
            f"SKILL.md `name` is {meta['name']!r}, expected {REQUIRED_NAME!r}"
        )

    if "description" not in meta or not meta["description"]:
        issues.append("SKILL.md frontmatter missing or empty `description`")
    elif _normalize_ws(meta["description"]) != REQUIRED_DESCRIPTION:
        issues.append("SKILL.md `description` does not match required metadata")

    for alias in REQUIRED_ALIASES:
        if alias not in text:
            issues.append(f"SKILL.md does not mention required alias: {alias!r}")

    return issues


_PATH_TABLE_ROW = re.compile(r"^\|\s*`([^`]+)`\s*\|", re.MULTILINE)


def referenced_paths_from_skill_md() -> list[Path]:
    """Pull `path-like` entries out of the markdown file table in SKILL.md."""
    text = _read_text(SKILL_MD)
    paths: list[Path] = []
    for match in _PATH_TABLE_ROW.finditer(text):
        raw = match.group(1).strip()
        if "/" not in raw and not raw.endswith(".md"):
            continue
        if raw.endswith("/"):
            paths.append(SKILL_DIR / raw.rstrip("/"))
        else:
            paths.append(SKILL_DIR / raw)
    return paths


def check_referenced_paths_exist() -> list[str]:
    issues: list[str] = []
    for path in referenced_paths_from_skill_md():
        if not path.exists():
            issues.append(
                f"SKILL.md file table references missing path: {_display_path(path)}"
            )
    return issues


_MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]+\]\(([^)]+)\)")
_EXTERNAL_LINK_PREFIXES = (
    "http://",
    "https://",
    "mailto:",
    "tel:",
    "#",
)


def iter_markdown_files() -> Iterable[Path]:
    yield from SKILL_DIR.rglob("*.md")


def _local_link_target(raw_target: str) -> str | None:
    target = raw_target.strip()
    if not target:
        return None
    if target.startswith("<") and ">" in target:
        target = target[1 : target.index(">")]
    else:
        # Strip optional Markdown title text in simple inline links.
        target = target.split(maxsplit=1)[0]
    if target.startswith(_EXTERNAL_LINK_PREFIXES):
        return None
    target = target.split("#", 1)[0]
    if not target:
        return None
    return unquote(target)


def check_markdown_links() -> list[str]:
    issues: list[str] = []
    for path in iter_markdown_files():
        text = _read_text(path)
        for match in _MARKDOWN_LINK.finditer(text):
            target = _local_link_target(match.group(1))
            if target is None:
                continue
            resolved = (path.parent / target).resolve()
            if not resolved.exists():
                issues.append(
                    f"{_display_path(path)} references missing local link: {target}"
                )
    return issues


def iter_xml_files() -> Iterable[Path]:
    for path in SKILL_DIR.rglob("*.xml"):
        yield path


def _manifest_child_text(root: ET.Element, tag: str) -> str:
    child = root.find(tag)
    return "" if child is None or child.text is None else _normalize_ws(child.text)


def _check_manifest_metadata(path: Path, root: ET.Element) -> list[str]:
    issues: list[str] = []
    rel = _display_path(path)

    if root.tag != "skill":
        issues.append(f"{rel}: root element is <{root.tag}>, expected <skill>")
    if root.get("id") != REQUIRED_NAME:
        issues.append(f"{rel}: root id is {root.get('id')!r}, expected {REQUIRED_NAME!r}")
    if _manifest_child_text(root, "name") != REQUIRED_NAME:
        issues.append(f"{rel}: <name> does not match {REQUIRED_NAME!r}")
    if _manifest_child_text(root, "description") != REQUIRED_DESCRIPTION:
        issues.append(f"{rel}: <description> does not match SKILL.md metadata")

    alias_values = {
        (alias_el.text or "").strip()
        for alias_el in root.findall("./aliases/alias")
        if (alias_el.text or "").strip()
    }
    for alias in REQUIRED_ALIASES:
        if alias not in alias_values:
            issues.append(f"{rel}: missing required alias {alias!r}")

    return issues


def check_xml_files() -> list[str]:
    issues: list[str] = []
    saw_any = False
    for path in iter_xml_files():
        saw_any = True
        rel = _display_path(path)
        try:
            raw = path.read_bytes()
        except OSError as exc:
            issues.append(f"{rel}: unreadable ({exc})")
            continue
        try:
            raw.decode("utf-8")
        except UnicodeDecodeError as exc:
            issues.append(f"{rel}: not valid UTF-8 ({exc})")
            continue
        try:
            tree = ET.parse(path)
        except ET.ParseError as exc:
            issues.append(f"{rel}: XML parse error: {exc}")
            continue
        root = tree.getroot()
        if root is None:
            issues.append(f"{rel}: no root element")
            continue

        if path == MANIFEST_XML:
            issues.extend(_check_manifest_metadata(path, root))

        seen_ids: dict[str, str] = {}
        for element in root.iter():
            value = (element.get("id") or "").strip()
            if not value:
                continue
            if value in seen_ids:
                issues.append(
                    f"{rel}: duplicate XML id {value!r} on <{element.tag}> "
                    f"(already used on <{seen_ids[value]}>)"
                )
            else:
                seen_ids[value] = element.tag

        # Validate any <file path="..."/> entries
        for file_el in root.iter("file"):
            ref = file_el.get("path")
            if not ref:
                issues.append(f"{rel}: <file> element missing path attribute")
                continue
            target = (path.parent / ref).resolve()
            if not target.exists():
                issues.append(
                    f"{rel}: <file path={ref!r}/> does not resolve on disk"
                )
        # Optional: check alias uniqueness if <aliases> present
        seen_aliases: set[str] = set()
        for alias_el in root.iter("alias"):
            value = (alias_el.text or "").strip()
            if not value:
                issues.append(f"{rel}: empty <alias/> element")
                continue
            if value in seen_aliases:
                issues.append(f"{rel}: duplicate alias {value!r}")
            seen_aliases.add(value)
    if not saw_any:
        # No XML is fine; skill is markdown-first. Just record an info note.
        pass
    return issues


def main() -> int:
    sections = [
        ("skill package layout", check_skill_layout),
        ("SKILL.md frontmatter / aliases", check_skill_md),
        ("SKILL.md file-table references", check_referenced_paths_exist),
        ("Markdown local links", check_markdown_links),
        ("XML well-formedness / manifest", check_xml_files),
    ]
    all_issues: list[str] = []
    for label, fn in sections:
        try:
            issues = fn()
        except ValidationError as exc:
            issues = [str(exc)]
        if issues:
            print(f"[FAIL] {label}")
            for issue in issues:
                print(f"  - {issue}")
            all_issues.extend(issues)
        else:
            print(f"[OK]   {label}")

    if all_issues:
        print(f"\nValidation failed with {len(all_issues)} issue(s).")
        return 1
    print("\nValidation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
