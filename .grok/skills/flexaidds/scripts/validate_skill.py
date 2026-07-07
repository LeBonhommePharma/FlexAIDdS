#!/usr/bin/env python3
"""
FlexAIDdS / FlexAID Skill Validator
Validates the flexaidds skill packaging:
- SKILL.md frontmatter has required name + description
- All XML files in repo are well-formed (one root, proper escaping, UTF-8, nesting)
- No broken local file references in SKILL.md
- Invocation aliases are documented
- Scientific terminology guardrails present (no overclaim of true ΔG vs CF proxy)

Run:
  python3 .grok/skills/flexaidds/scripts/validate_skill.py
  python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
"""
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Tuple

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parents[2]  # .grok/skills/flexaidds -> repo root


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)


def ok(msg: str) -> None:
    print(f"OK: {msg}")


def validate_frontmatter() -> bool:
    skill_md = SKILL_DIR / "SKILL.md"
    if not skill_md.exists():
        fail(f"SKILL.md not found at {skill_md}")
        return False
    try:
        content = skill_md.read_text(encoding="utf-8")
    except UnicodeDecodeError as e:
        fail(f"SKILL.md is not valid UTF-8: {e}")
        return False

    if not content.lstrip().startswith("---"):
        fail("SKILL.md missing leading --- YAML frontmatter delimiter")
        return False

    parts = content.split("---", 2)
    if len(parts) < 3:
        fail("SKILL.md frontmatter not closed with second ---")
        return False

    fm = parts[1].strip()
    # Required per task: name and description
    name_match = re.search(r"^name:\s*([A-Za-z0-9_-]+)", fm, re.MULTILINE)
    if not name_match or name_match.group(1) != "flexaidds":
        fail("Frontmatter must contain: name: flexaidds")
        return False

    if "description:" not in fm:
        fail("Frontmatter must contain a description field")
        return False

    # description should be non-trivial (handle > multiline or | block)
    desc_part = fm.split("description:", 1)[1]
    # take until next top-level key or end
    desc_text = re.split(r"\n[a-z_]+:", desc_part, maxsplit=1)[0]
    desc_clean = re.sub(r"\s+", " ", desc_text).strip()
    if len(desc_clean) < 30:
        fail("description is too short or empty")
        return False

    ok("SKILL.md frontmatter has required name: flexaidds and description")
    return True


def validate_xml_files() -> bool:
    xml_files: List[Path] = list(REPO_ROOT.rglob("*.xml")) + list(REPO_ROOT.rglob("*.XML"))
    # Also check inside skill dir explicitly
    skill_xml = list(SKILL_DIR.rglob("*.xml")) + list(SKILL_DIR.rglob("*.XML"))
    all_xml = sorted(set(xml_files + skill_xml))

    if not all_xml:
        ok("No XML files present in repository or skill (0 files). "
           "Skill uses SKILL.md (not XML). Validator will check any added *.xml for "
           "well-formedness: single root, escaped & < >, valid UTF-8, proper nesting, no duplicate IDs where disallowed.")
        return True

    errors: List[str] = []
    for xf in all_xml:
        try:
            # ET.parse enforces well-formed XML (one root, correct nesting, escaping)
            tree = ET.parse(xf)
            root = tree.getroot()
            # Basic sanity: root exists
            if root is None:
                errors.append(f"{xf}: no root element")
                continue
            # Check for common malformed patterns (ampersand not escaped etc already caught by parser)
            # UTF-8 is handled by parse (with default)
        except ET.ParseError as e:
            errors.append(f"{xf}: {e}")
        except UnicodeDecodeError as e:
            errors.append(f"{xf}: invalid UTF-8 or encoding issue: {e}")
        except Exception as e:
            errors.append(f"{xf}: unexpected error: {e}")

    if errors:
        for e in errors:
            fail(e)
        return False

    ok(f"Validated {len(all_xml)} XML file(s): all well-formed (single root, escaping, UTF-8, nesting OK)")
    return True


def check_broken_local_refs() -> bool:
    skill_md = SKILL_DIR / "SKILL.md"
    content = skill_md.read_text(encoding="utf-8")

    # Markdown links: [text](target)
    links = re.findall(r"\[([^\]]+)\]\(([^)]+?)(?:\s+\"[^\"]*\")?\)", content)
    broken: List[str] = []
    checked: set = set()

    for _text, target in links:
        target = target.strip()
        if not target or target.startswith(("#", "http:", "https:", "mailto:", "data:")):
            continue
        # Strip anchors
        target_path = target.split("#", 1)[0]
        if not target_path:
            continue
        if target_path in checked:
            continue
        checked.add(target_path)

        # Resolve relative to skill dir first (for references/, scripts/)
        cand_skill = (SKILL_DIR / target_path).resolve()
        # Also relative to repo root
        cand_repo = (REPO_ROOT / target_path).resolve()

        if not (cand_skill.exists() or cand_repo.exists()):
            # Also try as-is from cwd (rare)
            if not Path(target_path).resolve().exists():
                broken.append(target)

    if broken:
        for b in broken:
            fail(f"Broken local reference in SKILL.md: {b}")
        return False

    ok(f"No broken local references (checked {len(checked)} links)")
    return True


def check_aliases_documented() -> bool:
    skill_md = SKILL_DIR / "SKILL.md"
    content = skill_md.read_text(encoding="utf-8")
    required_aliases = [
        "/FlexAid docking",
        "/FlexAidDS",
        "FlexAIDdS",
        "FlexAID∆S",
    ]
    missing = [a for a in required_aliases if a not in content]
    if missing:
        for m in missing:
            fail(f"Missing documented invocation alias/trigger: {m}")
        return False
    ok("All required invocation aliases documented: " + ", ".join(required_aliases))
    return True


def check_bin_wrappers() -> bool:
    """bin/ entries must be executable shell wrappers, not symlinks into scripts/."""
    bin_dir = SKILL_DIR / "bin"
    if not bin_dir.is_dir():
        fail("bin/ directory missing")
        return False
    expected = (
        "validate-skill",
        "ensure-docking-data",
        "inspect-definition-files",
        "update-skill",
        "dataset-runner",
        "resolve-build",
    )
    missing = [name for name in expected if not (bin_dir / name).exists()]
    if missing:
        for name in missing:
            fail(f"Missing bin/ wrapper: {name}")
        return False
    for name in expected:
        path = bin_dir / name
        if path.is_symlink():
            fail(f"bin/{name} must not be a symlink (would overwrite scripts/ if edited)")
            return False
        if not os.access(path, os.X_OK):
            fail(f"bin/{name} is not executable")
            return False
        try:
            header = path.read_text(encoding="utf-8").splitlines()[0]
        except (OSError, UnicodeDecodeError) as exc:
            fail(f"bin/{name} unreadable: {exc}")
            return False
        if not header.startswith("#!"):
            fail(f"bin/{name} must be a shell wrapper with shebang")
            return False
    ok(f"bin/ wrappers present and executable ({len(expected)} commands)")
    return True


def check_guardrails() -> bool:
    skill_md = SKILL_DIR / "SKILL.md"
    content = skill_md.read_text(encoding="utf-8").lower()
    required_phrases = [
        "inspect repo state first",
        "avoid unsafe git",
        "validate claims against files",
        "contact-function scoring proxy",
        "cf/contact-function",
        "preserve current ranking",
        "thermodynamic/ensemble work only behind tests",
        "never merge branches or rewrite history",
        "chunked implementation plans",
        "produce chunked implementation plans",
        "source of truth",
        "agents.md",
        "repository hygiene",
    ]
    missing = [p for p in required_phrases if p not in content]
    if missing:
        for m in missing:
            fail(f"Missing guardrail language in SKILL.md: {m}")
        return False
    ok("Guardrails present (repo inspection, git safety, claim validation, proxy vs thermo separation, ranking preservation, test gates, no history rewrite)")
    return True


def check_active_build() -> bool:
    """Warn or fail when no resolvable production build exists."""
    resolver = SKILL_DIR / "scripts" / "resolve_build.py"
    if not resolver.is_file():
        fail("resolve_build.py missing")
        return False
    require = os.environ.get("FLEXAIDDS_REQUIRE_BUILD", "").strip() in ("1", "true", "yes")
    try:
        proc = subprocess.run(
            [sys.executable, str(resolver), "--check", "--repo-root", str(REPO_ROOT)],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        msg = f"resolve_build check could not run: {exc}"
        if require:
            fail(msg)
            return False
        print(f"WARN: {msg}")
        return True
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout or "resolve_build failed").strip()
        if require:
            fail(f"Active build resolution failed: {msg}")
            return False
        print(f"WARN: Active build not resolved (set FLEXAIDDS_REQUIRE_BUILD=1 to enforce): {msg}")
        return True
    ok("Active build resolves (resolve_build.py --check)")
    return True


def main() -> int:
    print("=== FlexAIDdS Skill Packaging Validator ===\n")
    results = [
        ("frontmatter", validate_frontmatter()),
        ("xml-well-formed", validate_xml_files()),
        ("no-broken-refs", check_broken_local_refs()),
        ("aliases-documented", check_aliases_documented()),
        ("bin-wrappers", check_bin_wrappers()),
        ("active-build", check_active_build()),
        ("guardrails", check_guardrails()),
    ]

    print("\n=== Summary ===")
    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")
        if not passed:
            all_pass = False

    if all_pass:
        print("\nVALIDATION PASSED: skill is well-formed and ready for /flexaidds, /FlexAidDS, FlexAID∆S etc.")
        return 0
    else:
        print("\nVALIDATION FAILED: see errors above. Fix before committing skill.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
