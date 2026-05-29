#!/usr/bin/env python3
"""
update_skill.py

Part of the flexaid-docking skill.

Built-in, safe, professional autoupdate for the flexaid-docking skill and its
sub-components (scripts, references, documentation, data matrices, bin shortcuts).

This makes it dramatically easier for users to keep the skill current without
manual rsync, copy, or re-zip steps.

Primary use cases:
- Update the active Grok skill installation (~/.grok/skills/flexaid-docking)
- Update a portable copy (the folder you share with Claude Code / Cursor / Aider)
- Refresh a skills-export/ tree before re-zipping for distribution

Safety & Design Principles (non-negotiable):
- Dry-run by default. Nothing is modified unless you explicitly request it.
- Pure symlinks in bin/ (never wrappers or duplicated logic).
- Always ends with the skill validator (unless --no-validate).
- Clear distinction between "I have a full FlexAIDδS checkout" vs "I only have the portable folder".
- Banners and shortcuts are purely informational. No scientific claim is ever valid
  without actually running the FlexAIDδS binary and the thermodynamic ledger code.
"""

import argparse
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

# =============================================================================
# CONFIGURATION
# =============================================================================

SKILL_NAME = "flexaid-docking"
VERSION = "2026-05"

# Files and directories that are part of the skill "subskills" surface
CORE_ITEMS = [
    "SKILL.md",
    "AGENTS.md",
    "CLAUDE.md",
    "PORTING.md",
    "QUICKSTART.md",
    "scripts/",
    "references/",
    "bin/",
    "data/README.md",   # We are conservative with data/ — only README by default
]

DATA_ITEMS = [
    "data/MC_st0r5.2_6.dat",
    "data/MC_10p_3.dat",
    "data/MC_5p_norm_P10_M2_2.dat",
]

# Relative to the script location
THIS_SCRIPT = Path(__file__).resolve()
SKILL_ROOT = THIS_SCRIPT.parent.parent   # .grok/skills/flexaid-docking or portable equivalent

# Common places a full FlexAIDδS checkout might live
CANONICAL_SEARCH_PATHS: List[Path] = [
    Path.home() / "FlexAIDdS",
    Path.home() / "flexaidds",
    Path.home() / "projects" / "FlexAIDdS",
    Path.home() / "work" / "FlexAIDdS",
    Path(__file__).resolve().parents[4],           # typical worktree depth
    Path(__file__).resolve().parents[5],
    Path.cwd().parent,
    Path.cwd().parent.parent,
]


# =============================================================================
# BANNER (Purely Informational — same style as ensure_docking_data.py)
# =============================================================================

def print_skill_banner(verbose: bool = False) -> None:
    """Print a clean, professional version banner. Informational only."""
    script = "update_skill.py"

    if verbose:
        print("+" + "-" * 60 + "+")
        print("|  flexaid-docking skill                                   |")
        print(f"|  Script   : {script:<44}|")
        print(f"|  Version  : {VERSION:<44}|")
        print("|                                                          |")
        print("|  Built-in autoupdate for skill + sub-components.         |")
        print("|  Shortcuts and banners are informational only.           |")
        print("|  No scientific claim is valid without running the code.  |")
        print("+" + "-" * 60 + "+")
    else:
        print(f"[flexaid-docking] {script}  v{VERSION}")


# =============================================================================
# SOURCE DETECTION
# =============================================================================

@dataclass
class SourceInfo:
    root: Path
    is_full_checkout: bool
    confidence: str


def find_canonical_source(explicit_source: Optional[Path] = None) -> Optional[SourceInfo]:
    """
    Locate the best source of truth for the skill files.
    Priority:
      1. Explicit --source
      2. FLEXAIDDS_ROOT environment variable
      3. Heuristic search for a full FlexAIDdS checkout containing .grok/skills/flexaid-docking
      4. If the current location already looks like a full checkout, use it
    """
    if explicit_source:
        p = explicit_source.resolve()
        if (p / ".grok/skills/flexaid-docking").exists() or (p / ".grok/skills/flexaid-docking/SKILL.md").exists():
            return SourceInfo(root=p, is_full_checkout=True, confidence="explicit")
        # Allow pointing directly at the skill dir inside a checkout
        if (p / "SKILL.md").exists() and p.name == "flexaid-docking":
            return SourceInfo(root=p.parent.parent.parent, is_full_checkout=True, confidence="explicit-skill-dir")
        print(f"[WARN] --source {explicit_source} does not look like a FlexAIDδS checkout.")
        return None

    # Environment variable (power users)
    env_root = Path.home() if "FLEXAIDDS_ROOT" not in __import__("os").environ else Path(__import__("os").environ["FLEXAIDDS_ROOT"])
    if env_root and (env_root / ".grok/skills/flexaid-docking/SKILL.md").exists():
        return SourceInfo(root=env_root, is_full_checkout=True, confidence="env")

    # Current working tree / parents (very common when user is inside the repo)
    for candidate in [Path.cwd()] + list(Path.cwd().parents):
        if (candidate / ".grok/skills/flexaid-docking/SKILL.md").exists():
            return SourceInfo(root=candidate, is_full_checkout=True, confidence="cwd-parent")

    # Heuristic search
    for base in CANONICAL_SEARCH_PATHS:
        if not base.exists():
            continue
        if (base / ".grok/skills/flexaid-docking/SKILL.md").exists():
            return SourceInfo(root=base, is_full_checkout=True, confidence="heuristic")

    # Last resort: maybe the user only has the portable skill and no full repo
    return None


def get_skill_target() -> Path:
    """Return the root of the skill installation we are currently running from."""
    return SKILL_ROOT


# =============================================================================
# UPDATE LOGIC
# =============================================================================

def compute_file_list(source_root: Path, include_data: bool = False) -> List[Path]:
    """Return list of relative paths that should be considered for update."""
    items = list(CORE_ITEMS)
    if include_data:
        items.extend(DATA_ITEMS)

    rel_paths: List[Path] = []
    for item in items:
        p = source_root / ".grok/skills/flexaid-docking" / item
        if p.exists():
            if p.is_dir():
                for f in p.rglob("*"):
                    if f.is_file():
                        rel_paths.append(f.relative_to(source_root / ".grok/skills/flexaid-docking"))
            else:
                rel_paths.append(Path(item))
    return sorted(set(rel_paths))


def sync_item(src: Path, dst: Path, dry_run: bool, use_links: bool = False) -> Tuple[bool, str]:
    """Copy or link a single file/dir. Returns (changed, action_description)."""
    if not src.exists():
        return False, "missing-in-source"

    if dst.exists():
        # Very simple mtime + size heuristic (good enough for this purpose)
        try:
            if dst.is_file() and src.is_file():
                if dst.stat().st_mtime >= src.stat().st_mtime and dst.stat().st_size == src.stat().st_size:
                    return False, "up-to-date"
        except Exception:
            pass

    action = "symlink" if use_links else "copy"
    if dry_run:
        return True, f"[dry-run] would {action}: {src} -> {dst}"

    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            if dst.is_dir():
                shutil.rmtree(dst)
            else:
                dst.unlink()

        if use_links and src.is_file():
            try:
                dst.symlink_to(src.resolve())
                return True, f"linked: {dst.name}"
            except OSError:
                shutil.copy2(src, dst)
                return True, f"copied (symlink failed): {dst.name}"
        else:
            if src.is_dir():
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)
            return True, f"updated: {dst}"
    except Exception as e:
        return False, f"ERROR: {e}"


def perform_update(
    source: SourceInfo,
    target: Path,
    dry_run: bool = True,
    include_data: bool = False,
    use_links: bool = False,
    verbose: bool = False,
) -> Tuple[bool, List[str]]:
    """
    Main sync routine.
    Returns (overall_success, list_of_human_readable_actions)
    """
    actions: List[str] = []
    source_skill_dir = source.root / ".grok/skills/flexaid-docking"
    if not source_skill_dir.exists():
        print(f"[ERROR] Source skill directory not found: {source_skill_dir}")
        return False, actions

    rel_files = compute_file_list(source.root, include_data=include_data)

    print(f"\nSource : {source.root} (confidence: {source.confidence})")
    print(f"Target : {target}")
    print(f"Mode   : {'DRY-RUN (no changes)' if dry_run else 'LIVE UPDATE'}")
    print(f"Items  : {len(rel_files)} files/dirs considered" + (" (+ data matrices)" if include_data else ""))
    print()

    success = True
    for rel in rel_files:
        src = source_skill_dir / rel
        dst = target / rel
        changed, msg = sync_item(src, dst, dry_run=dry_run, use_links=use_links)
        if changed or verbose:
            actions.append(msg)
            if verbose or "ERROR" in msg or "would" in msg:
                print(f"  {msg}")

    if not dry_run and actions:
        print(f"\n[SUCCESS] {len([a for a in actions if 'updated' in a or 'linked' in a or 'copied' in a])} item(s) refreshed.")

    return success, actions


def run_validator(target: Path, verbose: bool = False) -> bool:
    """Run the skill's own validator on the target location."""
    validator = target / "scripts" / "validate_skill.py"
    if not validator.exists():
        print("[WARN] Could not find validator in target. Skipping validation step.")
        return True

    print("\n--- Running skill validator on target ---")
    try:
        result = subprocess.run(
            [sys.executable, str(validator)],
            cwd=target,
            capture_output=True,
            text=True,
            timeout=60,
        )
        print(result.stdout)
        if result.returncode != 0:
            print(result.stderr)
            return False
        return "VALIDATION PASSED" in result.stdout
    except Exception as e:
        print(f"[ERROR] Validator failed to run: {e}")
        return False


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Safe, built-in autoupdate for the flexaid-docking skill and all sub-components.",
        epilog=(
            "Part of the flexaid-docking skill. "
            "Recommended: start with --dry-run. "
            "Requires a full FlexAIDδS checkout as source (or --source)."
        ),
    )
    parser.add_argument(
        "--source", "-s", type=Path,
        help="Path to a full FlexAIDδS checkout (the canonical source of truth)"
    )
    parser.add_argument(
        "--target", "-t", type=Path,
        help="Explicit target skill directory to update (defaults to the currently running skill location)"
    )
    parser.add_argument(
        "--export", action="store_true",
        help="Treat target as a skills-export tree (adjusts paths accordingly)"
    )
    parser.add_argument(
        "--data", action="store_true",
        help="Also update the bundled interaction matrices (use with care)"
    )
    parser.add_argument(
        "--link", action="store_true",
        help="Prefer symlinks instead of copies where possible (bin/ entries stay as symlinks)"
    )
    parser.add_argument(
        "--dry-run", "-n", action="store_true", default=True,
        help="Preview changes only (default, recommended)"
    )
    parser.add_argument(
        "--yes", "-y", "--force", dest="force", action="store_true",
        help="Actually perform the update (disables the safe dry-run default)"
    )
    parser.add_argument(
        "--no-validate", action="store_true",
        help="Skip running the validator after the update"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Show detailed file-by-file actions")
    parser.add_argument(
        "--check-only", action="store_true",
        help="Only report what would be different (implies --dry-run)"
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    print_skill_banner(verbose=args.verbose)

    # Resolve target
    target = args.target.resolve() if args.target else get_skill_target()

    # If user passed --force we turn off the default dry-run
    dry_run = False if args.force else True
    if args.check_only:
        dry_run = True

    print(f"\nCurrent skill location : {get_skill_target()}")
    print(f"Update target          : {target}")

    source_info = find_canonical_source(args.source)
    if not source_info:
        print("\n[ERROR] Could not locate a full FlexAIDδS checkout as source.")
        print("        Use --source /path/to/your/FlexAIDdS/checkout")
        print("        or set FLEXAIDDS_ROOT, or run this script from inside a checkout.")
        print("\nTip: The portable skill folder alone is not sufficient for updates.")
        print("     Keep a full clone of https://github.com/LeBonhommePharma/FlexAIDdS for the best experience.")
        return 2

    success, actions = perform_update(
        source=source_info,
        target=target,
        dry_run=dry_run,
        include_data=args.data,
        use_links=args.link,
        verbose=args.verbose,
    )

    if not success:
        return 1

    if dry_run and actions:
        print("\n[DRY-RUN] No files were modified. Re-run with --yes (or -y) to apply changes.")
        print("          Strongly recommended: inspect the listed actions first.")

    # Post-update validation (the skill demands this)
    if not args.no_validate and not dry_run:
        validator_ok = run_validator(target, verbose=args.verbose)
        if not validator_ok:
            print("\n[WARNING] Validator reported issues after update. Please investigate.")
            return 1

    if not dry_run:
        print("\n[SUCCESS] Skill update complete. Remember the mandatory ritual:")
        print("          git status + find + python3 .../validate_skill.py before any real work.")

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
