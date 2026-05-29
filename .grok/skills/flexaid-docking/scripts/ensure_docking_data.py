#!/usr/bin/env python3
"""
ensure_docking_data.py

Part of the flexaid-docking skill.

Ensures the critical precomputed interaction matrix files (MC_*.dat)
and definition files (*.def) required by the FlexAIDδS binary are available
in the location the binary expects them.

This is the primary, recommended tool for managing these otherwise-missing
runtime dependencies when working with FlexAIDδS through this skill.

Supports:
- Auto-discovery across common locations (including files bundled with the skill)
- Explicit copy from a known-good installation via --source
- Smart binary detection
- --check / --status mode
- --link mode for symlinks
- Clean professional output
"""

import argparse
import shutil
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# =============================================================================
# CONFIGURATION
# =============================================================================

SKILL_NAME = "flexaid-docking"
VERSION = "2026-05"

EXPECTED_MATRICES = [
    "MC_st0r5.2_6.dat",
]

# Definition files required by the FlexAIDδS binary (atom typing, nucleotides, etc.)
EXPECTED_DEF_FILES = [
    "AMINO.def",
    "AMINO8.def",
    "AMINO12.def",
    "AMINO26.def",
    "NUCLEOTIDES.def",
    "NUCLEOTIDES8.def",
    "NUCLEOTIDES12.def",
    "NUCLEOTIDES26.def",
]

# Additional critical runtime files commonly found in complete WRK/ data packs
# (rotamer libraries, scoring matrices, etc.). These are also required in the binary base path.
EXPECTED_EXTRA_FILES = [
    "Lovell_LIB.dat",
    "rotobs.lst",
    "SYBYL_emat.dat",
    "M6_cons_3.dat",
    "nrg_mat_BEST_011912.dat",
    "nrg_mat_BEST_012012.dat",
    "scr_bin.dat",
    "scr_mat.dat",
]

DEFAULT_SEARCH_PATHS: List[Path] = [
    # Bundled inside the skill (highest priority - makes the skill self-contained)
    Path(__file__).resolve().parent.parent / "data",
    Path.home() / ".flexaidds" / "data",
    Path("/usr/local/share/flexaidds/data"),
    Path("/opt/flexaidds/data"),
    Path.home() / "flexaidds-data",
    Path(__file__).resolve().parents[4] / "build",
]


# =============================================================================
# BANNER (Purely Informational)
# =============================================================================

def print_skill_banner(verbose: bool = False) -> None:
    """Print a clean, professional version banner. Informational only."""
    script = "ensure_docking_data.py"

    if verbose:
        print("+" + "-" * 60 + "+")
        print("|  flexaid-docking skill                                   |")
        print(f"|  Script   : {script:<44}|")
        print(f"|  Version  : {VERSION:<44}|")
        print("|                                                          |")
        print("|  This is part of the official skill.                     |")
        print("|  Shortcuts and banners are informational only.           |")
        print("|  No scientific claim is valid without running the code.  |")
        print("+" + "-" * 60 + "+")
    else:
        print(f"[flexaid-docking] {script}  v{VERSION}")


# =============================================================================
# CORE LOGIC
# =============================================================================

def find_matrix_files(search_roots: List[Path]) -> List[Path]:
    found = []
    for root in search_roots:
        if not root or not root.exists():
            continue
        for name in EXPECTED_MATRICES:
            direct = root / name
            if direct.is_file():
                found.append(direct)
            for candidate in root.rglob(name):
                if candidate.is_file():
                    found.append(candidate)
    seen = set()
    unique = []
    for p in found:
        if p not in seen:
            seen.add(p)
            unique.append(p)
    return unique


def find_def_files(search_roots: List[Path]) -> List[Path]:
    """Find *.def files (AMINO*, NUCLEOTIDES*) in the search roots."""
    found = []
    for root in search_roots:
        if not root or not root.exists():
            continue
        for name in EXPECTED_DEF_FILES:
            direct = root / name
            if direct.is_file():
                found.append(direct)
            for candidate in root.rglob(name):
                if candidate.is_file():
                    found.append(candidate)
    seen = set()
    unique = []
    for p in found:
        if p not in seen:
            seen.add(p)
            unique.append(p)
    return unique


def find_extra_files(search_roots: List[Path]) -> List[Path]:
    """Find additional runtime files (Lovell_LIB, rotobs, SYBYL_emat, scoring mats, etc.)."""
    found = []
    for root in search_roots:
        if not root or not root.exists():
            continue
        for name in EXPECTED_EXTRA_FILES:
            direct = root / name
            if direct.is_file():
                found.append(direct)
            for candidate in root.rglob(name):
                if candidate.is_file():
                    found.append(candidate)
    seen = set()
    unique = []
    for p in found:
        if p not in seen:
            seen.add(p)
            unique.append(p)
    return unique


def get_binary_base_path(binary: Optional[Path]) -> Path:
    if binary and binary.exists():
        return binary.parent.resolve()
    candidates = [
        Path.cwd(),
        Path(__file__).resolve().parents[4] / "build",
    ]
    for c in candidates:
        if (c / "FlexAIDδS").exists():
            return c.resolve()
    return Path.cwd()


def ensure_matrices(
    binary: Optional[Path] = None,
    source: Optional[Path] = None,
    dry_run: bool = False,
    use_links: bool = False,
    verbose: bool = False,
    check_only: bool = False,
) -> Tuple[bool, List[Path]]:
    binary_base = get_binary_base_path(binary)
    if verbose:
        print(f"[info] Binary base path: {binary_base}")

    search_roots = list(DEFAULT_SEARCH_PATHS)
    search_roots.append(binary_base)
    search_roots.append(binary_base.parent)

    if source:
        if not source.exists():
            print(f"[ERROR] --source path does not exist: {source}")
            return False, []
        search_roots = [source] + search_roots
        if verbose:
            print(f"[info] Using explicit source: {source}")

    found = find_matrix_files(search_roots)
    missing = [name for name in EXPECTED_MATRICES if not any(f.name == name for f in found)]

    if found:
        print(f"\nFound {len(found)} required matrix file(s):")
        for f in found:
            print(f"  + {f}")
    else:
        print("\n[ERROR] No required interaction matrix files found.")

    if missing:
        print(f"\nMissing: {', '.join(missing)}")

    if check_only:
        success = len(missing) == 0
        print(f"\n[check] {'READY' if success else 'DATA MISSING'}")
        return success, found

    if not missing:
        print("\n[SUCCESS] All required matrices are already available.")
        return True, found

    success = True
    for name in missing:
        candidates = [f for f in found if f.name == name]
        if not candidates:
            print(f"[ERROR] Cannot locate matrix: {name}")
            success = False
            continue

        src = candidates[0]
        dst = binary_base / name

        try:
            if dry_run:
                action = "symlink" if use_links else "copy"
                print(f"  [dry-run] Would {action}: {src} -> {dst}")
            else:
                dst.parent.mkdir(parents=True, exist_ok=True)
                if use_links:
                    try:
                        dst.symlink_to(src.resolve())
                        print(f"  Linked: {name}")
                    except OSError:
                        shutil.copy2(src, dst)
                        print(f"  Copied (symlink not supported): {name}")
                else:
                    shutil.copy2(src, dst)
                    print(f"  Copied: {name}")
        except Exception as e:
            print(f"  [ERROR] Failed to place {name}: {e}")
            success = False

    if success and not dry_run:
        print("\n[SUCCESS] All required matrices are now in place.")

    # --- Also ensure definition files (*.def) in the same base path ---
    def_found = find_def_files(search_roots)
    def_missing = [name for name in EXPECTED_DEF_FILES if not any(f.name == name for f in def_found)]

    if def_found:
        print(f"\nFound {len(def_found)} definition file(s):")
        for f in def_found:
            print(f"  + {f}")

    if def_missing:
        print(f"\nMissing definition files: {', '.join(def_missing)}")

    if check_only:
        overall_success = len(missing) == 0 and len(def_missing) == 0
        print(f"\n[check] {'READY' if overall_success else 'DATA MISSING'}")
        return overall_success, found + def_found

    if def_missing:
        for name in def_missing:
            candidates = [f for f in def_found if f.name == name]
            if not candidates:
                print(f"[ERROR] Cannot locate definition file: {name}")
                success = False
                continue

            src = candidates[0]
            dst = binary_base / name

            try:
                if dry_run:
                    action = "symlink" if use_links else "copy"
                    print(f"  [dry-run] Would {action}: {src} -> {dst}")
                else:
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    if use_links:
                        try:
                            dst.symlink_to(src.resolve())
                            print(f"  Linked: {name}")
                        except OSError:
                            shutil.copy2(src, dst)
                            print(f"  Copied (symlink not supported): {name}")
                    else:
                        shutil.copy2(src, dst)
                        print(f"  Copied: {name}")
            except Exception as e:
                print(f"  [ERROR] Failed to place {name}: {e}")
                success = False

        if success and not dry_run:
            print("\n[SUCCESS] All required definition files are now in place.")

    # --- Also ensure extra runtime files (Lovell_LIB, rotobs, scoring mats, etc.) ---
    extra_found = find_extra_files(search_roots)
    extra_missing = [name for name in EXPECTED_EXTRA_FILES if not any(f.name == name for f in extra_found)]

    if extra_found:
        print(f"\nFound {len(extra_found)} additional runtime file(s):")
        for f in extra_found:
            print(f"  + {f}")

    if extra_missing:
        print(f"\nMissing additional runtime files: {', '.join(extra_missing)}")

    if check_only:
        overall_success = (len(missing) == 0 and len(def_missing) == 0 and len(extra_missing) == 0)
        print(f"\n[check] {'READY' if overall_success else 'DATA MISSING'}")
        return overall_success, found + def_found + extra_found

    if extra_missing:
        for name in extra_missing:
            candidates = [f for f in extra_found if f.name == name]
            if not candidates:
                print(f"[ERROR] Cannot locate extra file: {name}")
                success = False
                continue

            src = candidates[0]
            dst = binary_base / name

            try:
                if dry_run:
                    action = "symlink" if use_links else "copy"
                    print(f"  [dry-run] Would {action}: {src} -> {dst}")
                else:
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    if use_links:
                        try:
                            dst.symlink_to(src.resolve())
                            print(f"  Linked: {name}")
                        except OSError:
                            shutil.copy2(src, dst)
                            print(f"  Copied (symlink not supported): {name}")
                    else:
                        shutil.copy2(src, dst)
                        print(f"  Copied: {name}")
            except Exception as e:
                print(f"  [ERROR] Failed to place {name}: {e}")
                success = False

        if success and not dry_run:
            print("\n[SUCCESS] All additional runtime files are now in place.")

    return success, found + def_found + extra_found


def print_definition_file_info(def_files: List[Path], verbose: bool = False) -> None:
    """Print useful diagnostic information about all runtime data files (matrices + defs + extras)."""
    print("\n=== FlexAIDδS Runtime Data Diagnostic ===")

    # Strong deduplication by filename (we only care about unique file names in the final report)
    seen_names = set()
    unique_files = []
    for f in def_files:
        if f.name not in seen_names:
            seen_names.add(f.name)
            unique_files.append(f)

    # Definition files
    amino_files = [f for f in unique_files if f.name.startswith("AMINO")]
    nucleotide_files = [f for f in unique_files if f.name.startswith("NUCLEOTIDES")]

    print(f"AMINO definition files found: {len(amino_files)}")
    for f in sorted(amino_files):
        print(f"  - {f.name}")

    print(f"NUCLEOTIDES definition files found: {len(nucleotide_files)}")
    for f in sorted(nucleotide_files):
        print(f"  - {f.name}")

    # Look for the preferred modern AMINO.def
    modern_amino = next((f for f in amino_files if f.name == "AMINO.def"), None)
    legacy_aminos = [f for f in amino_files if f.name != "AMINO.def"]

    if modern_amino:
        print("\n[OK] Modern AMINO.def (2011.12.08 recommended) is present.")
    elif legacy_aminos:
        print("\n[WARNING] Only legacy AMINO* variants found. Atom type numbers may be incompatible.")

    # Extra runtime files
    extra_files = [f for f in unique_files if f.name in EXPECTED_EXTRA_FILES]
    print(f"\nAdditional runtime files found: {len(extra_files)}")
    for f in sorted(extra_files, key=lambda x: x.name):
        print(f"  - {f.name}")

    missing_extra = [name for name in EXPECTED_EXTRA_FILES if not any(f.name == name for f in unique_files)]
    if missing_extra:
        print(f"  Missing: {', '.join(missing_extra)}")

    # FLEDIH summary
    print("\nSide-chain flexibility (FLEDIH) summary (from 2011 AMINO.def):")
    print("  ARG/LYS: 4 | GLN/GLU/MET: 3 | Many others: 1-2 | ALA/GLY/PRO: 0")

    print("\nAll these files (matrices + defs + extras) must be present in the binary base path for full functionality.")

    if verbose:
        print("\n[verbose] Full expected list in skill data/:")
        for name in EXPECTED_MATRICES + EXPECTED_DEF_FILES + EXPECTED_EXTRA_FILES:
            print(f"  - {name}")


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ensure critical FlexAIDδS runtime data (matrices + *.def + extra files like Lovell_LIB, rotobs, etc.) are available.",
        epilog="Part of the flexaid-docking skill. Run before any real docking task.",
    )
    parser.add_argument("--binary", "-b", type=Path, help="Path to FlexAIDδS binary")
    parser.add_argument("--source", "-s", dest="source", type=Path,
                        help="Path to a known-good installation to copy all runtime data from (matrices + defs + extras)")
    parser.add_argument("--link", action="store_true", help="Use symlinks instead of copies when possible")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Preview changes without modifying anything")
    parser.add_argument("--check", "--status", action="store_true", help="Only check, do not copy")
    parser.add_argument("--info", action="store_true", help="Print diagnostic information about found definition files (especially AMINO.def FLEDIH and variants)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show more details")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    print_skill_banner(verbose=args.verbose)

    success, found_defs = ensure_matrices(
        binary=args.binary,
        source=args.source,
        dry_run=args.dry_run,
        use_links=args.link,
        verbose=args.verbose,
        check_only=args.check,
    )

    if args.info:
        # Build search roots prioritizing the skill's own data/
        search_roots = list(DEFAULT_SEARCH_PATHS)
        if args.binary:
            base = get_binary_base_path(args.binary)
            search_roots = [base, base.parent] + search_roots
        if args.source:
            search_roots = [args.source] + search_roots

        all_found = find_def_files(search_roots) + find_extra_files(search_roots)
        # Deduplicate
        seen = set()
        unique = [f for f in all_found if not (f in seen or seen.add(f))]

        print_definition_file_info(unique, verbose=args.verbose)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())