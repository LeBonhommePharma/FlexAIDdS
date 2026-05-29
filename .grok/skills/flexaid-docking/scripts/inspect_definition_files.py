#!/usr/bin/env python3
"""
inspect_definition_files.py

Part of the flexaid-docking skill.

A small standalone helper to inspect and report on the critical runtime data files
that the FlexAIDδS binary requires in its base path:

- Interaction matrices (MC_*.dat)
- Definition files (AMINO*.def / NUCLEOTIDES*.def)
- Extra files (Lovell_LIB.dat, rotobs.lst, SYBYL_emat.dat, scoring support, etc.)

Useful for:
- Verifying completeness of your runtime data pack
- Understanding side-chain flexibility (FLEDIH) that will be sampled
- Detecting legacy vs modern definition variants
- Diagnosing "missing file" problems before expensive docking runs

Example usage:
    python3 .grok/skills/flexaid-docking/scripts/inspect_definition_files.py
    python3 .grok/skills/flexaid-docking/scripts/inspect_definition_files.py --binary /path/to/FlexAIDdS
    python3 .grok/skills/flexaid-docking/scripts/inspect_definition_files.py --source /path/to/good/install
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional

# Same lists as ensure_docking_data.py for consistency
EXPECTED_MATRICES = ["MC_st0r5.2_6.dat"]
EXPECTED_DEF_FILES = [
    "AMINO.def", "AMINO8.def", "AMINO12.def", "AMINO26.def",
    "NUCLEOTIDES.def", "NUCLEOTIDES8.def", "NUCLEOTIDES12.def", "NUCLEOTIDES26.def",
]
EXPECTED_EXTRA_FILES = [
    "Lovell_LIB.dat", "rotobs.lst", "SYBYL_emat.dat",
    "M6_cons_3.dat",
    "nrg_mat_BEST_011912.dat", "nrg_mat_BEST_012012.dat",
    "scr_bin.dat", "scr_mat.dat",
]

# Hardcoded high-value summary derived from the authoritative 2011 AMINO.def
FLEDIH_SUMMARY = {
    "ARG": 4, "LYS": 4,
    "GLN": 3, "GLU": 3, "MET": 3,
    "ASN": 2, "ASP": 2, "HIS": 2, "ILE": 2, "LEU": 2,
    "PHE": 2, "TRP": 2, "TYR": 2,
    "CYS": 1, "SER": 1, "THR": 1, "VAL": 1,
    "ALA": 0, "GLY": 0, "PRO": 0,
}

DEFAULT_SEARCH_PATHS = [
    Path(__file__).resolve().parent.parent / "data",
    Path.home() / ".flexaidds" / "data",
    Path("/usr/local/share/flexaidds/data"),
    Path("/opt/flexaidds/data"),
    Path.home() / "flexaidds-data",
    Path(__file__).resolve().parents[4] / "build",
    Path(__file__).resolve().parents[4] / "WRK",
]


def find_def_files(search_roots: List[Path]) -> List[Path]:
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
    for c in [Path.cwd(), Path(__file__).resolve().parents[4] / "build"]:
        if (c / "FlexAIDδS").exists():
            return c.resolve()
    return Path.cwd()


def should_use_light_mode(args) -> bool:
    """Automatically decide lightweight behavior (same logic as ensure_docking_data.py)."""
    if getattr(args, "quick", False):
        return True

    ci_env_vars = ("CI", "GITHUB_ACTIONS", "GITLAB_CI", "TRAVIS", "CIRCLECI", "JENKINS_URL")
    if any(os.environ.get(var) for var in ci_env_vars):
        return True

    if getattr(args, "info", False):   # if the user forced rich mode
        return False

    return False


def print_diagnostics(all_files: List[Path], verbose: bool = False) -> None:
    print("=== FlexAIDδS Runtime Data Inspector ===")
    print()

    # Strong deduplication by filename (clean report, no duplicate filenames)
    seen_names = set()
    unique_files = []
    for f in all_files:
        if f.name not in seen_names:
            seen_names.add(f.name)
            unique_files.append(f)

    matrices = [f for f in unique_files if f.name in EXPECTED_MATRICES]
    defs = [f for f in unique_files if f.name in EXPECTED_DEF_FILES]
    extras = [f for f in unique_files if f.name in EXPECTED_EXTRA_FILES]

    print(f"Matrices found: {len(matrices)}")
    for f in sorted(matrices, key=lambda x: x.name):
        print(f"  + {f.name}")

    amino = sorted([f for f in defs if f.name.startswith("AMINO")])
    nucl = sorted([f for f in defs if f.name.startswith("NUCLEOTIDES")])

    print(f"\nAMINO*.def files found: {len(amino)}")
    for f in amino:
        print(f"  + {f.name}")

    print(f"NUCLEOTIDES*.def files found: {len(nucl)}")
    for f in nucl:
        print(f"  + {f.name}")

    print(f"\nAdditional runtime files found: {len(extras)}")
    for f in sorted(extras, key=lambda x: x.name):
        print(f"  + {f.name}")

    modern = next((f for f in defs if f.name == "AMINO.def"), None)
    legacy = [f for f in defs if f.name.startswith("AMINO") and f.name != "AMINO.def"]

    print()
    if modern:
        print("[GOOD] Modern AMINO.def (2011) present.")
    elif legacy:
        print("[WARN] Only legacy AMINO variants — potential atom type mismatch.")

    print("\n--- Side-chain Flexibility (FLEDIH) ---")
    print("Controls which torsions the GA samples (from 2011 AMINO.def):")
    for res, count in sorted(FLEDIH_SUMMARY.items(), key=lambda x: -x[1]):
        print(f"  {res:3s}: {count} rotatable dihedral(s)")

    print("\nPerformance note: By default, a comprehensive inspection of all bundled runtime data files is performed. In CI pipelines and resource-limited environments, validation is automatically restricted to the minimum required set to reduce I/O and memory overhead.")

    print("\nTip: Use --verbose for search paths. Also available via ensure_docking_data.py --info")
    if verbose:
        print("\nSearch roots used:")
        for r in DEFAULT_SEARCH_PATHS:
            print(f"  {r}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Inspect all critical FlexAIDδS runtime data files (matrices + *.def + extras like Lovell_LIB, rotobs, etc.) and report on completeness + flexibility info.",
        epilog="Part of the flexaid-docking skill. Run this before important docking jobs."
    )
    parser.add_argument("--binary", "-b", type=Path, help="Path to FlexAIDδS binary (helps locate data)")
    parser.add_argument("--source", "-s", type=Path, help="Explicit directory to search for definition files")
    parser.add_argument("--quick", action="store_true", help="Force lightweight mode (normally auto-selected in CI / low-resource environments)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    search_roots = list(DEFAULT_SEARCH_PATHS)
    if args.binary:
        base = get_binary_base_path(args.binary)
        search_roots = [base, base.parent] + search_roots
    if args.source:
        search_roots = [args.source] + search_roots

    use_light = should_use_light_mode(args)

    if use_light:
        critical = ["MC_st0r5.2_6.dat", "AMINO.def"]
        all_found = []
        for root in search_roots:
            if root and root.exists():
                for name in critical:
                    p = root / name
                    if p.is_file():
                        all_found.append(p)
        if getattr(args, "quick", False) or not getattr(args, "info", False):
            print("Lightweight validation automatically enabled (CI or resource-limited environment detected).")
    else:
        # Single clean finder for all categories (rich --info style by default)
        all_found = []
        for root in search_roots:
            if root and root.exists():
                for name in EXPECTED_MATRICES + EXPECTED_DEF_FILES + EXPECTED_EXTRA_FILES:
                    p = root / name
                    if p.is_file():
                        all_found.append(p)

    # Deduplicate
    seen = set()
    unique = []
    for p in all_found:
        try:
            key = p.resolve()
        except Exception:
            key = p
        if key not in seen:
            seen.add(key)
            unique.append(p)

    print_diagnostics(unique, verbose=args.verbose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
