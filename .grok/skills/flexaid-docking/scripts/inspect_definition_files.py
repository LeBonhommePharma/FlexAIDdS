#!/usr/bin/env python3
"""
inspect_definition_files.py

Part of the flexaid-docking skill.

A small standalone helper to inspect and report on the critical definition files
(AMINO*.def and NUCLEOTIDES*.def) that the FlexAIDδS binary requires.

Useful for:
- Verifying which version of the definition files you have
- Understanding side-chain flexibility (FLEDIH) that will be sampled
- Diagnosing atom-typing or flexibility problems before a run
- Checking consistency between your data and the bundled skill data

Example usage:
    python3 .grok/skills/flexaid-docking/scripts/inspect_definition_files.py
    python3 .grok/skills/flexaid-docking/scripts/inspect_definition_files.py --binary /path/to/FlexAIDdS
    python3 .grok/skills/flexaid-docking/scripts/inspect_definition_files.py --source /path/to/good/install
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

# Same lists as ensure_docking_data.py for consistency
EXPECTED_DEF_FILES = [
    "AMINO.def", "AMINO8.def", "AMINO12.def", "AMINO26.def",
    "NUCLEOTIDES.def", "NUCLEOTIDES8.def", "NUCLEOTIDES12.def", "NUCLEOTIDES26.def",
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


def get_binary_base_path(binary: Optional[Path]) -> Path:
    if binary and binary.exists():
        return binary.parent.resolve()
    for c in [Path.cwd(), Path(__file__).resolve().parents[4] / "build"]:
        if (c / "FlexAIDδS").exists():
            return c.resolve()
    return Path.cwd()


def print_diagnostics(def_files: List[Path], verbose: bool = False) -> None:
    print("=== FlexAIDδS Definition Files Inspector ===")
    print()

    amino = sorted([f for f in def_files if f.name.startswith("AMINO")])
    nucl = sorted([f for f in def_files if f.name.startswith("NUCLEOTIDES")])

    print(f"AMINO*.def files found: {len(amino)}")
    for f in amino:
        print(f"  + {f}")

    print(f"\nNUCLEOTIDES*.def files found: {len(nucl)}")
    for f in nucl:
        print(f"  + {f}")

    modern = next((f for f in amino if f.name == "AMINO.def"), None)
    legacy = [f for f in amino if f.name != "AMINO.def"]

    print()
    if modern:
        print("[GOOD] Modern AMINO.def (2011.12.08) detected — recommended for current matrices.")
    elif legacy:
        print("[WARN] Only legacy AMINO* variants present. Atom type numbers differ from modern matrices.")
        print("       This frequently causes incorrect typing or scoring.")

    print("\n--- Side-chain Flexibility (FLEDIH dihedrals) ---")
    print("These control which torsions the GA will actually sample:")
    for res, count in sorted(FLEDIH_SUMMARY.items(), key=lambda x: -x[1]):
        print(f"  {res:3s}: {count} rotatable dihedral(s)")

    print("\nTip: Run with --verbose for full file paths and search roots.")
    if verbose:
        print("\nSearch roots used:")
        for r in DEFAULT_SEARCH_PATHS:
            print(f"  {r}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Inspect FlexAIDδS definition files (AMINO*.def / NUCLEOTIDES*.def) and report flexibility + compatibility info.",
        epilog="Part of the flexaid-docking skill. Run this before important docking jobs."
    )
    parser.add_argument("--binary", "-b", type=Path, help="Path to FlexAIDδS binary (helps locate data)")
    parser.add_argument("--source", "-s", type=Path, help="Explicit directory to search for definition files")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    search_roots = list(DEFAULT_SEARCH_PATHS)
    if args.binary:
        base = get_binary_base_path(args.binary)
        search_roots = [base, base.parent] + search_roots
    if args.source:
        search_roots = [args.source] + search_roots

    found = find_def_files(search_roots)
    print_diagnostics(found, verbose=args.verbose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
