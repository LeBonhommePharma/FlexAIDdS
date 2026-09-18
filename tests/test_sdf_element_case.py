"""Ligand SDF atom symbols must be canonical MDL case ("Cl", never "CL").

WHY THIS TEST EXISTS
--------------------
`DatasetRunner::extract_ligand` builds an MDL V2000 atom block from PDB/mmCIF
records.  The PDB element column (cols 77-78) is UPPERCASE by PDB spec, while
MDL V2000 atom symbols are CASE-SENSITIVE.  Copying the PDB casing verbatim
emitted "CL" / "BR", and the downstream chain failed silently:

    SDF "CL"
      -> ProcessLigand cannot perceive the symbol
      -> atom typed 39 (DUMMY), named "*"
      -> vcfunction.cpp classifies "Du" as NON-HEAVY
      -> the halogen is dropped from contact scoring entirely
      -> get_element(39) writes "Du" into every pose, which RDKit refuses
      -> validate_ligand_integrity.py counts one fewer heavy atom
         and still returns ok: True

Measured on Astex-85 (2026-09-17): 14/85 targets carry CL or BR, and the same
14 carried a Du atom in their elected pose.  Patched vs unpatched extractor on
identical inputs: all-caps symbols 14 -> 0, element identity preserved 85/85,
atom counts unchanged, and ProcessLigand types the halogen 39 -> 24.

The first test is the source contract.  The second is the data invariant that
would have caught this before a campaign ran -- point it at a prepared cache
via FLEXAIDDS_TEST_CACHE.
"""

from __future__ import annotations

import os
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
RUNNER = REPO / "LIB" / "DatasetRunner.cpp"

# Two-letter element symbols that are legal in a ligand and whose all-caps form
# is what the PDB column hands us.  A one-letter symbol (C/N/O/S/P/F/I) has no
# case ambiguity, which is why only these ever broke.
TWO_LETTER = ("Cl", "Br", "Se", "Mg", "Mn", "Fe", "Zn", "Ca", "Cu", "Ni", "Hg", "Cd", "Sr", "Co")


def _atom_symbols(sdf_path: pathlib.Path) -> list[str]:
    """Element symbols from a V2000 atom block, in file order."""
    lines = sdf_path.read_text(errors="replace").splitlines()
    if len(lines) < 4:
        return []
    try:
        natoms = int(lines[3][0:3])
    except ValueError:
        return []
    out: list[str] = []
    for ln in lines[4 : 4 + natoms]:
        parts = ln.split()
        if len(parts) >= 4:
            out.append(parts[3])
    return out


def test_extractor_canonicalises_the_element_symbol() -> None:
    """The writer must normalise case, not copy the PDB column verbatim."""
    src = RUNNER.read_text(encoding="utf-8")

    # The atom block is emitted here; the canonicalisation must sit before it.
    assert "std::string elem = atom.element;" in src, (
        "the SDF atom-block writer moved -- re-anchor this test"
    )
    assert "MDL V2000 atom symbols are CASE-SENSITIVE" in src, (
        "the rationale comment is gone; a future reader will re-introduce the bug"
    )
    # first char upper, remainder lower -- both halves required
    assert re.search(r"elem\[0\]\s*=\s*static_cast<char>\(std::toupper", src), (
        "first-character upper-casing missing from the SDF element write"
    )
    assert re.search(r"elem\[ci\]\s*=\s*static_cast<char>\(std::tolower", src), (
        "remainder lower-casing missing -- 'CL' would still reach the atom block"
    )


def test_the_guard_is_not_vacuous() -> None:
    """Deleting the fix must fail the contract test above, not pass it."""
    src = RUNNER.read_text(encoding="utf-8")
    without_fix = re.sub(
        r"\s*elem\[ci\]\s*=\s*static_cast<char>\(std::tolower[^\n]*\n", "\n", src
    )
    assert without_fix != src, "injection did not modify the source -- test is vacuous"
    assert not re.search(r"elem\[ci\]\s*=\s*static_cast<char>\(std::tolower", without_fix), (
        "the assertion used by the contract test still passes on a tree with the fix removed"
    )


def test_prepared_cache_has_no_uppercase_two_letter_symbols() -> None:
    """Data invariant: no prepared ligand SDF may carry an all-caps symbol.

    This is the gate that would have caught the defect before a campaign ran.
    Skipped unless a prepared cache is supplied.
    """
    cache = os.environ.get("FLEXAIDDS_TEST_CACHE")
    if not cache:
        pytest.skip("set FLEXAIDDS_TEST_CACHE to a prepared dataset cache to run this")
    root = pathlib.Path(cache)
    sdfs = sorted(root.glob("*/*/*_ligand.sdf")) or sorted(root.glob("*/*_ligand.sdf"))
    assert sdfs, f"no *_ligand.sdf under {root} -- nothing measured, so nothing proven"

    seen_two_letter = 0
    offenders: list[tuple[str, str]] = []
    for sdf in sdfs:
        for sym in _atom_symbols(sdf):
            if len(sym) == 2:
                seen_two_letter += 1
                if sym.isupper():
                    offenders.append((sdf.parent.name, sym))

    # Non-vacuity: a corpus with no two-letter symbols at all cannot detect the bug.
    assert seen_two_letter > 0, (
        f"{len(sdfs)} SDFs carry no two-letter element symbol; this corpus cannot "
        "witness the defect -- use one containing a halogenated ligand"
    )
    assert not offenders, (
        f"{len(set(t for t, _ in offenders))} target(s) emit an all-caps element symbol, "
        f"which ProcessLigand types as DUMMY: {sorted(set(offenders))[:12]}"
    )
