#!/usr/bin/env python3
"""Validate Astex Diverse apo receptors for residual cognate ligand atoms.

Canonical tree:
  benchmarks/astex_diverse/astex_diverse/<PDB>/

Prior audit finding: most ``*_apo.pdb`` files are byte-identical to the deposit
``*.pdb`` in the same directory. Ligands for this set are often already extracted
to ``*_ligand.sdf`` (via CIF), so identity alone is not always a failure — this
script reports per-target whether the cognate ligand residue / coordinates still
appear in the apo structure.

Default mode is **report only**. Optional ``--fix`` rewrites apo by stripping
matching ligand residue names (with ``.bak`` backup). Fix is intentionally
conservative: skipped for peptide-like titles that collide with standard amino
acids unless ``--force-peptide``.

Usage:
  python3 scripts/validate_astex_apo_strip.py
  python3 scripts/validate_astex_apo_strip.py --report benchmarks/datasets/astex_apo_strip_report.csv
  python3 scripts/validate_astex_apo_strip.py --fix --targets 1SQ5,1HP0,1G9V

Exit codes:
  0 — report written; no fail-status targets (or --allow-fail)
  1 — one or more targets have ligand residue atoms remaining in apo
  2 — usage / path errors
"""
from __future__ import annotations

import argparse
import csv
import filecmp
import hashlib
import shutil
import sys
from collections import Counter
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TREE = REPO_ROOT / "benchmarks" / "astex_diverse" / "astex_diverse"
DEFAULT_REPORT = REPO_ROOT / "benchmarks" / "datasets" / "astex_apo_strip_report.csv"

# Standard polymer + solvent + ions (kept in apo intentionally when not the ligand).
STD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
    "A", "C", "G", "U", "DA", "DC", "DG", "DT", "I",
    "HOH", "WAT", "H2O", "DOD",
    "CA", "MG", "ZN", "MN", "FE", "NA", "CL", "K", "CU", "NI", "CO",
    "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "BME",
    "HEM", "HEC", "HEA", "HEB", "NAD", "NAP", "NDP", "FAD", "FMN",
    "ATP", "ADP", "AMP", "GTP", "GDP", "SAH", "SAM",
}

PEPTIDE_LIKE = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


@dataclass
class TargetReport:
    pdb_id: str
    ligand_title: str
    apo_path: str
    deposit_path: str
    ligand_path: str
    apo_exists: bool
    deposit_exists: bool
    ligand_exists: bool
    apo_identical_to_deposit: bool
    deposit_sha256: str
    apo_sha256: str
    ligand_atoms_in_apo: int
    ligand_residues_in_apo: int
    ligand_chains_in_apo: str
    nonstd_hetatm_residues: str
    nonstd_hetatm_atoms: int
    coord_match_atoms: int
    status: str
    notes: str
    fix_applied: bool = False
    fix_atoms_removed: int = 0


def _sha256(path: Path) -> str:
    if not path.is_file():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sdf_title(path: Path) -> str:
    try:
        with path.open(encoding="utf-8", errors="replace") as f:
            line = f.readline()
        return line.strip()
    except OSError:
        return ""


def parse_sdf_coords(path: Path) -> list[tuple[float, float, float]]:
    """Parse V2000 SDF atom coordinates (first molecule only)."""
    coords: list[tuple[float, float, float]] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return coords
    if len(lines) < 4:
        return coords
    try:
        n_atoms = int(lines[3][0:3])
    except ValueError:
        return coords
    for i in range(4, 4 + n_atoms):
        if i >= len(lines):
            break
        parts = lines[i].split()
        if len(parts) < 4:
            continue
        try:
            coords.append((float(parts[0]), float(parts[1]), float(parts[2])))
        except ValueError:
            continue
    return coords


def parse_pdb_atoms(path: Path) -> list[dict]:
    atoms: list[dict] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return atoms
    for lineno, line in enumerate(text.splitlines(), 1):
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        res = line[17:20].strip()
        chain = line[21:22] if len(line) > 21 else " "
        atoms.append(
            {
                "line": line,
                "lineno": lineno,
                "record": line[:6].strip(),
                "res": res,
                "chain": chain,
                "xyz": (x, y, z),
            }
        )
    return atoms


def count_coord_matches(
    pdb_atoms: list[dict],
    lig_coords: list[tuple[float, float, float]],
    tol: float = 0.35,
) -> int:
    """Count PDB atoms within tol Å of any ligand SDF atom (nearest-neighbor)."""
    if not lig_coords or not pdb_atoms:
        return 0
    tol2 = tol * tol
    n = 0
    for a in pdb_atoms:
        x, y, z = a["xyz"]
        for lx, ly, lz in lig_coords:
            dx = x - lx
            dy = y - ly
            dz = z - lz
            if dx * dx + dy * dy + dz * dz <= tol2:
                n += 1
                break
    return n


def analyze_target(
    pdb_id: str,
    tree: Path,
    coord_tol: float,
) -> TargetReport:
    d = tree / pdb_id
    apo = d / f"{pdb_id}_apo.pdb"
    deposit = d / f"{pdb_id}.pdb"
    ligand = d / f"{pdb_id}_ligand.sdf"

    title = sdf_title(ligand) if ligand.is_file() else ""
    apo_atoms = parse_pdb_atoms(apo) if apo.is_file() else []
    lig_coords = parse_sdf_coords(ligand) if ligand.is_file() else []

    peptide_like = bool(title) and title in PEPTIDE_LIKE
    # Residue-name match is reliable only for non-polymer 3-letter ligand codes.
    # Peptide-like titles (ALA/TYR/…) collide with the protein sequence — use
    # coordinate matching exclusively for those cases.
    if peptide_like:
        lig_atoms: list[dict] = []
    else:
        lig_atoms = [a for a in apo_atoms if a["res"] == title and title]
    chains = sorted({a["chain"] for a in lig_atoms})
    res_keys = {(a["res"], a["chain"], a["line"][22:26].strip()) for a in lig_atoms}

    nonstd = Counter(
        a["res"]
        for a in apo_atoms
        if a["record"] == "HETATM" and a["res"] not in STD_RESIDUES
    )
    # Exclude the cognate ligand title from "other nonstd" noise when matched.
    nonstd_other = Counter(
        {r: n for r, n in nonstd.items() if r != title or not title}
    )

    identical = False
    if apo.is_file() and deposit.is_file():
        identical = filecmp.cmp(apo, deposit, shallow=False)

    coord_hits = count_coord_matches(apo_atoms, lig_coords, tol=coord_tol)

    notes: list[str] = []
    if not apo.is_file():
        notes.append("missing_apo")
    if not deposit.is_file():
        notes.append("missing_deposit")
    if not ligand.is_file():
        notes.append("missing_ligand_sdf")
    if identical:
        notes.append("apo_byte_identical_to_deposit")
    if peptide_like:
        notes.append("peptide_like_ligand_title")
    if coord_hits > 0 and len(lig_atoms) == 0:
        notes.append("coord_match_without_resname")

    n_lig = len(lig_atoms)
    # Status:
    #  fail — cognate non-peptide ligand residue still present, or strong coord match
    #  warn — weak coord matches, missing files, or peptide ambiguity without coords
    #  ok   — no residual ligand evidence (identity to deposit alone is not a fail when
    #         the cognate ligand lives only in CIF/SDF and is already absent from PDB)
    if not apo.is_file() or not ligand.is_file():
        status = "warn"
    elif n_lig > 0:
        status = "fail"
    elif coord_hits >= max(3, max(1, len(lig_coords) // 2)) and lig_coords:
        # Strong coordinate evidence ligand remains under another name / peptide
        status = "fail"
        notes.append("strong_coord_match")
    elif coord_hits > 0:
        status = "warn"
    elif peptide_like and not lig_coords:
        status = "warn"
        notes.append("peptide_needs_manual_review")
    elif identical:
        status = "ok"
        notes.append("identity_ok_if_ligand_absent_from_pdb")
    else:
        status = "ok"

    return TargetReport(
        pdb_id=pdb_id,
        ligand_title=title,
        apo_path=str(apo.relative_to(REPO_ROOT)) if apo.is_file() else str(apo),
        deposit_path=str(deposit.relative_to(REPO_ROOT)) if deposit.is_file() else str(deposit),
        ligand_path=str(ligand.relative_to(REPO_ROOT)) if ligand.is_file() else str(ligand),
        apo_exists=apo.is_file(),
        deposit_exists=deposit.is_file(),
        ligand_exists=ligand.is_file(),
        apo_identical_to_deposit=identical,
        deposit_sha256=_sha256(deposit),
        apo_sha256=_sha256(apo),
        ligand_atoms_in_apo=n_lig,
        ligand_residues_in_apo=len(res_keys),
        ligand_chains_in_apo=",".join(chains),
        nonstd_hetatm_residues=";".join(f"{r}:{n}" for r, n in sorted(nonstd_other.items())),
        nonstd_hetatm_atoms=sum(nonstd_other.values()),
        coord_match_atoms=coord_hits,
        status=status,
        notes=";".join(notes),
    )


def strip_ligand_residue(
    apo: Path,
    residue_name: str,
    backup: bool = True,
) -> int:
    """Remove ATOM/HETATM lines whose residue name matches residue_name.

    Returns number of atom lines removed. Writes ``apo`` in place; optional
    ``.bak`` beside the original.
    """
    if not residue_name:
        return 0
    text = apo.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines(keepends=True)
    kept: list[str] = []
    removed = 0
    for line in lines:
        if (line.startswith("ATOM  ") or line.startswith("HETATM")) and len(line) >= 20:
            res = line[17:20].strip()
            if res == residue_name:
                removed += 1
                continue
        kept.append(line)
    if removed == 0:
        return 0
    if backup:
        bak = apo.with_suffix(apo.suffix + ".bak")
        if not bak.exists():
            shutil.copy2(apo, bak)
    apo.write_text("".join(kept), encoding="utf-8")
    return removed


def discover_targets(tree: Path) -> list[str]:
    if not tree.is_dir():
        return []
    out = []
    for p in sorted(tree.iterdir()):
        if p.is_dir() and (p / f"{p.name}_apo.pdb").exists():
            out.append(p.name)
        elif p.is_dir() and (p / f"{p.name}.pdb").exists():
            out.append(p.name)
    return out


def write_csv(path: Path, rows: Iterable[TargetReport]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fieldnames = list(asdict(rows[0]).keys()) if rows else [f.name for f in TargetReport.__dataclass_fields__.values()]  # type: ignore[attr-defined]
    # dataclass fields order
    fieldnames = list(TargetReport.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--tree",
        type=Path,
        default=DEFAULT_TREE,
        help=f"Canonical Astex tree (default: {DEFAULT_TREE.relative_to(REPO_ROOT)})",
    )
    ap.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_REPORT,
        help=f"Output CSV path (default: {DEFAULT_REPORT.relative_to(REPO_ROOT)})",
    )
    ap.add_argument(
        "--targets",
        type=str,
        default="",
        help="Comma-separated PDB IDs (default: all under tree)",
    )
    ap.add_argument(
        "--coord-tol",
        type=float,
        default=0.35,
        help="Å tolerance for SDF↔PDB coordinate ligand match (default 0.35)",
    )
    ap.add_argument(
        "--fix",
        action="store_true",
        help="Strip cognate ligand residue name from apo (writes .bak; conservative)",
    )
    ap.add_argument(
        "--force-peptide",
        action="store_true",
        help="Allow --fix when ligand title is a standard amino-acid code",
    )
    ap.add_argument(
        "--allow-fail",
        action="store_true",
        help="Exit 0 even if some targets have status=fail",
    )
    ap.add_argument(
        "--summary-only",
        action="store_true",
        help="Print summary counts only (still writes CSV)",
    )
    args = ap.parse_args(argv)

    tree = args.tree if args.tree.is_absolute() else REPO_ROOT / args.tree
    report_path = args.report if args.report.is_absolute() else REPO_ROOT / args.report

    if not tree.is_dir():
        print(f"ERROR: canonical tree not found: {tree}", file=sys.stderr)
        return 2

    if args.targets.strip():
        targets = [t.strip().upper() for t in args.targets.split(",") if t.strip()]
    else:
        targets = discover_targets(tree)

    if not targets:
        print(f"ERROR: no targets under {tree}", file=sys.stderr)
        return 2

    rows: list[TargetReport] = []
    for pdb_id in targets:
        rep = analyze_target(pdb_id, tree, coord_tol=args.coord_tol)
        if args.fix and rep.status == "fail" and rep.ligand_title:
            if rep.ligand_title in PEPTIDE_LIKE and not args.force_peptide:
                rep.notes = (rep.notes + ";" if rep.notes else "") + "fix_skipped_peptide_like"
            else:
                apo = tree / pdb_id / f"{pdb_id}_apo.pdb"
                removed = strip_ligand_residue(apo, rep.ligand_title, backup=True)
                rep.fix_applied = removed > 0
                rep.fix_atoms_removed = removed
                # Re-analyze after fix
                rep = analyze_target(pdb_id, tree, coord_tol=args.coord_tol)
                rep.fix_applied = removed > 0
                rep.fix_atoms_removed = removed
                if removed:
                    rep.notes = (rep.notes + ";" if rep.notes else "") + f"stripped_{removed}_atoms"
        rows.append(rep)

    write_csv(report_path, rows)

    n = len(rows)
    n_fail = sum(1 for r in rows if r.status == "fail")
    n_warn = sum(1 for r in rows if r.status == "warn")
    n_ok = sum(1 for r in rows if r.status == "ok")
    n_ident = sum(1 for r in rows if r.apo_identical_to_deposit)
    n_lig = sum(1 for r in rows if r.ligand_atoms_in_apo > 0)
    n_fixed = sum(1 for r in rows if r.fix_applied)

    print(f"Astex apo strip validation — {n} targets")
    print(f"  tree:     {tree}")
    print(f"  report:   {report_path}")
    print(f"  status:   ok={n_ok}  warn={n_warn}  fail={n_fail}")
    print(f"  identical apo==deposit: {n_ident}/{n}")
    print(f"  ligand residue atoms remaining: {n_lig}/{n}")
    if args.fix:
        print(f"  fix applied: {n_fixed}")

    if not args.summary_only:
        # Show fails + a few warns
        shown = 0
        for r in rows:
            if r.status in ("fail", "warn") and shown < 20:
                print(
                    f"  [{r.status}] {r.pdb_id} lig={r.ligand_title!r} "
                    f"atoms={r.ligand_atoms_in_apo} coord={r.coord_match_atoms} "
                    f"ident={r.apo_identical_to_deposit} notes={r.notes}"
                )
                shown += 1
        if n_fail + n_warn > shown:
            print(f"  ... ({n_fail + n_warn - shown} more non-ok in CSV)")

    if n_fail and not args.allow_fail:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
