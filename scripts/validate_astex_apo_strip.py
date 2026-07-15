#!/usr/bin/env python3
"""Validate Astex Diverse apo receptors for residual cognate ligand atoms.

Canonical tree:
  benchmarks/astex_diverse/astex_diverse/<PDB>/

Prior audit finding: most ``*_apo.pdb`` files are byte-identical to the deposit
``*.pdb`` in the same directory. Ligands for this set are often already extracted
to ``*_ligand.sdf`` (via CIF), so identity alone is not always a failure — this
script reports per-target whether the cognate ligand residue / coordinates still
appear in the apo structure.

Default mode is **report only** (CSV + summary JSON). Optional ``--fix-dry-run``
describes strip operations without writing. Real strip requires ``--write`` and is
capped at 3 pilot targets unless ``--all-safe`` (and all strip candidates are
non-peptide, residue-name-safe).

Usage:
  python3 scripts/validate_astex_apo_strip.py
  python3 scripts/validate_astex_apo_strip.py --strict
  python3 scripts/validate_astex_apo_strip.py --fix-dry-run
  python3 scripts/validate_astex_apo_strip.py --write --targets 1SQ5,1HP0,1G9V
  python3 scripts/validate_astex_apo_strip.py \\
      --report benchmarks/datasets/astex_apo_strip_report.csv \\
      --summary-json benchmarks/datasets/astex_apo_strip_summary.json

Exit codes:
  0 — report written; no fail-status targets (or --allow-fail)
  1 — one or more targets have status=fail (default and under --strict)
  2 — usage / path errors
"""
from __future__ import annotations

import argparse
import csv
import filecmp
import hashlib
import json
import shutil
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TREE = REPO_ROOT / "benchmarks" / "astex_diverse" / "astex_diverse"
DEFAULT_REPORT = REPO_ROOT / "benchmarks" / "datasets" / "astex_apo_strip_report.csv"
DEFAULT_SUMMARY = REPO_ROOT / "benchmarks" / "datasets" / "astex_apo_strip_summary.json"

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

# Hard cap for real strip unless --all-safe.
DEFAULT_WRITE_PILOT_LIMIT = 3


@dataclass
class TargetReport:
    """Per-target apo-strip validation row.

    Field names keep both the mission contract names and prior CSV columns
    so existing consumers keep working.
    """

    pdb_id: str
    ligand_title: str
    ligand_source: str  # sdf | mol2 | none
    apo_path: str
    deposit_path: str
    ligand_path: str
    apo_exists: bool
    deposit_exists: bool
    ligand_exists: bool
    # Mission alias (byte-identical apo vs deposit PDB).
    identical_sha_to_deposit: bool
    apo_identical_to_deposit: bool
    deposit_sha256: str
    apo_sha256: str
    ligand_atoms_in_apo: int
    ligand_residues_in_apo: int
    ligand_chains_in_apo: str
    # Residue names of cognate ligand HETATM/ATOM hits in apo (semicolon list).
    ligand_hetatm_names: str
    nonstd_hetatm_residues: str
    nonstd_hetatm_atoms: int
    coord_match_atoms: int
    # pass|fail|warn  (pass == legacy ok)
    status: str
    pass_fail: str  # "pass" or "fail" (warn maps to pass for binary gate)
    notes: str
    fix_planned: bool = False
    fix_atoms_planned: int = 0
    fix_applied: bool = False
    fix_atoms_removed: int = 0


def _sha256_file(path: Path) -> str:
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


def mol2_title_and_coords(path: Path) -> tuple[str, list[tuple[float, float, float]]]:
    """Parse first @<TRIPOS>MOLECULE name and @<TRIPOS>ATOM coords from MOL2."""
    title = ""
    coords: list[tuple[float, float, float]] = []
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return title, coords
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if line.upper().startswith("@<TRIPOS>MOLECULE"):
            if i + 1 < n:
                title = lines[i + 1].strip()
            i += 1
            continue
        if line.upper().startswith("@<TRIPOS>ATOM"):
            i += 1
            while i < n and not lines[i].strip().startswith("@<TRIPOS>"):
                parts = lines[i].split()
                # atom_id name x y z type ...
                if len(parts) >= 5:
                    try:
                        coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                    except ValueError:
                        pass
                i += 1
            continue
        i += 1
    return title, coords


def resolve_ligand(d: Path, pdb_id: str) -> tuple[Path | None, str, list[tuple[float, float, float]], str]:
    """Return (path, title, coords, source) preferring SDF then MOL2."""
    sdf = d / f"{pdb_id}_ligand.sdf"
    mol2 = d / f"{pdb_id}_ligand.mol2"
    if sdf.is_file():
        return sdf, sdf_title(sdf), parse_sdf_coords(sdf), "sdf"
    if mol2.is_file():
        title, coords = mol2_title_and_coords(mol2)
        return mol2, title, coords, "mol2"
    for p in sorted(d.glob(f"{pdb_id}*.sdf")):
        return p, sdf_title(p), parse_sdf_coords(p), "sdf"
    for p in sorted(d.glob(f"{pdb_id}*.mol2")):
        title, coords = mol2_title_and_coords(p)
        return p, title, coords, "mol2"
    return None, "", [], "none"


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
) -> tuple[int, list[dict]]:
    """Count PDB atoms within tol Å of any ligand SDF/MOL2 atom.

    Returns (count, matching atom dicts).
    """
    if not lig_coords or not pdb_atoms:
        return 0, []
    tol2 = tol * tol
    hits: list[dict] = []
    for a in pdb_atoms:
        x, y, z = a["xyz"]
        for lx, ly, lz in lig_coords:
            dx = x - lx
            dy = y - ly
            dz = z - lz
            if dx * dx + dy * dy + dz * dz <= tol2:
                hits.append(a)
                break
    return len(hits), hits


def _rel(path: Path | None) -> str:
    if path is None:
        return ""
    try:
        if path.is_file():
            return str(path.relative_to(REPO_ROOT))
    except ValueError:
        pass
    return str(path)


def analyze_target(
    pdb_id: str,
    tree: Path,
    coord_tol: float,
) -> TargetReport:
    d = tree / pdb_id
    apo = d / f"{pdb_id}_apo.pdb"
    deposit = d / f"{pdb_id}.pdb"
    ligand_path, title, lig_coords, lig_source = resolve_ligand(d, pdb_id)

    apo_atoms = parse_pdb_atoms(apo) if apo.is_file() else []

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
    het_names = sorted({a["res"] for a in lig_atoms})

    nonstd = Counter(
        a["res"]
        for a in apo_atoms
        if a["record"] == "HETATM" and a["res"] not in STD_RESIDUES
    )
    nonstd_other = Counter({r: n for r, n in nonstd.items() if r != title or not title})

    identical = False
    if apo.is_file() and deposit.is_file():
        identical = filecmp.cmp(apo, deposit, shallow=False)

    coord_hits, coord_atoms = count_coord_matches(apo_atoms, lig_coords, tol=coord_tol)
    if coord_hits and not het_names:
        het_names = sorted({a["res"] for a in coord_atoms})

    notes: list[str] = []
    if not apo.is_file():
        notes.append("missing_apo")
    if not deposit.is_file():
        notes.append("missing_deposit")
    if ligand_path is None:
        notes.append("missing_ligand_sdf_mol2")
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
    #  pass — no residual ligand evidence (identity to deposit alone is not a fail when
    #         the cognate ligand lives only in CIF/SDF and is already absent from PDB)
    if not apo.is_file() or ligand_path is None:
        status = "warn"
    elif n_lig > 0:
        status = "fail"
    elif coord_hits >= max(3, max(1, len(lig_coords) // 2)) and lig_coords:
        status = "fail"
        notes.append("strong_coord_match")
    elif coord_hits > 0:
        status = "warn"
    elif peptide_like and not lig_coords:
        status = "warn"
        notes.append("peptide_needs_manual_review")
    elif identical:
        status = "pass"
        notes.append("identity_ok_if_ligand_absent_from_pdb")
    else:
        status = "pass"

    pass_fail = "fail" if status == "fail" else "pass"

    return TargetReport(
        pdb_id=pdb_id,
        ligand_title=title,
        ligand_source=lig_source,
        apo_path=_rel(apo) if apo.is_file() else str(apo),
        deposit_path=_rel(deposit) if deposit.is_file() else str(deposit),
        ligand_path=_rel(ligand_path) if ligand_path and ligand_path.is_file() else (
            str(ligand_path) if ligand_path else ""
        ),
        apo_exists=apo.is_file(),
        deposit_exists=deposit.is_file(),
        ligand_exists=ligand_path is not None and ligand_path.is_file(),
        identical_sha_to_deposit=identical,
        apo_identical_to_deposit=identical,
        deposit_sha256=_sha256_file(deposit),
        apo_sha256=_sha256_file(apo),
        ligand_atoms_in_apo=n_lig,
        ligand_residues_in_apo=len(res_keys),
        ligand_chains_in_apo=",".join(chains),
        ligand_hetatm_names=";".join(het_names),
        nonstd_hetatm_residues=";".join(f"{r}:{n}" for r, n in sorted(nonstd_other.items())),
        nonstd_hetatm_atoms=sum(nonstd_other.values()),
        coord_match_atoms=coord_hits,
        status=status,
        pass_fail=pass_fail,
        notes=";".join(notes),
    )


def strip_ligand_residue(
    apo: Path,
    residue_name: str,
    backup: bool = True,
    dry_run: bool = False,
) -> int:
    """Remove ATOM/HETATM lines whose residue name matches residue_name.

    Returns number of atom lines that would be / were removed.
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
    if removed == 0 or dry_run:
        return removed
    if backup:
        bak = apo.with_suffix(apo.suffix + ".bak")
        if not bak.exists():
            shutil.copy2(apo, bak)
    apo.write_text("".join(kept), encoding="utf-8")
    return removed


def plan_strip(rep: TargetReport) -> tuple[bool, int, str]:
    """Decide whether a conservative residue-name strip is planned.

    Returns (planned, estimated_atoms, reason).
    """
    if rep.status != "fail":
        return False, 0, "not_fail"
    if not rep.ligand_title:
        return False, 0, "no_ligand_title"
    if rep.ligand_title in PEPTIDE_LIKE:
        return False, 0, "peptide_like_requires_force"
    if rep.ligand_atoms_in_apo <= 0 and rep.coord_match_atoms <= 0:
        return False, 0, "no_atoms_to_strip"
    if rep.ligand_atoms_in_apo <= 0:
        return False, 0, "coord_only_needs_manual_review"
    return True, rep.ligand_atoms_in_apo, "strip_residue_name"


def discover_targets(tree: Path) -> list[str]:
    if not tree.is_dir():
        return []
    out = []
    for p in sorted(tree.iterdir()):
        if not p.is_dir():
            continue
        if (p / f"{p.name}_apo.pdb").exists() or (p / f"{p.name}.pdb").exists():
            out.append(p.name)
    return out


def write_csv(path: Path, rows: Iterable[TargetReport]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fieldnames = list(TargetReport.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def write_summary_json(path: Path, rows: list[TargetReport], tree: Path, extra: dict | None = None) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(rows)
    n_fail = sum(1 for r in rows if r.status == "fail")
    n_warn = sum(1 for r in rows if r.status == "warn")
    n_pass = sum(1 for r in rows if r.status in ("pass", "ok"))
    n_ident = sum(1 for r in rows if r.identical_sha_to_deposit)
    n_lig = sum(1 for r in rows if r.ligand_atoms_in_apo > 0)
    n_coord = sum(1 for r in rows if r.coord_match_atoms > 0)
    n_fixed = sum(1 for r in rows if r.fix_applied)
    n_planned = sum(1 for r in rows if r.fix_planned)

    fails = [r for r in rows if r.status == "fail"]

    def _worst_key(r: TargetReport) -> tuple:
        return (
            0 if r.status == "fail" else 1 if r.status == "warn" else 2,
            0 if not r.identical_sha_to_deposit else 1,
            -r.ligand_atoms_in_apo,
            -r.coord_match_atoms,
            -r.nonstd_hetatm_atoms,
            r.pdb_id,
        )

    candidates = [
        r
        for r in rows
        if r.status in ("fail", "warn")
        or r.ligand_atoms_in_apo > 0
        or r.coord_match_atoms > 0
        or not r.identical_sha_to_deposit
        or r.nonstd_hetatm_atoms >= 50
    ]
    ranked = sorted(candidates, key=_worst_key)
    worst = [
        {
            "pdb_id": r.pdb_id,
            "ligand_title": r.ligand_title,
            "status": r.status,
            "pass_fail": r.pass_fail,
            "identical_sha_to_deposit": r.identical_sha_to_deposit,
            "ligand_atoms_in_apo": r.ligand_atoms_in_apo,
            "ligand_hetatm_names": r.ligand_hetatm_names,
            "coord_match_atoms": r.coord_match_atoms,
            "nonstd_hetatm_residues": r.nonstd_hetatm_residues,
            "notes": r.notes,
        }
        for r in ranked[:15]
    ]
    non_id = [
        {
            "pdb_id": r.pdb_id,
            "ligand_title": r.ligand_title,
            "status": r.status,
            "identical_sha_to_deposit": r.identical_sha_to_deposit,
            "ligand_atoms_in_apo": r.ligand_atoms_in_apo,
            "coord_match_atoms": r.coord_match_atoms,
            "notes": r.notes,
        }
        for r in rows
        if not r.identical_sha_to_deposit
    ]

    try:
        tree_s = str(tree.relative_to(REPO_ROOT))
    except ValueError:
        tree_s = str(tree)

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tree": tree_s,
        "n_targets": n,
        "n_pass": n_pass,
        "n_warn": n_warn,
        "n_fail": n_fail,
        "n_identical_sha_to_deposit": n_ident,
        "n_ligand_atoms_remaining": n_lig,
        "n_coord_match": n_coord,
        "n_fix_planned": n_planned,
        "n_fix_applied": n_fixed,
        "fail_pdb_ids": [r.pdb_id for r in fails],
        "non_identical_apo_deposit": non_id,
        "worst_offenders": worst,
        "interpretation": (
            "status=fail means cognate ligand residue atoms (or strong coordinate "
            "match) remain in *_apo.pdb. Byte-identity of apo to deposit alone is "
            "NOT a fail when the cognate ligand is already absent from the deposit "
            "PDB (ligand lives in CIF/SDF only)."
        ),
    }
    if extra:
        summary.update(extra)
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
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
        "--summary-json",
        type=Path,
        default=DEFAULT_SUMMARY,
        help=f"Output summary JSON (default: {DEFAULT_SUMMARY.relative_to(REPO_ROOT)})",
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
        help="Å tolerance for SDF/MOL2↔PDB coordinate ligand match (default 0.35)",
    )
    ap.add_argument(
        "--fix-dry-run",
        action="store_true",
        help="Describe strip operations that would run; never write apo files",
    )
    ap.add_argument(
        "--write",
        action="store_true",
        help=(
            "Actually strip cognate ligand residue from apo (writes .bak). "
            f"Capped at {DEFAULT_WRITE_PILOT_LIMIT} targets unless --all-safe."
        ),
    )
    ap.add_argument(
        "--fix",
        action="store_true",
        help="Alias for --write (backward compatible)",
    )
    ap.add_argument(
        "--all-safe",
        action="store_true",
        help="Allow --write on all strip-safe fail targets (not just pilot ≤3)",
    )
    ap.add_argument(
        "--force-peptide",
        action="store_true",
        help="Allow strip when ligand title is a standard amino-acid code",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any target has status=fail (explicit fail gate)",
    )
    ap.add_argument(
        "--allow-fail",
        action="store_true",
        help="Exit 0 even if some targets have status=fail (overrides --strict)",
    )
    ap.add_argument(
        "--summary-only",
        action="store_true",
        help="Print summary counts only (still writes CSV + JSON)",
    )
    args = ap.parse_args(argv)

    do_write = bool(args.write or args.fix)
    if do_write and args.fix_dry_run:
        print("ERROR: use either --write/--fix or --fix-dry-run, not both", file=sys.stderr)
        return 2

    tree = args.tree if args.tree.is_absolute() else REPO_ROOT / args.tree
    report_path = args.report if args.report.is_absolute() else REPO_ROOT / args.report
    summary_path = (
        args.summary_json if args.summary_json.is_absolute() else REPO_ROOT / args.summary_json
    )

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
        planned, n_atoms, reason = plan_strip(rep)
        if rep.ligand_title in PEPTIDE_LIKE and rep.status == "fail" and args.force_peptide:
            planned, n_atoms, reason = True, rep.ligand_atoms_in_apo, "force_peptide"
        if args.fix_dry_run or do_write:
            if planned:
                rep.fix_planned = True
                rep.fix_atoms_planned = n_atoms
                rep.notes = (rep.notes + ";" if rep.notes else "") + f"strip_plan:{reason}:{n_atoms}"
            elif rep.status == "fail":
                rep.notes = (rep.notes + ";" if rep.notes else "") + f"strip_skip:{reason}"
        rows.append(rep)

    planned_rows = [r for r in rows if r.fix_planned]
    if do_write:
        write_targets = planned_rows
        if not args.all_safe and len(planned_rows) > DEFAULT_WRITE_PILOT_LIMIT:
            write_targets = planned_rows[:DEFAULT_WRITE_PILOT_LIMIT]
            print(
                f"NOTE: --write pilot cap: applying strip to {len(write_targets)}/"
                f"{len(planned_rows)} planned targets "
                f"(limit={DEFAULT_WRITE_PILOT_LIMIT}; pass --all-safe to strip all)",
                file=sys.stderr,
            )
        write_ids = {r.pdb_id for r in write_targets}
        new_rows: list[TargetReport] = []
        for rep in rows:
            if rep.pdb_id not in write_ids:
                new_rows.append(rep)
                continue
            if rep.ligand_title in PEPTIDE_LIKE and not args.force_peptide:
                rep.notes = (rep.notes + ";" if rep.notes else "") + "fix_skipped_peptide_like"
                new_rows.append(rep)
                continue
            apo = tree / rep.pdb_id / f"{rep.pdb_id}_apo.pdb"
            removed = strip_ligand_residue(apo, rep.ligand_title, backup=True, dry_run=False)
            re_rep = analyze_target(rep.pdb_id, tree, coord_tol=args.coord_tol)
            re_rep.fix_planned = True
            re_rep.fix_atoms_planned = rep.fix_atoms_planned
            re_rep.fix_applied = removed > 0
            re_rep.fix_atoms_removed = removed
            if removed:
                re_rep.notes = (
                    (re_rep.notes + ";" if re_rep.notes else "") + f"stripped_{removed}_atoms"
                )
            new_rows.append(re_rep)
        rows = new_rows

    write_csv(report_path, rows)
    try:
        csv_rel = str(report_path.relative_to(REPO_ROOT))
    except ValueError:
        csv_rel = str(report_path)
    summary = write_summary_json(
        summary_path,
        rows,
        tree,
        extra={
            "fix_dry_run": bool(args.fix_dry_run),
            "write": do_write,
            "all_safe": bool(args.all_safe),
            "strict": bool(args.strict),
            "csv_report": csv_rel,
        },
    )

    n = summary["n_targets"]
    n_fail = summary["n_fail"]
    n_warn = summary["n_warn"]
    n_pass = summary["n_pass"]
    n_ident = summary["n_identical_sha_to_deposit"]
    n_lig = summary["n_ligand_atoms_remaining"]
    n_fixed = summary["n_fix_applied"]
    n_planned = summary["n_fix_planned"]

    print(f"Astex apo strip validation — {n} targets")
    print(f"  tree:     {tree}")
    print(f"  report:   {report_path}")
    print(f"  summary:  {summary_path}")
    print(f"  status:   pass={n_pass}  warn={n_warn}  fail={n_fail}")
    print(f"  identical apo==deposit (SHA): {n_ident}/{n}")
    print(f"  ligand residue atoms remaining: {n_lig}/{n}")
    if args.fix_dry_run:
        print(f"  fix dry-run planned: {n_planned}")
        for r in rows:
            if r.fix_planned:
                print(
                    f"    would strip {r.pdb_id}: residue={r.ligand_title!r} "
                    f"atoms≈{r.fix_atoms_planned} apo={r.apo_path}"
                )
            elif r.status == "fail":
                print(f"    skip {r.pdb_id}: {r.notes}")
    if do_write:
        print(f"  fix applied: {n_fixed}")

    if not args.summary_only:
        shown = 0
        for r in rows:
            if r.status in ("fail", "warn") and shown < 20:
                print(
                    f"  [{r.status}/{r.pass_fail}] {r.pdb_id} lig={r.ligand_title!r} "
                    f"atoms={r.ligand_atoms_in_apo} het={r.ligand_hetatm_names!r} "
                    f"coord={r.coord_match_atoms} ident={r.identical_sha_to_deposit} "
                    f"notes={r.notes}"
                )
                shown += 1
        if n_fail + n_warn > shown:
            print(f"  ... ({n_fail + n_warn - shown} more non-pass in CSV)")
        non_id = [r for r in rows if not r.identical_sha_to_deposit]
        if non_id:
            print(f"  non-identical apo vs deposit ({len(non_id)}):")
            for r in non_id[:10]:
                print(
                    f"    {r.pdb_id} lig={r.ligand_title!r} status={r.status} "
                    f"atoms={r.ligand_atoms_in_apo} notes={r.notes}"
                )

    # Exit policy: fail on residual ligand unless --allow-fail.
    # --strict is the explicit mission gate (same outcome as default fail).
    if n_fail and not args.allow_fail:
        return 1
    if args.strict and n_fail and not args.allow_fail:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
