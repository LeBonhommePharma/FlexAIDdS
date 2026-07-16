#!/usr/bin/env python3
"""Ligand integrity checks for FlexAID prep and pose emission (pure Python).

Validates that LIG_ref / INI / pose PDBs have consistent heavy-atom topology:

  * heavy-atom count match (LIG_ref vs pose/INI)
  * optional max pairwise CONECT bond distance (default 3.0 Å)

Used by:
  * scripts/validate_ligand_integrity.py  (standalone CLI / canary)
  * scripts/generate_flexaid_inp.py       (prep gate after ProcessLigand)

Exit codes (CLI wrappers):
  0  OK
  2  usage / missing inputs
 10  missing LIG_ref
 11  missing pose/INI
 12  heavy-atom count mismatch
 13  CONECT bond length exceeds max
 14  parse failure (no heavy atoms)
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

# Sentinel elements / names treated as hydrogen-like (not heavy).
_H_LIKE = frozenset({"H", "D", "T", "DU", ""})

# FlexAID often writes serials ≥ 90000 for ligand atoms.
_LIGAND_SERIAL_MIN = 90000

DEFAULT_MAX_BOND_A = 3.0

# Exit codes shared with CLI
EXIT_OK = 0
EXIT_USAGE = 2
EXIT_MISSING_REF = 10
EXIT_MISSING_POSE = 11
EXIT_COUNT_MISMATCH = 12
EXIT_BOND_TOO_LONG = 13
EXIT_PARSE_FAIL = 14


@dataclass
class AtomRec:
    serial: int
    name: str
    resname: str
    x: float
    y: float
    z: float
    element: str
    record: str  # ATOM / HETATM


@dataclass
class BondBreach:
    a: int
    b: int
    distance: float


@dataclass
class IntegrityResult:
    ok: bool
    exit_code: int
    ref_path: str = ""
    pose_path: str = ""
    ref_heavy: int = 0
    pose_heavy: int = 0
    max_bond_a: Optional[float] = None
    max_bond_limit: float = DEFAULT_MAX_BOND_A
    bond_breaches: List[BondBreach] = field(default_factory=list)
    missing_serials_in_pose: List[int] = field(default_factory=list)
    extra_serials_in_pose: List[int] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        d = asdict(self)
        d["bond_breaches"] = [
            {"a": b.a, "b": b.b, "distance": round(b.distance, 4)} for b in self.bond_breaches
        ]
        return d


def _element_from_line(line: str) -> str:
    """PDB element (cols 77-78) or first alpha of atom name."""
    if len(line) >= 78:
        elem = line[76:78].strip().upper()
        if elem:
            return elem
    name = line[12:16] if len(line) >= 16 else ""
    for c in name:
        if c.isalpha():
            return c.upper()
    return ""


# Real elements that start with H (must not be classified as hydrogen).
_H_START_ELEMENTS = frozenset({"HE", "HF", "HG", "HO", "HS"})


def is_heavy_element(elem: str) -> bool:
    """True for non-hydrogen chemical elements (C, N, O, F, CL, BR, …)."""
    e = (elem or "").strip().upper()
    if not e or e in _H_LIKE:
        return False
    if e in _H_START_ELEMENTS:
        return True
    # Bare H / isotopes / numbered hydrogens
    if e == "H" or (e.startswith("H") and e[1:].isdigit()):
        return False
    # Single-letter non-H
    if len(e) == 1:
        return e != "H"
    # Two-letter element codes that are not H-named hydrogens (HA/HB from atom names)
    if len(e) == 2 and e[0] == "H" and e[1].isalpha() and e not in _H_START_ELEMENTS:
        return False
    return True


def parse_pdb_atoms(text: str) -> List[AtomRec]:
    atoms: List[AtomRec] = []
    for line in text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        try:
            serial = int(line[6:11])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        name = line[12:16].strip() if len(line) >= 16 else ""
        resname = line[17:20].strip() if len(line) >= 20 else ""
        elem = _element_from_line(line)
        atoms.append(
            AtomRec(
                serial=serial,
                name=name,
                resname=resname,
                x=x,
                y=y,
                z=z,
                element=elem,
                record=line[:6].strip(),
            )
        )
    return atoms


def parse_conect(text: str) -> List[Tuple[int, int]]:
    """Return undirected bond pairs from CONECT records (deduped a < b)."""
    seen: Set[Tuple[int, int]] = set()
    pairs: List[Tuple[int, int]] = []
    for line in text.splitlines():
        if not line.startswith("CONECT"):
            continue
        # Fixed-width 5-char serial fields after "CONECT"
        nums: List[int] = []
        body = line[6:]
        # Prefer fixed-width parse; fall back to whitespace split
        if body.strip():
            chunks = [body[i : i + 5] for i in range(0, len(body), 5)]
            for ch in chunks:
                s = ch.strip()
                if not s:
                    continue
                try:
                    nums.append(int(s))
                except ValueError:
                    # Whitespace-delimited fallback for nonstandard writers
                    nums = []
                    for tok in body.split():
                        try:
                            nums.append(int(tok))
                        except ValueError:
                            pass
                    break
        if len(nums) < 2:
            continue
        src = nums[0]
        for dst in nums[1:]:
            if src == dst:
                continue
            a, b = (src, dst) if src < dst else (dst, src)
            if (a, b) not in seen:
                seen.add((a, b))
                pairs.append((a, b))
    return pairs


def heavy_atoms(
    atoms: Sequence[AtomRec],
    *,
    ligand_only: bool = False,
    resname: Optional[str] = None,
    serial_min: Optional[int] = None,
) -> List[AtomRec]:
    out: List[AtomRec] = []
    for a in atoms:
        if not is_heavy_element(a.element):
            continue
        if ligand_only:
            high = serial_min is not None and a.serial >= serial_min
            res_ok = a.resname.upper() in {
                (resname or "").upper(),
                "LIG",
                "UNL",
                "GEO",
            }
            if not high and not res_ok:
                continue
        out.append(a)
    return out


def count_heavy_atoms_pdb(
    path: Path,
    *,
    ligand_only: bool = False,
    prefer_high_serial: bool = True,
) -> Tuple[int, List[AtomRec]]:
    """Count heavy atoms in a PDB file.

    For pose/INI files that contain receptor + ligand, set ligand_only=True
    (selects serial >= 90000 and/or resname LIG/GEO/UNL).
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    atoms = parse_pdb_atoms(text)
    if ligand_only:
        high = [a for a in atoms if a.serial >= _LIGAND_SERIAL_MIN and is_heavy_element(a.element)]
        if high and prefer_high_serial:
            return len(high), high
        lig_res = [
            a
            for a in atoms
            if a.resname.upper() in {"LIG", "UNL", "GEO"} and is_heavy_element(a.element)
        ]
        if lig_res:
            return len(lig_res), lig_res
        # Fallback: all HETATM heavy
        het = [
            a
            for a in atoms
            if a.record == "HETATM" and is_heavy_element(a.element)
        ]
        return len(het), het
    heavy = [a for a in atoms if is_heavy_element(a.element)]
    return len(heavy), heavy


def max_conect_bond_distance(
    path: Path,
    *,
    ligand_only: bool = True,
) -> Tuple[Optional[float], List[BondBreach], Dict[int, AtomRec]]:
    """Max distance among CONECT neighbors; list breaches is empty here.

    Returns (max_distance_or_None, [], atom_map).
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    atoms = parse_pdb_atoms(text)
    atom_map: Dict[int, AtomRec] = {a.serial: a for a in atoms if is_heavy_element(a.element)}
    if ligand_only:
        atom_map = {
            s: a
            for s, a in atom_map.items()
            if s >= _LIGAND_SERIAL_MIN or a.resname.upper() in {"LIG", "UNL", "GEO"}
        }
    pairs = parse_conect(text)
    max_d: Optional[float] = None
    for a, b in pairs:
        if a not in atom_map or b not in atom_map:
            continue
        aa, bb = atom_map[a], atom_map[b]
        d = math.sqrt(
            (aa.x - bb.x) ** 2 + (aa.y - bb.y) ** 2 + (aa.z - bb.z) ** 2
        )
        if max_d is None or d > max_d:
            max_d = d
    return max_d, [], atom_map


def find_bond_breaches(
    path: Path,
    max_bond: float = DEFAULT_MAX_BOND_A,
    *,
    ligand_only: bool = True,
) -> Tuple[Optional[float], List[BondBreach]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    atoms = parse_pdb_atoms(text)
    atom_map: Dict[int, AtomRec] = {a.serial: a for a in atoms}
    if ligand_only:
        atom_map = {
            s: a
            for s, a in atom_map.items()
            if is_heavy_element(a.element)
            and (
                s >= _LIGAND_SERIAL_MIN
                or a.resname.upper() in {"LIG", "UNL", "GEO"}
            )
        }
    else:
        atom_map = {s: a for s, a in atom_map.items() if is_heavy_element(a.element)}

    pairs = parse_conect(text)
    breaches: List[BondBreach] = []
    max_d: Optional[float] = None
    for a, b in pairs:
        if a not in atom_map or b not in atom_map:
            continue
        aa, bb = atom_map[a], atom_map[b]
        d = math.sqrt(
            (aa.x - bb.x) ** 2 + (aa.y - bb.y) ** 2 + (aa.z - bb.z) ** 2
        )
        if max_d is None or d > max_d:
            max_d = d
        if d > max_bond:
            breaches.append(BondBreach(a=a, b=b, distance=d))
    return max_d, breaches


def validate_ligand_integrity(
    ref_path: Path,
    pose_path: Optional[Path] = None,
    *,
    max_bond: float = DEFAULT_MAX_BOND_A,
    check_bonds: bool = True,
    pose_ligand_only: bool = True,
) -> IntegrityResult:
    """Compare LIG_ref vs pose/INI heavy-atom integrity.

    If pose_path is None, only validates that ref parses and (optionally)
    its own CONECT bond lengths are sane.
    """
    ref_path = Path(ref_path)
    messages: List[str] = []

    if not ref_path.is_file():
        return IntegrityResult(
            ok=False,
            exit_code=EXIT_MISSING_REF,
            ref_path=str(ref_path),
            messages=[f"missing LIG_ref: {ref_path}"],
        )

    try:
        ref_n, ref_atoms = count_heavy_atoms_pdb(ref_path, ligand_only=False)
    except OSError as e:
        return IntegrityResult(
            ok=False,
            exit_code=EXIT_PARSE_FAIL,
            ref_path=str(ref_path),
            messages=[f"cannot read LIG_ref: {e}"],
        )

    if ref_n <= 0:
        return IntegrityResult(
            ok=False,
            exit_code=EXIT_PARSE_FAIL,
            ref_path=str(ref_path),
            messages=["LIG_ref has zero heavy atoms"],
        )

    # Self-check bonds on ref if no pose
    check_path = pose_path if pose_path is not None else ref_path
    pose_n = ref_n
    pose_atoms: List[AtomRec] = ref_atoms
    missing: List[int] = []
    extra: List[int] = []

    if pose_path is not None:
        pose_path = Path(pose_path)
        if not pose_path.is_file():
            return IntegrityResult(
                ok=False,
                exit_code=EXIT_MISSING_POSE,
                ref_path=str(ref_path),
                pose_path=str(pose_path),
                ref_heavy=ref_n,
                messages=[f"missing pose/INI: {pose_path}"],
            )
        try:
            pose_n, pose_atoms = count_heavy_atoms_pdb(
                pose_path, ligand_only=pose_ligand_only
            )
        except OSError as e:
            return IntegrityResult(
                ok=False,
                exit_code=EXIT_PARSE_FAIL,
                ref_path=str(ref_path),
                pose_path=str(pose_path),
                ref_heavy=ref_n,
                messages=[f"cannot read pose: {e}"],
            )

        ref_serials = {a.serial for a in ref_atoms}
        pose_serials = {a.serial for a in pose_atoms}
        # Only report serial diffs when both use high-serial ligand numbering
        if ref_serials and min(ref_serials) >= _LIGAND_SERIAL_MIN:
            missing = sorted(ref_serials - pose_serials)
            extra = sorted(pose_serials - ref_serials)

        if pose_n != ref_n:
            msg = (
                f"heavy-atom count mismatch: LIG_ref={ref_n} pose={pose_n}"
            )
            if missing:
                msg += f" missing_serials={missing[:12]}"
            if extra:
                msg += f" extra_serials={extra[:12]}"
            messages.append(msg)
            return IntegrityResult(
                ok=False,
                exit_code=EXIT_COUNT_MISMATCH,
                ref_path=str(ref_path),
                pose_path=str(pose_path),
                ref_heavy=ref_n,
                pose_heavy=pose_n,
                missing_serials_in_pose=missing,
                extra_serials_in_pose=extra,
                messages=messages,
            )

    max_d: Optional[float] = None
    breaches: List[BondBreach] = []
    if check_bonds and check_path is not None and Path(check_path).is_file():
        max_d, breaches = find_bond_breaches(
            Path(check_path),
            max_bond=max_bond,
            # Pose/INI include receptor; LIG_ref is ligand-only — both OK with True.
            ligand_only=True,
        )
        if breaches:
            messages.append(
                f"CONECT bond(s) exceed {max_bond:.1f} Å "
                f"(max={max_d:.3f} if known): "
                + "; ".join(
                    f"{b.a}-{b.b}={b.distance:.3f}" for b in breaches[:8]
                )
            )
            return IntegrityResult(
                ok=False,
                exit_code=EXIT_BOND_TOO_LONG,
                ref_path=str(ref_path),
                pose_path=str(pose_path) if pose_path else "",
                ref_heavy=ref_n,
                pose_heavy=pose_n,
                max_bond_a=max_d,
                max_bond_limit=max_bond,
                bond_breaches=breaches,
                missing_serials_in_pose=missing,
                extra_serials_in_pose=extra,
                messages=messages,
            )

    messages.append(
        f"OK heavy_atoms={ref_n}"
        + (f" pose={pose_n}" if pose_path else "")
        + (f" max_conect_bond={max_d:.3f}Å" if max_d is not None else " (no CONECT bonds)")
    )
    return IntegrityResult(
        ok=True,
        exit_code=EXIT_OK,
        ref_path=str(ref_path),
        pose_path=str(pose_path) if pose_path else "",
        ref_heavy=ref_n,
        pose_heavy=pose_n,
        max_bond_a=max_d,
        max_bond_limit=max_bond,
        bond_breaches=[],
        missing_serials_in_pose=missing,
        extra_serials_in_pose=extra,
        messages=messages,
    )


def validate_work_dir(
    work: Path,
    *,
    max_bond: float = DEFAULT_MAX_BOND_A,
    check_bonds: bool = True,
    require_ini: bool = False,
) -> IntegrityResult:
    """Validate prep work tree: LIG_ref required; INI if present or required.

    Preflight contract:
      * After ProcessLigand: LIG_ref.pdb self-check (count + CONECT bonds).
      * After FlexAID start: * _INI.pdb ligand heavy count must match LIG_ref.
    """
    work = Path(work)
    ref = work / "LIG_ref.pdb"
    if not ref.is_file():
        # Also accept nested restart dirs? Prefer work-root LIG_ref.
        return IntegrityResult(
            ok=False,
            exit_code=EXIT_MISSING_REF,
            ref_path=str(ref),
            messages=[f"missing LIG_ref under {work}"],
        )

    # Self-check ref first
    base = validate_ligand_integrity(
        ref, None, max_bond=max_bond, check_bonds=check_bonds
    )
    if not base.ok:
        return base

    # Locate INI if any
    ini_candidates = [
        work / "INI.pdb",
        *sorted(work.glob("*_INI.pdb")),
        *sorted(work.glob("*_ini.pdb")),
    ]
    # Search one level of restart_* / r* for post-run canaries
    for sub in sorted(work.iterdir()) if work.is_dir() else []:
        if not sub.is_dir():
            continue
        if sub.name.startswith("restart_") or re.fullmatch(r"r\d+", sub.name):
            ini_candidates.extend(sorted(sub.glob("*_INI.pdb")))
            ini_candidates.append(sub / "INI.pdb")

    ini: Optional[Path] = next((p for p in ini_candidates if p.is_file()), None)
    if ini is None:
        if require_ini:
            return IntegrityResult(
                ok=False,
                exit_code=EXIT_MISSING_POSE,
                ref_path=str(ref),
                ref_heavy=base.ref_heavy,
                messages=[
                    "INI not found (post-FlexAID preflight). "
                    "Run after FlexAID emits *_INI.pdb, or pass --pose explicitly."
                ],
            )
        messages = list(base.messages)
        messages.append(
            "NOTE: no INI yet — prep gate passed on LIG_ref only. "
            "Re-run post-INI: python3 scripts/validate_ligand_integrity.py "
            f"--work {work} --require-ini"
        )
        base.messages = messages
        return base

    return validate_ligand_integrity(
        ref,
        ini,
        max_bond=max_bond,
        check_bonds=check_bonds,
        pose_ligand_only=True,
    )


def format_result(res: IntegrityResult) -> str:
    status = "PASS" if res.ok else "FAIL"
    lines = [f"{status} exit={res.exit_code}"]
    for m in res.messages:
        lines.append(f"  {m}")
    if res.ref_path:
        lines.append(f"  ref={res.ref_path} heavy={res.ref_heavy}")
    if res.pose_path:
        lines.append(f"  pose={res.pose_path} heavy={res.pose_heavy}")
    if res.max_bond_a is not None:
        lines.append(
            f"  max_conect_bond={res.max_bond_a:.3f} Å "
            f"(limit={res.max_bond_limit:.1f})"
        )
    return "\n".join(lines)
