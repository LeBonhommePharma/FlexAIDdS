#!/usr/bin/env python3
"""Clean apo TARGET for classic Astex redock (strip waters / optional metals).

Protocol (default redock):
  * Always strip HOH / WAT / H2O / DOD  (override: FLEXAIDDS_KEEP_HOH=1)
  * Strip common metals (MG, ZN, CA, …) (override: FLEXAIDDS_KEEP_METALS=1)
  * Do NOT strip HEM/cofactors by default (hemoglobin / heme proteins)
  * Clean orphan CONECT records that reference removed serials

EXCHET is often ON so waters may not score — still strip non-protocol HET for
clean apo per Astex redock practice (score pocket without crystal waters/metals
unless the metal is known-required and opt-in).

Usage:
  python3 scripts/clean_target_apo.py input_apo.pdb -o clean_apo.pdb
  python3 scripts/clean_target_apo.py input_apo.pdb --in-place
  python3 scripts/clean_target_apo.py input_apo.pdb -o out.pdb --keep-metals
  python3 scripts/clean_target_apo.py input_apo.pdb -o out.pdb --keep-hoh

Env:
  FLEXAIDDS_KEEP_HOH=1      keep waters
  FLEXAIDDS_KEEP_METALS=1   keep metal ions
  FLEXAIDDS_KEEP_METALS_NEAR_LIGAND=<Å>  keep metals within this distance of
      any ligand heavy atom (default 4.0 when ligand coords provided; 0=off)
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set, Tuple

WATER_RES = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3", "WAT2"})

# Common non-covalent metal ions (residue names and element symbols).
METAL_RES = frozenset(
    {
        "MG",
        "MN",
        "ZN",
        "CA",
        "NA",
        "K",
        "FE",
        "CU",
        "NI",
        "CO",
        "CD",
        "HG",
        "SR",
        "BA",
        "LI",
        "RB",
        "CS",
        "AL",
        "PT",
        "AU",
        "AG",
        "PB",
        "TL",
        "Y",
        "GD",
        "TB",
        "MO",
        "W",
        "V",
        "CR",
    }
)

# Keep these cofactors by default (not stripped).
KEEP_COFACTOR_RES = frozenset(
    {
        "HEM",
        "HEC",
        "HEA",
        "HEB",
        "NAD",
        "NAP",
        "NDP",
        "FAD",
        "FMN",
        "ATP",
        "ADP",
        "AMP",
        "GTP",
        "GDP",
        "SAH",
        "SAM",
    }
)


@dataclass
class CleanReport:
    input_path: str
    output_path: str
    atoms_in: int = 0
    atoms_out: int = 0
    waters_removed: int = 0
    metals_removed: int = 0
    metals_kept_near_ligand: int = 0
    other_removed: int = 0
    conect_removed: int = 0
    conect_rewritten: int = 0
    keep_hoh: bool = False
    keep_metals: bool = False
    metal_near_ligand_a: float = 0.0
    removed_resnames: List[str] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def _env_truthy(name: str) -> bool:
    v = os.environ.get(name, "").strip().lower()
    return v in {"1", "true", "yes", "on"}


def _resname(line: str) -> str:
    if len(line) < 20:
        return ""
    return line[17:20].strip().upper()


def _serial(line: str) -> Optional[int]:
    if len(line) < 11:
        return None
    try:
        return int(line[6:11])
    except ValueError:
        return None


def _element(line: str) -> str:
    if len(line) >= 78:
        e = line[76:78].strip().upper()
        if e:
            return e
    name = line[12:16] if len(line) >= 16 else ""
    for c in name:
        if c.isalpha():
            return c.upper()
    return ""


def _xyz(line: str) -> Optional[Tuple[float, float, float]]:
    if len(line) < 54:
        return None
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        return None


def load_ligand_heavy_xyz(ligand_path: Path) -> List[Tuple[float, float, float]]:
    """Heavy-atom coords from LIG_ref / ligand PDB (or SDF-like PDB columns)."""
    out: List[Tuple[float, float, float]] = []
    if not ligand_path.is_file():
        return out
    for line in ligand_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        elem = _element(line)
        if elem in {"H", "D"}:
            continue
        xyz = _xyz(line)
        if xyz:
            out.append(xyz)
    return out


def _min_dist_to_ligand(
    line: str, lig_xyz: Sequence[Tuple[float, float, float]]
) -> Optional[float]:
    xyz = _xyz(line)
    if xyz is None or not lig_xyz:
        return None
    return min(
        math.sqrt((xyz[0] - lx) ** 2 + (xyz[1] - ly) ** 2 + (xyz[2] - lz) ** 2)
        for lx, ly, lz in lig_xyz
    )


def should_strip_atom(
    line: str,
    *,
    keep_hoh: bool,
    keep_metals: bool,
    extra_strip_res: Optional[Set[str]] = None,
    lig_xyz: Optional[Sequence[Tuple[float, float, float]]] = None,
    metal_near_ligand_a: float = 0.0,
) -> str:
    """Return reason code if atom should be stripped, else empty string."""
    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
        return ""
    res = _resname(line)
    elem = _element(line)
    extra = extra_strip_res or set()

    if res in WATER_RES:
        return "" if keep_hoh else "water"

    # Cofactors (HEM, NAD, …) kept by default for redock unless explicitly extra-stripped.
    if res in KEEP_COFACTOR_RES and res not in extra:
        return ""

    if res in METAL_RES or (elem in METAL_RES and res in METAL_RES):
        if keep_metals:
            return ""
        # Keep catalytic / structural metals that coordinate the crystal ligand.
        if metal_near_ligand_a > 0.0 and lig_xyz:
            d = _min_dist_to_ligand(line, lig_xyz)
            if d is not None and d <= metal_near_ligand_a:
                return ""
        return "metal"

    if res in extra:
        return "extra"

    return ""


def parse_conect_serials(line: str) -> List[int]:
    if not line.startswith("CONECT"):
        return []
    body = line[6:]
    nums: List[int] = []
    chunks = [body[i : i + 5] for i in range(0, len(body), 5)]
    for ch in chunks:
        s = ch.strip()
        if not s:
            continue
        try:
            nums.append(int(s))
        except ValueError:
            nums = []
            for tok in body.split():
                try:
                    nums.append(int(tok))
                except ValueError:
                    pass
            break
    return nums


def rewrite_conect(line: str, kept_serials: Set[int]) -> Optional[str]:
    """Return rewritten CONECT line or None if source atom was removed / empty."""
    nums = parse_conect_serials(line)
    if not nums:
        return None
    src = nums[0]
    if src not in kept_serials:
        return None
    dsts = [n for n in nums[1:] if n in kept_serials]
    if not dsts:
        # Source kept but all neighbors gone — drop orphan CONECT
        return None
    out = f"CONECT{src:5d}" + "".join(f"{d:5d}" for d in dsts)
    return out + "\n" if line.endswith("\n") else out


def resolve_metal_near_ligand_a(explicit: Optional[float] = None) -> float:
    """Å cutoff for keeping metals near ligand; 0 disables. Default 4.0 via env/default."""
    if explicit is not None:
        return float(explicit)
    env = os.environ.get("FLEXAIDDS_KEEP_METALS_NEAR_LIGAND", "").strip()
    if env:
        try:
            return float(env)
        except ValueError:
            return 0.0
    # Default ON at 4.0 Å when ligand coords are supplied by callers.
    return 4.0


def clean_apo_pdb(
    text: str,
    *,
    keep_hoh: bool = False,
    keep_metals: bool = False,
    extra_strip_res: Optional[Iterable[str]] = None,
    lig_xyz: Optional[Sequence[Tuple[float, float, float]]] = None,
    metal_near_ligand_a: float = 0.0,
) -> tuple[str, CleanReport]:
    """Clean PDB text; return (new_text, report)."""
    extra = {r.strip().upper() for r in (extra_strip_res or []) if r.strip()}
    report = CleanReport(
        input_path="",
        output_path="",
        keep_hoh=keep_hoh,
        keep_metals=keep_metals,
        metal_near_ligand_a=float(metal_near_ligand_a or 0.0),
    )
    lines = text.splitlines(keepends=True)
    kept_lines: List[str] = []
    kept_serials: Set[int] = set()
    removed_res: Set[str] = set()
    atom_in = 0
    atom_out = 0

    # First pass: atom/hetatm filter
    pending_conect: List[str] = []
    other_tail: List[str] = []

    for line in lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            atom_in += 1
            reason = should_strip_atom(
                line,
                keep_hoh=keep_hoh,
                keep_metals=keep_metals,
                extra_strip_res=extra,
                lig_xyz=lig_xyz,
                metal_near_ligand_a=metal_near_ligand_a,
            )
            if reason:
                res = _resname(line)
                removed_res.add(res or "?")
                if reason == "water":
                    report.waters_removed += 1
                elif reason == "metal":
                    report.metals_removed += 1
                else:
                    report.other_removed += 1
                continue
            # Count metals kept only because they are ligand-proximal
            res = _resname(line)
            elem = _element(line)
            if (
                not keep_metals
                and metal_near_ligand_a > 0
                and lig_xyz
                and (res in METAL_RES or elem in METAL_RES)
            ):
                d = _min_dist_to_ligand(line, lig_xyz)
                if d is not None and d <= metal_near_ligand_a:
                    report.metals_kept_near_ligand += 1
            ser = _serial(line)
            if ser is not None:
                kept_serials.add(ser)
            kept_lines.append(line)
            atom_out += 1
        elif line.startswith("CONECT"):
            pending_conect.append(line)
        else:
            # Preserve non-atom header/footer; CONECT rewritten after
            if line.startswith("END"):
                other_tail.append(line)
            else:
                kept_lines.append(line)

    # CONECT cleanup
    for cl in pending_conect:
        new_cl = rewrite_conect(cl, kept_serials)
        if new_cl is None:
            report.conect_removed += 1
            continue
        # Detect rewrite vs identity
        if new_cl.rstrip("\n") != cl.rstrip("\n"):
            report.conect_rewritten += 1
        kept_lines.append(new_cl if new_cl.endswith("\n") else new_cl + "\n")

    kept_lines.extend(other_tail)
    # Ensure trailing newline
    out = "".join(kept_lines)
    if out and not out.endswith("\n"):
        out += "\n"

    report.atoms_in = atom_in
    report.atoms_out = atom_out
    report.removed_resnames = sorted(removed_res)
    report.messages.append(
        f"stripped waters={report.waters_removed} metals={report.metals_removed} "
        f"metals_kept_near_lig={report.metals_kept_near_ligand} "
        f"other={report.other_removed}; atoms {atom_in}->{atom_out}; "
        f"conect_removed={report.conect_removed}"
    )
    return out, report


def clean_apo_file(
    src: Path,
    dst: Path,
    *,
    keep_hoh: bool = False,
    keep_metals: bool = False,
    extra_strip_res: Optional[Iterable[str]] = None,
    ligand_ref: Optional[Path] = None,
    metal_near_ligand_a: Optional[float] = None,
) -> CleanReport:
    text = src.read_text(encoding="utf-8", errors="replace")
    lig_xyz: List[Tuple[float, float, float]] = []
    near_a = 0.0
    if ligand_ref is not None and Path(ligand_ref).is_file():
        lig_xyz = load_ligand_heavy_xyz(Path(ligand_ref))
        near_a = resolve_metal_near_ligand_a(metal_near_ligand_a)
    out, report = clean_apo_pdb(
        text,
        keep_hoh=keep_hoh,
        keep_metals=keep_metals,
        extra_strip_res=extra_strip_res,
        lig_xyz=lig_xyz or None,
        metal_near_ligand_a=near_a,
    )
    report.input_path = str(src)
    report.output_path = str(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(out, encoding="utf-8")
    return report


def resolve_keep_flags(
    keep_hoh: Optional[bool] = None,
    keep_metals: Optional[bool] = None,
) -> tuple[bool, bool]:
    """CLI flags override env; env defaults both False (strip)."""
    hoh = keep_hoh if keep_hoh is not None else _env_truthy("FLEXAIDDS_KEEP_HOH")
    metals = (
        keep_metals
        if keep_metals is not None
        else _env_truthy("FLEXAIDDS_KEEP_METALS")
    )
    return hoh, metals


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input", type=Path, help="Input apo PDB")
    ap.add_argument("-o", "--output", type=Path, default=None, help="Output PDB")
    ap.add_argument("--in-place", action="store_true", help="Overwrite input")
    ap.add_argument(
        "--keep-hoh",
        action="store_true",
        help="Keep waters (or set FLEXAIDDS_KEEP_HOH=1)",
    )
    ap.add_argument(
        "--keep-metals",
        action="store_true",
        help="Keep metal ions (or set FLEXAIDDS_KEEP_METALS=1)",
    )
    ap.add_argument(
        "--strip-res",
        default="",
        help="Comma-separated extra residue names to strip",
    )
    ap.add_argument("--dry-run", action="store_true", help="Report only, no write")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    src = args.input.expanduser().resolve()
    if not src.is_file():
        print(f"ERROR: input not found: {src}", file=sys.stderr)
        return 2

    if args.in_place:
        dst = src
    elif args.output:
        dst = args.output.expanduser().resolve()
    elif args.dry_run:
        dst = src  # report-only; no write
    else:
        print("ERROR: provide -o/--output or --in-place (or --dry-run)", file=sys.stderr)
        return 2

    keep_hoh, keep_metals = resolve_keep_flags(
        keep_hoh=True if args.keep_hoh else None,
        keep_metals=True if args.keep_metals else None,
    )
    if args.keep_hoh:
        keep_hoh = True
    if args.keep_metals:
        keep_metals = True

    extra = [r for r in args.strip_res.split(",") if r.strip()]

    text = src.read_text(encoding="utf-8", errors="replace")
    out, report = clean_apo_pdb(
        text, keep_hoh=keep_hoh, keep_metals=keep_metals, extra_strip_res=extra
    )
    report.input_path = str(src)
    report.output_path = str(dst)

    if not args.quiet:
        for m in report.messages:
            print(m)
        if report.removed_resnames:
            print(f"removed_resnames={','.join(report.removed_resnames)}")

    if args.dry_run:
        return 0

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(out, encoding="utf-8")
    if not args.quiet:
        print(f"wrote {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
