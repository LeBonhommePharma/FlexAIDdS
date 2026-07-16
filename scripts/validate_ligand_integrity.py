#!/usr/bin/env python3
"""Standalone ligand integrity gate (prep + post-INI / pose canary).

Checks LIG_ref heavy-atom count and optionally compares a pose/INI PDB.
Fails closed with stable exit codes (see scripts/ligand_integrity.py).

Usage:
  # Prep gate (after ProcessLigand; LIG_ref self-check)
  python3 scripts/validate_ligand_integrity.py --work $WORK/B0/1P62

  # Post-INI (after FlexAID emits *_INI.pdb)
  python3 scripts/validate_ligand_integrity.py --work $WORK/B0/1P62 --require-ini

  # Explicit paths
  python3 scripts/validate_ligand_integrity.py \\
      --ref LIG_ref.pdb --pose 1P62_0.pdb --max-bond 3.0

  # Pose-only count vs ref (skip bond check)
  python3 scripts/validate_ligand_integrity.py --ref LIG_ref.pdb --pose pose.pdb --no-bonds

Exit codes:
  0  OK
  2  usage
 10  missing LIG_ref
 11  missing pose/INI
 12  heavy-atom count mismatch
 13  CONECT bond > max
 14  parse failure

Preflight note:
  INI is written by FlexAID at run start, not by ProcessLigand. Prep can only
  self-check LIG_ref. Re-run with --require-ini after the first FlexAID launch
  (or pass --pose on any emitted pose PDB).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow `python3 scripts/validate_ligand_integrity.py` without package install.
_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from ligand_integrity import (  # noqa: E402
    DEFAULT_MAX_BOND_A,
    EXIT_USAGE,
    format_result,
    validate_ligand_integrity,
    validate_work_dir,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--work", type=Path, help="Prep work directory containing LIG_ref.pdb")
    ap.add_argument("--ref", type=Path, help="LIG_ref.pdb path")
    ap.add_argument("--pose", type=Path, help="Pose or INI PDB to compare")
    ap.add_argument(
        "--require-ini",
        action="store_true",
        help="With --work: fail if no *_INI.pdb is found",
    )
    ap.add_argument(
        "--max-bond",
        type=float,
        default=DEFAULT_MAX_BOND_A,
        help=f"Max allowed CONECT neighbor distance in Å (default {DEFAULT_MAX_BOND_A})",
    )
    ap.add_argument(
        "--no-bonds",
        action="store_true",
        help="Skip CONECT bond-distance check",
    )
    ap.add_argument(
        "--json",
        dest="json_out",
        type=Path,
        default=None,
        help="Optional JSON report path",
    )
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    check_bonds = not args.no_bonds

    if args.work:
        res = validate_work_dir(
            args.work.expanduser().resolve(),
            max_bond=args.max_bond,
            check_bonds=check_bonds,
            require_ini=args.require_ini,
        )
    elif args.ref:
        res = validate_ligand_integrity(
            args.ref.expanduser().resolve(),
            args.pose.expanduser().resolve() if args.pose else None,
            max_bond=args.max_bond,
            check_bonds=check_bonds,
            pose_ligand_only=True,
        )
    else:
        print("ERROR: provide --work or --ref", file=sys.stderr)
        return EXIT_USAGE

    if not args.quiet:
        print(format_result(res))

    if args.json_out:
        out = args.json_out.expanduser()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(res.as_dict(), indent=2) + "\n", encoding="utf-8")
        if not args.quiet:
            print(f"JSON written to {out}")

    return int(res.exit_code)


if __name__ == "__main__":
    sys.exit(main())
