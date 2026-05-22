"""truncate_chain.py — Write a nascent-chain PDB at a given length.

MVP: extended geometry (Cα atoms along +X axis at 3.8 Å spacing). Post-MVP can
plug in ESMFold or AlphaFold-derived φ/ψ angles per residue.

Used by ``scripts/run_dual_assembly_cotranslational.sh`` when the runner needs a
ligand PDB at chain length L_k for the cotranslational docking pipeline.

Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Three-letter codes for the canonical 20 amino acids (used in the PDB output).
_AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    "X": "UNK",
}

_CA_CA_DIST_A = 3.8  # canonical Cα–Cα separation in an extended chain


def write_extended_ca_chain(sequence: str, length: int, output_path: Path) -> Path:
    """Write a 1-atom-per-residue PDB at ``output_path`` for the first ``length``
    residues of ``sequence``. The Cα atoms lie along +X at 3.8 Å spacing.

    Returns ``output_path`` for convenience.
    """
    seq = sequence[:length]
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as fh:
        fh.write(f"REMARK 100  Synthetic nascent chain at L_k = {length}\n")
        fh.write(f"REMARK 100  Generator: flexaidds.truncate_chain (extended geometry)\n")
        for i, code in enumerate(seq, start=1):
            res3 = _AA3.get(code.upper(), "UNK")
            x = _CA_CA_DIST_A * (i - 1)
            fh.write(
                "ATOM  {atom_idx:>5d}  CA  {res3:<3s} A{res_idx:>4d}    "
                "{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           C\n".format(
                    atom_idx=i, res3=res3, res_idx=i, x=x, y=0.0, z=0.0,
                )
            )
        fh.write("TER\nEND\n")
    return output_path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sequence", required=True,
                   help="1-letter amino-acid sequence, OR @path/to.fasta")
    p.add_argument("--length", type=int, required=True,
                   help="Chain length L_k to emit")
    p.add_argument("--output", required=True, type=Path,
                   help="Destination PDB path")
    args = p.parse_args(argv)

    seq = args.sequence
    if seq.startswith("@"):
        seq_path = Path(seq[1:])
        seq = "".join(
            line.strip() for line in seq_path.read_text().splitlines()
            if not line.startswith(">")
        ).replace(" ", "")

    if not seq:
        print("ERROR: empty sequence", file=sys.stderr)
        return 2
    if args.length <= 0 or args.length > len(seq):
        print(f"ERROR: length must be in 1..{len(seq)}", file=sys.stderr)
        return 2

    write_extended_ca_chain(seq, args.length, args.output)
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
