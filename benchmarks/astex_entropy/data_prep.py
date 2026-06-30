"""data_prep.py - Prepare inputs for Vina, rDock, Boltz-2 from your Astex yamls.

Minimal, focused on native / non_native modes.
Outputs ready-to-use files under out_dir.
"""

from pathlib import Path
from typing import List, Dict
import yaml
import subprocess
import shutil

from .config import (
    ASTEX_DIVERSE_YAML, ASTEX_NONNATIVE_YAML,
    ASTEX_DIVERSE_DIR, ASTEX_NONNATIVE_DIR,
    DEFAULT_DATA_DIR,
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None

try:
    from openbabel import pybel
except ImportError:
    pybel = None


def load_targets(yaml_path: Path, mode: str = "native") -> List[Dict]:
    with open(yaml_path) as f:
        data = yaml.safe_load(f)
    targets = data.get("targets", [])
    # For non_native the yaml is lighter; real pairs from astex_nonnative/pairs or your loader.
    # Here we take the listed targets; pairing logic can be extended.
    return [{"id": t, "mode": mode} for t in targets]


def prepare_for_vina(target: Dict, out_dir: Path) -> Path:
    """Create pdbqt + simple config. Use existing structures."""
    tid = target["id"].upper()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Naive: use existing holo pdb as receptor, extract ligand or use sdf
    rec_pdb = ASTEX_DIVERSE_DIR / "structures" / f"{tid}.pdb"
    lig_sdf = ASTEX_DIVERSE_DIR / "data" / tid.lower() / f"{tid.lower()}.sdf"  # approximate

    if not rec_pdb.exists():
        # fallback for non_native or other layout
        rec_pdb = ASTEX_DIVERSE_DIR / "astex_diverse" / tid / f"{tid}_rec.pdb"  # example

    rec_pdbqt = out_dir / f"{tid}_rec.pdbqt"
    lig_pdbqt = out_dir / f"{tid}_lig.pdbqt"

    # Use obabel (openbabel) for conversion - minimal dep
    if pybel:
        # receptor
        mol = next(pybel.readfile("pdb", str(rec_pdb)))
        mol.write("pdbqt", str(rec_pdbqt), overwrite=True)
        # ligand - assume sdf or extract
        if lig_sdf.exists():
            mol = next(pybel.readfile("sdf", str(lig_sdf)))
            mol.write("pdbqt", str(lig_pdbqt), overwrite=True)
    else:
        # fallback rdkit + manual (very basic)
        if Chem and lig_sdf.exists():
            mol = Chem.MolFromMolFile(str(lig_sdf))
            if mol:
                # write simple pdbqt stub (real prep usually needs AutoDockTools)
                with open(lig_pdbqt, "w") as f:
                    f.write("@<TRIPOS>MOLECULE\n")  # stub, user should improve with mgltools if needed

    # Box from ligand or default ~20A
    cfg = out_dir / f"{tid}_vina.txt"
    with open(cfg, "w") as f:
        f.write(f"""center_x = 0
center_y = 0
center_z = 0
size_x = 22.5
size_y = 22.5
size_z = 22.5
num_modes = 10
""")
    return cfg


def prepare_for_rdock(target: Dict, out_dir: Path) -> Path:
    """Very basic .prm and cavity setup using existing data. Real .as grid needs rbcavity."""
    tid = target["id"].upper()
    out_dir.mkdir(parents=True, exist_ok=True)

    prm = out_dir / f"{tid}.prm"
    with open(prm, "w") as f:
        f.write(f"""RBT_PARAMETER_FILE_V1.00
TITLE {tid}
RECEPTOR_FILE {ASTEX_DIVERSE_DIR}/structures/{tid}.pdb
LIGAND_FILE {ASTEX_DIVERSE_DIR}/data/{tid.lower()}/{tid.lower()}.sdf
        """)
    return prm


def prepare_for_boltz(target: Dict, out_dir: Path) -> Path:
    """Boltz-2 input (yaml). Assumes boltz can take protein + ligand."""
    tid = target["id"].upper()
    out_dir.mkdir(parents=True, exist_ok=True)
    inp = out_dir / f"{tid}_boltz.yaml"
    # Minimal example - user will adjust sequences/ligands
    with open(inp, "w") as f:
        f.write(f"""protein: {tid}
ligand: {tid.lower()}
# Use your actual fasta/smiles or let boltz fetch
""")
    return inp


def data_prep(out_dir: Path = DEFAULT_DATA_DIR, modes: List[str] = None):
    """Main entry for data_prep command."""
    if modes is None:
        modes = ["native", "non_native"]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for mode in modes:
        yaml_path = ASTEX_DIVERSE_YAML if mode == "native" else ASTEX_NONNATIVE_YAML
        targets = load_targets(yaml_path, mode)

        mode_dir = out_dir / mode
        for t in targets:
            tdir = mode_dir / t["id"]
            prepare_for_vina(t, tdir / "vina")
            prepare_for_rdock(t, tdir / "rdock")
            prepare_for_boltz(t, tdir / "boltz")

    print(f"Data prep done under {out_dir}")
    print("Review the generated inputs and adjust boxes/ligands for your tools.")


if __name__ == "__main__":
    data_prep()
