from __future__ import annotations

import shutil
from pathlib import Path

from rdkit import Chem
from rdkit import RDLogger

from .io_utils import run_command

RDLogger.DisableLog("rdApp.error")
RDLogger.DisableLog("rdApp.warning")


def _normalise_counts_line(line: str) -> str:
    parts = line.split()
    if len(parts) < 3 or not parts[-1].upper().startswith("V2000"):
        return line
    try:
        atoms = int(parts[0])
        bonds = int(parts[1])
    except ValueError:
        return line
    return f"{atoms:>3}{bonds:>3}  0  0  0  0            999 V2000"


def normalise_v2000_counts_text(text: str) -> str:
    blocks: list[str] = []
    for block in text.split("$$$$"):
        lines = block.splitlines()
        if len(lines) >= 4:
            lines[3] = _normalise_counts_line(lines[3])
        normalised = "\n".join(lines).rstrip()
        if normalised:
            blocks.append(normalised)
    if not blocks:
        return text
    return "\n$$$$\n".join(blocks) + "\n$$$$\n"


def normalise_v2000_counts_file(source: str | Path, destination: str | Path | None = None) -> Path:
    source_path = Path(source)
    destination_path = Path(destination) if destination is not None else source_path
    destination_path.write_text(normalise_v2000_counts_text(source_path.read_text(errors="ignore")))
    return destination_path


def _supplier_first(path: str | Path, *, sanitize: bool) -> Chem.Mol | None:
    try:
        supplier = Chem.SDMolSupplier(str(path), sanitize=sanitize, removeHs=False, strictParsing=False)
        return supplier[0] if supplier and len(supplier) else None
    except Exception:
        return None


def _normalised_block_first(path: str | Path, *, sanitize: bool) -> Chem.Mol | None:
    text = normalise_v2000_counts_text(Path(path).read_text(errors="ignore"))
    first_block = text.split("$$$$", 1)[0].rstrip()
    if not first_block:
        return None
    try:
        return Chem.MolFromMolBlock(first_block + "\n", sanitize=sanitize, removeHs=False, strictParsing=False)
    except Exception:
        return None


def read_first_sdf_mol(
    path: str | Path,
    obabel: str | None = None,
    *,
    allow_unsanitized: bool = True,
) -> Chem.Mol:
    mol = _supplier_first(path, sanitize=True)
    if mol is not None:
        return mol

    mol = _normalised_block_first(path, sanitize=True)
    if mol is not None:
        return mol

    if allow_unsanitized:
        mol = _supplier_first(path, sanitize=False)
        if mol is not None:
            return mol
        mol = _normalised_block_first(path, sanitize=False)
        if mol is not None:
            return mol

    executable = obabel or shutil.which("obabel")
    if executable:
        repaired = Path(path).with_suffix(".rdkit.sdf")
        run_command(
            [executable, "-isdf", str(path), "-osdf", "-O", str(repaired)],
            log_path=repaired.with_suffix(".obabel.log"),
            check=False,
        )
        normalise_v2000_counts_file(repaired)
        mol = _supplier_first(repaired, sanitize=True)
        if mol is not None:
            return mol
        if allow_unsanitized:
            mol = _supplier_first(repaired, sanitize=False)
            if mol is not None:
                return mol

    raise RuntimeError(f"Could not read SDF molecule: {path}")


def smiles_from_sdf(path: str | Path, obabel: str | None = None) -> str:
    try:
        mol = read_first_sdf_mol(path, obabel, allow_unsanitized=False)
        return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    except Exception:
        pass

    executable = obabel or shutil.which("obabel")
    if not executable:
        raise RuntimeError(f"Could not derive ligand SMILES without OpenBabel: {path}")

    result = run_command(
        [executable, "-isdf", str(path), "-osmi"],
        log_path=Path(path).with_suffix(".smi.obabel.log"),
        check=False,
    )
    for line in result.stdout.splitlines():
        parts = line.strip().split()
        if parts:
            return parts[0]
    raise RuntimeError(f"Could not derive ligand SMILES: {path}")
