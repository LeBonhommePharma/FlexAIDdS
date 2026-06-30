"""runners.py - Thin wrappers to run Vina, rDock, Boltz-2.

Keep simple. Subprocess calls. User ensures binaries + inputs ready.
No over-abstraction.
"""

from pathlib import Path
import subprocess
import shutil
from typing import List, Optional

from .config import VINA_BIN, RDOCK_RBCAVITY, RDOCK_RBDOCK, BOLTZ_BIN


def run_vina(rec_pdbqt: Path, lig_pdbqt: Path, cfg: Path, out_sdf: Path, n_poses: int = 10) -> Path:
    """Run vina. Expects pdbqt inputs and config with center/size."""
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        VINA_BIN,
        "--receptor", str(rec_pdbqt),
        "--ligand", str(lig_pdbqt),
        "--config", str(cfg),
        "--num_modes", str(n_poses),
        "--out", str(out_sdf.with_suffix(".pdbqt")),
    ]
    subprocess.check_call(cmd)
    # Convert pdbqt -> sdf with obabel if available
    try:
        subprocess.check_call(["obabel", "-ipdbqt", str(out_sdf.with_suffix(".pdbqt")), "-osdf", "-O", str(out_sdf)])
    except Exception:
        shutil.copy(out_sdf.with_suffix(".pdbqt"), out_sdf)  # fallback
    return out_sdf


def run_rdock(prm: Path, lig: Path, out_dir: Path) -> Path:
    """Run rbcavity then rbdock. Very thin."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cavity = out_dir / "cavity.prm"
    subprocess.check_call([RDOCK_RBCAVITY, "-r", str(prm), "-was"])
    # rbdock
    out = out_dir / "rdock_out.sdf"
    cmd = [RDOCK_RBDOCK, "-i", str(lig), "-o", str(out.with_suffix("")), "-r", str(prm), "-p", "dock.prm", "-n", "10"]
    subprocess.check_call(cmd)
    return out


def run_boltz(input_yaml: Path, out_dir: Path) -> Path:
    """Run Boltz-2. Assumes it produces a pose file."""
    out_dir.mkdir(parents=True, exist_ok=True)
    # Example invocation (adjust to your boltz version)
    cmd = [BOLTZ_BIN, "predict", str(input_yaml), "--out_dir", str(out_dir)]
    subprocess.check_call(cmd)
    # Return first pose sdf if present
    for f in out_dir.rglob("*.sdf"):
        return f
    return out_dir / "predicted.sdf"
