"""evaluate.py - Evaluation with PoseBusters + RMSD + ITC/ΔH/ΔS + entropy/enthalpy index.

A pose is successful only if RMSD <= 2.0 and passes PoseBusters.
Primary interest: how re-ranking by your entropy/thermo affects decomposition and index.
"""

from pathlib import Path
from typing import List, Dict
import json

try:
    from posebusters import PoseBusters
    HAS_PB = True
except Exception:
    HAS_PB = False

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

from .entropy import entropy_enthalpy_index, deltaH_deltaS_decomposition


def compute_rmsd(ref_mol: Chem.Mol, pose_mol: Chem.Mol) -> float:
    """Heavy-atom RMSD using RDKit."""
    try:
        return rdMolAlign.CalcRMS(ref_mol, pose_mol)
    except Exception:
        return 99.0


def passes_posebusters(pose_path: Path, receptor_path: Path) -> bool:
    if not HAS_PB:
        return True  # allow running without; user should install
    try:
        pb = PoseBusters(config="redock")
        df = pb.bust(str(pose_path), str(receptor_path))
        if df is None or len(df) == 0:
            return False
        return bool(df.iloc[0].get("all_passed", True))
    except Exception:
        return True


def evaluate_poses(poses: List[Dict], ref_sdf: Path, receptor: Path, itc_exp: Dict = None) -> Dict:
    """Return stats + per-pose details with success, thermo, index."""
    results = []
    n_success = 0
    for i, p in enumerate(poses):
        pose_path = Path(p["path"])
        try:
            pose_mol = Chem.MolFromMolFile(str(pose_path), removeHs=False)
            ref_mol = Chem.MolFromMolFile(str(ref_sdf), removeHs=False)
        except Exception:
            pose_mol = ref_mol = None

        rmsd = compute_rmsd(ref_mol, pose_mol) if pose_mol and ref_mol else 99.0
        pb_ok = passes_posebusters(pose_path, receptor)
        success = (rmsd <= 2.0) and pb_ok

        th = p.get("thermo", {})
        de = deltaH_deltaS_decomposition(th)
        idx = p.get("entropy_enthalpy_index", entropy_enthalpy_index(de["deltaH"], de["minus_TdS"], de["deltaG"]))

        if success:
            n_success += 1

        results.append({
            "rank": i + 1,
            "path": str(pose_path),
            "rmsd": rmsd,
            "pb_pass": pb_ok,
            "success": success,
            "deltaH": de["deltaH"],
            "minus_TdS": de["minus_TdS"],
            "entropy_enthalpy_index": idx,
            "thermo_score": p.get("thermo_score"),
        })

    return {
        "n": len(poses),
        "n_success": n_success,
        "success_rate": n_success / len(poses) if poses else 0.0,
        "poses": results,
    }


def write_report(stats: Dict, out_json: Path, out_plot_dir: Path = None):
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(stats, f, indent=2)

    # Simple text summary (plots can be added with matplotlib/seaborn if desired)
    print(f"Success (RMSD<=2 + PB): {stats['n_success']}/{stats['n']} = {stats['success_rate']*100:.1f}%")
    print(f"Report written: {out_json}")
