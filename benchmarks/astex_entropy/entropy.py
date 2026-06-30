"""entropy.py - Wrappers for your existing entropy metrics + ITC/ΔH/ΔS + entropy/enthalpy index.

Uses your Shannon collapse, tENCoM, thermodynamic scoring (via python/flexaidds or binaries).
Focus: decomposition and index for re-ranking and reporting.

For ITC: hooks to your calibrate/fetch if experimental values available for targets.
"""

from pathlib import Path
from typing import Dict, Any, Optional
import subprocess
import json
import math

try:
    import yaml
except ImportError:
    yaml = None

# Your existing (import if available in the repo env)
try:
    import sys
    sys.path.insert(0, str(Path(__file__).parents[2] / "python"))
    from flexaidds import compute_shannon_entropy, run_shannon_thermo_stack  # type: ignore
    HAS_FLEX = True
except Exception:
    HAS_FLEX = False

from .config import TENCOM_BIN, CALIBRATE_ITC, FETCH_ITC


def compute_entropy_thermo(pose_file: Path, receptor_file: Path) -> Dict[str, float]:
    """Run your entropy stack on a pose.

    Returns dict with deltaG, deltaH, deltaS (or proxies), shannon, etc.
    Falls back to stubs if tools not in PATH.
    """
    # Prefer python bindings
    if HAS_FLEX:
        try:
            # Example using your stack (adapt signature to your API)
            res = run_shannon_thermo_stack(str(receptor_file), str(pose_file))
            return {
                "deltaG": getattr(res, "deltaG", 0.0),
                "shannonEntropy": getattr(res, "shannonEntropy", 0.0),
                "torsionalVibEntropy": getattr(res, "torsionalVibEntropy", 0.0),
                "entropyContribution": getattr(res, "entropyContribution", 0.0),
            }
        except Exception:
            pass

    # Fallback: call tENCoM if built
    try:
        out = subprocess.check_output([TENCOM_BIN, str(receptor_file), str(pose_file), "--json"], text=True)
        data = json.loads(out)
        return data
    except Exception:
        pass

    # Very simple stub for demo (replace with real call)
    return {"deltaG": -8.0, "shannonEntropy": 1.2, "torsionalVibEntropy": 0.3, "entropyContribution": 0.4}


def deltaH_deltaS_decomposition(thermo: Dict[str, float], T: float = 298.0) -> Dict[str, float]:
    """Extract or approximate ΔH and ΔS from your thermo result.

    In real use your thermodynamic scoring already gives them.
    Here we derive -TΔS and index from available fields.
    """
    dG = thermo.get("deltaG", -8.0)
    # If you have direct dH from scoring, prefer it. Stub uses contribution.
    # Real path: your FullThermoResult or tENCoM + shannon gives components.
    dH = thermo.get("deltaH", dG * 0.6)  # placeholder - replace with real
    TdS = thermo.get("entropyContribution", (dG - dH))  # approx -TΔS
    dS = TdS / T if T else 0.0

    return {"deltaH": dH, "deltaS": dS, "minus_T_deltaS": TdS, "deltaG": dG}


def entropy_enthalpy_index(dH: float, minus_TdS: float, dG: float) -> float:
    """Entropy/enthalpy index.

    Positive when entropy term dominates.
    Simple: (abs(minus_TdS) - abs(dH)) / (abs(dG) + 1e-6)
    """
    return (abs(minus_TdS) - abs(dH)) / (abs(dG) + 1e-6)


def load_itc_experimental(pdb_id: str) -> Optional[Dict[str, float]]:
    """Hook to your ITC data (fetch/calibrate) for ΔH/ΔS/ΔG experimental.

    Returns None if no data for this target.
    """
    # In real run you would call your scripts or load unified csv
    # Here a stub - extend with real call to calibrate_itc or itc_index.
    return None


def re_rank_poses(poses: list, receptor: Path, itc: Optional[Dict] = None) -> list:
    """Re-rank list of pose dicts using thermo + index.

    Each pose should have at least 'path' and original 'score'.
    Adds 'thermo', 'deltaH', 'minus_TdS', 'index', 'thermo_score'.
    """
    ranked = []
    for p in poses:
        th = compute_entropy_thermo(Path(p["path"]), receptor)
        decomp = deltaH_deltaS_decomposition(th)
        idx = entropy_enthalpy_index(decomp["deltaH"], decomp["minus_TdS"], decomp["deltaG"])

        # Example thermo-aware score: lower is better (can be -dG + weight*index etc.)
        thermo_score = decomp["deltaG"] - 0.5 * idx   # tune as you like

        p2 = dict(p)
        p2.update({
            "thermo": th,
            "deltaH": decomp["deltaH"],
            "minus_TdS": decomp["minus_TdS"],
            "entropy_enthalpy_index": idx,
            "thermo_score": thermo_score,
        })
        ranked.append(p2)

    # Sort by thermo_score (example: lower better)
    ranked.sort(key=lambda x: x["thermo_score"])
    return ranked
