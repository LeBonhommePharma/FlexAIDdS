"""Simple config and paths. No heavy abstraction."""

from pathlib import Path
import os

REPO_ROOT = Path(__file__).resolve().parents[2]

# Existing yamls (use yours)
ASTEX_DIVERSE_YAML = REPO_ROOT / "benchmarks" / "datasets" / "astex_diverse.yaml"
ASTEX_NONNATIVE_YAML = REPO_ROOT / "benchmarks" / "datasets" / "astex_nonnative.yaml"

# Data locations (existing in repo)
ASTEX_DIVERSE_DIR = REPO_ROOT / "benchmarks" / "astex_diverse"
ASTEX_NONNATIVE_DIR = REPO_ROOT / "benchmarks" / "astex_nonnative"

# Default output
DEFAULT_DATA_DIR = REPO_ROOT / "benchmarks" / "astex_entropy" / "data"
DEFAULT_RESULTS_DIR = REPO_ROOT / "results" / "astex_entropy"

# External tools (assume in PATH or set env)
VINA_BIN = os.environ.get("VINA_BIN", "vina")
RDOCK_RBCAVITY = os.environ.get("RDOCK_RBCAVITY", "rbcavity")
RDOCK_RBDOCK = os.environ.get("RDOCK_RBDOCK", "rbdock")
BOLTZ_BIN = os.environ.get("BOLTZ_BIN", "boltz")  # or python -m boltz ...

# Your entropy tools (for Shannon, tENCoM, thermo scoring / ITC deltaH deltaS)
# Prefer built binaries or python entry
TENCOM_BIN = os.environ.get("TENCOM_BIN", "tENCoM")
FLEXAIDDS_PYTHON = os.environ.get("FLEXAIDDS_PYTHON", "python -m flexaidds")  # if bindings available

# ITC / thermo data helpers (your existing)
CALIBRATE_ITC = REPO_ROOT / "benchmarks" / "calibrate_itc.py"
FETCH_ITC = REPO_ROOT / "benchmarks" / "fetch_itc_data.py"

# PoseBusters (pip install posebusters) - required for success = RMSD<=2 AND all_passed
# We import it where needed.

# Entropy/enthalpy index definition (simple, for re-ranking and reporting)
# index = (abs(TdS) - abs(dH)) / (abs(dG) + 1e-6)   positive favors entropy-driven
# Used together with your Shannon/tENCoM/thermo stack.
