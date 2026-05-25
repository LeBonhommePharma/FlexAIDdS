import os
import re
import sys
from pathlib import Path

from setuptools import Extension, setup

try:
    import pybind11
except ImportError as exc:
    raise SystemExit(
        "pybind11 is required to build flexaidds. "
        "Install it with `pip install pybind11` and retry."
    ) from exc

# --- P0: Wire source validator guard (python path) ---
# Runs early so developers doing `pip install -e .` (or equivalent) get
# immediate feedback if they added sources without wiring them.
try:
    import sys
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    from scripts.validate_sources import validate_sources

    # Respect env var for CI / strict local runs. Default is lenient for
    # normal developer `pip install -e` usage.
    strict = os.environ.get("FLEXAIDS_STRICT_SOURCE_VALIDATION", "0").lower() not in ("0", "false", "")
    validate_sources(root=str(repo_root), strict=strict)
except Exception as exc:
    # Never break the build due to the guard during normal development.
    print(f"[source-guard] Warning: validator skipped ({exc})", file=sys.stderr)

ROOT = Path(__file__).resolve().parent
LIB_DIR = ROOT.parent / "LIB"

# setuptools requires /-separated paths relative to setup.py directory
_rel_lib = "../LIB"

# Read version from __version__.py without importing the package (which
# would require _core to be built already).
_version_file = ROOT / "flexaidds" / "__version__.py"
_version_match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']',
                           _version_file.read_text(), re.MULTILINE)
_version = _version_match.group(1) if _version_match else "0.0.0"

_core_sources = [
    "flexaidds/_core.cpp",
    f"{_rel_lib}/statmech.cpp",
    f"{_rel_lib}/encom.cpp",
    f"{_rel_lib}/tENCoM/tencm.cpp",
    f"{_rel_lib}/ShannonThermoStack/ShannonThermoStack.cpp",
    f"{_rel_lib}/DiFT/DiFT.cpp",
]
_core_defs = []

# 256×256 soft contact matrix bindings (added when file exists)
_matrix_bindings = Path("bindings/bindings_matrix.cpp")
if _matrix_bindings.exists():
    _core_sources.append(str(_matrix_bindings))
    _core_defs.append(("FLEXAIDS_USE_256_MATRIX", "1"))

ext_modules = [
    Extension(
        "flexaidds._core",
        sources=_core_sources + [f"{_rel_lib}/fast_optics.cpp"],
        include_dirs=[
            str(LIB_DIR),
            str(LIB_DIR / "tENCoM"),
            str(LIB_DIR / "ShannonThermoStack"),
            pybind11.get_include(),
        ],
        define_macros=_core_defs + (
            [
                ("_CRT_SECURE_NO_WARNINGS", "1"),
                ("_USE_MATH_DEFINES", "1"),
                ("NOMINMAX", "1"),
            ]
            if os.name == "nt"
            else []
        ),
        language="c++",
        extra_compile_args=(
            ["/O2", "/std:c++latest", "/EHsc"]
            if os.name == "nt"
            else ["-std=c++26", "-O3"]
        ),
    ),
]

setup(
    name="flexaidds",
    version=_version,
    description="Python bindings for the FlexAID∆S thermodynamic core",
    author="Louis-Philippe Morency",
    packages=["flexaidds", "flexaidds.dataset_runner"],
    package_dir={"": "."},
    package_data={
        "flexaidds": ["py.typed"],
        "flexaidds.dataset_runner": ["datasets/*.yaml"],
    },
    ext_modules=ext_modules,
)
