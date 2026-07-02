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


def _find_eigen_include_dir():
    """Locate Eigen headers, mirroring the fallback chain used by the
    CMake build (cmake/FlexAIDDependencies.cmake): vendored submodule,
    then system package via pkg-config/Homebrew, then common install
    locations. Raises SystemExit with actionable instructions if none
    of these resolve, matching the CMake build's error message.
    """
    # 1. Vendored git submodule
    vendored = LIB_DIR / "vendor" / "eigen"
    if (vendored / "Eigen" / "Dense").exists():
        return str(vendored)

    # 2. pkg-config (works for apt/dnf libeigen3-dev and Homebrew alike)
    try:
        import subprocess
        out = subprocess.run(
            ["pkg-config", "--cflags-only-I", "eigen3"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        for token in out.split():
            if token.startswith("-I"):
                candidate = Path(token[2:])
                if (candidate / "Eigen" / "Dense").exists():
                    return str(candidate)
    except Exception:
        pass

    # 3. Homebrew (Apple Silicon + Intel prefixes), common Linux paths
    for candidate in (
        Path("/opt/homebrew/include/eigen3"),
        Path("/usr/local/include/eigen3"),
        Path("/usr/include/eigen3"),
    ):
        if (candidate / "Eigen" / "Dense").exists():
            return str(candidate)

    raise SystemExit(
        "Eigen3 headers not found. Fix with one of:\n"
        "  git submodule update --init LIB/vendor/eigen   (from repo root)\n"
        "  macOS:   brew install eigen\n"
        "  Linux:   sudo apt install libeigen3-dev\n"
        "  Windows: choco install eigen"
    )


_EIGEN_INCLUDE_DIR = _find_eigen_include_dir()

# Read version from __version__.py without importing the package (which
# would require _core to be built already).
_version_file = ROOT / "flexaidds" / "__version__.py"
_version_match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']',
                           _version_file.read_text(), re.MULTILINE)
_version = _version_match.group(1) if _version_match else "0.0.0"

# Keep in sync with python/CMakeLists.txt — statmech.cpp delegates log-sum-exp
# and other hot paths to UnifiedHardwareDispatch.
_core_sources = [
    "flexaidds/_core.cpp",
    f"{_rel_lib}/statmech.cpp",
    f"{_rel_lib}/TurboQuant.cpp",
    f"{_rel_lib}/UnifiedHardwareDispatch.cpp",
    f"{_rel_lib}/hardware_detect.cpp",
    f"{_rel_lib}/encom.cpp",
    f"{_rel_lib}/tENCoM/tencm.cpp",
    f"{_rel_lib}/ShannonThermoStack/ShannonThermoStack.cpp",
    f"{_rel_lib}/DiFT/DiFT.cpp",
    f"{_rel_lib}/fast_optics.cpp",
]
_core_defs = [("FLEXAIDS_HAS_EIGEN", "1")]

# 256×256 soft contact matrix bindings (added when file exists)
_matrix_bindings = Path("bindings/bindings_matrix.cpp")
if _matrix_bindings.exists():
    _core_sources.append(str(_matrix_bindings))
    _core_defs.append(("FLEXAIDS_USE_256_MATRIX", "1"))

ext_modules = [
    Extension(
        "flexaidds._core",
        sources=_core_sources,
        include_dirs=[
            str(LIB_DIR),
            str(LIB_DIR / "tENCoM"),
            str(LIB_DIR / "ShannonThermoStack"),
            _EIGEN_INCLUDE_DIR,
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
