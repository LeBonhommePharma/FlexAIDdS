import os
import sys
import warnings
from pathlib import Path

from setuptools import Extension, setup

# --- P0: Wire source validator guard (python path) ---
# Runs early so developers doing `pip install -e .` (or equivalent) get
# immediate feedback if they added sources without wiring them.
try:
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
    """Locate Eigen headers. Returns path or None.
    Does not raise so that `pip install` can succeed in pure-Python mode
    (the package has graceful fallbacks when the _core extension is absent).
    """
    # Allow explicit override (very useful for cibuildwheel / CI / conda)
    if os.environ.get("EIGEN_INCLUDE_DIR"):
        p = Path(os.environ["EIGEN_INCLUDE_DIR"])
        if (p / "Eigen" / "Dense").exists():
            return str(p)
        if p.name == "eigen" or p.name.startswith("eigen-"):
            # allow pointing at the extracted dir containing Eigen/
            if (p / "Eigen" / "Dense").exists():
                return str(p)

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

    # 3. Homebrew / common system / extracted locations (cibuildwheel friendly)
    for candidate in (
        Path("/opt/homebrew/include/eigen3"),
        Path("/usr/local/include/eigen3"),
        Path("/usr/include/eigen3"),
        Path("/usr/local/include/eigen"),   # sometimes extracted this way
        Path(os.environ.get("EIGEN_ROOT", "")) if os.environ.get("EIGEN_ROOT") else None,
    ):
        if candidate and (candidate / "Eigen" / "Dense").exists():
            return str(candidate)

    return None


_EIGEN_INCLUDE_DIR = _find_eigen_include_dir()

# Keep in sync with python/CMakeLists.txt and the main C++ build.
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

ext_modules = []
build_core = True

try:
    import pybind11
except ImportError:
    warnings.warn(
        "pybind11 not found. The accelerated _core extension will be skipped. "
        "Install with `pip install pybind11` (or use conda) for full performance. "
        "Pure-Python fallbacks will be used.",
        stacklevel=2,
    )
    build_core = False

if build_core and _EIGEN_INCLUDE_DIR is None:
    warnings.warn(
        "Eigen3 headers not found. The accelerated _core C++ extension will be "
        "skipped.\n"
        "  macOS:   brew install eigen\n"
        "  Linux:   sudo apt install libeigen3-dev\n"
        "  Windows: choco install eigen\n"
        "  Or vendor: git submodule update --init LIB/vendor/eigen\n\n"
        "The flexaidds package will still install and work in pure-Python fallback mode "
        "(results loading, models, pure StatMech, etc.).",
        stacklevel=2,
    )
    build_core = False

if build_core:
    try:
        ext = Extension(
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
                else ["-std=c++26", "-O3", "-ffast-math"]
            ),
        )
        ext_modules = [ext]
    except Exception as exc:  # pragma: no cover
        warnings.warn(f"Failed to prepare _core extension, installing without it: {exc}", stacklevel=2)
        ext_modules = []

# IMPORTANT: most metadata (name, version from dynamic, description, packages,
# package_data, scripts, etc.) lives in pyproject.toml for modern builds.
# We only pass what is needed for the compiled extension here.
setup(
    ext_modules=ext_modules,
)
