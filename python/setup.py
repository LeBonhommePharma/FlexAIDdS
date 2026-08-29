"""Build script for the flexaidds Python package.

The accelerated ``flexaidds._core`` extension is **optional**:
- Built when pybind11 + Eigen + C++ sources are available and compilation succeeds.
- Skipped cleanly (pure-Python fallback) when any of those are missing or the
  compile fails. ``pip install flexaidds`` must always succeed.

Important packaging constraints:
- ``setup()`` Extension sources must be relative paths under ``python/`` (never
  absolute, never ``../``). Modern setuptools rejects absolute paths even when
  assigned later in ``build_ext``.
- Monorepo ``../LIB`` sources are **staged** into the gitignored ``python/LIB/``
  tree at ``build_ext`` time so every source path stays relative to ``setup.py``.
- The custom ``sdist`` command stages the minimal LIB files into the *release
  tree* as ``LIB/`` so out-of-tree builds can still compile ``_core``.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist

# --- P0: Wire source validator guard (python path) ---
# Runs early so developers doing `pip install -e .` (or equivalent) get
# immediate feedback if they added sources without wiring them.
try:
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    from scripts.validate_sources import validate_sources

    # Respect env var for CI / strict local runs. Default is lenient for
    # normal developer `pip install -e` usage.
    strict = os.environ.get("FLEXAIDS_STRICT_SOURCE_VALIDATION", "0").lower() not in (
        "0",
        "false",
        "",
    )
    validate_sources(root=str(repo_root), strict=strict)
except Exception as exc:
    # Never break the build due to the guard during normal development.
    print(f"[source-guard] Warning: validator skipped ({exc})", file=sys.stderr)

ROOT = Path(__file__).resolve().parent

# Relative source paths (under LIB/) required to build flexaidds._core.
# Keep in sync with python/CMakeLists.txt and the main C++ build.
_CORE_LIB_REL_SOURCES = [
    "statmech.cpp",
    "TurboQuant.cpp",
    "UnifiedHardwareDispatch.cpp",
    "hardware_detect.cpp",
    "encom.cpp",
    "tENCoM/tencm.cpp",
    "ShannonThermoStack/ShannonThermoStack.cpp",
    "DiFT/DiFT.cpp",
    "fast_optics.cpp",
]

# Headers / extra files staged into the sdist so out-of-tree builds can compile.
# This is the FLOOR, not the whole set — the rest is derived below. A
# hand-maintained header list silently rots: this one was missing flexaid.h,
# flexaid_exception.h, EnvFlags.h, ThermodynamicEngine.h, log_sum_exp.h,
# metal_eval.h and simd_distance.h, so `pip install flexaidds` could never
# build flexaidds._core on ANY platform or compiler. It always fell back to
# pure Python while reporting only a warning:
#     LIB/ShannonThermoStack/../tENCoM/tencm.h:36:10:
#         fatal error: '../flexaid.h' file not found
# Verified on macOS arm64 / AppleClang 21 from a clean sdist.
_CORE_LIB_REL_HEADERS = [
    "statmech.h",
    "TurboQuant.h",
    "UnifiedHardwareDispatch.h",
    "hardware_detect.h",
    "encom.h",
    "fast_optics.hpp",
    "tENCoM/tencm.h",
    "ShannonThermoStack/ShannonThermoStack.h",
    "DiFT/DiFT.h",
]

# Quoted #include "..." — angle-bracket includes are system/third-party.
_INCLUDE_RE = re.compile(r'^[ \t]*#[ \t]*include[ \t]+"([^"]+)"', re.MULTILINE)

# Include roots handed to the compiler (see _prepare_core_extension): the LIB
# root plus the two subdirs. A header is resolved against the including file's
# directory first, then against each of these — same order the compiler uses.
_CORE_INCLUDE_ROOTS = ("", "tENCoM", "ShannonThermoStack")


def _transitive_local_headers(lib_dir: Path) -> List[str]:
    """Every header under *lib_dir* reachable from the _core sources.

    Derived by walking quoted includes rather than enumerated by hand, so
    adding an include to a core source cannot quietly break the sdist again.
    Includes that escape ``LIB/`` (e.g. ``../src/backends/webgpu/webgpu_eval.h``,
    reached only under a WebGPU build) are ignored: the pip path never enables
    those backends.
    """
    seen: set = set()
    found: set = set()
    stack: List[str] = list(_CORE_LIB_REL_SOURCES) + list(_CORE_LIB_REL_HEADERS)

    # The in-tree binding sources also reach into LIB/, but with monorepo-shaped
    # paths: bindings/bindings_matrix.cpp says #include "../../LIB/atom_typing_256.h",
    # which is correct from python/bindings/ in the repo and one level wrong in
    # the sdist. The compiler still finds it via -I<lib_dir>/tENCoM (…/LIB/tENCoM/
    # ../../LIB/x.h == …/LIB/x.h) — but only if the header was staged at all.
    # Seed from those files so it is.
    binding_seeds: List[Path] = [ROOT / "flexaidds" / "_core.cpp"]
    bindings_dir = ROOT / "bindings"
    if bindings_dir.is_dir():
        binding_seeds.extend(sorted(bindings_dir.glob("*.cpp")))
        binding_seeds.extend(sorted(bindings_dir.glob("*.h")))
    for seed in binding_seeds:
        if not seed.is_file():
            continue
        try:
            text = seed.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for match in _INCLUDE_RE.finditer(text):
            target = match.group(1).replace("\\", "/")
            # Keep only the part after the last LIB/ segment.
            if "LIB/" in target:
                rel = target.rsplit("LIB/", 1)[1]
                if (lib_dir / rel).is_file():
                    found.add(rel)
                    stack.append(rel)

    while stack:
        rel = stack.pop()
        if rel in seen:
            continue
        seen.add(rel)
        path = lib_dir / rel
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for match in _INCLUDE_RE.finditer(text):
            target = match.group(1)
            candidates = [os.path.normpath(os.path.join(os.path.dirname(rel), target))]
            candidates += [
                os.path.normpath(os.path.join(root, target))
                for root in _CORE_INCLUDE_ROOTS
            ]
            for cand in candidates:
                if cand.startswith(".."):
                    continue  # outside LIB/ — not part of the pip build
                if (lib_dir / cand).is_file():
                    found.add(cand)
                    stack.append(cand)
                    break

    return sorted(h for h in found if h.endswith((".h", ".hpp", ".hh", ".inc")))


def _skip_core_requested() -> bool:
    return os.environ.get("FLEXAIDDS_SKIP_CORE", "").lower() in ("1", "true", "yes")


def _resolve_lib_dir() -> Optional[Path]:
    """Locate LIB/ with C++ sources for the optional _core extension.

    Search order:
      1. FLEXAIDDS_LIB_DIR env override
      2. ./LIB           (sdist layout: sources staged beside setup.py)
      3. ../LIB          (monorepo: python/ next to LIB/)
    """
    candidates: List[Path] = []
    env = os.environ.get("FLEXAIDDS_LIB_DIR")
    if env:
        candidates.append(Path(env))
    # Prefer staged in-tree LIB (sdist) over monorepo parent so isolated builds
    # use the files that shipped with the sdist.
    candidates.extend(
        [
            ROOT / "LIB",
            ROOT.parent / "LIB",
        ]
    )
    for candidate in candidates:
        if candidate.is_dir() and (candidate / "statmech.cpp").is_file():
            return candidate.resolve()
    return None


def _find_eigen_include_dir() -> Optional[str]:
    """Locate Eigen headers. Returns path or None."""
    if os.environ.get("EIGEN_INCLUDE_DIR"):
        p = Path(os.environ["EIGEN_INCLUDE_DIR"])
        if (p / "Eigen" / "Dense").exists():
            return str(p)
        if (p.name == "eigen" or p.name.startswith("eigen-")) and (
            p / "Eigen" / "Dense"
        ).exists():
            return str(p)

    lib_dir = _resolve_lib_dir()
    if lib_dir is not None:
        vendored = lib_dir / "vendor" / "eigen"
        if (vendored / "Eigen" / "Dense").exists():
            return str(vendored)

    try:
        import subprocess

        out = subprocess.run(
            ["pkg-config", "--cflags-only-I", "eigen3"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        for token in out.split():
            if token.startswith("-I"):
                candidate = Path(token[2:])
                if (candidate / "Eigen" / "Dense").exists():
                    return str(candidate)
    except Exception:
        pass

    for candidate in (
        Path("/opt/homebrew/include/eigen3"),
        Path("/usr/local/include/eigen3"),
        Path("/usr/include/eigen3"),
        Path("/usr/local/include/eigen"),
        Path(os.environ["EIGEN_ROOT"]) if os.environ.get("EIGEN_ROOT") else None,
    ):
        if candidate and (candidate / "Eigen" / "Dense").exists():
            return str(candidate)

    return None


def _core_sources_available(lib_dir: Path) -> bool:
    for rel in _CORE_LIB_REL_SOURCES:
        if not (lib_dir / rel).is_file():
            return False
    if not (ROOT / "flexaidds" / "_core.cpp").is_file():
        return False
    return True


def _can_attempt_core() -> bool:
    if _skip_core_requested():
        return False
    try:
        import pybind11  # noqa: F401
    except ImportError:
        return False
    lib_dir = _resolve_lib_dir()
    if lib_dir is None or not _core_sources_available(lib_dir):
        return False
    if _find_eigen_include_dir() is None:
        return False
    return True


def _path_is_under_root(path: Path) -> bool:
    """True if *path* is the same as ROOT or a descendant (not via ``..``)."""
    try:
        path.resolve().relative_to(ROOT.resolve())
        return True
    except ValueError:
        return False


def _stage_lib_under_root(lib_dir: Path) -> Optional[Path]:
    """Ensure LIB sources live under ``python/`` as relative paths.

    setuptools rejects absolute paths and ``../`` in ``Extension.sources``.
    When the monorepo has sources at ``../LIB``, we copy the minimal set into
    the gitignored ``python/LIB/`` tree (same layout as the sdist).
    """
    if _path_is_under_root(lib_dir) and _core_sources_available(lib_dir):
        return lib_dir.resolve()

    dest = (ROOT / "LIB").resolve()
    dest.mkdir(parents=True, exist_ok=True)

    # Sources + headers + transitive headers from known subdirs.
    wanted = list(_CORE_LIB_REL_SOURCES) + list(_CORE_LIB_REL_HEADERS)
    # Everything the core sources actually #include, derived not enumerated.
    wanted += _transitive_local_headers(lib_dir)
    for sub, patterns in (
        ("tENCoM", ("*.h", "*.hpp", "*.hh")),
        ("ShannonThermoStack", ("*.h", "*.hpp", "*.hh")),
        ("DiFT", ("*.h", "*.hpp", "*.hh")),
    ):
        src_sub = lib_dir / sub
        if not src_sub.is_dir():
            continue
        for pattern in patterns:
            for src in src_sub.glob(pattern):
                wanted.append(f"{sub}/{src.name}")

    seen = set()
    staged = 0
    for rel in wanted:
        if rel in seen:
            continue
        seen.add(rel)
        src = lib_dir / rel
        if not src.is_file():
            continue
        dst = dest / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        # Copy when missing or source is newer (cheap mtime check).
        if (
            not dst.is_file()
            or src.stat().st_mtime > dst.stat().st_mtime
            or src.stat().st_size != dst.stat().st_size
        ):
            shutil.copy2(src, dst)
        staged += 1

    if not _core_sources_available(dest):
        return None
    return dest


def _relpath_from_root(path: Path) -> str:
    """Return a POSIX-style path relative to setup.py's directory."""
    rel = path.resolve().relative_to(ROOT.resolve())
    return rel.as_posix()


class optional_build_ext(_build_ext):
    """Build extensions, but never fail the whole install if compile breaks.

    Real LIB/ sources (monorepo ``../LIB`` or sdist ``LIB/``) are staged under
    ``python/`` here so every ``Extension.sources`` entry is a relative path
    (required by modern setuptools).

    Extensions that are skipped or fail to compile are removed from
    ``self.extensions`` so setuptools does not try to package a missing
    ``.so`` / ``.pyd`` (which would otherwise fail the wheel build).
    """

    def build_extensions(self) -> None:
        if not self.extensions:
            return
        requested = list(self.extensions)
        kept: List[Extension] = []
        for ext in requested:
            try:
                self.build_extension(ext)
            except Exception as exc:  # pragma: no cover - platform/compiler dependent
                warnings.warn(
                    f"Failed to build extension {ext.name} ({exc}). Skipping.",
                    stacklevel=2,
                )
                continue
            # Only keep extensions that actually produced a binary.
            try:
                built = Path(self.get_ext_fullpath(ext.name))
            except Exception:
                built = None
            if built is not None and built.is_file():
                kept.append(ext)
            else:
                self.announce(
                    f"skipping package of {ext.name} (no binary produced)",
                    level=2,
                )
        self.extensions = kept
        # When nothing was built, clear ext_modules so bdist_wheel emits a pure
        # ``py3-none-any`` wheel. Leaving empty Extension declarations makes
        # cibuildwheel/delocate tag a platform wheel with no binary and fail
        # (macOS: "Failed to find any binary with the required architecture").
        #
        # NOTE: clearing ``ext_modules`` here is necessary but NOT sufficient.
        # ``bdist_wheel.finalize_options()`` latches ``root_is_pure`` from
        # ``distribution.has_ext_modules()`` *before* ``build_ext`` runs, so by
        # the time we get here the wheel tag is already decided. Verified: on
        # gcc-11 (no C++26) with pybind11+Eigen present, the compile fails, the
        # fallback engages, and the wheel still comes out as
        # ``flexaidds-2.0.3-cp310-cp310-linux_aarch64.whl`` — a platform wheel
        # containing no binary, i.e. exactly the artifact this comment says it
        # prevents. ``bdist_wheel_optional_ext`` below closes that gap.
        if not kept:
            self.distribution.ext_modules = []
            if any(getattr(e, "name", "") == "flexaidds._core" for e in requested):
                warnings.warn(
                    "flexaidds._core was not built. Installing pure-Python fallback only. "
                    "Install a C++26 compiler + Eigen + pybind11 for the accelerated path.",
                    stacklevel=2,
                )

    def build_extension(self, ext: Extension) -> None:
        if ext.name == "flexaidds._core":
            if not self._prepare_core_extension(ext):
                self.announce(
                    "skipping flexaidds._core (sources/deps unavailable)",
                    level=2,
                )
                return
        try:
            super().build_extension(ext)
        except Exception as exc:  # pragma: no cover
            warnings.warn(
                f"Failed to build extension {ext.name}: {exc}. Skipping.",
                stacklevel=2,
            )

    def _prepare_core_extension(self, ext: Extension) -> bool:
        if _skip_core_requested():
            warnings.warn(
                "FLEXAIDDS_SKIP_CORE set — installing pure-Python package only.",
                stacklevel=2,
            )
            return False

        try:
            import pybind11
        except ImportError:
            warnings.warn(
                "pybind11 not found. The accelerated _core extension will be skipped.",
                stacklevel=2,
            )
            return False

        lib_src = _resolve_lib_dir()
        if lib_src is None or not _core_sources_available(lib_src):
            warnings.warn(
                "C++ sources for flexaidds._core not found. "
                "Installing pure-Python fallback only.",
                stacklevel=2,
            )
            return False

        lib_dir = _stage_lib_under_root(lib_src)
        if lib_dir is None:
            warnings.warn(
                "Failed to stage LIB sources under python/ for setuptools. "
                "Installing pure-Python fallback only.",
                stacklevel=2,
            )
            return False

        eigen = _find_eigen_include_dir()
        if eigen is None:
            warnings.warn(
                "Eigen3 headers not found. The accelerated _core extension will be "
                "skipped.\n"
                "  macOS:   brew install eigen\n"
                "  Linux:   sudo apt install libeigen3-dev\n"
                "  Or set:  EIGEN_INCLUDE_DIR=/path/to/eigen",
                stacklevel=2,
            )
            return False

        # Relative paths only — absolute / parent paths are rejected by setuptools.
        sources = ["flexaidds/_core.cpp"]
        for rel in _CORE_LIB_REL_SOURCES:
            sources.append(_relpath_from_root(lib_dir / rel))

        matrix_bindings = ROOT / "bindings" / "bindings_matrix.cpp"
        define_macros: List[Tuple[str, Optional[str]]] = list(ext.define_macros or [])
        if not any(k == "FLEXAIDS_HAS_EIGEN" for k, _ in define_macros):
            define_macros.append(("FLEXAIDS_HAS_EIGEN", "1"))
        if matrix_bindings.is_file():
            sources.append("bindings/bindings_matrix.cpp")
            if not any(k == "FLEXAIDS_USE_256_MATRIX" for k, _ in define_macros):
                define_macros.append(("FLEXAIDS_USE_256_MATRIX", "1"))

        if os.name == "nt":
            for macro in (
                ("_CRT_SECURE_NO_WARNINGS", "1"),
                ("_USE_MATH_DEFINES", "1"),
                ("NOMINMAX", "1"),
            ):
                if not any(k == macro[0] for k, _ in define_macros):
                    define_macros.append(macro)

        ext.sources = sources
        # include_dirs may be absolute (headers only; setuptools allows this).
        ext.include_dirs = [
            str(lib_dir),
            str(lib_dir / "tENCoM"),
            str(lib_dir / "ShannonThermoStack"),
            eigen,
            pybind11.get_include(),
        ]
        ext.define_macros = define_macros
        ext.language = "c++"
        if not ext.extra_compile_args:
            ext.extra_compile_args = (
                ["/O2", "/std:c++latest", "/EHsc"]
                if os.name == "nt"
                else ["-std=c++26", "-O3", "-ffast-math"]
            )
        return True


def _make_bdist_wheel_optional_ext():
    """bdist_wheel that re-decides purity *after* the extension build.

    ``bdist_wheel`` computes ``root_is_pure`` in ``finalize_options()``, which
    runs before ``build_ext``. For a package whose extension is genuinely
    optional that ordering is wrong: when the compile fails and we fall back to
    pure Python, the wheel is still stamped with an ABI + platform tag
    (``cp310-cp310-linux_aarch64``) despite containing no compiled object.

    Such a wheel is actively harmful on an index:
      * it is offered only to one interpreter/platform, so other users get
        nothing from it even though the payload is pure Python;
      * ``delocate``/``auditwheel`` reject it ("failed to find any binary with
        the required architecture"), breaking cibuildwheel; and
      * it silently misrepresents what is inside the artifact.

    We therefore run ``build`` first, then recompute the tag.
    """
    base = None
    try:  # setuptools >= 70 vendors it
        from setuptools.command.bdist_wheel import bdist_wheel as base
    except ImportError:
        try:
            from wheel.bdist_wheel import bdist_wheel as base  # type: ignore[no-redef]
        except ImportError:
            return None

    class bdist_wheel_optional_ext(base):  # type: ignore[misc,valid-type]
        def run(self):
            # Force the extension decision to happen before the tag is used.
            self.run_command("build")
            if not self.distribution.ext_modules:
                self.root_is_pure = True
            super().run()

        def get_tag(self):
            if not self.distribution.ext_modules:
                return ("py3", "none", "any")
            return super().get_tag()

    return bdist_wheel_optional_ext


class sdist_with_lib(_sdist):
    """Stage required LIB/ sources into the sdist so wheels can compile _core.

    Only the files needed by ``flexaidds._core`` are copied into the *release
    tree* (``base_dir/LIB``). Never write into the source checkout.
    """

    def make_release_tree(self, base_dir, files):  # type: ignore[no-untyped-def]
        super().make_release_tree(base_dir, files)

        # Prefer monorepo LIB/ so we do not re-stage a leftover python/LIB.
        monorepo_lib = ROOT.parent / "LIB"
        if monorepo_lib.is_dir() and (monorepo_lib / "statmech.cpp").is_file():
            lib_dir = monorepo_lib.resolve()
        else:
            lib_dir = _resolve_lib_dir()

        if lib_dir is None:
            self.announce(
                "sdist: LIB/ not found — sdist will be pure-Python only",
                level=2,
            )
            return

        dest_lib = Path(base_dir) / "LIB"
        staged = 0
        wanted = list(_CORE_LIB_REL_SOURCES) + list(_CORE_LIB_REL_HEADERS)
        # Everything the core sources actually #include, derived not enumerated.
        # Without this the sdist ships sources whose headers are absent, and the
        # extension can never compile from PyPI on any platform.
        wanted += _transitive_local_headers(lib_dir)

        for sub, patterns in (
            ("tENCoM", ("*.h", "*.hpp", "*.hh")),
            ("ShannonThermoStack", ("*.h", "*.hpp", "*.hh")),
            ("DiFT", ("*.h", "*.hpp", "*.hh")),
        ):
            src_sub = lib_dir / sub
            if not src_sub.is_dir():
                continue
            for pattern in patterns:
                for src in src_sub.glob(pattern):
                    wanted.append(f"{sub}/{src.name}")

        seen = set()
        for rel in wanted:
            if rel in seen:
                continue
            seen.add(rel)
            src = lib_dir / rel
            if not src.is_file():
                continue
            dst = dest_lib / rel
            self.mkpath(str(dst.parent))
            shutil.copy2(src, dst)
            staged += 1

        self.announce(f"sdist: staged {staged} LIB files for _core builds", level=2)


def _placeholder_extension_modules() -> List[Extension]:
    """Declare _core with only in-tree relative sources for setup()/sdist.

    Real sources are filled in by ``optional_build_ext`` at compile time.
    """
    if _skip_core_requested():
        warnings.warn(
            "FLEXAIDDS_SKIP_CORE set — installing pure-Python package only.",
            stacklevel=2,
        )
        return []

    if not _can_attempt_core():
        # Emit a single informative warning for the common missing-deps cases.
        if not (ROOT / "flexaidds" / "_core.cpp").is_file():
            return []
        try:
            import pybind11  # noqa: F401
        except ImportError:
            warnings.warn(
                "pybind11 not found. The accelerated _core extension will be skipped. "
                "Pure-Python fallbacks will be used.",
                stacklevel=2,
            )
            return []
        if _resolve_lib_dir() is None:
            warnings.warn(
                "C++ sources for flexaidds._core not found (expected LIB/ next to "
                "python/ or staged as python/LIB/ in the sdist). "
                "Installing pure-Python fallback only.",
                stacklevel=2,
            )
            return []
        if _find_eigen_include_dir() is None:
            warnings.warn(
                "Eigen3 headers not found. The accelerated _core C++ extension will be "
                "skipped. Pure-Python fallbacks will be used.",
                stacklevel=2,
            )
            return []
        return []

    # Placeholder only — absolute/parent paths are forbidden in setup().
    return [
        Extension(
            "flexaidds._core",
            sources=["flexaidds/_core.cpp"],
            language="c++",
        )
    ]


def _read_version() -> str:
    """Parse version from flexaidds/__version__.py without importing the package.

    Importing ``flexaidds`` pulls numpy (via ``__init__.py``). Build isolation
    does not install numpy, so ``import flexaidds.__version__`` fails and older
    tooling falls back to ``UNKNOWN-0.0.0`` wheels.
    """
    version_file = ROOT / "flexaidds" / "__version__.py"
    for line in version_file.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("__version__") and "=" in stripped:
            rhs = stripped.split("=", 1)[1].strip()
            return rhs.strip().strip("\"'")
    return "0.0.0"


# Metadata primarily lives in pyproject.toml. Explicit name/version here are a
# hard fallback for legacy pip/setuptools that do not fully apply PEP 621.
_CMDCLASS = {
    "build_ext": optional_build_ext,
    "sdist": sdist_with_lib,
}

_bdist_wheel_cls = _make_bdist_wheel_optional_ext()
if _bdist_wheel_cls is not None:
    _CMDCLASS["bdist_wheel"] = _bdist_wheel_cls

setup(
    name="flexaidds",
    version=_read_version(),
    ext_modules=_placeholder_extension_modules(),
    cmdclass=_CMDCLASS,
)
