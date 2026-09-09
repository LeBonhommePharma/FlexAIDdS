#!/usr/bin/env python3
"""Bind the landing install how-to to Formula/flexaidds.rb and pyproject.toml.

Reads the shipped files. Does not re-implement the support matrix.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOWTO = ROOT / "docs" / "INSTALL.md"
FORMULA = ROOT / "Formula" / "flexaidds.rb"
PYPROJECT = ROOT / "python" / "pyproject.toml"
SETUP_PY = ROOT / "python" / "setup.py"


def _section(text: str, heading: str) -> str:
    match = re.search(
        rf"^## {re.escape(heading)}\n(.*?)(?=^## |\Z)",
        text,
        flags=re.M | re.S,
    )
    assert match, f"missing how-to section {heading!r}"
    return match.group(1)


def _first_bash_fence(section: str) -> str:
    match = re.search(r"```bash\n(.*?)```", section, flags=re.S)
    assert match, "missing bash fence"
    return match.group(1)


def test_howto_homebrew_first_shot_is_head_and_matches_formula() -> None:
    howto = HOWTO.read_text(encoding="utf-8")
    formula = FORMULA.read_text(encoding="utf-8")
    brew_sec = _section(howto, "Engine: Homebrew (lowest friction on macOS)")
    fence = _first_bash_fence(brew_sec)

    assert "brew tap lebonhommepharma/flexaidds https://github.com/LeBonhommePharma/FlexAIDdS" in fence
    assert "brew trust --formula lebonhommepharma/flexaidds/flexaidds" in fence
    assert "brew install --HEAD lebonhommepharma/flexaidds/flexaidds" in fence
    assert "FlexAIDdS --help" in fence
    assert re.search(r"^brew install lebonhommepharma/flexaidds/flexaidds\s*$", fence, re.M) is None

    head = re.search(r'head "([^"]+)", branch: "([^"]+)"', formula)
    assert head, "formula missing head URL/branch"
    assert head.group(1) == "https://github.com/LeBonhommePharma/FlexAIDdS.git"
    assert head.group(2) == "main"
    url = re.search(r'^  url "([^"]+)"', formula, flags=re.M)
    assert url and "v2.0.3.tar.gz" in url.group(1)
    assert "unless build.head?" in formula
    assert "brew install --HEAD lebonhommepharma/flexaidds/flexaidds" in formula


def test_howto_pip_first_shot_is_not_unpublished_pypi() -> None:
    howto = HOWTO.read_text(encoding="utf-8")
    pyproject = PYPROJECT.read_text(encoding="utf-8")
    pip_sec = _section(howto, "Python: pip")
    fence = _first_bash_fence(pip_sec)

    assert "not on public PyPI" in pip_sec or "unpublished" in pip_sec.lower()
    assert "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python" in fence
    assert "FLEXAIDDS_SKIP_CORE=1" in fence
    assert "flexaidds --help" in fence
    assert re.search(r"^\s*pip install flexaidds\s*$", fence, flags=re.M) is None

    version = re.search(r'^version = "([^"]+)"', pyproject, flags=re.M)
    assert version, "pyproject.toml missing static [project] version"
    assert version.group(1) == "2.0.3"


def test_support_matrix_does_not_checkmark_broken_first_shots() -> None:
    howto = HOWTO.read_text(encoding="utf-8")
    matrix = _section(howto, "Support matrix")
    # Homebrew first-shot is --HEAD, not an unqualified formula ✅.
    assert "brew install --HEAD" in matrix
    assert "v2.0.3" in matrix and "provenance refuse" in matrix
    # PyPI pip is unpublished, not a working first-shot.
    pypi_rows = [
        line for line in matrix.splitlines() if "pip install flexaidds" in line and "PyPI" in line
    ]
    assert pypi_rows, "matrix must have an explicit PyPI pip row"
    for row in pypi_rows:
        assert "✅" not in row
        assert "unpublished" in row
    git_rows = [line for line in matrix.splitlines() if "git+subdirectory" in line]
    assert git_rows and all("✅" in row for row in git_rows)


def _homebrew_assertions_rb() -> Path:
    prefix = subprocess.check_output(["brew", "--prefix"], text=True).strip()
    path = Path(prefix) / "Library" / "Homebrew" / "formula_assertions.rb"
    assert path.is_file(), f"Homebrew assertions missing at {path}"
    return path


def test_formula_shell_output_matches_homebrew_api() -> None:
    """shell_output(cmd, result=0) — 2nd arg is exit status, not a CLI flag.

    Homebrew formula_assertions.rb is the shipped API. A Pathname cmd is
    exec'd with no argv, so `shell_output(bin/"FlexAIDdS", "--help")` never
    passes --help and fails brew test.
    """
    api = _homebrew_assertions_rb().read_text(encoding="utf-8")
    sig = re.search(
        r"def shell_output\(cmd,\s*result\s*=\s*0\)",
        api,
    )
    assert sig, "Homebrew shell_output(cmd, result=0) signature missing"
    assert "params(cmd: T.any(Pathname, String), result: Integer)" in api

    formula = FORMULA.read_text(encoding="utf-8")
    bad = re.search(
        r'shell_output\(\s*bin/"[^"]+"\s*,\s*"[^"]+"\s*\)',
        formula,
    )
    assert bad is None, (
        "Formula/flexaidds.rb passes a CLI flag as shell_output's result= "
        f"status argument: {bad.group(0)!r}"
    )
    assert 'shell_output("#{bin/"FlexAIDdS"} --help")' in formula


def test_setup_py_first_shot_has_no_package_owned_warning_prints() -> None:
    setup = SETUP_PY.read_text(encoding="utf-8")
    assert "[source-guard] Warning" not in setup
    skip_block = setup.split("if _skip_core_requested():", 1)[1][:400]
    assert "warnings.warn" not in skip_block


def test_setup_py_isolated_skip_core_is_silent(tmp_path: Path) -> None:
    """Drive shipped setup.py as an sdist would: no scripts/, SKIP_CORE=1."""
    dest = tmp_path / "python"
    shutil.copytree(
        ROOT / "python",
        dest,
        ignore=shutil.ignore_patterns(
            "dist", "build", "*.egg-info", "__pycache__", "LIB", ".pytest_cache"
        ),
    )
    env = os.environ.copy()
    env["FLEXAIDDS_SKIP_CORE"] = "1"
    env.pop("FORCE_COLOR", None)
    proc = subprocess.run(
        [sys.executable, str(dest / "setup.py"), "--name"],
        cwd=dest,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "flexaidds"
    assert "[source-guard] Warning" not in proc.stderr
    assert "UserWarning" not in proc.stderr
    assert "FLEXAIDDS_SKIP_CORE set" not in proc.stderr
