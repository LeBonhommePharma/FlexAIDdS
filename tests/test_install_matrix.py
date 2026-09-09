#!/usr/bin/env python3
"""Bind the landing install how-to to Formula/flexaidds.rb and pyproject.toml.

Reads the shipped files. Does not re-implement the support matrix.
"""

from __future__ import annotations

import ast
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOWTO = ROOT / "docs" / "INSTALL.md"
README = ROOT / "README.md"
FORMULA = ROOT / "Formula" / "flexaidds.rb"
PYPROJECT = ROOT / "python" / "pyproject.toml"
SETUP_PY = ROOT / "python" / "setup.py"
INSTALL_SH = ROOT / "scripts" / "install.sh"
UPDATE_SH = ROOT / "scripts" / "update.sh"
CURL_URL = (
    "https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/"
    "main/scripts/install.sh"
)


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


def test_curl_installer_is_documented_and_dry_run_safe() -> None:
    howto = HOWTO.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")
    script = INSTALL_SH.read_text(encoding="utf-8")
    assert INSTALL_SH.is_file()
    assert CURL_URL in howto
    assert CURL_URL in readme
    assert "brew install --HEAD" in script
    assert "subdirectory=python" in script
    assert "FLEXAIDDS_SKIP_CORE=1" in script
    # Must not present unpublished PyPI or refused stable brew as first-shot.
    assert not re.search(r"^\s*pip install flexaidds\s*$", script, flags=re.M)
    assert "brew install --HEAD" in script
    assert not re.search(
        r"brew install (?!--HEAD)lebonhommepharma/flexaidds/flexaidds",
        script,
    )

    proc = subprocess.run(
        ["bash", str(INSTALL_SH), "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    out = proc.stdout
    assert "brew install --HEAD" in out
    assert "subdirectory=python" in out
    assert "FLEXAIDDS_SKIP_CORE=1" in out
    assert "pip install flexaidds\n" not in out
    assert "done." in out


def test_unified_updater_dry_run() -> None:
    assert UPDATE_SH.is_file()
    howto = HOWTO.read_text(encoding="utf-8")
    assert "scripts/update.sh" in howto
    assert "python -m flexaidds --self-update" in howto
    proc = subprocess.run(
        ["bash", str(UPDATE_SH), "--dry-run"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr + proc.stdout
    out = proc.stdout
    assert "brew upgrade --fetch-HEAD" in out or "engine:" in out
    assert "pip install flexaidds\n" not in out or "unpublished" in out
    assert "done." in out
    script = UPDATE_SH.read_text(encoding="utf-8")
    assert "subdirectory=python" in script
    assert "brew install --HEAD" in script


def test_installer_uv_pipx_release_dry_runs() -> None:
    for extra in (["--python-only", "--uv"], ["--python-only", "--pipx"],
                  ["--from-release", "latest"]):
        proc = subprocess.run(
            ["bash", str(INSTALL_SH), "--dry-run", *extra],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, extra + [proc.stderr, proc.stdout]
    uv = subprocess.run(
        ["bash", str(INSTALL_SH), "--dry-run", "--python-only", "--uv"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    assert "uv tool install" in uv.stdout
    assert "subdirectory=python" in uv.stdout
    pipx = subprocess.run(
        ["bash", str(INSTALL_SH), "--dry-run", "--python-only", "--pipx"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    assert "pipx install" in pipx.stdout
    rel = subprocess.run(
        ["bash", str(INSTALL_SH), "--dry-run", "--from-release", "v2.2.0"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    assert "SHA256SUMS.txt" in rel.stdout


def test_from_archive_rejects_bad_hash_and_accepts_good(tmp_path: Path) -> None:
    dist = tmp_path / "dist"
    dist.mkdir()
    binary = dist / "FlexAIDdS"
    binary.write_text("#!/bin/sh\necho Usage: stub\n", encoding="utf-8")
    binary.chmod(0o755)
    archive = tmp_path / "flexaidds-macos.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(binary, arcname="dist/FlexAIDdS")
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    prefix = tmp_path / "prefix"

    bad = subprocess.run(
        ["bash", str(INSTALL_SH), "--from-archive", str(archive),
         "--sha256", "0" * 64, "--prefix", str(prefix)],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )
    assert bad.returncode != 0
    assert "SHA-256 mismatch" in bad.stderr + bad.stdout

    good = subprocess.run(
        ["bash", str(INSTALL_SH), "--from-archive", str(archive),
         "--sha256", digest, "--prefix", str(prefix)],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )
    assert good.returncode == 0, good.stderr + good.stdout
    installed = prefix / "bin" / "FlexAIDdS"
    assert installed.exists()


def test_extra_install_fronts_exist_and_parse() -> None:
    files = {
        ROOT / ".github" / "actions" / "setup-flexaidds" / "action.yml",
        ROOT / ".github" / "workflows" / "ghcr.yml",
        ROOT / ".github" / "workflows" / "homebrew-bottle.yml",
        ROOT / ".devcontainer" / "devcontainer.json",
        ROOT / ".devcontainer" / "Dockerfile",
        ROOT / "packaging" / "spack" / "package.py",
        ROOT / "packaging" / "easybuild" / "flexaidds-2.0.3.eb",
        ROOT / "packaging" / "modulefiles" / "flexaidds.lua",
        ROOT / "conda" / "conda-forge.md",
    }
    for path in files:
        assert path.is_file(), path
    action = (
        ROOT / ".github" / "actions" / "setup-flexaidds" / "action.yml"
    ).read_text(encoding="utf-8")
    assert "using: composite" in action
    assert "pip install" in action
    json.loads((ROOT / ".devcontainer" / "devcontainer.json").read_text(encoding="utf-8"))
    ast.parse((ROOT / "packaging" / "spack" / "package.py").read_text(encoding="utf-8"))
    ast.parse(
        (ROOT / "packaging" / "easybuild" / "flexaidds-2.0.3.eb").read_text(encoding="utf-8")
    )
    formula = FORMULA.read_text(encoding="utf-8")
    assert re.search(r"^\s*bottle do\b", formula, flags=re.M) is None
    howto = HOWTO.read_text(encoding="utf-8")
    assert "wget -qO-" in howto
    assert "uv tool install" in howto
    assert "pipx install" in howto
    assert "ghcr.io/lebonhommepharma/flexaidds" in howto
    release = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")
    assert "SHA256SUMS.txt" in release
    ver_py = (ROOT / "python" / "flexaidds" / "__version__.py").read_text(encoding="utf-8")
    ver = re.search(r'__version__\s*=\s*"([^"]+)"', ver_py)
    assert ver, "python/flexaidds/__version__.py missing __version__"
    meta = (ROOT / "conda" / "meta.yaml").read_text(encoding="utf-8")
    assert "load_file_regex" in meta
    assert "../python/flexaidds/__version__.py" in meta
    assert f'version = "{ver.group(1)}"' not in meta
    eb = (ROOT / "packaging" / "easybuild" / "flexaidds-2.0.3.eb").read_text(encoding="utf-8")
    assert f"version = '{ver.group(1)}'" in eb
    formula_sha = re.search(r'^  sha256 "([0-9a-f]{64})"', formula, flags=re.M)
    spack = (ROOT / "packaging" / "spack" / "package.py").read_text(encoding="utf-8")
    assert formula_sha and formula_sha.group(1) in spack
    assert formula_sha.group(1) in eb


def _pages_install_section() -> str:
    text = (ROOT / "site" / "FlexAIDdS" / "sections.jsx").read_text(encoding="utf-8")
    match = re.search(
        r"function InstallSection\(\) \{.*?\n\}\n\n// ─── BENCHMARKS",
        text,
        flags=re.S,
    )
    assert match, "InstallSection missing from site/FlexAIDdS/sections.jsx"
    return match.group(0)


def test_pages_install_shows_first_shot_and_extra_fronts() -> None:
    """Drive the shipped Pages #install source, not a restated command list."""
    section = _pages_install_section()
    required = (
        "scripts/install.sh",
        'cmd">wget</span> -qO-',
        "brew install</span> --HEAD lebonhommepharma/flexaidds/flexaidds",
        "subdirectory=python",
        "uv tool install",
        "pipx install",
        "scripts/update.sh",
        "python</span> -m flexaidds --self-update",
        "Two separate things",
        "flexaidds</span> Python",
        "--from-release latest",
        "SHA256SUMS.txt",
        "ghcr.io/lebonhommepharma/flexaidds",
        "Dockerfile.locked",
        "environment.yml",
        "setup-flexaidds",
        ".devcontainer/",
        "packaging/spack/package.py",
        "packaging/easybuild/flexaidds-2.0.3.eb",
        "packaging/modulefiles",
        "not on PyPI",
        "no bottle do",
        "not a conda-forge feedstock",
        "install-tabs",
        "code-box",
        "cmake-table",
        "binding-blurb",
    )
    missing = [token for token in required if token not in section]
    assert missing == [], f"Pages #install missing tokens: {missing}"
    assert 'maxWidth: "896px"' in section
    assert re.search(r"pip install flexaidds(?![^\n]*subdirectory=python)", section) is None
    assert "pip install flexaidds" not in section
    assert re.search(
        r"brew install(?![\s\S]{0,20}--HEAD)[\s\S]{0,80}lebonhommepharma/flexaidds/flexaidds",
        section,
    ) is None
    assert "conda install -c conda-forge flexaidds" not in section
    for tab in ("curl", "Homebrew", "Python", "Docker", "conda", "CI", "HPC", "CMake"):
        assert tab in section


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
