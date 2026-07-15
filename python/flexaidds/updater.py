"""Self-updater for FlexAID∆S.

Checks GitHub Releases API for new versions and optionally downloads/installs.
Uses only the standard library.
"""

from __future__ import annotations

import dataclasses
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Callable, Optional

from .__version__ import __version__, __version_info__

# Pre-release ordering for version comparison
_PRERELEASE_ORDER = {"alpha": 0, "beta": 1, "rc": 2}


@dataclasses.dataclass
class AssetInfo:
    """A downloadable release asset."""

    name: str
    url: str
    size: int
    content_type: str


@dataclasses.dataclass
class UpdateInfo:
    """Result of an update check."""

    latest_version: str
    current_version: str
    update_available: bool
    release_url: str
    release_notes: str
    assets: list[AssetInfo]
    published_at: str


def get_current_version() -> tuple[int, ...]:
    """Return numeric version tuple from __version_info__."""
    return tuple(v for v in __version_info__ if isinstance(v, int))


def _parse_version(tag: str) -> tuple[tuple[int, ...], int]:
    """Parse a version tag into (numeric_tuple, prerelease_rank).

    Examples:
        "v1.2.3"       → ((1, 2, 3), 3)       # release
        "1.0.0-alpha"  → ((1, 0, 0), 0)       # alpha
        "v2.1.0-rc.1"  → ((2, 1, 0), 2)       # release candidate
    """
    tag = tag.lstrip("v")

    # Split on hyphen to separate version from prerelease
    parts = tag.split("-", 1)
    nums = tuple(int(x) for x in re.findall(r"\d+", parts[0]))

    prerelease_rank = 3  # release (no suffix)
    if len(parts) > 1:
        suffix = parts[1].lower()
        for label, rank in _PRERELEASE_ORDER.items():
            if suffix.startswith(label):
                prerelease_rank = rank
                break

    return nums, prerelease_rank


def _version_newer(latest_tag: str, current_info: tuple) -> bool:
    """Return True if latest_tag represents a newer version than current_info."""
    latest_nums, latest_rank = _parse_version(latest_tag)

    # Build current comparable: first 3 numeric parts only (major, minor, patch)
    current_nums = tuple(v for v in current_info[:3] if isinstance(v, int))
    current_pre = 3  # default: release
    for v in current_info:
        if isinstance(v, str) and v in _PRERELEASE_ORDER:
            current_pre = _PRERELEASE_ORDER[v]
            break

    return (latest_nums, latest_rank) > (current_nums, current_pre)


def check_for_updates(repo: str = "LeBonhommePharma/FlexAIDdS") -> Optional[UpdateInfo]:
    """Check GitHub Releases for a newer version.

    Returns UpdateInfo if the API call succeeds, None on failure.
    """
    url = f"https://api.github.com/repos/{repo}/releases/latest"
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("User-Agent", "FlexAIDdS-updater")

    token = os.environ.get("GITHUB_TOKEN") or os.environ.get("GH_TOKEN")
    if token:
        req.add_header("Authorization", f"Bearer {token}")

    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
    except (urllib.error.URLError, urllib.error.HTTPError, OSError):
        return None

    tag = data.get("tag_name", "")
    if not tag:
        return None

    assets = [
        AssetInfo(
            name=a["name"],
            url=a["browser_download_url"],
            size=a["size"],
            content_type=a.get("content_type", ""),
        )
        for a in data.get("assets", [])
    ]

    return UpdateInfo(
        latest_version=tag,
        current_version=__version__,
        update_available=_version_newer(tag, __version_info__),
        release_url=data.get("html_url", ""),
        release_notes=data.get("body", ""),
        assets=assets,
        published_at=data.get("published_at", ""),
    )


def download_asset(
    asset: AssetInfo,
    dest_dir: str,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> Path:
    """Download a release asset to dest_dir.

    Args:
        asset: The asset to download.
        dest_dir: Directory to save the file.
        progress_callback: Optional callback(bytes_downloaded, total_bytes).

    Returns:
        Path to the downloaded file.
    """
    dest = Path(dest_dir) / asset.name
    req = urllib.request.Request(asset.url)
    req.add_header("User-Agent", "FlexAIDdS-updater")

    with urllib.request.urlopen(req, timeout=120) as resp:
        total = asset.size
        downloaded = 0
        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(8192)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if progress_callback:
                    progress_callback(downloaded, total)

    return dest


# Git install target used when the package is not (yet) on the default PyPI index.
_GIT_INSTALL_SPEC = (
    "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
)


def update_pip(version: str = "latest", *, index_url: Optional[str] = None) -> int:
    """Update the flexaidds Python package via pip.

    Args:
        version: Version to install, or "latest" for the newest.
        index_url: Optional extra index URL (e.g. TestPyPI). When the package
            is not yet on the default index, callers can pass TestPyPI or a
            direct VCS URL via ``version`` as ``git+https://...``.

    Returns:
        pip exit code.

    Notes:
        When ``version == "latest"`` and the default PyPI index does not yet
        host ``flexaidds``, this falls back to a direct GitHub VCS install so
        ``python -m flexaidds --self-update`` works before the first publish.
    """
    cmd = [sys.executable, "-m", "pip", "install", "--upgrade"]
    if index_url:
        cmd.extend(["--index-url", index_url])
    if version == "latest":
        cmd.append("flexaidds")
    elif version.startswith("git+") or version.startswith("https://") or version.startswith("http://"):
        cmd.append(version)
    else:
        cmd.append(f"flexaidds=={version.lstrip('v')}")

    result = subprocess.run(cmd, capture_output=False)
    # Pre-PyPI / missing-index fallback: try the GitHub VCS URL once.
    if (
        result.returncode != 0
        and version == "latest"
        and not index_url
    ):
        fallback = [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--upgrade",
            "--force-reinstall",
            _GIT_INSTALL_SPEC,
        ]
        result = subprocess.run(fallback, capture_output=False)
    return result.returncode


def select_asset_for_platform(assets: list[AssetInfo]) -> Optional[AssetInfo]:
    """Select the most appropriate binary asset for the current platform."""
    system = platform.system().lower()
    machine = platform.machine().lower()

    # Normalize architecture names
    arch_aliases = {
        "x86_64": ["x86_64", "amd64", "x64"],
        "aarch64": ["aarch64", "arm64"],
    }
    arch_names = arch_aliases.get(machine, [machine])

    for asset in assets:
        name_lower = asset.name.lower()
        if system in name_lower:
            if any(arch in name_lower for arch in arch_names):
                return asset

    # Fallback: match just the OS
    for asset in assets:
        name_lower = asset.name.lower()
        if system in name_lower:
            return asset

    return None
