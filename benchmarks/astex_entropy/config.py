from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml


PACKAGE_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = PACKAGE_DIR.parents[1]
DEFAULT_CONFIG_PATH = PACKAGE_DIR / "config.yaml"
LIVE_REPO_FALLBACK = Path("/Users/lp.more/Projects/FlexAIDdS")
ICLOUD_BENCHMARKS = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"


VALID_BUILD_PROFILES = {"lto", "metal"}


def _build_dir_from_profile(raw: dict[str, Any]) -> tuple[str, str]:
    profile = str(raw.get("build_profile", "lto")).strip().lower()
    if profile not in VALID_BUILD_PROFILES:
        raise ValueError(
            f"Unknown build_profile: {profile!r}; expected one of {sorted(VALID_BUILD_PROFILES)}"
        )
    return profile, f"build_{profile}"


def _repo_root_from_config(raw: dict[str, Any]) -> Path:
    configured = raw.get("repo_root") or os.environ.get("FLEXAIDDS_REPO_ROOT")
    if configured:
        path = Path(os.path.expandvars(str(configured))).expanduser()
        if path.exists():
            return path

    if (WORKSPACE_ROOT / "benchmarks" / "datasets").exists():
        return WORKSPACE_ROOT
    return LIVE_REPO_FALLBACK


def _expand(value: Any, context: dict[str, str]) -> Any:
    class _KeepUnknown(dict[str, str]):
        def __missing__(self, key: str) -> str:
            return "{" + key + "}"

    if isinstance(value, str):
        return os.path.expandvars(value.format_map(_KeepUnknown(context)))
    if isinstance(value, list):
        return [_expand(v, context) for v in value]
    if isinstance(value, dict):
        return {k: _expand(v, context) for k, v in value.items()}
    return value


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    path = Path(config_path) if config_path else DEFAULT_CONFIG_PATH
    with path.open() as fh:
        raw = yaml.safe_load(fh) or {}

    repo_root = _repo_root_from_config(raw)
    build_profile, build_dir = _build_dir_from_profile(raw)
    context = {
        "repo_root": str(repo_root),
        "workspace_root": str(WORKSPACE_ROOT),
        "package_dir": str(PACKAGE_DIR),
        "icloud_benchmarks": str(ICLOUD_BENCHMARKS),
        "build_profile": build_profile,
        "build_dir": build_dir,
    }
    raw_work_dir = _expand(raw.get("work_dir", "results/astex_entropy"), context)
    work_dir = (WORKSPACE_ROOT / str(raw_work_dir)).resolve()
    context["work_dir"] = str(work_dir)
    cfg = _expand(raw, context)
    cfg["repo_root"] = str(repo_root)
    cfg["workspace_root"] = str(WORKSPACE_ROOT)
    cfg["config_path"] = str(path.resolve())
    cfg["work_dir"] = str(work_dir)
    cfg["build_profile"] = build_profile
    cfg["build_dir"] = build_dir
    return cfg


def path_from_config(value: str | Path) -> Path:
    return Path(str(value)).expanduser().resolve()
