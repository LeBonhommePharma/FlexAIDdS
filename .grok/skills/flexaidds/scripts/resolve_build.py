#!/usr/bin/env python3
"""
resolve_build.py — Canonical FlexAIDdS build discovery with SHA pinning.

Rejects stale FLEXAIDDS_BUILD paths (missing binaries, outdated vs sources).
Pins the active engine SHA256 for reproducibility across agents and launches.

Usage:
  python3 resolve_build.py --check
  python3 resolve_build.py --json
  python3 resolve_build.py --export-shell
  python3 resolve_build.py --sync-env          # rewrite ~/.flexaidds_env (backup first)
  python3 resolve_build.py --write-pin         # persist resolved build as active pin
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

VERSION = "2026-07"

RUNTIME_DATA_FILES = ("MC_st0r5.2_6.dat", "AMINO.def", "NUCLEOTIDES.def")
SOURCE_FILES = ("LIB/benchmark_datasets.cpp", "LIB/DatasetRunner.cpp")
ENGINE_NAME = "FlexAIDdS"
RUNNER_NAME = "benchmark_datasets"

# Only reject known ephemeral FlexAIDdS bench trees under /tmp — not every /tmp path
# (pytest and CI use /tmp for legitimate ephemeral workspaces).
STALE_TMP_MARKERS = ("/tmp/flexaidds", "/private/tmp/flexaidds")


@dataclass(frozen=True)
class BuildResolution:
    build_dir: str
    engine_path: str
    runner_path: str
    engine_sha256: str
    runner_sha256: str
    git_head: str
    repo_root: str
    source: str
    fresh: bool
    pinned: bool
    rejected: tuple[tuple[str, str], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mtime(path: Path) -> int:
    try:
        return int(path.stat().st_mtime)
    except OSError:
        return 0


def discover_git_root(start: Path) -> Path | None:
    current = start.resolve()
    for _ in range(14):
        if (current / ".git").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    return None


def git_head(repo_root: Path | None) -> str:
    if not repo_root:
        return "unknown"
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=repo_root,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except (subprocess.CalledProcessError, OSError):
        return "unknown"


def source_mtime_max(repo_root: Path | None) -> int:
    if not repo_root:
        return 0
    latest = 0
    for rel in SOURCE_FILES:
        path = repo_root / rel
        if path.is_file():
            latest = max(latest, _mtime(path))
    return latest


def has_runtime_data(build_dir: Path) -> bool:
    return all((build_dir / name).is_file() and (build_dir / name).stat().st_size > 0 for name in RUNTIME_DATA_FILES)


def validate_build_dir(build_dir: Path) -> tuple[bool, str]:
    if not build_dir.is_dir():
        return False, "not a directory"
    engine = build_dir / ENGINE_NAME
    runner = build_dir / RUNNER_NAME
    if not engine.is_file() or not os.access(engine, os.X_OK):
        return False, f"missing or non-executable {ENGINE_NAME}"
    if not runner.is_file() or not os.access(runner, os.X_OK):
        return False, f"missing or non-executable {RUNNER_NAME}"
    if not has_runtime_data(build_dir):
        return False, "missing runtime data (MC_*.dat / *.def)"
    return True, "ok"


def is_fresh(build_dir: Path, repo_root: Path | None) -> bool:
    threshold = source_mtime_max(repo_root)
    if threshold == 0:
        return True
    engine_mt = _mtime(build_dir / ENGINE_NAME)
    runner_mt = _mtime(build_dir / RUNNER_NAME)
    return min(engine_mt, runner_mt) >= threshold


def is_stale_tmp_path(path: Path) -> bool:
    text = str(path.resolve()).lower()
    return any(marker in text for marker in STALE_TMP_MARKERS)


def active_manifest_path() -> Path:
    return Path.home() / ".flexaidds" / "active_build.json"


def load_pin_sha() -> str | None:
    env_pin = os.environ.get("FLEXAIDDS_ENGINE_SHA256", "").strip().lower()
    if env_pin:
        return env_pin
    pin_file = os.environ.get("FLEXAIDDS_PIN_FILE", "").strip()
    paths = []
    if pin_file:
        paths.append(Path(pin_file).expanduser())
    paths.append(active_manifest_path())
    repo = discover_git_root(Path(__file__))
    if repo:
        paths.append(repo / ".flexaidds" / "active_build.json")
    for path in paths:
        if not path.is_file():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        sha = str(data.get("engine_sha256", "")).strip().lower()
        if sha:
            return sha
    return None


def collect_candidates(repo_root: Path | None) -> list[Path]:
    seen: set[str] = set()
    ordered: list[Path] = []

    def add(path: Path | str | None) -> None:
        if not path:
            return
        resolved = Path(path).expanduser().resolve()
        key = str(resolved)
        if key in seen:
            return
        seen.add(key)
        ordered.append(resolved)

    for key in ("FLEXAIDDS_BUILD", "FLEXAIDDS_BUILD_DIR"):
        add(os.environ.get(key, "").strip() or None)

    if repo_root:
        add(repo_root / "build_lto")
        add(repo_root / "build")

    home = Path.home()
    add(home / "Projects" / "FlexAIDdS" / "build_lto")
    add(home / "Projects" / "FlexAIDdS" / "build")
    add(home / "FlexAIDdS" / "build_lto")
    add(home / "FlexAIDdS" / "build")

    return ordered


def resolve_build(
    repo_root: Path | None = None,
    *,
    require_fresh: bool = False,
    pin_sha: str | None = None,
    ignore_existing_pin: bool = False,
) -> BuildResolution:
    if repo_root is None:
        repo_root = discover_git_root(Path(__file__).resolve().parents[4])
    if ignore_existing_pin:
        pin = (pin_sha or "").strip().lower() or None
    else:
        pin = (pin_sha or load_pin_sha() or "").strip().lower() or None

    rejected: list[tuple[str, str]] = []
    valid: list[tuple[Path, bool, bool, str, str]] = []

    for candidate in collect_candidates(repo_root):
        ok, reason = validate_build_dir(candidate)
        if not ok:
            rejected.append((str(candidate), reason))
            continue

        engine = candidate / ENGINE_NAME
        runner = candidate / RUNNER_NAME
        engine_sha = _sha256(engine)
        runner_sha = _sha256(runner)
        fresh = is_fresh(candidate, repo_root)
        matches_pin = pin is None or engine_sha == pin

        if pin and not matches_pin:
            rejected.append((str(candidate), f"engine SHA mismatch (want {pin[:16]}…)"))
            continue

        if require_fresh and not fresh:
            rejected.append((str(candidate), "binaries older than benchmark source"))
            continue

        if is_stale_tmp_path(candidate) and not pin:
            rejected.append((str(candidate), "stale /tmp build tree (unpinned)"))
            continue

        valid.append((candidate, fresh, pin is not None, engine_sha, runner_sha))

    if not valid:
        lines = ["No valid FlexAIDdS build directory found."]
        if pin:
            lines.append(f"Pinned engine SHA256 required: {pin}")
        if rejected:
            lines.append("Rejected candidates:")
            for path, reason in rejected[:12]:
                lines.append(f"  - {path}: {reason}")
        raise SystemExit("\n".join(lines))

    # Prefer: pinned match > fresh > newest runner mtime
    valid.sort(key=lambda item: (item[2], item[1], _mtime(item[0] / RUNNER_NAME)), reverse=True)
    chosen, fresh, pinned, engine_sha, runner_sha = valid[0]

    return BuildResolution(
        build_dir=str(chosen),
        engine_path=str(chosen / ENGINE_NAME),
        runner_path=str(chosen / RUNNER_NAME),
        engine_sha256=engine_sha,
        runner_sha256=runner_sha,
        git_head=git_head(repo_root),
        repo_root=str(repo_root) if repo_root else "unknown",
        source="pin" if pinned else ("fresh-auto" if fresh else "fallback-auto"),
        fresh=fresh,
        pinned=pinned,
        rejected=tuple(rejected),
    )


def write_active_manifest(resolution: BuildResolution) -> Path:
    manifest_dir = Path.home() / ".flexaidds"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    path = active_manifest_path()
    payload = {
        "version": VERSION,
        "resolved_at": datetime.now(timezone.utc).isoformat(),
        **asdict(resolution),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def export_shell(resolution: BuildResolution) -> str:
    build = resolution.build_dir
    lines = [
        f"export FLEXAIDDS_BUILD={_shell_quote(build)}",
        f"export FLEXAIDDS_BINARY={_shell_quote(resolution.engine_path)}",
        f"export FLEXAIDDS_RUNNER={_shell_quote(resolution.runner_path)}",
        f"export FLEXAIDDS_ENGINE_SHA256={_shell_quote(resolution.engine_sha256)}",
        f"export FLEXAIDDS_RUNNER_SHA256={_shell_quote(resolution.runner_sha256)}",
        f"export PATH={_shell_quote(build + os.pathsep + os.environ.get('PATH', ''))}",
    ]
    bust = resolve_posebusters_bin(build)
    if bust:
        lines.insert(-1, f"export FLEXAIDDS_POSEBUSTERS_BIN={_shell_quote(bust)}")
    return "\n".join(lines)


def _shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def resolve_posebusters_bin(build_dir: str) -> str | None:
    """Locate the upstream PoseBusters `bust` CLI.

    claim_ready requires pb_backend == "bust_cli"; when `bust` cannot be found
    the engine falls back to the in-house NativePoseQC suite and claim_ready
    becomes unreachable. Exporting the resolved path here keeps runners from
    depending on whichever PATH the launcher happened to inherit.

    Mirrors the C++ lookup order in LIB/PoseBust/BustCli.cpp:resolve_bust_binary.
    """
    env = os.environ.get("FLEXAIDDS_POSEBUSTERS_BIN", "").strip()
    if env and os.access(env, os.X_OK):
        return env

    found = shutil.which("bust")
    if found:
        return found

    # The repo-local venv is the conventional install site for this project.
    candidates = []
    root = os.environ.get("FLEXAIDDS_ROOT", "").strip()
    if root:
        candidates.append(Path(root) / ".venv-posebusters" / "bin" / "bust")
    # build_dir is normally <repo>/build; walk up to the repo root.
    build_path = Path(build_dir)
    for parent in (build_path.parent, build_path.parent.parent):
        candidates.append(parent / ".venv-posebusters" / "bin" / "bust")

    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return None


def _parse_env_value(text: str, key: str) -> str | None:
    pattern = re.compile(rf"^{re.escape(key)}=(.+)$", re.MULTILINE)
    match = pattern.search(text)
    if not match:
        return None
    raw = match.group(1).strip()
    if raw.startswith('"') and raw.endswith('"'):
        return raw[1:-1]
    if raw.startswith("'") and raw.endswith("'"):
        return raw[1:-1]
    return raw


def sync_flexaidds_env(resolution: BuildResolution) -> Path:
    env_path = Path.home() / ".flexaidds_env"
    preserved: dict[str, str] = {}
    if env_path.is_file():
        old = env_path.read_text(encoding="utf-8")
        for key in (
            "FLEXAIDDS_RESULTS",
            "FLEXAIDDS_ICLOUD",
            "FLEXAIDDS_CACHE",
            "FLEXAIDDS_CACHE_ROOT",
            "FLEXAIDDS_RESULTS_ROOT",
            "OMP_NUM_THREADS",
        ):
            val = _parse_env_value(old, key)
            if val:
                preserved[key] = val
        backup = env_path.with_suffix(f".env.bak.{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        shutil.copy2(env_path, backup)

    results = preserved.get("FLEXAIDDS_RESULTS") or preserved.get("FLEXAIDDS_RESULTS_ROOT")
    if not results:
        # Local-first (AGENTS.md): live results on local APFS; iCloud is thin mirror only.
        local_root = os.environ.get("FLEXAIDDS_LOCAL_ROOT", "").strip()
        if local_root:
            results = str(Path(local_root).expanduser() / "results")
        else:
            results = str(Path.home() / "flexaidds_results" / "results")

    lines = [
        "# FlexAIDdS active environment — generated by resolve_build.py",
        f"# {datetime.now(timezone.utc).isoformat()}",
        f"export FLEXAIDDS_BUILD={_shell_quote(resolution.build_dir)}",
        f"export FLEXAIDDS_BINARY={_shell_quote(resolution.engine_path)}",
        f"export FLEXAIDDS_RUNNER={_shell_quote(resolution.runner_path)}",
        f"export FLEXAIDDS_ENGINE_SHA256={_shell_quote(resolution.engine_sha256)}",
        f"export FLEXAIDDS_RUNNER_SHA256={_shell_quote(resolution.runner_sha256)}",
        f"export FLEXAIDDS_RESULTS={_shell_quote(results)}",
        # claim_ready needs pb_backend == bust_cli; pin the CLI so campaigns do
        # not silently degrade to NativePoseQC when PATH lacks it.
        *(
            [f"export FLEXAIDDS_POSEBUSTERS_BIN={_shell_quote(bust)}"]
            if (bust := resolve_posebusters_bin(resolution.build_dir))
            else []
        ),
        # Quote only the build dir; leave :$PATH unquoted so the shell expands PATH.
        # Quoting '…:$PATH' freezes the literal string and wipes /usr/bin from PATH.
        f"export PATH={_shell_quote(resolution.build_dir)}:$PATH",
        "",
        "# Re-resolve after rebuilds: python3 .grok/skills/flexaidds/scripts/resolve_build.py --sync-env",
    ]
    for key, val in preserved.items():
        if key in {"FLEXAIDDS_RESULTS", "FLEXAIDDS_RESULTS_ROOT"}:
            continue
        lines.insert(-2, f"export {key}={_shell_quote(val)}")

    env_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_active_manifest(resolution)
    return env_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Resolve and pin FlexAIDdS build trees.")
    parser.add_argument("--repo-root", type=Path, help="FlexAIDdS git checkout root")
    parser.add_argument("--check", action="store_true", help="Exit 0 when a valid build resolves")
    parser.add_argument("--json", action="store_true", help="Print resolution JSON")
    parser.add_argument("--export-shell", action="store_true", help="Print shell export statements")
    parser.add_argument("--sync-env", action="store_true", help="Rewrite ~/.flexaidds_env from resolution")
    parser.add_argument("--write-pin", action="store_true", help="Write ~/.flexaidds/active_build.json pin")
    parser.add_argument("--require-fresh", action="store_true", help="Reject binaries older than sources")
    parser.add_argument("--pin-sha", help="Require exact engine SHA256")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    repo_root = args.repo_root
    if repo_root:
        repo_root = repo_root.expanduser().resolve()
    else:
        repo_root = discover_git_root(Path(__file__).resolve().parents[4])

    # --write-pin / --sync-env re-pin: ignore a stale existing pin so rebuilt
    # binaries can become the new pin without manual pin-file surgery.
    re_pin = bool((args.write_pin or args.sync_env) and not args.pin_sha)
    try:
        resolution = resolve_build(
            repo_root,
            require_fresh=args.require_fresh,
            pin_sha=args.pin_sha,
            ignore_existing_pin=re_pin,
        )
    except SystemExit as exc:
        print(exc, file=sys.stderr)
        return 1

    if args.write_pin or args.sync_env:
        write_active_manifest(resolution)

    if args.sync_env:
        path = sync_flexaidds_env(resolution)
        print(f"OK: synced {path}")
        print(f"     engine SHA256: {resolution.engine_sha256}")

    if args.export_shell:
        print(export_shell(resolution))

    if args.json:
        print(json.dumps(asdict(resolution), indent=2, sort_keys=True))

    if args.check and not (args.json or args.export_shell or args.sync_env):
        print("OK: build resolved")
        print(f"  build_dir:    {resolution.build_dir}")
        print(f"  engine_sha256: {resolution.engine_sha256}")
        print(f"  runner_sha256: {resolution.runner_sha256}")
        print(f"  source:       {resolution.source}  fresh={resolution.fresh}")

    if not any((args.check, args.json, args.export_shell, args.sync_env, args.write_pin)):
        print(json.dumps(asdict(resolution), indent=2, sort_keys=True))

    return 0


if __name__ == "__main__":
    sys.exit(main())