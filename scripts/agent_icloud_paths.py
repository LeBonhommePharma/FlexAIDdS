#!/usr/bin/env python3
"""agent_icloud_paths — pure path mapping for agent-home local↔iCloud backup/restore.

Hang-safe policy (mirrors AGENTS.md / docs/ICLOUD_BENCHMARK_STORAGE.md):
  - Live agent homes stay on local APFS under $HOME (never CloudDocs symlinks).
  - Durable mirror lands under $FLEXAIDDS_ICLOUD/agent_homes/ (or archive batches).
  - Default backup is thin/selective; large reinstallable runtimes are excluded.
  - No recursive find over Mobile Documents — mapping only.

Usage (library + CLI):
  python3 scripts/agent_icloud_paths.py --print-map
  python3 scripts/agent_icloud_paths.py --print-restore-cmds --archive-batch PATH
  python3 scripts/agent_icloud_paths.py --print-excludes claude_app

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

# Subdirectory under FLEXAIDDS_ICLOUD for the live agent-home mirror (not archive).
AGENT_HOMES_SUBDIR = "agent_homes"
HOME_DOTS_SUBDIR = "home_dots"
ARCHIVED_FROM_SSD_SUBDIR = "archived_from_ssd"

# Remote directory names match the 2026-07-25 archive_batch layout so restore
# and backup share one naming contract.
REMOTE_DOT_CLAUDE = "dot_claude"
REMOTE_APP_SUPPORT_CLAUDE = "Application_Support_Claude"
REMOTE_DOT_CLAUDE_SCIENCE = "dot_claude_science"
REMOTE_DOT_CODEX = "dot_codex"
REMOTE_DOT_GROK = "dot_grok"

# Named paths the user asked about (must be PRESENT or explicit SKIP).
NAMED_HOME_PATHS: tuple[str, ...] = (
    ".claude",
    ".codex",
    ".grok",
    ".venv",
    ".env",
)

# Meaningful home-dot files/dirs to mirror under home_dots/ (not ephemeral noise).
# Agent trees are handled separately under agent_homes/.
HOME_DOT_ALLOWLIST: tuple[str, ...] = (
    ".cursor",
    ".flexaidds_env",
    ".flexaidds_env.bak",
    ".gitconfig",
    ".zshrc",
    ".zprofile",
    ".bash_profile",
    ".alias_profile",
    ".npmrc",
    ".condarc",
    ".nanorc",
    ".claude.json",
    ".claude.json.backup",
    ".tcshrc",
    ".xonshrc",
)

# Ephemeral / skip prefixes under $HOME (never required durable mirrors).
HOME_DOT_SKIP_PREFIXES: tuple[str, ...] = (
    ".zcompdump",
    ".DS_Store",
    ".CFUserTextEncoding",
    ".bash_history",
    ".zsh_history",
    ".zhistory",
    ".lesshst",
    ".python_history",
    ".viminfo",
    ".wget-hsts",
    ".pymoltimestamp",
    ".flexaidds_campaign.pid",
    ".flexaidds_campaign.pids",
    ".passport_contact_sheet",
    ".ElvUI",
    ".elvui",
)

# Regenerable cache-like subpaths under live agent homes. Free only after the
# durable config mirror exists AND the path is listed as freeable (not the
# whole live home). These are intentionally NOT required on iCloud (thin
# backup excludes them); free is allowed only when the parent agent mirror
# verified OK and --free-regenerable is used, or when archive already holds
# an equivalent durable tree for that agent.
REGENERABLE_FREE_CANDIDATES: tuple[tuple[str, str], ...] = (
    # (agent_id, relative path under local home of that agent)
    ("claude", "cache"),
    ("claude", "telemetry"),
    ("claude", "debug"),
    ("claude", "shell-snapshots"),
    ("claude", "paste-cache"),
    ("claude", "file-history"),
    ("codex", "cache"),
    ("codex", "log"),
    ("codex", "node_repl"),
    ("codex", "browser"),
    ("codex", "attachments"),
    ("codex", "generated_images"),
    ("codex", "tmp"),
    ("grok", "logs"),
    ("grok", "marketplace-cache"),
    ("grok", "downloads"),
    ("grok", "upload_queue"),
    ("grok", "vendor"),
    ("grok", "bin"),
    ("grok", "migration"),
    ("claude_app", "Cache"),
    ("claude_app", "Code Cache"),
    ("claude_app", "GPUCache"),
    ("claude_app", "Crashpad"),
    ("claude_app", "DawnGraphiteCache"),
    ("claude_app", "DawnWebGPUCache"),
    ("claude_app", "blob_storage"),
    ("claude_app", "logs"),
)

# Thin excludes (rsync --exclude patterns). Default backup omits caches,
# reinstallable runtimes, and huge VM/blob trees. Secrets stay local unless
# the operator opts into --full (still never force-commits secrets to git).
DEFAULT_EXCLUDES: dict[str, tuple[str, ...]] = {
    "claude": (
        "cache/",
        "telemetry/",
        "shell-snapshots/",
        "debug/",
        "file-history/",
        "paste-cache/",
        "stats-cache.json",
        "gh-pr-status-cache.json",
        "mcp-needs-auth-cache.json",
        "*.lock",
    ),
    "claude_app": (
        "Cache/",
        "Code Cache/",
        "GPUCache/",
        "DawnGraphiteCache/",
        "DawnWebGPUCache/",
        "blob_storage/",
        "Crashpad/",
        "vm_bundles/",
        "Partitions/",
        "IndexedDB/",
        "File System/",
        "Service Worker/",
        "logs/",
        "Network Persistent State",
        "Cookies",
        "Cookies-journal",
        "DIPS",
        "DIPS-wal",
        "Local Storage/",
        "Session Storage/",
        "Shared Dictionary/",
        "SharedStorage/",
        "VideoDecodeStats/",
        "WebStorage/",
        "*.log",
    ),
    "claude_science": (
        "conda/",
        "runtime/",
        "r-libs/",
        "bin/",
        "sbx-bind-src/",
        "seed-assets/",
        "logs/",
        "*.lock",
    ),
    "codex": (
        "cache/",
        "log/",
        "*.sqlite-shm",
        "*.sqlite-wal",
        "node_repl/",
        "browser/",
        "attachments/",
        "generated_images/",
        "tmp/",
        "*.lock",
    ),
    "grok": (
        "logs/",
        "bin/",
        "vendor/",
        "marketplace-cache/",
        "downloads/",
        "upload_queue/",
        "migration/",
        "models_cache.json",
        "*.lock",
        "*.lock.*",
    ),
}


@dataclass(frozen=True)
class AgentMirror:
    """One local home directory ↔ remote mirror directory pair."""

    agent_id: str
    label: str
    local: Path
    remote_name: str
    excludes: tuple[str, ...] = field(default_factory=tuple)

    def remote(self, dest_root: Path) -> Path:
        return dest_root / self.remote_name


def default_icloud_root(home: Path | None = None) -> Path:
    """Resolve FLEXAIDDS_ICLOUD or the standard CloudDocs FlexAIDdS_benchmarks tree."""
    env = os.environ.get("FLEXAIDDS_ICLOUD")
    if env:
        return Path(env)
    h = home if home is not None else Path.home()
    return h / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"


def agent_homes_mirror_root(icloud: Path | None = None) -> Path:
    """Durable local→iCloud agent mirror root (not archive batches)."""
    root = icloud if icloud is not None else default_icloud_root()
    return root / AGENT_HOMES_SUBDIR


def home_dots_mirror_root(icloud: Path | None = None) -> Path:
    """Durable local→iCloud home-dot mirror root (shell/config files)."""
    root = icloud if icloud is not None else default_icloud_root()
    return root / HOME_DOTS_SUBDIR


def archive_batch_root(icloud: Path | None = None, batch: str | None = None) -> Path:
    """Path to archived_from_ssd[/batch] under iCloud."""
    root = icloud if icloud is not None else default_icloud_root()
    base = root / ARCHIVED_FROM_SSD_SUBDIR
    if batch:
        return base / batch
    return base


@dataclass(frozen=True)
class NamedPathStatus:
    """PRESENT / MISSING / SKIP status for a named or allowlisted home path."""

    name: str
    local: Path
    status: str  # PRESENT | MISSING | SKIP
    reason: str
    kind: str  # file | dir | missing
    remote_name: str | None = None


def inventory_named_home_paths(home: Path | None = None) -> list[NamedPathStatus]:
    """Status for user-named paths (.claude, .codex, .grok, .venv, .env)."""
    h = home if home is not None else Path.home()
    out: list[NamedPathStatus] = []
    for name in NAMED_HOME_PATHS:
        p = h / name
        if not p.exists():
            out.append(
                NamedPathStatus(
                    name=name,
                    local=p,
                    status="MISSING",
                    reason="absent locally; nothing to mirror",
                    kind="missing",
                    remote_name=_remote_name_for_home_dot(name),
                )
            )
            continue
        kind = "dir" if p.is_dir() else "file"
        out.append(
            NamedPathStatus(
                name=name,
                local=p,
                status="PRESENT",
                reason="exists locally",
                kind=kind,
                remote_name=_remote_name_for_home_dot(name),
            )
        )
    return out


def _remote_name_for_home_dot(name: str) -> str | None:
    mapping = {
        ".claude": REMOTE_DOT_CLAUDE,
        ".codex": REMOTE_DOT_CODEX,
        ".grok": REMOTE_DOT_GROK,
        ".venv": "dot_venv",
        ".env": "dot_env",
    }
    return mapping.get(name)


def _should_skip_home_dot(name: str) -> str | None:
    """Return skip reason if name is ephemeral noise, else None."""
    if name in (".", ".."):
        return "dot self"
    if name in {".claude", ".codex", ".grok", ".claude-science"}:
        return "handled by agent_homes (not home_dots)"
    for prefix in HOME_DOT_SKIP_PREFIXES:
        if name == prefix or name.startswith(prefix):
            return f"ephemeral skip prefix {prefix}"
    if name.endswith(".zwc"):
        return "zsh compiled dump"
    return None


def inventory_home_dots(home: Path | None = None) -> list[NamedPathStatus]:
    """Inventory allowlisted home dots + named missing; skip ephemeral noise.

    Does not walk beyond one-level names under $HOME (no recursive find).
    """
    h = home if home is not None else Path.home()
    out: list[NamedPathStatus] = []
    seen: set[str] = set()

    # Named paths first (including MISSING .venv/.env)
    for row in inventory_named_home_paths(home=h):
        out.append(row)
        seen.add(row.name)

    # Allowlisted configs
    for name in HOME_DOT_ALLOWLIST:
        if name in seen:
            continue
        seen.add(name)
        p = h / name
        remote = f"dot_{name[1:]}" if name.startswith(".") else name
        if not p.exists():
            out.append(
                NamedPathStatus(
                    name=name,
                    local=p,
                    status="MISSING",
                    reason="allowlisted but absent",
                    kind="missing",
                    remote_name=remote,
                )
            )
            continue
        kind = "dir" if p.is_dir() else "file"
        out.append(
            NamedPathStatus(
                name=name,
                local=p,
                status="PRESENT",
                reason="allowlisted config",
                kind=kind,
                remote_name=remote,
            )
        )

    # Other one-level dots: classify as SKIP with reason (honest accounting)
    try:
        names = sorted(os.listdir(h))
    except OSError as exc:
        out.append(
            NamedPathStatus(
                name=".",
                local=h,
                status="SKIP",
                reason=f"cannot list home: {exc}",
                kind="missing",
            )
        )
        return out

    for name in names:
        if not name.startswith(".") or name in seen:
            continue
        seen.add(name)
        skip = _should_skip_home_dot(name)
        if skip is None and name not in HOME_DOT_ALLOWLIST:
            # Unknown extra dot — skip by default (not force-upload)
            skip = "not in allowlist (not force-uploaded)"
        if skip is None:
            continue
        p = h / name
        kind = "dir" if p.is_dir() else "file" if p.exists() else "missing"
        out.append(
            NamedPathStatus(
                name=name,
                local=p,
                status="SKIP",
                reason=skip,
                kind=kind,
            )
        )
    return out


def build_home_dot_mirrors(
    home: Path | None = None,
    include_missing: bool = False,
) -> list[AgentMirror]:
    """Mirror specs for PRESENT allowlisted home dots (files or dirs).

    Missing named paths (.venv, .env) are not included unless include_missing
    (they still appear in inventory as MISSING — never invent empty uploads).
    """
    h = home if home is not None else Path.home()
    mirrors: list[AgentMirror] = []
    for name in HOME_DOT_ALLOWLIST:
        p = h / name
        if not p.exists():
            if include_missing:
                mirrors.append(
                    AgentMirror(
                        agent_id=f"home_{name.lstrip('.')}",
                        label=f"home dot {name} (missing)",
                        local=p,
                        remote_name=f"dot_{name[1:]}" if name.startswith(".") else name,
                        excludes=(),
                    )
                )
            continue
        mirrors.append(
            AgentMirror(
                agent_id=f"home_{name.lstrip('.')}",
                label=f"home dot {name}",
                local=p,
                remote_name=f"dot_{name[1:]}" if name.startswith(".") else name,
                excludes=(),
            )
        )
    return mirrors


def freeable_regenerable_paths(
    home: Path | None = None,
    agents: Sequence[str] | None = None,
) -> list[dict[str, str]]:
    """Local regenerable cache paths that may be freed after mirror verify.

    Does not free whole live agent roots. Returns dicts with agent_id, local,
    relative path — caller must verify durable agent mirror exists first.
    """
    mirrors = {
        m.agent_id: m
        for m in build_agent_mirrors(home=home, agents=agents, full=False)
    }
    wanted = set(mirrors.keys()) if agents is None else {
        a.strip().lower().replace("-", "_") for a in agents
    }
    rows: list[dict[str, str]] = []
    for agent_id, rel in REGENERABLE_FREE_CANDIDATES:
        if agent_id not in wanted:
            continue
        m = mirrors.get(agent_id)
        if m is None:
            continue
        local = m.local / rel
        rows.append(
            {
                "agent_id": agent_id,
                "relative": rel,
                "local": str(local),
                "agent_local_root": str(m.local),
                "remote_name": m.remote_name,
            }
        )
    return rows


def build_agent_mirrors(
    home: Path | None = None,
    agents: Sequence[str] | None = None,
    full: bool = False,
) -> list[AgentMirror]:
    """Build ordered mirror specs for the named agents (default: all)."""
    h = home if home is not None else Path.home()
    all_specs: list[AgentMirror] = [
        AgentMirror(
            agent_id="claude",
            label="Claude Code (~/.claude)",
            local=h / ".claude",
            remote_name=REMOTE_DOT_CLAUDE,
            excludes=() if full else DEFAULT_EXCLUDES["claude"],
        ),
        AgentMirror(
            agent_id="claude_app",
            label="Claude Desktop (Application Support)",
            local=h / "Library/Application Support/Claude",
            remote_name=REMOTE_APP_SUPPORT_CLAUDE,
            excludes=() if full else DEFAULT_EXCLUDES["claude_app"],
        ),
        AgentMirror(
            agent_id="claude_science",
            label="Claude Science (~/.claude-science, selective)",
            local=h / ".claude-science",
            remote_name=REMOTE_DOT_CLAUDE_SCIENCE,
            excludes=() if full else DEFAULT_EXCLUDES["claude_science"],
        ),
        AgentMirror(
            agent_id="codex",
            label="Codex (~/.codex)",
            local=h / ".codex",
            remote_name=REMOTE_DOT_CODEX,
            excludes=() if full else DEFAULT_EXCLUDES["codex"],
        ),
        AgentMirror(
            agent_id="grok",
            label="Grok (~/.grok)",
            local=h / ".grok",
            remote_name=REMOTE_DOT_GROK,
            excludes=() if full else DEFAULT_EXCLUDES["grok"],
        ),
    ]
    if agents is None:
        return all_specs
    wanted = {a.strip().lower().replace("-", "_") for a in agents if a.strip()}
    # Accept aliases
    aliases = {
        "claude_code": "claude",
        "dot_claude": "claude",
        "application_support_claude": "claude_app",
        "claude_desktop": "claude_app",
        "science": "claude_science",
        "claude-science": "claude_science",
    }
    normalized = {aliases.get(a, a) for a in wanted}
    return [m for m in all_specs if m.agent_id in normalized]


def print_map(
    mirrors: Iterable[AgentMirror],
    dest_root: Path,
    mode: str = "backup",
) -> list[dict[str, object]]:
    """Return serializable source→dest pairs for backup or restore."""
    rows: list[dict[str, object]] = []
    for m in mirrors:
        remote = m.remote(dest_root)
        if mode == "restore":
            src, dst = remote, m.local
        else:
            src, dst = m.local, remote
        rows.append(
            {
                "agent_id": m.agent_id,
                "label": m.label,
                "mode": mode,
                "source": str(src),
                "dest": str(dst),
                "excludes": list(m.excludes),
                "source_exists": src.exists() if mode == "backup" else None,
            }
        )
    return rows


def restore_rsync_commands(
    archive_batch: Path,
    home: Path | None = None,
    agents: Sequence[str] | None = None,
) -> list[str]:
    """Exact rsync lines matching the user's seed restore direction (iCloud→local).

    Seed pattern:
      A="…/archive_batch_…"
      rsync -a "$A/dot_claude/" "$HOME/.claude/"
      rsync -a "$A/Application_Support_Claude/" \\
          "$HOME/Library/Application Support/Claude/"
    """
    h = home if home is not None else Path.home()
    mirrors = build_agent_mirrors(home=h, agents=agents, full=True)
    cmds: list[str] = []
    # Prefer shell-friendly quoting; keep spaces in Application Support path.
    for m in mirrors:
        src = archive_batch / m.remote_name
        # Only emit if this remote name is part of the requested agents
        src_q = str(src)
        dst_q = str(m.local)
        cmds.append(f'rsync -a "{src_q}/" "{dst_q}/"')
    return cmds


def seed_restore_example(archive_batch: Path | str) -> str:
    """Human-readable seed restore block (Claude only) matching user intent."""
    a = str(archive_batch)
    return (
        f'A="{a}"\n'
        f'rsync -a "$A/{REMOTE_DOT_CLAUDE}/" "$HOME/.claude/"\n'
        f'rsync -a "$A/{REMOTE_APP_SUPPORT_CLAUDE}/" '
        f'"$HOME/Library/Application Support/Claude/"\n'
    )


def is_clouddocs_symlink(path: Path) -> bool:
    """True if path is a symlink whose target sits under Mobile Documents."""
    if not path.is_symlink():
        return False
    try:
        target = os.readlink(path)
    except OSError:
        return False
    return "Mobile Documents" in target or "com~apple~CloudDocs" in target


def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--print-map",
        action="store_true",
        help="Print JSON source/dest pairs for backup (local→iCloud).",
    )
    p.add_argument(
        "--print-restore-cmds",
        action="store_true",
        help="Print rsync restore commands (iCloud archive → local homes).",
    )
    p.add_argument(
        "--print-seed-restore",
        action="store_true",
        help="Print Claude-only seed restore block matching user archive commands.",
    )
    p.add_argument(
        "--print-excludes",
        metavar="AGENT",
        help="Print default thin excludes for one agent_id.",
    )
    p.add_argument(
        "--print-inventory",
        action="store_true",
        help="Print named + home-dot PRESENT/MISSING/SKIP inventory as JSON.",
    )
    p.add_argument(
        "--print-home-dots-map",
        action="store_true",
        help="Print JSON source/dest pairs for allowlisted home dots.",
    )
    p.add_argument(
        "--print-freeable",
        action="store_true",
        help="Print regenerable free candidates (verify mirror before free).",
    )
    p.add_argument(
        "--mode",
        choices=("backup", "restore"),
        default="backup",
        help="Direction for --print-map (default backup = local→iCloud).",
    )
    p.add_argument(
        "--agents",
        default="all",
        help="Comma-separated agent ids or 'all' "
        "(claude,claude_app,claude_science,codex,grok).",
    )
    p.add_argument(
        "--full",
        action="store_true",
        help="Include caches/conda/vm_bundles (no thin excludes).",
    )
    p.add_argument(
        "--icloud",
        default=None,
        help="Override FLEXAIDDS_ICLOUD root.",
    )
    p.add_argument(
        "--archive-batch",
        default=None,
        help="Archive batch path or name under archived_from_ssd/ for restore.",
    )
    p.add_argument(
        "--dest-root",
        default=None,
        help="Override remote root (default: $FLEXAIDDS_ICLOUD/agent_homes).",
    )
    p.add_argument("--home", default=None, help="Override $HOME for tests.")
    args = p.parse_args(list(argv) if argv is not None else None)

    home = Path(args.home) if args.home else Path.home()
    icloud = Path(args.icloud) if args.icloud else default_icloud_root(home)

    agent_list: Sequence[str] | None
    if args.agents.strip().lower() == "all":
        agent_list = None
    else:
        agent_list = [x.strip() for x in args.agents.split(",") if x.strip()]

    if args.print_excludes:
        key = args.print_excludes.strip().lower().replace("-", "_")
        if key not in DEFAULT_EXCLUDES:
            print(f"unknown agent: {key}", file=os.sys.stderr)
            return 2
        for ex in DEFAULT_EXCLUDES[key]:
            print(ex)
        return 0

    if args.print_inventory:
        rows = inventory_home_dots(home=home)
        print(
            json.dumps(
                {
                    "home": str(home),
                    "icloud": str(icloud),
                    "items": [
                        {
                            "name": r.name,
                            "local": str(r.local),
                            "status": r.status,
                            "reason": r.reason,
                            "kind": r.kind,
                            "remote_name": r.remote_name,
                        }
                        for r in rows
                    ],
                },
                indent=2,
            )
        )
        return 0

    if args.print_home_dots_map:
        dest = (
            Path(args.dest_root)
            if args.dest_root
            else home_dots_mirror_root(icloud)
        )
        mirrors = build_home_dot_mirrors(home=home)
        rows = print_map(mirrors, dest_root=dest, mode="backup")
        # Files use source path without trailing slash semantics in rsync
        print(
            json.dumps(
                {
                    "mode": "backup",
                    "dest_root": str(dest),
                    "icloud": str(icloud),
                    "pairs": rows,
                },
                indent=2,
            )
        )
        return 0

    if args.print_freeable:
        print(
            json.dumps(
                {
                    "freeable": freeable_regenerable_paths(
                        home=home, agents=agent_list
                    ),
                    "note": (
                        "Free only after durable agent_homes mirror for "
                        "that agent is verified non-empty; never free whole "
                        "live agent roots."
                    ),
                },
                indent=2,
            )
        )
        return 0

    # Resolve dest / archive root
    if args.archive_batch:
        ab = Path(args.archive_batch)
        if not ab.is_absolute():
            ab = archive_batch_root(icloud, batch=str(ab))
        dest_root = ab
    elif args.dest_root:
        dest_root = Path(args.dest_root)
    else:
        dest_root = agent_homes_mirror_root(icloud)

    if args.print_seed_restore:
        batch = dest_root if args.archive_batch else (
            archive_batch_root(icloud) / "archive_batch_LATEST"
        )
        if args.archive_batch:
            batch = dest_root
        print(seed_restore_example(batch), end="")
        return 0

    if args.print_restore_cmds:
        for line in restore_rsync_commands(
            archive_batch=dest_root, home=home, agents=agent_list
        ):
            print(line)
        return 0

    if args.print_map:
        mirrors = build_agent_mirrors(home=home, agents=agent_list, full=args.full)
        rows = print_map(mirrors, dest_root=dest_root, mode=args.mode)
        print(
            json.dumps(
                {
                    "mode": args.mode,
                    "dest_root": str(dest_root),
                    "icloud": str(icloud),
                    "full": args.full,
                    "pairs": rows,
                },
                indent=2,
            )
        )
        return 0

    p.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
