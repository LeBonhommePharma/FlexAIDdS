"""Unit tests for scripts/agent_icloud_paths.py (local↔iCloud agent home mapping).

Pure path logic — no real iCloud I/O, no mocked rsync-as-pass.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import agent_icloud_paths as aip  # noqa: E402


def test_default_icloud_root_uses_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    custom = tmp_path / "CloudDocs" / "FlexAIDdS_benchmarks"
    monkeypatch.setenv("FLEXAIDDS_ICLOUD", str(custom))
    assert aip.default_icloud_root() == custom


def test_default_icloud_root_falls_back_to_clouddocs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("FLEXAIDDS_ICLOUD", raising=False)
    got = aip.default_icloud_root(home=tmp_path)
    assert got == tmp_path / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
    assert "CloudDocs" in str(got)
    assert "FlexAIDdS_benchmarks" in str(got)


def test_agent_homes_mirror_root():
    icloud = Path("/tmp/fake_icloud")
    assert aip.agent_homes_mirror_root(icloud) == icloud / "agent_homes"


def test_build_agent_mirrors_all_five(tmp_path: Path):
    mirrors = aip.build_agent_mirrors(home=tmp_path)
    ids = [m.agent_id for m in mirrors]
    assert ids == ["claude", "claude_app", "claude_science", "codex", "grok"]
    by_id = {m.agent_id: m for m in mirrors}
    assert by_id["claude"].local == tmp_path / ".claude"
    assert by_id["claude"].remote_name == "dot_claude"
    assert by_id["claude_app"].local == tmp_path / "Library/Application Support/Claude"
    assert by_id["claude_app"].remote_name == "Application_Support_Claude"
    assert by_id["claude_science"].local == tmp_path / ".claude-science"
    assert by_id["codex"].local == tmp_path / ".codex"
    assert by_id["grok"].local == tmp_path / ".grok"


def test_thin_excludes_omit_conda_and_vm_bundles(tmp_path: Path):
    mirrors = aip.build_agent_mirrors(home=tmp_path, full=False)
    by_id = {m.agent_id: m for m in mirrors}
    assert any("conda/" in e for e in by_id["claude_science"].excludes)
    assert any("vm_bundles/" in e for e in by_id["claude_app"].excludes)
    assert any("cache/" in e for e in by_id["claude"].excludes)


def test_full_mode_clears_excludes(tmp_path: Path):
    mirrors = aip.build_agent_mirrors(home=tmp_path, full=True)
    for m in mirrors:
        assert m.excludes == ()


def test_filter_agents_aliases(tmp_path: Path):
    mirrors = aip.build_agent_mirrors(
        home=tmp_path, agents=["claude", "science", "codex"]
    )
    assert [m.agent_id for m in mirrors] == ["claude", "claude_science", "codex"]


def test_print_map_backup_direction(tmp_path: Path):
    dest = tmp_path / "icloud" / "agent_homes"
    mirrors = aip.build_agent_mirrors(home=tmp_path / "home", agents=["claude", "grok"])
    rows = aip.print_map(mirrors, dest_root=dest, mode="backup")
    assert len(rows) == 2
    assert rows[0]["source"] == str(tmp_path / "home" / ".claude")
    assert rows[0]["dest"] == str(dest / "dot_claude")
    assert rows[1]["source"] == str(tmp_path / "home" / ".grok")
    assert rows[1]["dest"] == str(dest / "dot_grok")


def test_print_map_restore_direction(tmp_path: Path):
    dest = tmp_path / "archive_batch"
    mirrors = aip.build_agent_mirrors(home=tmp_path / "home", agents=["claude_app"])
    rows = aip.print_map(mirrors, dest_root=dest, mode="restore")
    assert rows[0]["source"] == str(dest / "Application_Support_Claude")
    assert rows[0]["dest"] == str(
        tmp_path / "home" / "Library/Application Support/Claude"
    )


def test_restore_rsync_commands_match_seed_layout(tmp_path: Path):
    batch = tmp_path / "archive_batch_20260725T095624Z"
    cmds = aip.restore_rsync_commands(
        archive_batch=batch,
        home=tmp_path / "home",
        agents=["claude", "claude_app"],
    )
    assert len(cmds) == 2
    assert cmds[0] == (
        f'rsync -a "{batch}/dot_claude/" "{tmp_path / "home" / ".claude"}/"'
    )
    assert cmds[1] == (
        f'rsync -a "{batch}/Application_Support_Claude/" '
        f'"{tmp_path / "home" / "Library/Application Support/Claude"}/"'
    )


def test_seed_restore_example_contains_both_claude_paths(tmp_path: Path):
    batch = tmp_path / "archived_from_ssd" / "archive_batch_X"
    text = aip.seed_restore_example(batch)
    assert 'A="' in text
    assert "dot_claude" in text
    assert "Application_Support_Claude" in text
    assert "$HOME/.claude/" in text
    assert 'Application Support/Claude/"' in text
    assert "rsync -a" in text


def test_cli_print_map_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys):
    monkeypatch.setenv("FLEXAIDDS_ICLOUD", str(tmp_path / "icloud"))
    rc = aip.main(
        [
            "--print-map",
            "--home",
            str(tmp_path / "h"),
            "--agents",
            "claude,codex",
            "--icloud",
            str(tmp_path / "icloud"),
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    data = json.loads(out)
    assert data["mode"] == "backup"
    assert data["dest_root"] == str(tmp_path / "icloud" / "agent_homes")
    ids = [p["agent_id"] for p in data["pairs"]]
    assert ids == ["claude", "codex"]
    # Destination must be under CloudDocs-style icloud root, not live HOME
    for pair in data["pairs"]:
        assert str(tmp_path / "icloud") in pair["dest"]
        assert not pair["dest"].startswith(str(tmp_path / "h"))


def test_cli_print_seed_restore(tmp_path: Path, capsys):
    batch = tmp_path / "archive_batch_20260725T095624Z"
    rc = aip.main(["--print-seed-restore", "--archive-batch", str(batch)])
    assert rc == 0
    out = capsys.readouterr().out
    assert f'A="{batch}"' in out
    assert 'rsync -a "$A/dot_claude/" "$HOME/.claude/"' in out
    assert "Application_Support_Claude" in out


def test_is_clouddocs_symlink_false_for_real_dir(tmp_path: Path):
    d = tmp_path / ".claude"
    d.mkdir()
    assert aip.is_clouddocs_symlink(d) is False


def test_is_clouddocs_symlink_true_for_mobile_documents(tmp_path: Path):
    target = tmp_path / "Library/Mobile Documents/com~apple~CloudDocs/x"
    target.mkdir(parents=True)
    link = tmp_path / ".claude"
    link.symlink_to(target)
    assert aip.is_clouddocs_symlink(link) is True


def test_sync_script_exists_and_is_executable():
    script = SCRIPTS / "sync_agent_homes_to_icloud.sh"
    assert script.is_file()
    text = script.read_text(encoding="utf-8")
    assert "sync_agent_homes_to_icloud" in text
    assert "--dry-run" in text
    assert "--restore" in text
    assert "--backup" in text
    assert "dot_claude" in text
    assert "Application_Support_Claude" in text
    assert "claude_science" in text or "claude-science" in text
    assert "codex" in text
    assert "grok" in text
    # Safety: no find over CloudDocs
    assert "find " not in text or "Never" in text
    assert "gtimeout" in text or "timeout" in text


def test_inventory_named_paths_missing_venv_env(tmp_path: Path):
    (tmp_path / ".claude").mkdir()
    (tmp_path / ".codex").mkdir()
    (tmp_path / ".grok").mkdir()
    # no .venv / .env
    rows = {r.name: r for r in aip.inventory_named_home_paths(home=tmp_path)}
    assert rows[".claude"].status == "PRESENT"
    assert rows[".codex"].status == "PRESENT"
    assert rows[".grok"].status == "PRESENT"
    assert rows[".venv"].status == "MISSING"
    assert rows[".env"].status == "MISSING"
    assert "nothing to mirror" in rows[".venv"].reason


def test_home_dot_allowlist_mirrors_and_skips_zcompdump(tmp_path: Path):
    (tmp_path / ".zshrc").write_text("export X=1\n")
    (tmp_path / ".flexaidds_env").write_text("FLEXAIDDS_ROOT=/x\n")
    (tmp_path / ".zcompdump-foo").write_text("noise\n")
    (tmp_path / ".cursor").mkdir()
    mirrors = aip.build_home_dot_mirrors(home=tmp_path)
    names = {m.local.name for m in mirrors}
    assert ".zshrc" in names
    assert ".flexaidds_env" in names
    assert ".cursor" in names
    assert ".zcompdump-foo" not in names
    inv = aip.inventory_home_dots(home=tmp_path)
    skips = [r for r in inv if r.status == "SKIP" and "zcompdump" in r.name]
    assert skips
    dest = tmp_path / "icloud" / "home_dots"
    rows = aip.print_map(mirrors, dest_root=dest, mode="backup")
    assert any(r["dest"].endswith("dot_zshrc") for r in rows)


def test_freeable_never_includes_agent_root(tmp_path: Path):
    (tmp_path / ".claude" / "cache").mkdir(parents=True)
    (tmp_path / ".codex" / "cache").mkdir(parents=True)
    rows = aip.freeable_regenerable_paths(home=tmp_path, agents=["claude", "codex"])
    for row in rows:
        assert row["local"] != row["agent_local_root"]
        assert row["local"].startswith(row["agent_local_root"])
    assert any(r["relative"] == "cache" and r["agent_id"] == "claude" for r in rows)


def test_freeable_excludes_unique_install_and_user_trees(tmp_path: Path):
    """bin/vendor/attachments/generated_images/blob_storage must never be freeable."""
    for agent, sub in (
        (".grok", "bin"),
        (".grok", "vendor"),
        (".codex", "attachments"),
        (".codex", "generated_images"),
        (".claude", "cache"),
    ):
        (tmp_path / agent / sub).mkdir(parents=True, exist_ok=True)
    rows = aip.freeable_regenerable_paths(
        home=tmp_path, agents=["claude", "codex", "grok"]
    )
    relatives = {(r["agent_id"], r["relative"]) for r in rows}
    assert ("grok", "bin") not in relatives
    assert ("grok", "vendor") not in relatives
    assert ("codex", "attachments") not in relatives
    assert ("codex", "generated_images") not in relatives
    assert ("claude", "cache") in relatives
    for bad in aip.NEVER_FREE_RELATIVE:
        assert aip.is_never_free_relative(bad)


def test_resolve_free_proof_requires_path_level_duplicate(tmp_path: Path):
    """Parent agent_homes non-empty is not enough — relative path must exist."""
    homes = tmp_path / "agent_homes"
    arch = tmp_path / "archive"
    # Parent has other files but NOT cache/
    (homes / "dot_claude" / "settings.json").parent.mkdir(parents=True)
    (homes / "dot_claude" / "settings.json").write_text("{}")
    assert (
        aip.resolve_free_proof(
            "cache", "dot_claude", agent_homes_root=homes, archive_root=arch
        )
        is None
    )
    # Path-level proof on agent_homes
    (homes / "dot_claude" / "cache").mkdir()
    (homes / "dot_claude" / "cache" / "x").write_text("1")
    proof = aip.resolve_free_proof(
        "cache", "dot_claude", agent_homes_root=homes, archive_root=arch
    )
    assert proof is not None
    assert proof == homes / "dot_claude" / "cache"
    # Never free bin even if remote has it
    (homes / "dot_grok" / "bin").mkdir(parents=True)
    (homes / "dot_grok" / "bin" / "x").write_text("1")
    assert (
        aip.resolve_free_proof(
            "bin", "dot_grok", agent_homes_root=homes, archive_root=None
        )
        is None
    )


def test_freeable_require_remote_proof_filters(tmp_path: Path):
    home = tmp_path / "home"
    (home / ".claude" / "cache").mkdir(parents=True)
    (home / ".codex" / "cache").mkdir(parents=True)
    homes = tmp_path / "agent_homes"
    # Only claude/cache mirrored
    (homes / "dot_claude" / "cache").mkdir(parents=True)
    (homes / "dot_claude" / "cache" / "a").write_text("1")
    rows = aip.freeable_regenerable_paths(
        home=home,
        agents=["claude", "codex"],
        agent_homes_root=homes,
        archive_root=tmp_path / "empty_arch",
        require_remote_proof=True,
    )
    assert len(rows) == 1
    assert rows[0]["agent_id"] == "claude"
    assert rows[0]["relative"] == "cache"
    assert rows[0]["proof_path"].endswith("dot_claude/cache")


def test_freeable_archive_proof_counts(tmp_path: Path):
    home = tmp_path / "home"
    (home / ".claude" / "telemetry").mkdir(parents=True)
    homes = tmp_path / "agent_homes"
    homes.mkdir()
    arch = tmp_path / "archive"
    (arch / "dot_claude" / "telemetry").mkdir(parents=True)
    (arch / "dot_claude" / "telemetry" / "t").write_text("t")
    proof = aip.resolve_free_proof(
        "telemetry",
        "dot_claude",
        agent_homes_root=homes,
        archive_root=arch,
    )
    assert proof == arch / "dot_claude" / "telemetry"


def test_home_dots_and_safe_free_scripts_exist():
    for name in (
        "sync_home_dots_to_icloud.sh",
        "safe_free_verified_icloud_duplicates.sh",
    ):
        script = SCRIPTS / name
        assert script.is_file(), name
        text = script.read_text(encoding="utf-8")
        assert "FLEXAIDDS_ICLOUD" in text
        assert "dry-run" in text or "--dry-run" in text
    free_sh = (SCRIPTS / "safe_free_verified_icloud_duplicates.sh").read_text(
        encoding="utf-8"
    )
    assert "path_level_proof" in free_sh
    assert "never-free" in free_sh or "REFUSE never-free" in free_sh
    assert "bin|vendor|attachments" in free_sh
