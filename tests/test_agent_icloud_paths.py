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
