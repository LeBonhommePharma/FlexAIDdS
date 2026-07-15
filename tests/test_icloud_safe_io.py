"""Unit tests for scripts/icloud_safe_io.py (local-first / thin-iCloud I/O).

No real iCloud required — uses tmp dirs and path-string CloudDocs markers.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import icloud_safe_io as iso  # noqa: E402


def test_is_clouddocs_detects_mobile_documents():
    cloud = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/foo/bar.dat"
    assert iso.is_clouddocs(cloud) is True
    assert iso.is_clouddocs(str(cloud)) is True


def test_is_clouddocs_false_for_local_paths(tmp_path: Path):
    local = tmp_path / "flexaidds_results" / "campaigns" / "x" / "result.csv"
    assert iso.is_clouddocs(local) is False
    assert iso.is_clouddocs("/Users/someone/Projects/FlexAIDdS/WRK/file.dat") is False
    assert iso.is_local_apfs(local) is True


def test_safe_md5_local_file(tmp_path: Path):
    p = tmp_path / "matrix.dat"
    data = b"MC_st0r5.2_6 mock matrix payload\n" * 10
    p.write_bytes(data)
    expected = hashlib.md5(data).hexdigest()
    got = iso.safe_md5(p)
    assert got == expected


def test_materialize_local_returns_same_path(tmp_path: Path):
    p = tmp_path / "local_only.dat"
    p.write_text("hello local")
    out = iso.materialize(p)
    assert out is not None
    assert out.resolve() == p.resolve()


def test_materialize_missing_local_returns_none(tmp_path: Path):
    p = tmp_path / "does_not_exist.dat"
    assert iso.materialize(p) is None


def test_prefer_local_campaign_prefers_local_with_result_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    local_root = tmp_path / "flexaidds_results"
    campaigns = local_root / "campaigns"
    camp = campaigns / "C0_test_campaign"
    target = camp / "1ABC"
    target.mkdir(parents=True)
    (target / "result.csv").write_text("target,rmsd\n1ABC,1.0\n")

    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(local_root))
    # Reset module cache of env via direct args
    chosen = iso.prefer_local_campaign(
        "C0_test_campaign",
        local_campaigns=campaigns,
        icloud_results=tmp_path / "icloud" / "results" / "campaigns",
    )
    assert chosen == camp


def test_prefer_local_campaign_falls_back_to_empty_local_dir(
    tmp_path: Path,
):
    campaigns = tmp_path / "campaigns"
    camp = campaigns / "empty_camp"
    camp.mkdir(parents=True)
    chosen = iso.prefer_local_campaign(
        "empty_camp",
        local_campaigns=campaigns,
        icloud_results=tmp_path / "icloud",
    )
    assert chosen == camp


def test_prefer_local_campaign_missing_returns_local_path(tmp_path: Path):
    campaigns = tmp_path / "campaigns"
    campaigns.mkdir()
    chosen = iso.prefer_local_campaign(
        "never_created",
        local_campaigns=campaigns,
    )
    assert chosen == campaigns / "never_created"


def test_safe_glob_result_csvs_one_level(tmp_path: Path):
    camp = tmp_path / "camp"
    (camp / "t1").mkdir(parents=True)
    (camp / "t2").mkdir(parents=True)
    (camp / "t1" / "result.csv").write_text("a\n")
    (camp / "t2" / "result.csv").write_text("b\n")
    # nested should NOT be found (one-level only)
    deep = camp / "t1" / "nested"
    deep.mkdir()
    (deep / "result.csv").write_text("deep\n")
    # incomplete archive skipped
    inc = camp / "t3_incomplete_backup"
    inc.mkdir()
    (inc / "result.csv").write_text("skip\n")

    found = iso.safe_glob_result_csvs(camp)
    names = sorted(p.parent.name for p in found)
    assert names == ["t1", "t2"]


def test_local_root_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(tmp_path / "custom_local"))
    assert iso.local_root() == tmp_path / "custom_local"


def test_pin_cache_dir_under_local(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    monkeypatch.setenv("FLEXAIDDS_LOCAL_ROOT", str(tmp_path / "lr"))
    monkeypatch.delenv("FLEXAIDDS_PIN_CACHE", raising=False)
    d = iso.pin_cache_dir()
    assert d == tmp_path / "lr" / "pins" / "materialize"
    assert d.is_dir()
